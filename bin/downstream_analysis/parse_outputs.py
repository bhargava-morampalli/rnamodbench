#!/usr/bin/env python3
"""
Output Parser Module for RNA Modification Detection Tools

Parses outputs from 9 different modification detection tools into a
standardized format for downstream analysis.

Supported tools:
- TOMBO: CSV with p-values and KS statistics
- Yanocomp: BED + JSON with FDR and KS statistics
- Nanocompore: TSV with multiple p-value types
- xPore: diffmod.table with modification rates and p-values
- ELIGOS: TXT with p-values and odds ratios
- EpiNano: CSV with z-scores and error differences
- DiffErr: BED with FDR, p-values, and odds ratios
- DRUMMER: summary.txt with p-values and odds ratios
- JACUSA2: BED with scores

Output format (standardized DataFrame):
    - tool: str (tool name)
    - reference: str (reference/chromosome name)
    - position: int (1-based position)
    - score: float (primary score for ranking)
    - score_type: str (pvalue, fdr, zscore, probability, score)
    - pvalue: float (p-value if available, else NaN)
    - replicate: str (replicate identifier)
    - sample_id: str (full sample identifier)
    - coverage: str (coverage level, e.g., "100x", "500x", or NaN if not applicable)
"""

import os
import json
import logging
import re
from pathlib import Path
from typing import Optional, Dict, List, Tuple, Union

import pandas as pd
import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

DIR_CAPABLE_TOOLS = {"nanocompore", "xpore", "eligos", "epinano", "drummer"}
SUPPORTED_TOOLS = [
    "tombo",
    "yanocomp",
    "nanocompore",
    "xpore",
    "eligos",
    "epinano",
    "differr",
    "drummer",
    "jacusa2",
    "nanodoc",
]


def _empty_dataframe() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "tool",
            "reference",
            "position",
            "score",
            "score_type",
            "pvalue",
            "replicate",
            "sample_id",
            "coverage",
            "_imputed",
        ]
    )


def _safe_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _find_column(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    """Find best matching column by exact or case-insensitive name."""
    lower_map = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c in df.columns:
            return c
        if c.lower() in lower_map:
            return lower_map[c.lower()]
    return None


def _extract_reference_from_path(filepath: Path) -> str:
    """Best-effort reference extraction from filename/path tokens."""
    tokens = [filepath.stem] + [p.name for p in list(filepath.parents)[:4]]
    for token in tokens:
        token_l = token.lower()
        if "16s" in token_l:
            return "16s"
        if "23s" in token_l:
            return "23s"
        if "5s" in token_l:
            return "5s"
    return "unknown"


def _extract_replicate_from_path(filepath: Path) -> str:
    """Extract replicate from full path (filename + parent directories)."""
    patterns = [
        r"rep(?:licate)?[_-]?(\d+)",
        r"(?:^|[_-])r(\d+)(?:$|[_-])",
    ]

    candidates = [filepath.stem]
    candidates.extend([p.name for p in list(filepath.parents)[:6]])

    for text in candidates:
        for pat in patterns:
            m = re.search(pat, text, flags=re.IGNORECASE)
            if m:
                return f"rep{m.group(1)}"
    return "rep1"


def _extract_coverage_from_path(filepath: Path) -> Optional[str]:
    for part in filepath.parts:
        m = re.search(r"(?:coverage[_-]?)?(\d+)[xX]", part)
        if m:
            return f"{m.group(1)}x"
        if re.fullmatch(r"\d+", part):
            return f"{part}x"
    return None


def _read_auto_delim(filepath: Path) -> pd.DataFrame:
    """Read CSV/TSV where delimiter may vary."""
    try:
        return pd.read_csv(filepath, sep=None, engine="python")
    except Exception:
        # fallback to comma first, then tab
        try:
            return pd.read_csv(filepath)
        except Exception:
            return pd.read_csv(filepath, sep="\t")


def _benjamini_hochberg(pvalues: pd.Series) -> pd.Series:
    """Return Benjamini-Hochberg adjusted p-values with NaNs preserved."""
    vals = pd.to_numeric(pvalues, errors="coerce").astype(float)
    arr = vals.to_numpy(copy=True)
    valid = np.isfinite(arr)
    if not valid.any():
        return pd.Series(np.nan, index=pvalues.index, dtype=float)

    clipped = np.clip(arr[valid], 0.0, 1.0)
    order = np.argsort(clipped)
    ranked = clipped[order]
    n = ranked.size
    ranks = np.arange(1, n + 1, dtype=float)
    adjusted = ranked * n / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)

    out = np.full(arr.shape, np.nan, dtype=float)
    valid_idx = np.flatnonzero(valid)
    out[valid_idx[order]] = adjusted
    return pd.Series(out, index=pvalues.index, dtype=float)


def _parse_epinano_chr_pos(value: str) -> Tuple[str, Optional[int]]:
    """Parse chr_pos like '16s_88_rrsE 100 A +' -> (reference, position)."""
    if pd.isna(value):
        return "unknown", None

    tokens = str(value).strip().split()
    if len(tokens) < 2:
        return "unknown", None

    reference = tokens[0]
    try:
        pos = int(tokens[1])
    except Exception:
        pos = None
    return reference, pos


def _finalize_result(
    df: pd.DataFrame,
    tool: str,
    filepath: Path,
    sample_id: Optional[str],
    replicate: Optional[str],
    coverage: Optional[str],
) -> pd.DataFrame:
    if df.empty:
        return _empty_dataframe()

    result = df.copy()
    result["tool"] = tool
    result["replicate"] = replicate or _extract_replicate_from_path(filepath)
    result["sample_id"] = sample_id or filepath.stem
    result["coverage"] = coverage if coverage is not None else _extract_coverage_from_path(filepath)
    result["_imputed"] = False

    # enforce required columns
    for col in ["reference", "position", "score", "score_type", "pvalue"]:
        if col not in result.columns:
            result[col] = np.nan

    result["reference"] = result["reference"].astype(str)
    result["position"] = _safe_numeric(result["position"]).astype("Int64")
    result["score"] = _safe_numeric(result["score"])
    result["pvalue"] = _safe_numeric(result["pvalue"])

    # drop rows without parseable position/reference
    result = result[result["reference"].notna()]
    result = result[result["reference"].astype(str).str.lower() != "unknown"]
    result = result[result["position"].notna()]

    if result.empty:
        return _empty_dataframe()

    # cast positions back to int after filtering
    result["position"] = result["position"].astype(int)

    keep_cols = [
        "tool",
        "reference",
        "position",
        "score",
        "score_type",
        "pvalue",
        "replicate",
        "sample_id",
        "coverage",
        "_imputed",
    ]

    # keep extra context columns if present
    extra_cols = [
        c
        for c in [
            "bed_score_rounded",
            "strand",
            "odds_ratio",
            "g_stat",
            "g_pval",
            "g_fdr",
            "g_pval_neglog10",
            "g_fdr_neglog10",
            "or_pval",
            "or_fdr",
            "log_odds_ratio",
            "diff_mod_rate",
            "raw_pvalue",
            "bh_fdr",
            "delta_mis",
            "delta_sum_err",
            "z_scores",
            "z_score_prediction",
            "max.G_padj",
            "OR_padj",
            "G_padj",
        ]
        if c in result.columns
    ]

    return result[keep_cols + extra_cols].reset_index(drop=True)


# -----------------------------------------------------------------------------
# Tool-specific parsers
# -----------------------------------------------------------------------------

def parse_tombo(
    filepath: Union[str, Path],
    sample_id: str = None,
    replicate: str = None
) -> pd.DataFrame:
    """
    Parse TOMBO CSV output.

    TOMBO outputs from tombo_stats.LevelStats contain:
    - pos: position
    - stat: KS statistic or other test statistic
    - valid_cov: valid coverage
    - p-value or similar statistics depending on test type

    Args:
        filepath: Path to TOMBO CSV file
        sample_id: Sample identifier
        replicate: Replicate identifier

    Returns:
        Standardized DataFrame
    """
    filepath = Path(filepath)
    logger.info(f"Parsing TOMBO output: {filepath}")

    if not filepath.exists():
        logger.warning(f"File not found: {filepath}")
        return _empty_dataframe()

    try:
        df = pd.read_csv(filepath)

        if df.empty:
            logger.warning(f"Empty TOMBO file: {filepath}")
            return _empty_dataframe()

        # Identify column names (TOMBO can have different column names)
        pos_col = _find_column(df, ['pos', 'position', 'start'])
        stat_col = _find_column(df, ['stat', 'statistic', 'ks_stat', 'p_value', 'pvalue'])

        if pos_col is None:
            logger.error(f"Could not find position column in TOMBO output: {df.columns.tolist()}")
            return _empty_dataframe()

        # Extract reference name from filename if not in data
        ref_name = _extract_reference_from_filename(filepath)

        result = pd.DataFrame({
            'tool': 'tombo',
            'reference': ref_name,
            'position': df[pos_col].astype(int),
            'score': df[stat_col].astype(float) if stat_col else np.nan,
            'score_type': 'pvalue',
            'pvalue': df[stat_col].astype(float) if stat_col else np.nan,
            'replicate': replicate or _extract_replicate_from_filename(filepath),
            'sample_id': sample_id or filepath.stem
        })

        logger.info(f"Parsed {len(result)} positions from TOMBO")
        return result

    except Exception as e:
        logger.error(f"Error parsing TOMBO file {filepath}: {e}")
        return _empty_dataframe()


def parse_yanocomp(
    bed_filepath: Union[str, Path],
    json_filepath: Optional[Union[str, Path]] = None,
    sample_id: str = None,
    replicate: str = None
) -> pd.DataFrame:
    """
    Parse Yanocomp BED output (and optionally JSON predictions).

    Yanocomp outputs:
    - BED file: chr, start, end, name, score (FDR), strand
    - JSON file: detailed predictions with KS statistics

    Args:
        bed_filepath: Path to Yanocomp BED file
        json_filepath: Optional path to JSON predictions file
        sample_id: Sample identifier
        replicate: Replicate identifier

    Returns:
        Standardized DataFrame
    """
    bed_filepath = Path(bed_filepath)
    logger.info(f"Parsing Yanocomp output: {bed_filepath}")

    if not bed_filepath.exists():
        logger.warning(f"File not found: {bed_filepath}")
        return _empty_dataframe()

    try:
        # Parse BED file (no header)
        df = pd.read_csv(
            bed_filepath,
            sep='\t',
            header=None,
            names=['chr', 'start', 'end', 'name', 'fdr', 'strand'],
            comment='#'
        )

        if df.empty:
            logger.warning(f"Empty Yanocomp BED file: {bed_filepath}")
            return _empty_dataframe()

        # Try to parse JSON for additional info (KS statistic)
        ks_stats = {}
        if json_filepath:
            json_filepath = Path(json_filepath)
            if json_filepath.exists():
                try:
                    with open(json_filepath, 'r') as f:
                        json_data = json.load(f)
                    # Extract KS statistics from JSON
                    for key, value in json_data.items():
                        if isinstance(value, dict) and 'ks_stat' in value:
                            ks_stats[key] = value['ks_stat']
                except Exception as e:
                    logger.warning(f"Could not parse Yanocomp JSON: {e}")

        result = pd.DataFrame({
            'tool': 'yanocomp',
            'reference': df['chr'],
            'position': (df['start'] + 1).astype(int),  # Convert to 1-based
            'score': df['fdr'].astype(float),
            'score_type': 'fdr',
            'pvalue': np.nan,  # FDR is not p-value
            'replicate': replicate or _extract_replicate_from_filename(bed_filepath),
            'sample_id': sample_id or bed_filepath.stem
        })
        return result

    except Exception as e:
        logger.error(f"Error parsing Yanocomp file {bed_filepath}: {e}")
        return _empty_dataframe()


def parse_nanocompore(
    filepath: Union[str, Path],
    sample_id: str = None,
    replicate: str = None,
    pvalue_col: str = 'GMM_logit_pvalue'
) -> pd.DataFrame:
    filepath = Path(filepath)

    if filepath.is_dir():
        candidates = [
            filepath / "outnanocompore_results.tsv",
            filepath / "outSampComp_results.tsv",
            filepath / "sampcomp_results.tsv",
            filepath / "results.tsv",
        ]
        filepath = next((p for p in candidates if p.exists()), filepath / "outnanocompore_results.tsv")

    logger.info("Parsing Nanocompore output: %s", filepath)
    if not filepath.exists():
        return _empty_dataframe()

    try:
        df = pd.read_csv(filepath, sep="\t")
        if df.empty:
            return _empty_dataframe()

        pos_col = _find_column(df, ["genomicPos", "pos", "ref_pos", "position"])
        raw_pos_col = _find_column(df, ["pos", "ref_pos", "position"])
        ref_col = _find_column(df, ["ref_id", "reference", "chr", "chrom", "contig"])
        pval_col = _find_column(df, [pvalue_col, "GMM_logit_pvalue", "KS_dwell_pvalue", "pvalue"])

        if pos_col is None:
            return _empty_dataframe()

        pos_vals = _safe_numeric(df[pos_col])
        if raw_pos_col is not None:
            fallback_pos = _safe_numeric(df[raw_pos_col])
            if pos_col.lower() == "genomicpos":
                pos_vals = pos_vals.where(pos_vals.notna(), fallback_pos) + 2
            else:
                pos_vals = pos_vals + 2
        else:
            pos_vals = pos_vals + 2

        out = pd.DataFrame(
            {
                "reference": df[ref_col] if ref_col else _extract_reference_from_path(filepath),
                "position": pos_vals,
                "score": _safe_numeric(df[pval_col]) if pval_col else np.nan,
                "score_type": "pvalue",
                "pvalue": _safe_numeric(df[pval_col]) if pval_col else np.nan,
            }
        )
        return _finalize_result(out, "nanocompore", filepath, sample_id, replicate, None)
    except Exception as exc:
        logger.error("Error parsing Nanocompore file %s: %s", filepath, exc)
        return _empty_dataframe()


def parse_xpore(
    filepath: Union[str, Path],
    sample_id: str = None,
    replicate: str = None
) -> pd.DataFrame:
    filepath = Path(filepath)

    if filepath.is_dir():
        diffmod_file = filepath / "diffmod.table"
        if not diffmod_file.exists():
            tables = sorted(filepath.glob("*.table"))
            if tables:
                diffmod_file = tables[0]
        filepath = diffmod_file

    logger.info("Parsing xPore output: %s", filepath)
    if not filepath.exists():
        return _empty_dataframe()

    try:
        df = _read_auto_delim(filepath)
        if df.empty:
            return _empty_dataframe()

        pos_col = _find_column(df, ["position", "pos"])
        ref_col = _find_column(df, ["id", "ref_id", "reference", "chr"])
        pval_cols = [c for c in df.columns if c.startswith("pval_")]
        pval_col = pval_cols[0] if pval_cols else _find_column(df, ["pvalue", "p_value"])

        if pos_col is None:
            return _empty_dataframe()

        raw_pvalue = _safe_numeric(df[pval_col]) if pval_col else pd.Series(np.nan, index=df.index)
        bh_fdr = _benjamini_hochberg(raw_pvalue)

        out = pd.DataFrame(
            {
                "reference": df[ref_col] if ref_col else _extract_reference_from_path(filepath),
                "position": _safe_numeric(df[pos_col]) + 2,
                "score": bh_fdr,
                "score_type": "fdr",
                "pvalue": raw_pvalue,
                "raw_pvalue": raw_pvalue,
                "bh_fdr": bh_fdr,
            }
        )

        rate_cols = [c for c in df.columns if "diff_mod_rate" in c]
        if rate_cols:
            out["diff_mod_rate"] = _safe_numeric(df[rate_cols[0]])

        return _finalize_result(out, "xpore", filepath, sample_id, replicate, None)
    except Exception as exc:
        logger.error("Error parsing xPore file %s: %s", filepath, exc)
        return _empty_dataframe()


def parse_eligos(
    filepath: Union[str, Path],
    sample_id: str = None,
    replicate: str = None
) -> pd.DataFrame:
    filepath = Path(filepath)

    if filepath.is_dir():
        candidates = sorted(filepath.glob("*_combine.txt"))
        if not candidates:
            candidates = sorted(filepath.glob("*result*.txt"))
        filepath = candidates[0] if candidates else filepath / "test_paired_diff_mod_result.txt"

    logger.info("Parsing ELIGOS output: %s", filepath)
    if not filepath.exists():
        return _empty_dataframe()

    try:
        df = pd.read_csv(filepath, sep="\t")
        if df.empty:
            return _empty_dataframe()

        pos_col = _find_column(df, ["start_loc", "start", "pos", "position"])
        ref_col = _find_column(df, ["chrom", "chr", "reference", "contig"])
        score_col = _find_column(df, ["adjPval", "adj_pval", "pval", "pvalue", "p_value"])

        if pos_col is None:
            return _empty_dataframe()

        out = pd.DataFrame(
            {
                "reference": df[ref_col] if ref_col else _extract_reference_from_path(filepath),
                "position": _safe_numeric(df[pos_col]) + 1,
                "score": _safe_numeric(df[score_col]) if score_col else np.nan,
                "score_type": "pvalue",
                "pvalue": _safe_numeric(df[score_col]) if score_col else np.nan,
            }
        )

        odds_col = _find_column(df, ["oddR", "odds_ratio", "oddsratio"])
        if odds_col:
            out["odds_ratio"] = _safe_numeric(df[odds_col])

        return _finalize_result(out, "eligos", filepath, sample_id, replicate, None)
    except Exception as exc:
        logger.error("Error parsing ELIGOS file %s: %s", filepath, exc)
        return _empty_dataframe()


def parse_epinano(
    filepath: Union[str, Path],
    sample_id: str = None,
    replicate: str = None
) -> pd.DataFrame:
    filepath = Path(filepath)

    if filepath.is_dir():
        preferred = sorted(filepath.glob("*delta-sum_err*.prediction.csv"))
        if not preferred:
            preferred = sorted(filepath.glob("*sumerr*.prediction.csv"))
        if not preferred:
            preferred = sorted(filepath.glob("*mismatch*.prediction.csv"))
        if preferred:
            filepath = preferred[0]
        else:
            diff_err_files = sorted(filepath.glob("*diff_err*.csv"))
            if diff_err_files:
                filepath = diff_err_files[0]
            else:
                per_site = [filepath / "native_per_site.csv", filepath / "ivt_per_site.csv"]
                existing_per_site = next((p for p in per_site if p.exists()), None)
                if existing_per_site:
                    filepath = existing_per_site
                else:
                    any_csv = sorted(filepath.glob("*.csv"))
                    filepath = any_csv[0] if any_csv else filepath / "native_per_site.csv"

    logger.info("Parsing EpiNano output: %s", filepath)
    if not filepath.exists():
        return _empty_dataframe()

    try:
        df = pd.read_csv(filepath)
        if df.empty:
            return _empty_dataframe()

        filename_l = filepath.name.lower()
        if "prediction" in filename_l and "chr_pos" in df.columns:
            parsed = df["chr_pos"].apply(_parse_epinano_chr_pos)
            ref_series = parsed.apply(lambda x: x[0])
            pos_series = parsed.apply(lambda x: x[1])

            score_col = _find_column(
                df,
                ["delta_sum_err", "delta_mis", "sum_err", "z_scores", "z_score", "zscore"],
            )
            score_type = "error"
            if score_col and "z" in score_col.lower():
                score_type = "zscore"

            out = pd.DataFrame(
                {
                    "reference": ref_series,
                    "position": _safe_numeric(pos_series),
                    "score": _safe_numeric(df[score_col]) if score_col else np.nan,
                    "score_type": score_type,
                    "pvalue": np.nan,
                    "delta_sum_err": _safe_numeric(df["delta_sum_err"]) if "delta_sum_err" in df.columns else np.nan,
                    "delta_mis": _safe_numeric(df["delta_mis"]) if "delta_mis" in df.columns else np.nan,
                    "z_scores": _safe_numeric(df["z_scores"]) if "z_scores" in df.columns else np.nan,
                    "z_score_prediction": df["z_score_prediction"].astype(str) if "z_score_prediction" in df.columns else np.nan,
                }
            )
        else:
            ref_col = _find_column(df, ["#Ref", "Ref", "X.Ref", "reference", "chr"])
            pos_col = _find_column(df, ["pos", "position", "start"])
            score_col = _find_column(df, ["sum_err", "mis", "error", "z_scores"])

            if ref_col is None or pos_col is None:
                return _empty_dataframe()

            score_type = "error"
            if score_col and "z" in score_col.lower():
                score_type = "zscore"

            out = pd.DataFrame(
                {
                    "reference": df[ref_col].astype(str),
                    "position": _safe_numeric(df[pos_col]) + 1,
                    "score": _safe_numeric(df[score_col]) if score_col else np.nan,
                    "score_type": score_type,
                    "pvalue": np.nan,
                }
            )

        return _finalize_result(out, "epinano", filepath, sample_id, replicate, None)
    except Exception as exc:
        logger.error("Error parsing EpiNano file %s: %s", filepath, exc)
        return _empty_dataframe()


def parse_differr(
    filepath: Union[str, Path],
    sample_id: str = None,
    replicate: str = None,
    score_field: str = "g_fdr_neglog10",
) -> pd.DataFrame:
    """
    Parse DiffErr BED output.

    DiffErr outputs BED format with columns:
    chr, start, end, odds_ratio, g_stat, pval, fdr

    Args:
        filepath: Path to DiffErr BED file
        sample_id: Sample identifier
        replicate: Replicate identifier

    Returns:
        Standardized DataFrame
    """
    filepath = Path(filepath)
    logger.info(f"Parsing DiffErr output: {filepath}")

    if not filepath.exists():
        logger.warning(f"File not found: {filepath}")
        return _empty_dataframe()

    try:
        if score_field not in {"g_fdr_neglog10", "g_stat"}:
            raise ValueError(
                f"Unsupported DiffErr score field: {score_field}. "
                f"Expected one of ['g_fdr_neglog10', 'g_stat']"
            )

        df = pd.read_csv(filepath, sep="\t", comment="#", header=None)
        if df.empty:
            return _empty_dataframe()

        if str(df.iloc[0, 0]).lower() in {"chr", "chrom", "reference"}:
            df = df.iloc[1:].reset_index(drop=True)
            if df.empty:
                return _empty_dataframe()

        n_cols = df.shape[1]

        if n_cols >= 14:
            score_col = 7 if score_field == "g_stat" else 9
            score_type = "g_stat" if score_field == "g_stat" else "neglog10_fdr"
            out = pd.DataFrame(
                {
                    "reference": df.iloc[:, 0].astype(str),
                    "position": _safe_numeric(df.iloc[:, 1]) + 1,
                    "score": _safe_numeric(df.iloc[:, score_col]),
                    "score_type": score_type,
                    "pvalue": np.nan,
                    "bed_score_rounded": _safe_numeric(df.iloc[:, 4]),
                    "strand": df.iloc[:, 5].astype(str),
                    "log_odds_ratio": _safe_numeric(df.iloc[:, 6]),
                    "g_stat": _safe_numeric(df.iloc[:, 7]),
                    "g_pval_neglog10": _safe_numeric(df.iloc[:, 8]),
                    "g_fdr_neglog10": _safe_numeric(df.iloc[:, 9]),
                }
            )
        elif n_cols >= 7:
            score_col = 4 if score_field == "g_stat" else 6
            score_type = "g_stat" if score_field == "g_stat" else "neglog10_fdr"
            out = pd.DataFrame(
                {
                    "reference": df.iloc[:, 0].astype(str),
                    "position": _safe_numeric(df.iloc[:, 1]) + 1,
                    "score": _safe_numeric(df.iloc[:, score_col]),
                    "score_type": score_type,
                    "pvalue": np.nan,
                    "odds_ratio": _safe_numeric(df.iloc[:, 3]),
                    "g_stat": _safe_numeric(df.iloc[:, 4]),
                    "g_pval": _safe_numeric(df.iloc[:, 5]),
                    "g_fdr": _safe_numeric(df.iloc[:, 6]),
                }
            )
        else:
            return _empty_dataframe()

        return _finalize_result(out, "differr", filepath, sample_id, replicate, None)
    except Exception as e:
        logger.error(f"Error parsing DiffErr file {filepath}: {e}")
        return _empty_dataframe()


def parse_drummer(
    filepath: Union[str, Path],
    sample_id: str = None,
    replicate: str = None
) -> pd.DataFrame:
    """
    Parse DRUMMER summary.txt output.

    DRUMMER outputs summary.txt with columns:
    transcript_id, chrom, ref_base, position, read_depth,
    base_fractions, odds_ratio, pval, motif, g_test, genomic_pos

    Args:
        filepath: Path to DRUMMER directory or summary.txt file
        sample_id: Sample identifier
        replicate: Replicate identifier

    Returns:
        Standardized DataFrame
    """
    filepath = Path(filepath)

    # Handle directory input
    if filepath.is_dir():
        summary_file = filepath / 'summary.txt'
        if not summary_file.exists():
            txt_files = list(filepath.glob('*.txt'))
            if txt_files:
                summary_file = txt_files[0]
        filepath = summary_file

    logger.info(f"Parsing DRUMMER output: {filepath}")

    if not filepath.exists():
        logger.warning(f"File not found: {filepath}")
        return _empty_dataframe()

    try:
        df = pd.read_csv(filepath, sep='\t', comment='#')

        if df.empty:
            logger.warning(f"Empty DRUMMER file: {filepath}")
            return _empty_dataframe()

        pos_col = _find_column(df, ["transcript_pos", "position", "pos", "genomic_pos"])
        ref_col = _find_column(df, ["transcript_id", "chrom", "chr", "reference"])
        score_col = _find_column(df, ["max.G_padj", "G_padj", "OR_padj", "pval", "pvalue", "p_value"])

        if pos_col is None:
            return _empty_dataframe()

        out = pd.DataFrame(
            {
                "reference": df[ref_col] if ref_col else _extract_reference_from_path(filepath),
                "position": _safe_numeric(df[pos_col]),
                "score": _safe_numeric(df[score_col]) if score_col else np.nan,
                "score_type": "pvalue",
                "pvalue": _safe_numeric(df[score_col]) if score_col else np.nan,
            }
        )

        for extra in ["max.G_padj", "OR_padj", "G_padj", "odds_ratio"]:
            col = _find_column(df, [extra])
            if col:
                out[extra] = _safe_numeric(df[col])

        return _finalize_result(out, "drummer", filepath, sample_id, replicate, None)
    except Exception as e:
        logger.error(f"Error parsing DRUMMER file {filepath}: {e}")
        return _empty_dataframe()


def parse_jacusa2(
    filepath: Union[str, Path],
    sample_id: str = None,
    replicate: str = None
) -> pd.DataFrame:
    filepath = Path(filepath)
    logger.info("Parsing JACUSA2 output: %s", filepath)

    if not filepath.exists():
        return _empty_dataframe()

    try:
        data_lines: List[str] = []
        with open(filepath, "r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("##"):
                    continue
                if line.startswith("#"):
                    continue
                if line.strip():
                    data_lines.append(line.rstrip("\n"))

        if not data_lines:
            return _empty_dataframe()

        records = [ln.split("\t") for ln in data_lines]
        n_cols = len(records[0])

        if n_cols >= 11:
            cols = [
                "contig",
                "start",
                "end",
                "name",
                "score",
                "strand",
                "bases11",
                "bases21",
                "info",
                "filter",
                "ref",
            ]
        else:
            cols = [
                "contig",
                "start",
                "end",
                "name",
                "score",
                "strand",
                "info",
                "filter",
                "ref",
            ]

        df = pd.DataFrame(records, columns=cols[:n_cols])

        out = pd.DataFrame(
            {
                "reference": df["contig"].astype(str),
                "position": _safe_numeric(df["start"]) + 1,
                "score": _safe_numeric(df["score"]),
                "score_type": "score",
                "pvalue": np.nan,
            }
        )

        return _finalize_result(out, "jacusa2", filepath, sample_id, replicate, None)
    except Exception as exc:
        logger.error("Error parsing JACUSA2 file %s: %s", filepath, exc)
        return _empty_dataframe()


def parse_nanodoc(
    filepath: Union[str, Path], sample_id: str = None, replicate: str = None
) -> pd.DataFrame:
    filepath = Path(filepath)
    logger.info("Parsing nanoDoc output: %s", filepath)

    if not filepath.exists():
        return _empty_dataframe()

    try:
        col_names = [
            "pos", "kmer", "depth_tgt", "depth_ref",
            "med_current", "mad_current", "med_currentR", "mad_currentR",
            "current_ratio", "scoreSide1", "scoreSide2", "scoreTotal",
        ]
        df = pd.read_csv(filepath, sep="\t", header=None, names=col_names)
        if df.empty:
            return _empty_dataframe()

        # Extract replicate from nanodoc filename pattern: {cov}x_ndoc_{ref}_{rep}.txt
        if replicate is None:
            m = re.search(r"_(\d+)$", filepath.stem)
            if m:
                replicate = f"rep{m.group(1)}"

        out = pd.DataFrame(
            {
                "reference": _extract_reference_from_path(filepath),
                "position": _safe_numeric(df["pos"]),
                "score": _safe_numeric(df["scoreTotal"]),
                "score_type": "score",
                "pvalue": np.nan,
            }
        )

        return _finalize_result(out, "nanodoc", filepath, sample_id, replicate, None)
    except Exception as exc:
        logger.error("Error parsing nanoDoc file %s: %s", filepath, exc)
        return _empty_dataframe()


def parse_tool_output(
    tool: str,
    filepath: Union[str, Path],
    sample_id: str = None,
    replicate: str = None,
    coverage: str = None,
    differr_score_field: str = "g_fdr_neglog10",
    **kwargs
) -> pd.DataFrame:
    """
    Parse output from any supported tool.

    Args:
        tool: Tool name (tombo, yanocomp, nanocompore, xpore, eligos,
              epinano, differr, drummer, jacusa2)
        filepath: Path to output file or directory
        sample_id: Sample identifier
        replicate: Replicate identifier
        coverage: Coverage level (e.g., "100x"). If None, extracted from path.
        **kwargs: Additional tool-specific arguments

    Returns:
        Standardized DataFrame
    """
    parsers = {
        "tombo": parse_tombo,
        "yanocomp": parse_yanocomp,
        "nanocompore": parse_nanocompore,
        "xpore": parse_xpore,
        "eligos": parse_eligos,
        "epinano": parse_epinano,
        "differr": parse_differr,
        "drummer": parse_drummer,
        "jacusa2": parse_jacusa2,
        "nanodoc": parse_nanodoc,
    }

    tool_lower = tool.lower()
    if tool_lower not in parsers:
        raise ValueError(f"Unknown tool: {tool}. Supported tools: {list(parsers.keys())}")

    parser_kwargs = dict(kwargs)
    if tool_lower == "differr":
        parser_kwargs["score_field"] = differr_score_field

    df = parsers[tool_lower](filepath, sample_id=sample_id, replicate=replicate, **parser_kwargs)

    # Add coverage column
    if not df.empty:
        if coverage is None:
            coverage = _extract_coverage_from_path(Path(filepath))
        df['coverage'] = coverage

        # Add _imputed column (False for all natively parsed positions)
        df['_imputed'] = False

        # Ensure score column uses NaN for invalid/missing values
        df['score'] = pd.to_numeric(df['score'], errors='coerce')
        if 'pvalue' in df.columns:
            df['pvalue'] = pd.to_numeric(df['pvalue'], errors='coerce')

    return df


# -----------------------------------------------------------------------------
# Loader
# -----------------------------------------------------------------------------

def _discover_tool_inputs(tool: str, tool_dir: Path) -> List[Path]:
    """Return canonical input paths for each tool."""
    inputs: List[Path] = []

    if not tool_dir.exists():
        return inputs

    if tool == "tombo":
        # Support both flat (tool/*.csv) and nested (tool/<ref>/*.csv) layouts.
        inputs.extend(sorted(tool_dir.glob("*.csv")))
        inputs.extend(sorted(tool_dir.glob("*/*.csv")))
    elif tool == "yanocomp":
        inputs.extend(sorted(tool_dir.glob("*/*.bed")))
    elif tool == "nanocompore":
        inputs.extend(sorted(tool_dir.glob("*_sampcomp")))
    elif tool == "xpore":
        # Prefer explicit diffmod.table files
        tables = sorted(tool_dir.glob("*_diffmod/diffmod.table"))
        if tables:
            inputs.extend(tables)
        else:
            inputs.extend(sorted(tool_dir.glob("*_diffmod")))
    elif tool == "eligos":
        inputs.extend(sorted(tool_dir.glob("*/*_eligos")))
    elif tool == "epinano":
        inputs.extend(sorted(tool_dir.glob("*/*_epinano*")))
    elif tool == "differr":
        inputs.extend(sorted(tool_dir.glob("*/*.bed")))
    elif tool == "drummer":
        inputs.extend(sorted(tool_dir.glob("*/*_drummer")))
    elif tool == "jacusa2":
        inputs.extend(sorted(tool_dir.glob("*/*.bed")))
    elif tool == "nanodoc":
        inputs.extend(sorted(tool_dir.glob("*.txt")))

    # Generic fallback
    if not inputs:
        for item in sorted(tool_dir.iterdir()):
            if item.is_file() and item.suffix in {".csv", ".tsv", ".txt", ".bed", ".table"}:
                inputs.append(item)
            elif item.is_dir() and tool in DIR_CAPABLE_TOOLS:
                inputs.append(item)

    return inputs


def load_all_tool_outputs(
    output_dir: Union[str, Path],
    tools: Optional[List[str]] = None,
    coverage_dirs: bool = False,
    differr_score_field: str = "g_fdr_neglog10",
) -> Dict[str, pd.DataFrame]:
    """
    Load outputs from all tools in a standard output directory structure.

    Expected directory structure (standard):
    output_dir/
    ├── tombo/
    ├── yanocomp/{reference}/
    ├── nanocompore/
    ├── xpore/
    ├── eligos/{reference}/
    ├── epinano/{reference}/
    ├── differr/{reference}/
    ├── drummer/{reference}/
    └── jacusa2/{reference}/

    Or with coverage directories (coverage_dirs=True):
    output_dir/
    ├── 100x/
    │   ├── tombo/
    │   ├── yanocomp/
    │   └── ...
    ├── 500x/
    │   ├── tombo/
    │   └── ...
    └── 1000x/
        └── ...

    Args:
        output_dir: Path to modifications output directory
        tools: List of tools to load (default: all)
        coverage_dirs: If True, expect {coverage}/tool/ structure

    Returns:
        Dictionary mapping tool names to DataFrames (with coverage column)
    """
    output_dir = Path(output_dir)

    if tools is None:
        tools = ['tombo', 'yanocomp', 'nanocompore', 'xpore', 'eligos',
                 'epinano', 'differr', 'drummer', 'jacusa2']

    if coverage_dirs:
        return _load_with_coverage_dirs(output_dir, tools, differr_score_field=differr_score_field)
    return _load_standard_structure(output_dir, tools, differr_score_field=differr_score_field)


def _load_standard_structure(
    output_dir: Path,
    tools: List[str],
    differr_score_field: str = "g_fdr_neglog10",
) -> Dict[str, pd.DataFrame]:
    results: Dict[str, pd.DataFrame] = {}

    for tool in tools:
        tool_dir = output_dir / tool
        if not tool_dir.exists():
            logger.warning("Tool directory not found: %s", tool_dir)
            results[tool] = _empty_dataframe()
            continue

        inputs = _discover_tool_inputs(tool, tool_dir)
        dfs: List[pd.DataFrame] = []

        for inp in inputs:
            try:
                df = parse_tool_output(tool, inp, differr_score_field=differr_score_field)
                if not df.empty:
                    dfs.append(df)
            except Exception as exc:
                logger.warning("Failed to parse %s for %s: %s", inp, tool, exc)

        results[tool] = pd.concat(dfs, ignore_index=True) if dfs else _empty_dataframe()
        logger.info("Loaded %d rows for %s", len(results[tool]), tool)

    return results


def _load_with_coverage_dirs(
    output_dir: Path,
    tools: List[str],
    differr_score_field: str = "g_fdr_neglog10",
) -> Dict[str, pd.DataFrame]:
    """
    Load tool outputs from coverage-structured directories.

    Expected structure: output_dir/{coverage}/tool/output.csv
    """
    results: Dict[str, List[pd.DataFrame]] = {tool: [] for tool in tools}

    # Find coverage directories (match patterns like 100x, 500, coverage_100x)
    coverage_pattern = re.compile(r'^(?:coverage[_-]?)?(\d+)[xX]?$')

    for coverage_dir in sorted(output_dir.iterdir()):
        if not coverage_dir.is_dir():
            continue

        # Check if this is a coverage directory
        match = coverage_pattern.match(coverage_dir.name)
        if not match:
            continue

        coverage = f"{match.group(1)}x"

        for tool in tools:
            tool_dir = coverage_dir / tool
            if not tool_dir.exists():
                continue

            for inp in _discover_tool_inputs(tool, tool_dir):
                try:
                    df = parse_tool_output(
                        tool,
                        inp,
                        coverage=coverage,
                        differr_score_field=differr_score_field,
                    )
                    if not df.empty:
                        results[tool].append(df)
                except Exception as exc:
                    logger.warning("Failed to parse %s (%s): %s", inp, coverage, exc)

    # Concatenate all dataframes for each tool
    final_results: Dict[str, pd.DataFrame] = {}
    for tool in tools:
        final_results[tool] = (
            pd.concat(results[tool], ignore_index=True) if results[tool] else _empty_dataframe()
        )
        logger.info("Loaded %d rows for %s across coverage dirs", len(final_results[tool]), tool)

    return final_results


# =============================================================================
# Helper Functions
# =============================================================================

def _empty_dataframe() -> pd.DataFrame:
    """Return an empty DataFrame with standard columns."""
    return pd.DataFrame(columns=[
        'tool', 'reference', 'position', 'score', 'score_type',
        'pvalue', 'replicate', 'sample_id', 'coverage', '_imputed'
    ])


def _find_column(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    """Find the first matching column from a list of candidates."""
    for col in candidates:
        # Exact match
        if col in df.columns:
            return col
        # Case-insensitive match
        for df_col in df.columns:
            if df_col.lower() == col.lower():
                return df_col
    return None


def _extract_reference_from_filename(filepath: Path) -> str:
    """Extract reference name from filepath."""
    # Common patterns: 16s_rep1_tombo.csv, 23s_sample1.bed
    name = filepath.stem.lower()
    if '16s' in name:
        return '16s'
    elif '23s' in name:
        return '23s'
    elif '5s' in name:
        return '5s'
    # Try parent directory
    parent = filepath.parent.name.lower()
    if '16s' in parent:
        return '16s'
    elif '23s' in parent:
        return '23s'
    elif '5s' in parent:
        return '5s'
    return 'unknown'


def _extract_replicate_from_filename(filepath: Path) -> str:
    """Extract replicate identifier from filepath."""
    import re
    name = filepath.stem
    # Match patterns like rep1, rep_1, replicate1, R1
    match = re.search(r'rep(?:licate)?[_-]?(\d+)|R(\d+)', name, re.IGNORECASE)
    if match:
        num = match.group(1) or match.group(2)
        return f"rep{num}"
    return "rep1"  # Default


def _extract_coverage_from_path(filepath: Path) -> Optional[str]:
    """
    Extract coverage level from filepath.

    Looks for patterns like:
    - {coverage}/tool/output.csv (e.g., 100x/tombo/results.csv)
    - coverage_100x/tool/output.csv
    - 500X/tool/output.csv

    Args:
        filepath: Path to file

    Returns:
        Coverage string (e.g., "100x") or None if not found
    """
    import re

    # Check each parent directory for coverage pattern
    parts = filepath.parts
    for part in parts:
        # Match patterns like: 100x, 500X, coverage_100x, cov100
        match = re.search(r'(?:coverage[_-]?)?(\d+)[xX]', part, re.IGNORECASE)
        if match:
            return f"{match.group(1)}x"
        # Also match plain numbers that could be coverage (e.g., "100", "500")
        if re.match(r'^\d+$', part):
            return f"{part}x"

    return None


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Parse RNA modification detection tool outputs'
    )
    parser.add_argument('tool', help='Tool name')
    parser.add_argument('filepath', help='Path to output file or directory')
    parser.add_argument('--output', '-o', help='Output CSV file')
    parser.add_argument('--sample-id', help='Sample identifier')
    parser.add_argument('--replicate', help='Replicate identifier')

    args = parser.parse_args()

    df = parse_tool_output(
        args.tool,
        args.filepath,
        sample_id=args.sample_id,
        replicate=args.replicate
    )

    if args.output:
        df.to_csv(args.output, index=False)
        print(f"Saved {len(df)} rows to {args.output}")
    else:
        print(df.to_string())
