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
from pathlib import Path
from typing import Optional, Dict, List, Union

import pandas as pd
import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Tools whose parsers support directory inputs (have is_dir() handling)
DIR_CAPABLE_TOOLS = {'nanocompore', 'xpore', 'eligos', 'epinano', 'drummer'}


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
            'score_type': 'statistic',
            'pvalue': df[stat_col].astype(float) if stat_col and 'p' in stat_col.lower() else np.nan,
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

        logger.info(f"Parsed {len(result)} positions from Yanocomp")
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
    """
    Parse Nanocompore sampcomp results TSV.

    Nanocompore outputs outSampComp_results.tsv with columns:
    - ref_id, ref_pos (or pos), ref_kmer
    - Various p-value columns: GMM_logit_pvalue, KS_pvalue, MW_pvalue, etc.
    - Shift statistics

    Args:
        filepath: Path to Nanocompore results TSV
        sample_id: Sample identifier
        replicate: Replicate identifier
        pvalue_col: Which p-value column to use as primary score

    Returns:
        Standardized DataFrame
    """
    filepath = Path(filepath)

    # Handle directory input (find results file)
    if filepath.is_dir():
        results_file = filepath / 'outSampComp_results.tsv'
        if not results_file.exists():
            # Try alternative names
            for name in ['outSampComp_results.tsv', 'sampcomp_results.tsv', 'results.tsv']:
                candidate = filepath / name
                if candidate.exists():
                    results_file = candidate
                    break
        filepath = results_file

    logger.info(f"Parsing Nanocompore output: {filepath}")

    if not filepath.exists():
        logger.warning(f"File not found: {filepath}")
        return _empty_dataframe()

    try:
        df = pd.read_csv(filepath, sep='\t')

        if df.empty:
            logger.warning(f"Empty Nanocompore file: {filepath}")
            return _empty_dataframe()

        # Find position column
        pos_col = _find_column(df, ['ref_pos', 'pos', 'position'])
        ref_col = _find_column(df, ['ref_id', 'ref', 'reference', 'chr', 'chrom'])

        if pos_col is None:
            logger.error(f"Could not find position column in Nanocompore: {df.columns.tolist()}")
            return _empty_dataframe()

        # Find best p-value column
        pval_col = _find_column(df, [pvalue_col, 'GMM_logit_pvalue', 'KS_pvalue', 'pvalue'])

        result = pd.DataFrame({
            'tool': 'nanocompore',
            'reference': df[ref_col] if ref_col else _extract_reference_from_filename(filepath),
            'position': df[pos_col].astype(int),
            'score': df[pval_col].astype(float) if pval_col else np.nan,
            'score_type': 'pvalue',
            'pvalue': df[pval_col].astype(float) if pval_col else np.nan,
            'replicate': replicate or _extract_replicate_from_filename(filepath),
            'sample_id': sample_id or filepath.stem
        })

        logger.info(f"Parsed {len(result)} positions from Nanocompore")
        return result

    except Exception as e:
        logger.error(f"Error parsing Nanocompore file {filepath}: {e}")
        return _empty_dataframe()


def parse_xpore(
    filepath: Union[str, Path],
    sample_id: str = None,
    replicate: str = None
) -> pd.DataFrame:
    """
    Parse xPore diffmod.table output.

    xPore outputs diffmod.table with columns:
    - id, position, kmer
    - diff_mod_rate_*: differential modification rate
    - pval_*: p-values
    - z_score_*: z-scores

    Args:
        filepath: Path to xPore diffmod.table or directory
        sample_id: Sample identifier
        replicate: Replicate identifier

    Returns:
        Standardized DataFrame
    """
    filepath = Path(filepath)

    # Handle directory input
    if filepath.is_dir():
        diffmod_file = filepath / 'diffmod.table'
        if not diffmod_file.exists():
            # Try to find any .table file
            table_files = list(filepath.glob('*.table'))
            if table_files:
                diffmod_file = table_files[0]
        filepath = diffmod_file

    logger.info(f"Parsing xPore output: {filepath}")

    if not filepath.exists():
        logger.warning(f"File not found: {filepath}")
        return _empty_dataframe()

    try:
        df = pd.read_csv(filepath, sep='\t')

        if df.empty:
            logger.warning(f"Empty xPore file: {filepath}")
            return _empty_dataframe()

        # Find relevant columns
        pos_col = _find_column(df, ['position', 'pos'])
        ref_col = _find_column(df, ['id', 'ref_id', 'reference'])

        # Find p-value column (xPore uses pval_* naming)
        pval_cols = [c for c in df.columns if c.startswith('pval_')]
        pval_col = pval_cols[0] if pval_cols else _find_column(df, ['pvalue', 'p_value'])

        # Find diff_mod_rate column
        rate_cols = [c for c in df.columns if 'diff_mod_rate' in c]
        rate_col = rate_cols[0] if rate_cols else None

        result = pd.DataFrame({
            'tool': 'xpore',
            'reference': df[ref_col] if ref_col else _extract_reference_from_filename(filepath),
            'position': df[pos_col].astype(int) if pos_col else np.nan,
            'score': df[pval_col].astype(float) if pval_col else np.nan,
            'score_type': 'pvalue',
            'pvalue': df[pval_col].astype(float) if pval_col else np.nan,
            'replicate': replicate or _extract_replicate_from_filename(filepath),
            'sample_id': sample_id or filepath.stem
        })

        # Add diff_mod_rate as additional column if available
        if rate_col:
            result['diff_mod_rate'] = df[rate_col].astype(float)

        logger.info(f"Parsed {len(result)} positions from xPore")
        return result

    except Exception as e:
        logger.error(f"Error parsing xPore file {filepath}: {e}")
        return _empty_dataframe()


def parse_eligos(
    filepath: Union[str, Path],
    sample_id: str = None,
    replicate: str = None
) -> pd.DataFrame:
    """
    Parse ELIGOS pair_diff_mod output.

    ELIGOS outputs test_paired_diff_mod_result.txt with columns:
    - chrom, start_loc, end_loc, strand
    - pval, adj_pval, oddR (odds ratio)
    - ESB (error-specific bias)

    Args:
        filepath: Path to ELIGOS result file or directory
        sample_id: Sample identifier
        replicate: Replicate identifier

    Returns:
        Standardized DataFrame
    """
    filepath = Path(filepath)

    # Handle directory input
    if filepath.is_dir():
        result_file = filepath / 'test_paired_diff_mod_result.txt'
        if not result_file.exists():
            txt_files = list(filepath.glob('*result*.txt'))
            if txt_files:
                result_file = txt_files[0]
        filepath = result_file

    logger.info(f"Parsing ELIGOS output: {filepath}")

    if not filepath.exists():
        logger.warning(f"File not found: {filepath}")
        return _empty_dataframe()

    try:
        df = pd.read_csv(filepath, sep='\t')

        if df.empty:
            logger.warning(f"Empty ELIGOS file: {filepath}")
            return _empty_dataframe()

        # Find relevant columns
        pos_col = _find_column(df, ['start_loc', 'start', 'pos', 'position'])
        ref_col = _find_column(df, ['chrom', 'chr', 'reference'])
        pval_col = _find_column(df, ['pval', 'pvalue', 'p_value', 'adj_pval'])
        odds_col = _find_column(df, ['oddR', 'odds_ratio', 'oddsratio'])

        result = pd.DataFrame({
            'tool': 'eligos',
            'reference': df[ref_col] if ref_col else _extract_reference_from_filename(filepath),
            'position': (df[pos_col] + 1).astype(int) if pos_col else np.nan,  # Convert to 1-based
            'score': df[pval_col].astype(float) if pval_col else np.nan,
            'score_type': 'pvalue',
            'pvalue': df[pval_col].astype(float) if pval_col else np.nan,
            'replicate': replicate or _extract_replicate_from_filename(filepath),
            'sample_id': sample_id or filepath.stem
        })

        # Add odds ratio if available
        if odds_col:
            result['odds_ratio'] = df[odds_col].astype(float)

        logger.info(f"Parsed {len(result)} positions from ELIGOS")
        return result

    except Exception as e:
        logger.error(f"Error parsing ELIGOS file {filepath}: {e}")
        return _empty_dataframe()


def parse_epinano(
    filepath: Union[str, Path],
    sample_id: str = None,
    replicate: str = None
) -> pd.DataFrame:
    """
    Parse EpiNano error output.

    EpiNano outputs include:
    - native_per_site.csv / ivt_per_site.csv: raw per-site error data
    - *_diff_err.csv: differential error analysis with z-scores

    Args:
        filepath: Path to EpiNano output file or directory
        sample_id: Sample identifier
        replicate: Replicate identifier

    Returns:
        Standardized DataFrame
    """
    filepath = Path(filepath)

    # Handle directory input - look for diff_err file first
    if filepath.is_dir():
        diff_files = list(filepath.glob('*diff_err*.csv'))
        if diff_files:
            filepath = diff_files[0]
        else:
            # Fall back to native per_site file
            per_site_files = list(filepath.glob('native_per_site.csv'))
            if per_site_files:
                filepath = per_site_files[0]

    logger.info(f"Parsing EpiNano output: {filepath}")

    if not filepath.exists():
        logger.warning(f"File not found: {filepath}")
        return _empty_dataframe()

    try:
        # Read file, handling comment lines
        df = pd.read_csv(filepath, comment='#')

        if df.empty:
            logger.warning(f"Empty EpiNano file: {filepath}")
            return _empty_dataframe()

        # Find relevant columns
        pos_col = _find_column(df, ['pos', 'position', 'start'])
        ref_col = _find_column(df, ['X.Ref', '#Ref', 'Ref', 'chr', 'chrom', 'reference'])

        # Look for z-score or error columns
        zscore_col = _find_column(df, ['z_score', 'zscore', 'z-score', 'SumErr_z'])
        error_col = _find_column(df, ['sum_err', 'sumerr', 'diff_err', 'mis', 'error'])

        # Use z-score if available, otherwise error difference
        if zscore_col:
            score_values = df[zscore_col].astype(float)
            score_type = 'zscore'
        elif error_col:
            score_values = df[error_col].astype(float)
            score_type = 'error'
        else:
            score_values = np.nan
            score_type = 'unknown'

        result = pd.DataFrame({
            'tool': 'epinano',
            'reference': df[ref_col] if ref_col else _extract_reference_from_filename(filepath),
            'position': df[pos_col].astype(int) if pos_col else np.nan,
            'score': score_values,
            'score_type': score_type,
            'pvalue': np.nan,  # EpiNano uses z-scores, not p-values
            'replicate': replicate or _extract_replicate_from_filename(filepath),
            'sample_id': sample_id or filepath.stem
        })

        logger.info(f"Parsed {len(result)} positions from EpiNano")
        return result

    except Exception as e:
        logger.error(f"Error parsing EpiNano file {filepath}: {e}")
        return _empty_dataframe()


def parse_differr(
    filepath: Union[str, Path],
    sample_id: str = None,
    replicate: str = None
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
        # Read BED file with expected columns
        df = pd.read_csv(
            filepath,
            sep='\t',
            comment='#',
            names=['chr', 'start', 'end', 'odds_ratio', 'g_stat', 'pval', 'fdr'],
            header=None
        )

        # Check if first row is header and skip if needed
        if df.iloc[0]['chr'] == 'chr':
            df = df.iloc[1:].reset_index(drop=True)

        if df.empty:
            logger.warning(f"Empty DiffErr file: {filepath}")
            return _empty_dataframe()

        result = pd.DataFrame({
            'tool': 'differr',
            'reference': df['chr'],
            'position': (df['start'].astype(float) + 1).astype(int),  # Convert to 1-based
            'score': pd.to_numeric(df['fdr'], errors='coerce'),
            'score_type': 'fdr',
            'pvalue': pd.to_numeric(df['pval'], errors='coerce'),
            'replicate': replicate or _extract_replicate_from_filename(filepath),
            'sample_id': sample_id or filepath.stem
        })

        # Add odds ratio
        result['odds_ratio'] = pd.to_numeric(df['odds_ratio'], errors='coerce')

        logger.info(f"Parsed {len(result)} positions from DiffErr")
        return result

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

        # Find relevant columns
        pos_col = _find_column(df, ['position', 'pos', 'genomic_pos'])
        ref_col = _find_column(df, ['chrom', 'chr', 'transcript_id', 'reference'])
        pval_col = _find_column(df, ['pval', 'pvalue', 'p_value'])
        odds_col = _find_column(df, ['odds_ratio', 'oddsratio', 'oddR'])

        result = pd.DataFrame({
            'tool': 'drummer',
            'reference': df[ref_col] if ref_col else _extract_reference_from_filename(filepath),
            'position': df[pos_col].astype(int) if pos_col else np.nan,
            'score': df[pval_col].astype(float) if pval_col else np.nan,
            'score_type': 'pvalue',
            'pvalue': df[pval_col].astype(float) if pval_col else np.nan,
            'replicate': replicate or _extract_replicate_from_filename(filepath),
            'sample_id': sample_id or filepath.stem
        })

        # Add odds ratio if available
        if odds_col:
            result['odds_ratio'] = df[odds_col].astype(float)

        logger.info(f"Parsed {len(result)} positions from DRUMMER")
        return result

    except Exception as e:
        logger.error(f"Error parsing DRUMMER file {filepath}: {e}")
        return _empty_dataframe()


def parse_jacusa2(
    filepath: Union[str, Path],
    sample_id: str = None,
    replicate: str = None
) -> pd.DataFrame:
    """
    Parse JACUSA2 BED output.

    JACUSA2 outputs BED format with columns:
    #contig, start, end, name, score, strand, info, filter, ref

    Args:
        filepath: Path to JACUSA2 BED file
        sample_id: Sample identifier
        replicate: Replicate identifier

    Returns:
        Standardized DataFrame
    """
    filepath = Path(filepath)
    logger.info(f"Parsing JACUSA2 output: {filepath}")

    if not filepath.exists():
        logger.warning(f"File not found: {filepath}")
        return _empty_dataframe()

    try:
        # Read file, skipping JACUSA2 header lines
        with open(filepath, 'r') as f:
            lines = f.readlines()

        # Find data lines (skip ## comments and header)
        data_lines = []
        header_idx = None
        for i, line in enumerate(lines):
            if line.startswith('##'):
                continue
            if line.startswith('#contig') or line.startswith('#'):
                header_idx = i
                continue
            data_lines.append(line.strip())

        if not data_lines:
            logger.warning(f"No data in JACUSA2 file: {filepath}")
            return _empty_dataframe()

        # Parse as TSV
        from io import StringIO
        df = pd.read_csv(
            StringIO('\n'.join(data_lines)),
            sep='\t',
            names=['contig', 'start', 'end', 'name', 'score', 'strand', 'info', 'filter', 'ref'],
            header=None
        )

        if df.empty:
            logger.warning(f"Empty JACUSA2 file: {filepath}")
            return _empty_dataframe()

        result = pd.DataFrame({
            'tool': 'jacusa2',
            'reference': df['contig'],
            'position': (df['start'].astype(int) + 1),  # Convert to 1-based
            'score': pd.to_numeric(df['score'], errors='coerce'),
            'score_type': 'score',
            'pvalue': np.nan,  # JACUSA2 uses score, not p-value
            'replicate': replicate or _extract_replicate_from_filename(filepath),
            'sample_id': sample_id or filepath.stem
        })

        logger.info(f"Parsed {len(result)} positions from JACUSA2")
        return result

    except Exception as e:
        logger.error(f"Error parsing JACUSA2 file {filepath}: {e}")
        return _empty_dataframe()


def parse_tool_output(
    tool: str,
    filepath: Union[str, Path],
    sample_id: str = None,
    replicate: str = None,
    coverage: str = None,
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
        'tombo': parse_tombo,
        'yanocomp': parse_yanocomp,
        'nanocompore': parse_nanocompore,
        'xpore': parse_xpore,
        'eligos': parse_eligos,
        'epinano': parse_epinano,
        'differr': parse_differr,
        'drummer': parse_drummer,
        'jacusa2': parse_jacusa2,
    }

    tool_lower = tool.lower()
    if tool_lower not in parsers:
        raise ValueError(f"Unknown tool: {tool}. Supported tools: {list(parsers.keys())}")

    df = parsers[tool_lower](filepath, sample_id=sample_id, replicate=replicate, **kwargs)

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


def load_all_tool_outputs(
    output_dir: Union[str, Path],
    tools: Optional[List[str]] = None,
    coverage_dirs: bool = False
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

    results = {}

    if coverage_dirs:
        # Handle {coverage}/tool/output.csv structure
        results = _load_with_coverage_dirs(output_dir, tools)
    else:
        # Standard structure: tool/{reference}/output.csv
        results = _load_standard_structure(output_dir, tools)

    return results


def _load_standard_structure(
    output_dir: Path,
    tools: List[str]
) -> Dict[str, pd.DataFrame]:
    """Load tool outputs from standard directory structure."""
    results = {}

    for tool in tools:
        tool_dir = output_dir / tool
        if not tool_dir.exists():
            logger.warning(f"Tool directory not found: {tool_dir}")
            continue

        # Collect all output files for this tool
        tool_dfs = []

        # Handle nested directory structure (tool/{reference}/)
        for subdir in tool_dir.iterdir():
            if subdir.is_dir():
                # Look for output files in subdirectory
                for filepath in subdir.glob('*'):
                    if filepath.is_file() and filepath.suffix in ['.csv', '.tsv', '.txt', '.bed']:
                        try:
                            df = parse_tool_output(tool, filepath)
                            if not df.empty:
                                tool_dfs.append(df)
                        except Exception as e:
                            logger.warning(f"Failed to parse {filepath}: {e}")
                    elif filepath.is_dir() and tool in DIR_CAPABLE_TOOLS:
                        # Directory outputs (e.g. ELIGOS, EPINANO, DRUMMER)
                        try:
                            df = parse_tool_output(tool, filepath)
                            if not df.empty:
                                tool_dfs.append(df)
                        except Exception as e:
                            logger.warning(f"Failed to parse directory {filepath}: {e}")
            elif subdir.is_file() and subdir.suffix in ['.csv', '.tsv', '.txt', '.bed']:
                # Direct file in tool directory
                try:
                    df = parse_tool_output(tool, subdir)
                    if not df.empty:
                        tool_dfs.append(df)
                except Exception as e:
                    logger.warning(f"Failed to parse {subdir}: {e}")

        if tool_dfs:
            results[tool] = pd.concat(tool_dfs, ignore_index=True)
            logger.info(f"Loaded {len(results[tool])} total positions for {tool}")
        else:
            results[tool] = _empty_dataframe()
            logger.warning(f"No data loaded for {tool}")

    return results


def _load_with_coverage_dirs(
    output_dir: Path,
    tools: List[str]
) -> Dict[str, pd.DataFrame]:
    """
    Load tool outputs from coverage-structured directories.

    Expected structure: output_dir/{coverage}/tool/output.csv
    """
    import re

    results = {tool: [] for tool in tools}

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
        logger.info(f"Processing coverage level: {coverage}")

        for tool in tools:
            tool_dir = coverage_dir / tool
            if not tool_dir.exists():
                continue

            # Load files from tool directory
            for item in tool_dir.iterdir():
                if item.is_dir():
                    # Nested structure (e.g., tool/{reference}/)
                    for filepath in item.glob('*'):
                        if filepath.is_file() and filepath.suffix in ['.csv', '.tsv', '.txt', '.bed']:
                            try:
                                df = parse_tool_output(tool, filepath, coverage=coverage)
                                if not df.empty:
                                    results[tool].append(df)
                            except Exception as e:
                                logger.warning(f"Failed to parse {filepath}: {e}")
                        elif filepath.is_dir() and tool in DIR_CAPABLE_TOOLS:
                            try:
                                df = parse_tool_output(tool, filepath, coverage=coverage)
                                if not df.empty:
                                    results[tool].append(df)
                            except Exception as e:
                                logger.warning(f"Failed to parse directory {filepath}: {e}")
                elif item.is_file() and item.suffix in ['.csv', '.tsv', '.txt', '.bed']:
                    try:
                        df = parse_tool_output(tool, item, coverage=coverage)
                        if not df.empty:
                            results[tool].append(df)
                    except Exception as e:
                        logger.warning(f"Failed to parse {item}: {e}")

    # Concatenate all dataframes for each tool
    final_results = {}
    for tool in tools:
        if results[tool]:
            final_results[tool] = pd.concat(results[tool], ignore_index=True)
            n_coverages = final_results[tool]['coverage'].nunique()
            logger.info(f"Loaded {len(final_results[tool])} positions for {tool} across {n_coverages} coverage levels")
        else:
            final_results[tool] = _empty_dataframe()
            logger.warning(f"No data loaded for {tool}")

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
