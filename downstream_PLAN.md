## Downstream Analysis Recovery Plan (Full Correctness, Surgical)

### Summary
Audit of `bin/downstream_analysis/` against real pipeline outputs (`results_100x_allfixed` and `results_100x_allparams_fixcheck`) shows the downstream stack is currently not reliable for benchmarking:
- It is **not wired** into the main workflow (`workflows/rnamodbench.nf` has no downstream call).
- The parser loads only a subset correctly; in real data, several tools are dropped or misparsed.
- Core comparison/metrics logic uses `position` only (not `reference+position`), which can produce incorrect overlap/metrics.
- There are runtime crashes/incompatibilities in the analysis scripts.

This plan fixes correctness first with minimal, targeted changes and no overengineering.

### Findings (Severity-Ordered)

1. **Critical: parser/schema mismatch drops or corrupts multiple tools**
- `parse_yanocomp` assumes 6-column BED and misreads current 18-column output (`bin/downstream_analysis/parse_outputs.py:148`).
- `parse_xpore` never loads `.table` because loader ignores `.table` suffix; it ends up parsing empty `models/` dirs (`bin/downstream_analysis/parse_outputs.py:869`, `bin/downstream_analysis/parse_outputs.py:941`).
- `parse_eligos` expects `test_paired_diff_mod_result.txt`, but current output is `*_combine.txt` (`bin/downstream_analysis/parse_outputs.py:381`).
- `parse_epinano` uses `comment='#'`, which drops `#Ref` header and breaks dataframe construction (`bin/downstream_analysis/parse_outputs.py:471`).
- `parse_differr` hardcodes 7 columns but current output has 14 (`bin/downstream_analysis/parse_outputs.py:543`).
- `parse_drummer` expects legacy columns; current output uses `transcript_pos` / `OR_padj` (`bin/downstream_analysis/parse_outputs.py:626`).
- `parse_jacusa2` assumes 9 fields; current files have 11 (`bin/downstream_analysis/parse_outputs.py:705`).
- `_extract_replicate_from_filename` collapses many directory-based tools to `rep1` (`bin/downstream_analysis/parse_outputs.py:1023`).

2. **Critical: cross-reference collisions in metrics/comparison**
- Multiple modules compare sets of `position` only, ignoring reference identity:
  - `bin/downstream_analysis/benchmark_metrics.py:183`
  - `bin/downstream_analysis/replicate_analysis.py:186`
  - `bin/downstream_analysis/tool_comparison.py:145`
  - `bin/downstream_analysis/data_quality.py:486`
  - `bin/downstream_analysis/position_standardization.py:204`
- This can merge 16S and 23S calls at same numeric positions and invalidate AUROC/AUPRC and concordance.

3. **High: runtime breakages in analysis orchestration**
- `visualization.py` uses package-relative imports that fail when called from `run_analysis.py` script execution path (`bin/downstream_analysis/visualization.py:83`, `bin/downstream_analysis/visualization.py:867`).
- `run_analysis.py` writes a DataFrame directly as text (`bin/downstream_analysis/run_analysis.py:176`, `bin/downstream_analysis/run_analysis.py:178`; source function returns DataFrame at `bin/downstream_analysis/data_quality.py:573`).
- Coverage analysis references non-existent metric fields (`bin/downstream_analysis/coverage_analysis.py:180`).

4. **High: downstream module logic under-collects inputs**
- `DOWNSTREAM_ANALYSIS` takes only the first matched artifact per tool due `break` in each loop (`modules/local/downstream_analysis/main.nf:61` onward).

5. **Medium: downstream stage is not invoked**
- No include/call in `workflows/rnamodbench.nf`, so scripts are effectively inactive in pipeline execution.

---

## Implementation Plan

### Phase 1: Parser Compatibility (No metric semantics change yet)
Files: `bin/downstream_analysis/parse_outputs.py`

1. Update loader rules to parse canonical artifacts only.
- Treat xpore `diffmod.table` as valid input.
- For nanocompore, parse result table only (`outnanocompore_results.tsv` / `outSampComp_results.tsv`), not shift-stats for ranking rows.
- For ELIGOS, support `*_combine.txt`.
- For EpiNano, prefer prediction outputs (`*mismatch*.prediction.csv`, `*sumerr*.prediction.csv`), fallback to per-site CSV.
- Keep backward-compatible fallbacks where possible.

2. Update tool-specific parsers to current real formats.
- Yanocomp: parse extended BED, extract reference from k-mer name field, derive central position.
- xPore: support directory and file paths robustly.
- EpiNano: stop stripping header with comment mode; parse `chr_pos` into `(reference, position)` when using prediction files.
- DiffErr: support 14-column format and robust numeric coercion.
- DRUMMER: map modern columns (`transcript_pos`, adjusted p-value fields).
- JACUSA2: parse 11-column data layout after metadata lines.

3. Fix replicate/sample extraction for directory-based outputs.
- Add path-aware replicate extraction from parent directory names (not only filename stem).

### Phase 2: Full Correctness for Comparisons/Metrics (Reference-aware keys)
Files:
- `bin/downstream_analysis/benchmark_metrics.py`
- `bin/downstream_analysis/replicate_analysis.py`
- `bin/downstream_analysis/tool_comparison.py`
- `bin/downstream_analysis/data_quality.py`
- `bin/downstream_analysis/position_standardization.py`

1. Replace position-only set logic with `(reference, position)` keys.
- All overlap, Jaccard, consensus, and label prep paths use reference-aware tuples.
- Keep output tables readable by writing both `reference` and `position` columns explicitly.

2. Fix standardization to avoid reference collisions.
- Replace `position_ref_map[position]` with tuple-safe mapping.
- Union/intersection across replicates operates on `(reference, position)`.

3. Keep behavior surgical.
- No new statistical methods.
- No change to caller thresholds/tool execution.
- Only correctness of downstream bookkeeping and matching.

### Phase 3: Runtime Stability Fixes
Files:
- `bin/downstream_analysis/run_analysis.py`
- `bin/downstream_analysis/visualization.py`
- `bin/downstream_analysis/coverage_analysis.py`

1. Make imports script-safe and package-safe.
- Use dual import fallback (`from .x` then `from x`) in visualization/comparison paths used by script execution.

2. Fix known orchestration bugs.
- Save quality summary DataFrame as CSV and write a text summary string separately.
- Fix coverage metrics field mapping to `f1_optimal`, `precision_optimal`, `recall_optimal`.

3. Keep failure handling non-fatal where intended.
- Preserve warning-and-continue behavior.

### Phase 4: Downstream Module Wiring (Optional, default off)
Files:
- `modules/local/downstream_analysis/main.nf`
- `workflows/rnamodbench.nf`
- `nextflow.config` / schema/docs as needed

1. Make module consume full outputs, not first match.
- Remove first-hit `break` behavior.
- Prefer passing a single `--input-dir` style directory to analysis rather than one-path-per-tool CLI flags.

2. Wire invocation behind a default-off flag.
- Avoid changing existing successful pipeline behavior unless explicitly enabled.

---

## Public APIs / Interfaces / Types

### Core phases (1–3)
- No pipeline user-facing API changes required.
- Internal downstream dataframe semantics change:
  - set/overlap logic becomes `reference+position` aware.

### Optional phase (4)
- Add optional pipeline param(s), default disabled, e.g.:
  - `--run_downstream` (bool, default `false`)
  - `--ground_truth` (path, optional)

---

## Test Cases and Scenarios

1. **Parser regression on real outputs**
- Input: `results_100x_allfixed/modifications`
- Input: `results_100x_allparams_fixcheck/modifications`
- Expected: non-empty parsed outputs for all active tools with artifacts present (especially xpore, eligos, epinano, differr, jacusa2).

2. **Replicate extraction correctness**
- For directory-based tools, parsed replicates should reflect `rep1/rep2/rep3`, not all `rep1`.

3. **Reference-aware collision test**
- Synthetic dataset with same numeric position on two references.
- Expected: no cross-reference merging in metrics/comparison/consensus.

4. **Run-analysis smoke tests**
- Script mode without ground truth: completes, writes quality/comparison artifacts.
- Script mode with synthetic ground truth: completes PR/ROC/metrics path without import/type errors.

5. **Coverage analysis unit test**
- Confirm no attribute errors and non-NaN metric fields for valid input.

6. **(Optional phase 4) workflow integration test**
- With downstream disabled: no change to existing pipeline outputs/behavior.
- With downstream enabled: downstream output directory produced.

---

## Assumptions and Defaults

1. Scope is limited to downstream analysis correctness and activation readiness.
2. No changes to upstream caller execution logic or thresholds.
3. No refactor to a new framework/package layout; keep current file structure.
4. Optional workflow wiring remains default-off to avoid collateral regressions.
5. Existing successful main pipeline behavior is treated as baseline and must remain intact.
