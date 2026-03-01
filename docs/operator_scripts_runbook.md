# Operator Scripts Runbook

## 1) Purpose And Audience
This runbook is the canonical operator guide for running and maintaining the core RNAModBench execution scripts:

1. `generate_samplesheets.sh`
2. `run_all_coverages.sh`
3. `bin/downstream_analysis/collate_runs.py`
4. `bin/downstream_analysis/run_analysis.py`

Use this document when you need reproducible command history, clear success/failure criteria, and operational recovery steps.

## 2) Prerequisites And Environment Assumptions

### Required working directory
All commands are written assuming you run from repo root:

```bash
cd /home/bmorampa/rnamodbench
```

### Pipeline runtime prerequisites
1. `nextflow` installed and on `PATH`
2. Container runtime/profile support for your chosen profile (`singularity` in current batch script)
3. Input mapping file present: `references.csv`
4. Ground truth file present: `ground_truth_mod_positions.csv`

### Data path assumptions used by helper scripts
`generate_samplesheets.sh` currently expects FASTQ and FAST5 at:
1. `/home/bmorampa/covbench/results_all_tools_16S/filtlong`
2. `/home/bmorampa/covbench/results_all_tools_23S/filtlong`
3. `/home/bmorampa/k12_native_fast5`
4. `/home/bmorampa/k12_ivt_fast5`

If your data lives elsewhere, update `generate_samplesheets.sh` before running.

### Python environment for downstream scripts
`bin/downstream_analysis/run_analysis.py` requires Python packages including `scikit-learn`.

Local environment file:
`modules/local/downstream_analysis/environment.yml`

`bin/downstream_analysis/collate_runs.py` only requires common Python data stack (`pandas`) and can be run standalone.

## 3) Repository Conventions And Path Map

### Core directories
1. `samplesheets/`: generated per-coverage samplesheets (`samplesheet_5x.csv`, etc.)
2. `covbench_results/`: one output directory per coverage (`results_5x`, `results_10x`, ...)
3. `covbench_logs/`: one Nextflow batch log per coverage (`run_5x.log`, ...)
4. `covbench_results/results_*/downstream_analysis/collation/`: per-run collated downstream tables

### Canonical multi-run collation output
Recommended merged output location:
`covbench_results/collated_multi_cov/`

Typical files:
1. `metrics_long.csv`
2. `metrics_summary_long.csv`
3. `window_metrics_long.csv`
4. `window_metrics_summary_long.csv`
5. `lag_metrics_long.csv`
6. `lag_metrics_summary_long.csv`
7. `collate_report.md`
8. `collate_report.json`

## 4) Script Matrix

| Script | When To Use | Required Inputs | Key Outputs | Success Signals |
| --- | --- | --- | --- | --- |
| `generate_samplesheets.sh` | Build all coverage samplesheets and validate FASTQ paths before pipeline run | FASTQ trees + FAST5 roots hardcoded in script | `samplesheets/samplesheet_<COV>.csv` | `All FASTQ paths verified.` and exit `0` |
| `run_all_coverages.sh` | Run Nextflow sequentially across all or selected coverages | `samplesheets/`, `references.csv`, `ground_truth_mod_positions.csv` | `covbench_results/results_<COV>/` + `covbench_logs/run_<COV>.log` | Final batch summary + each desired coverage has downstream `metrics_long.csv` |
| `bin/downstream_analysis/collate_runs.py` | Merge multiple downstream collation outputs into one analysis dataset | Runs via `--runs` and/or discovery via `--runs-root` | Merged CSVs + `collate_report.md/json` | CLI summary with expected run/table counts and report files |
| `bin/downstream_analysis/run_analysis.py` | Re-run downstream analysis for one run (or custom tool inputs) without rerunning full pipeline | Modifications dir or per-tool paths; optional ground truth | `downstream_analysis/` with `collation/`, `metadata/`, plots, metrics | Process completes without Python/runtime error and writes `collation/*.csv` |

## 5) End-To-End Workflow Playbooks

### Playbook A: First full multi-coverage run

```bash
cd /home/bmorampa/rnamodbench

# 1) Build and validate samplesheets
./generate_samplesheets.sh

# 2) Run all configured coverages
./run_all_coverages.sh

# 3) Collate all completed runs
python bin/downstream_analysis/collate_runs.py \
  --runs-root covbench_results \
  --run-glob 'results_*/downstream_analysis' \
  --output covbench_results/collated_multi_cov
```

### Playbook B: Rerun selected coverages only

```bash
cd /home/bmorampa/rnamodbench

# Re-run only selected coverages; completed ones are auto-skipped
./run_all_coverages.sh 50x 100x 500x
```

### Playbook C: Resume after failure

1. Inspect failed coverage log: `covbench_logs/run_<COV>.log`
2. Fix the underlying issue (missing data, environment, disk, etc.)
3. Re-run batch script with same coverage list or all coverages:

```bash
cd /home/bmorampa/rnamodbench
./run_all_coverages.sh
```

The script skips already-completed outputs (where downstream `metrics_long.csv` exists).

### Playbook D: Collate all coverages (discovery mode)

```bash
cd /home/bmorampa/rnamodbench
python bin/downstream_analysis/collate_runs.py \
  --runs-root covbench_results \
  --run-glob 'results_*/downstream_analysis' \
  --output covbench_results/collated_multi_cov
```

### Playbook E: Collate explicit subset and strict QA mode

```bash
cd /home/bmorampa/rnamodbench

# Explicit subset
python bin/downstream_analysis/collate_runs.py \
  --runs \
    covbench_results/results_5x/downstream_analysis \
    covbench_results/results_10x/downstream_analysis \
    covbench_results/results_100x/downstream_analysis \
  --output covbench_results/collated_subset

# Strict mode (non-zero exit on soft issues)
python bin/downstream_analysis/collate_runs.py \
  --runs-root covbench_results \
  --run-glob 'results_*/downstream_analysis' \
  --output covbench_results/collated_multi_cov \
  --strict
```

### Playbook F: Re-run standalone downstream analysis for one run

```bash
cd /home/bmorampa/rnamodbench
python bin/downstream_analysis/run_analysis.py \
  --input-dir covbench_results/results_100x/modifications \
  --ground-truth ground_truth_mod_positions.csv \
  --references-csv references.csv \
  --output-dir covbench_results/results_100x/downstream_analysis_rerun \
  --run-id run_100x_rerun \
  --coverage-label 100x \
  --quality-label base \
  --expected-replicates 3 \
  --min-replicates 2
```

## 6) Per-Script Reference Blocks

### `generate_samplesheets.sh`

### When to use
Before running multi-coverage Nextflow batches, to generate and validate all `samplesheet_<coverage>.csv` files.

### Inputs required
1. FASTQ files in configured 16S and 23S trees
2. FAST5 root dirs for native and IVT

### Run command
```bash
cd /home/bmorampa/rnamodbench
./generate_samplesheets.sh
```

### Expected outputs
1. `samplesheets/samplesheet_5x.csv` ... `samplesheets/samplesheet_1000x.csv`
2. One header + 12 data rows per file (13 lines total)

### How to verify success
```bash
cd /home/bmorampa/rnamodbench
ls samplesheets/samplesheet_*.csv | wc -l
wc -l samplesheets/samplesheet_100x.csv
```
Expected:
1. 25 samplesheet files
2. 13 lines per file

### If it fails, do this
1. Read missing file lines printed by script (`ERROR: missing ...`)
2. Fix path roots in `generate_samplesheets.sh` if your storage layout changed
3. Re-run script and ensure it ends with `All FASTQ paths verified.`

### `run_all_coverages.sh`

### When to use
To run the full RNAModBench pipeline sequentially across all coverages or a selected subset.

### Inputs required
1. `samplesheets/samplesheet_<COV>.csv`
2. `ground_truth_mod_positions.csv`
3. `references.csv`
4. Nextflow and configured profile runtime

### Run command
```bash
cd /home/bmorampa/rnamodbench

# all configured coverages
./run_all_coverages.sh

# subset rerun
./run_all_coverages.sh 50x 100x 500x
```

### Defaults
1. Coverage list: `5x ... 1000x` (defined in script)
2. Results dir: `covbench_results/`
3. Logs dir: `covbench_logs/`
4. Nextflow profile: `singularity`
5. Work dir flag: `-w /mnt/nvme_work`

### Exit behavior
1. Exits non-zero on preflight failure (missing ground truth/samplesheets)
2. Stops batch at first coverage execution failure
3. Completed runs are skipped automatically on rerun

### Expected outputs
1. `covbench_results/results_<COV>/`
2. `covbench_logs/run_<COV>.log`
3. For successful completed run: `downstream_analysis/collation/metrics_long.csv`

### How to verify success
```bash
cd /home/bmorampa/rnamodbench
ls covbench_logs/run_*.log | wc -l
ls covbench_results/results_*/downstream_analysis/collation/metrics_long.csv | wc -l
```

### If it fails, do this
1. Open failed run log: `covbench_logs/run_<COV>.log`
2. Fix root cause
3. Re-run script; completed coverages will be skipped

### `bin/downstream_analysis/collate_runs.py`

### When to use
After multiple runs are complete, to produce one merged downstream table set plus diagnostics.

### Inputs required
1. Either explicit run paths (`--runs`) or discovery root (`--runs-root`)
2. Run directories containing `downstream_analysis/collation/*.csv`

### Run command
```bash
cd /home/bmorampa/rnamodbench

# discovery collation
python bin/downstream_analysis/collate_runs.py \
  --runs-root covbench_results \
  --run-glob 'results_*/downstream_analysis' \
  --output covbench_results/collated_multi_cov

# quiet mode + custom report prefix
python bin/downstream_analysis/collate_runs.py \
  --runs-root covbench_results \
  --run-glob 'results_*/downstream_analysis' \
  --output covbench_results/collated_multi_cov \
  --report-prefix covbench_collate \
  --quiet
```

### Optional args
1. `--tables <table...>` for selected table collation
2. `--strict` to return non-zero on soft issues
3. `--quiet` to suppress per-issue stderr details

### Exit behavior
1. `0`: complete (warnings may still exist)
2. `1`: strict mode + soft issues detected
3. `2`: hard failure (for example no usable runs)

### Expected outputs
1. Merged CSV tables in output dir
2. `<report_prefix>.md`
3. `<report_prefix>.json`

### How to verify success
```bash
cd /home/bmorampa/rnamodbench
python bin/downstream_analysis/collate_runs.py \
  --runs-root covbench_results \
  --run-glob 'results_*/downstream_analysis' \
  --output covbench_results/collated_multi_cov

test -f covbench_results/collated_multi_cov/collate_report.md
```

### If it fails, do this
1. Read CLI summary and `collate_report.md`
2. Inspect `collate_report.json` for issue codes (`missing_table_file`, `empty_file`, `csv_parse_error`, etc.)
3. Re-run missing failed coverages or regenerate downstream tables

### `bin/downstream_analysis/run_analysis.py`

### When to use
To re-run downstream analysis independently for a run without rerunning mapping/caller steps.

### Inputs required
1. `--output-dir` (required)
2. One of:
   - `--input-dir <modifications_dir>`
   - tool-specific args (`--tombo`, `--yanocomp`, `--nanocompore`, `--xpore`, `--eligos`, `--epinano`, `--differr`, `--drummer`, `--jacusa2`)
3. Optional but recommended:
   - `--ground-truth ground_truth_mod_positions.csv`
   - `--references-csv references.csv`

### Run command
```bash
cd /home/bmorampa/rnamodbench
python bin/downstream_analysis/run_analysis.py \
  --input-dir covbench_results/results_100x/modifications \
  --ground-truth ground_truth_mod_positions.csv \
  --references-csv references.csv \
  --output-dir covbench_results/results_100x/downstream_analysis_rerun \
  --run-id run_100x_rerun \
  --coverage-label 100x \
  --quality-label base \
  --expected-replicates 3 \
  --min-replicates 2
```

### Optional args
1. `--references` subset filter
2. `--threshold` score threshold
3. `--coverage-dirs` for `{coverage}/tool/` layout
4. `--verbose`

### Defaults
1. `--min-replicates 2`
2. `--expected-replicates 3`

### Exit behavior
1. Exits non-zero if input directory missing
2. Exits non-zero if no tool outputs loaded
3. Warns if ground truth file is missing but can still run without benchmark labels

### Expected outputs
Under output directory:
1. `metadata/` (`run_metadata.json`, quality reports)
2. `collation/` (`metrics_long.csv`, summary tables, window/lag tables)
3. by-reference analysis folders and visualizations

### How to verify success
```bash
cd /home/bmorampa/rnamodbench
ls covbench_results/results_100x/downstream_analysis_rerun/collation/
```

### If it fails, do this
1. If `ModuleNotFoundError: sklearn`, activate/install downstream-analysis dependencies
2. Confirm `--input-dir` exists and contains expected tool output structure
3. Re-run with `--verbose` for additional diagnostics

## 7) Verification Checklist Commands

Run this full checklist from repo root:

```bash
cd /home/bmorampa/rnamodbench

# CLI surface checks
python bin/downstream_analysis/collate_runs.py --help
python bin/downstream_analysis/run_analysis.py --help

# Script syntax checks
bash -n generate_samplesheets.sh
bash -n run_all_coverages.sh

# Discovery collation smoke test
python bin/downstream_analysis/collate_runs.py \
  --runs-root covbench_results \
  --run-glob 'results_*/downstream_analysis' \
  --output /tmp/collated_multi_cov_doccheck \
  --quiet

# Strict exit semantics smoke test
tmpd=$(mktemp -d)
mkdir -p "$tmpd/run_bad/collation"
python bin/downstream_analysis/collate_runs.py \
  --runs "$tmpd/run_bad" \
  --output "$tmpd/out" \
  --tables metrics_long.csv \
  --strict \
  --quiet
rc=$?
echo "strict_mode_exit=$rc"  # expected 1
```

## 8) Failure Recovery And Troubleshooting

### `run_all_coverages.sh` stops on one coverage
1. Open `covbench_logs/run_<COV>.log`
2. Identify failing process and root cause
3. Fix environment/data issue
4. Re-run `./run_all_coverages.sh` (completed runs skipped)

### `generate_samplesheets.sh` exits with missing FASTQ errors
1. Verify data roots in script match actual storage
2. Validate FASTQ naming pattern for each coverage and replicate
3. Re-run and confirm zero missing paths

### `collate_runs.py` reports missing/empty tables
1. Open generated `collate_report.md`
2. Inspect affected run/table paths
3. Re-run failed coverages or downstream analysis
4. Re-run collation

### `run_analysis.py` dependency error (`sklearn` missing)
1. Use environment from `modules/local/downstream_analysis/environment.yml`
2. Re-run command in that environment

### `run_analysis.py` reports no tool outputs loaded
1. Confirm input structure under `--input-dir`
2. If custom paths are used, provide explicit tool arguments

## 9) Audit And Handoff Checklist
Before sharing or archiving a completed batch, capture:

1. Commands used (exact shell history block)
2. Script versions/commit hash (`git rev-parse HEAD`)
3. Batch logs (`covbench_logs/run_*.log`)
4. Collation reports (`collate_report.md`, `collate_report.json`)
5. Final merged tables (`covbench_results/collated_multi_cov/*.csv`)
6. Per-run downstream metadata (`results_*/downstream_analysis/metadata/*`)

Suggested manifest command:

```bash
cd /home/bmorampa/rnamodbench
{
  echo "commit: $(git rev-parse HEAD)"
  echo "date: $(date -Iseconds)"
  echo "runs:"
  ls -1 covbench_results/results_* | sed 's#^#  - #' 
} > covbench_results/collated_multi_cov/run_manifest.txt
```

## 10) Known Pitfalls And Anti-Patterns
1. Running commands outside repo root and mixing relative paths
2. Editing data-root paths in one script but not documenting the change
3. Treating `run_all_coverages.sh` completion as sufficient without checking downstream `metrics_long.csv`
4. Ignoring `collate_report.md/json` and assuming merged outputs are complete
5. Running `run_analysis.py` in a Python env missing `scikit-learn`
6. Updating examples in multiple docs without updating this canonical runbook
