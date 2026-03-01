# Tool Availability Reporting

This branch adds backfillable tool-availability and error reporting for completed
`results_*x` run directories.

## Scope

- additive reporting only
- no tool command changes
- no downstream metric changes
- no notebook-facing CSV schema changes
- no dependency on Nextflow `work/` directories

The reporting code reconstructs run-level availability from published artifacts:

- `pipeline_info/execution_trace_*.txt`
- `logs/`
- `modifications/`

Because deleted `work/` directories are not required, historical runs can be
backfilled without rerunning the benchmark.

## Per-run outputs

Running the report on a completed run writes these files into `pipeline_info/`:

- `process_status.tsv`
- `log_events.tsv`
- `tool_availability_per_run.tsv`
- `error_summary.csv`
- `error_summary.html`

## Collated outputs

When `collate_runs.py` is run over downstream directories whose parent runs have
per-run reporting files, it also writes:

- `tool_availability_matrix.tsv`
- `tool_failure_summary.tsv`
- `tool_log_events.tsv`
- `availability_report.md`

These outputs are additive and do not replace the six benchmark collation CSVs.

## Usage

Single completed run:

```bash
python /home/bmorampa/rnamodbench/bin/generate_error_report.py \
  --run-dir /home/bmorampa/rnamodbench/covbench_results/results_1000x
```

Backfill all historical coverage runs:

```bash
python /home/bmorampa/rnamodbench/bin/backfill_run_reports.py \
  --runs-root /home/bmorampa/rnamodbench/covbench_results \
  --run-glob 'results_*x'
```

Collate benchmark outputs and the new reporting artifacts into a separate output
directory:

```bash
python /home/bmorampa/rnamodbench/bin/downstream_analysis/collate_runs.py \
  --runs-root /home/bmorampa/rnamodbench/covbench_results \
  --run-glob 'results_*/downstream_20260301_parserfix_rernaeval' \
  --output /tmp/rnamodbench_tool_availability_collate_test
```

## Status model

`tool_availability_per_run.tsv` uses:

- `nf_status`: `COMPLETED`, `FAILED`, `ABORTED`, `PARTIAL`, `MISSING_TRACE`
- `output_state`: `parsed_nonempty`, `parsed_empty`, `raw_present_unparsed`,
  `final_output_missing`, `no_observed_artifacts`
- `reason_source`: `log`, `trace`, `output_scan`, `unknown`
- `diagnostic_completeness`: `trace_and_log`, `trace_only`, `log_only`,
  `output_only`, `limited`

This is designed to distinguish:

- process failure
- successful process with empty output
- successful process with missing published output
- output present but not parser-readable
- successful process with usable output

## Limits

Without historical `work/` directories, the report cannot reliably recover:

- `.command.err` / `.command.out`
- unpublished stderr
- full shell command state for failed tasks

When that information is unavailable, the report records inferred reasons and
their source explicitly instead of pretending the diagnosis is complete.
