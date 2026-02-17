#!/bin/bash
# Run the rnamodbench pipeline sequentially for all coverage levels.
# Cleans work directory after each successful run to save disk space.
# Skips coverages that already have completed results.
#
# Usage:
#   ./run_all_coverages.sh              # run all coverages
#   ./run_all_coverages.sh 50x 100x     # run only specific coverages

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

SAMPLESHEET_DIR="${SCRIPT_DIR}/samplesheets"
RESULTS_BASE="${SCRIPT_DIR}/covbench_results"
LOG_DIR="${SCRIPT_DIR}/covbench_logs"
GROUND_TRUTH="${SCRIPT_DIR}/ground_truth_mod_positions.csv"
REFERENCES="${SCRIPT_DIR}/references.csv"

mkdir -p "$RESULTS_BASE" "$LOG_DIR"

ALL_COVERAGES=(5x 10x 15x 20x 25x 30x 35x 40x 45x 50x 55x 60x 65x 70x 75x 80x 85x 90x 95x 100x 150x 200x 500x 750x 1000x)

# Use command-line args if provided, otherwise run all
if [[ $# -gt 0 ]]; then
    COVERAGES=("$@")
else
    COVERAGES=("${ALL_COVERAGES[@]}")
fi

# Pre-flight checks
if [[ ! -f "$GROUND_TRUTH" ]]; then
    echo "ERROR: Ground truth file not found: $GROUND_TRUTH"
    echo "Copy it from /tmp/ground_truth_mod_positions.csv first."
    exit 1
fi

if [[ ! -d "$SAMPLESHEET_DIR" ]]; then
    echo "ERROR: Samplesheet directory not found: $SAMPLESHEET_DIR"
    echo "Run ./generate_samplesheets.sh first."
    exit 1
fi

echo "========================================"
echo "  RNAModBench Multi-Coverage Batch Run"
echo "========================================"
echo "Coverages: ${COVERAGES[*]}"
echo "Results:   ${RESULTS_BASE}/"
echo "Logs:      ${LOG_DIR}/"
echo "Started:   $(date)"
echo "========================================"
echo ""

completed=0
failed=0
skipped=0

for COV in "${COVERAGES[@]}"; do
    SAMPLESHEET="${SAMPLESHEET_DIR}/samplesheet_${COV}.csv"
    OUTDIR="${RESULTS_BASE}/results_${COV}"
    LOGFILE="${LOG_DIR}/run_${COV}.log"
    RUN_NAME="covbench_${COV}"

    echo "--- ${COV} ---"

    # Check samplesheet exists
    if [[ ! -f "$SAMPLESHEET" ]]; then
        echo "  SKIP: samplesheet not found (${SAMPLESHEET})"
        skipped=$((skipped + 1))
        continue
    fi

    # Skip if already completed (downstream_analysis dir exists with metrics)
    if [[ -f "${OUTDIR}/downstream_analysis/collation/metrics_long.csv" ]]; then
        echo "  SKIP: already completed (${OUTDIR})"
        skipped=$((skipped + 1))
        continue
    fi

    echo "  Starting at $(date)"
    echo "  Log: ${LOGFILE}"

    # Run Nextflow
    if nextflow run main.nf \
        --input "$SAMPLESHEET" \
        --references "$REFERENCES" \
        --outdir "$OUTDIR" \
        --run_downstream true \
        --ground_truth "$GROUND_TRUTH" \
        --downstream_run_id "run_${COV}" \
        --downstream_coverage_label "${COV}" \
        --downstream_quality_label "base" \
        -profile singularity \
        -name "$RUN_NAME" \
        > "$LOGFILE" 2>&1; then

        echo "  DONE at $(date)"
        completed=$((completed + 1))

        # Clean work directory to save disk space
        echo "  Cleaning work directory..."
        rm -rf work/
        echo "  Cleaned."
    else
        echo "  FAILED at $(date) — see ${LOGFILE}"
        failed=$((failed + 1))
        echo ""
        echo "Stopping batch due to failure on ${COV}."
        echo "To resume, fix the issue and re-run this script."
        echo "Completed runs will be skipped automatically."
        break
    fi

    echo ""
done

echo "========================================"
echo "  Batch Complete: $(date)"
echo "  Completed: ${completed}"
echo "  Skipped:   ${skipped}"
echo "  Failed:    ${failed}"
echo "========================================"
