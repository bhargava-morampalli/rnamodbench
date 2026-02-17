#!/bin/bash
# Generate samplesheets for all coverage levels
# Each samplesheet has 12 data rows: 3 reps × 2 types × 2 targets

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${SCRIPT_DIR}/samplesheets"
mkdir -p "$OUTDIR"

DATA_16S="/home/bmorampa/covbench/results_all_tools_16S/filtlong"
DATA_23S="/home/bmorampa/covbench/results_all_tools_23S/filtlong"
NATIVE_FAST5="/home/bmorampa/k12_native_fast5"
IVT_FAST5="/home/bmorampa/k12_ivt_fast5"

COVERAGES=(5x 10x 15x 20x 25x 30x 35x 40x 45x 50x 55x 60x 65x 70x 75x 80x 85x 90x 95x 100x 150x 200x 500x 750x 1000x)

errors=0

for COV in "${COVERAGES[@]}"; do
    OUTFILE="${OUTDIR}/samplesheet_${COV}.csv"
    echo "sample,fastq,type,replicate,fast5_dir,target" > "$OUTFILE"

    # Collect all fastq paths to verify
    all_ok=true

    for REP in 1 2 3; do
        # Native 16S — lowercase 16s in filename
        FASTQ="${DATA_16S}/native/${COV}/native_16s_${REP}_${COV}_filtlong.fastq.gz"
        if [[ ! -f "$FASTQ" ]]; then
            echo "ERROR: missing $FASTQ"
            all_ok=false
            errors=$((errors + 1))
        fi
        echo "native_16s_rep${REP},${FASTQ},native,rep${REP},${NATIVE_FAST5},16s" >> "$OUTFILE"

        # IVT 16S — uppercase 16S in filename
        FASTQ="${DATA_16S}/ivt/${COV}/ivt_16S_${REP}_${COV}_filtlong.fastq.gz"
        if [[ ! -f "$FASTQ" ]]; then
            echo "ERROR: missing $FASTQ"
            all_ok=false
            errors=$((errors + 1))
        fi
        echo "ivt_16s_rep${REP},${FASTQ},ivt,rep${REP},${IVT_FAST5},16s" >> "$OUTFILE"

        # Native 23S
        FASTQ="${DATA_23S}/native/${COV}/native_23S_${REP}_${COV}_filtlong.fastq.gz"
        if [[ ! -f "$FASTQ" ]]; then
            echo "ERROR: missing $FASTQ"
            all_ok=false
            errors=$((errors + 1))
        fi
        echo "native_23s_rep${REP},${FASTQ},native,rep${REP},${NATIVE_FAST5},23s" >> "$OUTFILE"

        # IVT 23S
        FASTQ="${DATA_23S}/ivt/${COV}/ivt_23S_${REP}_${COV}_filtlong.fastq.gz"
        if [[ ! -f "$FASTQ" ]]; then
            echo "ERROR: missing $FASTQ"
            all_ok=false
            errors=$((errors + 1))
        fi
        echo "ivt_23s_rep${REP},${FASTQ},ivt,rep${REP},${IVT_FAST5},23s" >> "$OUTFILE"
    done

    if $all_ok; then
        echo "OK: ${OUTFILE} ($(wc -l < "$OUTFILE") lines)"
    else
        echo "WARN: ${OUTFILE} has missing files"
    fi
done

echo ""
echo "Generated ${#COVERAGES[@]} samplesheets in ${OUTDIR}/"
if [[ $errors -gt 0 ]]; then
    echo "WARNING: ${errors} missing FASTQ file(s) detected!"
    exit 1
else
    echo "All FASTQ paths verified."
fi
