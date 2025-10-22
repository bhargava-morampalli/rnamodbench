#!/bin/bash

echo "========================================"
echo "RNA Modifications Pipeline - Quick Setup"
echo "========================================"
echo ""

# Step 1: Create directories
echo "Step 1: Creating test data directories..."
mkdir -p test_data/fastq
mkdir -p test_data/fast5_native
mkdir -p test_data/fast5_ivt
mkdir -p test_data/references
echo "✅ Directories created"
echo ""

# Step 2: Instructions for data
echo "Step 2: Copy your test data"
echo "----------------------------"
echo "Run these commands (replace paths with your actual data):"
echo ""
echo "# Copy FASTQ files:"
echo "cp /path/to/your/native.fastq.gz test_data/fastq/native_rep1.fastq.gz"
echo "cp /path/to/your/ivt.fastq.gz test_data/fastq/ivt_rep1.fastq.gz"
echo ""
echo "# Copy FAST5 directories:"
echo "cp -r /path/to/native_fast5/* test_data/fast5_native/"
echo "cp -r /path/to/ivt_fast5/* test_data/fast5_ivt/"
echo ""
echo "# Copy reference files:"
echo "cp /path/to/16s.fasta test_data/references/16s.fasta"
echo "cp /path/to/23s.fasta test_data/references/23s.fasta"
echo ""
read -p "Press ENTER when you've copied your data..."

# Step 3: Check references
echo ""
echo "Step 3: Checking reference chromosome names..."
if [ -f "test_data/references/16s.fasta" ]; then
    echo "16S reference header:"
    grep ">" test_data/references/16s.fasta
    echo ""
    read -p "Should this be '>16s_88_rrsE'? (y/n): " fix16s
    if [ "$fix16s" = "y" ]; then
        read -p "Enter current header (without >): " old_header
        sed -i "s/>$old_header/>16s_88_rrsE/" test_data/references/16s.fasta
        echo "✅ Fixed 16S header"
    fi
else
    echo "⚠️  16S reference not found at test_data/references/16s.fasta"
fi

if [ -f "test_data/references/23s.fasta" ]; then
    echo "23S reference header:"
    grep ">" test_data/references/23s.fasta
    echo ""
    read -p "Should this be '>23s_78_rrlB'? (y/n): " fix23s
    if [ "$fix23s" = "y" ]; then
        read -p "Enter current header (without >): " old_header
        sed -i "s/>$old_header/>23s_78_rrlB/" test_data/references/23s.fasta
        echo "✅ Fixed 23S header"
    fi
else
    echo "⚠️  23S reference not found at test_data/references/23s.fasta"
fi

# Step 4: Create samplesheet
echo ""
echo "Step 4: Creating samplesheet..."
cat > test_data/samplesheet.csv << 'SAMPLESHEET'
sample,fastq,type,replicate,fast5_dir
native_rep1,test_data/fastq/native_rep1.fastq.gz,native,rep1,test_data/fast5_native
ivt_rep1,test_data/fastq/ivt_rep1.fastq.gz,ivt,rep1,test_data/fast5_ivt
SAMPLESHEET

echo "✅ Samplesheet created at test_data/samplesheet.csv"
echo ""
cat test_data/samplesheet.csv
echo ""

# Step 5: Verify files
echo "Step 5: Verifying files..."
echo ""
echo "FASTQ files:"
ls -lh test_data/fastq/*.fastq.gz 2>/dev/null || echo "⚠️  No FASTQ files found"
echo ""
echo "FAST5 directories:"
ls -d test_data/fast5_* 2>/dev/null || echo "⚠️  No FAST5 directories found"
echo ""
echo "Reference files:"
ls -lh test_data/references/*.fasta 2>/dev/null || echo "⚠️  No reference files found"
echo ""

# Step 6: Ready to test
echo "========================================"
echo "Setup Complete!"
echo "========================================"
echo ""
echo "Next steps:"
echo "1. Run stub test:"
echo "   nextflow run main.nf -profile test,singularity --input test_data/samplesheet.csv --ref_16s test_data/references/16s.fasta --ref_23s test_data/references/23s.fasta --outdir results_stub -stub-run"
echo ""
echo "2. Run real test:"
echo "   nextflow run main.nf -profile test,singularity --input test_data/samplesheet.csv --ref_16s test_data/references/16s.fasta --ref_23s test_data/references/23s.fasta --outdir results_test -resume"
echo ""
echo "See TESTING_GUIDE.md for detailed instructions!"
