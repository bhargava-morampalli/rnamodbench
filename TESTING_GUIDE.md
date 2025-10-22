# Complete Testing Guide for RNA Modifications Pipeline

## Step-by-Step Instructions

---

## STEP 1: Prepare Your Test Data

### 1.1 Create a Test Data Directory

```bash
cd /home/bmorampa/new_improved
mkdir -p test_data/fastq
mkdir -p test_data/fast5_native
mkdir -p test_data/fast5_ivt
mkdir -p test_data/references
```

### 1.2 Copy Your Test Data

Copy your smallest dataset for testing:

```bash
# Copy FASTQ files (native and IVT samples)
cp /path/to/your/native_sample.fastq.gz test_data/fastq/native_rep1.fastq.gz
cp /path/to/your/ivt_sample.fastq.gz test_data/fastq/ivt_rep1.fastq.gz

# Copy FAST5 directories
cp -r /path/to/your/native_fast5_dir/* test_data/fast5_native/
cp -r /path/to/your/ivt_fast5_dir/* test_data/fast5_ivt/

# Copy reference files
cp /path/to/16s_reference.fasta test_data/references/16s.fasta
cp /path/to/23s_reference.fasta test_data/references/23s.fasta
```

**Important:** Start with a SMALL subset (e.g., 1000-5000 reads) for initial testing!

---

## STEP 2: Verify Reference FASTA Files

### 2.1 Check Chromosome Names

Your reference FASTA files MUST have these EXACT header names:

```bash
# Check 16S reference
grep ">" test_data/references/16s.fasta
# Expected output: >16s_88_rrsE

# Check 23S reference
grep ">" test_data/references/23s.fasta
# Expected output: >23s_78_rrlB
```

### 2.2 Fix Chromosome Names If Needed

If your headers are different, fix them:

```bash
# For 16S (replace YOUR_CURRENT_HEADER with actual header)
sed -i 's/>YOUR_CURRENT_HEADER/>16s_88_rrsE/' test_data/references/16s.fasta

# For 23S (replace YOUR_CURRENT_HEADER with actual header)
sed -i 's/>YOUR_CURRENT_HEADER/>23s_78_rrlB/' test_data/references/23s.fasta
```

### 2.3 Index References

```bash
# Index for minimap2 (optional but recommended)
samtools faidx test_data/references/16s.fasta
samtools faidx test_data/references/23s.fasta
```

---

## STEP 3: Create Samplesheet

### 3.1 Create Samplesheet File

```bash
cat > test_data/samplesheet.csv << 'SAMPLESHEET'
sample,fastq,type,replicate,fast5_dir
native_rep1,test_data/fastq/native_rep1.fastq.gz,native,rep1,test_data/fast5_native
ivt_rep1,test_data/fastq/ivt_rep1.fastq.gz,ivt,rep1,test_data/fast5_ivt
SAMPLESHEET
```

### 3.2 Verify Samplesheet Format

```bash
# Check the samplesheet
cat test_data/samplesheet.csv

# Verify files exist
ls -lh test_data/fastq/
ls -d test_data/fast5_*
```

**Samplesheet Requirements:**
- **Header:** Must be exactly: `sample,fastq,type,replicate,fast5_dir`
- **sample:** Unique sample name (no spaces)
- **fastq:** Path to FASTQ file (can be .fastq, .fastq.gz, .fq, .fq.gz)
- **type:** Must be either `native` or `ivt`
- **replicate:** Replicate identifier (e.g., rep1, rep2)
- **fast5_dir:** Path to directory containing FAST5 files

---

## STEP 4: Configure Containers

### 4.1 Check Your Existing Containers

```bash
ls -lh /home/bmorampa/containers/*.sif
```

You have these containers:
- ✅ samtools_1.16.1--h6899075_1.sif
- ✅ minimap2_2.24--h7132678_1.sif
- ✅ f5c_1.1--h0326b38_1.sif
- ✅ ont-fast5-api_4.1.0--pyhdfd78af_0.sif
- ✅ tombo_new.sif (or ont-tombo_1.5.1--py36r36h39af1c6_2.sif)
- ✅ seqkit_2.3.1--h9ee0642_0.sif
- ✅ nanoplot_1.40.2--pyhdfd78af_0.sif
- ✅ xpore_2.1--pyh5e36f6f_0.sif
- ✅ yanocomp.sif

### 4.2 Update nextflow.config for Singularity

Add this to your `nextflow.config` (if not already present):

```bash
cat >> nextflow.config << 'CONFIG'

// Singularity configuration
singularity {
    enabled = true
    autoMounts = true
    cacheDir = '/home/bmorampa/containers/'
}

profiles {
    singularity {
        singularity.enabled = true
        singularity.autoMounts = true
    }
    
    test {
        params {
            max_cpus = 4
            max_memory = '16.GB'
            max_time = '6.h'
        }
    }
}
CONFIG
```

---

## STEP 5: Run Stub Test (Dry Run)

This tests that all modules are properly configured WITHOUT actually running them.

```bash
cd /home/bmorampa/new_improved

# Run stub test
nextflow run main.nf \
    -profile test,singularity \
    --input test_data/samplesheet.csv \
    --ref_16s test_data/references/16s.fasta \
    --ref_23s test_data/references/23s.fasta \
    --outdir results_stub \
    -stub-run \
    -resume
```

**What to look for:**
- ✅ No syntax errors
- ✅ All modules initialize
- ✅ Creates stub output files
- ✅ No container errors

**If it fails:** Check the error messages carefully. Common issues:
- Missing containers
- Incorrect paths in samplesheet
- Syntax errors in modules

---

## STEP 6: Run Small Test (Real Run)

Once stub test passes, run with actual data:

```bash
cd /home/bmorampa/new_improved

# Run pipeline with small dataset
nextflow run main.nf \
    -profile test,singularity \
    --input test_data/samplesheet.csv \
    --ref_16s test_data/references/16s.fasta \
    --ref_23s test_data/references/23s.fasta \
    --outdir results_test \
    -resume \
    -with-report report.html \
    -with-timeline timeline.html \
    -with-trace
```

**Monitoring the run:**
```bash
# In another terminal, watch progress
watch -n 5 'ls -lh work/ | tail -20'

# Check logs
tail -f .nextflow.log
```

---

## STEP 7: Troubleshooting Common Issues

### Issue 1: "Cannot find fast5 files"

**Fix:**
```bash
# Verify FAST5 files exist
ls test_data/fast5_native/*.fast5 | head -5
ls test_data/fast5_ivt/*.fast5 | head -5

# Check permissions
chmod -R 755 test_data/fast5_*
```

### Issue 2: "Container not found" or "Image file does not exist"

**Fix:**
```bash
# Check if container paths are correct
ls -lh /home/bmorampa/containers/*.sif

# The pipeline will try to pull from depot.galaxyproject.org
# If you want to use local containers, modify the module container paths
```

### Issue 3: "Tombo resquiggle fails"

**Possible causes:**
- Wrong chromosome names in reference
- FAST5 files not in correct format
- Not enough reads mapped

**Fix:**
```bash
# Verify chromosome names match
grep ">" test_data/references/16s.fasta
grep ">" test_data/references/23s.fasta

# Check BAM file has reads
samtools view -c results_test/mapping_rrna/sample_sorted.bam
```

### Issue 4: "Out of memory"

**Fix:**
```bash
# Reduce resource requirements in conf/base.config
# Or run with fewer reads in your test dataset
```

---

## STEP 8: Verify Output Files

### 8.1 Check Directory Structure

```bash
tree -L 2 results_test/
```

**Expected directories:**
```
results_test/
├── input_check/
├── mapping_rrna/
├── qc_stats/
├── signal_processing/
├── modification_calling/
└── pipeline_info/
```

### 8.2 Check Key Output Files

```bash
# Mapping outputs
ls results_test/mapping_rrna/

# QC outputs
ls results_test/qc_stats/

# Tombo outputs (CSV files with statistics)
ls results_test/modification_calling/tombo_*/*.csv

# Yanocomp outputs (BED files)
ls results_test/modification_calling/yanocomp_*/*.bed

# Xpore outputs
ls results_test/modification_calling/xpore_*/
```

### 8.3 Inspect Tombo CSV Output

```bash
# Check first few lines of Tombo output
head results_test/modification_calling/tombo_*/16s_rep1.csv

# Expected columns from get_reg_stats():
# - pos (position)
# - coverage
# - stat (test statistic)
# - possibly others depending on Tombo version
```

---

## STEP 9: Quick Validation Checks

### 9.1 Check Read Counts

```bash
# Count reads in FASTQ
zcat test_data/fastq/native_rep1.fastq.gz | echo $((`wc -l`/4))

# Count mapped reads
samtools view -c results_test/mapping_rrna/*native*16s*sorted.bam
```

### 9.2 Check Coverage

```bash
# Check depth files
head results_test/qc_stats/*depth*.txt
```

### 9.3 Check Modification Calls

```bash
# Tombo results
wc -l results_test/modification_calling/tombo_*/*.csv

# Look for positions with high test statistics
cat results_test/modification_calling/tombo_*/16s*.csv | \
    awk -F',' 'NR>1 {if($3>10) print}' | head -10
```

---

## STEP 10: Clean Up and Re-run

### 10.1 Clean Work Directory (Saves Space)

```bash
# After successful run
nextflow clean -f -k

# Or remove old work directories manually
rm -rf work/
```

### 10.2 Resume Failed Run

```bash
# Resume from where it failed
nextflow run main.nf \
    -profile test,singularity \
    --input test_data/samplesheet.csv \
    --ref_16s test_data/references/16s.fasta \
    --ref_23s test_data/references/23s.fasta \
    --outdir results_test \
    -resume
```

---

## STEP 11: Run with Full Dataset

Once small test works:

```bash
# Create full samplesheet
cat > samplesheet_full.csv << 'FULL'
sample,fastq,type,replicate,fast5_dir
native_rep1,/path/to/full/native_rep1.fastq.gz,native,rep1,/path/to/fast5/native_rep1
native_rep2,/path/to/full/native_rep2.fastq.gz,native,rep2,/path/to/fast5/native_rep2
ivt_rep1,/path/to/full/ivt_rep1.fastq.gz,ivt,rep1,/path/to/fast5/ivt_rep1
ivt_rep2,/path/to/full/ivt_rep2.fastq.gz,ivt,rep2,/path/to/fast5/ivt_rep2
FULL

# Run full pipeline
nextflow run main.nf \
    -profile singularity \
    --input samplesheet_full.csv \
    --ref_16s /path/to/16s.fasta \
    --ref_23s /path/to/23s.fasta \
    --outdir results_full \
    -resume \
    -with-report report_full.html \
    -with-timeline timeline_full.html
```

---

## Resource Requirements

### Minimum for Testing:
- **CPUs:** 4 cores
- **Memory:** 16 GB RAM
- **Storage:** 50 GB (for small test with ~5000 reads)

### Recommended for Full Run:
- **CPUs:** 16+ cores
- **Memory:** 64+ GB RAM
- **Storage:** 500 GB - 1 TB (depending on dataset size)

---

## Quick Reference Commands

```bash
# Check pipeline status
nextflow log

# Kill running pipeline
pkill -f nextflow

# View execution report
firefox report.html

# Check specific process logs
cat .nextflow.log | grep ERROR

# Find work directory for failed process
nextflow log last -f workdir,name,status | grep FAILED
cd <workdir>
cat .command.err
```

---

## Expected Runtime (Small Test Dataset ~5000 reads)

| Stage | Approximate Time |
|-------|-----------------|
| Input Check | 1-2 min |
| Mapping | 5-10 min |
| QC Stats | 2-5 min |
| Signal Processing | 10-20 min |
| Tombo Resquiggle | 20-30 min |
| Modification Calling | 10-20 min |
| **Total** | **~1-2 hours** |

Full dataset times will scale roughly linearly with read count.

---

## Contact & Help

If you encounter issues:
1. Check `.nextflow.log` for errors
2. Look at process-specific logs in `work/` directory
3. Verify input file formats and paths
4. Check container availability

---

**Good luck with testing! 🚀**
