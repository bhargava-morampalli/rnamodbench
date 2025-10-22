# Quick Commands to Run Tests

## Step 1: Stub Test (Dry Run - Fast!)

```bash
cd /home/bmorampa/new_improved

nextflow run main.nf \
    -profile test,singularity \
    --input test_data/samplesheet.csv \
    --ref_16s test_data/references/16S.fa \
    --ref_23s test_data/references/23S.fa \
    --outdir results_stub \
    -stub-run
```

## Step 2: Real Test (If stub test passes)

```bash
cd /home/bmorampa/new_improved

nextflow run main.nf \
    -profile test,singularity \
    --input test_data/samplesheet.csv \
    --ref_16s test_data/references/16S.fa \
    --ref_23s test_data/references/23S.fa \
    --outdir results_test \
    -resume \
    -with-report report.html \
    -with-timeline timeline.html
```

## Check for Errors

```bash
# Check the log file
tail -100 .nextflow.log

# Or search for errors
grep ERROR .nextflow.log
```

## Files Created:
- ✅ conf/test.config
- ✅ conf/test_full.config
- ✅ lib/WorkflowMain.groovy
- ✅ lib/NfcoreTemplate.groovy
- ✅ lib/NfcoreSchema.groovy
- ✅ lib/WorkflowRnamodifications.groovy
- ✅ lib/Utils.groovy
- ✅ workflows/rnamodifications.nf (added softwareVersionsToYAML function)

Now try the stub test!
