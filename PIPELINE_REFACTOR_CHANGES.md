# Pipeline Refactor: Dynamic Target-Based Data Flow

## Summary

This document describes the refactoring of the RNA modifications pipeline from a hardcoded dual-mapping approach (16S/23S) to a dynamic, target-based architecture that supports arbitrary rRNA targets.

## Problem with Previous Implementation

The previous implementation had these issues:

1. **Dual Mapping Overhead**: Every FASTQ was mapped to BOTH 16S and 23S references, regardless of content
2. **Double Remapping**: Extracted reads were remapped again to their respective references
3. **Hardcoded Targets**: All downstream code explicitly handled only `16s` and `23s` targets
4. **Poor Extensibility**: Adding new targets (5S, 5.8S, 18S, 28S) required extensive code changes
5. **Computational Waste**: ~2-4x mapping overhead if data was already target-specific

## New Implementation

### Key Changes

1. **Single Direct Mapping**: Each sample is mapped ONCE to its designated reference
2. **Dynamic Target Support**: Supports any rRNA target, not just 16S/23S
3. **Target Specified in Samplesheet**: Users declare target in the input CSV
4. **References CSV**: Flexible mapping of targets to reference files
5. **Unified Channels**: Downstream processing uses unified channels with dynamic grouping

---

## Files Changed

### 1. Input Schema (`assets/schema_input.json`)

**Change**: Added `target` column

```json
"target": {
    "type": "string",
    "pattern": "^[a-zA-Z0-9_]+[sS]?$",
    "errorMessage": "Target rRNA must be provided (e.g., 16s, 23s, 5s, 18s, 28s)",
    "meta": ["rrna"]
}
```

**New samplesheet format**:
```csv
sample,fastq,type,replicate,fast5_dir,target
sample1_16s,sample1_16s.fastq,native,rep1,/path/to/fast5,16s
sample1_23s,sample1_23s.fastq,native,rep1,/path/to/fast5,23s
sample2_16s,sample2_16s.fastq,ivt,rep1,/path/to/fast5,16s
sample2_23s,sample2_23s.fastq,ivt,rep1,/path/to/fast5,23s
```

### 2. References Schema (`assets/schema_references.json`) - NEW FILE

**Purpose**: Validates the references CSV file

```json
{
    "items": {
        "properties": {
            "target": { "type": "string" },
            "reference": { "type": "string", "format": "file-path" }
        },
        "required": ["target", "reference"]
    }
}
```

**References CSV format**:
```csv
target,reference
16s,/path/to/16s_reference.fa
23s,/path/to/23s_reference.fa
5s,/path/to/5s_reference.fa
```

### 3. Nextflow Schema (`nextflow_schema.json`)

**Change**: Added `references` parameter, deprecated `ref_16s`/`ref_23s`

```json
"references": {
    "type": "string",
    "format": "file-path",
    "pattern": "^\\S+\\.csv$",
    "description": "Path to CSV file mapping target names to reference FASTA files"
}
```

### 4. Config (`nextflow.config`)

**Change**: Added `references` parameter

```groovy
references = null  // CSV file mapping targets to reference FASTA files
ref_16s    = null  // [DEPRECATED]
ref_23s    = null  // [DEPRECATED]
```

### 5. Main Workflow (`workflows/rnamodbench.nf`)

**Changes**:
- Parses references CSV into a channel and map
- Parses samplesheet with `target` column (sets `meta.rrna`)
- Validates that all targets have corresponding references
- Passes unified channels to subworkflows
- Backwards compatibility: supports legacy `--ref_16s`/`--ref_23s` parameters

**Key new logic**:
```groovy
// Parse references CSV
ch_references = Channel.fromPath(params.references)
    .splitCsv(header: true)
    .map { row -> [ row.target.toLowerCase(), file(row.reference) ] }

// Parse samplesheet - target becomes meta.rrna
meta.rrna = row.target.toLowerCase()

// Pass unified channels
MAPPING_RRNA(reads, ch_references)
QC_STATS(MAPPING_RRNA.out.mapped_bams)
```

### 6. Mapping Subworkflow (`subworkflows/local/mapping_rrna.nf`)

**Complete rewrite**. Previous version:
- 184 lines
- 4 hardcoded MINIMAP2 aliases
- 4 hardcoded REMAP aliases
- Explicit 16s/23s branching

New version:
- 93 lines
- Single MINIMAP2 call
- Dynamic reference selection via map lookup
- No remapping step (not needed with pre-separated data)

**Key new logic**:
```groovy
// Convert references to map
ch_ref_map = references.collect().map { ref_list ->
    def ref_map = [:]
    ref_list.each { target, ref_file -> ref_map[target] = ref_file }
    ref_map
}

// Combine reads with matching reference
ch_reads_with_ref = reads.combine(ch_ref_map).map { meta, fastq, ref_map ->
    def reference = ref_map[meta.rrna.toLowerCase()]
    [ meta, fastq, reference ]
}

// Single mapping call
MINIMAP2_ALIGN(...)
```

### 7. QC Stats Subworkflow (`subworkflows/local/qc_stats.nf`)

**Simplified**. Previous version:
- 4 separate input channels (16s_native, 16s_ivt, 23s_native, 23s_ivt)

New version:
- Single unified `bams` channel
- Processes all BAMs regardless of target type

### 8. Prepare Signal Data (`subworkflows/local/prepare_signal_data.nf`)

**Simplified**. Previous version:
- 8 separate input channels
- 4 separate output channels with explicit branching

New version:
- 2 unified input channels (mapped_fastqs, mapped_bams)
- 2 unified output channels (single_fast5, f5c_ready)
- No explicit target branching

### 9. Signal Processing (`subworkflows/local/signal_processing.nf`)

**Changes**:
- Takes unified `f5c_ready` channel
- Takes `ref_map` for dynamic reference selection
- Outputs unified channels (consumers filter by `meta.rrna`)

**Key new logic**:
```groovy
// Dynamic reference selection
ch_inputs = f5c_ready.combine(ref_map).map { meta, fastq, bam, bai, fast5, refs ->
    def reference = refs[meta.rrna.toLowerCase()]
    [ meta, fastq, bam, bai, fast5, reference ]
}
```

### 10. Modification Calling (`subworkflows/local/modification_calling.nf`)

**Changes**:
- Takes unified channels for all inputs
- Takes `ref_map` for dynamic reference selection
- Groups by `${meta.rrna}_${meta.replicate}` (works for any target)
- Pairs native/IVT samples dynamically

**Key grouping logic** (works for any target):
```groovy
ch_bam_paired = bams
    .map { meta, bam, bai ->
        def key = "${meta.rrna}_${meta.replicate}"
        [ key, meta, bam, bai, meta.type ]
    }
    .groupTuple()
    .filter { ... has_native && has_ivt ... }
    .combine(ref_map)
    .map { ... select reference by meta.rrna ... }
```

---

## Migration Guide

### From Old to New Format

**Old command**:
```bash
nextflow run main.nf \
    --input old_samplesheet.csv \
    --ref_16s /path/to/16s.fa \
    --ref_23s /path/to/23s.fa \
    --outdir results
```

**New command**:
```bash
nextflow run main.nf \
    --input new_samplesheet.csv \
    --references references.csv \
    --outdir results
```

**Old samplesheet** (`old_samplesheet.csv`):
```csv
sample,fastq,type,replicate,fast5_dir
sample1,sample1.fastq,native,rep1,/path/to/fast5
sample2,sample2.fastq,ivt,rep1,/path/to/fast5
```

**New samplesheet** (`new_samplesheet.csv`):
```csv
sample,fastq,type,replicate,fast5_dir,target
sample1_16s,sample1_16s.fastq,native,rep1,/path/to/fast5,16s
sample1_23s,sample1_23s.fastq,native,rep1,/path/to/fast5,23s
sample2_16s,sample2_16s.fastq,ivt,rep1,/path/to/fast5,16s
sample2_23s,sample2_23s.fastq,ivt,rep1,/path/to/fast5,23s
```

**References CSV** (`references.csv`):
```csv
target,reference
16s,/path/to/16s.fa
23s,/path/to/23s.fa
```

---

## Benefits

| Aspect | Before | After |
|--------|--------|-------|
| Mapping operations | 4x per sample | 1x per sample |
| Supported targets | Only 16s, 23s | Any target |
| Adding new target | Major code changes | Just add to CSV |
| Code complexity | High (explicit branching) | Lower (dynamic grouping) |
| Input flexibility | Mixed data assumed | Pre-separated data |

---

## Input Validation

The pipeline performs two critical validations at startup:

### 1. Target-Reference Validation
Ensures every target in the samplesheet has a corresponding reference in the references CSV.

### 2. Native/IVT Pair Validation
Ensures every `target + replicate` combination has BOTH native AND IVT samples.

**Why this matters:** All downstream modification calling tools (Tombo, Yanocomp, Nanocompore, Xpore, ELIGOS, EpiNano, DiffErr, DRUMMER, JACUSA2) require comparison between native and IVT samples. Without both, no analysis can be performed.

**Example error if IVT is missing:**
```
ERROR: Target+replicate '16s_rep1' is missing ivt sample(s).
All modification calling tools require both native and IVT samples for comparison.
```

This validation prevents wasted computation by failing early with a clear error message.

---

## Backwards Compatibility

The pipeline maintains backwards compatibility:
- Legacy `--ref_16s` and `--ref_23s` parameters still work
- A deprecation warning is shown encouraging migration to `--references`
- Old parameters are converted internally to the new format

---

## Example: Adding 5S rRNA Support

With the new implementation, adding 5S support requires NO code changes:

1. Add to `references.csv`:
```csv
target,reference
16s,/path/to/16s.fa
23s,/path/to/23s.fa
5s,/path/to/5s.fa
```

2. Add samples to samplesheet:
```csv
sample,fastq,type,replicate,fast5_dir,target
sample1_5s,sample1_5s.fastq,native,rep1,/path/to/fast5,5s
sample2_5s,sample2_5s.fastq,ivt,rep1,/path/to/fast5,5s
```

3. Run pipeline - 5S samples will be processed automatically!

---

## nf-test Updates

All nf-tests have been updated to work with the new channel structures:

### Test Files Updated

1. **`tests/data/references.csv`** - NEW: Test references file
2. **`workflows/tests/data/samplesheet.csv`** - Updated with `target` column
3. **`subworkflows/local/mapping_rrna/tests/main.nf.test`** - Updated inputs
4. **`subworkflows/local/qc_stats/tests/main.nf.test`** - Updated inputs
5. **`subworkflows/local/prepare_signal_data/tests/main.nf.test`** - Updated inputs
6. **`subworkflows/local/signal_processing/tests/main.nf.test`** - Updated inputs
7. **`subworkflows/local/modification_calling/tests/main.nf.test`** - Updated inputs
8. **`workflows/tests/main.nf.test`** - Uses `references` param
9. **`tests/main.nf.test`** - Updated samplesheet and `references` param

### nf-test Configuration

Changed `nf-test.config` to use conda profile instead of docker:
```groovy
profile = "test,conda"
```

### Test Status

- **MAPPING_RRNA**: PASSES
- **QC_STATS**: PASSES
- **SIGNAL_PROCESSING**: PASSES
- **PREPARE_SIGNAL_DATA**: PASSES
- **MODIFICATION_CALLING**: Requires conda environment fixes (eligos2 dependency)

---

## Testing Recommendations

1. Test with existing 16S/23S data using new format
2. Test backwards compatibility with legacy parameters
3. Test with a third target type (e.g., 5S) to verify extensibility
4. Verify all modification calling tools produce correct output groupings
