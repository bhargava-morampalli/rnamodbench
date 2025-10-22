# Implementation Guide for nf-core/rnamodifications Pipeline

This document provides a comprehensive guide for completing and deploying the RNA modifications detection pipeline.

## Overview

This pipeline follows the latest nf-core standards (DSL2) with improved modularity, scalability, and maintainability. The structure has been completely reorganized from your original pipeline.

## Directory Structure

```
new_improved/
├── main.nf                          # Main pipeline entry point
├── nextflow.config                  # Main configuration file
├── nextflow_schema.json             # Parameter schema for validation
├── README.md                        # Pipeline documentation
├── MODULES_README.md                # Module implementation guide
├── IMPLEMENTATION_GUIDE.md          # This file
│
├── workflows/
│   └── rnamodifications.nf         # Main workflow logic
│
├── subworkflows/local/
│   ├── input_check.nf              # Input validation subworkflow
│   ├── mapping_rrna.nf             # Read mapping subworkflow
│   ├── qc_stats.nf                 # QC statistics subworkflow
│   ├── prepare_signal_data.nf      # Signal preparation subworkflow
│   ├── signal_processing.nf        # Signal processing subworkflow
│   └── modification_calling.nf     # Modification detection subworkflow
│
├── modules/local/
│   ├── samplesheet_check.nf        # ✓ Sample sheet validation
│   ├── minimap2_align.nf           # ✓ Read alignment
│   ├── samtools_view.nf            # ✓ SAM to BAM conversion
│   ├── samtools_sort.nf            # TO IMPLEMENT
│   ├── samtools_index.nf           # TO IMPLEMENT
│   ├── samtools_flagstat.nf        # TO IMPLEMENT
│   ├── samtools_depth.nf           # TO IMPLEMENT
│   ├── extract_mapped_reads.nf     # TO IMPLEMENT
│   ├── extract_read_ids.nf         # TO IMPLEMENT
│   ├── fast5_subset.nf             # TO IMPLEMENT
│   ├── multi_to_single_fast5.nf    # TO IMPLEMENT
│   ├── tombo_resquiggle.nf         # ✓ Signal resquiggling
│   ├── tombo_detect_modifications.nf # TO IMPLEMENT
│   ├── tombo_text_output.nf        # TO IMPLEMENT
│   ├── f5c_index.nf                # TO IMPLEMENT
│   ├── f5c_eventalign.nf           # ✓ Event alignment
│   ├── yanocomp_prepare.nf         # TO IMPLEMENT
│   ├── yanocomp_analysis.nf        # ✓ Modification detection
│   ├── xpore_dataprep.nf           # TO IMPLEMENT
│   ├── xpore_diffmod.nf            # TO IMPLEMENT
│   ├── nanoplot_bam.nf             # TO IMPLEMENT
│   └── coverage_plot.nf            # TO IMPLEMENT
│
├── conf/
│   ├── base.config                 # Base resource configuration
│   ├── modules.config              # Module-specific configuration
│   ├── test.config                 # TO CREATE: Test configuration
│   └── test_full.config            # TO CREATE: Full test configuration
│
├── bin/
│   └── check_samplesheet.py        # ✓ Sample sheet validation script
│
└── assets/
    ├── samplesheet.csv             # ✓ Example sample sheet
    └── schema_input.json           # ✓ Input schema
```

## Key Improvements Over Original Pipeline

### 1. Modular Architecture
- **Subworkflows**: Logical grouping of related processes
- **Reusable Modules**: Each tool is a standalone module
- **Clear Separation**: Workflows, subworkflows, and modules are clearly separated

### 2. Metadata Handling
- **meta Map**: All samples carry metadata (id, type, rrna, replicate)
- **Channel Operations**: Uses `.branch()`, `.join()`, `.groupTuple()` for complex operations
- **Type Safety**: Proper channel typing throughout

### 3. Configuration Management
- **Profile-based**: Easy switching between conda, docker, singularity
- **Resource Labels**: Automatic resource allocation (process_single/low/medium/high)
- **Module-specific configs**: Easy customization via conf/modules.config

### 4. Parameter Validation
- **nextflow_schema.json**: Automatic parameter validation
- **Input Schema**: Validates sample sheet structure
- **Error Handling**: Clear error messages for invalid inputs

### 5. Reproducibility
- **Version Tracking**: All tools report versions
- **Container Support**: Full support for Docker, Singularity, Podman
- **Conda Environments**: Alternative to containers

## Implementation Steps

### Phase 1: Complete Missing Modules (Priority)

Implement the remaining modules following the template in `MODULES_README.md`:

1. **SAM/BAM Processing** (High Priority)
   - samtools_sort.nf
   - samtools_index.nf
   - samtools_flagstat.nf
   - samtools_depth.nf

2. **Read/FAST5 Processing** (High Priority)
   - extract_mapped_reads.nf
   - extract_read_ids.nf
   - fast5_subset.nf (may need ont-fast5-api)
   - multi_to_single_fast5.nf (needs ont-fast5-api)

3. **Signal Processing** (High Priority)
   - f5c_index.nf
   - tombo_detect_modifications.nf
   - tombo_text_output.nf

4. **Modification Calling** (High Priority)
   - yanocomp_prepare.nf
   - xpore_dataprep.nf
   - xpore_diffmod.nf

5. **QC Modules** (Medium Priority)
   - nanoplot_bam.nf
   - coverage_plot.nf

### Phase 2: Add Missing Helper Classes

Create these files in `lib/` directory:

```groovy
// lib/WorkflowMain.groovy
class WorkflowMain {
    public static void initialise(workflow, params, log) {
        // Parameter validation
        // Help message
    }
}

// lib/WorkflowRnamodifications.groovy
class WorkflowRnamodifications {
    public static void initialise(params, log) {
        // Pipeline-specific initialization
    }
}

// lib/NfcoreSchema.groovy
class NfcoreSchema {
    public static Map paramsSummaryMap(workflow, params) {
        // Generate parameter summary
    }
}

// lib/NfcoreTemplate.groovy
class NfcoreTemplate {
    public static void email(workflow, params, summary_params, projectDir, log) {
        // Email notification
    }
    public static void summary(workflow, params, log) {
        // Pipeline summary
    }
}
```

### Phase 3: Testing

1. **Create Test Data**
   ```bash
   mkdir -p test_data/{native,ivt}/{fastq,fast5}
   mkdir -p test_data/references
   ```

2. **Create Test Config**
   Create `conf/test.config`:
   ```groovy
   params {
       config_profile_name = 'Test profile'
       config_profile_description = 'Minimal test dataset'
       max_cpus = 2
       max_memory = '6.GB'
       max_time = '6.h'

       input = 'test_data/samplesheet_test.csv'
       ref_16s = 'test_data/references/16S.fa'
       ref_23s = 'test_data/references/23S.fa'
   }
   ```

3. **Run Test**
   ```bash
   nextflow run main.nf -profile test,singularity
   ```

### Phase 4: Container Management

#### Option 1: Use Existing Containers
Most tools are available in biocontainers. Update module definitions:

```groovy
conda "bioconda::tool=version"
container "${ workflow.containerEngine == 'singularity' ?
    'https://depot.galaxyproject.org/singularity/tool:version' :
    'biocontainers/tool:version' }"
```

#### Option 2: Build Custom Containers
For tools not in biocontainers (yanocomp, custom scripts):

1. Create Dockerfile:
   ```dockerfile
   FROM continuumio/miniconda3:latest
   RUN conda install -c conda-forge -c bioconda yanocomp
   ```

2. Build and push:
   ```bash
   docker build -t your-registry/yanocomp:latest .
   docker push your-registry/yanocomp:latest
   ```

3. Convert to Singularity:
   ```bash
   singularity build yanocomp.sif docker://your-registry/yanocomp:latest
   ```

### Phase 5: Documentation

1. **Update README.md** with:
   - Actual container paths
   - Real test data locations
   - Example outputs
   - Troubleshooting section

2. **Create CHANGELOG.md**:
   ```markdown
   # Changelog

   ## v1.0.0 - Initial Release
   - Complete pipeline implementation
   - Support for tombo, yanocomp, xpore
   - Comprehensive QC metrics
   ```

3. **Add CITATIONS.md**:
   List all tools with proper citations

## Migration from Old Pipeline

### Key Differences

| Aspect | Old Pipeline | New Pipeline |
|--------|-------------|--------------|
| Structure | Flat, all in main.nf | Modular (workflows/subworkflows/modules) |
| Channels | Named channels | Meta-map based channels |
| Configuration | Hardcoded in processes | Centralized in conf/ |
| Parameters | Scattered | Validated via schema |
| Reusability | Low | High (modules can be reused) |
| Scalability | Limited | Excellent |

### Data Flow Comparison

**Old Pipeline:**
```
main.nf (273 lines)
├── Manual channel creation
├── Process definitions inline
├── Hardcoded paths
└── Manual file handling
```

**New Pipeline:**
```
main.nf (57 lines)
└── RNAMODIFICATIONS workflow
    ├── INPUT_CHECK subworkflow
    ├── MAPPING_RRNA subworkflow
    ├── QC_STATS subworkflow
    ├── PREPARE_SIGNAL_DATA subworkflow
    ├── SIGNAL_PROCESSING subworkflow
    └── MODIFICATION_CALLING subworkflow
```

## Usage Examples

### Basic Run
```bash
nextflow run /home/bmorampa/new_improved \
    -profile singularity \
    --input samplesheet.csv \
    --ref_16s references/k12_16S.fa \
    --ref_23s references/k12_23S.fa \
    --outdir results
```

### With Custom Parameters
```bash
nextflow run /home/bmorampa/new_improved \
    -profile singularity \
    --input samplesheet.csv \
    --ref_16s references/k12_16S.fa \
    --ref_23s references/k12_23S.fa \
    --minimap2_args '-ax splice -uf -k14' \
    --yanocomp_fdr_threshold 0.05 \
    --max_cpus 32 \
    --max_memory 256.GB \
    --outdir results
```

### Resume a Failed Run
```bash
nextflow run /home/bmorampa/new_improved \
    -profile singularity \
    --input samplesheet.csv \
    --ref_16s references/k12_16S.fa \
    --ref_23s references/k12_23S.fa \
    --outdir results \
    -resume
```

## Troubleshooting

### Common Issues

1. **Container Not Found**
   - Ensure singularity/docker is installed
   - Check container paths in modules
   - Try pulling containers manually

2. **Memory Issues**
   - Adjust `--max_memory` parameter
   - Check `conf/base.config` resource definitions
   - Use appropriate process labels

3. **FAST5 Processing Errors**
   - Verify FAST5 directory structure
   - Check ont-fast5-api installation
   - Ensure read IDs match between FASTQ and FAST5

4. **Channel Errors**
   - Verify meta map structure
   - Check channel branching logic
   - Use `.view()` to debug channels

## Next Steps

1. **Implement Remaining Modules** (See Phase 1)
2. **Create Helper Classes** (See Phase 2)
3. **Test Pipeline** (See Phase 3)
4. **Build/Collect Containers** (See Phase 4)
5. **Update Documentation** (See Phase 5)

## Additional Tools to Consider

Once the core pipeline is working, consider adding:

1. **Nanocompore**: Alternative to yanocomp
2. **EpiNano**: Another modification detection tool
3. **DRUMMER**: Deep learning-based modification detection
4. **m6Anet**: m6A-specific detection
5. **MINES**: Multi-instance learning approach
6. **Nanom6A**: Direct m6A detection

Each can be added as a new module in `modules/local/` and integrated via `modification_calling.nf` subworkflow.

## Support and Resources

- **Nextflow Documentation**: https://www.nextflow.io/docs/latest/
- **nf-core Guidelines**: https://nf-co.re/docs/
- **nf-core Modules**: https://nf-co.re/modules
- **nf-core Tools**: https://nf-co.re/tools

## Summary

This new pipeline structure provides:
- ✅ Modular, maintainable code
- ✅ Proper metadata handling
- ✅ Flexible configuration
- ✅ Resource management
- ✅ Parameter validation
- ✅ Container support
- ✅ Comprehensive documentation
- ✅ Scalability for adding new tools

The foundation is solid. Focus on implementing the remaining modules following the provided templates, and you'll have a production-ready, nf-core compliant pipeline!
