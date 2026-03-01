# Quick Start Guide

For current day-to-day operator workflows and copy/paste script commands, use the
canonical runbook: [`docs/operator_scripts_runbook.md`](docs/operator_scripts_runbook.md).

## What Has Been Created

A complete rnamodbench compliant pipeline structure for RNA modification detection has been created in `/home/bmorampa/rnamodbench/`.

## Pipeline Architecture

### ✅ Completed Components

1. **Main Entry Point** (`main.nf`)
   - Clean DSL2 structure
   - Workflow orchestration

2. **Main Workflow** (`workflows/rnamodbench.nf`)
   - Modular workflow design
   - Six subworkflows for logical grouping

3. **Six Subworkflows** (`subworkflows/local/`)
   - `input_check.nf` - Input validation and sample parsing
   - `mapping_rrna.nf` - Read mapping to 16S/23S rRNA
   - `qc_stats.nf` - Quality control and statistics
   - `prepare_signal_data.nf` - FAST5 processing and tombo resquiggle
   - `signal_processing.nf` - F5C eventalign
   - `modification_calling.nf` - Multi-tool modification detection

4. **Sample Modules** (`modules/local/`)
   - SAMPLESHEET_CHECK
   - MINIMAP2_ALIGN
   - SAMTOOLS_VIEW
   - TOMBO_RESQUIGGLE
   - F5C_EVENTALIGN
   - YANOCOMP_ANALYSIS

5. **Configuration Files**
   - `nextflow.config` - Main config with profiles
   - `conf/base.config` - Resource management
   - `conf/modules.config` - Module-specific settings
   - `nextflow_schema.json` - Parameter validation schema

6. **Documentation**
   - `README.md` - Complete pipeline documentation
   - `MODULES_README.md` - Module implementation guide
   - `IMPLEMENTATION_GUIDE.md` - Detailed implementation steps
   - `QUICK_START.md` - This file

7. **Helper Scripts**
   - `bin/check_samplesheet.py` - Sample sheet validation

8. **Example Files**
   - `assets/samplesheet.csv` - Example sample sheet
   - `assets/schema_input.json` - Input validation schema

## What Still Needs to Be Done

### 1. Implement Remaining Modules (~20 modules)

Create these files in `modules/local/` following the template in `MODULES_README.md`:

**High Priority:**
- [ ] samtools_sort.nf
- [ ] samtools_index.nf
- [ ] samtools_flagstat.nf
- [ ] samtools_depth.nf
- [ ] extract_mapped_reads.nf
- [ ] extract_read_ids.nf
- [ ] fast5_subset.nf
- [ ] multi_to_single_fast5.nf
- [ ] f5c_index.nf
- [ ] tombo_detect_modifications.nf
- [ ] tombo_text_output.nf
- [ ] yanocomp_prepare.nf
- [ ] xpore_dataprep.nf
- [ ] xpore_diffmod.nf

**Medium Priority:**
- [ ] nanoplot_bam.nf
- [ ] coverage_plot.nf

### 2. Create Helper Classes

Create `lib/` directory with:
- [ ] WorkflowMain.groovy
- [ ] WorkflowRnamodbench.groovy
- [ ] SchemaUtils.groovy
- [ ] TemplateUtils.groovy

### 3. Setup Test Data

- [ ] Create small test dataset
- [ ] Create `conf/test.config`
- [ ] Test pipeline end-to-end

### 4. Container Management

- [ ] Verify/update container paths for all modules
- [ ] Build custom containers for tools not in biocontainers
- [ ] Test with Singularity and Docker

## How to Complete the Pipeline

### Step 1: Implement a Module (Example: SAMTOOLS_SORT)

```bash
# Create the module file
nano /home/bmorampa/rnamodbench/modules/local/samtools_sort.nf
```

Copy this template:
```groovy
process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: bam
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools sort \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.sorted.bam \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
```

### Step 2: Test Individual Module

```bash
# Create a simple test script
cat > test_module.nf << 'EOF'
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SAMTOOLS_SORT } from './modules/local/samtools_sort'

workflow {
    Channel
        .fromPath('test.bam')
        .map { file ->
            def meta = [:]
            meta.id = 'test'
            [ meta, file ]
        }
        .set { ch_bam }

    SAMTOOLS_SORT(ch_bam)
}
EOF

# Run test
nextflow run test_module.nf -profile singularity
```

### Step 3: Repeat for All Modules

Follow the templates in `MODULES_README.md` for each module.

### Step 4: Test Complete Pipeline

```bash
# 1. Prepare your samplesheet
cat > samplesheet.csv << EOF
sample,fastq,type,replicate,fast5_dir
native_rep1,/path/to/native_rep1.fastq,native,rep1,/path/to/native_fast5
ivt_rep1,/path/to/ivt_rep1.fastq,ivt,rep1,/path/to/ivt_fast5
EOF

# 2. Run the pipeline
nextflow run /home/bmorampa/rnamodbench/main.nf \\
    -profile singularity \\
    --input samplesheet.csv \\
    --ref_16s /path/to/16S.fa \\
    --ref_23s /path/to/23S.fa \\
    --outdir results

# 3. Check outputs
ls -R results/
```

## Directory Structure Overview

```
rnamodbench/
├── main.nf                      # ✅ Entry point
├── nextflow.config              # ✅ Main config
├── nextflow_schema.json         # ✅ Parameter schema
├── README.md                    # ✅ Documentation
├── MODULES_README.md            # ✅ Module guide
├── IMPLEMENTATION_GUIDE.md      # ✅ Implementation details
├── QUICK_START.md              # ✅ This file
│
├── workflows/
│   └── rnamodbench.nf     # ✅ Main workflow
│
├── subworkflows/local/
│   ├── input_check.nf          # ✅ Complete
│   ├── mapping_rrna.nf         # ✅ Complete
│   ├── qc_stats.nf             # ✅ Complete
│   ├── prepare_signal_data.nf  # ✅ Complete
│   ├── signal_processing.nf    # ✅ Complete
│   └── modification_calling.nf # ✅ Complete
│
├── modules/local/
│   ├── samplesheet_check.nf    # ✅ Complete
│   ├── minimap2_align.nf       # ✅ Complete
│   ├── samtools_view.nf        # ✅ Complete
│   ├── tombo_resquiggle.nf     # ✅ Complete
│   ├── f5c_eventalign.nf       # ✅ Complete
│   ├── yanocomp_analysis.nf    # ✅ Complete
│   └── [16 more to create]     # ⚠️ TODO
│
├── conf/
│   ├── base.config             # ✅ Complete
│   └── modules.config          # ✅ Complete
│
├── bin/
│   └── check_samplesheet.py    # ✅ Complete
│
└── assets/
    ├── samplesheet.csv         # ✅ Example
    └── schema_input.json       # ✅ Schema
```

## Key Improvements vs Original Pipeline

| Feature | Old | New |
|---------|-----|-----|
| Lines in main.nf | 273 | 57 |
| Modularity | Low | High |
| Reusability | No | Yes |
| Configuration | Hardcoded | Flexible |
| Testing | Difficult | Easy |
| Maintenance | Hard | Easy |
| Scalability | Limited | Excellent |
| rnamodbench compliance | No | Yes |

## Expected Timeline

- **Phase 1** (2-3 days): Implement remaining modules
- **Phase 2** (1 day): Create helper classes
- **Phase 3** (1-2 days): Testing and debugging
- **Phase 4** (1 day): Container setup
- **Phase 5** (1 day): Documentation updates

**Total**: ~1 week for a complete, production-ready pipeline

## Getting Help

1. **Check existing documentation**:
   - `README.md` - Pipeline usage
   - `MODULES_README.md` - Module templates
   - `IMPLEMENTATION_GUIDE.md` - Detailed steps

2. **Nextflow resources**:
   - https://www.nextflow.io/docs/latest/
   - https://github.com/bhargava-morampalli/rnamodbench/docs/

3. **Module examples**:
   - Check `modules/local/` for completed examples
   - Look at rnamodbench modules: https://github.com/bhargava-morampalli/rnamodbench/modules

## Quick Reference Commands

```bash
# Run pipeline
nextflow run /home/bmorampa/rnamodbench/main.nf \\
    -profile singularity \\
    --input samplesheet.csv \\
    --ref_16s refs/16S.fa \\
    --ref_23s refs/23S.fa \\
    --outdir results

# Resume failed run
nextflow run ... -resume

# Clean work directory
nextflow clean -f

# Check pipeline syntax
nextflow run /home/bmorampa/rnamodbench/main.nf --help

# View DAG
nextflow run /home/bmorampa/rnamodbench/main.nf -with-dag flowchart.html
```

## Success Criteria

Your pipeline is complete when:
- ✅ All modules are implemented
- ✅ Test run completes successfully
- ✅ All outputs are generated correctly
- ✅ Documentation is updated
- ✅ Version tracking works
- ✅ Resume functionality works
- ✅ Resource management is optimal

## Next Steps

1. Start implementing modules from the high-priority list
2. Test each module individually as you create it
3. Once all modules are done, test the complete pipeline
4. Optimize resource usage
5. Add more modification detection tools as needed

Good luck with your implementation! The hard architectural work is done - now it's mostly following the module template for the remaining processes.
