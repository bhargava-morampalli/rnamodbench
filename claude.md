# CLAUDE.md - nf-core Nextflow Pipeline Development Rules

> **Last Updated:** January 2026  
> **Purpose:** Rules and constraints for Claude Code when developing, modifying, or evaluating nf-core compliant Nextflow pipelines

---

## CRITICAL: Version Requirements (January 2026)

Before making ANY changes, verify you are using current versions:

| Component | Current Version | Notes |
|-----------|-----------------|-------|
| nf-core/tools | **3.5.1** | Commands require `nf-core pipelines` prefix |
| nf-schema | **2.6.1** | Replaces deprecated nf-validation |
| Nextflow | **≥25.04.0** (minimum), **25.10.0** (recommended) | Required for nf-schema 2.6.x |
| nf-test | **0.9.3** | pytest-workflow is deprecated |
| Python | **3.10+** | Python 3.9 is EOL |

### DEPRECATED - DO NOT USE:
- ❌ `nf-validation` plugin → Use `nf-schema`
- ❌ `Channel.fromSamplesheet()` → Use `Channel.fromList(samplesheetToList())`
- ❌ `params.validationXxx` → Use `validation { }` config scope
- ❌ `definitions` in JSON schema → Use `$defs`
- ❌ JSON Schema draft-07 → Use draft 2020-12
- ❌ pytest-workflow → Use nf-test
- ❌ `check_max()` function → Use `resourceLimits`
- ❌ `nf-core create` → Use `nf-core pipelines create`
- ❌ Gitpod → Use GitHub Codespaces with devcontainers

---

## 1. nf-schema Plugin (CURRENT STANDARD)

### 1.1 Installation in nextflow.config

```groovy
plugins {
    id 'nf-schema@2.6.1'
}

validation {
    parametersSchema     = "nextflow_schema.json"
    monochromeLogs       = false
    lenientMode          = false
    showHiddenParams     = false
    
    help {
        enabled = true
        command = "nextflow run <pipeline> --input samplesheet.csv --outdir results"
    }
    
    summary {
        beforeText = ""
        afterText  = ""
    }
    
    logging {
        unrecognisedParams = "warn"  // Options: skip, debug, info, warn, error
        unrecognisedHeaders = "warn"
    }
}
```

### 1.2 Usage in main.nf

```groovy
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// Validate parameters
validateParameters()

// Print parameter summary
log.info paramsSummaryLog(workflow)

// Parse samplesheet - THIS IS THE CORRECT WAY
ch_input = Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))
```

### 1.3 JSON Schema Format (draft 2020-12)

**nextflow_schema.json:**
```json
{
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "https://raw.githubusercontent.com/nf-core/<pipeline>/master/nextflow_schema.json",
    "title": "Pipeline Parameters",
    "type": "object",
    "$defs": {
        "input_output_options": {
            "title": "Input/output options",
            "type": "object",
            "required": ["input", "outdir"],
            "properties": {
                "input": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "schema": "assets/schema_input.json",
                    "description": "Path to samplesheet"
                },
                "outdir": {
                    "type": "string",
                    "format": "directory-path",
                    "description": "Output directory"
                }
            }
        }
    },
    "allOf": [
        {"$ref": "#/$defs/input_output_options"}
    ]
}
```

**Sample sheet schema (assets/schema_input.json):**
```json
{
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "type": "array",
    "items": {
        "type": "object",
        "required": ["sample", "fastq_1"],
        "properties": {
            "sample": {
                "type": "string",
                "pattern": "^\\S+$"
            },
            "fastq_1": {
                "type": "string",
                "format": "file-path",
                "exists": true
            },
            "fastq_2": {
                "type": "string",
                "format": "file-path",
                "exists": true
            }
        }
    },
    "uniqueEntries": ["sample"],
    "dependentRequired": {
        "fastq_2": ["fastq_1"]
    }
}
```

---

## 2. Pipeline Directory Structure

```
nf-core-<pipeline>/
├── main.nf                      # Entry point - REQUIRED
├── nextflow.config              # Central config - REQUIRED
├── nextflow_schema.json         # Parameter schema - REQUIRED
├── modules.json                 # Tracks installed modules
├── .nf-core.yml                 # nf-core tools configuration
├── README.md                    # Pipeline documentation - REQUIRED
├── LICENSE                      # MIT License - REQUIRED
├── CHANGELOG.md                 # Version history - REQUIRED
├── CITATIONS.md                 # How to cite
├── CODE_OF_CONDUCT.md           # Community guidelines
│
├── assets/
│   ├── schema_input.json        # Samplesheet schema
│   └── multiqc_config.yml       # MultiQC configuration
│
├── conf/
│   ├── base.config              # Default resource allocations
│   ├── modules.config           # Module-specific ext.args
│   ├── test.config              # Minimal test configuration
│   └── test_full.config         # Full-size test data
│
├── docs/
│   ├── usage.md                 # How to run - REQUIRED
│   └── output.md                # Output documentation - REQUIRED
│
├── modules/
│   ├── local/                   # Pipeline-specific modules
│   │   └── <module>/
│   │       ├── main.nf
│   │       ├── meta.yml
│   │       ├── environment.yml
│   │       └── tests/
│   └── nf-core/                 # Installed community modules
│
├── subworkflows/
│   ├── local/
│   └── nf-core/
│
├── workflows/
│   └── <pipeline>.nf            # Main workflow logic
│
├── tests/
│   ├── main.nf.test             # Pipeline-level nf-test
│   └── nextflow.config          # Test configuration
│
└── .github/
    └── workflows/               # CI/CD workflows
```

---

## 3. Module Development Rules

### 3.1 Module Structure

Each module MUST contain:
```
modules/<scope>/<tool>/
├── main.nf              # Process definition
├── meta.yml             # Documentation (YAML)
├── environment.yml      # Conda dependencies
└── tests/
    ├── main.nf.test     # nf-test file
    ├── main.nf.test.snap # Snapshot file
    └── nextflow.config  # Test config (optional)
```

### 3.2 Process Template (2026 Standard)

```groovy
process TOOL_SUBTOOL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community.wave.seqera.io/library/<tool>:<version>' :
        'community.wave.seqera.io/library/<tool>:<version>' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai, optional: true
    tuple val("${task.process}"), val('<tool>'), eval('<tool> --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    <tool> \\
        $args \\
        --threads $task.cpus \\
        --input $reads \\
        --output ${prefix}.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bai
    """
}
```

### 3.3 Version Emission with Topic Channels (NEW in 2025)

For NEW modules, use topic output qualifier:
```groovy
output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('samtools'), 
          eval('samtools --version | head -1 | sed "s/samtools //"'), 
          emit: versions, topic: versions
```

For modules with templates or complex versioning, continue using versions.yml but add topic:
```groovy
output:
    path "versions.yml", emit: versions, topic: versions
```

### 3.4 Naming Conventions

| Element | Convention | Example |
|---------|------------|---------|
| Directory | lowercase, no punctuation | `modules/nf-core/bwa/mem/` |
| Process name | UPPERCASE with underscores | `process BWA_MEM` |
| Parameters | snake_case | `params.output_format` |
| Functions | camelCase | `def getPrefix()` |
| Channels | snake_case lowercase | `ch_sorted_bam` |
| Output files | `${prefix}.<ext>` | `${prefix}.bam` |

### 3.5 Module Input/Output Rules

**MUST:**
- All input files defined in input channels
- Named extensions for ALL outputs: `path "*.txt", emit: txt`
- One output channel per file type
- Stub block for every module
- `$task.ext.args` for optional command arguments
- Use `task.cpus` for multi-threading

**MUST NOT:**
- Hardcode custom meta fields (only `meta.id` and `meta.single_end` allowed)
- Modify the `when` statement (use `process.ext.when` in config)
- Use `params` directly in modules (pass via channels or ext.args)

### 3.6 Resource Labels

```groovy
// In conf/base.config
process {
    withLabel:process_single { cpus = 1;  memory = 6.GB;  time = 4.h }
    withLabel:process_low    { cpus = 2;  memory = 12.GB; time = 4.h }
    withLabel:process_medium { cpus = 6;  memory = 36.GB; time = 8.h }
    withLabel:process_high   { cpus = 12; memory = 72.GB; time = 16.h }
    withLabel:process_long   { time = 20.h }
    withLabel:process_gpu    { accelerator = 1 }
}
```

---

## 4. Configuration Patterns

### 4.1 conf/modules.config

```groovy
process {
    withName: 'FASTQC' {
        ext.args   = '--quiet'
        ext.prefix = { "${meta.id}" }
        publishDir = [
            path: { "${params.outdir}/fastqc" },
            mode: params.publish_dir_mode,
            pattern: "*.{html,zip}"
        ]
    }
    
    withName: 'MULTIQC' {
        ext.args   = { params.multiqc_title ? "--title \"${params.multiqc_title}\"" : '' }
        publishDir = [
            path: { "${params.outdir}/multiqc" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
}
```

### 4.2 Standard Profiles

```groovy
profiles {
    docker {
        docker.enabled         = true
        docker.userEmulation   = true
        conda.enabled          = false
        singularity.enabled    = false
    }
    singularity {
        singularity.enabled    = true
        singularity.autoMounts = true
        conda.enabled          = false
        docker.enabled         = false
    }
    conda {
        conda.enabled          = true
        docker.enabled         = false
        singularity.enabled    = false
    }
    arm64 {
        docker.runOptions      = '-u $(id -u):$(id -g) --platform=linux/arm64'
    }
    test {
        includeConfig 'conf/test.config'
    }
    test_full {
        includeConfig 'conf/test_full.config'
    }
}
```

---

## 5. Testing with nf-test

### 5.1 Module Test Template

```groovy
nextflow_process {
    name "Test Process TOOL_NAME"
    script "../main.nf"
    process "TOOL_NAME"
    
    tag "modules"
    tag "modules_nfcore"
    tag "tool_name"
    
    config "./nextflow.config"

    test("sarscov2 - bam") {
        when {
            process {
                """
                input[0] = [
                    [ id:'test', single_end:false ],
                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam', checkIfExists: true)
                ]
                """
            }
        }
        then {
            assertAll(
                { assert process.success },
                { assert snapshot(process.out).match() }
            )
        }
    }
    
    test("sarscov2 - bam - stub") {
        options "-stub"
        when {
            process {
                """
                input[0] = [
                    [ id:'test', single_end:false ],
                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam', checkIfExists: true)
                ]
                """
            }
        }
        then {
            assertAll(
                { assert process.success },
                { assert snapshot(process.out).match() }
            )
        }
    }
}
```

### 5.2 Pipeline-Level Test

```groovy
// tests/main.nf.test
nextflow_pipeline {
    name "Test pipeline"
    script "../main.nf"

    test("Test with test profile") {
        when {
            params {
                outdir = "$outputDir"
            }
        }
        then {
            assertAll(
                { assert workflow.success },
                { assert snapshot(
                    path("$outputDir/multiqc/multiqc_report.html"),
                    path("$outputDir/pipeline_info/")
                ).match() }
            )
        }
    }
}
```

### 5.3 Test Configuration (tests/nextflow.config)

```groovy
process {
    withName: 'TOOL_NAME' {
        ext.args = params.module_args ?: ''
    }
}
```

### 5.4 nf-test Commands

```bash
# Run all tests
nf-test test

# Run specific test
nf-test test tests/modules/tool/main.nf.test

# Update snapshots
nf-test test --update-snapshot

# Run with specific profile
nf-test test --profile docker
```

---

## 6. meta.yml Documentation Structure (2026 Format)

```yaml
name: "tool_subtool"
description: Brief description of what the module does
keywords:
  - keyword1
  - keyword2
  - alignment  # Domain-specific
tools:
  - tool_name:
      description: |
        Full description of the tool
      homepage: https://tool-homepage.org
      documentation: https://tool-docs.org
      tool_dev_url: https://github.com/tool/repo
      doi: "10.xxxx/xxxxx"
      licence: ["MIT"]
      identifier: "biotools:tool_name"
      args_id: "$args"
  - samtools:  # If piped
      description: SAMtools for BAM processing
      homepage: http://www.htslib.org/
      identifier: "biotools:samtools"
      args_id: "$args2"

input:
  - - meta:
        type: map
        description: |
          Groovy Map containing sample information
          e.g. [ id:'sample1', single_end:false ]
    - reads:
        type: file
        description: Input FASTQ files
        pattern: "*.{fastq,fq,fastq.gz,fq.gz}"
        ontology: "http://edamontology.org/format_1930"
  - - reference:
        type: file
        description: Reference genome
        pattern: "*.{fa,fasta,fa.gz}"

output:
  - bam:
      - meta:
          type: map
          description: Sample metadata
      - "*.bam":
          type: file
          description: Aligned BAM file
          pattern: "*.bam"
          ontology: "http://edamontology.org/format_2572"
  - versions:
      - versions.yml:
          type: file
          description: Software versions
          pattern: "versions.yml"

authors:
  - "@github_username"
maintainers:
  - "@github_username"
```

---

## 7. Container Specifications

### 7.1 Container Directive Format

```groovy
conda "${moduleDir}/environment.yml"
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://community.wave.seqera.io/library/tool:1.0.0--abcdef123' :
    'community.wave.seqera.io/library/tool:1.0.0--abcdef123' }"
```

### 7.2 Container Registries (Priority Order)

1. `community.wave.seqera.io/library/` - Seqera Containers (preferred)
2. `quay.io/biocontainers/` - BioContainers
3. `biocontainers/` - Docker Hub BioContainers
4. `depot.galaxyproject.org/singularity/` - Galaxy Singularity

### 7.3 environment.yml Format

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - bioconda::samtools=1.18
  - bioconda::bwa=0.7.17
```

**Rules:**
- Pin to channel AND version: `bioconda::tool=1.0.0`
- DO NOT pin to build number (varies by platform)
- No `latest` or `dev` tags

---

## 8. Linting and Code Quality

### 8.1 Running Linting

```bash
# Lint entire pipeline
nf-core pipelines lint

# Lint for release
nf-core pipelines lint --release

# Lint specific module
nf-core modules lint <tool_name>

# Fix meta.yml issues automatically
nf-core modules lint --fix
```

### 8.2 .nf-core.yml Configuration

```yaml
repository_type: pipeline
nf_core_version: "3.5.1"

lint:
  # Disable specific checks
  actions_awsfulltest: false
  pipeline_todos: false
  
  # Skip specific files
  files_unchanged:
    - assets/email_template.html
    - CODE_OF_CONDUCT.md
  files_exist:
    - CODE_OF_CONDUCT.md

# Template features (set during creation)
template:
  name: mypipeline
  org: nf-core
  skip:
    - igenomes
```

### 8.3 Code Formatting

- Use Prettier for Markdown, YAML, JSON
- Use EditorConfig for whitespace
- Follow "Harshil Alignment" for Nextflow code

**Harshil Alignment Example:**
```groovy
process EXAMPLE {
    input:
    tuple val(meta), path(reads)
    path  reference
    val   mode

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tool \\
        --input  $reads \\
        --ref    $reference \\
        --mode   $mode \\
        --output ${prefix}.bam \\
        $args
    """
}
```

---

## 9. Release Process

### 9.1 Version Bump

```bash
nf-core pipelines bump-version 1.0.0
```

### 9.2 Pre-Release Checklist

1. All tests pass: `nf-test test`
2. Linting passes: `nf-core pipelines lint --release`
3. CHANGELOG.md updated with version and date
4. Schema validated: `nf-core pipelines schema build`
5. Documentation complete (usage.md, output.md)

### 9.3 Semantic Versioning

- **Major** (1.0→2.0): Breaking changes (renamed params, changed outputs)
- **Minor** (1.0→1.1): New features, template updates
- **Patch** (1.0.0→1.0.1): Bug fixes, documentation

---

## 10. Common Mistakes to Avoid

### ❌ WRONG - Using deprecated nf-validation
```groovy
plugins {
    id 'nf-validation@1.1.3'  // DEPRECATED
}
```

### ✅ CORRECT - Using nf-schema
```groovy
plugins {
    id 'nf-schema@2.6.1'
}
```

---

### ❌ WRONG - Old samplesheet parsing
```groovy
Channel.fromSamplesheet("input")  // DEPRECATED
```

### ✅ CORRECT - New samplesheet parsing
```groovy
Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))
```

---

### ❌ WRONG - Old schema format
```json
{
    "$schema": "http://json-schema.org/draft-07/schema",
    "definitions": { }
}
```

### ✅ CORRECT - New schema format
```json
{
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$defs": { }
}
```

---

### ❌ WRONG - Old validation params
```groovy
params.validationFailUnrecognisedParams = true
params.validationShowHiddenParams = false
```

### ✅ CORRECT - New validation config
```groovy
validation {
    logging.unrecognisedParams = "error"
    showHiddenParams = false
}
```

---

### ❌ WRONG - versions.yml only
```groovy
output:
    path "versions.yml", emit: versions
```

### ✅ CORRECT - Topic channel for versions (new modules)
```groovy
output:
    tuple val("${task.process}"), val('tool'), 
          eval('tool --version'), emit: versions, topic: versions
```

---

### ❌ WRONG - Old nf-core commands
```bash
nf-core create
nf-core lint
```

### ✅ CORRECT - New command structure
```bash
nf-core pipelines create
nf-core pipelines lint
```

---

## 11. Quick Reference Commands

```bash
# Pipeline creation and management
nf-core pipelines create                    # Create new pipeline
nf-core pipelines lint                      # Lint pipeline
nf-core pipelines lint --release            # Lint for release
nf-core pipelines schema build              # Build/update schema
nf-core pipelines bump-version 1.0.0        # Bump version
nf-core pipelines sync                      # Sync with template

# Module management
nf-core modules list remote                 # List available modules
nf-core modules install <tool>              # Install module
nf-core modules update --all                # Update all modules
nf-core modules create                      # Create new module
nf-core modules lint <tool>                 # Lint module
nf-core modules lint --fix                  # Auto-fix meta.yml

# Subworkflow management
nf-core subworkflows install <name>         # Install subworkflow
nf-core subworkflows create                 # Create subworkflow

# Testing
nf-test test                                # Run all tests
nf-test test --update-snapshot              # Update snapshots
nf-test test --profile docker               # Test with profile

# Running pipeline
nextflow run . -profile test,docker --outdir results
nextflow run . -profile test,singularity --outdir results
```

---

## 12. Documentation Links

Always verify against these authoritative sources:

- **nf-core docs**: https://nf-co.re/docs
- **nf-schema docs**: https://nextflow-io.github.io/nf-schema/2.6.1/
- **nf-schema migration**: https://nextflow-io.github.io/nf-schema/latest/migration_guide/
- **Nextflow docs**: https://www.nextflow.io/docs/latest/index.html
- **nf-test docs**: https://www.nf-test.com/
- **Module specifications**: https://nf-co.re/docs/guidelines/components/modules

---

## 13. Validation Checklist

Before submitting any pipeline code:

- [ ] Using nf-schema@2.6.1 (NOT nf-validation)
- [ ] JSON schemas use draft 2020-12 with `$defs`
- [ ] Samplesheet parsed with `samplesheetToList()`
- [ ] All modules have stub blocks
- [ ] nf-test tests exist and pass
- [ ] Linting passes: `nf-core pipelines lint`
- [ ] Resource labels applied to all processes
- [ ] Containers pinned to specific versions
- [ ] meta.yml follows 2026 structure with ontologies
- [ ] CHANGELOG.md updated
- [ ] Topic channels used for version reporting (new modules)
