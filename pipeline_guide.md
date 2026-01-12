# Pipeline Guide: Data Flow & Usage

## 1. Data Flow Diagram

The following flowchart illustrates how the pipeline processes **3 Replicates** of Native and IVT data for RNA modification detection.

```mermaid
graph TD
    %% Inputs
    Input[Samplesheet CSV] --> |sample, type, replicate, fastq, fast5_dir| Parseinput(Parse Input)
    
    %% MAPPING_RRNA Subworkflow
    subgraph "1. MAPPING_RRNA"
        Parseinput --> |FASTQ| Minimap[MINIMAP2_ALIGN]
        Minimap --> |BAM| SamSplit[SAMTOOLS_VIEW (Split 16S/23S)]
        SamSplit --> |BAM 16S| Bam16[Sort & Index 16S]
        SamSplit --> |BAM 23S| Bam23[Sort & Index 23S]
    end

    %% PREPARE_SIGNAL_DATA Subworkflow
    subgraph "2. PREPARE_SIGNAL_DATA"
        Parseinput --> |FAST5 Directory| Fast5Subset[FAST5_SUBSET (Extract Read IDs)]
        Bam16 --> |BAM/BAI| Join16[Join FASTQ+BAM+FAST5]
        Bam23 --> |BAM/BAI| Join23[Join FASTQ+BAM+FAST5]
        
        Fast5Subset --> |Multi FAST5| Multi[MULTI_TO_SINGLE_FAST5]
        
        Fast5Subset --> |Multi FAST5| Join16
        Fast5Subset --> |Multi FAST5| Join23
    end

    %% SIGNAL_PROCESSING Subworkflow
    subgraph "3. SIGNAL_PROCESSING (F5C & Tombo Paths)"
        Join16 --> |Indexed| Eventalign16[F5C_INDEX & EVENTALIGN]
        Join23 --> |Indexed| Eventalign23[F5C_INDEX & EVENTALIGN]
        Multi --> |Single FAST5| Resquiggle[TOMBO_RESQUIGGLE]
    end

    %% MODIFICATION_CALLING Subworkflow
    subgraph "4. MODIFICATION_CALLING (Grouped by Replicate)"
        Eventalign16 --> |Eventalign Txt| ModCall16
        Eventalign23 --> |Eventalign Txt| ModCall23
        Resquiggle --> |Resquiggled FAST5| ModCallTombo
        
        subgraph "Modification Callers"
            ModCall16 --> Nanocompore
            ModCall16 --> Xpore
            ModCall16 --> Yanocomp
            ModCallTombo[Tombo Grouping] --> Detect[TOMBO_DETECT_MODIFICATIONS]
            Detect --> Text[TOMBO_TEXT_OUTPUT]
        end
    end

    %% Outputs
    Nanocompore --> Results[Modification Tables]
    Xpore --> Results
    Yanocomp --> Results
    Text --> Results
```

## 2. Test Run Guide

To perform a test run, you need to provide a **Samplesheet CSV**.

### A. Constructing the Samplesheet

Create a file named `samplesheet.csv` with the following columns:

| Column | Description | Example |
|--------|-------------|---------|
| `sample` | Unique sample ID | `native_rep1` |
| `type` | Condition type (`native` or `ivt`) | `native` |
| `replicate` | Replicate number (1, 2, 3) | `1` |
| `fastq` | Path to FASTQ file (gzipped) | `/data/native_rep1.fastq.gz` |
| `fast5_dir` | Path to directory containing raw FAST5s | `/data/native_rep1_fast5/` |

**Example Content (3 Replicates):**
```csv
sample,type,replicate,fastq,fast5_dir
native_rep1,native,1,/path/to/native1.fastq.gz,/path/to/native1/fast5/
native_rep2,native,2,/path/to/native2.fastq.gz,/path/to/native2/fast5/
native_rep3,native,3,/path/to/native3.fastq.gz,/path/to/native3/fast5/
ivt_rep1,ivt,1,/path/to/ivt1.fastq.gz,/path/to/ivt1/fast5/
ivt_rep2,ivt,2,/path/to/ivt2.fastq.gz,/path/to/ivt2/fast5/
ivt_rep3,ivt,3,/path/to/ivt3.fastq.gz,/path/to/ivt3/fast5/
```

### B. Running the Pipeline

Use the following command to run the pipeline. 

> **Note**: For a quick test, you can use the `-stub` profile to skip massive computation and just verify data flow. For real processing, use `-profile docker` (or singularity).

**Command:**
```bash
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --ref_16s /path/to/16S_ref.fasta \
    --ref_23s /path/to/23S_ref.fasta
```

### C. Troubleshooting Nanocompore

If `NANOCOMPORE` fails, ensure:
1.  **Docker is running**: Nanocompore requires a container.
2.  **Memory**: It can be memory intensive. Try increasing memory in `conf/base.config`.
3.  **Data Quality**: Ensure enough reads align to the references.
