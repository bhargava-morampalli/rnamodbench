# Modification Calling Tools: Command Comparison

Comparison between `pipeline_all_method.md` (reference document) and the current Nextflow pipeline implementation.

---

## 1. TOMBO (Signal-level native vs IVT comparison)

### Commands

**Reference (`pipeline_all_method.md`):**
```bash
tombo detect_modifications level_sample_compare \
    --fast5-basedirs wt/single/ \
    --alternate-fast5-basedirs ivt/single/ \
    --statistics-file-basename sample.level_compare_sample
```

**Current Pipeline (`modules/local/tombo_detect_modifications/main.nf`):**
```bash
tombo detect_modifications level_sample_compare \
    --fast5-basedirs ${native_dir} \
    --alternate-fast5-basedirs ${ivt_dir} \
    --statistics-file-basename ${prefix} \
    --store-p-value \
    --minimum-test-reads 1 \
    --statistic-type ks \
    --processes ${task.cpus}
```

### Differences

| Parameter | Reference | Pipeline | Impact |
|-----------|-----------|----------|--------|
| `--store-p-value` | Missing | Added | Pipeline saves p-values for downstream analysis (needed for ROC curves) |
| `--minimum-test-reads` | Missing (default=50) | `1` | Pipeline uses lower threshold - detects more positions but may be noisier |
| `--statistic-type` | Missing (default=ks) | `ks` | Same - both use Kolmogorov-Smirnov test |
| `--processes` | Missing | Added | Pipeline uses parallelization |

### Summary
Pipeline adds parameters for ROC curve generation and uses more lenient read threshold.

---

## 2. NANOCOMPORE (Eventalign-based sample comparison)

### Commands

**Reference (`pipeline_all_method.md`):**
```bash
nanocompore sampcomp \
    --file_list1 ./wt/eventalign_collapsed_reads_tsv/out_eventalign_collapse.tsv \
    --file_list2 ./ivt/eventalign_collapsed_reads_tsv/out_eventalign_collapse.tsv \
    --label1 wt \
    --label2 ivt \
    --fasta SINV_Toto1101.fa \
    --outpath ./nanocompore_results
```

**Current Pipeline (`modules/local/nanocompore_sampcomp/main.nf`):**
```bash
nanocompore sampcomp \
    --file_list1 ${cond1_files} \
    --file_list2 ${cond2_files} \
    --label1 ${condition1_label} \
    --label2 ${condition2_label} \
    --fasta $fasta \
    --outpath ${prefix} \
    --nthreads ${task.cpus} \
    --sequence_context 2 \
    --pvalue_thr 0.01 \
    --logit
```

### Differences

| Parameter | Reference | Pipeline | Impact |
|-----------|-----------|----------|--------|
| `--nthreads` | Missing | Added | Pipeline uses parallelization |
| `--sequence_context` | Missing (default=0) | `2` | Pipeline uses 5-mer context (2 bases on each side) for better accuracy |
| `--pvalue_thr` | Missing (default=0.05) | `0.01` | Pipeline uses stricter p-value threshold |
| `--logit` | Missing | Added | Pipeline uses logistic regression model (recommended for better performance) |

### Summary
Pipeline uses enhanced parameters recommended by benchmarking studies for improved accuracy.

---

## 3. XPORE (Differential modification detection)

### Commands

**Reference (`pipeline_all_method.md`):**
```bash
xpore dataprep --eventalign ivt.eventalign.txt --out_dir ivt_dataprep
xpore dataprep --eventalign wt.eventalign.txt --out_dir wt_dataprep

# Manual YAML config creation:
# data:
#     IVT:
#         rep1: ./IVT/IVT_dataprep/
#     wt:
#         rep1: ./wt/wt_dataprep/
# out: ./xpore_results

xpore diffmod --config xpore.yml
```

**Current Pipeline (`modules/local/xpore_diffmod/main.nf`):**
```bash
# Auto-generates YAML config
cat > config.yml << YAML
out: ${out_dir}
data:
  native:
    rep1: ${native_dir}
  ivt:
    rep1: ${ivt_dir}
YAML

xpore diffmod --config config.yml --n_processes ${task.cpus}
```

### Differences

| Parameter | Reference | Pipeline | Impact |
|-----------|-----------|----------|--------|
| `--n_processes` | Missing | Added | Pipeline uses parallelization |
| Config creation | Manual | Automatic | Pipeline auto-generates config from native/ivt directory inputs |

### Summary
Functionally equivalent - same core command, pipeline just automates config generation and adds parallelization.

---

## 4. ELIGOS (Alignment-based modification detection)

### Commands

**Reference (`pipeline_all_method.md`):**
```bash
eligos2 pair_diff_mod \
    -tbam wt.bam \
    -cbam ivt.bam \
    -reg gene.bed \
    -ref SINV_Toto1101.fa \
    -t 16 \
    --pval 1 \
    --oddR 0 \
    --esb 0 \
    -o eligos2_results
```

**Current Pipeline (`modules/local/eligos_pair_diff_mod/main.nf`):**
```bash
# Auto-generates BED file from reference
awk -v OFS='\t' '{print $1, 0, $2, $1, 0, "+"}' ${reference}.fai > regions.bed

eligos2 pair_diff_mod \
    -tbam $native_bam \
    -cbam $ivt_bam \
    -reg regions.bed \
    -ref $reference \
    -t $task.cpus \
    $args \
    -o $prefix
```

### Differences

| Parameter | Reference | Pipeline | Impact |
|-----------|-----------|----------|--------|
| `-reg` | `gene.bed` (user-provided) | Auto-generated from FASTA | Pipeline covers entire reference automatically; reference requires manual BED file |
| `--pval 1` | Explicit | Via `$args` | **Not hardcoded** - must be passed via config if needed |
| `--oddR 0` | Explicit | Via `$args` | **Not hardcoded** - must be passed via config |
| `--esb 0` | Explicit | Via `$args` | **Not hardcoded** - must be passed via config |

### Summary
**Key difference:** Reference document sets `--pval 1 --oddR 0 --esb 0` to output ALL positions for ROC analysis. Pipeline requires these to be passed via `task.ext.args` in modules.config.

**Recommendation:** Add these defaults to pipeline config for ROC curve compatibility.

---

## 5. EPINANO (Error pattern-based detection)

### Commands

**Reference (`pipeline_all_method.md`):**
```bash
# Variants extraction
python path/EpiNano/Epinano_Variants.py \
    -n 20 \
    -R SINV_Toto1101.fa \
    -b ivt.bam \
    -s path/EpiNano/misc/sam2tsv.jam \
    --type t

python path/EpiNano/Epinano_Variants.py \
    -n 20 \
    -R SINV_Toto1101.fa \
    -b wt.bam \
    -s path/EpiNano/misc/sam2tsv.jam \
    --type t

# Differential error (mismatch)
Rscript path/EpiNano/Epinano_DiffErr.R \
    -k ivt/ivt.plus_strand.per.site.csv \
    -w wt/wt.plus_strand.per.site.csv \
    -o epinano_mismatch_ \
    -f mis \
    -d 0.1 \
    -p
```

**Current Pipeline (`modules/local/epinano_error/main.nf`):**
```bash
# Variants extraction
python EpiNano/Epinano_Variants.py \
    -r $reference \
    -b $native_bam \
    -c $task.cpus \
    -o native_variants

python EpiNano/Epinano_Variants.py \
    -r $reference \
    -b $ivt_bam \
    -c $task.cpus \
    -o ivt_variants

# Differential error (sum_err)
Rscript EpiNano/Epinano_DiffErr.R \
    -k "ivt_fixed.csv" \
    -w "native_fixed.csv" \
    -t $zscore_threshold \
    -o ${prefix}/${meta.id} \
    -c $coverage_threshold \
    -f sum_err \
    -d $error_threshold
```

### Differences

| Parameter | Reference | Pipeline | Impact |
|-----------|-----------|----------|--------|
| `-n 20` vs `-c` | `-n 20` (threads) | `-c $task.cpus` | Same purpose, different flag name (newer API) |
| `-R` vs `-r` | `-R` (uppercase) | `-r` (lowercase) | API change - newer EpiNano uses lowercase |
| `-s sam2tsv.jar` | Required | Not used | **Newer EpiNano doesn't need sam2tsv** |
| `--type t` | Explicit | Not used | Pipeline uses default behavior |
| `-f mis` vs `-f sum_err` | `mis` (mismatch only) | `sum_err` | **Different metric** - pipeline uses combined errors |
| `-p` | Plot output | Not used | Pipeline skips plot generation |
| `-t` (z-score) | Not used | Added | Pipeline adds z-score threshold filter |
| `-c` (coverage) | Not used | Added | Pipeline adds coverage threshold filter |

### Summary
**Key difference:** Reference uses `-f mis` (mismatch only), pipeline uses `-f sum_err` (sum of mismatch + insertion + deletion rates). The `sum_err` metric may be more sensitive but could also be noisier.

**API changes:** Pipeline uses newer EpiNano API that doesn't require sam2tsv.jar.

---

## 6. DIFFERR (Differential error rate analysis)

### Commands

**Reference (`pipeline_all_method.md`):**
```bash
python main.py \
    -b wt/wt.bam \
    -a ivt.bam \
    -r SINV_Toto1101.fa \
    -o differr.bed
```

**Current Pipeline (`modules/local/differr/main.nf`):**
```bash
python -m differr \
    -a $ivt_bam \
    -b $native_bam \
    -r $reference \
    -o ${prefix}.bed \
    --fdr ${fdr_threshold}
```

### Differences

| Parameter | Reference | Pipeline | Impact |
|-----------|-----------|----------|--------|
| Invocation | `python main.py` | `python -m differr` | Pipeline uses module invocation (more robust, works after pip install) |
| `-a` and `-b` order | `-b wt -a ivt` | `-a ivt -b native` | Same assignment (IVT=alternate, native=reference) |
| `--fdr` | Not used | Default `1.0` | Pipeline outputs ALL positions for ROC analysis |

### Summary
Pipeline adds `--fdr 1.0` to ensure all positions are output for downstream ROC curve generation. The module invocation style is more portable.

---

## 7. DRUMMER (Statistical modification detection)

### Commands

**Reference (`pipeline_all_method.md`):**
```bash
python DRUMMER.py \
    -r SINV_Toto1101.fa \
    -t ivt.bam \
    -c wt.bam \
    -p 1 \
    -o DURMMER_result \
    -a exome \
    -n NC_001547.1 \
    -m True
```

**Current Pipeline (`modules/local/drummer/main.nf`):**
```bash
python DRUMMER/DRUMMER.py \
    -r $reference \
    -n $ref_name \
    -c $ivt_bam \
    -t $native_bam \
    -o ${prefix} \
    -a exome \
    -p ${pval_threshold}
```

### Differences

| Parameter | Reference | Pipeline | Impact |
|-----------|-----------|----------|--------|
| `-t` and `-c` | `-t ivt -c wt` | `-t native -c ivt` | **Reference has them REVERSED** - Pipeline is correct |
| `-n` | Hardcoded `NC_001547.1` | Auto-extracted from FASTA | Pipeline automatically gets reference name from .fai |
| `-m True` | Explicit | Not used | Missing multiprocessing flag in pipeline |
| `-p` | `1` | Default `1.0` | Same - output all positions |

### Summary
**Critical issue in reference document:** The `-t` (treatment) and `-c` (control) parameters are **inverted** in `pipeline_all_method.md`.
- Treatment (`-t`) should be the **modified** sample (native/wt)
- Control (`-c`) should be the **unmodified** sample (ivt)

The pipeline has the **correct** assignment. The reference document appears to have an error.

**Missing in pipeline:** The `-m True` flag for multiprocessing.

---

## Overall Summary Table

| Tool | Command Match | Key Differences | Action Needed |
|------|---------------|-----------------|---------------|
| **Tombo** | ~90% | Pipeline adds `--store-p-value`, `--minimum-test-reads 1` | None - pipeline enhancements are beneficial |
| **Nanocompore** | ~80% | Pipeline adds `--sequence_context 2`, `--logit`, stricter p-value | None - recommended parameters |
| **Xpore** | ~100% | Functionally identical, pipeline auto-generates config | None |
| **ELIGOS** | ~70% | Pipeline missing `--pval 1 --oddR 0 --esb 0` defaults | Consider adding to modules.config |
| **EpiNano** | ~60% | Different error metric (`sum_err` vs `mis`), newer API | Decide which metric is preferred |
| **DiffErr** | ~90% | Pipeline adds `--fdr 1.0` for full output | None |
| **DRUMMER** | ~80% | Reference has **inverted** treatment/control; pipeline is correct | Reference doc has error |

---

## Recommendations

1. **ELIGOS**: Add default args `--pval 1 --oddR 0 --esb 0` to `conf/modules.config` for full position output
2. **EpiNano**: Confirm whether `sum_err` or `mis` is the preferred metric for your analysis
3. **DRUMMER**: Consider adding `-m True` for multiprocessing support
4. **Reference doc**: Correct the DRUMMER command to swap `-t` and `-c` parameters

---

*Generated: 2026-01-12*
*Comparison between: `pipeline_all_method.md` and Nextflow pipeline modules*
