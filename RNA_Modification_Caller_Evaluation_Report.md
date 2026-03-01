# Repository-Derived Report on RNA Modification Caller Evaluation

## Scope and Caveats

This report is derived strictly from the repository contents in `pipeline.nf`, `pipeline.conf`, `Scripts/postprocessing.R`, `Scripts/statistical_analysis.R`, `Scripts/plottingScript.R`, and the caller/support Dockerfiles under `Docker/`. It does not import column semantics from external tool documentation unless the repository itself makes the meaning explicit. Where the repository selects anonymous columns (`V1`, `V2`, and so on), the report keeps those indices and labels any biological meaning as repo-inferred rather than externally validated (`Scripts/postprocessing.R:18-646`).

The caller layer is the primary scope. Preprocessing, orchestration, and runtime-support tools are included only to explain how native outputs are produced, how coordinates are represented before postprocessing, and how operational figures such as running time are assembled (`pipeline.nf:50-1170`; `pipeline.conf:74-99`; `Scripts/plottingScript.R:1228-1419`).

Verification performed for this report:

- Every benchmarked caller listed below appears in `pipeline.nf` and in `Scripts/postprocessing.R` (`pipeline.nf:401-1170`; `Scripts/postprocessing.R:25-605`).
- Every native output file named below is matched to the pipeline process that writes it and the postprocessing branch that reads it (`pipeline.nf:401-1170`; `Scripts/postprocessing.R:25-605`).
- Every downstream statistic is tied to an explicit overlap, threshold, or PRROC computation in `Scripts/statistical_analysis.R` or to a plotting/stat-summary block in `Scripts/plottingScript.R` (`Scripts/statistical_analysis.R:125-390`; `Scripts/plottingScript.R:116-1419`).
- A repository-wide scan found descriptive statistics, PR curves, overlap heatmaps, and box/bar/histogram plots, but no formal inferential test calls such as `t.test`, `wilcox`, `fisher.test`, `chisq.test`, `kruskal.test`, `cor.test`, `lm`, or `glm`; the evaluation implemented here is descriptive rather than inferential (`Scripts/statistical_analysis.R:19-22,125-390`; `Scripts/plottingScript.R:117-258,528-1419`).

## Pipeline Overview

The benchmark pipeline is a Nextflow workflow that starts from a sample sheet and raw direct RNA FAST5 directories, converts multi-read FAST5 to single-read FAST5, extracts FASTQ, aligns reads with minimap2, performs signal-level preprocessing with Tombo and Nanopolish where needed, runs multiple m6A callers, standardizes each caller output into a BED-like schema, and then performs overlap-based benchmarking and figure generation (`pipeline.nf:50-1170`; `README.md:1-43`).

The repository wires fourteen RNA modification callers into the workflow:

- Single-condition or test-condition callers: `DENA`, `EpiNano-SVM`, `m6Anet`, `MINES`, `Nanom6A` (`pipeline.nf:401-442,568-683,1028-1078`).
- Pairwise or differential callers: `DiffErr`, `DRUMMER`, `ELIGOS`, `EpiNano-Error`, `Nanocompore`, `NanoDoc`, `Tombo` sample comparison, `xPore`, `Yanocomp` (`pipeline.nf:351-398,451-534,685-809,817-1170`).

The benchmark is intentionally permissive at execution time and stricter at evaluation time. `pipeline.conf` defaults `threshold = "relaxed"` for postprocessing, while `Scripts/statistical_analysis.R` applies a separate fixed default threshold table when deriving default-threshold performance metrics and PR rankings (`pipeline.conf:65-99`; `Scripts/postprocessing.R:619-633`; `Scripts/statistical_analysis.R:27-32,176-234`).

The reference space is not uniform before postprocessing:

- Genome-native or genome-like outputs: `DiffErr`, `ELIGOS`, `EpiNano-Error`, `EpiNano-SVM`, `NanoDoc`, `Tombo` comparison; `DRUMMER` is treated as direct coordinates in postprocessing but its exact native coordinate semantics remain ambiguous from the repo alone (`Scripts/postprocessing.R:25-40,117-135,181-195,277-338,365-395,437-519`; `pipeline.nf:451-809,817-858`).
- Transcriptome-native outputs requiring transcript-to-genome lift-over: `DENA`, `MINES`, `m6Anet`, and `Tombo` comparison explicitly use `transcriptToGenome()` in postprocessing (`Scripts/postprocessing.R:42-116,203-276,437-519,527-605`).
- Mixed coordinate handling: `Nanocompore` uses a column already named `genomicPos`; `xPore` keeps the numeric position from `diffmod.table` and only looks up chromosome/strand by gene ID; `Nanom6A` parses chromosome from the TSV and looks up strand by gene ID (`Scripts/postprocessing.R:159-179,340-364,396-436`).

`NanoDoc` is wired into the pipeline but disabled by default in `pipeline.conf`, so it is part of the repository benchmark design but not part of the default executable configuration (`pipeline.conf:82-99`).

## Caller Inventory

### `caller_inventory`

| tool | pipeline_process | execution_mode | comparison_type | reference_space_native | native_output_file | container_source | version_inferred_from_repo | default_evaluation_threshold_direction | default_threshold_used_in_statistics | binary_only_or_scored | notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DENA | `dena` | test condition only (`condition1`) | single-condition RRACH candidate scoring | transcript coordinates lifted to genome in postprocessing | `${test_condition}/dena/dena_label.tsv` | `Docker/dena-23_02_2022/Dockerfile` plus `pipeline.conf` `bproject/dena:v1` | dated snapshot `23_02_2022` with bundled `DENA` repo and model archive | maximize | `0.1` | scored | Uses RRACH candidate extraction, LSTM prediction, then postprocessing filters `V6 > threshold` and carries `Mod.Ratio` to column 6 (`pipeline.nf:568-603`; `Scripts/postprocessing.R:42-109`; `pipeline.conf:197-203`). |
| DiffErr | `differr` | two-condition | pairwise genome BAM differential calling | genome-native BED-like output | `differr/differrOut.bed` | `Docker/differr-0.2/Dockerfile` plus `pipeline.conf` `bproject/differr:v1` | `0.2` | minimize | `0.05` | scored | Native run uses permissive `params.differrFDR` at execution; postprocessing does not re-filter and converts native `-log10(FDR)` to `FDR` in standardized column 6 (`pipeline.nf:451-473`; `Scripts/postprocessing.R:25-40`; `pipeline.conf:47-49,176-182`). |
| DRUMMER | `drummer` | two-condition | pairwise comparison by chromosome | repo-treated as direct coordinates from `transcript_id` and `position`; exact semantic space is unclear from repo alone | `drummer/**/multiple_comp.txt` | `Docker/drummer-28_02_2022/Dockerfile` plus `pipeline.conf` `bproject/drummer:v1` | dated snapshot `28_02_2022` | minimize | `0.05` in `statistical_analysis.R`; `0.01` in `plottingScript.R` | scored | Output is already filtered at run time by `params.drummerPval`; postprocessing carries `max.G_padj` into column 6 without additional filtering. This is one of the repo’s explicit threshold inconsistencies (`pipeline.nf:817-852`; `Scripts/postprocessing.R:117-135`; `Scripts/statistical_analysis.R:30-32`; `Scripts/plottingScript.R:71-74`). |
| ELIGOS | `eligos` | two-condition | pairwise differential comparison on genome alignments | genome-native | `eligos/merged/minimap.sortG.1_vs_minimap.sortG.2_on_genome_combine.txt` | `Docker/eligos-2.1.0/Dockerfile` plus `pipeline.conf` `bproject/eligos:v1` | `2.1.0` | minimize | `0.0001` | scored | Downstream consumes only the merged-replicate comparison, not the split-replicate directory; postprocessing filters on adjusted p-value and odds ratio and writes adjusted p-value to column 6 (`pipeline.nf:488-529,1169`; `Scripts/postprocessing.R:181-195`; `pipeline.conf:183-189`). |
| EpiNano-Error | `epinanoError` | two-condition | pairwise error-profile differential analysis | genome-native plus/minus strand CSVs | `epinanoError/plus/diffErr.delta-sum_err.prediction.csv` and/or `epinanoError/minus/diffErr.delta-sum_err.prediction.csv` | `Docker/epinano-1.2/Dockerfile` plus `pipeline.conf` `bproject/epinano:v1` | directory `1.2`, tarball `Epinano1.2.1` | maximize | `0.1` | scored | Postprocessing keeps rows labeled `mod` by native z-score prediction and carries numeric `delta_sum_err` into standardized column 6 (`pipeline.nf:685-767`; `Scripts/postprocessing.R:277-304`; `pipeline.conf:53-57,211-217`; `Docker/epinano-1.2/Dockerfile:1-34`). |
| EpiNano-SVM | `epinanoSVM` | test condition only (`condition1`) | single-condition classifier on per-site variant features | genome-native | `${test_condition}/epinanoSVM/plus_mod_prediction...csv` and/or `minus_mod_prediction...csv` | `Docker/epinano-1.2/Dockerfile` plus `pipeline.conf` `bproject/epinano:v1` | directory `1.2`, tarball `Epinano1.2.1` | maximize | `0.5` | scored | Uses only `condition1` outputs; postprocessing filters on probability and RRACH k-mer membership and writes `ProbM` to column 6 (`pipeline.nf:610-683`; `Scripts/postprocessing.R:306-338`; `pipeline.conf:204-210`; `Docker/epinano-1.2/Dockerfile:1-34`). |
| m6Anet | `m6anet2` | test condition only (`condition1`) | single-condition transcriptome ML inference | transcript coordinates lifted to genome in postprocessing | `m6anet/data.result.csv` | `Docker/m6anet-v1.1.0/Dockerfile` plus `pipeline.conf` `bproject/m6anet:v1` | `v1.1.0` | maximize | `0.9` | scored | Inference is run only on `condition1`; postprocessing filters on probability, lifts transcript positions to genome, and writes `Prob_mod` to column 6 (`pipeline.nf:1028-1078`; `Scripts/postprocessing.R:527-598`; `pipeline.conf:95-99`; `Docker/m6anet-v1.1.0/Dockerfile:1-10`). |
| MINES | `mines` | test condition only (`condition1`) | single-condition calling from Tombo de novo signal summaries | transcript coordinates lifted to genome in postprocessing | `${test_condition}/mines/m6A_output_filename.bed` | `Docker/mines-23_02_2022/Dockerfile` plus `pipeline.conf` `bproject/mines:v1` | dated snapshot `23_02_2022` | binary | `N/A` | binary-only | Native column 7 is carried transiently but dropped from the standardized output, so downstream treats MINES as binary overlap-only with no general score column (`pipeline.nf:536-566`; `Scripts/postprocessing.R:203-269`; `pipeline.conf:190-196`). |
| Nanocompore | `nanocompore2` | two-condition | pairwise transcriptome eventalign comparison | mixed but already expressed as `genomicPos` in native TSV | `nanocompore/outnanocompore_results.tsv` | `Docker/nanocompore-1.0.3/Dockerfile` plus `pipeline.conf` `bproject/nanocompore:v1` | `1.0.3` | minimize | `0.01` | scored | Postprocessing filters on `GMM_logit_pvalue` and `abs(Logit_LOR) > 0.5`, then writes p-value to column 6 (`pipeline.nf:959-1026`; `Scripts/postprocessing.R:159-179`; `pipeline.conf:93-99`). |
| NanoDoc | `nanodoc` | two-condition | pairwise raw-signal comparison | genome-window text outputs, chromosome parsed from filename | `nanodoc/nanoDoc_results_*.txt` | `Docker/nanodoc-28_02_2022/Dockerfile` plus `pipeline.conf` `bproject/nanodoc:v2` | dated snapshot `28_02_2022` | maximize | `0.02` | scored | Present in the benchmark code path but disabled by default in `pipeline.conf`; postprocessing uses column 12 as score and assigns strand `*` (`pipeline.nf:775-809`; `Scripts/postprocessing.R:365-395`; `pipeline.conf:88-90,218-220`; `Docker/nanodoc-28_02_2022/Dockerfile:1-31`). |
| Nanom6A | `nanom6a` | test condition only (`condition1`) | single-condition site prediction across multiple probability cutoffs | mixed: chromosome parsed from TSV, strand looked up from gene BED | `${test_condition}/nanom6a/result_final/ratio.*.tsv` | `Docker/nanom6a-22_10_2021/Dockerfile` plus `pipeline.conf` `bproject/nanom6a:v2` | dated release `2021_10_22` | binary for default metrics; derived ranking for PR | `0.5` for default metrics, `max_thr` across files for PR | binary-plus-derived-score | Multiple probability-specific files are generated from `params.nanom6AP`. Postprocessing drops the numeric mod-ratio from the standardized BEDs, so default metrics are binary and PR ranking is reconstructed later from the highest threshold per bin (`pipeline.nf:401-428`; `Scripts/postprocessing.R:396-436`; `Scripts/statistical_analysis.R:150-156,257-315`; `pipeline.conf:41-42`). |
| Tombo (sample comparison) | `tombo3` / postprocessing key `tomboComparison` | two-condition | pairwise sample-level raw-signal comparison | transcript coordinates lifted to genome in postprocessing | `tomboComparison/sample.level_samp_comp_detect.statistic.plus.wig` | `Docker/tombo-1.5.1/Dockerfile` plus `pipeline.conf` `bproject/tombo:v4` | `1.5.1` | minimize | `0.05` | scored | Only the sample-comparison branch is benchmarked directly. Tombo de novo (`tombo2`) is upstream for MINES but not itself standardized as a benchmark caller (`pipeline.nf:372-398,1145-1158`; `Scripts/postprocessing.R:437-519`; `pipeline.conf:148-167`; `Docker/tombo-1.5.1/Dockerfile:1-17`). |
| xPore | `xpore2` | two-condition | pairwise transcriptome eventalign comparison | mixed: numeric position from native table plus chromosome/strand lookup by gene ID | `xpore/diffmod.table` | `Docker/xpore-2.0/Dockerfile` plus `pipeline.conf` `bproject/xpore:v1` | `2.0` | minimize | `0.05` | scored | Postprocessing recomputes BH-adjusted FDR from native column 5, filters by that FDR, and uses gene BED only for chromosome and strand lookup (`pipeline.nf:900-957`; `Scripts/postprocessing.R:340-357`; `pipeline.conf:92-93`; `Docker/xpore-2.0/Dockerfile:1-10`). |
| Yanocomp | `yanocomp2` | two-condition | pairwise GMM test on eventalign-derived HDF5 | genome BED output | `yanocomp/yanocomp_output.bed` | `Docker/yanocomp-0.2/Dockerfile` plus `pipeline.conf` `bproject/yanocomp:v1` | `0.2` | minimize | `0.05` | scored | Native run uses permissive `params.yanocompFDR`; postprocessing carries native column 9 into standardized column 6 without further filtering (`pipeline.nf:1080-1140`; `Scripts/postprocessing.R:143-157`; `pipeline.conf:44-45,246-251`; `Docker/yanocomp-0.2/Dockerfile:1-18`). |

## Native Output Columns Consumed Per Caller

### `native_to_standardized_columns`

| tool | native_file | native_columns_read | native_column_names_if_available | native_columns_used_for_filtering | native_columns_used_for_coordinate_projection | native_columns_preserved_into_standardized_output | standardized_output_columns_written | standardized_columns_used_in_downstream_analysis | unused_native_columns_after_postprocessing | transformations_applied |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DENA | `dena_label.tsv` | `V1`, `V2`, `V6` via `data_dena[, c(1,2,2,6,6)]` | no header; repo treats `V1` as transcript ID, `V2` as transcript position, `V6` as mod ratio (`Scripts/postprocessing.R:44-48`) | `V6 > filtering_parameter` (`Scripts/postprocessing.R:48-49`) | `V1` and `V2` feed `IRanges` and `transcriptToGenome()` (`Scripts/postprocessing.R:50-92`) | `V6` survives as `Mod.Ratio`; coordinates and strand come from lift-over (`Scripts/postprocessing.R:105-109`) | `Chr`, `Start`, `End`, `Strand`, `Status`, `Mod.Ratio` | `Chr`, `Start`, `End`, `Strand` for overlaps; `Mod.Ratio` as column 6 score in stats and plotting; `Status` not used as an analysis variable (`Scripts/statistical_analysis.R:169-177`; `Scripts/plottingScript.R:145-172,191-207`) | All native columns except selected `V1`, `V2`, `V6` are discarded | Filters positives, lifts transcript sites to genome, deduplicates lift-over collisions, writes numeric ratio to standardized column 6 (`Scripts/postprocessing.R:48-109`). |
| DiffErr | `differrOut.bed` | `V1`, `V2`, `V3`, `V5`, `V6` (`Scripts/postprocessing.R:27-33`) | no header; repo treats `V5` as `-log10(FDR)` and `V6` as strand | none in postprocessing; row set is taken as-is (`Scripts/postprocessing.R:25-40`) | `V1`, `V2`, `V3`, `V6` become direct coordinates and strand | `V5` becomes standardized `FDR`; `V1`, `V2`, `V3`, `V6` become coordinates/strand | `Chr`, `Start`, `End`, `Strand`, `Status`, `FDR` | `Chr`, `Start`, `End`, `Strand`, `FDR` used downstream; `Status` only labels output (`Scripts/statistical_analysis.R:173-234`; `Scripts/plottingScript.R:145-172`) | Any native columns other than 1,2,3,5,6 | Converts start to 1-based by adding 1 to `V2`, keeps end unchanged, sets `Status="Mod"` for all rows, converts `10^(-V5)` into FDR (`Scripts/postprocessing.R:28-33`). |
| DRUMMER | `multiple_comp.txt` files under `drummer/DRUMMER/**/` | `transcript_id`, `position`, `max.G_padj` (`Scripts/postprocessing.R:117-135`) | header present | none in postprocessing (`Scripts/postprocessing.R:126-135`) | `transcript_id` and `position` are used directly; no lift-over | `max.G_padj` survives as `Padj` | `Chr`, `Start`, `End`, `Strand`, `Status`, `Padj` | `Chr`, `Start`, `End`, `Strand`, `Padj`; `Status` is not downstream analytical input (`Scripts/statistical_analysis.R:173-234`; `Scripts/plottingScript.R:145-172`) | All other DRUMMER columns | Concatenates all chromosome-level result files, duplicates `position` for `Start`/`End`, assigns `Strand="*"`, sets `Status="Mod"` for all rows (`Scripts/postprocessing.R:117-135`). |
| ELIGOS | `minimap.sortG.1_vs_minimap.sortG.2_on_genome_combine.txt` | columns `1:4`, `16`, `18` via `data_eligos[, c(1:4,18,16,18)]` (`Scripts/postprocessing.R:183-194`) | header present; uses `start_loc`, `end_loc`, `adjPval`, `oddR` | `adjPval < filtering_parameter` and `oddR > 1.2` (`Scripts/postprocessing.R:189-192`) | direct genome coordinates from columns 1-4 | `adjPval` survives as `Padj` | `Chr`, `Start`, `End`, `Strand`, `Status`, `Padj` | `Chr`, `Start`, `End`, `Strand`, `Padj`; `Status` only marks retained rows (`Scripts/statistical_analysis.R:173-234`; `Scripts/plottingScript.R:145-172`) | All non-selected columns; `oddR` is used only for filtering and then discarded | Adds 1 to `start_loc`, filters by adjusted p-value and odds ratio, keeps only modified rows, writes adjusted p-value to column 6 (`Scripts/postprocessing.R:187-195`). |
| EpiNano-Error | `plus/diffErr.delta-sum_err.prediction.csv` and/or `minus/diffErr.delta-sum_err.prediction.csv` | `chr_pos`, `z_score_prediction`, `delta_sum_err` (`Scripts/postprocessing.R:279-298`) | header present | `z_score_prediction == "mod"` (`Scripts/postprocessing.R:294-298`) | coordinates and strand are parsed from `chr_pos` text field (`Scripts/postprocessing.R:290-295`) | `delta_sum_err` survives as standardized column 6 | `Chr`, `Start`, `End`, `Strand`, `Status`, `Delta sum err` | `Chr`, `Start`, `End`, `Strand`, `Delta sum err`; `Status` is not reused analytically (`Scripts/statistical_analysis.R:173-234`; `Scripts/plottingScript.R:145-172`) | Any additional native EpiNano columns | Merges plus/minus files if both exist, parses `chr_pos`, filters by native z-score label, writes numeric delta error to column 6 (`Scripts/postprocessing.R:277-304`). |
| EpiNano-SVM | `plus_mod_prediction...csv` and/or `minus_mod_prediction...csv` | `V1`, `V2`, `V3`, `V4`, `V28` (`Scripts/postprocessing.R:306-332`) | no header; repo treats `V1` as k-mer, `V2` as position string, `V3` as chromosome, `V4` as strand, `V28` as probability | `V28 > filtering_parameter` and `V1 %in% rrach` (`Scripts/postprocessing.R:319-330`) | `V3` direct chromosome, start parsed from `V2`, strand from `V4` | `V28` survives as `ProbM` | `Chr`, `Start`, `End`, `Strand`, `Status`, `ProbM` | `Chr`, `Start`, `End`, `Strand`, `ProbM`; `Status` only records the filter result (`Scripts/statistical_analysis.R:173-234`; `Scripts/plottingScript.R:145-172`) | All non-selected native columns; `Kmer` is filter-only and discarded | Parses position from `V2`, shifts start/end by `+2`, filters to RRACH motifs and probability-positive rows, writes probability to column 6 (`Scripts/postprocessing.R:319-332`). |
| m6Anet | `data.result.csv` | native columns 1, 2, 4 via `data_m6anet[,1]`, `[,2]`, `[,4]` (`Scripts/postprocessing.R:529-537`) | header present in file load, but the code accesses by column index rather than explicit names | native probability column 4 `> filtering_parameter` (`Scripts/postprocessing.R:533-537`) | transcript ID and position feed `transcriptToGenome()` (`Scripts/postprocessing.R:538-581`) | probability survives as `Prob_mod` | `Chr`, `Start`, `End`, `Strand`, `Status`, `Prob_mod` | `Chr`, `Start`, `End`, `Strand`, `Prob_mod` drive overlaps and score-based evaluation; `Status` is not a separate downstream variable (`Scripts/statistical_analysis.R:173-234`; `Scripts/plottingScript.R:145-172`) | Other m6Anet result columns | Adds `+2` to position before lift-over, filters on probability, deduplicates lift-over collisions, writes probability to column 6 (`Scripts/postprocessing.R:529-598`). |
| MINES | `m6A_output_filename.bed` | `V1`, `V2`, `V3`, `V7` via `data_mines[, c(1:3,7,7)]` (`Scripts/postprocessing.R:203-209`) | no header; repo treats `V1` as transcript ID and `V7` as a modification ratio-like value | none in postprocessing (`Scripts/postprocessing.R:203-269`) | `V1`, `V2`, `V3` feed `transcriptToGenome()` (`Scripts/postprocessing.R:210-254`) | native `V7` is carried transiently but dropped before final write | `Chr`, `Start`, `End`, `Strand`, `Status` | only `Chr`, `Start`, `End`, `Strand` participate downstream; there is no standardized score column for the scored-tool path (`Scripts/statistical_analysis.R:239-255`; `Scripts/plottingScript.R:174-207`) | All other native MINES columns, including the numeric `V7` after transient use | Adds 1 to start, lifts transcript coordinates to genome, deduplicates collisions, assigns `Status="Mod"`, but drops the native score from the final BED (`Scripts/postprocessing.R:205-269`). |
| Nanocompore | `outnanocompore_results.tsv` | native columns `2`, `3`, `5`, `7`, `13` via `[, c(2,3,3,5,7,13,7)]` (`Scripts/postprocessing.R:159-173`) | header present; code uses `genomicPos`, `GMM_logit_pvalue`, `Logit_LOR` | `GMM_logit_pvalue < filtering_parameter` and `abs(Logit_LOR) > 0.5` (`Scripts/postprocessing.R:165-170`) | direct coordinates from `genomicPos`, direct strand from column 5 | p-value survives as standardized `Pvalue`; `Logit_LOR` is filter-only | `Chr`, `Start`, `End`, `Strand`, `Status`, `Pvalue` | `Chr`, `Start`, `End`, `Strand`, `Pvalue`; `Status` only flags retained rows (`Scripts/statistical_analysis.R:173-234`; `Scripts/plottingScript.R:145-172`) | All non-selected native columns; `Logit_LOR` after filter | Adds 2 to `genomicPos`, converts `NC` log-odds to `NA`, filters by p-value and effect-size magnitude, writes p-value to column 6 (`Scripts/postprocessing.R:161-173`). |
| NanoDoc | `nanoDoc_results_*.txt` | `V1`, `V12`, filename `x` via `data_nanodoc[, c(1,1,12,12,13)]` (`Scripts/postprocessing.R:374-386`) | no header; repo treats `V12` as a score and extracts chromosome from filename | `V12 > filtering_parameter` (`Scripts/postprocessing.R:375-379`) | `V1` direct start/end, chromosome parsed from filename, strand forced to `*` | `V12` survives as `Score` | `Chr`, `Start`, `End`, `Strand`, `Status`, `Score` | `Chr`, `Start`, `End`, `Strand`, `Score`; `Status` is only the filter label (`Scripts/statistical_analysis.R:173-234`; `Scripts/plottingScript.R:145-172`) | All other NanoDoc text columns | Concatenates all region text files, filters on column 12, extracts chromosome from filename regex, writes score to column 6, and writes nothing if no positive rows survive (`Scripts/postprocessing.R:365-395`). |
| Nanom6A | `ratio.*.tsv` | first field plus per-cell parsed subfields; the code extracts `geneID`, `chr`, `start`, `mod_ratio` from tokenized strings (`Scripts/postprocessing.R:396-426`) | no fixed header contract; parsing is custom | none in postprocessing beyond non-empty cell checks (`Scripts/postprocessing.R:402-418`) | `chr` is parsed directly from field 1, `start` from each populated cell; strand is looked up by `geneID` in gene BED (`Scripts/postprocessing.R:404-423`) | no numeric score survives; mod ratio is parsed but dropped | `Chr`, `Start`, `End`, `Strand`, `Status` | only `Chr`, `Start`, `End`, `Strand` are used downstream from standardized BEDs; default metrics use the `ratio 0.5` file as binary hits and PR ranking is reconstructed across files later (`Scripts/statistical_analysis.R:150-156,257-315`; `Scripts/plottingScript.R:62-64`) | Parsed `geneID` and `mod_ratio` after strand lookup | Produces one standardized BED per probability file, keeps only coordinates/strand/status, and drops the parsed `mod_ratio` from the final BED (`Scripts/postprocessing.R:396-436`). |
| Tombo (sample comparison) | `sample.level_samp_comp_detect.statistic.plus.wig` | two parsed columns from WIG body; `transcriptID`, `start/end`, and native p-value are reconstructed row-wise (`Scripts/postprocessing.R:437-519`) | header exists, but the code parses rows manually rather than by named columns | `10^(-native_pvalue) < filtering_parameter` (`Scripts/postprocessing.R:463-465`) | transcript ID and position feed `transcriptToGenome()` (`Scripts/postprocessing.R:466-509`) | native p-value survives as standardized `Pvalue` | `Chr`, `Start`, `End`, `Strand`, `Status`, `Pvalue` | `Chr`, `Start`, `End`, `Strand`, `Pvalue`; `Status` is not used separately downstream (`Scripts/statistical_analysis.R:173-234`; `Scripts/plottingScript.R:145-172`) | Other WIG/header content outside the parsed transcript-position-pvalue pattern | Parses transcript IDs from non-numeric lines, parses p-values from numeric lines, filters by p-value, lifts to genome, and adds 1 to final end coordinate (`Scripts/postprocessing.R:439-519`). |
| xPore | `diffmod.table` | native columns 1, 2, 5 (`Scripts/postprocessing.R:342-357`) | header present but code accesses by index; column 5 is treated as p-value | `BH-adjusted native column 5 < filtering_parameter` (`Scripts/postprocessing.R:346-350`) | position comes directly from native column 2; chromosome and strand are looked up from gene BED using column 1 gene ID (`Scripts/postprocessing.R:343-356`) | BH-adjusted FDR survives as `FDR` | `Chr`, `Start`, `End`, `Strand`, `Status`, `FDR` | `Chr`, `Start`, `End`, `Strand`, `FDR`; `Status` only records the filter result (`Scripts/statistical_analysis.R:173-234`; `Scripts/plottingScript.R:145-172`) | Other native xpore columns | Recomputes BH FDR from native p-values, adds `+2` to position, looks up chromosome and strand from gene BED, writes FDR to column 6 (`Scripts/postprocessing.R:342-357`). |
| Yanocomp | `yanocomp_output.bed` | `V1`, `V2`, `V3`, `V6`, `V9` via `data_yanocomp[, c(1:3,6,9,9)]` (`Scripts/postprocessing.R:143-151`) | no header; repo treats `V6` as strand and `V9` as score | none in postprocessing (`Scripts/postprocessing.R:143-157`) | direct coordinates and strand | `V9` survives as `Score` | `Chr`, `Start`, `End`, `Strand`, `Status`, `Score` | `Chr`, `Start`, `End`, `Strand`, `Score`; `Status` only labels retained rows (`Scripts/statistical_analysis.R:173-234`; `Scripts/plottingScript.R:145-172`) | Any non-selected native BED columns | Adds 2 to start, forces end equal to shifted start, sets `Status="Mod"` for all rows, writes native score to column 6 (`Scripts/postprocessing.R:145-151`). |

### Critical Answer: Exact Output Column Used Downstream Per Tool

This is the exact per-tool answer to “which output column is actually used downstream for PR curves and score-based evaluation?” For the score-bearing tools, `Scripts/statistical_analysis.R` always reads standardized output column 6 for PR curves, threshold sweeps, and default-threshold binarization (`Scripts/statistical_analysis.R:173-234`). `Scripts/plottingScript.R` also uses standardized column 6 for score histograms and default-threshold filtering when the BED has more than five columns (`Scripts/plottingScript.R:154-173`).

| tool | native output file | native output column used downstream | native column number from tool output | standardized output column used downstream | downstream analyses using that column | exact note |
| --- | --- | --- | --- | --- | --- | --- |
| DENA | `dena_label.tsv` | `V6` (`Mod.Ratio` after postprocessing) | `6` | standardized column `6` = `Mod.Ratio` | PR curves, manual threshold sweep, default-threshold recall/precision/F1, score histogram, default-threshold filtering in plotting | Native `V6` is the only numeric confidence column preserved into downstream analysis (`Scripts/postprocessing.R:44-49,105-109`; `Scripts/statistical_analysis.R:176-234`). |
| DiffErr | `differrOut.bed` | native `V5` (`-log10(FDR)`), transformed to `FDR` | `5` | standardized column `6` = `FDR` | PR curves, manual threshold sweep, default-threshold metrics, score histogram, default-threshold plotting filter | Downstream does not use native column 5 directly; it uses standardized column 6, which is `10^(-V5)` (`Scripts/postprocessing.R:27-33`; `Scripts/statistical_analysis.R:176-234`). |
| DRUMMER | `multiple_comp.txt` | `max.G_padj` | numeric position not inferable from repo alone; accessed by header name | standardized column `6` = `Padj` | PR curves, manual threshold sweep, default-threshold metrics, score histogram, default-threshold plotting filter | The repo uses the header name `max.G_padj`, not a numeric index. Numeric column position may depend on DRUMMER output layout (`Scripts/postprocessing.R:126-135`). |
| ELIGOS | `minimap.sortG.1_vs_minimap.sortG.2_on_genome_combine.txt` | `adjPval` | `18` | standardized column `6` = `Padj` | PR curves, manual threshold sweep, default-threshold metrics, score histogram, default-threshold plotting filter | Native column `16` (`oddR`) is filter-only; native column `18` (`adjPval`) is the one carried forward into downstream scoring (`Scripts/postprocessing.R:183-195`). |
| EpiNano-Error | `diffErr.delta-sum_err.prediction.csv` | `delta_sum_err` | numeric position not inferable from repo alone; accessed by header name | standardized column `6` = `Delta sum err` | PR curves, manual threshold sweep, default-threshold metrics, score histogram, default-threshold plotting filter | Native `z_score_prediction` is used only to retain `mod` rows; the downstream score is `delta_sum_err` (`Scripts/postprocessing.R:290-298`). |
| EpiNano-SVM | `plus_mod_prediction...csv` / `minus_mod_prediction...csv` | `V28` (`ProbM`) | `28` | standardized column `6` = `ProbM` | PR curves, manual threshold sweep, default-threshold metrics, score histogram, default-threshold plotting filter | Native `V1` k-mer is filter-only for RRACH membership; the downstream score is `V28` (`Scripts/postprocessing.R:319-332`). |
| m6Anet | `data.result.csv` | native column 4 (`Prob_mod`) | `4` | standardized column `6` = `Prob_mod` | PR curves, manual threshold sweep, default-threshold metrics, score histogram, default-threshold plotting filter | The code accesses this by numeric index, not by named column (`Scripts/postprocessing.R:529-537,594-597`). |
| MINES | `m6A_output_filename.bed` | none retained downstream | native `V7` is parsed but dropped before write | no standardized score column | binary overlap only; recall, precision, F1 only | MINES does not contribute a score column to PR curves or score histograms (`Scripts/postprocessing.R:205-209,265-269`; `Scripts/statistical_analysis.R:239-255`). |
| Nanocompore | `outnanocompore_results.tsv` | `GMM_logit_pvalue` | `7` as accessed in the repo subset | standardized column `6` = `Pvalue` | PR curves, manual threshold sweep, default-threshold metrics, score histogram, default-threshold plotting filter | `Logit_LOR` is also used, but only as a filter; the score carried into downstream analysis is `GMM_logit_pvalue` (`Scripts/postprocessing.R:161-173`). |
| NanoDoc | `nanoDoc_results_*.txt` | `V12` (`Score`) | `12` | standardized column `6` = `Score` | PR curves, manual threshold sweep, default-threshold metrics, score histogram, default-threshold plotting filter | Native filename is used only to recover chromosome; score comes from `V12` (`Scripts/postprocessing.R:374-385`). |
| Nanom6A | `ratio.*.tsv` | no native numeric column is retained in standardized BEDs | parsed `mod_ratio` exists transiently but is not written | no standardized score column; PR uses filename-derived threshold, default metrics use binary calls from the `0.5` file | binary default recall/precision/F1 from `0.5`; PR curve from derived `max_thr` across multiple files | This repo does not use a Nanom6A output column for PR curves. It uses the threshold encoded in each filename as the derived ranking score (`Scripts/postprocessing.R:410-426`; `Scripts/statistical_analysis.R:281-315`). |
| Tombo (sample comparison) | `sample.level_samp_comp_detect.statistic.plus.wig` | parsed p-value from the second parsed field in the row-wise reader | repo-parsed as column `2` in the WIG table loop | standardized column `6` = `Pvalue` | PR curves, manual threshold sweep, default-threshold metrics, score histogram, default-threshold plotting filter | The repo parses WIG rows manually, so “column 2” is the parser’s second field rather than a stable named output column (`Scripts/postprocessing.R:443-456,513-516`). |
| xPore | `diffmod.table` | native column 5 (p-value, then BH-adjusted) | `5` | standardized column `6` = `FDR` | PR curves, manual threshold sweep, default-threshold metrics, score histogram, default-threshold plotting filter | Downstream uses the BH-adjusted version written to standardized column 6, not the raw native p-value directly (`Scripts/postprocessing.R:342-357`). |
| Yanocomp | `yanocomp_output.bed` | `V9` (`Score`) | `9` | standardized column `6` = `Score` | PR curves, manual threshold sweep, default-threshold metrics, score histogram, default-threshold plotting filter | Native `V6` is strand; native `V9` is the downstream score (`Scripts/postprocessing.R:145-151`). |

## Standardized Output Schema and Column Usage

The standardized outputs written to `output_bed_files/` are BED-like tables with a common leading schema:

- `Chr`
- `Start`
- `End`
- `Strand`
- `Status`
- tool-specific score column in column 6 when the tool remains score-bearing after postprocessing (`Scripts/postprocessing.R:25-605`).

The common schema is not perfectly uniform across callers:

- Score-bearing standardized outputs: `DENA`, `DiffErr`, `DRUMMER`, `ELIGOS`, `EpiNano-Error`, `EpiNano-SVM`, `m6Anet`, `Nanocompore`, `NanoDoc`, `Tombo`, `xPore`, `Yanocomp` (`Scripts/postprocessing.R:25-195,277-395,437-598`).
- Binary-only standardized outputs: `MINES`, `Nanom6A` (`Scripts/postprocessing.R:203-269,396-436`).

Downstream usage of standardized columns is split across two scripts:

- `Scripts/statistical_analysis.R` reads each BED with `header = TRUE`, converts the first four standardized columns into `GRanges`, computes overlaps against genome bins, and, for scored tools, reads column 6 as the operative confidence parameter (`Scripts/statistical_analysis.R:168-177`).
- `Scripts/plottingScript.R` uses columns 1-4 to compute peak widths, metagene placement, gene overlaps, and filtered subsets; when a file has more than five columns, it treats column 6 as the score to histogram and threshold by the default threshold table (`Scripts/plottingScript.R:145-173,180-207`).

Three distinct notions of “used columns” exist in this repository:

1. Native columns used only to decide whether a row survives postprocessing.
   Example: `Nanocompore` uses `Logit_LOR` only for the filter `abs(Logit_LOR) > 0.5` and then discards it (`Scripts/postprocessing.R:165-172`).
2. Native columns used to populate the standardized BED-like file.
   Example: `DiffErr` converts native `V5` into standardized `FDR` and carries that forward as column 6 (`Scripts/postprocessing.R:27-33`).
3. Standardized columns actually consumed downstream.
   Example: in `statistical_analysis.R`, only standardized column 6 participates in score aggregation and threshold sweeps for scored tools; `Status` does not drive those analyses (`Scripts/statistical_analysis.R:173-234`).

Actual downstream column usage is:

- `Chr`, `Start`, `End`, `Strand`: used to build `GRanges` and compute overlaps with gene bins, RRACH motifs, reference peaks, genes, and metagene coordinates (`Scripts/statistical_analysis.R:168-172`; `Scripts/plottingScript.R:126-207,420-521`).
- Standardized column 6 when present: used as the score/confidence parameter in PR curve generation, default-threshold binarization, histogram plotting, and default-threshold filtering for metagene and gene-level summaries (`Scripts/statistical_analysis.R:176-234`; `Scripts/plottingScript.R:154-173`).
- `Status`: primarily records postprocessing retention and is not used as the analytical decision variable in the core benchmarking code because the score-bearing tools are re-thresholded downstream from column 6 (`Scripts/statistical_analysis.R:176-234`; `Scripts/plottingScript.R:154-173`).

Special downstream handling:

- `MINES` is binary-only because postprocessing drops its native score before writing the standardized output; the statistical script therefore routes it to the `MINES` branch that computes only overlap-based recall, precision, and F1 (`Scripts/postprocessing.R:265-269`; `Scripts/statistical_analysis.R:239-255`).
- `Nanom6A` is also treated as binary in the standardized BEDs because postprocessing discards the parsed `mod_ratio`; the statistical script reconstructs a derived ranking (`max_thr`) across the multiple threshold-specific files for PR analysis, while default metrics use only the `0.5` output (`Scripts/postprocessing.R:412-426`; `Scripts/statistical_analysis.R:150-156,257-315`; `Scripts/plottingScript.R:62-64`).

## Downstream Analyses

### `downstream_analysis_map`

| analysis_name | script | input_object | columns_used | comparison | statistic_computed | formula_or_rule | output_artifact | importance_for_tool_evaluation | formal_test_present |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Gene-window binning | `Scripts/statistical_analysis.R` | `genesbed` -> `genesBins` | gene BED columns `chr`, `start`, `end`, `name`, `strand` | each gene split into fixed windows of length `binLength` | fixed-width analysis bins | per gene, split coordinates into windows and discard bins whose width is not exactly `w` (`Scripts/statistical_analysis.R:42-88`) | in-memory `genesBins`; saved indirectly in `Results_window_*bp.rda` | foundational; defines the site-to-bin benchmark resolution | no |
| Reference-set labeling | `Scripts/statistical_analysis.R` | `peaks_bed`, `genesBins_par` | peak `chr/start/end/strand`; bin coordinates | gold-standard peak bins vs all bins | binary reference label | `hitsMatrix[subjectHits(peaks_overlap), "Reference_set"] <- 1` (`Scripts/statistical_analysis.R:147-148`) | `hitsMatrix`, `overlapMatrix`, `Results_window_*bp.rda` | foundational; all recall/precision/F1 and PR analyses compare tool predictions against these labels | no |
| RRACH subset restriction | `Scripts/statistical_analysis.R` | genome FASTA motifs, `genesBins`, reference peaks | motif coordinates, bin coordinates, peak overlaps | RRACH-bearing bins and RRACH-bearing peak bins vs full set | subset selection | identify `RRACH` / `DGTYY` motif positions, keep bins overlapping motifs, and keep reference-positive motif bins (`Scripts/statistical_analysis.R:34-40,90-101`) | `Results_window_*bp_RRACH.rda`; RRACH performance TSV/PDF suffixes | important for motif-bias analysis, not a general headline metric | no |
| Optional high-coverage restriction | `Scripts/statistical_analysis.R` | optional `highcov_bed_file` | high-coverage BED coordinates | high-coverage bins vs full bins | subset selection | overlaps high-coverage regions with bins, then with reference peaks and RRACH motifs (`Scripts/statistical_analysis.R:103-123,377-390`) | `Results_window_*bp_highcov.rda` and `_highcov_RRACH.rda` when input exists | important for robustness under coverage filtering, but not wired by default pipeline invocation | no |
| Per-bin score aggregation | `Scripts/statistical_analysis.R` | standardized BEDs, `findOverlaps` output | columns 1-4 and standardized column 6 | multiple calls from one tool within one bin | aggregated bin score | maximize column 6 for `DENA`, `EpiNano-Error`, `EpiNano-SVM`, `NanoDoc`, `m6Anet`; minimize then sign-flip for `DiffErr`, `DRUMMER`, `Yanocomp`, `Nanocompore`, `ELIGOS`, `xPore`, `Tombo` (`Scripts/statistical_analysis.R:27-32,173-195`) | `hitsMatrix` numeric scores | critical for translating site-level calls into comparable bin-level rankings | no |
| Precision-recall curve generation | `Scripts/statistical_analysis.R` | `positive`, `negative` bin-score vectors | aggregated bin scores and reference labels | score distribution in reference-positive vs reference-negative bins | PR curve and PR AUC via PRROC | `pr.curve(scores.class0 = positive, scores.class1 = negative, curve = TRUE, rand.compute = TRUE)` (`Scripts/statistical_analysis.R:202-215`) | per-tool `*_PRcurve*.Rdata` and `*.pdf`; summary PR PDF | primary ranking metric when a tool has a usable score column | no |
| Manual threshold sweep | `Scripts/statistical_analysis.R` | aggregated bin scores | bin scores and reference labels | thresholds across one tool’s score range | recall, precision, F1 | threshold grid = 100 evenly spaced values between min/max score plus default threshold; `recall = TP/(TP+FN)`, `precision = TP/(TP+FP)`, `F1 = 2RP/(R+P)` (`Scripts/statistical_analysis.R:219-234`) | `Performances*_window_*bp.tsv` | primary default-threshold and threshold-trace summary | no |
| Default-threshold binarization | `Scripts/statistical_analysis.R` | `hitsMatrix` and `threshold_default` | aggregated bin score column per tool | predicted-positive bins at repo default threshold vs reference-positive bins | default recall, precision, F1 and binary overlap state | `pred_pos_def <- which(hitsMatrix[,x] >= default_thr)` after maximizing or sign-flipping (`Scripts/statistical_analysis.R:176-200`) | `overlapMatrix`; last row or appended default entry in `Performances` | crucial for “out-of-the-box” tool behavior | no |
| Binary-only MINES evaluation | `Scripts/statistical_analysis.R` | MINES standardized BED | columns 1-4 only | binary MINES-positive bins vs reference-positive bins | recall, precision, F1 | counts overlaps only; no score sweep (`Scripts/statistical_analysis.R:239-255`) | one-row `Performances[[x]]` entry for MINES | necessary because MINES has no retained score column | no |
| Multi-threshold Nanom6A ranking | `Scripts/statistical_analysis.R` | multiple Nanom6A BEDs from different `ratio.*.tsv` files | columns 1-4 from each BED; threshold encoded in filename | highest threshold retaining a positive bin vs reference labels | derived score `max_thr`, PR curve, default `0.5` recall/precision/F1 | per bin, choose the highest Nanom6A threshold file in which the bin is positive; use that numeric threshold as ranking score (`Scripts/statistical_analysis.R:150-156,257-315`) | `Nanom6A_PRcurve*.Rdata`, `Nanom6A_PRcurve*.pdf`, `Performances` entry for `0.5` | essential for making Nanom6A comparable to scored tools despite binary BEDs | no |
| Tool-overlap heatmap | `Scripts/statistical_analysis.R` | `overlapMatrix` | default-threshold binary call matrix | tool `i` positives vs tool `j` positives | asymmetric overlap fraction | `| positives_i ∩ positives_j | / | positives_i |` (`Scripts/statistical_analysis.R:346-358`) | `Tools_overlap_default_par*_window_*bp.pdf` | useful for agreement structure and redundancy, but not for correctness on its own | no |
| Gene-level secondary evaluation | `Scripts/plottingScript.R` | standardized BEDs, gene GRanges | columns 1-4 and optionally column 6 for filtering | tool-positive genes vs reference-positive genes | gene-level recall, precision, F1 | overlap filtered tool peaks with genes, then compare gene sets with reference-set genes (`Scripts/plottingScript.R:126-207`) | console summaries; informs metagene preparation | secondary; less specific than site/bin-level benchmarking | no |
| RRACH motif-stratified recall | `Scripts/plottingScript.R` | saved `Results_window_50bp.rda`, motif-split bins | `hitsMatrix` plus RRACH motif assignment | reference-positive bins stratified by RRACH motif vs tool positives | per-motif recall, plotted as “accuracy” | `recall = TP/(TP+FN)` within bins assigned to each motif (`Scripts/plottingScript.R:420-540`) | `RRACH_accuracy.Rdata`, `RRACH_accuracy.pdf`, `RRACHsAccuracies.pdf` | useful for motif bias profiling; not a general accuracy metric | no |
| TP/FN/FP sequence-context diagnostics | `Scripts/plottingScript.R` | `hitsMatrix`, extracted sequences, RNAfold outputs | TP/FN/FP bin sequences | TP vs FN vs FP vs control sequences for each tool | GC content, free energy, Shannon entropy distributions | GC fraction from FASTA, RNAfold free energy, Shannon entropy; plotted as boxplots (`Scripts/plottingScript.R:571-813`) | `FP_bins.RData`, `TP_bins.RData`, `FN_bins.RData`, `sequenceFeatures.pdf` | useful for diagnosing sequence/context bias, not for headline ranking | no |
| Runtime summaries | `Scripts/plottingScript.R` | hard-coded CPU-hour vectors from Nextflow Tower | per-tool runtime totals split by preprocessing/tool execution | tools within each dataset | CPU hours | stacked sums of hard-coded preprocessing and execution times (`Scripts/plottingScript.R:1228-1300`) | `Running_time_*.pdf` | operational comparison only | no |
| Downsampling robustness | `Scripts/plottingScript.R` | hard-coded HEK293T metric vectors | tool, fraction of dataset, metric value | 25/50/75/100% subsets | F1, PR AUC, hit count, normalized hit count | plot metric vs fraction; normalized hits divide each subset’s hit count by the full-data hit count (`Scripts/plottingScript.R:1303-1419`) | `HEK293T_downsampling_*.pdf` | operational robustness after accuracy is understood | no |

## Visualization Audit

### `visualization_audit`

| figure_name | script | input_data | x_axis | y_axis | grouping | statistic_shown | what_is_compared | formal_statistical_test | interpretation_value | caveats |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Per-tool peak-size histograms | `Scripts/plottingScript.R` | standardized BEDs after optional default filtering | peak size in nucleotides | count of peaks | one tool per panel | empirical distribution of `End - Start` | size distribution of called sites within a tool | none | descriptive check for whether a tool emits narrow sites or broad regions | Depends on whether the BED has already been filtered by default thresholds and on `proper_bed_flag = 0` (`Scripts/plottingScript.R:136-179`). |
| Per-tool score distributions with default threshold line | `Scripts/plottingScript.R` | standardized BEDs with column 6 present | tool-specific score/FDR/P-value/probability | count of peaks | one tool per panel | histogram of standardized column 6 plus default threshold line | score distribution relative to the repo’s default cutoff | none | important for understanding threshold placement and hit inflation | Uses hard-coded default thresholds from `plottingScript.R`, which are inconsistent for DRUMMER (`Scripts/plottingScript.R:154-173`). |
| Reference-set metagene plot | `Scripts/plottingScript.R` | reference-set BED written to `/tmp` and transcript database | normalized metagene position | density from Guitar | reference set only | metagene density | placement of gold-standard sites across transcript structure | none | biological baseline for comparing caller positional bias | Dataset-specific reference files are external inputs, not generated in this repo (`Scripts/plottingScript.R:95-123`). |
| Per-tool metagene plots | `Scripts/plottingScript.R` | filtered per-tool BEDs expanded by `tol` | normalized metagene position | density from Guitar | one tool per plot | metagene density | positional distribution of one tool’s calls | none | useful for biological plausibility and transcript-region bias | Downsamples large BEDs to `max_filt_sites` and expands intervals by `tol = 100`, which changes apparent positional support (`Scripts/plottingScript.R:183-237`). |
| Combined metagene plot across tools | `Scripts/plottingScript.R` | `p3[[1]]$data` from Guitar for all tools | metagene coordinate `x` | density | tool | density curves | positional distributions of all tools on the same axes | none | useful side-by-side view of positional bias | Purely descriptive; no statistical comparison or CI is shown (`Scripts/plottingScript.R:240-260`). |
| Tool hit-count barplots by dataset | `Scripts/plottingScript.R` | hard-coded `num_filt_peaks_*` vectors | tool | number of hits on log10 scale | dataset | number of retained hits per tool | relative call volume across tools within each dataset | none | useful hit-inflation sanity check | Counts are hard-coded manuscript/reporting values rather than recomputed from current outputs (`Scripts/plottingScript.R:262-389`). |
| RRACH motif-stratified “accuracy” barplot | `Scripts/plottingScript.R` | `hitsMatrix` and motif-split RRACH-positive bins | tool | value labeled `accuracy` | RRACH motif variant | actually recall within motif-stratified reference-positive bins | how well each tool recovers positives for each RRACH variant | none | useful for motif bias profiling | The plotted quantity is recall, not standard accuracy (`Scripts/plottingScript.R:498-540`). |
| RRACH common-vs-uncommon boxplots | `Scripts/plottingScript.R` | saved `RRACH_accuracy.Rdata` | tool groups | motif-stratified recall values | common vs uncommon RRACH groups | distribution of motif-specific recall values | motif frequency classes within each tool | none | descriptive summary of whether tools do better on common motifs | No inferential test is applied even though the plot suggests a two-group comparison (`Scripts/plottingScript.R:542-566`). |
| TP/FN/FP sequence-feature boxplots | `Scripts/plottingScript.R` | TP/FN/FP FASTAs, RNAfold outputs, entropy function | tool blocks | GC %, free energy, entropy | TP vs FN vs FP vs control | distribution summaries as boxplots | sequence-context properties of true positives, false negatives, and false positives | none | valuable for diagnosing systematic biases in which bins each tool misses or overcalls | Motif enrichment itself is offloaded externally to XSTREME online and is not reproducibly implemented in the repo (`Scripts/plottingScript.R:571-813`). |
| Default-parameter metric barplots | `Scripts/plottingScript.R` | hard-coded recall, precision, F1 arrays by dataset | tool | metric value | metric (`Recall`, `Precision`, `F1 score`) | default-threshold recall, precision, F1 | default out-of-the-box behavior of tools | none | one of the most important summary plots for practical use | Values are hard-coded from prior outputs rather than recomputed live; dataset-specific reference sets differ (`Scripts/plottingScript.R:818-1071`). |
| Summary PR curve plots | `Scripts/statistical_analysis.R` and `Scripts/plottingScript.R` | saved PRROC objects and hard-coded default recall/precision points | recall | precision | tool plus random baseline | PR curves, implicit PR AUC, default operating-point markers | score-based ranking performance of tools | none | primary comparative visualization for scored tools | AUC values are not exported to TSV by `statistical_analysis.R`; later downsampling plots use hard-coded AUC values read from prior PR plots (`Scripts/statistical_analysis.R:209-215,321-335`; `Scripts/plottingScript.R:1072-1225,1331-1355`). |
| Tool-overlap heatmap at default thresholds | `Scripts/statistical_analysis.R` | `overlapMatrix` at default thresholds | tool `j` | tool `i` | tool pair | asymmetric overlap fraction | how much of tool `i`’s default-positive set is shared with tool `j` | none | useful for agreement structure and redundancy | This is not Jaccard or symmetric concordance; interpretation must be directional (`Scripts/statistical_analysis.R:346-358`). |
| Running-time stacked barplots | `Scripts/plottingScript.R` | hard-coded CPU-hour vectors from Nextflow Tower | tool | CPU hours | preprocessing vs tool execution | stacked runtime totals | operational cost of running each caller | none | important for deployment trade-offs after accuracy is understood | Times are manually entered, dataset-specific, and not reproducibly recomputed from current Nextflow runs (`Scripts/plottingScript.R:1228-1300`). |
| HEK293T downsampling plots | `Scripts/plottingScript.R` | hard-coded F1, AUC, and hit-count vectors for 25/50/75/100% subsets | fraction of dataset | F1, AUC, number of hits, normalized hits | tool | trend lines over downsampling fractions | robustness of each tool to reduced data volume | none | useful operational robustness diagnostic | Only HEK293T is covered; all values are hard-coded; normalized hit plots divide by each tool’s full-data hit count (`Scripts/plottingScript.R:1303-1419`). |

## Which Statistics Matter Most

1. **PR AUC from PR curves is the strongest ranking statistic when a tool exposes a usable score.**  
   In this repository, scored tools are reduced to per-bin numeric scores and passed to `PRROC::pr.curve`, which evaluates how well score rank separates reference-positive bins from reference-negative bins (`Scripts/statistical_analysis.R:173-215`). This matters more than a single cutoff because it evaluates the entire ranking behavior, not just one operating point. It is the best headline statistic for `DENA`, `DiffErr`, `DRUMMER`, `ELIGOS`, `EpiNano-Error`, `EpiNano-SVM`, `m6Anet`, `Nanocompore`, `NanoDoc`, `Tombo`, `xPore`, and `Yanocomp`. For `Nanom6A`, the PR ranking is reconstructed from `max_thr`, so its PR curve is derived rather than native (`Scripts/statistical_analysis.R:281-315`).

2. **Recall, precision, and F1 at the default threshold are the most important practical deployment statistics.**  
   The benchmark explicitly computes these from the default threshold table after aggregating calls to bins (`Scripts/statistical_analysis.R:176-234`). These values answer the operational question, “What happens if I run the tool and use the repository’s default evaluation cutoff?” The grouped barplots in `plottingScript.R` are therefore a strong practical summary even though they are manuscript-style hard-coded reproductions (`Scripts/plottingScript.R:818-1071`).

3. **Hit-count inflation is an important behavioral diagnostic, but not a correctness metric by itself.**  
   The repository repeatedly visualizes number of calls because tools differ by orders of magnitude in call volume (`Scripts/plottingScript.R:262-389,1400-1419`). A tool with extremely high hit count may have high recall only because it overcalls; a tool with very low hit count may have high precision but poor recall. Call volume should always be interpreted together with PR AUC or default-threshold recall/precision/F1.

4. **The tool-overlap heatmap is useful for agreement structure and redundancy, not for correctness.**  
   The overlap statistic is directional and measures how much of tool `i`’s positive set is contained in tool `j`’s positive set (`Scripts/statistical_analysis.R:346-358`). High overlap can mean shared biological signal or shared bias. It does not tell whether either tool is correct with respect to the reference set.

5. **Gene-level recall/precision/F1 in `plottingScript.R` are secondary summaries.**  
   These metrics collapse site-level calls to genes by overlap with gene coordinates (`Scripts/plottingScript.R:126-207`). They are easier to interpret biologically but lose resolution relative to the bin-level benchmark. They are useful if the downstream use case is gene prioritization rather than exact-site detection.

6. **Metagene distributions are biological plausibility metrics, not primary accuracy metrics.**  
   The Guitar-based plots test whether a tool’s calls accumulate in expected transcript regions (`Scripts/plottingScript.R:116-260`). They are valuable for spotting gross positional bias or unexpected enrichment patterns, but a visually plausible metagene profile does not imply high precision or recall.

7. **RRACH motif-stratified recall is useful for bias profiling, not overall ranking.**  
   The repository’s RRACH figure measures recall within reference-positive motif subsets, even though the y-axis is labeled “Accuracy” (`Scripts/plottingScript.R:498-540`). This can reveal that a tool works better on common RRACH variants than on rare ones, but it should not replace PR AUC or default-threshold recall/precision/F1 as the main evaluation statistic.

8. **TP/FN/FP sequence-context features are diagnostic, not headline metrics.**  
   GC content, RNAfold free energy, and Shannon entropy highlight which sequence contexts generate false positives or false negatives (`Scripts/plottingScript.R:571-813`). These are valuable when choosing or improving a caller, but they explain errors rather than summarize overall performance.

9. **Running time and downsampling robustness matter after accuracy is characterized.**  
   Runtime and robustness plots show whether a tool is operationally viable and whether it degrades gracefully at lower depth (`Scripts/plottingScript.R:1228-1419`). They matter in production settings, but they do not replace accuracy metrics.

## Cross-Repo Inconsistencies and Reproducibility Caveats

- `plottingScript.R` is partly manual and not fully pipeline-reproducible. Large sections hard-code counts, recall/precision/F1 values, AUC values, and runtime totals instead of recomputing them from the current run outputs (`Scripts/plottingScript.R:262-389,818-1419`).
- No formal inferential hypothesis tests are implemented for tool-to-tool comparisons. The repository computes descriptive overlaps, PR curves, and distribution plots only (`Scripts/statistical_analysis.R:125-390`; `Scripts/plottingScript.R:117-258,528-1419`).
- The RRACH “accuracy” plot is actually motif-stratified recall over reference-positive bins: the code calculates `TP / (TP + FN)` and never computes `TN`-based accuracy (`Scripts/plottingScript.R:504-518`).
- `DRUMMER` default thresholds are inconsistent across scripts: `Scripts/postprocessing.R` and `Scripts/statistical_analysis.R` use `0.05`, while `Scripts/plottingScript.R` hard-codes `0.01` in its default-threshold table (`Scripts/postprocessing.R:619-620`; `Scripts/statistical_analysis.R:30-32`; `Scripts/plottingScript.R:71-74`).
- High-coverage analyses exist in `Scripts/statistical_analysis.R`, but the default Nextflow `postprocessing` process does not pass a `highcov_bed_file`, so those branches are dormant unless the script is run separately with additional inputs (`Scripts/statistical_analysis.R:103-123,377-390`; `pipeline.nf:1169-1170`).
- `NanoDoc` is disabled by default in `pipeline.conf`, even though both postprocessing and plotting code include it in the benchmark design (`pipeline.conf:88-90`; `Scripts/postprocessing.R:365-395`; `Scripts/plottingScript.R:67-74,244-246`).
- Several figures rely on external artifacts rather than on pipeline-generated outputs alone:
  - dataset-specific gold-standard BEDs (`Scripts/plottingScript.R:39-57`);
  - Nextflow Tower runtime totals (`Scripts/plottingScript.R:1230,1248,1266,1284`);
  - online MEME-suite/XSTREME motif enrichment outside the repo (`Scripts/plottingScript.R:603-604`).
- The benchmark saves PR curve objects and plots, but `Scripts/statistical_analysis.R` does not export numeric PR AUC values to TSV. Later AUC summaries in `plottingScript.R` are hard-coded from prior figure outputs, not recomputed from the saved PR objects (`Scripts/statistical_analysis.R:209-215,301-305,317-335`; `Scripts/plottingScript.R:1331-1355`).
- `MINES` loses its numeric native score during postprocessing and is forced into binary-only downstream evaluation (`Scripts/postprocessing.R:265-269`; `Scripts/statistical_analysis.R:239-255`).
- `Nanom6A` also loses its native parsed `mod_ratio` in the standardized BEDs, so its default metrics are binary and its PR ranking is reconstructed indirectly from filename thresholds rather than from a native numeric output column (`Scripts/postprocessing.R:412-426`; `Scripts/statistical_analysis.R:281-315`).
- `ELIGOS` computes a `rows_eligos` completeness mask but then immediately reassigns `eligos <- data_eligos[, c(1:4,18,16,18)]`, so the intended NA-row filtering is not actually preserved in the downstream object (`Scripts/postprocessing.R:183-186`).
- The `MINES` branch in `Scripts/statistical_analysis.R` writes `threshold = thresholds` even though `thresholds` is only defined in the scored-tool branch, so the threshold column for MINES is inherited from prior loop state rather than being explicitly defined for MINES (`Scripts/statistical_analysis.R:219-255`).
- The summary PR curve plot is guarded by `if (length(negative) > 0)` after the per-file loop, which means the control variable comes from the loop’s last score-bearing tool rather than being evaluated independently for each tool or for the PR-curve list as a whole (`Scripts/statistical_analysis.R:321-337`).

## Appendix: Support/Preprocessing Tools

The repository uses the following support stack to create caller inputs, standardize outputs, and generate evaluation figures:

| component | role in benchmark | repo evidence |
| --- | --- | --- |
| Nextflow | orchestration of the full benchmark workflow | `pipeline.nf`, `README.md`, `pipeline.conf` (`README.md:1-43`; `pipeline.nf:1-1177`) |
| Singularity | container runtime for the Nextflow workflow | `pipeline.conf` singularity block (`pipeline.conf:102-106`) |
| `ont-fast5-api` / `multi_to_single_fast5` | convert multi-read FAST5 into single-read FAST5 directories | `pipeline.nf` `multi2single`; `Docker/ont_fast5_api-4.0.0/Dockerfile` (`pipeline.nf:102-138`; `Docker/ont_fast5_api-4.0.0/Dockerfile:1-17`) |
| Poretools | extract FASTQ from single-read FAST5 | `pipeline.nf` `fastq`; `Docker/poretools-0.6.0/Dockerfile` (`pipeline.nf:140-169`; `Docker/poretools-0.6.0/Dockerfile:1-8`) |
| minimap2 + samtools | transcriptome and genome alignments and BAM merging | `pipeline.nf` `minimap2`/`minimap2Merge`; `Docker/minimap2-2.24.0/Dockerfile` (`pipeline.nf:171-300`; `Docker/minimap2-2.24.0/Dockerfile:1-12`) |
| Tombo resquiggle and Tombo de novo | resquiggling plus de novo statistics used upstream of MINES and for Tombo comparison | `pipeline.nf` `tombo1`, `tombo2`, `tombo3`; `Docker/tombo-1.5.1/Dockerfile` (`pipeline.nf:302-398`; `Docker/tombo-1.5.1/Dockerfile:1-17`) |
| Nanopolish | transcriptome eventalign generation for `xPore`, `Nanocompore`, `Yanocomp`, and `m6Anet` | `pipeline.nf` `nanopolish1`; `Docker/nanopolish-0.8.4/Dockerfile` (`pipeline.nf:860-898`; `Docker/nanopolish-0.8.4/Dockerfile:1-17`) |
| Picard | sequence-dictionary creation for Nanom6A and EpiNano workflows | `pipeline.nf` `nanom6a`; `Docker/nanom6a-22_10_2021/Dockerfile` (`pipeline.nf:420-421`; `Docker/nanom6a-22_10_2021/Dockerfile:1-21`) |
| `wig2bed` | convert Tombo WIG outputs for MINES input | `pipeline.nf` `mines` (`pipeline.nf:553-557`) |
| `Scripts/postprocessing.R` | normalize heterogeneous caller outputs into standardized BED-like files | `pipeline.nf` postprocessing invocation; script itself (`pipeline.nf:1142-1170`; `Scripts/postprocessing.R:18-646`) |
| `Scripts/statistical_analysis.R` | bin-level benchmarking, PR curves, overlap matrices, RRACH/high-coverage subsets | `pipeline.nf` postprocessing invocation; script itself (`pipeline.nf:1169-1170`; `Scripts/statistical_analysis.R:1-390`) |
| `Guitar`, `PRROC`, `pheatmap`, `ggplot2` | metagene plots, PR curves, overlap heatmaps, barplots | `Scripts/statistical_analysis.R` and `Scripts/plottingScript.R` imports and usage (`Scripts/statistical_analysis.R:14-22,209-215,357-358`; `Scripts/plottingScript.R:5-17,123-258`) |
| ViennaRNA `RNAfold` | sequence free-energy calculation for TP/FN/FP diagnostics | `Scripts/plottingScript.R:656-759` |
| MEME-suite / XSTREME | external motif-enrichment step for exported FASTAs | `Scripts/plottingScript.R:603-604` |
| Nextflow Tower | source of manually reported CPU-hour runtime values in plotting script | `pipeline.conf:108-112`; `Scripts/plottingScript.R:1230-1299` |

### Exact Caller Commands Used in the Pipeline

Below are the exact caller invocations used in `pipeline.nf`, stripped down to the commands that actually produce caller outputs.

#### DENA

```bash
python3 /DENA/step4_predict/LSTM_extract.py get_pos --fasta transcriptome.fa --motif 'RRACH' --output ${params.resultsDir}/${condition1}/dena/candidate_predict_pos.txt
python3 /DENA/step4_predict/LSTM_extract.py predict --fast5 ${params.resultsDir}/${condition1} --corr_grp ${params.tombo_slot} --bam minimap.filt.sort.1.bam --sites ${params.resultsDir}/${condition1}/dena/candidate_predict_pos.txt --label "dena_label" --windows 2 2 --processes ${task.cpus}
python3 /DENA/step4_predict/LSTM_predict.py -i ${params.resultsDir}/${condition1}/dena/ -m /DENA/denaModels/ -o ${params.resultsDir}/${condition1}/dena/ -p "dena_label"
```

Source: [pipeline.nf](/Users/bmorampa/NanOlympicsMod/pipeline.nf#L584C1).

#### DiffErr

```bash
differr -p ${task.cpus} \
$(for file in minimap.sortG.2*.bam; do echo -a $file; done) \
$(for file in minimap.sortG.1*.bam; do echo -b $file; done) \
-r genome.fa -o ${params.resultsDir}/differr/differrOut.bed \
-f ${params.differrFDR}
```

Source: [pipeline.nf](/Users/bmorampa/NanOlympicsMod/pipeline.nf#L467C1).

#### DRUMMER

```bash
while read -r line;
do
  python3 DRUMMER.py -r ${params.genome_fasta} \
  -t ${params.resultsDir}/drummer/minimap.sortG.2.*.bam \
  -c ${params.resultsDir}/drummer/minimap.sortG.1.*.bam \
  -o ${params.resultsDir}/drummer/DRUMMER/$line/ \
  -a exome \
  -p ${params.drummerPval} \
  -n $line \
  -m true ;
done < ${params.resultsDir}/drummer/chromosomes.txt || true
```

Source: [pipeline.nf](/Users/bmorampa/NanOlympicsMod/pipeline.nf#L842C1).

#### ELIGOS

```bash
/eligos2-v2.1.0/eligos2 pair_diff_mod -t ${task.cpus} \
-tbam minimap.sortG.1.bam \
-cbam minimap.sortG.2.bam \
-reg genome.bed \
-ref genome.fa \
-o ${params.resultsDir}/eligos/merged/
```

```bash
/eligos2-v2.1.0/eligos2 pair_diff_mod -t ${task.cpus} \
-tbam minimap.sortG.1.*.bam \
-cbam minimap.sortG.2.*.bam \
-reg genome.bed \
-ref genome.fa \
-o ${params.resultsDir}/eligos/
```

Only the merged output is postprocessed downstream. Source: [pipeline.nf](/Users/bmorampa/NanOlympicsMod/pipeline.nf#L504C1).

#### EpiNano-Error

```bash
/usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Variants.py -R genome.fa -b minimap.sort.1.bam -s /EpiNano-Epinano1.2.1/misc/sam2tsv.jar --type g -n ${task.cpus}
/usr/bin/python3 /EpiNano-Epinano1.2.1/misc/Epinano_sumErr.py --kmer 0 --file minimap.sort.1.plus*site.csv --out minimap.sort.1.plus.sumErrOut.csv
/usr/bin/python3 /EpiNano-Epinano1.2.1/misc/Epinano_sumErr.py --kmer 0 --file minimap.sort.1.minus*site.csv --out minimap.sort.1.minus.sumErrOut.csv
/usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Variants.py -R genome.fa -b minimap.sort.2.bam -s /EpiNano-Epinano1.2.1/misc/sam2tsv.jar --type g -n ${task.cpus}
/usr/bin/python3 /EpiNano-Epinano1.2.1/misc/Epinano_sumErr.py --kmer 0 --file minimap.sort.2.plus*site.csv --out minimap.sort.2.plus.sumErrOut.csv
/usr/bin/python3 /EpiNano-Epinano1.2.1/misc/Epinano_sumErr.py --kmer 0 --file minimap.sort.2.minus*site.csv --out minimap.sort.2.minus.sumErrOut.csv
/bin/miniconda3/bin/Rscript /EpiNano-Epinano1.2.1/Epinano_DiffErr.R -k minimap.sort.2.plus.sumErrOut.csv -w minimap.sort.1.plus.sumErrOut.csv -d ${params.epinanoErrorSumErr} -t 3 -p -o diffErr -f sum_err
/bin/miniconda3/bin/Rscript /EpiNano-Epinano1.2.1/Epinano_DiffErr.R -k minimap.sort.2.minus.sumErrOut.csv -w minimap.sort.1.minus.sumErrOut.csv -d ${params.epinanoErrorSumErr} -t 3 -p -o diffErr -f sum_err
```

Source: [pipeline.nf](/Users/bmorampa/NanOlympicsMod/pipeline.nf#L703C1).

#### EpiNano-SVM

```bash
/usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Variants.py -n ${task.cpus} -T g -R genome.fa -b minimap.sort.1.bam -s /EpiNano-Epinano1.2.1/misc/sam2tsv.jar
/usr/bin/python3 /EpiNano-Epinano1.2.1/misc/Slide_Variants.py minimap.sort.1.plus_strand.per.site.csv 5
/usr/bin/python3 /EpiNano-Epinano1.2.1/misc/Slide_Variants.py minimap.sort.1.minus_strand.per.site.csv 5
/usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Predict.py --model /EpiNano-Epinano1.2.1/models/rrach.q3.mis3.del3.linear.dump --predict minimap.sort.1.plus_strand.per.site.5mer.csv --columns 8,13,23 --out_prefix plus_mod_prediction
/usr/bin/python3 /EpiNano-Epinano1.2.1/Epinano_Predict.py --model /EpiNano-Epinano1.2.1/models/rrach.q3.mis3.del3.linear.dump --predict minimap.sort.1.minus_strand.per.site.5mer.csv --columns 8,13,23 --out_prefix minus_mod_prediction
```

Source: [pipeline.nf](/Users/bmorampa/NanOlympicsMod/pipeline.nf#L631C1).

#### m6Anet

```bash
m6anet-dataprep --eventalign ${params.resultsDir}/${condition}/${sample}/nanopolish/transcriptome/eventalign_readIndex.txt \
--out_dir ${params.resultsDir}/${condition}/${sample}/m6anet --n_processes ${task.cpus}

m6anet-run_inference --input_dir $preprocessing_dirs --out_dir ${params.resultsDir}/m6anet --infer_mod_rate --n_processes ${task.cpus}
zcat ${params.resultsDir}/m6anet/data.result.csv.gz > ${params.resultsDir}/m6anet/data.result.csv
```

Source: [pipeline.nf](/Users/bmorampa/NanOlympicsMod/pipeline.nf#L1038C1) and [pipeline.nf](/Users/bmorampa/NanOlympicsMod/pipeline.nf#L1068C1).

#### MINES

```bash
wig2bed < ${params.resultsDir}/${condition1}/tomboDenovo/output_filename.fraction_modified_reads.plus.wig > ${params.resultsDir}/${condition1}/mines/output_filename.fraction_modified_reads.plus.wig.bed
python3 /MINES/cDNA_MINES.py --fraction_modified ${params.resultsDir}/${condition1}/mines/output_filename.fraction_modified_reads.plus.wig.bed --coverage ${params.resultsDir}/${condition1}/tomboDenovo/output_filename.coverage.plus.bedgraph --output ${params.resultsDir}/${condition1}/mines/m6A_output_filename.bed --ref transcriptome.fa --kmer_models /MINES/Final_Models/names.txt
```

Source: [pipeline.nf](/Users/bmorampa/NanOlympicsMod/pipeline.nf#L553C1).

#### Nanocompore

```bash
nanocompore eventalign_collapse --eventalign ${params.resultsDir}/${condition}/${sample}/nanopolish/transcriptome/eventalign_readName.txt --nthreads ${task.cpus} --outpath ${params.resultsDir}/${condition}/${sample}/nanocompore/ --overwrite

nanocompore sampcomp --file_list1 "${f1[*]}" --file_list2 "${f2[*]}" \
--label1 ${condition1} \
--label2 ${condition2} \
--fasta transcriptome.fa \
--bed transcriptome.bed \
--outpath ${params.resultsDir}/nanocompore/ \
--allow_warnings \
--logit \
--nthreads ${task.cpus} \
--overwrite
```

Source: [pipeline.nf](/Users/bmorampa/NanOlympicsMod/pipeline.nf#L968C1) and [pipeline.nf](/Users/bmorampa/NanOlympicsMod/pipeline.nf#L1005C1).

#### NanoDoc

```bash
/bin/miniconda3/bin/python /nanoDoc/src/nanoDoc.py formatfile -i ${params.resultsDir}/${condition1}/ -o . -r genome.fa -t ${task.cpus}
/bin/miniconda3/bin/python /nanoDoc/src/nanoDoc.py formatfile -i ${params.resultsDir}/${condition2}/ -o . -r genome.fa -t ${task.cpus}
cat genome.bed | while read line; do chr=$(echo $line | cut -d' ' -f1); start=$(echo $line | cut -d' ' -f2); end=$(echo $line | cut -d' ' -f3); /bin/miniconda3/bin/python /nanoDoc/src/nanoDoc.py analysis -w /nanoDoc/weight5mer/ -p /nanoDoc/param20.txt -r genome.fa -rraw ${params.resultsDir}/nanodoc/${condition2}_output/ -traw ${params.resultsDir}/nanodoc/${condition1}_output/ -chrom $chr --start $start --end $end -o "nanoDoc_results_"$chr"_"$start"_"$end".txt"; done
```

Source: [pipeline.nf](/Users/bmorampa/NanOlympicsMod/pipeline.nf#L789C1).

#### Nanom6A

```bash
/nanom6A_2021_10_22/bin/extract_raw_and_feature_fast --cpu=${task.cpus} --fl=${params.resultsDir}/${condition1}/nanom6a/files.txt -o ${params.resultsDir}/${condition1}/nanom6a/result --clip=10 --basecall_group ${params.tombo_slot} --basecall_subgroup ${params.tombo_subslot}
for prob in ${params.nanom6AP}; do /nanom6A_2021_10_22/bin/predict_sites --cpu ${task.cpus} -i ${params.resultsDir}/${condition1}/nanom6a/result -o ${params.resultsDir}/${condition1}/nanom6a/result_final -r transcriptome.fa -g genome.fa -b ${params.genes2transcripts} --model /nanom6A_2021_10_22/bin/model/ --proba $prob; done
```

Source: [pipeline.nf](/Users/bmorampa/NanOlympicsMod/pipeline.nf#L425C1).

#### Tombo (sample comparison)

```bash
/bin/miniconda3/bin/tombo detect_modifications level_sample_compare \
--fast5-basedirs ${params.resultsDir}/${condition1}/ \
--alternate-fast5-basedirs ${params.resultsDir}/${condition2}/ \
--minimum-test-reads 50 \
--processes ${task.cpus} --statistics-file-basename ${params.resultsDir}/tomboComparison/sample.level_samp_comp_detect \
--store-p-value

/bin/miniconda3/bin/tombo text_output browser_files --statistics-filename ${params.resultsDir}/tomboComparison/sample.level_samp_comp_detect.tombo.stats \
--browser-file-basename ${params.resultsDir}/tomboComparison/sample.level_samp_comp_detect --file-types statistic
```

Source: [pipeline.nf](/Users/bmorampa/NanOlympicsMod/pipeline.nf#L381C1).

#### xPore

```bash
xpore dataprep --eventalign ${params.resultsDir}/${condition}/${sample}/nanopolish/transcriptome/eventalign_readIndex.txt --out_dir ${params.resultsDir}/${condition}/${sample}/xpore --gtf_or_gff genome.gtf --transcript_fasta transcriptome.fa --genome

xpore diffmod --config ${params.resultsDir}/xpore/xpore.yaml --n_processes ${task.cpus}
xpore postprocessing --diffmod_dir ${params.resultsDir}/xpore/
```

Source: [pipeline.nf](/Users/bmorampa/NanOlympicsMod/pipeline.nf#L913C1) and [pipeline.nf](/Users/bmorampa/NanOlympicsMod/pipeline.nf#L939C1).

#### Yanocomp

```bash
/bin/miniconda3/envs/yanocomp/bin/yanocomp prep -p ${task.cpus} -e ${params.resultsDir}/${condition}/${sample}/nanopolish/transcriptome/eventalign_readName.txt -h ${params.resultsDir}/${condition}/${sample}/yanocomp/transcriptome/output.hdf5
/bin/miniconda3/envs/yanocomp/bin/yanocomp prep -p ${task.cpus} -e ${params.resultsDir}/${condition}/${sample}/nanopolish/transcriptome/eventalign_readName.txt -h ${params.resultsDir}/${condition}/${sample}/yanocomp/genome/output.hdf5 -g genome.gtf

/bin/miniconda3/envs/yanocomp/bin/yanocomp gmmtest \
$(for file in outputG.1.*.hdf5; do echo -c $file; done) \
$(for file in outputG.2.*.hdf5; do echo -t $file; done) \
-p ${task.cpus} \
-o ${params.resultsDir}/yanocomp/yanocomp_output.bed \
-s ${params.resultsDir}/yanocomp/yanocomp_output.json.gzip \
-f ${params.yanocompFDR}
```

Source: [pipeline.nf](/Users/bmorampa/NanOlympicsMod/pipeline.nf#L1091C1) and [pipeline.nf](/Users/bmorampa/NanOlympicsMod/pipeline.nf#L1126C1).
