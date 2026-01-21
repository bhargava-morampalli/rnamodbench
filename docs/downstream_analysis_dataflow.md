# Downstream Analysis - Data Flow

## High-Level Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              INPUTS                                          │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐          │
│  │  TOMBO   │ │ YANOCOMP │ │NANOCOMPORE│ │  XPORE   │ │  ELIGOS  │          │
│  │  .csv    │ │  .bed    │ │   .tsv   │ │  .table  │ │   .txt   │          │
│  └────┬─────┘ └────┬─────┘ └────┬─────┘ └────┬─────┘ └────┬─────┘          │
│       │            │            │            │            │                  │
│  ┌────┴─────┐ ┌────┴─────┐ ┌────┴─────┐ ┌────┴─────┐                        │
│  │ EPINANO  │ │ DIFFERR  │ │ DRUMMER  │ │ JACUSA2  │                        │
│  │  .csv    │ │  .bed    │ │  .txt    │ │  .bed    │                        │
│  └────┬─────┘ └────┬─────┘ └────┬─────┘ └────┬─────┘                        │
│       │            │            │            │                               │
│       └────────────┴─────┬──────┴────────────┘                               │
│                          │                                                   │
│                          ▼                                                   │
│                 ┌─────────────────┐      ┌─────────────────┐                │
│                 │  GROUND TRUTH   │      │   PARAMETERS    │                │
│                 │  (known_sites)  │      │  - threshold    │                │
│                 │     .csv        │      │  - min_reps     │                │
│                 └────────┬────────┘      │  - reference    │                │
│                          │               └────────┬────────┘                │
└──────────────────────────┼────────────────────────┼─────────────────────────┘
                           │                        │
                           ▼                        ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         PARSE & STANDARDIZE                                  │
│                         (parse_outputs.py)                                   │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│    Raw Tool Outputs ──────────────────────────► Standardized DataFrame       │
│                                                                              │
│    TOMBO CSV:                    ┌─────────────────────────────────────┐    │
│    pos,stat,valid_cov            │  tool      │ reference │ position   │    │
│    142,0.003,45          ──►     │  tombo     │ 16s       │ 142        │    │
│                                  │  score     │ score_type│ replicate  │    │
│    DRUMMER summary.txt:          │  0.003     │ pvalue    │ rep1       │    │
│    chrom,position,pval           └─────────────────────────────────────┘    │
│    16s,1402,0.001        ──►     (Same standardized format)                 │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                            ANALYSIS PIPELINE                                 │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│    ┌─────────────────────────────────────────────────────────────────┐      │
│    │                                                                 │      │
│    │    Standardized      ┌──────────────┐      Ground Truth        │      │
│    │    Tool Outputs ────►│  BENCHMARK   │◄──── (known sites)       │      │
│    │                      │   METRICS    │                          │      │
│    │                      └──────┬───────┘                          │      │
│    │                             │                                  │      │
│    │                             ▼                                  │      │
│    │                      ┌──────────────┐                          │      │
│    │                      │   AUPRC      │                          │      │
│    │                      │   AUROC      │                          │      │
│    │                      │   F1, P, R   │                          │      │
│    │                      └──────────────┘                          │      │
│    │                                                                 │      │
│    └─────────────────────────────────────────────────────────────────┘      │
│                                                                              │
│    ┌─────────────────────────────────────────────────────────────────┐      │
│    │                                                                 │      │
│    │    Tool Output       ┌──────────────┐                          │      │
│    │    (per replicate)──►│  REPLICATE   │                          │      │
│    │                      │   ANALYSIS   │                          │      │
│    │                      └──────┬───────┘                          │      │
│    │                             │                                  │      │
│    │                             ▼                                  │      │
│    │                      ┌──────────────┐                          │      │
│    │                      │  Consensus   │                          │      │
│    │                      │  Jaccard     │                          │      │
│    │                      │  Correlation │                          │      │
│    │                      └──────────────┘                          │      │
│    │                                                                 │      │
│    └─────────────────────────────────────────────────────────────────┘      │
│                                                                              │
│    ┌─────────────────────────────────────────────────────────────────┐      │
│    │                                                                 │      │
│    │    All Tool          ┌──────────────┐                          │      │
│    │    Outputs ─────────►│    TOOL      │                          │      │
│    │                      │  COMPARISON  │                          │      │
│    │                      └──────┬───────┘                          │      │
│    │                             │                                  │      │
│    │                             ▼                                  │      │
│    │                      ┌──────────────┐                          │      │
│    │                      │  UpSet data  │                          │      │
│    │                      │  Venn data   │                          │      │
│    │                      │  Rankings    │                          │      │
│    │                      └──────────────┘                          │      │
│    │                                                                 │      │
│    └─────────────────────────────────────────────────────────────────┘      │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                              OUTPUTS                                         │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Detailed Data Flow

### Stage 1: Input Parsing

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ TOOL RAW OUTPUTS                                                             │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│ TOMBO (CSV)              DRUMMER (summary.txt)        DIFFERR (BED)         │
│ ┌──────────────────┐     ┌─────────────────────┐     ┌──────────────────┐   │
│ │pos,stat,valid_cov│     │chrom,pos,pval,oddR  │     │chr,start,end,fdr │   │
│ │142,0.003,45      │     │16s,1402,0.001,3.2   │     │16s,141,142,0.01  │   │
│ │143,0.012,42      │     │16s,1518,0.05,1.8    │     │16s,1401,1402,0.02│   │
│ └────────┬─────────┘     └──────────┬──────────┘     └────────┬─────────┘   │
│          │                          │                          │             │
│          │    parse_tombo()         │    parse_drummer()       │             │
│          │                          │                          │             │
│          ▼                          ▼                          ▼             │
│ ┌────────────────────────────────────────────────────────────────────────┐  │
│ │                    STANDARDIZED DATAFRAME                               │  │
│ ├────────┬───────────┬──────────┬────────┬────────────┬─────────────────┤  │
│ │ tool   │ reference │ position │ score  │ score_type │ replicate       │  │
│ ├────────┼───────────┼──────────┼────────┼────────────┼─────────────────┤  │
│ │ tombo  │ 16s       │ 142      │ 0.003  │ pvalue     │ rep1            │  │
│ │ tombo  │ 16s       │ 143      │ 0.012  │ pvalue     │ rep1            │  │
│ │ drummer│ 16s       │ 1402     │ 0.001  │ pvalue     │ rep1            │  │
│ │ drummer│ 16s       │ 1518     │ 0.05   │ pvalue     │ rep1            │  │
│ │ differr│ 16s       │ 142      │ 0.01   │ fdr        │ rep1            │  │
│ │ differr│ 16s       │ 1402     │ 0.02   │ fdr        │ rep1            │  │
│ └────────┴───────────┴──────────┴────────┴────────────┴─────────────────┘  │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

### Stage 2: Benchmark Metrics Calculation

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ BENCHMARK METRICS FLOW                                                       │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  INPUTS:                                                                     │
│  ┌─────────────────────┐    ┌─────────────────────────────────┐             │
│  │ Tool Output         │    │ Ground Truth (known_sites.csv)  │             │
│  │ (standardized)      │    │                                 │             │
│  │                     │    │ reference,position,mod_type     │             │
│  │ position │ score    │    │ 16s,142,m6A                     │             │
│  │ 142      │ 0.003    │    │ 16s,1402,m6A                    │             │
│  │ 143      │ 0.012    │    │ 16s,1518,Nm                     │             │
│  │ 500      │ 0.8      │    │ 23s,2445,m6A                    │             │
│  │ 1402     │ 0.001    │    │                                 │             │
│  └──────────┬──────────┘    └───────────────┬─────────────────┘             │
│             │                               │                                │
│             └───────────────┬───────────────┘                                │
│                             │                                                │
│                             ▼                                                │
│             ┌───────────────────────────────┐                                │
│             │    LABEL ASSIGNMENT           │                                │
│             │                               │                                │
│             │  For each tool position:      │                                │
│             │  - Match to ground truth?     │                                │
│             │    YES → label = 1 (positive) │                                │
│             │    NO  → label = 0 (negative) │                                │
│             └───────────────┬───────────────┘                                │
│                             │                                                │
│                             ▼                                                │
│             ┌───────────────────────────────┐                                │
│             │  position │ score  │ label    │                                │
│             │  142      │ 0.003  │ 1 (GT)   │  ◄── In ground truth           │
│             │  143      │ 0.012  │ 0        │  ◄── Not in ground truth       │
│             │  500      │ 0.8    │ 0        │                                │
│             │  1402     │ 0.001  │ 1 (GT)   │  ◄── In ground truth           │
│             └───────────────┬───────────────┘                                │
│                             │                                                │
│                             ▼                                                │
│             ┌───────────────────────────────┐                                │
│             │    SCORE TRANSFORMATION       │                                │
│             │                               │                                │
│             │  p-value/FDR: -log10(score)   │                                │
│             │  (lower p-value = higher rank)│                                │
│             │                               │                                │
│             │  z-score: abs(score)          │                                │
│             └───────────────┬───────────────┘                                │
│                             │                                                │
│                             ▼                                                │
│             ┌───────────────────────────────┐                                │
│             │    METRICS CALCULATION        │                                │
│             │    (sklearn)                  │                                │
│             └───────────────┬───────────────┘                                │
│                             │                                                │
│                             ▼                                                │
│  OUTPUT:    ┌───────────────────────────────────────────────────────────┐   │
│             │ MetricsResult                                              │   │
│             │ ─────────────────────────────────────────────────────────  │   │
│             │ • tool: "drummer"                                          │   │
│             │ • reference: "16s"                                         │   │
│             │ • auprc: 0.847        ◄── Area under PR curve              │   │
│             │ • auroc: 0.912        ◄── Area under ROC curve             │   │
│             │ • f1_optimal: 0.78    ◄── Best F1 score                    │   │
│             │ • precision_optimal: 0.82                                  │   │
│             │ • recall_optimal: 0.75                                     │   │
│             │ • optimal_threshold: 2.5  ◄── -log10(pval) threshold       │   │
│             │ • n_true_positives: 45                                     │   │
│             │ • n_false_positives: 10                                    │   │
│             │ • n_false_negatives: 15                                    │   │
│             │ • n_ground_truth: 60                                       │   │
│             └───────────────────────────────────────────────────────────┘   │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

**What the metrics mean:**
| Metric | Meaning | Good Value |
|--------|---------|------------|
| **AUPRC** | Overall ranking quality (how well high scores correspond to true modifications) | >0.5 (higher = better) |
| **AUROC** | Discrimination ability (separating true from false) | >0.5 (higher = better) |
| **F1** | Balance of precision and recall at optimal threshold | >0.7 (higher = better) |
| **Precision** | What fraction of calls are correct? | Higher = fewer false positives |
| **Recall** | What fraction of true sites found? | Higher = fewer missed sites |

---

### Stage 3: Replicate Analysis

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ REPLICATE ANALYSIS FLOW                                                      │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  INPUT: Tool output with 3 replicates                                        │
│                                                                              │
│  ┌────────────────┐  ┌────────────────┐  ┌────────────────┐                 │
│  │ Replicate 1    │  │ Replicate 2    │  │ Replicate 3    │                 │
│  │ ──────────     │  │ ──────────     │  │ ──────────     │                 │
│  │ pos: 142  ✓    │  │ pos: 142  ✓    │  │ pos: 142  ✓    │                 │
│  │ pos: 143       │  │ pos: 143       │  │ pos: 200       │                 │
│  │ pos: 500  ✓    │  │ pos: 500  ✓    │  │ pos: 500  ✓    │                 │
│  │ pos: 1402 ✓    │  │ pos: 1402 ✓    │  │ pos: 1402 ✓    │                 │
│  │ pos: 800       │  │ pos: 900       │  │ pos: 1000      │                 │
│  └───────┬────────┘  └───────┬────────┘  └───────┬────────┘                 │
│          │                   │                   │                          │
│          └───────────────────┼───────────────────┘                          │
│                              │                                               │
│                              ▼                                               │
│          ┌───────────────────────────────────────┐                          │
│          │         CONSENSUS CALLING             │                          │
│          │         (min_replicates=2)            │                          │
│          └───────────────────┬───────────────────┘                          │
│                              │                                               │
│                              ▼                                               │
│          ┌───────────────────────────────────────────────────────────┐      │
│          │ Consensus Sites (in ≥2 replicates)                        │      │
│          │ ─────────────────────────────────────────────────────     │      │
│          │ position │ n_replicates │ mean_score │ replicates         │      │
│          │ 142      │ 3            │ 0.004      │ rep1,rep2,rep3     │      │
│          │ 143      │ 2            │ 0.015      │ rep1,rep2          │      │
│          │ 500      │ 3            │ 0.002      │ rep1,rep2,rep3     │      │
│          │ 1402     │ 3            │ 0.001      │ rep1,rep2,rep3     │      │
│          └───────────────────────────────────────────────────────────┘      │
│                                                                              │
│          ┌───────────────────────────────────────┐                          │
│          │         CONCORDANCE METRICS           │                          │
│          └───────────────────┬───────────────────┘                          │
│                              │                                               │
│                              ▼                                               │
│          ┌───────────────────────────────────────────────────────────┐      │
│          │ Pairwise Jaccard Similarity                               │      │
│          │ ─────────────────────────────────────────────────────     │      │
│          │                                                           │      │
│          │ Jaccard = |A ∩ B| / |A ∪ B|                               │      │
│          │                                                           │      │
│          │ Rep1 vs Rep2: 0.75  (3 shared / 4 union)                  │      │
│          │ Rep1 vs Rep3: 0.60  (3 shared / 5 union)                  │      │
│          │ Rep2 vs Rep3: 0.67  (3 shared / 4.5 union)                │      │
│          │                                                           │      │
│          │ Mean Jaccard: 0.67 ± 0.08                                 │      │
│          └───────────────────────────────────────────────────────────┘      │
│                                                                              │
│  OUTPUT:                                                                     │
│  ┌───────────────────────────────────────────────────────────────────────┐  │
│  │ ReplicateStats                                                         │  │
│  │ ───────────────────────────────────────────────────────────────────    │  │
│  │ • tool: "drummer"                                                      │  │
│  │ • n_replicates: 3                                                      │  │
│  │ • n_sites_per_replicate: [5, 5, 5]                                     │  │
│  │ • n_consensus_sites: 4     ◄── Sites in ≥2 replicates                  │  │
│  │ • mean_jaccard: 0.67       ◄── Reproducibility score                   │  │
│  │ • std_jaccard: 0.08                                                    │  │
│  └───────────────────────────────────────────────────────────────────────┘  │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

**What the replicate metrics mean:**
| Metric | Meaning | Good Value |
|--------|---------|------------|
| **n_consensus_sites** | Sites reproducibly detected | Higher = more reliable |
| **mean_jaccard** | How similar are replicate calls? | >0.6 = good reproducibility |

---

### Stage 4: Tool Comparison

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ TOOL COMPARISON FLOW                                                         │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  INPUT: Significant positions from each tool (at threshold)                  │
│                                                                              │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐        │
│  │   TOMBO     │  │  DRUMMER    │  │   ELIGOS    │  │   XPORE     │        │
│  │   ─────     │  │   ──────    │  │   ──────    │  │   ─────     │        │
│  │   pos: 142  │  │   pos: 142  │  │   pos: 142  │  │   pos: 200  │        │
│  │   pos: 500  │  │   pos: 500  │  │   pos: 600  │  │   pos: 500  │        │
│  │   pos: 1402 │  │   pos: 1402 │  │   pos: 1402 │  │   pos: 1402 │        │
│  │   pos: 800  │  │   pos: 900  │  │             │  │   pos: 1000 │        │
│  └──────┬──────┘  └──────┬──────┘  └──────┬──────┘  └──────┬──────┘        │
│         │                │                │                │                 │
│         └────────────────┴────────────────┴────────────────┘                 │
│                                    │                                         │
│                                    ▼                                         │
│                    ┌───────────────────────────────┐                         │
│                    │      SET OPERATIONS           │                         │
│                    └───────────────┬───────────────┘                         │
│                                    │                                         │
│                                    ▼                                         │
│  ┌───────────────────────────────────────────────────────────────────────┐  │
│  │                                                                        │  │
│  │    TOMBO ────────┐                                                     │  │
│  │    {142,500,     │                                                     │  │
│  │     1402,800}    │    UNION: {142,200,500,600,800,900,1000,1402}       │  │
│  │                  ├──► = 8 unique positions                             │  │
│  │    DRUMMER ──────┤                                                     │  │
│  │    {142,500,     │    INTERSECTION: {142, 500, 1402}                   │  │
│  │     1402,900}    │    = 3 positions detected by ALL tools              │  │
│  │                  │                                                     │  │
│  │    ELIGOS ───────┤                                                     │  │
│  │    {142,600,     │                                                     │  │
│  │     1402}        │                                                     │  │
│  │                  │                                                     │  │
│  │    XPORE ────────┘                                                     │  │
│  │    {200,500,                                                           │  │
│  │     1402,1000}                                                         │  │
│  │                                                                        │  │
│  └───────────────────────────────────────────────────────────────────────┘  │
│                                    │                                         │
│                                    ▼                                         │
│  ┌───────────────────────────────────────────────────────────────────────┐  │
│  │  UpSet Plot Data                                                       │  │
│  │  ────────────────────────────────────────────────────────              │  │
│  │  position │ TOMBO │ DRUMMER │ ELIGOS │ XPORE                           │  │
│  │  142      │   1   │    1    │   1    │   0     ◄── 3 tools             │  │
│  │  200      │   0   │    0    │   0    │   1     ◄── 1 tool (xpore only) │  │
│  │  500      │   1   │    1    │   0    │   1     ◄── 3 tools             │  │
│  │  600      │   0   │    0    │   1    │   0     ◄── 1 tool (eligos only)│  │
│  │  800      │   1   │    0    │   0    │   0     ◄── 1 tool (tombo only) │  │
│  │  900      │   0   │    1    │   0    │   0     ◄── 1 tool              │  │
│  │  1000     │   0   │    0    │   0    │   1     ◄── 1 tool              │  │
│  │  1402     │   1   │    1    │   1    │   1     ◄── ALL 4 tools         │  │
│  └───────────────────────────────────────────────────────────────────────┘  │
│                                    │                                         │
│                                    ▼                                         │
│  OUTPUT:                                                                     │
│  ┌───────────────────────────────────────────────────────────────────────┐  │
│  │ ToolComparisonResult                                                   │  │
│  │ ───────────────────────────────────────────────────────────────────    │  │
│  │ • n_tools: 4                                                           │  │
│  │ • n_union: 8              ◄── Total unique positions                   │  │
│  │ • n_intersection: 1       ◄── Positions found by ALL tools             │  │
│  │ • positions_by_n_tools:                                                │  │
│  │   - 1 tool:  4 positions  ◄── Tool-specific calls                      │  │
│  │   - 2 tools: 0 positions                                               │  │
│  │   - 3 tools: 3 positions  ◄── High-confidence calls                    │  │
│  │   - 4 tools: 1 position   ◄── Highest-confidence                       │  │
│  │ • pairwise_jaccard:                                                    │  │
│  │   - TOMBO-DRUMMER: 0.60                                                │  │
│  │   - TOMBO-ELIGOS: 0.40                                                 │  │
│  │   - DRUMMER-XPORE: 0.50                                                │  │
│  │   ...                                                                  │  │
│  └───────────────────────────────────────────────────────────────────────┘  │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

### Stage 5: Visualization & Output

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ VISUALIZATION OUTPUTS                                                        │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │  PR CURVES (pr_curves.png)                                          │    │
│  │  ─────────────────────────                                          │    │
│  │                                                                     │    │
│  │  Precision │    ╭──── DRUMMER (AUPRC=0.85)                          │    │
│  │      1.0   │   ╱  ╲                                                 │    │
│  │            │  ╱    ╲──── TOMBO (AUPRC=0.72)                         │    │
│  │      0.5   │ ╱      ╲                                               │    │
│  │            │╱        ╲──── ELIGOS (AUPRC=0.65)                      │    │
│  │      0.0   └────────────────────────────────                        │    │
│  │            0.0       0.5        1.0                                 │    │
│  │                    Recall                                           │    │
│  │                                                                     │    │
│  │  INTERPRETATION: Higher curve = better tool                         │    │
│  │  DRUMMER has best precision at all recall levels                    │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
│                                                                              │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │  UPSET PLOT (upset_plot.png)                                        │    │
│  │  ───────────────────────────                                        │    │
│  │                                                                     │    │
│  │  Count │                                                            │    │
│  │    4   │  ████                  ◄── 4 positions by 1 tool only      │    │
│  │    3   │        ████            ◄── 3 positions by 3 tools          │    │
│  │    1   │              ████      ◄── 1 position by all 4 tools       │    │
│  │        └─────────────────────                                       │    │
│  │           ●         ●    ●      TOMBO                               │    │
│  │           ●         ●    ●      DRUMMER                             │    │
│  │                     ●    ●      ELIGOS                              │    │
│  │                          ●      XPORE                               │    │
│  │                                                                     │    │
│  │  INTERPRETATION: Most positions detected by only 1 tool             │    │
│  │  Only 1 position detected by all tools (high confidence)            │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
│                                                                              │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │  CORRELATION HEATMAP (tool_correlation.png)                         │    │
│  │  ──────────────────────────────────────────                         │    │
│  │                                                                     │    │
│  │           TOMBO  DRUMMER  ELIGOS  XPORE                             │    │
│  │  TOMBO    1.00    0.60     0.40   0.45                              │    │
│  │  DRUMMER  0.60    1.00     0.55   0.50                              │    │
│  │  ELIGOS   0.40    0.55     1.00   0.35                              │    │
│  │  XPORE    0.45    0.50     0.35   1.00                              │    │
│  │                                                                     │    │
│  │  INTERPRETATION: Values = Jaccard similarity                        │    │
│  │  TOMBO-DRUMMER have highest agreement (0.60)                        │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Complete Output Directory Structure

```
downstream_analysis/
│
├── metrics/
│   ├── all_metrics.csv              # All tools × all references
│   │   ┌─────────────────────────────────────────────────────────┐
│   │   │ tool,reference,auprc,auroc,f1_optimal,precision,recall  │
│   │   │ drummer,16s,0.847,0.912,0.78,0.82,0.75                  │
│   │   │ tombo,16s,0.723,0.856,0.65,0.70,0.61                    │
│   │   │ eligos,16s,0.651,0.801,0.58,0.62,0.55                   │
│   │   └─────────────────────────────────────────────────────────┘
│   │
│   └── tool_ranking.csv             # Tools ranked by AUPRC
│       ┌───────────────────────────────────────────┐
│       │ rank,tool,auprc,auroc,f1                  │
│       │ 1,drummer,0.847,0.912,0.78                │
│       │ 2,tombo,0.723,0.856,0.65                  │
│       │ 3,eligos,0.651,0.801,0.58                 │
│       └───────────────────────────────────────────┘
│
├── replicate_analysis/
│   ├── concordance.csv              # Jaccard per tool
│   │   ┌─────────────────────────────────────────────────────────┐
│   │   │ tool,n_replicates,mean_jaccard,std_jaccard,n_consensus  │
│   │   │ drummer,3,0.67,0.08,45                                  │
│   │   │ tombo,3,0.58,0.12,38                                    │
│   │   └─────────────────────────────────────────────────────────┘
│   │
│   ├── drummer_consensus.csv        # Consensus sites per tool
│   │   ┌───────────────────────────────────────────────────────┐
│   │   │ position,n_replicates,mean_score,replicates           │
│   │   │ 142,3,0.004,rep1,rep2,rep3                            │
│   │   │ 1402,3,0.001,rep1,rep2,rep3                           │
│   │   └───────────────────────────────────────────────────────┘
│   │
│   └── tombo_consensus.csv
│
├── tool_comparison/
│   ├── summary.csv                  # High-level comparison
│   │   ┌─────────────────────────────┐
│   │   │ metric,value                │
│   │   │ n_tools,9                   │
│   │   │ n_union,523                 │
│   │   │ n_intersection,12           │
│   │   │ n_drummer,89                │
│   │   │ n_tombo,76                  │
│   │   └─────────────────────────────┘
│   │
│   ├── upset_data.csv               # For UpSet plot
│   ├── pairwise_agreement.csv       # Tool-tool Jaccard
│   ├── consensus_2tools.csv         # Sites by ≥2 tools
│   ├── consensus_3tools.csv         # Sites by ≥3 tools
│   └── consensus_9tools.csv         # Sites by ALL tools
│
├── visualizations/
│   ├── pr_curves.png                # Precision-Recall curves
│   ├── roc_curves.png               # ROC curves
│   ├── score_distributions.png      # Score histograms per tool
│   ├── upset_plot.png               # Tool agreement UpSet
│   ├── venn_diagram.png             # 3-tool Venn
│   ├── metrics_comparison.png       # Bar chart of AUPRC
│   ├── replicate_concordance.png    # Jaccard bar chart
│   └── tool_correlation.png         # Heatmap
│
└── report/
    ├── report.html                  # Comprehensive HTML report
    ├── pr_curves.png
    ├── roc_curves.png
    └── ...
```

---

## Quick Reference: What Each Output Tells You

| File | Question Answered |
|------|-------------------|
| `all_metrics.csv` | "Which tool is most accurate?" (compare AUPRC) |
| `tool_ranking.csv` | "What's the ranking of tools?" |
| `concordance.csv` | "Which tool is most reproducible across replicates?" |
| `*_consensus.csv` | "Which sites are reliably detected?" |
| `pairwise_agreement.csv` | "Which tools agree with each other?" |
| `consensus_Ntools.csv` | "What are the highest-confidence modification sites?" |
| `pr_curves.png` | "How does precision trade off with recall?" |
| `upset_plot.png` | "How do tool calls overlap?" |
| `tool_correlation.png` | "Are tools detecting similar things?" |
