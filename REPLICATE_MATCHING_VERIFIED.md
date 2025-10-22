# Replicate Matching Logic - VERIFIED ✅

## Your Question:
"When comparing native and ivt → it should be comparison of native rep1 data with ivt rep1 data and native rep2 with ivt rep2 and native rep3 with ivt rep3"

## Answer: ✅ YES, THIS IS EXACTLY HOW IT WORKS!

---

## How the Grouping Works:

All three modification calling tools (Tombo, Yanocomp, Xpore) use the **SAME grouping key**:

```groovy
def key = "${meta.rrna}_${meta.replicate}"
```

### For Your 3 Replicates:

#### 16S rRNA Comparisons:
| Native Sample | IVT Sample | Key Used | Result |
|---------------|------------|----------|---------|
| native_rep1 (16s) | ivt_rep1 (16s) | `16s_rep1` | ✅ Compared together |
| native_rep2 (16s) | ivt_rep2 (16s) | `16s_rep2` | ✅ Compared together |
| native_rep3 (16s) | ivt_rep3 (16s) | `16s_rep3` | ✅ Compared together |

#### 23S rRNA Comparisons:
| Native Sample | IVT Sample | Key Used | Result |
|---------------|------------|----------|---------|
| native_rep1 (23s) | ivt_rep1 (23s) | `23s_rep1` | ✅ Compared together |
| native_rep2 (23s) | ivt_rep2 (23s) | `23s_rep2` | ✅ Compared together |
| native_rep3 (23s) | ivt_rep3 (23s) | `23s_rep3` | ✅ Compared together |

---

## Total Comparisons Made:

You will get **6 comparison results**:
1. 16s_rep1: native_rep1 vs ivt_rep1
2. 16s_rep2: native_rep2 vs ivt_rep2
3. 16s_rep3: native_rep3 vs ivt_rep3
4. 23s_rep1: native_rep1 vs ivt_rep1
5. 23s_rep2: native_rep2 vs ivt_rep2
6. 23s_rep3: native_rep3 vs ivt_rep3

Each comparison will be performed by:
- ✅ Tombo (CSV output)
- ✅ Yanocomp (BED output)
- ✅ Xpore (table output)

---

## Code Evidence:

### 1. TOMBO (prepare_signal_data.nf, line 100):
```groovy
TOMBO_RESQUIGGLE.out.resquiggled
    .map { meta, fast5 ->
        def key = "${meta.rrna}_${meta.replicate}"  // ← Groups by rRNA + replicate
        [ key, meta.type, fast5 ]
    }
    .groupTuple()  // ← Groups samples with same key together
```

### 2. YANOCOMP (modification_calling.nf, line 54):
```groovy
YANOCOMP_PREPARE.out.hdf5
    .map { meta, hdf5 ->
        def key = "${meta.rrna}_${meta.replicate}"  // ← Same grouping
        [ key, meta.type, hdf5 ]
    }
    .groupTuple()
```

### 3. XPORE (modification_calling.nf, line 88):
```groovy
XPORE_DATAPREP.out.dataprep
    .map { meta, dataprep_dir ->
        def key = "${meta.rrna}_${meta.replicate}"  // ← Same grouping
        [ key, meta.type, dataprep_dir ]
    }
    .groupTuple()
```

---

## What `.groupTuple()` Does:

The `.groupTuple()` function groups all items with the **same key** together.

**Example for 16s_rep1:**
```
Input after .map:
[ "16s_rep1", "native", native_fast5_rep1 ]
[ "16s_rep1", "ivt",    ivt_fast5_rep1 ]

After .groupTuple():
[ "16s_rep1", ["native", "ivt"], [native_fast5_rep1, ivt_fast5_rep1] ]
                   ↑                         ↑
            types are grouped        files are grouped
```

Then the code finds the native and ivt within the group:
```groovy
def native_idx = types.findIndexOf { it == 'native' }
def ivt_idx = types.findIndexOf { it == 'ivt' }
```

And compares them: `files[native_idx]` vs `files[ivt_idx]`

---

## Expected Output Files:

### Tombo CSV files (in results/modification_calling/):
```
tombo_16s_rep1.csv   ← native_rep1 vs ivt_rep1 for 16S
tombo_16s_rep2.csv   ← native_rep2 vs ivt_rep2 for 16S
tombo_16s_rep3.csv   ← native_rep3 vs ivt_rep3 for 16S
tombo_23s_rep1.csv   ← native_rep1 vs ivt_rep1 for 23S
tombo_23s_rep2.csv   ← native_rep2 vs ivt_rep2 for 23S
tombo_23s_rep3.csv   ← native_rep3 vs ivt_rep3 for 23S
```

### Yanocomp BED files:
```
yanocomp_16s_rep1.bed   ← native_rep1 vs ivt_rep1 for 16S
yanocomp_16s_rep2.bed   ← native_rep2 vs ivt_rep2 for 16S
yanocomp_16s_rep3.bed   ← native_rep3 vs ivt_rep3 for 16S
yanocomp_23s_rep1.bed   ← native_rep1 vs ivt_rep1 for 23S
yanocomp_23s_rep2.bed   ← native_rep2 vs ivt_rep2 for 23S
yanocomp_23s_rep3.bed   ← native_rep3 vs ivt_rep3 for 23S
```

### Xpore tables (similar naming):
```
xpore_16s_rep1/diffmod.table
xpore_16s_rep2/diffmod.table
xpore_16s_rep3/diffmod.table
xpore_23s_rep1/diffmod.table
xpore_23s_rep2/diffmod.table
xpore_23s_rep3/diffmod.table
```

---

## ✅ CONCLUSION:

**Your concern is already handled correctly!**

The pipeline will **automatically**:
1. Group samples by rRNA type (16s/23s) AND replicate (rep1/rep2/rep3)
2. Compare matching replicates: native_rep1 vs ivt_rep1, etc.
3. Generate separate results for each replicate comparison

**No changes needed** - the logic is already exactly what you want! 🎉

---

## Note About Shared FAST5 Directories:

Your setup with shared FAST5 directories is also fine:
- `/home/bmorampa/k12_native_fast5` (shared by all native replicates)
- `/home/bmorampa/k12_ivt_fast5` (shared by all ivt replicates)

The `FAST5_SUBSET` module will:
1. Extract read IDs from each replicate's mapped FASTQ
2. Find only those specific reads in the shared FAST5 directory
3. Create replicate-specific FAST5 subsets

So each replicate will get its own subset of FAST5 files extracted from the shared directory.

---

**Everything is configured correctly for your 3-replicate experimental design!** ✅
