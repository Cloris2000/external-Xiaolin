# Nextflow Channel Operators: Detailed Explanation

## Overview

Nextflow channels are asynchronous queues that connect processes together. Channel operators transform and combine data streams to create the combinations needed for parallel processing. This document explains how `cross`, `join`, and `groupTuple` are used in the GWAS pipeline.

---

## 1. The `cross` Operator

### Purpose
Creates the **Cartesian product** of two channels - every item from channel A paired with every item from channel B.

### Pipeline Usage

```groovy
// Create channels for cohorts and cell types
def cohort_ch = Channel.fromList(params.cohorts)  // ["ROSMAP", "Mayo", "MSBB"]
def celltype_ch = Channel.fromList(cell_types)     // ["Astrocyte", "Microglia", ...]

// Create all combinations
def cohort_celltype_ch = cohort_ch.cross(celltype_ch)
```

### What This Creates

If you have:
- **Cohorts**: `["ROSMAP", "Mayo"]`
- **Cell Types**: `["Astrocyte", "Microglia"]`

The `cross` operator produces:
```
[ROSMAP, Astrocyte]
[ROSMAP, Microglia]
[Mayo, Astrocyte]
[Mayo, Microglia]
```

### Visual Representation

```
cohort_ch:     [ROSMAP] ──┐
                          ├──→ cross ──→ [ROSMAP, Astrocyte]
celltype_ch:   [Astrocyte]─┘              [ROSMAP, Microglia]
                                        [Mayo, Astrocyte]
                                        [Mayo, Microglia]
```

### Why We Need This

Instead of writing nested loops:
```groovy
// BAD: Sequential processing
for (cohort in cohorts) {
    for (cell_type in cell_types) {
        run_gwas(cohort, cell_type)  // Runs one at a time
    }
}
```

We get parallel execution:
```groovy
// GOOD: Parallel processing
cohort_ch.cross(celltype_ch).map { cohort, cell_type ->
    run_gwas(cohort, cell_type)  // All combinations run in parallel!
}
```

**Result**: 2 cohorts × 2 cell types = 4 GWAS runs that execute **simultaneously** instead of sequentially.

---

## 2. The `join` Operator

### Purpose
**Merges** two channels based on a **matching key**, similar to a SQL JOIN operation.

### Pipeline Usage

```groovy
// Step 1: Create a keyed channel from pred_list files
def pred_list_keyed = regenie_step1_ch.out.pred_list
    .map { pred_list_file ->
        // Extract cohort and cell_type from filename
        def basename = pred_list_file.toString().split('/').last()
        def parts = basename.replace('_step1_pred.list', '').split('_')
        def cohort = parts[0]
        def cell_type = parts.drop(1).join('_')
        
        // Create key-value pair: [key, value]
        [[cohort, cell_type], pred_list_file]
    }

// Step 2: Create keyed channel from cohort_celltype combinations
def cohort_celltype_keyed = cohort_celltype_ch
    .map { cohort, cell_type -> 
        [[cohort, cell_type], [cohort, cell_type]] 
    }

// Step 3: Join them on the key
def regenie_step2_input = cohort_celltype_keyed
    .join(pred_list_keyed, by: 0)  // Join on first element (the key)
```

### What This Does

**Before Join:**

`cohort_celltype_keyed`:
```
[[ROSMAP, Astrocyte], [ROSMAP, Astrocyte]]
[[ROSMAP, Microglia], [ROSMAP, Microglia]]
[[Mayo, Astrocyte], [Mayo, Astrocyte]]
```

`pred_list_keyed`:
```
[[ROSMAP, Astrocyte], /path/to/ROSMAP_Astrocyte_step1_pred.list]
[[ROSMAP, Microglia], /path/to/ROSMAP_Microglia_step1_pred.list]
[[Mayo, Astrocyte], /path/to/Mayo_Astrocyte_step1_pred.list]
```

**After Join (by: 0)**:

```
[[ROSMAP, Astrocyte], [ROSMAP, Astrocyte], /path/to/ROSMAP_Astrocyte_step1_pred.list]
[[ROSMAP, Microglia], [ROSMAP, Microglia], /path/to/ROSMAP_Microglia_step1_pred.list]
[[Mayo, Astrocyte], [Mayo, Astrocyte], /path/to/Mayo_Astrocyte_step1_pred.list]
```

### Visual Representation

```
cohort_celltype_keyed:  [key: [ROSMAP, Astrocyte], value: [ROSMAP, Astrocyte]]
                              │
                              ├─── join (match on key) ────┐
                              │                            │
pred_list_keyed:              [key: [ROSMAP, Astrocyte], value: pred_list_file]
                                                           │
                                                           ↓
Result: [ROSMAP, Astrocyte, pred_list_file]
```

### Why We Need This

**Problem**: Step 2 needs both:
1. The cohort and cell_type information
2. The pred_list file from Step 1

**Solution**: Join ensures each Step 2 process gets the **correct** pred_list file that matches its cohort and cell_type.

**Without join**, you'd have to manually match files, which is error-prone:
```groovy
// BAD: Manual matching
def pred_list = find_file_for(cohort, cell_type)  // Error-prone!
```

**With join**, Nextflow automatically matches based on the key:
```groovy
// GOOD: Automatic matching
.join(pred_list_keyed, by: 0)  // Nextflow handles matching!
```

---

## 3. The `groupTuple` Operator

### Purpose
**Groups** items from a channel by a key, collecting all values with the same key into a tuple.

### Pipeline Usage

```groovy
// Step 1: Extract cell_type from each raw_p file and create key-value pairs
def raw_p_by_celltype = raw_p_ch.out.raw_p_file
    .map { raw_p_file ->
        def file_name = raw_p_file.toString()
        def basename = file_name.split('/').last()
        def parts = basename.replace('_step2.regenie.raw_p', '').split('_')
        def cohort = parts[0]
        def cell_type_raw = parts.drop(1).join('_')
        def cell_type = cell_type_mapping.find { k, v -> v == cell_type_raw }?.key ?: cell_type_raw
        
        // Create key-value pair
        [cell_type, raw_p_file]
    }
    .groupTuple(by: 0)  // Group by first element (cell_type)
```

### What This Creates

**Before groupTuple:**

```
[Astrocyte, /path/to/ROSMAP_Astrocyte_step2.regenie.raw_p]
[Astrocyte, /path/to/Mayo_Astrocyte_step2.regenie.raw_p]
[Astrocyte, /path/to/MSBB_Astrocyte_step2.regenie.raw_p]
[Microglia, /path/to/ROSMAP_Microglia_step2.regenie.raw_p]
[Microglia, /path/to/Mayo_Microglia_step2.regenie.raw_p]
```

**After groupTuple(by: 0):**

```
[Astrocyte, [/path/to/ROSMAP_Astrocyte_step2.regenie.raw_p,
             /path/to/Mayo_Astrocyte_step2.regenie.raw_p,
             /path/to/MSBB_Astrocyte_step2.regenie.raw_p]]
[Microglia, [/path/to/ROSMAP_Microglia_step2.regenie.raw_p,
             /path/to/Mayo_Microglia_step2.regenie.raw_p]]
```

### Visual Representation

```
Input Stream:
[Astrocyte, file1] ──┐
[Astrocyte, file2] ──┼──→ groupTuple ──→ [Astrocyte, [file1, file2, file3]]
[Astrocyte, file3] ──┘
[Microglia, file4] ──┐
[Microglia, file5] ──┼──→ groupTuple ──→ [Microglia, [file4, file5]]
```

### Why We Need This

**Problem**: METAL meta-analysis needs **all** raw_p files for a cell type from **all** cohorts, but they arrive as separate items.

**Solution**: `groupTuple` collects them into a single tuple per cell type.

**Without groupTuple**, you'd have to manually collect files:
```groovy
// BAD: Manual collection
def astrocyte_files = []
for (file in all_files) {
    if (file.contains("Astrocyte")) {
        astrocyte_files.add(file)
    }
}
// But when do you know you have all files? Timing issues!
```

**With groupTuple**, Nextflow automatically waits and groups:
```groovy
// GOOD: Automatic grouping
.groupTuple(by: 0)  // Nextflow waits for all items, then groups!
```

**Key Point**: `groupTuple` **waits** until all items with the same key have been emitted before creating the grouped tuple. This ensures meta-analysis has all cohort files.

---

## 4. Complete Data Flow Example

Let's trace a complete example through the pipeline:

### Starting Point

```groovy
params.cohorts = ["ROSMAP", "Mayo"]
cell_types = ["Astrocyte", "Microglia"]
```

### Stage 1: Create Combinations

```groovy
def cohort_ch = Channel.fromList(["ROSMAP", "Mayo"])
def celltype_ch = Channel.fromList(["Astrocyte", "Microglia"])
def cohort_celltype_ch = cohort_ch.cross(celltype_ch)
```

**Result**: 4 items
```
[ROSMAP, Astrocyte]
[ROSMAP, Microglia]
[Mayo, Astrocyte]
[Mayo, Microglia]
```

### Stage 2: Regenie Step 1 (runs in parallel)

Each combination triggers a Regenie Step 1 process:
- `ROSMAP_Astrocyte_step1` → produces `ROSMAP_Astrocyte_step1_pred.list`
- `ROSMAP_Microglia_step1` → produces `ROSMAP_Microglia_step1_pred.list`
- `Mayo_Astrocyte_step1` → produces `Mayo_Astrocyte_step1_pred.list`
- `Mayo_Microglia_step1` → produces `Mayo_Microglia_step1_pred.list`

### Stage 3: Join for Step 2

```groovy
// pred_list_keyed contains:
[[ROSMAP, Astrocyte], ROSMAP_Astrocyte_step1_pred.list]
[[ROSMAP, Microglia], ROSMAP_Microglia_step1_pred.list]
[[Mayo, Astrocyte], Mayo_Astrocyte_step1_pred.list]
[[Mayo, Microglia], Mayo_Microglia_step1_pred.list]

// Join with cohort_celltype_ch
.join(pred_list_keyed, by: 0)
```

**Result**: Each Step 2 process gets the correct pred_list file
```
[ROSMAP, Astrocyte, ROSMAP_Astrocyte_step1_pred.list]
[ROSMAP, Microglia, ROSMAP_Microglia_step1_pred.list]
[Mayo, Astrocyte, Mayo_Astrocyte_step1_pred.list]
[Mayo, Microglia, Mayo_Microglia_step1_pred.list]
```

### Stage 4: Regenie Step 2 (runs in parallel)

Each produces a `.regenie` file, then converted to `.raw_p`:
```
ROSMAP_Astrocyte_step2.regenie.raw_p
ROSMAP_Microglia_step2.regenie.raw_p
Mayo_Astrocyte_step2.regenie.raw_p
Mayo_Microglia_step2.regenie.raw_p
```

### Stage 5: Group by Cell Type

```groovy
// After mapping to extract cell_type:
[Astrocyte, ROSMAP_Astrocyte_step2.regenie.raw_p]
[Astrocyte, Mayo_Astrocyte_step2.regenie.raw_p]
[Microglia, ROSMAP_Microglia_step2.regenie.raw_p]
[Microglia, Mayo_Microglia_step2.regenie.raw_p]

// After groupTuple(by: 0):
[Astrocyte, [ROSMAP_Astrocyte_step2.regenie.raw_p, Mayo_Astrocyte_step2.regenie.raw_p]]
[Microglia, [ROSMAP_Microglia_step2.regenie.raw_p, Mayo_Microglia_step2.regenie.raw_p]]
```

### Stage 6: Meta-Analysis (runs in parallel)

Each cell type gets one meta-analysis process with all cohort files:
- `Astrocyte_meta_analysis` uses both ROSMAP and Mayo files
- `Microglia_meta_analysis` uses both ROSMAP and Mayo files

---

## 5. Key Benefits of This Approach

### 5.1 Automatic Parallelization

**Without channels**: You'd write nested loops and process sequentially
```groovy
for (cohort in cohorts) {
    for (cell_type in cell_types) {
        run_gwas(cohort, cell_type)  // Waits for each to finish
    }
}
// Total time: sum of all processing times
```

**With channels**: All combinations process simultaneously
```groovy
cohort_ch.cross(celltype_ch).map { cohort, cell_type ->
    run_gwas(cohort, cell_type)  // All run at once!
}
// Total time: max of processing times (for same resource level)
```

### 5.2 Automatic Dependency Management

Nextflow automatically ensures:
- Step 2 waits for Step 1 to complete
- Meta-analysis waits for all cohort files
- PRS calculation waits for PRS-cs to complete

**You don't write**: `if (step1_done) { run_step2() }`
**Nextflow handles**: Dependencies automatically via channel connections

### 5.3 Lazy Evaluation

Channels are **lazy** - they only process items as they're needed:
- If Step 1 fails for one combination, Step 2 doesn't try to process it
- If you only need results for one cell type, other cell types aren't processed
- Resources aren't wasted on unnecessary computations

### 5.4 Scalability

The same code works for:
- 1 cohort × 19 cell types = 19 combinations
- 4 cohorts × 19 cell types = 76 combinations
- 10 cohorts × 19 cell types = 190 combinations

**No code changes needed** - just update the configuration!

---

## 6. Common Patterns in the Pipeline

### Pattern 1: Create All Combinations
```groovy
channel1.cross(channel2)  // Every item from channel1 × every item from channel2
```

### Pattern 2: Match Related Items
```groovy
channel1.join(channel2, by: 0)  // Match items with same key
```

### Pattern 3: Collect Related Items
```groovy
channel.groupTuple(by: 0)  // Group items with same key into tuple
```

### Pattern 4: Transform Items
```groovy
channel.map { item ->
    // Extract information, create new structure
    [key, value]
}
```

---

## 7. Summary

| Operator | Purpose | Example Use Case |
|----------|---------|------------------|
| `cross` | Cartesian product | Create all cohort × cell_type combinations |
| `join` | Match by key | Connect Step 1 output to Step 2 input |
| `groupTuple` | Collect by key | Gather all cohort files for meta-analysis |
| `map` | Transform items | Extract information, create key-value pairs |

**Key Insight**: These operators enable the pipeline to automatically:
1. **Create** all necessary combinations
2. **Match** related items together
3. **Group** items that need to be processed together
4. **Execute** everything in parallel where possible

This is why the pipeline can handle 4 cohorts × 19 cell types = 76 GWAS runs efficiently, with automatic dependency management and parallel execution!

