# CPU Requirements Fix

## Problem

The pipeline failed with this error:
```
ERROR ~ Error executing process > 'COMBINED_PIPELINE:GWAS_PIPELINE:REGENIE_STEP1'
Caused by:
  Process requirement exceeds available CPUs -- req: 24; avail: 12
```

**Root cause**: The configs were set for 24 CPUs, but your system only has 12 CPUs available.

## Solution Applied to ROSMAP Config

✅ **Fixed** `nextflow.config.combined.rosmap`:
1. Changed `regenie_threads` from 24 → 12
2. Added process-specific CPU limits
3. Fixed `conda_env` warning

## Apply Same Fix to All Other Cohorts

You need to update the other 7 cohort configs too. Here's how:

### Quick Fix Script

```bash
cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

# Update all configs
for cohort in mayo msbb gtex nabec cmc_mssm cmc_penn cmc_pitt; do
    echo "Fixing nextflow.config.combined.${cohort}..."
    
    # Replace regenie_threads 24 with 12
    sed -i 's/regenie_threads = 24/regenie_threads = 12  \/\/ Reduced to match available CPUs/' \
        nextflow.config.combined.${cohort}
    
    # Add conda_env if missing
    if ! grep -q "conda_env" nextflow.config.combined.${cohort}; then
        sed -i '/study = /a\    \n    // Conda environment (optional)\n    conda_env = "test"' \
            nextflow.config.combined.${cohort}
    fi
done

echo "All configs updated!"
```

### Or Manual Fix

For each config file (`mayo`, `msbb`, `gtex`, `nabec`, `cmc_mssm`, `cmc_penn`, `cmc_pitt`):

**1. Change regenie_threads:**
```groovy
// OLD
regenie_threads = 24

// NEW
regenie_threads = 12  // Reduced to match available CPUs
```

**2. Add conda_env (after study parameter):**
```groovy
cohorts = ["Mayo"]
study = "Mayo"

// Add this:
conda_env = "test"
```

**3. Add process configuration at the end** (optional but recommended):
```groovy
// At the end of the file, after cohort_configs
process {
    withName: 'COMBINED_PIPELINE:GWAS_PIPELINE:REGENIE_STEP1' {
        cpus = 12
        memory = '48 GB'
    }
    
    withName: 'COMBINED_PIPELINE:GWAS_PIPELINE:REGENIE_STEP2' {
        cpus = 12
        memory = '48 GB'
    }
}
```

## Verification

After fixing all configs, verify with:

```bash
# Check all configs have regenie_threads = 12
grep -n "regenie_threads" nextflow.config.combined.*

# Should show:
# nextflow.config.combined.cmc_mssm:    regenie_threads = 12
# nextflow.config.combined.cmc_penn:    regenie_threads = 12
# ... etc
```

## Retry Pipeline

After fixing the config:

```bash
# Resume the pipeline
nextflow run combined_pipeline.nf \
    -c nextflow.config.combined.rosmap \
    -resume
```

The `-resume` flag will skip completed steps and continue from where it failed.

## For Systems with Different CPU Counts

| System CPUs | Set regenie_threads to |
|-------------|------------------------|
| 8 CPUs | 8 |
| 12 CPUs | 12 |
| 16 CPUs | 16 |
| 24+ CPUs | 24 |

**Rule**: Set `regenie_threads` to match (or be less than) your available CPUs.

## Check Your System's CPU Count

```bash
# Linux
nproc

# Or
lscpu | grep "^CPU(s):"

# Or check in your current session
echo "Available CPUs: $(nproc)"
```

## Performance Impact

**Regenie with 12 CPUs vs 24 CPUs**:
- **Speed**: ~1.5-2x slower per GWAS analysis
- **Total time**: Still acceptable (maybe add 4-8 hours per cohort)
- **Still much faster than sequential cohort execution**

## Alternative: Use SLURM

If your system has SLURM, you can request nodes with more CPUs:

```bash
#!/bin/bash
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G

nextflow run combined_pipeline.nf -c config -resume
```

Then you can keep `regenie_threads = 24` in the configs.

