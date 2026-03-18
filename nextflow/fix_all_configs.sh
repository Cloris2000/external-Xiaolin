#!/bin/bash

# Fix software paths and CPU settings in all cohort configs
# Only bcftools path needs to be changed (user's location)
# Also reduce CPU requirements from 24 to 12

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

echo "Fixing all cohort configs..."
echo ""

for cohort in mayo msbb gtex nabec cmc_mssm cmc_penn cmc_pitt; do
    config_file="nextflow.config.combined.${cohort}"
    
    if [ ! -f "$config_file" ]; then
        echo "❌ Config not found: $config_file"
        continue
    fi
    
    echo "Fixing $config_file..."
    
    # 1. Fix bcftools path (user's location)
    sed -i 's|bcftools_path = "/external/rprshnas01/kcni/mwainberg/software/bcftools"|bcftools_path = "/nethome/kcni/xzhou/software/bcftools"|g' "$config_file"
    
    # 2. Fix CPU requirements (24 -> 12)
    sed -i 's/regenie_threads = 24/regenie_threads = 12  \/\/ Reduced to match available CPUs/' "$config_file"
    
    # 3. Add conda_env if missing (to avoid warning)
    if ! grep -q "conda_env" "$config_file"; then
        sed -i '/study = /a\    \n    // Conda environment (optional)\n    conda_env = "test"' "$config_file"
    fi
    
    echo "  ✅ Fixed bcftools path"
    echo "  ✅ Fixed CPU requirements (24 → 12)"
    echo "  ✅ Added conda_env parameter"
    echo ""
done

echo "All configs updated!"
echo ""
echo "Summary of software paths:"
echo "  plink2:   /external/rprshnas01/kcni/mwainberg/software/plink2"
echo "  bcftools: /nethome/kcni/xzhou/software/bcftools"
echo "  regenie:  /external/rprshnas01/kcni/mwainberg/software/regenie"
echo "  METAL:    /external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/METAL/generic-metal/executables"
echo ""
echo "Verification:"
grep "bcftools_path" nextflow.config.combined.* | head -8
echo ""
grep "regenie_threads" nextflow.config.combined.* | head -8

