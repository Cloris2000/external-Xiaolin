# Combined Pipeline - Quick Reference Card

## One-Line Commands

### Run Single Cohort
```bash
nextflow run combined_pipeline.nf -c nextflow.config.combined.cmc_mssm -resume
```

### Run All Cohorts
```bash
./run_all_cohorts_combined.sh
```

### Run Specific Cohorts
```bash
./run_all_cohorts_combined.sh ROSMAP Mayo MSBB
```

## Pipeline Flow
```
GENOTYPING_QC_PIPELINE
         ↓
    (generates .psam + pca.csv)
         ↓
PHENOTYPE_PIPELINE
         ↓
    (uses .psam + pca.csv)
         ↓
GWAS_PIPELINE
         ↓
    Results
    
All stages run SEQUENTIALLY (not parallel)
```

## Files You Need to Create

For each cohort, create: `nextflow.config.combined.cohort_name`

Copy from template:
```bash
cp nextflow.config.combined.cmc_mssm nextflow.config.combined.my_cohort
```

## Required Parameters

```groovy
params {
    // General
    cohorts = ["MY_COHORT"]
    study = "MY_COHORT"
    output_dir = "${projectDir}/results/MY_COHORT"
    
    // Phenotype Pipeline Inputs
    count_matrix_file = "/path/to/expression.tsv"
    metadata_file = "/path/to/metadata.csv"
    tissue_filter = "DLPFC"
    col_sample_id_for_matching = "individualID"
    
    // Genotyping QC Inputs
    vcf_dir = "/path/to/vcf"
    vcf_pattern = "*.vcf.gz"
    plink_path = "/path/to/plink2"
    bcftools_path = "/path/to/bcftools"
    
    // GWAS Configuration
    regenie_path = "/path/to/regenie"
    metal_path = "/path/to/metal"
    
    cohort_configs = [
        "MY_COHORT": [
            pgen_file: "${output_dir}/MY_COHORT.QC.final",
            prune_in_file: "${output_dir}/MY_COHORT.QC.prune.in",
            pheno_file: "${output_dir}/phenotypes_RINT.txt",
            covar_file: "${output_dir}/covariates.txt"
        ]
    ]
}
```

## Output Files

```
results/MY_COHORT/
├── phenotypes_RINT.txt              # Phenotypes for 19 cell types
├── covariates.txt                   # Covariates (age, sex, PCs, etc.)
├── MY_COHORT.QC.final.{pgen,psam,pvar}  # QC'd genotypes
├── pca.csv                          # Genetic PCs
├── regenie_step1/                   # GWAS step 1 (19 files)
├── regenie_step2/                   # GWAS step 2 (19 files)
└── combined_pipeline_report.html    # Execution report
```

## Troubleshooting

| Error | Solution |
|-------|----------|
| File not found | Check paths in config file |
| Out of memory | Increase memory in process config |
| GWAS starts too early | Should not happen - check you're using combined_pipeline.nf |
| Resume doesn't work | Try `nextflow clean -f` then rerun |

## Monitoring

```bash
# View running pipeline
tail -f .nextflow.log

# View reports (after completion)
firefox results/MY_COHORT/combined_pipeline_report.html
firefox results/MY_COHORT/combined_pipeline_timeline.html

# Check specific process logs
ls -lh work/*/*/.command.log
```

## Runtime Estimates

| Stage | Time |
|-------|------|
| Phenotype Pipeline | 2-6h |
| Genotyping QC | 4-12h |
| GWAS Pipeline | 8-24h |
| **Total per cohort** | **14-42h** |

## Common Options

```bash
# Resume from last successful step
-resume

# Use specific config
-c nextflow.config.combined.my_cohort

# Generate reports
-with-report report.html
-with-timeline timeline.html
-with-dag flowchart.png

# Dry run (show what would run)
-preview
```

## Cell Types Analyzed

19 cell types per cohort:
- Astrocyte, Endothelial, IT, L4.IT, L5.ET
- L5.6.IT.Car3, L5.6.NP, L6.CT, L6b, LAMP5
- Microglia, OPC, Oligodendrocyte, PAX6, PVALB
- Pericyte, SST, VIP, VLMC

## Documentation Files

| File | Purpose |
|------|---------|
| `COMBINED_PIPELINE_SUMMARY.md` | Overview and quick start |
| `COMBINED_PIPELINE_SETUP_GUIDE.md` | Detailed setup instructions |
| `COMBINED_PIPELINE_README.md` | Complete user guide |
| `COMBINED_PIPELINE_QUICKREF.md` | This reference card |

## Support Files

| File | Purpose |
|------|---------|
| `combined_pipeline.nf` | Main pipeline |
| `nextflow.config.combined` | Base config template |
| `nextflow.config.combined.cmc_mssm` | Example config |
| `nextflow.config.combined.rosmap` | Example config |
| `run_all_cohorts_combined.sh` | Run all cohorts script |

## Dependencies

The combined pipeline uses these sub-workflows:
- `phenotype_pipeline.nf` → PHENOTYPE_PIPELINE
- `genotyping_qc_pipeline.nf` → GENOTYPING_QC_PIPELINE  
- `gwas_pipeline.nf` → GWAS_PIPELINE

All original pipelines still work independently!

