# 1. Pre-process metadata and count matrices (run once manually)
python3 scripts/preprocess_amp_ad_diverse.py

# 2. Per-sub-cohort phenotype + GWAS (can run in parallel after step 1)
nextflow run combined_pipeline_v2.nf -c nextflow.config.combined.amp_ad_columbia -w work/AMP_AD_Columbia -ansi-log false > logs/AMP_AD_Columbia/nf.log 2>&1 &

nextflow run combined_pipeline_v2.nf -c nextflow.config.combined.amp_ad_mssm -w work/AMP_AD_MtSinai -ansi-log false > logs/AMP_AD_MtSinai/nf.log 2>&1 &

nextflow run combined_pipeline_v2.nf -c nextflow.config.combined.amp_ad_rush -w work/AMP_AD_Rush -ansi-log false > logs/AMP_AD_Rush/nf.log 2>&1 &

tail -f logs/AMP_AD_Columbia/nf.log logs/AMP_AD_MtSinai/nf.log logs/AMP_AD_Rush/nf.log

# 3. All-cohorts (15-cohort) meta-analysis
nextflow run standalone_meta_analysis.nf -c nextflow.config.standalone_meta_all_cohorts -resume