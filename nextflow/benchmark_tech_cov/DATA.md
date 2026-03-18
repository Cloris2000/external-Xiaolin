# Technical covariate benchmark – raw data manifest

Use this list to obtain or point to the raw data needed to run the tech-covariate vs. MGP accuracy benchmark. Paths are the ones used on the lab server; your PI can replace them with local paths.

---

## Ground truth (snRNA-derived cell type proportions)

| Cohort     | Description                                      | Current path |
|-----------|---------------------------------------------------|--------------|
| CMC_MSSM  | PsychENCODE label-transferred snRNA proportions   | `/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/PsychEncode_label_transferred_snCTP.csv` |
| ROSMAP    | ROSMAP single-nucleus proportions                 | `/external/rprshnas01/netdata_kcni/stlab/cross_cohort_MGPs/rosmap_single_nuc_proportions.csv` |

---

## Bulk RNA-seq inputs (count matrix + metadata per cohort)

| Cohort    | File          | Current path |
|----------|---------------|--------------|
| CMC_MSSM | Count matrix  | `/nethome/kcni/xzhou/GWAS_tut/CMC_reimputed/CMC_MSSM_count_matrix.csv` |
| CMC_MSSM | Metadata      | `/nethome/kcni/xzhou/GWAS_tut/CMC_reimputed/CMC_MSSM_metadata.csv` |
| ROSMAP   | Count matrix  | `/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_DLPFC_batch_all.csv` |
| ROSMAP   | Metadata      | `/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_combined_metrics.csv` |

---

## Reference data (markers, HGNC)

| File          | Current path |
|---------------|--------------|
| MGP marker list | `/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/new_MTGnCgG_lfct2.5_Publication.csv` |
| HGNC mapping  | `/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/hgnc_complete_set.txt` |

---

## Shared tech-covariate lists (in this folder)

- `benchmark_configs/tech_cov_sets/shared_cmc_mssm.txt` — PMI, RIN, platform, libraryPrep, runType, readLength, isStranded
- `benchmark_configs/tech_cov_sets/shared_rosmap.txt` — pmi, RIN, platform, libraryPrep, runType, readLength, isStranded

No separate download needed.
