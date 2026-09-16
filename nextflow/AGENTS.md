# Agent rules

## SCC directory is read-only

`/project/rrg-shreejoy/zhoux156/Xiaolin/SCC/` is a transferred snapshot. **Do not modify it.**

| Allowed | Not allowed |
|---|---|
| Read `Xiaolin/SCC/nextflow/{data,data_input,results}` | Write, edit, delete, or overwrite anything under `Xiaolin/SCC/` |
| Edit this repo: `external-Xiaolin/nextflow/` | Save figures or intermediates into `Xiaolin/SCC/nextflow/manuscript_figure/` |

**Write outputs here:** `/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/`

**Figure 2 outputs:** `/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/manuscript_figure/`

If a SLURM compute node cannot write to `/project`, write under `/scratch/zhoux156/` and copy the results into `external-Xiaolin/nextflow/`. Never copy results into SCC.
