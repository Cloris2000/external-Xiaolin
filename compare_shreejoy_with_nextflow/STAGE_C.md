# Stage C — deferred genome-wide job

Stage A and B are done. The six figures and `SUMMARY.md` cover the four
analysis stages at two loci. This file is the review of what is still missing
and the job that would fill it. **It has not been submitted.** Per
`CLAUDE.md`, nothing goes to Slurm until you say so.

## What Stage B already answers

- Stage 1: deconvolution agreement, and that 520 of 687 hits come from traits
  with no portable counterpart.
- Stage 2: TMEM106B cancellation / glial-shadow / SST attenuation.
- Stage 3: 77% high-I2 on my side; his own ancestry meta recovers the same two
  loci as his mega.
- Stage 4: coloc agrees at TMEM106B when the disease file is shared; GRN is
  missing because stage 2 never found it.
- GRN effect estimates: same sign, same SE, ~4x smaller on my SST.

## What it cannot answer

1. **My genomic-control lambda per trait.** His is median 1.0026, 0 of 101
   outside [0.95, 1.05] (`REPORT.md`). 687 hits at 77% heterogeneity invites
   the calibration question. Nothing in Stage B answers it.
2. **Genome-wide z-score concordance** between my 19 classes and his mapped
   supertypes. Everything above is at TMEM106B and GRN.

Both need one pass over the 19 `*.annotated.tsv` files (~18 GB) plus his
per-chromosome `.regenie` files. That is a compute-node job, not a login-node
one.

## Proposed job (shown, not submitted)

A single-node script that, per trait, streams the annotated meta file and
writes:

- `data/my_lambda.tsv` — genomic-control lambda from median chi-square
  (`lambda = median(qchisq(1-P, 1), na.rm) / qchisq(0.5, 1)`), plus n
  variants and n with P < 5e-8. Completes fig2 panel (c) against his
  `lambda/lambda_by_trait.csv`.
- `data/z_concordance.tsv` — for each (my class, his mapped supertype) pair,
  Pearson r of z-scores on variants matched by chr:pos:REF:ALT after applying
  the two empirical offsets. Restricted to autosomes already in both files.

Sketch (not a submission):

```bash
# login node: inspect only
# compute node, if you approve:
#   8 hours, 1 node, no --mem (Trillium linear select)
#   RSCRIPT or python3 streaming awk + python
#   write only under compare_shreejoy_with_nextflow/data/
```

I will write `07_stage_c.sbatch` and show it to you if you want this run.
I will not submit it without that go-ahead.

## Recommendation

Worth doing. The lambda is the one number that would change how 687 hits are
read, and it is cheap relative to a full re-analysis. Concordance is secondary
but would show whether TMEM106B/GRN are typical or special.
