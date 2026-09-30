# Sample overlap across the 15-cohort meta-analysis

Checks whether the same individual entered the meta-analysis under more than one
cohort. Built 2026-09-17.

## What "a cohort's samples" means here

The post-QC genotype sample list actually analyzed by REGENIE:

    results/{cohort}/*.QC.final.psam

3,620 samples across the 15 cohorts in `nextflow.config.standalone_meta_15cohorts`.
This is deliberately the `.psam`, not an upstream manifest, so the answer reflects
who survived QC and entered the GWAS.

## Evidence used

1. **Genotype identity** — `/project/rrg-shreejoy/Public_datasets/_registry/genotype_matches.csv`,
   pairwise matches over 5,401 common autosomal SNPs. Identity called at
   `ibs0_rate < 0.005` (the registry README documents 0.005–0.015 as empty and
   0.02–0.03 as first-degree relatives, so the threshold is unambiguous).
   Samples are unioned into connected components = persons.

2. **Identifier crosswalk (ROSMAP family only)** — the registry never
   systematically compared ROSMAP WGS against ROSMAP array, so that pair was
   resolved by ID instead:

       ROSMAP        SM-*          -> biospecimen metadata      -> individualID (R*)
       ROSMAP_array  MAP*/ROS*     -> samples_reheader_map.txt  -> 11AD* specimenID
                                   -> biospecimen metadata      -> individualID
                     MAP*/ROS*     -> projid                    -> individualID
       AMP_AD_Rush   R*            -> already individualID

   All 796 + 162 + 133 IDs resolve (100%).

## Results

134 distinct individuals appear in 2+ cohorts; 145 redundant samples
(~4.0% of 3,620).

| cohort A | cohort B | shared |
|---|---|---|
| GVEX | Mayo | 46 |
| GVEX | NABEC | 40 |
| CMC_MSSM | MSBB | 38 |
| Mayo | NABEC | 15 |
| AMP_AD_Rush | ROSMAP_array | 6 |
| NABEC | NIMH_HBCC_h650 | 5 |
| NABEC | NIMH_HBCC_1M | 3 |
| AMP_AD_Mayo | GVEX | 2 |
| AMP_AD_Rush | ROSMAP | 1 |

Per-individual detail: `FINAL_overlapping_individuals.tsv`
(`dup_id`, `n_cohorts`, `cohort`, `analyzed_sample_id`).

GVEX/Mayo/NABEC form 3-way groups — the Banner Sun Health donor pool. The most
exposed cohorts are NABEC (24.8%), Mayo (19.5%), GVEX (19.5%), CMC_MSSM (15.7%)
and MSBB (15.1%).

## Verified independent

- **The 3 NIMH_HBCC platform cohorts do not share individuals.** The registry has
  zero HBCC-vs-HBCC edges, so this was untested there. Resolved indirectly
  (`hbcc_bridge.py`): 290 HBCC samples match external PsychAD/Multiome WGS
  anchors, and every anchor maps to exactly one HBCC sample — no anchor pulls in
  two platforms. 257 distinct samples, matching the known HBCC↔PsychAD figure.
- **ROSMAP vs ROSMAP_array = 0 overlap**, by two independent methods (ID
  crosswalk above, and 26 ROSMAP-internal genotype rows, none of which pair an
  analyzed WGS sample with an analyzed array sample). Consistent with
  ROSMAP_array being the TOPMed-imputed samples *not* in the WGS pipeline.
- **The 3 CMC sites are mutually disjoint** at person level (separate brain banks).
- **Mayo vs AMP_AD_Mayo = 0, by ancestry construction.** These share the Mayo/Emory
  site name but are separate studies: Mayo (AMP-AD 1.0) is 257/257 EUR, while
  AMP-AD Diverse: Mayo/Emory was recruited specifically for ancestry diversity
  (52 AFR + 178 AMR, zero EUR). The strata do not intersect, so no individual can
  be in both. They also differ in build (b37 vs GRCh38) and use different WGS
  pools (349 vs the shared 743-sample set). Note this is the *structural* reason;
  the fact that their numeric IDs do not collide is NOT itself evidence, since
  unrelated numbering systems would not collide even for a shared donor.

## Caveats — cohorts that could NOT be fully checked

Absence of evidence is not evidence of absence: `genotype_matches.csv` only
contains pairs that *matched*, so a cohort missing from it was either never
compared or compared and matched nothing. Only 9 of 105 cohort pairs have direct
genotype evidence.

- **CMC_PITT — unresolvable.** No representation in the genotype panel, and its
  `.psam` IDs (`0_PITT_516`+) use different numbering from the registry
  (`CMC_PITT_014`–`170`), so 0/161 map. Cannot be checked by genotype or by ID.
- **CMC_PENN — genotype-blind.** No panel representation; 80/93 map to registry
  persons by ID, and none collide with another GWAS cohort. Lower risk than PITT
  but not genotype-verified.
- **GTEx_v10 — genotype-blind.** No panel representation; all 312 map to registry
  persons, none shared with another GWAS cohort. GTEx donors are independently
  recruited, so prior risk is low.
- **AMP_AD_Mayo (0.9% in panel) and ROSMAP (2.5%)** have thin genotype coverage.
  ROSMAP is covered by the ID crosswalk. AMP_AD_Mayo is not, but is resolved by
  the ancestry argument above rather than by genotype.

### Ancestry context for the remaining gap

Ancestry strata bound where overlap is even possible
(`docs/ancestry_specific/ancestry_sample_assignments.tsv`):

| cohort | AFR | AMR | EAS | EUR | OTHER | UNKNOWN |
|---|---|---|---|---|---|---|
| Mayo | 0 | 0 | 0 | 257 | 0 | 0 |
| AMP_AD_Mayo | 52 | 178 | 0 | 0 | 0 | 0 |
| CMC_MSSM | 10 | 5 | 0 | 230 | 0 | 0 |
| CMC_PENN | 13 | 0 | 0 | 80 | 1 | 0 |
| CMC_PITT | 27 | 0 | 0 | 139 | 0 | 0 |
| MSBB | 28 | 17 | 1 | 205 | 1 | 0 |

This does **not** relieve the CMC_PITT gap: PITT is 139 EUR / 27 AFR, overlapping
the same strata as CMC_MSSM (230 EUR) and MSBB (205 EUR), which are exactly the
cohorts already shown to share 38 individuals. CMC_PITT sits in the highest-risk
stratum and remains the one genuinely unchecked cohort.

## Related check: ROSMAP genotype↔phenotype ID join (verified OK)

`nextflow.config.combined.rosmap` sets `fid_method = 'biospec_specimen'` with
`biospec_assay_filter = 'wholeGenomeSeq'`, so `pheno_prep.R` rewrites phenotype
IDs into the WGS `SM-*` specimenID space. Verified on the actual output
(`idjoin.py`): `phenotypes_RINT.txt` has 852 rows — 817 `SM-*` and 35 `R*` — and
all 796 analyzed genotype samples have a matching phenotype row. Nothing is
silently dropped; the 35 `R*` rows are donors with no WGS genotype and are
correctly excluded by REGENIE. The crosswalk is wired and working.

## Reproducing

    python3 final_report.py    # headline numbers + FINAL_overlapping_individuals.tsv
    python3 coverage.py        # per-cohort panel representation / blind spots
    python3 hbcc_bridge.py     # HBCC cross-platform test
    python3 rosmap_xwalk2.py   # ROSMAP-family ID crosswalk

Read-only; nothing writes into the results tree.
