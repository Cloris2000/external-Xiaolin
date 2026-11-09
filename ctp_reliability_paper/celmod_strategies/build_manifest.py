"""Build strategy_pairs.tsv: Shreejoy's retention manifest + the new pairs and recipes.

Existing rows are copied unchanged from
  /scratch/shreejoy/composition/wt_gtex/manifests/retention_pairs.tsv
(the manifest that built results_2026-10-21c/retention, git 134bb7b), so every
existing trainer/target is fitted exactly as before. LBP and ROSMAP_FC/TC rows are
kept although they score nothing: they are part of the common gene space, and
dropping them would change every existing fit.

Added (counts from preflight/pairs_matched.tsv, release 2026-10-21c):
  ROSMAP x CUIMC2                       236  both
  ROSMAP_PCC x CUIMC2                   164  both
  NIMH_HBCC x Catching_2026_HBCC         89  both  (needs MIN_FIT 89; decision 2026-11-04)
  LIBD_szControl_HBCC x Catching_2026_HBCC 63 target
  NABEC x Catching_2026_NABEC           198  both
  ROSMAP consensus Green+Mathys+Multiome2025_DLPFC+CUIMC2 and its leave-one-out
    variants (the G+M+Mu one is the existing production trio, already declared)
  NIMH_HBCC consensus PsychAD_HBCC+Catching_2026_HBCC

Adds a `class` column for strategy S3, by DIAGNOSIS (decision 2026-11-04), from the
bulk cohort's resolved persons (registry covariates.tsv dx, dx_vocabulary family):
  adult_psych  if >= 20% psychiatric-family dx, or (< 5% dementia-family dx and
               median age < 65)   [HBCC, LIBD_HBCC, CommonMind, NABEC]
  aging_AD     if >= 20% dementia-family dx  [ROSMAP*, MSBB*, Diverse*]
  other        >= 20% movement-family dx (LBP; scores nothing), or none of the above
Region/lobe of the new snRNA rows come from the registry datasets.tsv (CUIMC2 and
Catching_2026 are DLPFC); bulk fields are copied from the cohort's existing rows,
except NABEC (registry bulk_NABEC region FCTX, frontal).
"""
import os
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
SRC = Path("/scratch/shreejoy/composition/wt_gtex/manifests/retention_pairs.tsv")
REL = Path("/project/rrg-shreejoy/brain-data-registry/releases/2026-10-21c")
PF = HERE / "preflight"

man = pd.read_csv(SRC, sep="\t", dtype=str, keep_default_na=False)
cols = list(man.columns)
matched = pd.read_csv(PF / "pairs_matched.tsv", sep="\t")
nm = {(r.bulk, r.snrna): int(r.n_matched) for r in matched.itertuples()
      if pd.notna(r.n_matched)}


def bulk_fields(cohort):
    r = man[(man.bulk_cohort == cohort) & (man.role != "excluded")].iloc[0]
    return dict(family=r.family, bulk_region=r.bulk_region, bulk_lobe=r.bulk_lobe)


def row(bulk, snrna, role, study, note, snrna_region="DLPFC", snrna_lobe="frontal",
        bf=None):
    bf = bf or bulk_fields(bulk)
    n = "" if "+" in snrna else str(nm[(bulk, snrna)])
    return dict(bulk_cohort=bulk, snrna=snrna, role=role, **bf,
                snrna_region=snrna_region, snrna_lobe=snrna_lobe, study=study,
                n_matched=n, note=note)


TAG = "Added 2026-11-04 (ctp_reliability_paper/celmod_strategies)."
new = [
    row("ROSMAP", "CUIMC2", "both", "CUIMC2_ROSMAP", f"{TAG} CUIMC2 joined the snRNA pairs 2026-10-29."),
    row("ROSMAP_PCC", "CUIMC2", "both", "CUIMC2_ROSMAP", TAG),
    row("NIMH_HBCC", "Catching_2026_HBCC", "both", "Catching_2026",
        f"{TAG} 89 matched; trains only with MIN_FIT lowered to 89 (decision 2026-11-04)."),
    row("LIBD_szControl_HBCC", "Catching_2026_HBCC", "target", "Catching_2026", TAG),
    row("NABEC", "Catching_2026_NABEC", "both", "Catching_2026", TAG,
        bf=dict(family="NABEC", bulk_region="FCTX", bulk_lobe="frontal")),
]
FULL = ["Green", "Mathys", "Multiome2025_DLPFC", "CUIMC2"]
recipes = ["+".join(FULL)] + ["+".join(m for m in FULL if m != drop) for drop in FULL]
# leave-one-out variants of the production trio (Green+Multiome2025_DLPFC already declared)
recipes += ["Green+Mathys", "Mathys+Multiome2025_DLPFC"]
have = set(man.snrna)
for rc in recipes:
    if rc in have:  # G+M+Mu = production trio, already declared
        continue
    new.append(row("ROSMAP", rc, "trainer", "consensus",
                   f"{TAG} DLPFC-set consensus / leave-one-out variant, so every "
                   f"ROSMAP scorer has a recipe that excludes it.", snrna_region="PFC"))
new.append(row("NIMH_HBCC", "PsychAD_HBCC+Catching_2026_HBCC", "trainer", "consensus",
               f"{TAG} HBCC consensus.", snrna_region="PFC"))

out = pd.concat([man, pd.DataFrame(new)[cols]], ignore_index=True)
dup = out.duplicated(["bulk_cohort", "snrna"])
assert not dup.any(), out[dup]

# ---- class by diagnosis --------------------------------------------------------
cc = pd.read_csv(PF / "cohort_class.tsv", sep="\t")  # median_age per cohort
vocab = pd.read_csv(REL / "dx_vocabulary.tsv", sep="\t").set_index("dx").family
cov = pd.read_csv(REL / "covariates.tsv", sep="\t", usecols=["person_id", "dx"],
                  dtype=str).set_index("person_id")
# persons per bulk cohort: re-derived by preflight; stored there only as counts, so
# recompute fractions from the preflight person sets written below
persons = pd.read_csv(PF / "bulk_persons.tsv", sep="\t", dtype=str)
fam = persons.join(cov, on="person_id").assign(family=lambda d: d.dx.map(vocab))
frac = (fam.dropna(subset=["family"]).groupby("cohort").family
        .value_counts(normalize=True).unstack(fill_value=0))
cls = {}
for c in sorted(set(out.bulk_cohort)):
    f = frac.loc[c] if c in frac.index else pd.Series(dtype=float)
    psy, dem = f.get("psychiatric", 0.0), f.get("dementia", 0.0)
    mov = f.get("movement", 0.0)
    age = cc.set_index("cohort").median_age.get(c, float("nan"))
    if mov >= 0.20:  # LBP (Parkinson's DBS cohorts): neither class
        cls[c] = "other"
    elif psy >= 0.20 or (dem < 0.05 and age < 65):
        cls[c] = "adult_psych"
    elif dem >= 0.20:
        cls[c] = "aging_AD"
    else:
        cls[c] = "other"
    print(f"{c:22s} psych {psy:5.2f} dementia {dem:5.2f} movement {mov:5.2f} "
          f"median_age {age:6.1f} -> {cls[c]}")
out["class"] = out.bulk_cohort.map(cls)
out.to_csv(HERE / "strategy_pairs.tsv", sep="\t", index=False)
print(f"\nwrote {HERE/'strategy_pairs.tsv'}: {len(man)} existing + {len(new)} new rows")
print(out.loc[len(man):, ["bulk_cohort", "snrna", "role", "n_matched", "class"]].to_string(index=False))
