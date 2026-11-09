"""Pre-flight for the CelMod training-strategy comparison (final figure).

Lightweight, login-node safe: reads only parquet indexes/schemas, the snRNA count
tables and registry TSVs. Fits nothing.

Reports, for the plan at ~/.claude/plans/hi-for-the-ctp-dapper-bumblebee.md step 0:
  1. matched people for every candidate (bulk, snRNA) pair, existing and new, using
     Shreejoy's own id translation (ctcomp.bulk.to_persons, release 2026-10-21c)
  2. within-bulk overlap between snRNA datasets (what the person hold-out removes
     when one trains and the other scores)
  3. the common gene space with and without the new cohorts
  4. per bulk cohort: median age and dx mix of its resolved persons -> S3 class
     (rule fixed before any scoring: median age >= 65 -> aging_AD, else adult_psych)

Run (login node is fine):
  export CTCOMP_REPO=/scratch/shreejoy/composition/wt_gtex REGISTRY_RELEASE=2026-10-21c
  source $CTCOMP_REPO/env.sh; unset LABCORE_READS
  $PY preflight.py
Writes: celmod_strategies/preflight/*.tsv (next to this script).
"""
import os
import sys
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

REPO = Path(os.environ["CTCOMP_REPO"])
sys.path.insert(0, str(REPO / "src"))
from ctcomp.bulk import to_persons, matrix_stem, bank_of, PART_OF  # noqa: E402

RES = Path(os.environ["CTCOMP_RESULTS"])
COUNTS = RES / "counts"
BULK = Path(os.environ["CTCOMP_BULK"])
REL = Path(os.environ["REGISTRY_RELEASES"]) / os.environ["REGISTRY_RELEASE"]
OUT = Path(__file__).resolve().parent / "preflight"
OUT.mkdir(exist_ok=True)
LEVEL = "supertype"
AGE_CUT = 65

MAN = pd.read_csv(REPO / "manifests" / "retention_pairs.tsv", sep="\t", dtype=str,
                  keep_default_na=False)
MAN = MAN[(MAN.role != "excluded") & ~MAN.snrna.str.contains(r"\+")]
NEW = [("ROSMAP", "CUIMC2"), ("ROSMAP_PCC", "CUIMC2"),
       ("NIMH_HBCC", "Catching_2026_HBCC"), ("LIBD_szControl_HBCC", "Catching_2026_HBCC"),
       ("NABEC", "Catching_2026_NABEC")]
PAIRS = [(r.bulk_cohort, r.snrna, "existing") for r in MAN.itertuples()] + \
        [(b, s, "new") for b, s in NEW]
NEW_COHORTS = {"NABEC", "LIBD_szControl_HBCC"} - set(MAN.bulk_cohort)


def bulk_path(c, v):
    return BULK / f"{matrix_stem(c)}_{v}.parquet"


_persons = {}


def bulk_persons(c):
    """person ids of a bulk cohort's raw matrix, via Shreejoy's translation.
    Bank-split stems (CommonMind_MSSM) are filtered the way read_target does."""
    if c not in _persons:
        idx = pd.read_parquet(bulk_path(c, "raw"), columns=[]).index.astype(str)
        p, rep = to_persons(idx, matrix_stem(c))
        s = set(p.dropna())
        if c in PART_OF:  # one bank of a split matrix, as read_target filters it
            banks = bank_of(matrix_stem(c))
            s = {x for x in s if banks.get(x) == PART_OF[c][1]}
        _persons[c] = s
        print(f"  bulk {c}: {rep['samples']} samples -> {len(s)} persons", flush=True)
    return _persons[c]


_sn = {}


def sn_persons(s):
    if s not in _sn:
        _sn[s] = set(pd.read_parquet(COUNTS / f"{s}__{LEVEL}.parquet", columns=[]).index.astype(str))
    return _sn[s]


def main():
    # 1. matched people per pair
    rows = []
    for b, s, kind in PAIRS:
        try:
            n = len(bulk_persons(b) & sn_persons(s))
            err = ""
        except Exception as e:  # report, do not guess
            n, err = None, f"{type(e).__name__}: {e}"
        rows.append(dict(bulk=b, snrna=s, kind=kind, n_matched=n,
                         n_bulk=len(_persons.get(b, ())), n_snrna=len(_sn.get(s, ())),
                         role_by_floor=("both" if n and n >= 90 else
                                        "target" if n and n >= 50 else "none"),
                         error=err))
    P = pd.DataFrame(rows)
    P.to_csv(OUT / "pairs_matched.tsv", sep="\t", index=False)
    print("\n== matched people per pair (fit floor 90, score floor 50) ==")
    print(P.drop(columns="error").to_string(index=False))
    if P.error.any():
        print(P[P.error != ""][["bulk", "snrna", "error"]].to_string(index=False))

    # 2. within-bulk overlap of snRNA datasets (the hold-out casualties)
    ov = []
    for b, g in P[P.n_matched.fillna(0) > 0].groupby("bulk"):
        ss = list(g.snrna)
        for i, s1 in enumerate(ss):
            for s2 in ss[i + 1:]:
                a = bulk_persons(b) & sn_persons(s1)
                c = bulk_persons(b) & sn_persons(s2)
                ov.append(dict(bulk=b, a=s1, b=s2, n_a=len(a), n_b=len(c),
                               shared=len(a & c), b_only=len(c - a), a_only=len(a - c)))
    O = pd.DataFrame(ov)
    O.to_csv(OUT / "within_bulk_overlap.tsv", sep="\t", index=False)
    print("\n== within-bulk snRNA overlap (a trains -> b scored on b_only; and reverse) ==")
    print(O.to_string(index=False))

    # 3. common gene space
    cohorts = sorted(set(P.bulk))
    genes = {}
    for c in cohorts:
        sets = [set(pq.read_schema(bulk_path(c, v)).names) for v in ("light", "raw")]
        genes[c] = set.intersection(*sets)
    old = set.intersection(*[genes[c] for c in cohorts if c not in NEW_COHORTS])
    allg = set.intersection(*genes.values())
    print(f"\n== gene space == existing cohorts: {len(old)}; with {sorted(NEW_COHORTS)}: {len(allg)}")
    for c in sorted(NEW_COHORTS):
        print(f"  {c}: own {len(genes[c])}, shared with existing {len(old & genes[c])}")
    pd.DataFrame([dict(cohort=c, n_genes=len(genes[c]), in_existing_space=len(old & genes[c]))
                  for c in cohorts]).to_csv(OUT / "gene_space.tsv", sep="\t", index=False)

    # 4. cohort class from registry covariates
    cov = pd.read_csv(REL / "covariates.tsv", sep="\t", dtype={"person_id": str},
                      usecols=["person_id", "age", "dx"], low_memory=False).set_index("person_id")
    cls = []
    for c in cohorts:
        ps = sorted(bulk_persons(c))
        cv = cov.reindex(ps)
        dx = cv.dx.fillna("NA").value_counts()
        med = cv.age.median()
        cls.append(dict(cohort=c, n_persons=len(ps), n_age=int(cv.age.notna().sum()),
                        median_age=med, q25=cv.age.quantile(.25), q75=cv.age.quantile(.75),
                        cls=("aging_AD" if med >= AGE_CUT else "adult_psych"),
                        dx_top=", ".join(f"{k} {v}" for k, v in dx.head(5).items())))
    C = pd.DataFrame(cls)
    C.to_csv(OUT / "cohort_class.tsv", sep="\t", index=False)
    pd.DataFrame([dict(cohort=c, person_id=p) for c in cohorts for p in sorted(bulk_persons(c))]
                 ).to_csv(OUT / "bulk_persons.tsv", sep="\t", index=False)
    print(f"\n== cohort class (median age >= {AGE_CUT} -> aging_AD) ==")
    print(C.to_string(index=False))

    # snRNA-side class too (for reference; S3 classes trainers by their bulk cohort)
    srow = []
    for s in sorted(set(P.snrna)):
        cv = cov.reindex(sorted(sn_persons(s)))
        srow.append(dict(snrna=s, n=len(cv), median_age=cv.age.median(),
                         dx_top=", ".join(f"{k} {v}" for k, v in
                                          cv.dx.fillna("NA").value_counts().head(4).items())))
    S = pd.DataFrame(srow)
    S.to_csv(OUT / "snrna_age_dx.tsv", sep="\t", index=False)
    print("\n== snRNA datasets: age / dx ==")
    print(S.to_string(index=False))


if __name__ == "__main__":
    main()
