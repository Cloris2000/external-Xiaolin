"""Assign CelMod training strategies to evaluation units and score them on common donors.

Input: the trainer x target matrix written by run_matrix.py ($STRATEGY_OUT, default
/scratch/zhoux156/ctp_reliability_paper/celmod_strategies/matrix): scores.csv,
pred/<trainer>__<bulk>.parquet, oof/<target>.parquet, plus strategy_pairs.tsv.
Output: ../figure9/data/fig9_*.tsv. Lightweight (no fitting); login node is fine.

Evaluation unit = (target bulk b, scorer snRNA s). Ground truth is CLR(s) of b's own
people (decision 2026-11-04):
  * cross  - a same-bulk trainer that does NOT contain s is scorable on b x s after the
             person hold-out (>= 50 donors). If b has any cross unit, only cross units
             of b are used.
  * oof    - otherwise (b has one snRNA, or every cross pairing is below the floor):
             one unit per scorable s, S1 = the same-dataset out-of-fold fit (flagged).

Strategies (fixed before scoring, plan step 3):
  S1 paired   same-bulk trainers excluding s. The declared one with the most members
              (a consensus) if any has > 1 member, else the per-cell-type MEAN r over
              the eligible single-dataset trainers (no selection on the target).
              oof units: the out-of-fold diagonal.
  S2 global   one recipe for every unit, chosen leave-one-cohort-family-out: best mean
              (over target bulk cohorts in the other families) of the unit median Core r,
              among candidates covering >= 80% of those cohorts. A recipe containing s
              is replaced by its declared leave-s-out variant; none -> undefined.
  S3 matched  as S2, candidates and selection units restricted to the unit's class.
  prod        current production trio ROSMAP Green+Mathys+Multiome2025_DLPFC (variant rule).
  oracle      best eligible trainer per unit, from scores.csv (optimistic; not for claims).
Selection reads scores.csv (each trainer's own hold-out donors). The comparison is then
RECOMPUTED on one donor set per unit: s's matched donors in b minus every person any of
the unit's S1/S2/S3/prod trainers trained on. Per-cell-type Pearson r on that set.

Run:
  export CTCOMP_REPO=/scratch/shreejoy/composition/wt_gtex REGISTRY_RELEASE=2026-10-21c
  source $CTCOMP_REPO/env.sh; unset LABCORE_READS; $PY strategy_scores.py
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
MAT = Path(os.environ.get("STRATEGY_OUT",
                          "/scratch/zhoux156/ctp_reliability_paper/celmod_strategies/matrix"))
FIG = HERE.parent / "figure9" / "data"
FIG.mkdir(parents=True, exist_ok=True)
COMMON = HERE.parent / "common" / "data"
REPO = Path(os.environ["CTCOMP_REPO"])
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "stages" / "composition"))
from ctcomp.bulk import to_persons, matrix_stem, bank_of, PART_OF  # noqa: E402
import retention_matrix as RM  # noqa: E402  (read-only use: load_counts)

MIN_SCORE = 50
COVER = 0.80
# METRIC=spearman (user 2026-11-10): selection reads scores_spearman.csv (scores_spearman.py;
# same rows/donors as scores.csv, Spearman rho) and the common-donor step uses Spearman.
# Outputs get a "_spearman" suffix; the default (pearson) run and its files are unchanged.
METRIC = os.environ.get("METRIC", "pearson")
assert METRIC in ("pearson", "spearman"), METRIC
SUF = "" if METRIC == "pearson" else "_spearman"


def corr(a, y):
    return np.corrcoef(a, y)[0, 1] if METRIC == "pearson" else stats.spearmanr(a, y).correlation
PROD = "ROSMAP__Green+Mathys+Multiome2025_DLPFC"
# The paper's snRNA datasets (Figs 1-4, 8 pairs, plus the matched-bulk Catching_NABEC).
# Multiome2025_STG is excluded as scorer and trainer member: a second region of the
# Multiome2025_DLPFC donors (259/273 shared in ROSMAP), not in the paper's pairs, and
# as a scorer the person hold-out leaves ~6 donors (dry test 2026-11-04).
PAPER_SN = {"Green", "Mathys", "Multiome2025_DLPFC", "CUIMC2", "PsychAD_MSBB",
            "PsychAD_HBCC", "PsychAD_ROSMAP", "BrainSCOPE_CMC", "Catching_2026_HBCC",
            "Catching_2026_NABEC"}
# One bulk matrix per cohort, the region closest to DLPFC (PI request 2026-11-04):
# ROSMAP DLPFC (not PCC); MSBB FP = BM10 granular prefrontal, adjacent to DLPFC BA9/46
# (not IFG BM44 dysgranular / STG temporal; chosen on anatomy, not accuracy); Diverse
# Rush DLPFC (not STG); NIMH HBCC bulk (not LIBD's bulk of the same HBCC donors). The
# dropped matrices are the same people re-measured, so as separate cohorts they were
# pseudo-replicates in the selection and the cohort-level sign tests.
REP_BULK = ["ROSMAP", "MSBB_FP", "Diverse_Rush2_DLPFC", "CommonMind_MSSM", "NIMH_HBCC", "NABEC"]


def norm(ct):
    return ct.replace("/", "_").replace(" ", "_")


# ---- inputs -------------------------------------------------------------------
man = pd.read_csv(HERE / "strategy_pairs.tsv", sep="\t", dtype=str, keep_default_na=False)
man = man[man.role != "excluded"]
man = man[man.snrna.map(lambda s: set(s.split("+")) <= PAPER_SN)]
# REP_BULK is applied only in the figure (panel b); the strategy analysis (panels c, d)
# keeps every region (user decision 2026-11-04, restoring the previous panel c).
man["pair"] = man.bulk_cohort + "__" + man.snrna
TRAINERS = man[man.role.isin(["both", "trainer"])].set_index("pair")
cls_of = dict(zip(man.bulk_cohort, man["class"]))
fam_of = dict(zip(man.bulk_cohort, man.family))
members = {p: r.snrna.split("+") for p, r in TRAINERS.iterrows()}
bulk_of = {p: r.bulk_cohort for p, r in TRAINERS.iterrows()}

ci = pd.read_csv(COMMON / "reliability_ci_summary.tsv", sep="\t")
ci = ci[ci.level == "supertype"].assign(ct=lambda d: d.cell_type.map(norm))
CORE = set(ci.loc[ci.ci_lo > 0, "ct"])
ALL = set(ci.ct)

S = pd.read_csv(MAT / ("scores.csv" if METRIC == "pearson" else "scores_spearman.csv"))
print(f"metric: {METRIC}")
assert CORE <= set(S.cell_type), sorted(CORE - set(S.cell_type))
S = S[S.trainer.isin(set(man.pair)) & S.target.isin(set(man.pair))]
S["core"] = S.cell_type.isin(CORE)
S = S[S.cell_type.isin(ALL)]
print(f"scores: {len(S)} rows, {S.trainer.nunique()} trainers, {S.target.nunique()} targets; "
      f"Core {len(CORE)} of {len(ALL)} supertypes")

# unit-level median Core r per (trainer, target) from scores.csv (selection + oracle)
cell = (S[S.core].groupby(["trainer", "target", "diag"])
        .agg(med_core=("r_transfer", "median"), n=("n", "median")).reset_index())


def variant(R, s):
    """R itself if it does not contain s, else its declared leave-s-out variant."""
    if s not in members[R]:
        return R
    rem = [m for m in members[R] if m != s]
    if not rem:
        return None
    cand = f"{bulk_of[R]}__{'+'.join(rem)}"
    if cand in TRAINERS.index:
        return cand
    for p in TRAINERS.index:  # declared in another member order
        if bulk_of[p] == bulk_of[R] and sorted(members[p]) == sorted(rem):
            return p
    return None


# ---- units --------------------------------------------------------------------
tg = man[man.role.isin(["both", "target"]) & ~man.snrna.str.contains(r"\+")]
scorable = set(S.target)  # target pairs with any row (>= MIN_SCORE after hold-out)
units = []
for b, g in tg.groupby("bulk_cohort"):
    if cls_of.get(b) == "other":
        continue
    # both kinds are built; after the common-donor step a cohort keeps its cross units
    # if any survive the 50-donor floor, else its oof units
    for s in g.snrna:
        u = f"{b}__{s}"
        if u not in scorable:
            continue
        el = [p for p in TRAINERS.index if bulk_of[p] == b and s not in members[p]]
        ok = [p for p in el if ((cell.trainer == p) & (cell.target == u)).any()]
        base = dict(bulk=b, scorer=s, target=u, family=fam_of[b], cls=cls_of[b])
        if ok:
            units.append(dict(unit=f"{b}|{s}|cross", kind="cross", s1_pool=ok, **base))
        if ((cell.trainer == u) & (cell.target == u) & cell.diag).any():
            units.append(dict(unit=f"{b}|{s}|oof", kind="oof", s1_pool=[], **base))
U = pd.DataFrame(units)
print(f"\n{len(U)} candidate units over {U.bulk.nunique()} bulk cohorts "
      f"({(U.kind == 'cross').sum()} cross, {(U.kind == 'oof').sum()} oof)")


def s1_trainers(row):
    if row.kind == "oof":
        return [row.target]  # the diagonal
    pool = row.s1_pool
    multi = [p for p in pool if len(members[p]) > 1]
    if multi:
        k = max(len(members[p]) for p in multi)
        top = [p for p in multi if len(members[p]) == k]
        assert len(top) == 1, top
        return top
    return pool


U["S1"] = U.apply(s1_trainers, axis=1)


def unit_score(R, row):
    """median Core r of recipe R (variant-adjusted) on this unit, from scores.csv."""
    v = variant(R, row.scorer)
    if v is None:
        return np.nan, None
    m = cell[(cell.trainer == v) & (cell.target == row.target) & ~cell.diag]
    return (float(m.med_core.iloc[0]), v) if len(m) else (np.nan, v)


def select(cands, units_sel):
    """best candidate by mean over bulk cohorts of mean unit score; coverage >= COVER.
    Selection uses one unit per (bulk, scorer): cross if the cohort has any, else oof."""
    has_cross = set(units_sel.loc[units_sel.kind == "cross", "bulk"])
    units_sel = units_sel[(units_sel.kind == "cross") | ~units_sel.bulk.isin(has_cross)]
    nb = units_sel.bulk.nunique()
    best, rows = None, []
    for R in cands:
        sc = units_sel.apply(lambda r: unit_score(R, r)[0], axis=1)
        per_bulk = sc.groupby(units_sel.bulk.values).mean().dropna()
        cov = len(per_bulk) / nb if nb else 0
        crit = per_bulk.mean() if cov >= COVER else np.nan
        rows.append(dict(candidate=R, coverage=cov, criterion=crit))
    T = pd.DataFrame(rows).sort_values("criterion", ascending=False)
    ok = T.dropna(subset=["criterion"])
    return list(ok.candidate), T


sel_rows, choice = [], {}  # choice: ranked candidate list
for strat in ("S2", "S3"):
    for fam in sorted(U.family.unique()):
        for c in (sorted(U.cls.unique()) if strat == "S3" else [None]):
            if strat == "S3" and not (U[(U.family == fam)].cls == c).any():
                continue
            us = U[U.family != fam]
            cands = list(TRAINERS.index)
            if strat == "S3":
                us = us[us.cls == c]
                cands = [p for p in cands if cls_of[bulk_of[p]] == c]
            ranked, T = select(cands, us)
            choice[(strat, fam, c)] = ranked
            T = T.assign(strategy=strat, held_out_family=fam, cls=c,
                         rank=range(1, len(T) + 1))
            sel_rows.append(T)
            print(f"{strat} held-out {fam:8s} class {str(c):12s} -> "
                  f"{ranked[:2] if ranked else None}")
SEL = pd.concat(sel_rows)
SEL.to_csv(FIG / f"fig9_selection{SUF}.tsv", sep="\t", index=False)

# full-data choices (all families) for the transfer-only recommendation
glob_all = select(list(TRAINERS.index), U)[0][:1]
cls_all = {c: select([p for p in TRAINERS.index if cls_of[bulk_of[p]] == c], U[U.cls == c])[0][:1]
           for c in sorted(U.cls.unique())}
print(f"\nall-data choice: S2 {glob_all}; S3 {cls_all}")


def defined(R, row):
    """R's variant for this unit exists and was scored on it (availability, not value)."""
    v = variant(R, row.scorer)
    return v if v and ((cell.trainer == v) & (cell.target == row.target) & ~cell.diag).any() else None


def assign(row):
    """Rule fixed before scoring: the top-ranked selection candidate that is DEFINED for
    the unit (its leave-scorer-out variant exists and has a scores.csv row); the rank
    used is recorded."""
    out = {}
    for strat in ("S2", "S3"):
        ranked = choice[(strat, row.family, row.cls if strat == "S3" else None)]
        pick = next(((i + 1, defined(R, row)) for i, R in enumerate(ranked) if defined(R, row)),
                    (None, None))
        out[strat] = [pick[1]] if pick[1] else []
        out[f"{strat}_rank_used"] = pick[0]
    v = defined(PROD, row)
    out["prod"] = [v] if v else []
    sc = cell[(cell.target == row.target) & ~cell.diag]
    sc = sc[sc.trainer.map(lambda p: p in members and row.scorer not in members[p])]
    out["oracle_trainer"] = sc.loc[sc.med_core.idxmax(), "trainer"] if len(sc) else None
    out["oracle_scores_csv"] = sc.med_core.max() if len(sc) else np.nan
    return pd.Series(out)


U = pd.concat([U, U.apply(assign, axis=1)], axis=1)

# ---- common-donor recomputation -------------------------------------------------
_light = {}


def src_persons(b):
    """persons of b's SOURCE (light) matrix, as RM.load_bulk/read_target resolve them."""
    if b not in _light:
        idx = pd.read_parquet(Path(os.environ["CTCOMP_BULK"]) / f"{matrix_stem(b)}_light.parquet",
                              columns=[]).index.astype(str)
        p, _ = to_persons(idx, matrix_stem(b))
        s = set(p.dropna())
        if b in PART_OF:
            banks = bank_of(matrix_stem(b))
            s = {x for x in s if banks.get(x) == PART_OF[b][1]}
        _light[b] = s
    return _light[b]


labels = RM.load_counts(sorted(set(U.scorer) | {m for p in TRAINERS.index for m in members[p]}))


def train_people(p):
    lab = set().union(*[set(labels[m].index) for m in members[p]])
    return src_persons(bulk_of[p]) & lab


_pred = {}


def pred(p, b, diag):
    k = (p, b, diag)
    if k not in _pred:
        f = MAT / "oof" / f"{p}.parquet" if diag else MAT / "pred" / f"{p}__{b}.parquet"
        _pred[k] = pd.read_parquet(f) if f.exists() else None
    return _pred[k]


ct_rows, unit_rows = [], []
for row in U.itertuples():
    Y = labels[row.scorer]
    base = pred(row.target, row.bulk, True)  # target's own OOF baseline donors
    if base is None:
        continue
    D = set(base.index) & set(Y.index)
    strat_tr = {"S1": row.S1, "S2": row.S2, "S3": row.S3, "prod": row.prod}
    for st, trs in strat_tr.items():
        for p in trs:
            if not (row.kind == "oof" and st == "S1"):
                D -= train_people(p)
    D = sorted(D)
    rec = dict(unit=row.unit, bulk=row.bulk, scorer=row.scorer, kind=row.kind, cls=row.cls,
               family=row.family, n_common=len(D),
               S1_trainers=";".join(row.S1), S2_trainer=";".join(row.S2),
               S3_trainer=";".join(row.S3), prod_trainer=";".join(row.prod),
               S2_rank_used=row.S2_rank_used, S3_rank_used=row.S3_rank_used,
               oracle_trainer=row.oracle_trainer, oracle_med_core_scores_csv=row.oracle_scores_csv)
    if len(D) < MIN_SCORE:
        rec["note"] = f"common donors {len(D)} < {MIN_SCORE}; unit dropped"
        unit_rows.append(rec)
        continue
    for st, trs in strat_tr.items():
        if not trs:
            continue
        per = []
        for p in trs:
            P = pred(p, row.bulk, row.kind == "oof" and st == "S1")
            if P is None:
                continue
            r = {}
            for ct in P.columns.intersection(Y.columns):
                if ct not in ALL:
                    continue
                a, y = P[ct].reindex(D).values.astype(float), Y[ct].reindex(D).values.astype(float)
                ok = np.isfinite(a) & np.isfinite(y)
                if ok.sum() >= MIN_SCORE and a[ok].std() > 1e-9 and y[ok].std() > 1e-9:
                    r[ct] = corr(a[ok], y[ok])
            per.append(pd.Series(r, name=p))
        if not per:
            continue
        R = pd.concat(per, axis=1).mean(axis=1)  # S1 over singles: per-type mean
        for ct, v in R.items():
            ct_rows.append(dict(unit=row.unit, strategy=st, cell_type=ct, core=ct in CORE, r=v))
        rec[f"{st}_med_core"] = R[R.index.isin(CORE)].median()
        rec[f"{st}_med_all"] = R.median()
    unit_rows.append(rec)

CT = pd.DataFrame(ct_rows)
UN = pd.DataFrame(unit_rows)
# per cohort: keep cross units if any survived the common-donor floor, else oof units
UN["scored"] = UN.get("note").isna() if "note" in UN else True
live_cross = set(UN.loc[UN.scored & (UN.kind == "cross"), "bulk"])
UN["used"] = UN.scored & ((UN.kind == "cross") | ~UN.bulk.isin(live_cross))
CT = CT[CT.unit.isin(set(UN.loc[UN.used, "unit"]))]

# paired Wilcoxon per unit over Core cell types
W = CT[CT.core].pivot_table(index=["unit", "cell_type"], columns="strategy", values="r")
for a, b in (("S2", "S1"), ("S3", "S1"), ("S3", "S2"), ("S2", "prod"), ("S3", "prod")):
    ps = {}
    for u, g in W.groupby(level=0):
        if a in g and b in g:
            d = (g[a] - g[b]).dropna()
            ps[u] = stats.wilcoxon(d).pvalue if len(d) >= 10 and (d != 0).any() else np.nan
    UN[f"p_{a}_vs_{b}"] = UN.unit.map(ps)

CT.to_csv(FIG / f"fig9_celltype_r{SUF}.tsv", sep="\t", index=False)
UN.to_csv(FIG / f"fig9_units{SUF}.tsv", sep="\t", index=False)

# across bulk cohorts: mean over a cohort's units, then sign test
B = UN[UN.used].groupby(["bulk", "cls", "kind"])[
    [c for c in UN.columns if c.endswith("_med_core")]].mean()
summ = []
for a, b in (("S2", "S1"), ("S3", "S1"), ("S3", "S2"), ("S2", "prod"), ("S3", "prod")):
    for kinds in (("cross",), ("cross", "oof")):
        d = (B[f"{a}_med_core"] - B[f"{b}_med_core"])[B.index.get_level_values("kind").isin(kinds)].dropna()
        ties = int((d.abs() < 1e-12).sum())
        d = d[d.abs() >= 1e-12]  # same trainer under both strategies: tie, not a loss
        if len(d) or ties:
            w = int((d > 0).sum())
            summ.append(dict(comparison=f"{a} - {b}", units="+".join(kinds),
                             n_bulk=len(d) + ties, ties=ties, wins=w, losses=len(d) - w,
                             median_diff=d.median() if len(d) else 0.0,
                             sign_p=stats.binomtest(w, len(d)).pvalue if len(d) else np.nan))
SUM = pd.DataFrame(summ)
SUM.to_csv(FIG / f"fig9_summary{SUF}.tsv", sep="\t", index=False)
B.reset_index().to_csv(FIG / f"fig9_by_bulk{SUF}.tsv", sep="\t", index=False)

# trainer x target matrix (median Core r; diag = same-dataset OOF)
M = cell.assign(cv=False)
# Consensus x member-dataset cells have no donors after the hold-out; they are filled
# from cons_oof.py (5-fold out-of-fold, drawn circled like the diagonal). Heatmap ONLY:
# `cell` above (selection, oracle, strategies) does not see these rows.
cf = MAT / "cons_oof_scores.csv"
if cf.exists() and METRIC == "pearson":   # Pearson r_cv; the Spearman CV cells are in heatmap_spearman.py
    C = pd.read_csv(cf)
    C = C[C.cell_type.isin(CORE)]
    C = (C.groupby(["trainer", "target"]).agg(med_core=("r_cv", "median"), n=("n", "median"))
         .reset_index().assign(diag=False, cv=True))
    have = set(zip(M.trainer, M.target))
    C = C[[(a, b) not in have for a, b in zip(C.trainer, C.target)]]
    M = pd.concat([M, C], ignore_index=True)
    print(f"heatmap: added {len(C)} consensus out-of-fold cells from {cf.name}")
M = M.assign(trainer_cls=M.trainer.map(lambda p: cls_of.get(p.split("__")[0])),
             target_cls=M.target.map(lambda p: cls_of.get(p.split("__")[0])))
M.to_csv(FIG / f"fig9_trainer_target{SUF}.tsv", sep="\t", index=False)

pd.DataFrame([dict(strategy="S2", cls="all", recipe=";".join(glob_all))] +
             [dict(strategy="S3", cls=c, recipe=";".join(r)) for c, r in cls_all.items()]
             ).to_csv(FIG / f"fig9_recommended{SUF}.tsv", sep="\t", index=False)

pd.set_option("display.width", 250)
print("\n== units ==")
print(UN[["unit", "used", "cls", "n_common", "S1_med_core", "S2_med_core", "S3_med_core",
          "prod_med_core", "oracle_med_core_scores_csv", "S2_rank_used", "S3_rank_used"]]
      .round(3).to_string(index=False))
print("\n== across bulk cohorts ==")
print(SUM.round(4).to_string(index=False))
