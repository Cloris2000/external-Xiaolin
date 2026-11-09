"""Figure 8 inputs: does dissection (cortical) depth explain per-donor discordance?

Question (user, 2026-10-19): for the five independent matched-donor pairs used
in Figure 1 (common/reliability_stats.py), is a donor's cross-dataset CLR
discordance larger when the two datasets' tissue blocks sampled different
cortical depths? This is DISSECTION depth (which layers the block covered),
inferred from the excitatory supertypes present, via the cortdepth package
(/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/cortical-depth). It is not
sequencing depth; the earlier UMI-based version of this figure
(fig8_donor_depth.tsv, fig8_pair_donor.tsv) is superseded and left in place.

Depth discordance per shared donor (read-only, the PI's run):
  /scratch/shreejoy/nserc2027_aim3_fig/discordant/discordant_pairs_checked.csv
  = `cortdepth discordant --vintage 2026-10-13`, every person with >= 200
  excitatory nuclei in two datasets, scored on EXCITATORY nuclei only:
    wasserstein  earth mover's distance between the two samples' excitatory
                 depth profiles, in depth units (0 = pia, 1 = white-matter
                 boundary; 0.1 = an average shift of a tenth of the cortex)
    layer_tv     total-variation distance between L2/3, L4, L5, L6 shares
    mean_shift   B's mean excitatory depth minus A's (sign: B deeper > 0)
  plus the PI's checks (discordant_search.py): frontal, not the Mathys /
  Mathys_2019 duplicate, no label collapse (>= 2 effective excitatory
  supertypes on both sides), sex from XIST/chrY concordant with each other and
  the registry, identity not merely assumed; `pass` = all of these.
  Dataset orientation (dataset_a, dataset_b) matches PAIRS below for all five
  pairs (checked 2026-10-20), so mean_shift keeps the B-minus-A sign.

Donor ids: the registry keys people as BDR*; the wp3 caches behind Figures 1-2
use TLP*. releases/2026-10-13/legacy_ids.tsv maps TLP -> BDR one-to-one. The
registry's DFC-2026 labels ARE the wp3 cache labels (nucleus counts match
exactly for all five datasets, checked 2026-10-20), so the CLR discordance is
computed from the wp3 caches exactly as before and joined on the bridged id.

Composition discordance per shared donor, three flavours:
  inh     inhibitory supertypes only (CLR within the inhibitory set). The
          HEADLINE (user decision 2026-10-23): the depth score is built from
          excitatory composition, so excitatory types must not sit on both
          axes, and glial proportions (oligodendrocytes in white matter) are
          themselves largely a function of sampling depth, so glia are left
          out too. Same design as Figure 3 (depth from excitatory nuclei,
          agreement from inhibitory).
  nonexc  inhibitory + glial supertypes (CLR within each set), the previous
          headline, kept as a sensitivity column.
  all     every supertype (neuronal compartment = excitatory + inhibitory;
          glia), sensitivity only; partly circular.
  Each is the mean over that donor's cell types of |residual CLR_A - residual
  CLR_B|, residual = CLR minus the cell type's mean over the pair's shared
  donors in that dataset (so a systematic between-dataset offset is not
  counted). Supertype
  class (Excitatory / Inhibitory / Glia) comes from taxonomy_DFC_2026.tsv
  class_id, so one supertype is always in one set.

Output (data/ next to this script):
  fig8_depth_pair_donor.tsv   one row per (pair, donor) with depth scores on
                              both sides: TLP donor, BDR person, discordance
                              (inh, nonexc, all) with the number of types
                              behind each, wasserstein, layer_tv, mean_shift,
                              n_exc_a/b, exc_diversity_a/b, the PI's flags, pass
Run (login node is fine; the PsychAD cache needs ~6 GB):
  PYTHONNOUSERSITE=1 /home/zhoux156/miniforge3/envs/brisc/bin/python fig8_data.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "common"))
import reliability_stats as rs  # noqa: E402

TAX = rs.TAX
DEPTH = Path("/scratch/shreejoy/nserc2027_aim3_fig/discordant/discordant_pairs_checked.csv")
OUT = Path(__file__).resolve().parent / "data"

# (dataset a, dataset b, label): the eight pairs of common/reliability_stats.py
# (2026-10-29; were five, from the wp3 caches with a TLP -> BDR bridge). Count
# tables are the registry composition tables keyed by person_id, and the dataset
# ids are the registry dataset_ids the depth run uses, so no id bridge is needed.
# Orientation (dataset_a, dataset_b) matches the depth run for all eight
# (checked 2026-11-02; asserted below).
PAIRS = rs.PAIRS
PSEUDO = 0.5
DEPTH_COLS = ["wasserstein", "layer_tv", "mean_shift", "n_exc_a", "n_exc_b",
              "exc_diversity_a", "exc_diversity_b", "frontal", "duplicate",
              "label_collapse", "sex_ok", "identity_assumed", "pass"]


def clr(counts):
    p = counts + PSEUDO
    p = p / p.sum(axis=1, keepdims=True)
    L = np.log(p)
    return L - L.mean(axis=1, keepdims=True)


def pair_donor_discordance(tA, tB, groups):
    """Per shared donor, mean |residual CLR_A - residual CLR_B| over the supertypes
    of `groups` ({set name: supertypes}); CLR is taken within each set, residual =
    CLR minus the type's mean over the pair's shared donors (those kept for the
    set) in that dataset. Donors with no
    nuclei of a set in either dataset contribute no types from that set."""
    shared = sorted(set(tA.index) & set(tB.index))
    tot = pd.Series(0.0, index=shared)
    cnt = pd.Series(0, index=shared)
    for types in groups.values():
        types = [t for t in types if t in tA.columns and t in tB.columns]
        if len(types) < 2:
            continue
        dA = tA.reindex(index=shared, columns=types).fillna(0).values
        dB = tB.reindex(index=shared, columns=types).fillna(0).values
        keep = (dA.sum(1) > 0) & (dB.sum(1) > 0)
        if keep.sum() < 2:
            continue
        yA, yB = clr(dA[keep]), clr(dB[keep])
        diff = np.abs((yA - yA.mean(0)) - (yB - yB.mean(0)))
        donors = np.array(shared)[keep]
        tot.loc[donors] += diff.sum(1)
        cnt.loc[donors] += diff.shape[1]
    return pd.DataFrame({"donor": shared, "disc": (tot / cnt.replace(0, np.nan)).values,
                         "n_types": cnt.values})


def main():
    OUT.mkdir(exist_ok=True)
    tax = pd.read_csv(TAX, sep="\t")
    assert set(tax.class_id) <= {"Excitatory", "Inhibitory", "Glia"}, set(tax.class_id)
    by_class = {k: set(g.supertype_label) for k, g in tax.groupby("class_id")}
    groups = {
        "inh": {"inh": by_class["Inhibitory"]},
        "nonexc": {"inh": by_class["Inhibitory"], "glia": by_class["Glia"]},
        "all": {"neu": by_class["Excitatory"] | by_class["Inhibitory"], "glia": by_class["Glia"]},
    }

    depth = pd.read_csv(DEPTH)
    out, cache = [], {}
    for a, b, label in PAIRS:
        reg_a, reg_b = a, b
        for n in (a, b):
            if n not in cache:
                c = rs.load(n)
                cache[n] = c.loc[:, c.sum(axis=0) > 0]   # types seen, as the former crosstab
                print(f"{n}: {len(cache[n])} donors")
        tA, tB = cache[a], cache[b]
        d = None
        for flavour, g in groups.items():
            x = pair_donor_discordance(tA, tB, g).rename(
                columns={"disc": f"disc_{flavour}", "n_types": f"n_types_{flavour}"})
            d = x if d is None else d.merge(x, on="donor")
        d.insert(0, "pair", label)
        d["person"] = d.donor   # registry person_id already
        n_shared, n_mapped = len(d), int(d.person.notna().sum())

        sc = depth[(depth.dataset_a == reg_a) & (depth.dataset_b == reg_b)]
        assert len(sc), f"{label}: no rows for {reg_a} x {reg_b} in {DEPTH.name}"
        assert not depth[(depth.dataset_a == reg_b) & (depth.dataset_b == reg_a)].shape[0], \
            f"{label}: reversed orientation present in {DEPTH.name}"
        d = d.merge(sc[["person", *DEPTH_COLS]], on="person", how="inner")
        print(f"  {label:30s} shared {n_shared}, with depth scores "
              f"{len(d)}, pass {int(d['pass'].sum())}; W median {d.wasserstein.median():.3f}")
        out.append(d)

        later = {x for p in PAIRS[PAIRS.index((a, b, label)) + 1:] for x in p[:2]}
        for n in list(cache):
            if n not in later:
                del cache[n]

    res = pd.concat(out, ignore_index=True)
    res.to_csv(OUT / "fig8_depth_pair_donor.tsv", sep="\t", index=False)
    print(f"wrote {OUT}/fig8_depth_pair_donor.tsv ({len(res)} rows, {int(res['pass'].sum())} pass)")


if __name__ == "__main__":
    main()
