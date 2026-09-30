#!/usr/bin/env python3
"""For each of the 15 GWAS cohorts, how much of it is represented in the
genotype-matching panel at all? A cohort whose samples never appear in
genotype_matches.csv as either member of ANY pair was either (a) never
genotype-compared, or (b) compared and matched nothing. Those two cases look
identical in the matches table, so we quantify representation to flag blind spots.
"""
import csv, glob, os
from collections import defaultdict

REG = "/project/rrg-shreejoy/Public_datasets/_registry"
RES = "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/results"
COHORTS = ["ROSMAP","ROSMAP_array","Mayo","MSBB","CMC_MSSM","CMC_PENN","CMC_PITT",
           "GTEx_v10","NABEC","NIMH_HBCC_1M","NIMH_HBCC_h650","NIMH_HBCC_Omni5M",
           "GVEX","AMP_AD_Rush","AMP_AD_Mayo"]

cohort_ids = {}
for c in COHORTS:
    f = sorted(glob.glob(os.path.join(RES, c, "*.psam")))[0]
    with open(f) as fh:
        cohort_ids[c] = {l.split("\t")[1].strip() for l in fh if not l.startswith("#")}

matched_ids = set()
src_of = defaultdict(set)
with open(os.path.join(REG, "genotype_matches.csv")) as fh:
    for r in csv.DictReader(fh):
        matched_ids.add(r["sample_a"]); matched_ids.add(r["sample_b"])
        src_of[r["sample_a"]].add(r["source_a"])
        src_of[r["sample_b"]].add(r["source_b"])

print(f"{'cohort':<20} {'N':>5} {'in_panel':>9} {'pct':>6}   registry_sources_seen")
print("-" * 95)
for c in COHORTS:
    ids = cohort_ids[c]
    hit = ids & matched_ids
    srcs = sorted({s for i in hit for s in src_of[i]})
    pct = 100.0 * len(hit) / len(ids) if ids else 0
    print(f"{c:<20} {len(ids):>5} {len(hit):>9} {pct:>5.1f}%   {', '.join(srcs) if srcs else '-- NONE (blind spot) --'}")
