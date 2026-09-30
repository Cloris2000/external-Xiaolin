#!/usr/bin/env python3
"""Detect overlapping individuals across the 15 GWAS cohorts actually analyzed.

Ground truth for "who was analyzed" = post-QC .psam IDs per cohort.
Identity evidence = /project/rrg-shreejoy/Public_datasets/_registry/genotype_matches.csv
(pairwise genotype matches, ibs0_rate < 0.005 == same person).

Method: build a graph where nodes are (cohort_source, sample_id) and edges are
called genotype matches; connected components = persons. Then any component
touching >=2 of the 15 GWAS cohorts is a duplicated individual.
"""
import csv, sys, glob, os
from collections import defaultdict

REG = "/project/rrg-shreejoy/Public_datasets/_registry"
RES = "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/results"

COHORTS = ["ROSMAP","ROSMAP_array","Mayo","MSBB","CMC_MSSM","CMC_PENN","CMC_PITT",
           "GTEx_v10","NABEC","NIMH_HBCC_1M","NIMH_HBCC_h650","NIMH_HBCC_Omni5M",
           "GVEX","AMP_AD_Rush","AMP_AD_Mayo"]

# ---- 1. load the analyzed sample list per GWAS cohort -----------------------
cohort_ids = {}
for c in COHORTS:
    f = sorted(glob.glob(os.path.join(RES, c, "*.psam")))[0]
    ids = []
    with open(f) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            ids.append(line.split("\t")[1].strip())
    cohort_ids[c] = set(ids)

# id -> list of GWAS cohorts that analyzed that literal id string
id_to_cohorts = defaultdict(set)
for c, ids in cohort_ids.items():
    for i in ids:
        id_to_cohorts[i].add(c)

# ---- 2. union-find over genotype match edges -------------------------------
parent = {}
def find(x):
    parent.setdefault(x, x)
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[ra] = rb

IBS0_THRESH = 0.005
edges = 0
with open(os.path.join(REG, "genotype_matches.csv")) as fh:
    for row in csv.DictReader(fh):
        if float(row["ibs0_rate"]) < IBS0_THRESH:
            union(row["sample_a"], row["sample_b"])
            edges += 1

# every analyzed id is a node even if it has no edge
for i in id_to_cohorts:
    find(i)

# ---- 3. collapse to persons, keep those spanning >=2 GWAS cohorts ----------
comp = defaultdict(set)   # root -> set of (cohort, id) actually analyzed
for i, cs in id_to_cohorts.items():
    r = find(i)
    for c in cs:
        comp[r].add((c, i))

dups = {r: m for r, m in comp.items() if len({c for c, _ in m}) >= 2}

pair_counts = defaultdict(set)
for r, m in dups.items():
    cs = sorted({c for c, _ in m})
    for a in range(len(cs)):
        for b in range(a + 1, len(cs)):
            pair_counts[(cs[a], cs[b])].add(r)

# ---- 4. report -------------------------------------------------------------
print(f"genotype-match edges used (ibs0<{IBS0_THRESH}): {edges}")
print(f"total analyzed samples across 15 cohorts: {sum(len(v) for v in cohort_ids.values())}")
print(f"distinct persons spanning >=2 GWAS cohorts: {len(dups)}")
print(f"redundant samples (double-counted): "
      f"{sum(len({c for c,_ in m}) - 1 for m in dups.values())}")
print()
print("=== cohort-pair overlap counts (persons shared) ===")
for (a, b), rs in sorted(pair_counts.items(), key=lambda kv: -len(kv[1])):
    print(f"{len(rs):5d}  {a} <-> {b}")

with open("/tmp/claude-3127300/overlap_pairs.tsv", "w") as out:
    out.write("cohort_a\tcohort_b\tn_shared_persons\n")
    for (a, b), rs in sorted(pair_counts.items(), key=lambda kv: -len(kv[1])):
        out.write(f"{a}\t{b}\t{len(rs)}\n")

with open("/tmp/claude-3127300/overlap_persons.tsv", "w") as out:
    out.write("person_component\tn_cohorts\tcohort\tanalyzed_sample_id\n")
    for n, (r, m) in enumerate(sorted(dups.items(), key=lambda kv: -len({c for c, _ in kv[1]})), 1):
        ncoh = len({c for c, _ in m})
        for c, i in sorted(m):
            out.write(f"DUP{n:05d}\t{ncoh}\t{c}\t{i}\n")
print("\nwrote /tmp/claude-3127300/overlap_pairs.tsv and overlap_persons.tsv")
