#!/usr/bin/env python3
"""Final consolidated overlap report for the 15-cohort meta-analysis."""
import csv, glob, os
from collections import defaultdict

REG = "/project/rrg-shreejoy/Public_datasets/_registry"
RES = "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/results"
COHORTS = ["ROSMAP","ROSMAP_array","Mayo","MSBB","CMC_MSSM","CMC_PENN","CMC_PITT",
           "GTEx_v10","NABEC","NIMH_HBCC_1M","NIMH_HBCC_h650","NIMH_HBCC_Omni5M",
           "GVEX","AMP_AD_Rush","AMP_AD_Mayo"]

ids = {}
for c in COHORTS:
    f = sorted(glob.glob(os.path.join(RES, c, "*.psam")))[0]
    ids[c] = {l.split("\t")[1].strip() for l in open(f) if not l.startswith("#")}

owner = defaultdict(set)
for c in COHORTS:
    for i in ids[c]:
        owner[i].add(c)

parent = {}
def find(x):
    parent.setdefault(x, x)
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: parent[ra] = rb

with open(os.path.join(REG, "genotype_matches.csv")) as fh:
    for r in csv.DictReader(fh):
        if float(r["ibs0_rate"]) < 0.005:
            union(r["sample_a"], r["sample_b"])

# add the ID-crosswalk evidence for the ROSMAP family (individualID namespace)
META = "/project/rrg-shreejoy/ROSMAP/Metadata/ROSMAP_biospecimen_metadata.csv"
HARM = "/project/rrg-shreejoy/ROSMAP/Metadata/RNAseq_Harmonization_ROSMAP_combined_metadata.csv"
RMAP = "/project/rrg-shreejoy/ROSMAP/Genotype/ROSMAP_array_sample_maps/samples_reheader_map.txt"
spec2ind = {}
for r in csv.DictReader(open(META)):
    if r["specimenID"].strip() and r["individualID"].strip():
        spec2ind[r["specimenID"].strip()] = r["individualID"].strip()
proj2ind = {}
for r in csv.DictReader(open(HARM)):
    p, i = r.get("projid", "").strip(), r.get("individualID", "").strip()
    if p and i: proj2ind[p.zfill(8)] = i
new2spec = {}
for line in open(RMAP):
    q = line.split()
    if len(q) == 2: new2spec[q[1]] = q[0].split("_")[0]

def to_ind(i):
    if i in spec2ind: return spec2ind[i]
    if i in new2spec and new2spec[i] in spec2ind: return spec2ind[new2spec[i]]
    if i[:3] in ("MAP","ROS") and i[3:].isdigit(): return proj2ind.get(i[3:].zfill(8))
    if i.startswith("R"): return i
    return None

ind_group = defaultdict(list)
for c in ("ROSMAP","ROSMAP_array","AMP_AD_Rush"):
    for i in ids[c]:
        v = to_ind(i)
        if v: ind_group[v].append(i)
for v, g in ind_group.items():
    for x in g[1:]:
        union(g[0], x)

comp = defaultdict(set)
for i, cs in owner.items():
    for c in cs:
        comp[find(i)].add((c, i))
dups = {r: m for r, m in comp.items() if len({c for c, _ in m}) >= 2}

redundant = defaultdict(int)   # cohort -> samples that duplicate another cohort
pair = defaultdict(set)
for r, m in dups.items():
    cs = sorted({c for c, _ in m})
    for a in range(len(cs)):
        for b in range(a+1, len(cs)):
            pair[(cs[a], cs[b])].add(r)
    # attribute the redundancy to every cohort but the largest contributor
    for c in cs[1:]:
        redundant[c] += 1

print("=" * 78)
print("OVERLAPPING INDIVIDUALS ACROSS THE 15-COHORT META-ANALYSIS")
print("=" * 78)
tot = sum(len(v) for v in ids.values())
print(f"\nAnalyzed samples (post-QC .psam) across 15 cohorts : {tot}")
print(f"Distinct individuals appearing in >=2 cohorts      : {len(dups)}")
print(f"Redundant (double-counted) samples                 : "
      f"{sum(len({c for c,_ in m})-1 for m in dups.values())}")

print("\n--- cohort pairs sharing individuals ---")
print(f"{'cohort A':<18}{'cohort B':<18}{'shared':>7}")
for (a, b), rs in sorted(pair.items(), key=lambda kv: -len(kv[1])):
    print(f"{a:<18}{b:<18}{len(rs):>7}")

print("\n--- per-cohort exposure ---")
print(f"{'cohort':<20}{'N':>6}{'dup':>6}{'pct':>8}")
for c in COHORTS:
    d = sum(1 for m in dups.values() if c in {x for x, _ in m})
    print(f"{c:<20}{len(ids[c]):>6}{d:>6}{100.0*d/len(ids[c]):>7.1f}%")

with open("/tmp/claude-3127300/FINAL_overlapping_individuals.tsv", "w") as out:
    out.write("dup_id\tn_cohorts\tcohort\tanalyzed_sample_id\n")
    for n, (r, m) in enumerate(sorted(dups.items(),
                               key=lambda kv: (-len({c for c,_ in kv[1]}), sorted(kv[1]))), 1):
        k = len({c for c, _ in m})
        for c, i in sorted(m):
            out.write(f"DUP{n:05d}\t{k}\t{c}\t{i}\n")
print("\nwrote /tmp/claude-3127300/FINAL_overlapping_individuals.tsv")
