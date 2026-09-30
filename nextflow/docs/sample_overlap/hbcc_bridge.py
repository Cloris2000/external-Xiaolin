#!/usr/bin/env python3
"""Indirect test for HBCC cross-platform duplication.

The registry never compared HBCC chip-vs-chip, but it DID compare each HBCC
sample against PsychAD_WGS / Multiome_wgs. If two HBCC samples on different
platforms both match the SAME external WGS sample, they are the same person.
That bridges the untested pair without new genotype work.
"""
import csv, glob, os
from collections import defaultdict

REG = "/project/rrg-shreejoy/Public_datasets/_registry"
RES = "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/results"
HB = ["NIMH_HBCC_1M", "NIMH_HBCC_h650", "NIMH_HBCC_Omni5M"]

ids, owner = {}, {}
for c in HB:
    f = sorted(glob.glob(os.path.join(RES, c, "*.psam")))[0]
    with open(f) as fh:
        ids[c] = {l.split("\t")[1].strip() for l in fh if not l.startswith("#")}
    for i in ids[c]:
        owner[i] = c

# external anchor -> set of (platform, hbcc_id) matching it
anchor = defaultdict(set)
with open(os.path.join(REG, "genotype_matches.csv")) as fh:
    for r in csv.DictReader(fh):
        if float(r["ibs0_rate"]) >= 0.005:
            continue
        a, b, sa, sb = r["sample_a"], r["sample_b"], r["source_a"], r["source_b"]
        for hb_id, hb_src, ex_id, ex_src in ((a, sa, b, sb), (b, sb, a, sa)):
            if hb_id in owner and ex_id not in owner:
                anchor[(ex_src, ex_id)].add((owner[hb_id], hb_id))

bridged = {k: v for k, v in anchor.items() if len({p for p, _ in v}) >= 2}
print(f"external WGS anchors matching HBCC samples: {len(anchor)}")
print(f"anchors bridging >=2 DIFFERENT HBCC platforms: {len(bridged)}")
pair = defaultdict(set)
for k, v in bridged.items():
    ps = sorted({p for p, _ in v})
    for i in range(len(ps)):
        for j in range(i + 1, len(ps)):
            pair[(ps[i], ps[j])].add(k)
print()
for (a, b), ks in sorted(pair.items(), key=lambda kv: -len(kv[1])):
    print(f"{len(ks):4d}  {a} <-> {b}  (same person, bridged via external WGS)")
if bridged:
    print("\nexamples:")
    for k, v in list(sorted(bridged.items()))[:8]:
        print(f"  anchor {k[1]} ({k[0]}): " + ", ".join(f"{p}:{i}" for p, i in sorted(v)))
