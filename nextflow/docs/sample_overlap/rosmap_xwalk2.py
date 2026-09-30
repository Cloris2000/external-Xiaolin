#!/usr/bin/env python3
"""Close the ROSMAP WGS vs array question using the pipeline's own reheader map.

Chain: ROSMAP_array .psam id (MAP*/ROS*) <- samples_reheader_map.txt <- array
specimenID (e.g. 11AD39717) -> ROSMAP_biospecimen_metadata.csv -> individualID (R*).
ROSMAP WGS .psam ids (SM-*) map to individualID directly via the same metadata.
Both cohorts then live in one namespace and can be intersected exactly.
"""
import csv, glob, os

META = "/project/rrg-shreejoy/ROSMAP/Metadata/ROSMAP_biospecimen_metadata.csv"
HARM = "/project/rrg-shreejoy/ROSMAP/Metadata/RNAseq_Harmonization_ROSMAP_combined_metadata.csv"
RMAP = "/project/rrg-shreejoy/ROSMAP/Genotype/ROSMAP_array_sample_maps/samples_reheader_map.txt"
RES = "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/results"

def psam(c):
    f = sorted(glob.glob(os.path.join(RES, c, "*.psam")))[0]
    return {l.split("\t")[1].strip() for l in open(f) if not l.startswith("#")}

spec2ind = {}
with open(META) as fh:
    for r in csv.DictReader(fh):
        s, i = r["specimenID"].strip(), r["individualID"].strip()
        if s and i:
            spec2ind[s] = i

# projid (zero-padded) -> individualID, from the harmonization table
proj2ind = {}
with open(HARM) as fh:
    for r in csv.DictReader(fh):
        p, i = r.get("projid", "").strip(), r.get("individualID", "").strip()
        if p and i:
            proj2ind[p.zfill(8)] = i

# array .psam id -> array specimenID  (reverse the reheader map)
new2spec = {}
with open(RMAP) as fh:
    for line in fh:
        parts = line.split()
        if len(parts) == 2:
            old, new = parts
            new2spec[new] = old.split("_")[0]

rosmap, array, rush = psam("ROSMAP"), psam("ROSMAP_array"), psam("AMP_AD_Rush")

def to_ind(ids):
    out, unres = {}, []
    for i in ids:
        ind = None
        if i in spec2ind:                       # direct specimenID (WGS SM-*)
            ind = spec2ind[i]
        elif i in new2spec and new2spec[i] in spec2ind:   # array via reheader
            ind = spec2ind[new2spec[i]]
        elif i[:3] in ("MAP", "ROS") and i[3:].isdigit():  # study-prefixed projid
            ind = proj2ind.get(i[3:].zfill(8))
        elif i.startswith("R"):                 # already an individualID
            ind = i
        if ind:
            out[i] = ind
        else:
            unres.append(i)
    return out, unres

w, wu = to_ind(rosmap)
a, au = to_ind(array)
r, ru = to_ind(rush)
for n, ids, m, u in (("ROSMAP(WGS)", rosmap, w, wu), ("ROSMAP_array", array, a, au),
                     ("AMP_AD_Rush", rush, r, ru)):
    print(f"{n:<14} N={len(ids):4d}  resolved={len(m):4d}  unresolved={len(u):4d}"
          f"  distinct individuals={len(set(m.values()))}")

W, A, R = set(w.values()), set(a.values()), set(r.values())
print("\n=== ROSMAP-family individual-level overlap ===")
print(f"ROSMAP(WGS) n ROSMAP_array : {len(W & A)}")
print(f"ROSMAP(WGS) n AMP_AD_Rush  : {len(W & R)}")
print(f"ROSMAP_array n AMP_AD_Rush : {len(A & R)}")
print(f"all three                  : {len(W & A & R)}")

rev = {}
for src, m in (("ROSMAP", w), ("ROSMAP_array", a), ("AMP_AD_Rush", r)):
    for sid, ind in m.items():
        rev.setdefault(ind, []).append((src, sid))
with open("/tmp/claude-3127300/rosmap_family_overlap.tsv", "w") as out:
    out.write("individualID\tn_cohorts\tcohort\tanalyzed_sample_id\n")
    for ind in sorted((W & A) | (W & R) | (A & R)):
        e = rev[ind]
        for src, sid in sorted(e):
            out.write(f"{ind}\t{len({s for s,_ in e})}\t{src}\t{sid}\n")
print("\nwrote /tmp/claude-3127300/rosmap_family_overlap.tsv")
if au:
    print(f"\nunresolved ROSMAP_array examples: {au[:8]}")
