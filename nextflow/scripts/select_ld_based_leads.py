#!/usr/bin/env python3
"""
Ancestry-aware LD-based lead selection for pooled 15-cohort cell-type-proportion
(CTP) GWAS, plus effect extraction, allele harmonization, rsID annotation and
heterogeneity / sensitivity statistics for the "Generalizability" figure.

This script is additive: it does NOT overwrite any existing script and writes all
outputs under results/meta_sensitivity/generalizability/.

Pipeline (see the project plan for full spec):
  1. Reference-free variant normalization of pooled METAL markers.
  2. Genome-wide-significance filter (pooled fixed-effect P <= 5e-8).
  3. Ancestry-aware LD clumping (max r2 across in-sample AFR + 1000G EUR panels;
     500 kb distance fallback when no LD estimate is available).
  4. Cross-cell-type shared-locus map.
  5. rsID annotation from local hg19 dbSNP (exact chr:pos:REF:ALT).
  6. Effect extraction + allele harmonization (pooled / EUR / AFR / AMR / cohorts).
  7. Cross-cohort and cross-ancestry heterogeneity + random-effects + LOO.
  8. QC report + plot-ready TSVs.

Discovery uses ONLY the pooled fixed-effect P value. Heterogeneity, ancestry
availability and effect-direction concordance never influence lead selection.
No SNP-over-indel preference: the most significant variant in a clump is the lead
whether it is a SNP or an indel.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import math
import os
import shutil
import subprocess
import sys
from collections import defaultdict, OrderedDict
from pathlib import Path

# ────────────────────────────────────────────────────────────────────────────
# Static configuration (detected, real paths only)
# ────────────────────────────────────────────────────────────────────────────
PROJECT = Path("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow")

CELL_TYPES = [
    "Astrocyte", "Endothelial", "IT", "L4.IT", "L5.6.IT.Car3", "L5.6.NP",
    "L5.ET", "L6.CT", "L6b", "LAMP5", "Microglia", "OPC", "Oligodendrocyte",
    "PAX6", "PVALB", "Pericyte", "SST", "VIP", "VLMC",
]

# 15 pooled-meta cohorts (base result dirs hold their REGENIE step2 output)
COHORTS_15 = [
    "ROSMAP", "ROSMAP_array", "Mayo", "MSBB", "CMC_MSSM", "CMC_PENN", "CMC_PITT",
    "GTEx_v10", "NABEC", "GVEX", "NIMH_HBCC_1M", "NIMH_HBCC_h650",
    "NIMH_HBCC_Omni5M", "AMP_AD_Rush", "AMP_AD_Mayo",
]
# Cohort -> broad ancestry of the cohort's contribution to the POOLED meta
# (used only for annotation; several cohorts are ancestry-mixed).
COHORT_ANCESTRY = {
    "ROSMAP": "EUR", "ROSMAP_array": "EUR", "Mayo": "EUR", "MSBB": "EUR",
    "CMC_MSSM": "EUR", "CMC_PENN": "EUR", "CMC_PITT": "EUR", "GTEx_v10": "EUR",
    "NABEC": "EUR", "GVEX": "EUR", "NIMH_HBCC_1M": "AFR",
    "NIMH_HBCC_h650": "AFR", "NIMH_HBCC_Omni5M": "EUR", "AMP_AD_Rush": "AFR",
    "AMP_AD_Mayo": "AMR",
}

BROAD_CLASS = {
    "Astrocyte": "Non-neuronal", "Oligodendrocyte": "Non-neuronal",
    "OPC": "Non-neuronal", "Microglia": "Non-neuronal",
    "Endothelial": "Non-neuronal", "Pericyte": "Non-neuronal",
    "VLMC": "Non-neuronal",
    "VIP": "Inhibitory", "SST": "Inhibitory", "PVALB": "Inhibitory",
    "LAMP5": "Inhibitory", "PAX6": "Inhibitory",
    "IT": "Excitatory", "L4.IT": "Excitatory", "L5.6.IT.Car3": "Excitatory",
    "L5.6.NP": "Excitatory", "L5.ET": "Excitatory", "L6.CT": "Excitatory",
    "L6b": "Excitatory",
}

POOLED_DIR = PROJECT / "results/meta_analysis_15cohorts"
EUR_DIR = PROJECT / "results/meta_sensitivity/ancestry_specific_EUR"
AFR_DIR = PROJECT / "results/meta_sensitivity/ancestry_specific_AFR"
AMR_REGENIE = PROJECT / "results/AMP_AD_Mayo_AMR/regenie_step2/AMP_AD_Mayo_AMR_{ct}_step2_{ct}.regenie"

# LD panels
EUR_1000G_PREFIX = "/external/rprshnas01/kcni/mwainberg/ldsc/1000G_Phase3_plinkfiles/1000G.EUR.QC.{chrom}"
AFR_PGEN_DIRS = [
    PROJECT / "results/NIMH_HBCC_1M_AFR/CMC_HBCC.QC.final",
    PROJECT / "results/NIMH_HBCC_h650_AFR/CMC_HBCC.QC.final",
]

DBSNP_VCF = Path("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/data/dbsnp_human_hg37_rsID/common_all_20180423.vcf.gz")
DBSNP_BUILD = "dbSNP151_GRCh37"
LEAD_RSIDS_FALLBACK = PROJECT / "results/meta_sensitivity/ancestry_lead_effects/lead_rsids.tsv"

# Tools
PLINK2 = "/external/rprshnas01/kcni/mwainberg/software/plink2"
PLINK1 = shutil.which("plink") or "/nethome/kcni/xzhou/.anaconda3/bin/plink"
CONDA_BCFTOOLS_BIN = "/nethome/kcni/xzhou/.anaconda3/envs/bcftools_env/bin"
BCFTOOLS = f"{CONDA_BCFTOOLS_BIN}/bcftools"
TABIX = f"{CONDA_BCFTOOLS_BIN}/tabix"
BGZIP = f"{CONDA_BCFTOOLS_BIN}/bgzip"

OUT_DIR = PROJECT / "results/meta_sensitivity/generalizability"
LD_DIR = OUT_DIR / "ld_panels"
GENOME_BUILD = "hg19/GRCh37"
SQRT2 = math.sqrt(2.0)


# ────────────────────────────────────────────────────────────────────────────
# Small helpers
# ────────────────────────────────────────────────────────────────────────────
def to_float(v):
    if v is None:
        return None
    s = str(v).strip()
    if s == "" or s.upper() in ("NA", "NAN", ".", "-NAN", "INF", "-INF"):
        return None
    try:
        f = float(s)
        return None if (math.isnan(f) or math.isinf(f)) else f
    except ValueError:
        return None


def norm_chrom(c):
    return str(c).replace("chr", "").replace("Chr", "").strip()


def parse_marker(marker):
    """chr{c}:{pos}:{A}:{B} -> (chrom_no_chr, pos_int, A_up, B_up) or None."""
    parts = marker.strip().split(":")
    if len(parts) < 4:
        return None
    chrom = norm_chrom(parts[0])
    try:
        pos = int(parts[1])
    except ValueError:
        return None
    a = parts[2].upper()
    b = parts[3].upper()
    return chrom, pos, a, b


def variant_type(ref, alt):
    if "*" in (ref, alt) or "<" in alt or "<" in ref:
        return "symbolic_or_spanning"
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    if len(ref) > len(alt):
        return "deletion"
    if len(alt) > len(ref):
        return "insertion"
    return "mnv"


def is_palindromic(ref, alt):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    if len(ref) != 1 or len(alt) != 1:
        return False
    return comp.get(ref) == alt


def norm_sf(z):
    return math.erfc(abs(z) / SQRT2)


def norm_key(chrom, pos, a, b):
    """Order-independent allele key for cross-file matching."""
    return (norm_chrom(chrom), int(pos), frozenset((a.upper(), b.upper())))


def find_pooled_tbl(meta_dir, ct):
    """Prefer the METAL ANALYZE-HETEROGENEITY primary output (…1.tbl)."""
    cands = [p for p in meta_dir.glob(f"{ct}_meta_analysis_*.tbl") if p.stat().st_size > 0]
    if not cands:
        return None
    cands.sort(key=lambda p: (0 if p.name.endswith("1.tbl") else 1, -p.stat().st_size))
    return cands[0]


def sibling_tbl(primary):
    """The non-'1' sibling .tbl for validation, if present."""
    if primary.name.endswith("1.tbl"):
        sib = primary.with_name(primary.name[:-5] + ".tbl")
        return sib if sib.exists() else None
    return None


def write_tsv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def run(cmd, **kw):
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          text=True, **kw)


# ────────────────────────────────────────────────────────────────────────────
# Stage 1-2 : normalization + GW filter (reads pooled METAL .tbl)
# ────────────────────────────────────────────────────────────────────────────
TBL_COLS = ["MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE", "MinFreq",
            "MaxFreq", "Effect", "StdErr", "P-value", "Direction",
            "HetISq", "HetChiSq", "HetDf", "HetPVal"]


def read_tbl_header_idx(fh):
    header = fh.readline().decode("utf-8", "replace").rstrip("\n").split("\t")
    return {c: header.index(c) for c in TBL_COLS if c in header}, header


def scan_pooled_gws(ct, cfg, norm_audit, warnings):
    """Return list of GW-significant, normalized lead-candidate dicts for a cell type.

    Uses awk to prefilter rows by the P-value column (streaming C), then parses the
    few genome-wide-significant rows in Python.
    """
    tbl = find_pooled_tbl(POOLED_DIR, ct)
    if tbl is None:
        warnings.append(f"[{ct}] no pooled .tbl found")
        return [], 0
    with open(tbl, "rb") as fh:
        idx, header = read_tbl_header_idx(fh)
    het_ok = all(c in idx for c in ("HetISq", "HetPVal", "HetDf", "Direction"))
    if not het_ok:
        warnings.append(f"[{ct}] {tbl.name}: ANALYZE HETEROGENEITY columns missing")
    pcol = idx["P-value"] + 1  # awk is 1-based
    awk_prog = (f'NR>1 && ${pcol}!="" && (${pcol}+0)<={cfg.p_thresh}')
    r = run(["awk", "-F", "\t", awk_prog, str(tbl)])
    seen = {}
    gws = []
    n_below = 0
    for line in r.stdout.splitlines():
        parts = [p.encode() for p in line.rstrip("\n").split("\t")]
        maxi = max(idx.values())
        if len(parts) <= maxi:
            continue
        p = to_float(parts[idx["P-value"]].decode())
        if p is None or p > cfg.p_thresh:
            continue
        n_below += 1
        marker = parts[idx["MarkerName"]].decode().strip()
        pm = parse_marker(marker)
        if pm is None:
            norm_audit.append(dict(cell_type=ct, original_marker=marker,
                normalized_variant="", chromosome="", position="", REF="",
                ALT="", variant_type="", genome_build=GENOME_BUILD,
                normalization_status="unparseable_marker", notes="split<4"))
            continue
        chrom, pos, a_ref, a_alt = pm
        a1 = parts[idx["Allele1"]].decode().strip().upper()  # effect allele
        a2 = parts[idx["Allele2"]].decode().strip().upper()
        eff = to_float(parts[idx["Effect"]].decode())
        se = to_float(parts[idx["StdErr"]].decode())
        freq1 = to_float(parts[idx["Freq1"]].decode())
        direction = parts[idx["Direction"]].decode().strip() if "Direction" in idx else ""
        hetisq = to_float(parts[idx["HetISq"]].decode()) if "HetISq" in idx else None
        hetp = to_float(parts[idx["HetPVal"]].decode()) if "HetPVal" in idx else None
        hetdf = to_float(parts[idx["HetDf"]].decode()) if "HetDf" in idx else None
        n_cohorts = sum(1 for c in direction if c in "+-")
        vt = variant_type(a_ref, a_alt)
        normalized = f"chr{chrom}:{pos}:{a_ref}:{a_alt}"

        # required-field validity
        invalid = []
        if eff is None: invalid.append("beta")
        if se is None or se <= 0: invalid.append("se")
        if p is None: invalid.append("p")
        if not a1 or not a2: invalid.append("alleles")
        status = "reference_free_ok" if not invalid else "invalid:" + ",".join(invalid)
        notes = "REF unverified (no FASTA); indels not left-aligned"
        if "*" in (a_ref, a_alt):
            notes += "; spanning/symbolic allele"

        # duplicate / conflict detection within cell type
        if normalized in seen:
            prev = seen[normalized]
            if abs((prev["beta"] or 0) - (eff or 0)) > 1e-6 or \
               {prev["effect_allele"], prev["other_allele"]} != {a1, a2}:
                raise SystemExit(
                    f"CONFLICT in {ct}: normalized variant {normalized} has "
                    f"conflicting alleles/effects between duplicate rows.")
            norm_audit.append(dict(cell_type=ct, original_marker=marker,
                normalized_variant=normalized, chromosome=chrom, position=pos,
                REF=a_ref, ALT=a_alt, variant_type=vt, genome_build=GENOME_BUILD,
                normalization_status="duplicate_dropped", notes=notes))
            continue

        rec = dict(cell_type=ct, marker=marker, normalized=normalized,
                   chrom=chrom, pos=pos, ref=a_ref, alt=a_alt, vtype=vt,
                   effect_allele=a1, other_allele=a2, beta=eff, se=se,
                   freq1=freq1, p=p, direction=direction, n_cohorts=n_cohorts,
                   het_isq=hetisq, het_p=hetp, het_df=hetdf, invalid=invalid)
        seen[normalized] = rec
        norm_audit.append(dict(cell_type=ct, original_marker=marker,
            normalized_variant=normalized, chromosome=chrom, position=pos,
            REF=a_ref, ALT=a_alt, variant_type=vt, genome_build=GENOME_BUILD,
            normalization_status=status, notes=notes))
        if not invalid and n_cohorts >= cfg.min_cohorts:
            gws.append(rec)
    return gws, n_below


# ────────────────────────────────────────────────────────────────────────────
# Stage 3 : LD panels + ancestry-aware clumping
# ────────────────────────────────────────────────────────────────────────────
def build_afr_panel(warnings):
    """Merge the two hg19 in-sample AFR pgen filesets into one PLINK1 bed panel."""
    LD_DIR.mkdir(parents=True, exist_ok=True)
    merged = LD_DIR / "AFR_insample"
    if (merged.with_suffix(".bed")).exists() and (merged.with_suffix(".bim")).exists():
        return merged
    beds = []
    for i, pfx in enumerate(AFR_PGEN_DIRS):
        if not Path(str(pfx) + ".pgen").exists():
            warnings.append(f"AFR panel: missing {pfx}.pgen")
            continue
        out = LD_DIR / f"afr_part{i}"
        r = run([PLINK2, "--pfile", str(pfx), "--max-alleles", "2",
                 "--snps-only", "just-acgt", "--autosome",
                 "--set-all-var-ids", "chr@:#:$r:$a", "--new-id-max-allele-len", "200",
                 "--make-bed", "--out", str(out)])
        if not (out.with_suffix(".bed")).exists():
            warnings.append(f"AFR panel: plink2 make-bed failed for {pfx}: {r.stderr[-300:]}")
            continue
        beds.append(out)
    if not beds:
        warnings.append("AFR panel: no parts built")
        return None
    if len(beds) == 1:
        for ext in (".bed", ".bim", ".fam"):
            shutil.copy(str(beds[0]) + ext, str(merged) + ext)
        return merged
    mlist = LD_DIR / "afr_merge_list.txt"
    mlist.write_text("\n".join(str(b) for b in beds[1:]) + "\n")
    r = run([PLINK1, "--bfile", str(beds[0]), "--merge-list", str(mlist),
             "--make-bed", "--allow-no-sex", "--out", str(merged)])
    if not (merged.with_suffix(".bed")).exists():
        # merge failed (multiallelic mismatch) -> fall back to first part
        warnings.append(f"AFR panel merge failed ({r.stderr[-200:]}); using first part only")
        for ext in (".bed", ".bim", ".fam"):
            shutil.copy(str(beds[0]) + ext, str(merged) + ext)
    return merged


def index_bim_positions(bim_path, want_pos):
    """Map norm_key -> panel variant ID for positions we care about."""
    out = {}
    if not Path(bim_path).exists():
        return out
    with open(bim_path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 6:
                continue
            chrom, vid, _, pos, a1, a2 = f[0], f[1], f[2], f[3], f[4], f[5]
            try:
                ip = int(pos)
            except ValueError:
                continue
            if (norm_chrom(chrom), ip) not in want_pos:
                continue
            out[norm_key(chrom, ip, a1, a2)] = vid
    return out


def plink_pairwise_r2(bfile, chrom, panel_ids, extra_args=None):
    """Return {(idA,idB): r2} within 500kb for the given panel IDs on one chrom."""
    if not panel_ids:
        return {}
    tmp = LD_DIR / f"tmp_extract_{os.getpid()}.txt"
    tmp.write_text("\n".join(sorted(set(panel_ids))) + "\n")
    out = LD_DIR / f"tmp_ld_{os.getpid()}"
    cmd = [PLINK1, "--bfile", str(bfile), "--extract", str(tmp),
           "--r2", "--ld-window-kb", "500", "--ld-window", "999999",
           "--ld-window-r2", "0", "--allow-no-sex", "--out", str(out)]
    if chrom is not None:
        cmd += ["--chr", str(chrom)]
    if extra_args:
        cmd += extra_args
    run(cmd)
    res = {}
    ldf = out.with_suffix(".ld")
    if ldf.exists():
        with open(ldf) as fh:
            hdr = fh.readline().split()
            try:
                ia, ib, ir = hdr.index("SNP_A"), hdr.index("SNP_B"), hdr.index("R2")
            except ValueError:
                ia, ib, ir = 2, 5, 6
            for line in fh:
                p = line.split()
                if len(p) <= ir:
                    continue
                r2 = to_float(p[ir])
                if r2 is None:
                    continue
                res[(p[ia], p[ib])] = r2
                res[(p[ib], p[ia])] = r2
        ldf.unlink()
    for junk in (tmp, out.with_suffix(".log"), out.with_suffix(".nosex")):
        if Path(junk).exists():
            Path(junk).unlink()
    return res


def compute_r2_tables(ct_gws, cfg, warnings):
    """For each cell type, build per-panel r2 dicts keyed by normalized variant pairs.

    Returns: r2_by_ct[ct][panel][(normA,normB)] = r2, and panel availability sets.
    Panels: 'EUR' (1000G per chrom), 'AFR' (in-sample merged).
    """
    r2_by_ct = defaultdict(lambda: {"EUR": {}, "AFR": {}})
    panel_has = defaultdict(lambda: {"EUR": set(), "AFR": set()})  # normalized ids present

    afr_bfile = build_afr_panel(warnings) if cfg.use_ld else None

    for ct, gws in ct_gws.items():
        if not gws or not cfg.use_ld:
            continue
        by_chrom = defaultdict(list)
        for g in gws:
            by_chrom[g["chrom"]].append(g)
        for chrom, group in by_chrom.items():
            want_pos = {(g["chrom"], g["pos"]) for g in group}
            # ---- EUR 1000G (per chrom bed) ----
            eur_bim = EUR_1000G_PREFIX.format(chrom=chrom) + ".bim"
            eur_map = index_bim_positions(eur_bim, want_pos)
            id2norm_eur = {}
            for g in group:
                pid = eur_map.get(norm_key(g["chrom"], g["pos"], g["ref"], g["alt"]))
                if pid:
                    id2norm_eur[pid] = g["normalized"]
                    panel_has[ct]["EUR"].add(g["normalized"])
            if len(id2norm_eur) >= 2:
                eur_bfile = EUR_1000G_PREFIX.format(chrom=chrom)
                pr = plink_pairwise_r2(eur_bfile, chrom, list(id2norm_eur.keys()))
                for (ida, idb), r2 in pr.items():
                    na, nb = id2norm_eur.get(ida), id2norm_eur.get(idb)
                    if na and nb:
                        r2_by_ct[ct]["EUR"][(na, nb)] = r2
            # ---- AFR in-sample (merged bed) ----
            if afr_bfile is not None:
                afr_map = index_bim_positions(str(afr_bfile) + ".bim", want_pos)
                id2norm_afr = {}
                for g in group:
                    pid = afr_map.get(norm_key(g["chrom"], g["pos"], g["ref"], g["alt"]))
                    if pid:
                        id2norm_afr[pid] = g["normalized"]
                        panel_has[ct]["AFR"].add(g["normalized"])
                if len(id2norm_afr) >= 2:
                    pr = plink_pairwise_r2(afr_bfile, chrom, list(id2norm_afr.keys()))
                    for (ida, idb), r2 in pr.items():
                        na, nb = id2norm_afr.get(ida), id2norm_afr.get(idb)
                        if na and nb:
                            r2_by_ct[ct]["AFR"][(na, nb)] = r2
    return r2_by_ct, panel_has


def clump_cell_type(ct, gws, r2_tables, panel_has, cfg):
    """Greedy ancestry-aware clumping. Returns (leads, membership_rows)."""
    ordered = sorted(gws, key=lambda g: (g["p"], g["normalized"]))
    assigned = {}          # normalized -> lead normalized
    leads = []
    membership = []
    eur_r2 = r2_tables.get(ct, {}).get("EUR", {})
    afr_r2 = r2_tables.get(ct, {}).get("AFR", {})
    for idx in ordered:
        if idx["normalized"] in assigned:
            continue
        # new index variant (lead)
        lead = idx
        leads.append(lead)
        clump_members = [lead]
        assigned[lead["normalized"]] = lead["normalized"]
        lo = lead["pos"] - cfg.clump_kb * 1000
        hi = lead["pos"] + cfg.clump_kb * 1000
        for cand in ordered:
            if cand["normalized"] in assigned:
                continue
            if cand["chrom"] != lead["chrom"]:
                continue
            if not (lo <= cand["pos"] <= hi):
                continue
            key = (lead["normalized"], cand["normalized"])
            r2s = {}
            if key in eur_r2:
                r2s["EUR"] = eur_r2[key]
            if key in afr_r2:
                r2s["AFR"] = afr_r2[key]
            method = None
            max_r2 = None
            if r2s:
                max_r2 = max(r2s.values())
                method = "LD"
                assign = max_r2 >= cfg.clump_r2
            else:
                # no LD estimate in any panel -> distance fallback
                method = "distance_fallback"
                assign = True
                max_r2 = None
            if assign:
                assigned[cand["normalized"]] = lead["normalized"]
                clump_members.append(cand)
                membership.append(dict(
                    cell_type=ct, gws_variant=cand["normalized"],
                    assigned_lead=lead["normalized"], pooled_p=cand["p"],
                    EUR_r2=r2s.get("EUR"), AFR_r2=r2s.get("AFR"),
                    maximum_r2=max_r2, assignment_method=method,
                    distance_fallback_used=(method == "distance_fallback")))
        # membership row for the lead itself
        membership.append(dict(
            cell_type=ct, gws_variant=lead["normalized"],
            assigned_lead=lead["normalized"], pooled_p=lead["p"],
            EUR_r2=1.0, AFR_r2=1.0, maximum_r2=1.0,
            assignment_method="index", distance_fallback_used=False))
        # annotate lead with clump extent + neighbour info
        positions = [m["pos"] for m in clump_members]
        lead["clump_start"] = min(positions)
        lead["clump_end"] = max(positions)
        lead["n_gws_in_clump"] = len(clump_members)
        # nearest secondary GWS variant (2nd best by p in clump)
        others = [m for m in clump_members if m["normalized"] != lead["normalized"]]
        if others:
            sec = min(others, key=lambda m: m["p"])
            lead["nearest_secondary"] = sec["normalized"]
        else:
            lead["nearest_secondary"] = ""
        # per-lead LD panel bookkeeping
        panels = []
        if lead["normalized"] in panel_has.get(ct, {}).get("EUR", set()):
            panels.append("EUR")
        if lead["normalized"] in panel_has.get(ct, {}).get("AFR", set()):
            panels.append("AFR")
        lead["ld_panels_available"] = ";".join(panels) if panels else "none"
        # was any assignment in this clump a distance fallback?
        fb = any(m.get("assignment_method") == "distance_fallback"
                 for m in membership if m["assigned_lead"] == lead["normalized"])
        lead["distance_fallback_used"] = fb
        lead["clumping_method"] = "LD+distance" if lead["ld_panels_available"] != "none" else "distance_only"
    return leads, membership


# ────────────────────────────────────────────────────────────────────────────
# Cell-type display order (manuscript class blocks: Non-neuronal, Inhib, Excit)
# ────────────────────────────────────────────────────────────────────────────
CT_ORDER = [
    "Microglia", "Endothelial", "Oligodendrocyte", "Pericyte", "VLMC",
    "Astrocyte", "OPC",
    "LAMP5", "VIP", "PAX6", "PVALB", "SST",
    "IT", "L5.6.IT.Car3", "L5.6.NP", "L5.ET", "L6.CT", "L4.IT", "L6b",
]
CT_RANK = {c: i for i, c in enumerate(CT_ORDER)}
CLASS_ORDER = {"Non-neuronal": 0, "Inhibitory": 1, "Excitatory": 2}

RG = shutil.which("rg") or "rg"
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


# ────────────────────────────────────────────────────────────────────────────
# Stage 5 : rsID annotation from local hg19 dbSNP (exact chr:pos:REF:ALT)
# ────────────────────────────────────────────────────────────────────────────
def ensure_dbsnp_index(warnings):
    """Return a tabix-queryable dbSNP path (symlink + .tbi under our out dir)."""
    if not DBSNP_VCF.exists():
        warnings.append(f"dbSNP VCF missing: {DBSNP_VCF}")
        return None
    local = LD_DIR / "dbsnp_hg19.vcf.gz"
    LD_DIR.mkdir(parents=True, exist_ok=True)
    if not local.exists():
        try:
            local.symlink_to(DBSNP_VCF)
        except OSError:
            local = DBSNP_VCF  # fall back to original location
    tbi = Path(str(local) + ".tbi")
    if not tbi.exists():
        r = run([TABIX, "-p", "vcf", str(local)])
        if not tbi.exists():
            warnings.append(f"tabix index build failed: {r.stderr[-300:]}")
            return None
    return local


def load_rsid_fallback():
    m = {}
    if LEAD_RSIDS_FALLBACK.exists():
        for row in csv.DictReader(open(LEAD_RSIDS_FALLBACK), delimiter="\t"):
            mk = row.get("marker") or row.get("normalized") or ""
            rs = row.get("rsid") or row.get("rsID") or ""
            if mk and rs and rs.startswith("rs"):
                m[norm_chrom(mk.split(":")[0]) + ":" + mk.split(":")[1]] = rs
    return m


def annotate_rsids(leads, dbsnp_path, fallback, warnings):
    audit = []
    for lead in leads:
        chrom, pos = lead["chrom"], lead["pos"]
        ref, alt = lead["ref"], lead["alt"]
        rsid, status, source = "", "no_rsID", "dbSNP"
        if dbsnp_path is not None:
            reg = f"{chrom}:{pos}-{pos}"
            r = run([TABIX, str(dbsnp_path), reg])
            hits = []
            for line in r.stdout.splitlines():
                f = line.split("\t")
                if len(f) < 5:
                    continue
                d_ref = f[3].upper()
                d_alts = [a.upper() for a in f[4].split(",")]
                d_rs = f[2]
                # exact match: our allele pair == dbSNP {REF, one ALT} (order-free)
                if {ref, alt} == {d_ref, alt} and alt in d_alts and ref == d_ref:
                    hits.append(d_rs)
                elif ref == d_ref and alt in d_alts:
                    hits.append(d_rs)
                elif {ref, alt} == set([d_ref] + d_alts) or \
                     (ref in ([d_ref] + d_alts) and alt in ([d_ref] + d_alts) and ref != alt):
                    hits.append(d_rs)
            hits = [h for h in dict.fromkeys(hits) if h.startswith("rs")]
            if len(hits) == 1:
                rsid, status = hits[0], "exact_match"
            elif len(hits) > 1:
                rsid, status = hits[0], "multiple_rsIDs"
            elif r.stdout.strip():
                status = "allele_mismatch"
        if not rsid:
            fb = fallback.get(f"{chrom}:{pos}")
            if fb:
                rsid, status, source = fb, "exact_match", "lead_rsids.tsv"
        lead["rsid"] = rsid
        lead["rsid_status"] = status
        audit.append(dict(normalized_variant=lead["normalized"], chromosome=chrom,
            position=pos, REF=ref, ALT=alt, rsID=rsid, mapping_status=status,
            mapping_source=source, dbSNP_build=DBSNP_BUILD, genome_build=GENOME_BUILD,
            notes=("reference-free alleles; matched allele set at position"
                   if status != "no_rsID" else "no matching dbSNP record")))
    return audit


# ────────────────────────────────────────────────────────────────────────────
# Stage 6 : effect extraction + allele harmonization
# ────────────────────────────────────────────────────────────────────────────
def harmonize(eff_allele, oth_allele, beta, se, p, freq, tgt_eff, tgt_oth,
              tgt_freq=None):
    """Align an estimate to the target effect allele. Returns dict with status."""
    if beta is None:
        return dict(beta=None, se=se, p=p, freq=freq, status="missing")
    ea, oa = eff_allele.upper(), oth_allele.upper()
    te, to = tgt_eff.upper(), tgt_oth.upper()
    palindromic = is_palindromic(te, to)
    # direct / reverse (order) match
    if {ea, oa} == {te, to}:
        if ea == te:
            b, fr, st = beta, freq, "ok"
        else:
            b = -beta
            fr = (1 - freq) if freq is not None else None
            st = "flipped"
    else:
        # strand complement (SNPs only)
        eac, oac = COMP.get(ea), COMP.get(oa)
        if eac and oac and {eac, oac} == {te, to}:
            if eac == te:
                b, fr, st = beta, freq, "strand_ok"
            else:
                b = -beta
                fr = (1 - freq) if freq is not None else None
                st = "strand_flipped"
        else:
            return dict(beta=None, se=se, p=p, freq=freq, status="allele_mismatch")
    # palindrome resolution via allele frequency
    if palindromic:
        if fr is None or tgt_freq is None or 0.4 <= min(fr, 1 - fr) <= 0.5:
            return dict(beta=None, se=se, p=p, freq=fr,
                        status="excluded_palindromic_ambiguous")
        # if oriented freq is on the opposite side of 0.5 vs target -> flip
        if (fr - 0.5) * (tgt_freq - 0.5) < 0:
            b = -b
            fr = 1 - fr
            st = st + "+freqflip"
    return dict(beta=b, se=se, p=p, freq=fr, status=st)


def rg_extract(path, patterns):
    """Return matching data lines from a text file using ripgrep (fixed strings)."""
    if not Path(path).exists() or not patterns:
        return []
    cmd = [RG, "-N", "--no-heading"]
    for pat in patterns:
        cmd += ["-e", pat]
    cmd.append(str(path))
    r = run(cmd)
    return r.stdout.splitlines()


def lookup_metal_tbl(tbl, leads):
    """Map norm_key -> record (effect=Allele1) for the given leads in a METAL tbl."""
    out = {}
    if tbl is None or not Path(tbl).exists():
        return out
    pats = sorted({f"chr{l['chrom']}:{l['pos']}:" for l in leads})
    for line in rg_extract(tbl, pats):
        f = line.split("\t")
        if len(f) < 10:
            continue
        pm = parse_marker(f[0])
        if pm is None:
            continue
        chrom, pos, a, b = pm
        rec = dict(eff=f[1].upper(), oth=f[2].upper(), freq=to_float(f[3]),
                   beta=to_float(f[7]), se=to_float(f[8]), p=to_float(f[9]))
        out[norm_key(chrom, pos, a, b)] = rec
    return out


def parse_regenie_header(path):
    with open(path) as fh:
        hdr = fh.readline().split()
    return {c: i for i, c in enumerate(hdr)}


def lookup_regenie(path, leads):
    out = {}
    if not Path(path).exists():
        return out
    idx = parse_regenie_header(path)
    need = ["ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P"]
    if not all(c in idx for c in need):
        return out
    pats = sorted({f"chr{l['chrom']}:{l['pos']}:" for l in leads})
    for line in rg_extract(path, pats):
        f = line.split()
        if len(f) <= idx["LOG10P"]:
            continue
        pm = parse_marker(f[idx["ID"]])
        if pm is None:
            continue
        chrom, pos, a, b = pm
        log10p = to_float(f[idx["LOG10P"]])
        rec = dict(eff=f[idx["ALLELE1"]].upper(), oth=f[idx["ALLELE0"]].upper(),
                   freq=to_float(f[idx["A1FREQ"]]), beta=to_float(f[idx["BETA"]]),
                   se=to_float(f[idx["SE"]]),
                   p=(10 ** (-log10p) if log10p is not None else None),
                   n=(int(float(f[idx["N"]])) if to_float(f[idx["N"]]) else None))
        out[norm_key(chrom, pos, a, b)] = rec
    return out


# ────────────────────────────────────────────────────────────────────────────
# Stage 7 : heterogeneity / sensitivity
# ────────────────────────────────────────────────────────────────────────────
def fixed_effect(pairs):
    w = [1 / (s * s) for _, s in pairs]
    sw = sum(w)
    b = sum(wi * bi for (bi, _), wi in zip(pairs, w)) / sw
    se = math.sqrt(1 / sw)
    return b, se, norm_sf(b / se)


def dl_random_effects(pairs):
    k = len(pairs)
    w = [1 / (s * s) for _, s in pairs]
    sw = sum(w)
    bfe = sum(wi * bi for (bi, _), wi in zip(pairs, w)) / sw
    Q = sum(wi * (bi - bfe) ** 2 for (bi, _), wi in zip(pairs, w))
    df = k - 1
    if df > 0:
        c = sw - sum(wi * wi for wi in w) / sw
        tau2 = max(0.0, (Q - df) / c) if c > 0 else 0.0
        i2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
        hetp = norm_sf(0) if Q <= 0 else _chisq_sf(Q, df)
    else:
        tau2, i2, hetp = 0.0, None, None
    wr = [1 / (s * s + tau2) for _, s in pairs]
    swr = sum(wr)
    br = sum(wi * bi for (bi, _), wi in zip(pairs, wr)) / swr
    ser = math.sqrt(1 / swr)
    return dict(beta=br, se=ser, p=norm_sf(br / ser), tau2=tau2, i2=i2,
                Q=Q, df=df, het_p=hetp)


def _chisq_sf(x, k):
    """Survival function of chi-square via regularized upper incomplete gamma."""
    if x <= 0:
        return 1.0
    a = k / 2.0
    xx = x / 2.0
    # use series/continued fraction for regularized gamma Q(a,x)
    if xx < a + 1:
        # series for P then Q=1-P
        term = 1.0 / a
        summ = term
        n = a
        for _ in range(500):
            n += 1
            term *= xx / n
            summ += term
            if term < summ * 1e-12:
                break
        p = summ * math.exp(-xx + a * math.log(xx) - math.lgamma(a))
        return max(0.0, 1.0 - p)
    else:
        b = xx + 1 - a
        c = 1e300
        d = 1.0 / b
        h = d
        for i in range(1, 500):
            an = -i * (i - a)
            b += 2
            d = an * d + b
            if abs(d) < 1e-30:
                d = 1e-30
            c = b + an / c
            if abs(c) < 1e-30:
                c = 1e-30
            d = 1.0 / d
            dl = d * c
            h *= dl
            if abs(dl - 1) < 1e-12:
                break
        return h * math.exp(-xx + a * math.log(xx) - math.lgamma(a))


def cross_cohort_stats(pairs_named, pooled_beta):
    """pairs_named = list of (cohort, beta, se) harmonized. Returns dict."""
    pairs = [(b, s) for _, b, s in pairs_named]
    k = len(pairs)
    if k == 0:
        return dict(k=0)
    fe_b, fe_se, fe_p = fixed_effect(pairs)
    re = dl_random_effects(pairs) if k >= 2 else dict(beta=fe_b, se=fe_se, p=fe_p,
        tau2=0, i2=None, Q=0, df=0, het_p=None)
    ref_sign = math.copysign(1, pooled_beta if pooled_beta else fe_b)
    concord = sum(1 for b, _ in pairs if math.copysign(1, b) == ref_sign)
    # leave-one-out
    loo_signs_stable = True
    loo_bmin = loo_bmax = None
    if k >= 2:
        for i in range(k):
            sub = pairs[:i] + pairs[i + 1:]
            lb, _, _ = fixed_effect(sub)
            loo_bmin = lb if loo_bmin is None else min(loo_bmin, lb)
            loo_bmax = lb if loo_bmax is None else max(loo_bmax, lb)
            if math.copysign(1, lb) != ref_sign:
                loo_signs_stable = False
    return dict(k=k, fe_beta=fe_b, fe_se=fe_se, fe_p=fe_p,
                re_beta=re["beta"], re_se=re["se"], re_p=re["p"],
                i2=re["i2"], Q=re["Q"], het_p=re["het_p"],
                pct_concordant=100.0 * concord / k,
                loo_sign_stable=loo_signs_stable,
                loo_beta_min=loo_bmin, loo_beta_max=loo_bmax)


def cross_ancestry_i2(estimates):
    """estimates = list of (beta, se) across ancestry metas. Returns (i2, Q, hetp, k)."""
    est = [(b, s) for b, s in estimates if b is not None and s and s > 0]
    k = len(est)
    if k < 2:
        return None, None, None, k
    w = [1 / (s * s) for _, s in est]
    bbar = sum(wi * bi for (bi, _), wi in zip(est, w)) / sum(w)
    Q = sum(wi * (bi - bbar) ** 2 for (bi, _), wi in zip(est, w))
    i2 = max(0.0, (Q - (k - 1)) / Q * 100) if Q > 0 else 0.0
    hetp = _chisq_sf(Q, k - 1)
    return i2, Q, hetp, k


# ────────────────────────────────────────────────────────────────────────────
# Stage 4 : cross-cell-type shared-locus map
# ────────────────────────────────────────────────────────────────────────────
def build_locus_map(all_leads, cfg):
    """Group leads across cell types into shared regions by physical overlap."""
    rows = []
    by_chrom = defaultdict(list)
    for lead in all_leads:
        by_chrom[lead["chrom"]].append(lead)
    locus_id = 0
    lead_to_locus = {}
    for chrom, leads in by_chrom.items():
        leads.sort(key=lambda l: l["pos"])
        clusters = []
        for lead in leads:
            placed = False
            for cl in clusters:
                if abs(lead["pos"] - cl["center"]) <= cfg.clump_kb * 1000:
                    cl["members"].append(lead)
                    cl["center"] = sum(m["pos"] for m in cl["members"]) / len(cl["members"])
                    placed = True
                    break
            if not placed:
                clusters.append(dict(center=lead["pos"], members=[lead]))
        for cl in clusters:
            locus_id += 1
            positions = [m["pos"] for m in cl["members"]]
            cts = sorted({m["cell_type"] for m in cl["members"]})
            for m in cl["members"]:
                lead_to_locus[(m["cell_type"], m["normalized"])] = locus_id
                rel = "shared_multi_celltype" if len(cts) > 1 else "single_celltype"
                rows.append(dict(
                    shared_locus_id=f"L{locus_id:03d}", cell_type=m["cell_type"],
                    lead_variant=m["normalized"], rsID=m.get("rsid", ""),
                    chromosome=chrom, position=m["pos"],
                    locus_start=min(positions), locus_end=max(positions),
                    relationship_to_other_leads=rel,
                    LD_evidence="physical_overlap(<=%dkb)" % cfg.clump_kb,
                    notes=("cell types at locus: " + ",".join(cts))))
    return rows, lead_to_locus, locus_id


# ────────────────────────────────────────────────────────────────────────────
# main
# ────────────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--p-thresh", type=float, default=5e-8)
    p.add_argument("--clump-r2", type=float, default=0.1)
    p.add_argument("--clump-kb", type=int, default=500)
    p.add_argument("--min-cohorts", type=int, default=2)
    p.add_argument("--min-total-n", type=int, default=None)
    p.add_argument("--reconciliation", choices=["max"], default="max",
                   help="LD reconciliation across ancestry panels")
    p.add_argument("--no-ld", action="store_true",
                   help="Distance-only clumping (skip LD panels)")
    p.add_argument("--cell-types", nargs="+", default=None)
    return p.parse_args()


def main():
    args = parse_args()
    args.use_ld = not args.no_ld
    args.min_cohorts = args.min_cohorts
    cfg = args
    cell_types = args.cell_types or CELL_TYPES
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LD_DIR.mkdir(parents=True, exist_ok=True)
    warnings = []
    norm_audit = []

    print(f"[cfg] p<={cfg.p_thresh:g} r2>={cfg.clump_r2} kb={cfg.clump_kb} "
          f"min_cohorts={cfg.min_cohorts} LD={'on' if cfg.use_ld else 'off'}")

    # Stage 1-2 : normalize + GW filter (parallel awk prefilter across cell types)
    ct_gws = OrderedDict((ct, []) for ct in cell_types)
    n_below = {}

    def _scan(ct):
        local_audit, local_warn = [], []
        gws, nb = scan_pooled_gws(ct, cfg, local_audit, local_warn)
        return ct, gws, nb, local_audit, local_warn

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        for ct, gws, nb, la, lw in ex.map(_scan, cell_types):
            ct_gws[ct] = gws
            n_below[ct] = nb
            norm_audit.extend(la)
            warnings.extend(lw)
            print(f"[{ct}] GW-sig variants (n_cohorts>={cfg.min_cohorts}): "
                  f"{len(gws)}  (raw P<=thr: {nb})", flush=True)

    # Stage 3 : LD panels + clumping
    r2_tables, panel_has = compute_r2_tables(ct_gws, cfg, warnings)
    all_leads = []
    membership_rows = []
    for ct in cell_types:
        leads, mem = clump_cell_type(ct, ct_gws[ct], r2_tables, panel_has, cfg)
        for l in leads:
            l["cell_type"] = ct
            l["broad_class"] = BROAD_CLASS[ct]
        all_leads.extend(leads)
        membership_rows.extend(mem)

    print(f"\n[leads] total lead SNP-cell-type associations: {len(all_leads)}")

    # Stage 5 : rsID
    dbsnp = ensure_dbsnp_index(warnings)
    fallback = load_rsid_fallback()
    rsid_audit = annotate_rsids(all_leads, dbsnp, fallback, warnings)

    # Stage 4 : cross-cell-type locus map
    locus_rows, lead_to_locus, n_loci = build_locus_map(all_leads, cfg)

    # Stage 6-7 : effect extraction, harmonization, heterogeneity
    harmon_audit = []
    effects_long = []
    lead_annot = []
    panelA_rows = []
    panelB_rows = []

    leads_by_ct = defaultdict(list)
    for lead in all_leads:
        leads_by_ct[lead["cell_type"]].append(lead)

    for ct_i, (ct, leads) in enumerate(leads_by_ct.items(), 1):
        print(f"[extract {ct_i}/{len(leads_by_ct)}] {ct}: {len(leads)} leads",
              flush=True)
        # parallel I/O: ancestry metas, AMR, and 15 cohort REGENIE files at once
        tasks = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex:
            tasks["__EUR__"] = ex.submit(lookup_metal_tbl, find_pooled_tbl(EUR_DIR, ct), leads)
            tasks["__AFR__"] = ex.submit(lookup_metal_tbl, find_pooled_tbl(AFR_DIR, ct), leads)
            tasks["__AMR__"] = ex.submit(lookup_regenie, str(AMR_REGENIE).format(ct=ct), leads)
            for coh in COHORTS_15:
                rgpath = PROJECT / f"results/{coh}/regenie_step2/{coh}_{ct}_step2_{ct}.regenie"
                tasks[coh] = ex.submit(lookup_regenie, rgpath, leads)
            results = {k: v.result() for k, v in tasks.items()}
        eur_map = results["__EUR__"]
        afr_map = results["__AFR__"]
        amr_map = results["__AMR__"]
        cohort_maps = {coh: results[coh] for coh in COHORTS_15}

        for lead in leads:
            key = norm_key(lead["chrom"], lead["pos"], lead["ref"], lead["alt"])
            te, to = lead["effect_allele"], lead["other_allele"]
            tfreq = lead["freq1"]

            # pooled (self)
            strat = OrderedDict()
            strat["Pooled"] = dict(beta=lead["beta"], se=lead["se"], p=lead["p"],
                                   freq=lead["freq1"], status="ok")
            # ancestry metas + AMR single cohort
            for name, mp in (("EUR", eur_map), ("AFR", afr_map)):
                rec = mp.get(key)
                if rec is None:
                    strat[name] = dict(beta=None, se=None, p=None, freq=None,
                                       status="missing")
                else:
                    strat[name] = harmonize(rec["eff"], rec["oth"], rec["beta"],
                        rec["se"], rec["p"], rec["freq"], te, to, tfreq)
            amr = amr_map.get(key)
            if amr is None:
                strat["AMR"] = dict(beta=None, se=None, p=None, freq=None,
                                    status="missing")
            else:
                strat["AMR"] = harmonize(amr["eff"], amr["oth"], amr["beta"],
                    amr["se"], amr["p"], amr["freq"], te, to, tfreq)

            for sname, s in strat.items():
                if s["status"] not in ("ok", "flipped", "strand_ok",
                                       "strand_flipped", "missing") \
                        and not s["status"].startswith(("ok", "flipped")):
                    harmon_audit.append(dict(cell_type=ct, lead=lead["normalized"],
                        target=sname, effect_allele=te, other_allele=to,
                        status=s["status"], notes="stratum-level"))
                b = s["beta"]
                lo = (b - 1.96 * s["se"]) if (b is not None and s["se"]) else None
                hi = (b + 1.96 * s["se"]) if (b is not None and s["se"]) else None
                effects_long.append(dict(cell_type=ct, lead_variant=lead["normalized"],
                    rsID=lead.get("rsid", ""), stratum=sname, level="ancestry",
                    beta=b, se=s["se"], ci_lo=lo, ci_hi=hi, p=s["p"],
                    freq_effect=s["freq"], status=s["status"]))

            # per-cohort harmonized
            coh_named = []
            for coh in COHORTS_15:
                rec = cohort_maps[coh].get(key)
                if rec is None:
                    h = dict(beta=None, se=None, p=None, freq=None, status="missing")
                    n_c = None
                else:
                    h = harmonize(rec["eff"], rec["oth"], rec["beta"], rec["se"],
                                  rec["p"], rec["freq"], te, to, tfreq)
                    n_c = rec.get("n")
                    if h["status"] not in ("ok", "flipped", "strand_ok",
                                           "strand_flipped"):
                        harmon_audit.append(dict(cell_type=ct, lead=lead["normalized"],
                            target=coh, effect_allele=te, other_allele=to,
                            status=h["status"], notes="cohort-level"))
                b = h["beta"]
                effects_long.append(dict(cell_type=ct, lead_variant=lead["normalized"],
                    rsID=lead.get("rsid", ""), stratum=coh, level="cohort",
                    beta=b, se=h["se"],
                    ci_lo=(b - 1.96 * h["se"]) if (b is not None and h["se"]) else None,
                    ci_hi=(b + 1.96 * h["se"]) if (b is not None and h["se"]) else None,
                    p=h["p"], freq_effect=h["freq"], status=h["status"]))
                panelB_rows.append(dict(cell_type=ct, lead_variant=lead["normalized"],
                    rsID=lead.get("rsid", ""), cohort=coh,
                    cohort_ancestry=COHORT_ANCESTRY[coh], beta=b, status=h["status"]))
                if b is not None and h["se"] and h["status"] in (
                        "ok", "flipped", "strand_ok", "strand_flipped") \
                        or (b is not None and h["se"] and h["status"].startswith(
                            ("ok", "flipped"))):
                    coh_named.append((coh, b, h["se"], n_c))

            # heterogeneity across cohorts
            cc = cross_cohort_stats([(c, b, s) for c, b, s, _ in coh_named],
                                    lead["beta"])
            total_n = sum(n for *_ , n in coh_named if n) or None
            # recomputed FE vs METAL
            fe_vs_metal = None
            if cc.get("k", 0) >= 1 and lead["beta"] is not None:
                fe_vs_metal = abs(cc["fe_beta"] - lead["beta"])
            # cross-ancestry I2 (EUR, AFR, AMR estimates)
            anc_est = []
            for nm in ("EUR", "AFR", "AMR"):
                s = strat[nm]
                if s["beta"] is not None and s["se"]:
                    anc_est.append((s["beta"], s["se"]))
            anc_i2, anc_Q, anc_hetp, anc_k = cross_ancestry_i2(anc_est)

            lead_annot.append(dict(
                cell_type=ct, broad_class=lead["broad_class"],
                lead_variant=lead["normalized"], rsID=lead.get("rsid", ""),
                rsid_status=lead.get("rsid_status", ""),
                shared_locus_id="L%03d" % lead_to_locus.get((ct, lead["normalized"]), 0),
                variant_type=lead["vtype"],
                pooled_beta=lead["beta"], pooled_se=lead["se"], pooled_p=lead["p"],
                pooled_ci_lo=lead["beta"] - 1.96 * lead["se"],
                pooled_ci_hi=lead["beta"] + 1.96 * lead["se"],
                metal_het_isq=lead["het_isq"], metal_het_p=lead["het_p"],
                n_contributing_cohorts=lead["n_cohorts"], total_N=total_n,
                cohort_k_extracted=cc.get("k", 0),
                cohort_i2=cc.get("i2"), cohort_het_p=cc.get("het_p"),
                pct_cohorts_concordant=cc.get("pct_concordant"),
                loo_sign_stable=cc.get("loo_sign_stable"),
                loo_beta_min=cc.get("loo_beta_min"), loo_beta_max=cc.get("loo_beta_max"),
                re_beta=cc.get("re_beta"), re_se=cc.get("re_se"), re_p=cc.get("re_p"),
                fe_recomputed_beta=cc.get("fe_beta"),
                fe_vs_metal_abs_diff=fe_vs_metal,
                ancestry_i2=anc_i2, ancestry_het_p=anc_hetp,
                n_ancestry_groups=anc_k,
                cross_ancestry_assessable=("yes" if anc_k >= 2 else "no"),
                ld_panels_available=lead.get("ld_panels_available", "none"),
                clumping_method=lead.get("clumping_method", ""),
                distance_fallback_used=lead.get("distance_fallback_used", False),
                n_gws_in_clump=lead.get("n_gws_in_clump", 1)))

            panelA_rows.append(dict(
                cell_type=ct, broad_class=lead["broad_class"],
                lead_variant=lead["normalized"], rsID=lead.get("rsid", ""),
                pooled_beta=strat["Pooled"]["beta"], pooled_se=strat["Pooled"]["se"],
                EUR_beta=strat["EUR"]["beta"], EUR_se=strat["EUR"]["se"],
                AFR_beta=strat["AFR"]["beta"], AFR_se=strat["AFR"]["se"],
                AMR_beta=strat["AMR"]["beta"], AMR_se=strat["AMR"]["se"],
                ancestry_i2=anc_i2, ancestry_het_p=anc_hetp,
                n_ancestry_groups=anc_k,
                cross_ancestry_assessable=("yes" if anc_k >= 2 else "no")))

    # ── ordering for figure rows ────────────────────────────────────────────
    def row_sort_key(r):
        return (CLASS_ORDER.get(r["broad_class"], 9),
                CT_RANK.get(r["cell_type"], 99), r["lead_variant"])
    lead_annot.sort(key=row_sort_key)
    for i, r in enumerate(lead_annot):
        r["row_order"] = i
    order_lookup = {(r["cell_type"], r["lead_variant"]): r["row_order"]
                    for r in lead_annot}
    for r in panelA_rows:
        r["row_order"] = order_lookup.get((r["cell_type"], r["lead_variant"]), 999)
    for r in panelB_rows:
        r["row_order"] = order_lookup.get((r["cell_type"], r["lead_variant"]), 999)
    panelA_rows.sort(key=lambda r: r["row_order"])
    panelB_rows.sort(key=lambda r: (r["row_order"], COHORTS_15.index(r["cohort"])
                                    if r["cohort"] in COHORTS_15 else 99))

    # lead-selection audit table
    sel_rows = []
    for lead in sorted(all_leads, key=lambda l: (CLASS_ORDER.get(l["broad_class"], 9),
                       CT_RANK.get(l["cell_type"], 99), l["p"])):
        m = next((a for a in lead_annot if a["cell_type"] == lead["cell_type"]
                  and a["lead_variant"] == lead["normalized"]), {})
        sel_rows.append(dict(
            cell_type=lead["cell_type"], broad_cell_class=lead["broad_class"],
            lead_variant=lead["normalized"], chromosome=lead["chrom"],
            position=lead["pos"], REF=lead["ref"], ALT=lead["alt"],
            variant_type=lead["vtype"], rsID=lead.get("rsid", ""),
            pooled_beta=lead["beta"], pooled_se=lead["se"], pooled_p=lead["p"],
            total_N=m.get("total_N"),
            number_of_contributing_cohorts=lead["n_cohorts"],
            number_of_GWS_variants_in_clump=lead.get("n_gws_in_clump", 1),
            clump_start=lead.get("clump_start"), clump_end=lead.get("clump_end"),
            nearest_secondary_GWS_variant=lead.get("nearest_secondary", ""),
            EUR_r2="", AFR_r2="", AMR_LAT_r2="",
            maximum_r2="", LD_panels_available=lead.get("ld_panels_available", "none"),
            clumping_method=lead.get("clumping_method", ""),
            distance_fallback_used=lead.get("distance_fallback_used", False),
            QC_notes=("indel_lead" if lead["vtype"] != "SNP" else "")))

    # ── write outputs ───────────────────────────────────────────────────────
    write_tsv(OUT_DIR / "variant_normalization_audit.tsv", norm_audit,
              ["cell_type", "original_marker", "normalized_variant", "chromosome",
               "position", "REF", "ALT", "variant_type", "genome_build",
               "normalization_status", "notes"])
    write_tsv(OUT_DIR / "ld_based_lead_selection.tsv", sel_rows,
              ["cell_type", "broad_cell_class", "lead_variant", "chromosome",
               "position", "REF", "ALT", "variant_type", "rsID", "pooled_beta",
               "pooled_se", "pooled_p", "total_N", "number_of_contributing_cohorts",
               "number_of_GWS_variants_in_clump", "clump_start", "clump_end",
               "nearest_secondary_GWS_variant", "EUR_r2", "AFR_r2", "AMR_LAT_r2",
               "maximum_r2", "LD_panels_available", "clumping_method",
               "distance_fallback_used", "QC_notes"])
    write_tsv(OUT_DIR / "ld_clump_membership.tsv", membership_rows,
              ["cell_type", "gws_variant", "assigned_lead", "pooled_p", "EUR_r2",
               "AFR_r2", "maximum_r2", "assignment_method", "distance_fallback_used"])
    write_tsv(OUT_DIR / "cross_celltype_locus_map.tsv", locus_rows,
              ["shared_locus_id", "cell_type", "lead_variant", "rsID", "chromosome",
               "position", "locus_start", "locus_end", "relationship_to_other_leads",
               "LD_evidence", "notes"])
    write_tsv(OUT_DIR / "rsid_mapping_audit.tsv", rsid_audit,
              ["normalized_variant", "chromosome", "position", "REF", "ALT", "rsID",
               "mapping_status", "mapping_source", "dbSNP_build", "genome_build",
               "notes"])
    write_tsv(OUT_DIR / "allele_harmonization_audit.tsv", harmon_audit,
              ["cell_type", "lead", "target", "effect_allele", "other_allele",
               "status", "notes"])
    write_tsv(OUT_DIR / "lead_effects_long.tsv", effects_long,
              ["cell_type", "lead_variant", "rsID", "stratum", "level", "beta",
               "se", "ci_lo", "ci_hi", "p", "freq_effect", "status"])
    write_tsv(OUT_DIR / "lead_heterogeneity_sensitivity.tsv", lead_annot,
              list(lead_annot[0].keys()) if lead_annot else ["cell_type"])
    write_tsv(OUT_DIR / "generalizability_panelA.tsv", panelA_rows,
              ["row_order", "cell_type", "broad_class", "lead_variant", "rsID",
               "pooled_beta", "pooled_se", "EUR_beta", "EUR_se", "AFR_beta",
               "AFR_se", "AMR_beta", "AMR_se", "ancestry_i2", "ancestry_het_p",
               "n_ancestry_groups", "cross_ancestry_assessable"])
    write_tsv(OUT_DIR / "generalizability_panelB_matrix.tsv", panelB_rows,
              ["row_order", "cell_type", "lead_variant", "rsID", "cohort",
               "cohort_ancestry", "beta", "status"])

    # ── QC report ───────────────────────────────────────────────────────────
    n_assoc = len(all_leads)
    n_cts = len({l["cell_type"] for l in all_leads})
    n_snp = sum(1 for l in all_leads if l["vtype"] == "SNP")
    n_indel = n_assoc - n_snp
    n_ld = sum(1 for m in membership_rows if m["assignment_method"] == "LD")
    n_fb = sum(1 for m in membership_rows if m["assignment_method"] == "distance_fallback")
    multi = defaultdict(int)
    for l in all_leads:
        multi[l["cell_type"]] += 1
    multi_ct = {c: n for c, n in multi.items() if n > 1}
    zero_ct = [c for c in cell_types if multi.get(c, 0) == 0]
    shared_loci = len({r["shared_locus_id"] for r in locus_rows
                       if r["relationship_to_other_leads"] == "shared_multi_celltype"})
    no_rsid = [a["normalized_variant"] for a in rsid_audit
               if a["mapping_status"] in ("no_rsID", "allele_mismatch", "unresolved")]
    allele_mm = [a for a in harmon_audit if a["status"] == "allele_mismatch"]
    palindromic_excl = [a for a in harmon_audit
                        if a["status"].startswith("excluded_palindromic")]
    single_anc = [a for a in lead_annot if a["cross_ancestry_assessable"] == "no"]
    under_cohort = [a for a in lead_annot
                    if (a["n_contributing_cohorts"] or 0) < cfg.min_cohorts]
    fe_diffs = [a["fe_vs_metal_abs_diff"] for a in lead_annot
                if a["fe_vs_metal_abs_diff"] is not None]
    max_fe_diff = max(fe_diffs) if fe_diffs else None
    nonfinite_ci = [a for a in lead_annot
                    if not (math.isfinite(a["pooled_ci_lo"]) and
                            math.isfinite(a["pooled_ci_hi"]))]
    dup_rows = len(effects_long) - len({(e["cell_type"], e["lead_variant"], e["stratum"])
                                        for e in effects_long})

    qc = []
    qc.append("Generalizability lead-selection QC report")
    qc.append("=" * 60)
    qc.append(f"Genome build: {GENOME_BUILD}")
    qc.append(f"Config: p<={cfg.p_thresh:g}, clump_r2>={cfg.clump_r2}, "
              f"clump_kb={cfg.clump_kb}, min_cohorts={cfg.min_cohorts}, "
              f"LD={'on' if cfg.use_ld else 'off'}, reconciliation={cfg.reconciliation}")
    qc.append("")
    qc.append("LD panels used:")
    qc.append(f"  EUR: 1000G Phase3 EUR (per-chrom bed, hg19)")
    qc.append(f"  AFR: in-sample NIMH_HBCC_1M_AFR + NIMH_HBCC_h650_AFR (merged bed)")
    qc.append(f"  AMR: no hg19 in-sample panel (AMP_AD genotypes are hg38) -> distance fallback")
    qc.append("")
    qc.append(f"Lead SNP-cell-type associations: {n_assoc}")
    qc.append(f"Previous expected count: 19  ->  new count: {n_assoc}")
    if n_assoc != 19:
        qc.append("  Count changed vs prior distance+SNP-preference pipeline because: "
                  "(a) LD r2>=0.1 can split nearby independent signals that pure "
                  "distance clumping merged; (b) no SNP-over-indel replacement; "
                  "(c) min_cohorts filter. See ld_clump_membership.tsv.")
    qc.append(f"Contributing cell types: {n_cts}; cell types with no GW lead: "
              f"{len(zero_ct)} ({', '.join(zero_ct) if zero_ct else 'none'})")
    qc.append(f"Cell types with >1 lead: {multi_ct if multi_ct else 'none'}")
    qc.append(f"Distinct cross-cell-type genomic regions: {n_loci}; "
              f"shared by >1 cell type: {shared_loci}")
    qc.append(f"SNP leads: {n_snp}; indel leads: {n_indel}")
    qc.append(f"Clump assignments via LD: {n_ld}; via distance fallback: {n_fb}")
    qc.append("")
    qc.append(f"rsID mapped (exact/multiple): {n_assoc - len(no_rsid)}/{n_assoc}; "
              f"without reliable rsID: {len(no_rsid)}")
    qc.append(f"Allele mismatches remaining: {len(allele_mm)}; "
              f"palindromic excluded: {len(palindromic_excl)}")
    qc.append(f"Leads without cross-ancestry evaluation (<2 ancestries): {len(single_anc)}")
    qc.append(f"Leads under min cohorts ({cfg.min_cohorts}): {len(under_cohort)}")
    qc.append(f"Duplicate cell-type/variant/stratum rows: {dup_rows}")
    qc.append(f"Non-finite pooled CIs: {len(nonfinite_ci)}")
    qc.append(f"Max |recomputed FE beta - METAL beta|: "
              f"{max_fe_diff:.4g}" if max_fe_diff is not None else
              "Max |recomputed FE - METAL|: n/a")
    qc.append("")
    if warnings:
        qc.append("WARNINGS requiring manual review:")
        for w in warnings:
            qc.append(f"  - {w}")
    else:
        qc.append("No warnings.")
    (OUT_DIR / "generalizability_figure_QC_report.txt").write_text("\n".join(qc) + "\n")

    # ── console summary ─────────────────────────────────────────────────────
    print("\n" + "\n".join(qc[-25:] if len(qc) > 25 else qc))
    print(f"\nOutputs written under {OUT_DIR}")


if __name__ == "__main__":
    main()
