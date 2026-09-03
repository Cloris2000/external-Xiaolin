#!/usr/bin/env python3
"""Resolve chr:pos:ref:alt lead markers (hg19) to dbSNP rsIDs via Ensembl GRCh37 REST.
Writes a cache TSV: marker<TAB>rsid  (rsid = '.' if none matched)."""
import csv, json, time, urllib.request
from pathlib import Path

ROOT = Path("/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow")
LEADS = ROOT / "results/meta_sensitivity/ancestry_lead_effects/ancestry_leads.tsv"
OUT   = ROOT / "results/meta_sensitivity/ancestry_lead_effects/lead_rsids.tsv"

def norm(a):  # normalise allele string for comparison
    return a.upper()

def query(chrom, pos):
    url = (f"https://grch37.rest.ensembl.org/overlap/region/human/"
           f"{chrom}:{pos}-{pos}?feature=variation;content-type=application/json")
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "rsid-lookup"})
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            time.sleep(2 * (attempt + 1))
    return []

markers = {}
for r in csv.DictReader(open(LEADS), delimiter='\t'):
    markers[r['marker']] = (r['chrom'], int(r['pos']),
                            norm(r['effect_allele']), norm(r['other_allele']))

cache = {}
for mk, (chrom, pos, ea, oa) in markers.items():
    hits = query(chrom, pos)
    rsid = '.'
    best = None
    for h in hits:
        if not str(h.get('id', '')).startswith('rs'):
            continue
        alleles = {norm(a) for a in h.get('alleles', [])}
        if h.get('start') != pos:
            continue
        # prefer variant whose allele set contains both our alleles
        if ea in alleles and oa in alleles:
            best = h['id']; break
        if best is None and (ea in alleles or oa in alleles):
            best = h['id']
        if best is None:
            best = h['id']  # fallback: any rs at this exact position
    if best:
        rsid = best
    cache[mk] = rsid
    print(f"{mk}\t{rsid}", flush=True)
    time.sleep(0.15)

with open(OUT, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['marker', 'rsid'])
    for mk, rs in cache.items():
        w.writerow([mk, rs])
print(f"\nWrote {OUT}")
