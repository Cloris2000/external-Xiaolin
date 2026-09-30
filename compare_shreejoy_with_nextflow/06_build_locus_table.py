#!/usr/bin/env python3
"""Stage A / step 6

Write the definition of all 13 of Shreejoy's loci, so the extractors and the
figures all read the same windows from one place.

Positions are hg38, taken from the lead variant of his own locus table. My side
has to be extracted on a deliberately wide hg19 window, because without a chain
file the hg19 position is not known until 08_align_all_coords.py derives the
offset. Four megabases is comfortably wider than any hg19-to-hg38 shift on these
chromosomes (the largest here is chr17 at 1.92 Mb) and costs nothing but a
slightly longer awk filter.

Read-only. Stdlib only.
"""

import csv
from collections import defaultdict
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data"
HITS = Path("/scratch/shreejoy/ctpgwas/results/loci/annotated_hits.csv")

PLOT_HALF_WINDOW = 350_000      # what the figures show, hg38
SEARCH_HALF_WINDOW = 4_000_000  # what to pull from my hg19 files


def main() -> int:
    DATA.mkdir(parents=True, exist_ok=True)
    by_locus: dict[str, list[dict]] = defaultdict(list)
    with HITS.open() as fh:
        for r in csv.DictReader(fh):
            by_locus[r["locus"]].append(r)

    rows = []
    for locus, hits in by_locus.items():
        best = min(hits, key=lambda r: float(r["P"]))
        pos = int(float(best["GENPOS"]))
        rows.append({
            "locus": locus,
            "nearest_gene": best["nearest_gene"],
            "CHROM": int(best["CHROM"]),
            "lead_variant": best["variant"],
            "lead_pos_hg38": pos,
            "best_P": best["P"],
            "best_trait": best["trait"],
            "n_arms": len({h["arm"] for h in hits}),
            "n_traits": len({h["trait"] for h in hits}),
            "plot_start_hg38": pos - PLOT_HALF_WINDOW,
            "plot_end_hg38": pos + PLOT_HALF_WINDOW,
            "search_start_hg19": max(1, pos - SEARCH_HALF_WINDOW),
            "search_end_hg19": pos + SEARCH_HALF_WINDOW,
        })

    rows.sort(key=lambda r: float(r["best_P"]))
    out = DATA / "loci.tsv"
    with out.open("w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    chroms = sorted({r["CHROM"] for r in rows})
    print(f"wrote {out}: {len(rows)} loci on {len(chroms)} chromosomes "
          f"{chroms}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
