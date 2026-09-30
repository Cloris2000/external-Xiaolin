#!/usr/bin/env python3
"""Stage A / step 4

Put Xiaolin's hg19 locus slices onto hg38 so the two pipelines can be plotted on
one axis and joined variant by variant.

There is no liftover chain file and no pyliftover on this system, and installing
either needs approval. Inside a window this small, with no assembly change
between builds, the hg19->hg38 map is a pure translation, so the offset is
recovered from the data instead:

  1. Index both slices by (REF, ALT) -- these are build-invariant for a SNV.
  2. Vote for every candidate offset implied by a REF/ALT match whose allele
     frequencies also agree. There are only four SNV allele combinations, so
     without the frequency condition each variant would vote against thousands
     of others and bury the true offset in noise.
  3. Take the modal offset, then require that it accounts for nearly all of my
     variants that have any counterpart at all, and that it beats the runner-up
     offset by a wide margin.

If the modal offset does not dominate, the window is not a clean translation and
the script fails rather than emitting a silently wrong join.

Writes:
  data/coord_offsets.tsv   the derived offset per chromosome, with evidence
  data/my_loci_hg38.tsv    my_loci.tsv with a POS_hg38 column added

Read-only against the pipelines. Stdlib only.
"""

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data"

# Share of my matchable variants the winning offset must account for. This is
# deliberately not near 1.0: the two pipelines have different variant universes
# (8.58M against 9.52M), so some of my variants have no counterpart at all, and
# with only four SNV allele classes in a 700 kb window a fair number of those
# will still find a same-alleles, similar-frequency partner at a random offset.
MIN_DOMINANCE = 0.75
# The winning offset must beat the runner-up by at least this factor. This is
# the statistic that actually discriminates: a genuine translation puts
# thousands of variants on one offset and a handful on each spurious one.
MIN_SHARPNESS = 50.0
# Allele frequency must agree to this tolerance for a pair to count as a match.
MAX_FREQ_DIFF = 0.05
# Builds differ by a few Mb at most in these windows; ignore wilder candidates.
MAX_ABS_OFFSET = 5_000_000

# Independent anchor. Gene spans from the build annotations, used only to check
# that the empirically derived offset lands where the assemblies say it should.
# TMEM106B  hg19 chr7:12,250,867   hg38 chr7:12,211,240
# GRN       hg19 chr17:42,422,491  hg38 chr17:44,345,086
EXPECTED_OFFSET = {7: 12_211_240 - 12_250_867, 17: 44_345_086 - 42_422_491}
MAX_ANCHOR_DEVIATION = 2_000


def load_mine() -> list[dict]:
    with (DATA / "my_loci.tsv").open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_his_positions() -> dict[int, dict[tuple[str, str], list[tuple[int, float]]]]:
    """chrom -> (REF, ALT) -> [(pos, freq), ...], from the pooled arm only."""
    out: dict[int, dict[tuple[str, str], list[tuple[int, float]]]] = defaultdict(
        lambda: defaultdict(list))
    seen: set[tuple[int, int]] = set()
    with (DATA / "his_loci.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["arm"] != "pooled":
                continue
            chrom, pos = int(r["CHROM"]), int(r["POS_hg38"])
            if (chrom, pos) in seen:
                continue           # one row per variant, not per trait
            seen.add((chrom, pos))
            try:
                freq = float(r["A1FREQ"])
            except ValueError:
                continue
            out[chrom][(r["REF"].upper(), r["ALT"].upper())].append((pos, freq))
    return out


def derive_offset(mine: list[dict], his: dict, chrom: int) -> tuple[int, dict]:
    # One row per variant on my side; my_loci.tsv repeats each variant per trait.
    rows: dict[int, dict] = {}
    for r in mine:
        if int(r["CHROM"]) == chrom:
            rows.setdefault(int(r["POS_hg19"]), r)

    votes: Counter[int] = Counter()
    matchable = 0
    for pos19, r in rows.items():
        key = (r["REF"].upper(), r["ALT"].upper())
        try:
            freq19 = float(r["FREQ_ALT"])
        except ValueError:
            continue
        hit = False
        for pos38, freq38 in his[chrom].get(key, ()):
            off = pos38 - pos19
            if abs(off) > MAX_ABS_OFFSET:
                continue
            if abs(freq19 - freq38) > MAX_FREQ_DIFF:
                continue
            votes[off] += 1
            hit = True
        matchable += hit

    if not votes:
        sys.exit(f"FATAL: chr{chrom} has no REF/ALT matches between builds")

    ranked = votes.most_common(2)
    offset, n_modal = ranked[0]
    runner_up, n_runner = ranked[1] if len(ranked) > 1 else ("NA", 0)

    dominance = n_modal / matchable if matchable else 0.0
    sharpness = n_modal / n_runner if n_runner else float("inf")

    anchor = EXPECTED_OFFSET.get(chrom)
    return offset, {
        "chrom": chrom,
        "offset": offset,
        "n_supporting": n_modal,
        "n_matchable_variants": matchable,
        "dominance": round(dominance, 4),
        "runner_up_offset": runner_up,
        "runner_up_n": n_runner,
        "sharpness": round(sharpness, 1) if n_runner else "inf",
        "expected_offset": anchor if anchor is not None else "NA",
        "anchor_deviation_bp": (abs(offset - anchor)
                                if anchor is not None else "NA"),
    }


def validate(ev: dict) -> list[str]:
    """Reasons this offset should not be trusted. Empty means it is fine."""
    chrom, offset = ev["chrom"], ev["offset"]
    problems = []
    if ev["dominance"] < MIN_DOMINANCE:
        problems.append(
            f"chr{chrom} offset {offset} accounts for only "
            f"{ev['dominance']:.1%} of matchable variants (floor "
            f"{MIN_DOMINANCE:.0%})")
    if ev["sharpness"] != "inf" and ev["sharpness"] < MIN_SHARPNESS:
        problems.append(
            f"chr{chrom} offset {offset} beats runner-up "
            f"{ev['runner_up_offset']} by only {ev['sharpness']}x "
            f"(floor {MIN_SHARPNESS:.0f}x)")
    if ev["anchor_deviation_bp"] != "NA" and \
            ev["anchor_deviation_bp"] > MAX_ANCHOR_DEVIATION:
        problems.append(
            f"chr{chrom} offset {offset} is {ev['anchor_deviation_bp']} bp from "
            f"the gene-annotation offset {ev['expected_offset']}")
    return problems


def main() -> int:
    mine = load_mine()
    his = load_his_positions()

    offsets: dict[int, int] = {}
    evidence: list[dict] = []
    for chrom in sorted({int(r["CHROM"]) for r in mine}):
        off, ev = derive_offset(mine, his, chrom)
        offsets[chrom] = off
        evidence.append(ev)
        print(f"  chr{chrom}: hg19 -> hg38 offset {off:+d} "
              f"({ev['dominance']:.1%} of {ev['n_matchable_variants']} "
              f"matchable variants, {ev['sharpness']}x the runner-up, "
              f"{ev['anchor_deviation_bp']} bp from the gene annotation)",
              file=sys.stderr)

    # Always record the evidence, including when it fails, so a rejected offset
    # can be inspected rather than guessed at.
    with (DATA / "coord_offsets.tsv").open("w") as out:
        w = csv.DictWriter(out, fieldnames=list(evidence[0]), delimiter="\t")
        w.writeheader()
        w.writerows(evidence)

    problems = [p for ev in evidence for p in validate(ev)]
    if problems:
        sys.exit("FATAL: coordinate alignment rejected\n  "
                 + "\n  ".join(problems)
                 + f"\nEvidence written to {DATA}/coord_offsets.tsv")

    src_fields = list(mine[0])
    with (DATA / "my_loci_hg38.tsv").open("w") as out:
        w = csv.DictWriter(out, fieldnames=src_fields + ["POS_hg38"],
                           delimiter="\t")
        w.writeheader()
        for r in mine:
            r["POS_hg38"] = int(r["POS_hg19"]) + offsets[int(r["CHROM"])]
            w.writerow(r)

    print(f"wrote {DATA}/coord_offsets.tsv and my_loci_hg38.tsv "
          f"({len(mine)} rows)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
