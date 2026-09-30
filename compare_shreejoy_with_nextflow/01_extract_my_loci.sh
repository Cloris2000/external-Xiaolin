#!/bin/bash
# Stage A / step 1
#
# Pull the TMEM106B and GRN windows out of the 15-cohort METAL meta-analysis.
#
# Reads the .annotated.tsv rather than the raw .tbl: METAL lowercases alleles and
# reports Effect against its own Allele1, which is not consistently REF or ALT.
# The annotated file carries BETA_ALT / FREQ_ALT, already oriented to ALT, which
# is the same orientation as REGENIE's ALLELE1 in Shreejoy's ID scheme
# (chr:pos:REF:ALT with ALLELE0=REF, ALLELE1=ALT). Without this the effect signs
# would not be comparable between the two pipelines.
#
# Read-only. Roughly 5 minutes of sequential I/O over ~25 GB.

set -euo pipefail
source "$(dirname "$0")/config.sh"

OUT="${DATA_DIR}/my_loci.tsv"
TMP="${OUT}.partial"

printf 'trait\tvariant\tCHROM\tPOS_hg19\tREF\tALT\tBETA_ALT\tSE\tP\tFREQ_ALT\tHetISq\tHetPVal\tDIRECTION_ALT\n' > "${TMP}"

n_files=0
for f in "${MINE_META}"/*.annotated.tsv; do
    trait="$(basename "$f")"
    trait="${trait%%_meta_analysis_*}"
    n_files=$((n_files + 1))
    echo "  [${n_files}] ${trait}" >&2

    # Columns in .annotated.tsv:
    #   9 StdErr  10 P-value  12 HetISq  15 HetPVal
    #  16 CHROM  17 POS  18 REF  19 ALT  23 BETA_ALT  24 FREQ_ALT  27 DIRECTION_ALT
    awk -F'\t' -v OFS='\t' -v trait="${trait}" '
        NR == 1 {
            # Fail loudly if the column layout is not what this script assumes.
            if ($16 != "CHROM" || $17 != "POS" || $23 != "BETA_ALT" || $24 != "FREQ_ALT")
                { print "FATAL: unexpected column layout in " FILENAME > "/dev/stderr"; exit 1 }
            next
        }
        ($16 == 7  && $17 >= 11900000 && $17 <= 12600000) ||
        ($16 == 17 && $17 >= 42100000 && $17 <= 42800000) {
            print trait, $1, $16, $17, $18, $19, $23, $9, $10, $24, $12, $15, $27
        }
    ' "$f" >> "${TMP}"
done

mv "${TMP}" "${OUT}"
echo "wrote ${OUT}: $(( $(wc -l < "${OUT}") - 1 )) rows from ${n_files} traits" >&2
