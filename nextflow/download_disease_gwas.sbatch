#!/bin/bash
#SBATCH --job-name=dl_disease_gwas
#SBATCH --array=1-8%8              # 8 datasets, all at once
#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/coloc/dl_gwas_%A_%a.out
#SBATCH --error=logs/coloc/dl_gwas_%A_%a.err

# Download all disease GWAS in parallel (one SLURM task per dataset).
# Usage: sbatch download_disease_gwas.sbatch
# After all tasks complete, check logs and then run: sbatch run_coloc.sbatch

set -euo pipefail

PROJECT_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow"
OUTDIR="$PROJECT_DIR/results/coloc/disease_gwas"
mkdir -p logs/coloc "$OUTDIR"

export PATH="/nethome/kcni/xzhou/.anaconda3/bin:$PATH"

# ---------------------------------------------------------------------------
# Dataset definitions: NAME | RAW_FILE | URL | N | Ncase | Ncont | FORMAT
# ---------------------------------------------------------------------------
# Formats: harmonised (GWAS Catalog harmonised TSV), pgc (CHR/BP/A1/A2/BETA/SE/P),
#          pgc_or (CHR/BP/A1/A2/OR/SE/P), mdd (CHR/BP/A1/A2/OR/SE/P with Nca/Nco)

declare -A NAMES RAWS URLS NS NCASES NCONTS FMTS

NAMES[1]="AD_Bellenguez2022"
RAWS[1]="raw.tsv.gz"
URLS[1]="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/harmonised/35379992-GCST90027158-MONDO_0004975.h.tsv.gz"
NS[1]=788989; NCASES[1]=111326; NCONTS[1]=677663; FMTS[1]="harmonised"

NAMES[2]="PD_Nalls2019"
RAWS[2]="raw.tsv.gz"
URLS[2]="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009325/harmonised/GCST009325.h.tsv.gz"
NS[2]=1474097; NCASES[2]=37688; NCONTS[2]=1436409; FMTS[2]="harmonised_pd"

NAMES[3]="SCZ_Trubetskoy2022"
RAWS[3]="raw.tsv.gz"
URLS[3]="https://ndownloader.figshare.com/files/34517828"
NS[3]=161405; NCASES[3]=67390; NCONTS[3]=94015; FMTS[3]="pgc_vcftsv"

NAMES[4]="BD_Mullins2021"
RAWS[4]="raw.tsv.gz"
URLS[4]="https://ndownloader.figshare.com/files/26603681"
NS[4]=413466; NCASES[4]=41917; NCONTS[4]=371549; FMTS[4]="pgc_vcftsv"

NAMES[5]="MDD_Wray2018"
RAWS[5]="raw.gz"
URLS[5]="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005839/MDD2018_ex23andMe.gz"
NS[5]=172953; NCASES[5]=59851; NCONTS[5]=113154; FMTS[5]="mdd_or"

NAMES[6]="LBD_Chia2021"
RAWS[6]="raw.tsv.gz"
URLS[6]="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90001001-GCST90002000/GCST90001390/harmonised/33589841-GCST90001390-EFO_0006792.h.tsv.gz"
NS[6]=8392; NCASES[6]=2591; NCONTS[6]=5801; FMTS[6]="harmonised"

NAMES[7]="ALS_vanRheenen2021"
RAWS[7]="raw.tsv.gz"
URLS[7]="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90012001-GCST90013000/GCST90012201/harmonised/32915819-GCST90012201-EFO_0005680.h.tsv.gz"
NS[7]=138086; NCASES[7]=27205; NCONTS[7]=110881; FMTS[7]="harmonised"

# Task 8: FTD — skip (manual download required); just write README
NAMES[8]="FTD_vanderLee2019"
RAWS[8]=""; URLS[8]=""; NS[8]=0; NCASES[8]=0; NCONTS[8]=0; FMTS[8]="manual"

# ---------------------------------------------------------------------------
IDX=${SLURM_ARRAY_TASK_ID}
NAME="${NAMES[$IDX]}"
RAW="${RAWS[$IDX]}"
URL="${URLS[$IDX]}"
N="${NS[$IDX]}"
NCASE="${NCASES[$IDX]}"
NCONT="${NCONTS[$IDX]}"
FMT="${FMTS[$IDX]}"
DIR="$OUTDIR/$NAME"

mkdir -p "$DIR"
echo "=== Task $IDX: $NAME ==="
echo "URL: $URL"
echo "Format: $FMT"

# ---- FTD: manual only ----
if [ "$FMT" = "manual" ]; then
  cat > "$DIR/README.txt" <<'FTDREADME'
FTD/FTLD GWAS - Manual Download Required
==========================================
Recommended: van der Lee et al. 2019 (GCST009258).
Request access from: https://www.niagads.org/ or contact the FTD consortium.
Once downloaded, standardize to:
  snp (chr:pos:REF:ALT), chr, pos, ref, alt, beta, se, p, N, Ncase, Ncont
and save as: FTD_vanderLee2019_hg19.tsv in this directory.
FTDREADME
  echo "FTD: manual download required. README written."
  exit 0
fi

RAW_PATH="$DIR/$RAW"

# ---- Download (skip if already non-empty) ----
echo "Downloading (resuming if partial)..."
if [[ "$URL" == *figshare* || "$URL" == *ndownloader* ]]; then
  # figshare: pre-signed URLs expire in 10s; always start fresh (no resume)
  [ -s "$RAW_PATH" ] || \
    wget -q -L --max-redirect=5 --tries=5 --timeout=600 --waitretry=30 -O "$RAW_PATH" "$URL"
else
  # EBI FTP and others: use -c to resume partial downloads
  wget -q -c -L --max-redirect=5 --tries=5 --timeout=600 --waitretry=30 -O "$RAW_PATH" "$URL"
fi

if [ ! -s "$RAW_PATH" ]; then
  echo "ERROR: Download failed or file is empty: $RAW_PATH"
  echo "Manual download: $URL"
  exit 1
fi

# Verify gzip integrity before standardizing (catches truncated downloads)
if [[ "$RAW_PATH" == *.gz ]]; then
  if ! gzip -t "$RAW_PATH" 2>/dev/null; then
    echo "WARN: $RAW_PATH is truncated — retrying download..."
    rm -f "$RAW_PATH"
    wget -q -L --max-redirect=5 --tries=5 --timeout=1800 --waitretry=30 -O "$RAW_PATH" "$URL"
    if ! gzip -t "$RAW_PATH" 2>/dev/null; then
      echo "ERROR: Download still truncated after retry. Aborting."
      exit 1
    fi
  fi
fi

SIZE=$(du -sh "$RAW_PATH" | cut -f1)
echo "Downloaded: $RAW_PATH ($SIZE)"

# ---- Standardize ----
OUT_TSV="$DIR/${NAME}_hg19.tsv"

python3 - <<PYEOF
import gzip, csv, math, sys

infile  = "$RAW_PATH"
outfile = "$OUT_TSV"
fmt     = "$FMT"
N, Ncase, Ncont = $N, $NCASE, $NCONT

opener = gzip.open if infile.endswith(".gz") else open
n_written = 0

with opener(infile, "rt") as fin, open(outfile, "w") as fout:
    fout.write("snp\tchr\tpos\tref\talt\tbeta\tse\tp\tN\tNcase\tNcont\n")
    # Skip VCF-style ## meta-info lines; the actual header starts with #CHROM or a field name
    header_line = fin.readline()
    while header_line.startswith("##"):
        header_line = fin.readline()
    # Strip leading '#' from header (e.g. #CHROM -> CHROM)
    header_line = header_line.lstrip("#").rstrip("\n")
    fieldnames = header_line.split("\t")
    reader = csv.DictReader(fin, fieldnames=fieldnames, delimiter="\t")

    for row in reader:
        try:
            if fmt in ("harmonised", "harmonised_pd"):
                chrom = str(row.get("chromosome", "")).lstrip("chr")
                pos   = int(row.get("base_pair_location", ""))
                ref   = row.get("effect_allele", "").upper()
                alt   = row.get("other_allele", "").upper()
                beta  = row.get("beta", "")
                se    = row.get("standard_error", "")
                pval  = row.get("p_value", "")
                # PD harmonised file has per-SNP N_cases/N_controls
                if fmt == "harmonised_pd":
                    nc  = row.get("N_cases", "")
                    nct = row.get("N_controls", "")
                    if nc and nct:
                        try:
                            Ncase = int(float(nc))
                            Ncont = int(float(nct))
                            N     = Ncase + Ncont
                        except ValueError:
                            pass

            elif fmt == "pgc":
                chrom = str(row.get("CHR", row.get("CHROM", ""))).lstrip("chr")
                pos   = int(row.get("BP", row.get("POS", "")))
                ref   = row.get("A1", row.get("REF", "")).upper()
                alt   = row.get("A2", row.get("ALT", "")).upper()
                beta  = row.get("BETA", "")
                if not beta:
                    or_ = float(row.get("OR", 0))
                    beta = str(math.log(or_)) if or_ > 0 else ""
                se    = row.get("SE", "")
                pval  = row.get("P", "")

            elif fmt == "pgc_vcftsv":
                # PGC vcf.tsv.gz format: CHROM (or #CHROM) POS A1 A2 BETA SE PVAL NCAS NCON
                chrom = str(row.get("CHROM", row.get("#CHROM", ""))).lstrip("chr")
                pos   = int(row.get("POS", "") or 0)
                ref   = (row.get("A1") or "").upper()
                alt   = (row.get("A2") or "").upper()
                beta  = row.get("BETA") or ""
                se    = row.get("SE") or ""
                pval  = row.get("PVAL") or ""
                nc    = row.get("NCAS") or ""
                nct   = row.get("NCON") or ""
                if nc and nct:
                    try:
                        Ncase = int(float(nc))
                        Ncont = int(float(nct))
                        N     = Ncase + Ncont
                    except (ValueError, TypeError):
                        pass

            elif fmt == "mdd_or":
                # PGC MDD: CHR SNP BP A1 A2 OR SE P Nca Nco
                chrom = str(row.get("CHR") or "").lstrip("chr")
                pos   = int(row.get("BP") or 0)
                ref   = (row.get("A1") or "").upper()
                alt   = (row.get("A2") or "").upper()
                or_   = float(row.get("OR") or 0)
                beta  = str(math.log(or_)) if or_ > 0 else ""
                se    = row.get("SE") or ""
                pval  = row.get("P") or ""
                nc    = row.get("Nca") or ""
                nct   = row.get("Nco") or ""
                if nc and nct:
                    try:
                        Ncase = int(float(nc))
                        Ncont = int(float(nct))
                        N     = Ncase + Ncont
                    except (ValueError, TypeError):
                        pass

            else:
                continue

            if not all([chrom, pos, ref, alt, beta, se, pval]):
                continue

            snp_id = f"chr{chrom}:{pos}:{ref}:{alt}"
            fout.write(f"{snp_id}\t{chrom}\t{pos}\t{ref}\t{alt}\t{beta}\t{se}\t{pval}\t{N}\t{Ncase}\t{Ncont}\n")
            n_written += 1

        except (IndexError, ValueError, KeyError):
            continue

print(f"Written {n_written:,} variants -> {outfile}")
if n_written == 0:
    sys.exit(1)
PYEOF

echo "=== Done: $NAME ==="
