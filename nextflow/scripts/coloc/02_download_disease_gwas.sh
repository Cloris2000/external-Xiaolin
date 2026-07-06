#!/usr/bin/env bash
# Download and standardize disease GWAS summary statistics (hg19/GRCh37).
#
# Output directory: results/coloc/disease_gwas/
# Each disease gets a standardized TSV:
#   snp (chr:pos:REF:ALT), chr, pos, ref, alt, beta, se, p, N, Ncase, Ncont
#
# Run from the nextflow project root:
#   bash scripts/coloc/02_download_disease_gwas.sh
#
# URLs verified June 2026. All harmonised files use the GWAS Catalog
# harmonised format: chromosome, base_pair_location, effect_allele,
# other_allele, beta, standard_error, p_value, ...

set -euo pipefail

OUTDIR="results/coloc/disease_gwas"
mkdir -p "$OUTDIR"

echo "=== Disease GWAS download script ==="
echo "Output: $OUTDIR"
echo ""

# ---------------------------------------------------------------------------
# Helper: download with wget, following redirects, retrying once
# ---------------------------------------------------------------------------
dl() {
  local url="$1" out="$2"
  wget -q -L --tries=2 --timeout=120 -O "$out" "$url" && \
    [ -s "$out" ] && return 0
  # If file is empty, remove it and report failure
  rm -f "$out"
  echo "  WARN: Download failed for $url"
  return 1
}

# ---------------------------------------------------------------------------
# Helper: standardize GWAS Catalog harmonised format to common TSV
# Harmonised columns: chromosome, base_pair_location, effect_allele,
#   other_allele, beta, standard_error, p_value [, N_cases, N_controls]
# ---------------------------------------------------------------------------
standardize_harmonised() {
  local input="$1" output="$2" N="$3" Ncase="$4" Ncont="$5"
  python3 - <<PYEOF
import gzip, csv, sys

infile  = "$input"
outfile = "$output"
N, Ncase, Ncont = $N, $Ncase, $Ncont

opener = gzip.open if infile.endswith(".gz") else open
with opener(infile, "rt") as fin, open(outfile, "w") as fout:
    fout.write("snp\tchr\tpos\tref\talt\tbeta\tse\tp\tN\tNcase\tNcont\n")
    reader = csv.DictReader(fin, delimiter="\t")
    n_written = 0
    for row in reader:
        try:
            chrom = str(row.get("chromosome", row.get("hm_chrom",""))).lstrip("chr")
            pos   = int(row.get("base_pair_location", row.get("hm_pos","")))
            ref   = row.get("effect_allele", row.get("hm_effect_allele","")).upper()
            alt   = row.get("other_allele",  row.get("hm_other_allele","")).upper()
            beta  = row.get("beta", row.get("hm_beta",""))
            se    = row.get("standard_error","")
            pval  = row.get("p_value","")
            if not all([chrom, pos, ref, alt, beta, se, pval]):
                continue
            snp_id = f"chr{chrom}:{pos}:{ref}:{alt}"
            fout.write(f"{snp_id}\t{chrom}\t{pos}\t{ref}\t{alt}\t{beta}\t{se}\t{pval}\t{N}\t{Ncase}\t{Ncont}\n")
            n_written += 1
        except (IndexError, ValueError, KeyError):
            continue
print(f"  Written {n_written:,} variants -> {outfile}")
PYEOF
}

# ===========================================================================
# 1. ALZHEIMER'S DISEASE — Bellenguez et al. 2022 (Nature Genetics)
#    GWAS Catalog: GCST90027158  (harmonised hg38→hg37 by EBI)
#    N=788,989 (111,326 cases, 677,663 controls), European
#    PMID: 35379992
# ===========================================================================
echo "--- AD: Bellenguez 2022 ---"
AD_DIR="$OUTDIR/AD_Bellenguez2022"
mkdir -p "$AD_DIR"
AD_RAW="$AD_DIR/raw.tsv.gz"
AD_OUT="$AD_DIR/AD_Bellenguez2022_hg19.tsv"
AD_URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/harmonised/35379992-GCST90027158-MONDO_0004975.h.tsv.gz"

if [ ! -s "$AD_RAW" ]; then
  echo "Downloading Bellenguez 2022 AD (~1.1 GB)..."
  dl "$AD_URL" "$AD_RAW" || true
fi
if [ -s "$AD_RAW" ]; then
  standardize_harmonised "$AD_RAW" "$AD_OUT" 788989 111326 677663
  echo "AD Bellenguez2022 written to $AD_OUT"
else
  echo "  SKIP: $AD_RAW not downloaded. Get manually from:"
  echo "  $AD_URL"
fi

# ===========================================================================
# 2. PARKINSON'S DISEASE — Nalls et al. 2019 (Lancet Neurology)
#    GWAS Catalog: GCST009325  (harmonised)
#    N=~1,474,097 (37,688 cases, 1,436,409 controls, proxy-case design)
#    PMID: 31701892
# ===========================================================================
echo "--- PD: Nalls 2019 ---"
PD_DIR="$OUTDIR/PD_Nalls2019"
mkdir -p "$PD_DIR"
PD_RAW="$PD_DIR/raw.tsv.gz"
PD_OUT="$PD_DIR/PD_Nalls2019_hg19.tsv"
PD_URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009325/harmonised/GCST009325.h.tsv.gz"

if [ ! -s "$PD_RAW" ]; then
  echo "Downloading Nalls 2019 PD (~490 MB)..."
  dl "$PD_URL" "$PD_RAW" || true
fi
if [ -s "$PD_RAW" ]; then
  python3 - <<'PDPY'
import gzip, csv
infile  = "results/coloc/disease_gwas/PD_Nalls2019/raw.tsv.gz"
outfile = "results/coloc/disease_gwas/PD_Nalls2019/PD_Nalls2019_hg19.tsv"
Ncase, Ncont = 37688, 1436409
N = Ncase + Ncont
n_written = 0
with gzip.open(infile, "rt") as fin, open(outfile, "w") as fout:
    fout.write("snp\tchr\tpos\tref\talt\tbeta\tse\tp\tN\tNcase\tNcont\n")
    reader = csv.DictReader(fin, delimiter="\t")
    for row in reader:
        try:
            chrom = str(row.get("chromosome","")).lstrip("chr")
            pos   = int(row.get("base_pair_location",""))
            ref   = row.get("effect_allele","").upper()
            alt   = row.get("other_allele","").upper()
            beta  = row.get("beta","")
            se    = row.get("standard_error","")
            pval  = row.get("p_value","")
            if not all([chrom, pos, ref, alt, beta, se, pval]):
                continue
            # Use N_cases/N_controls from file if available, else use totals
            nc = row.get("N_cases","")
            nct = row.get("N_controls","")
            ncase_v = nc if nc else Ncase
            ncont_v = nct if nct else Ncont
            n_v = int(ncase_v) + int(ncont_v) if nc else N
            snp_id = f"chr{chrom}:{pos}:{ref}:{alt}"
            fout.write(f"{snp_id}\t{chrom}\t{pos}\t{ref}\t{alt}\t{beta}\t{se}\t{pval}\t{n_v}\t{ncase_v}\t{ncont_v}\n")
            n_written += 1
        except (ValueError, KeyError):
            continue
print(f"  Written {n_written:,} variants -> {outfile}")
print(f"PD Nalls2019 written to {outfile}")
PDPY
else
  echo "  SKIP: $PD_RAW not downloaded. Get manually from:"
  echo "  $PD_URL"
fi

# ===========================================================================
# 3. SCHIZOPHRENIA — Trubetskoy et al. 2022 (PGC SCZ wave 3, Nature)
#    PGC download: https://pgc.unc.edu/for-researchers/download-results/
#    Figshare: https://figshare.com/ndownloader/files/34517861
#    N=161,405 (67,390 cases, 94,015 controls)
#    PMID: 35396579
# ===========================================================================
echo "--- SCZ: Trubetskoy 2022 (PGC wave 3) ---"
SCZ_DIR="$OUTDIR/SCZ_Trubetskoy2022"
mkdir -p "$SCZ_DIR"
SCZ_RAW="$SCZ_DIR/raw.tsv.gz"
SCZ_OUT="$SCZ_DIR/SCZ_Trubetskoy2022_hg19.tsv"
SCZ_URL="https://figshare.com/ndownloader/files/34517861"

if [ ! -s "$SCZ_RAW" ]; then
  echo "Downloading PGC SCZ wave 3 (~200 MB, via figshare)..."
  # Figshare uses async download generation; use curl with longer timeout
  curl -sL --max-time 300 --retry 2 -o "$SCZ_RAW" "$SCZ_URL" && \
    [ -s "$SCZ_RAW" ] || { rm -f "$SCZ_RAW"; echo "  WARN: figshare download failed."; }
fi
if [ -s "$SCZ_RAW" ]; then
  python3 - <<'SCZPY'
import gzip, csv
infile  = "results/coloc/disease_gwas/SCZ_Trubetskoy2022/raw.tsv.gz"
outfile = "results/coloc/disease_gwas/SCZ_Trubetskoy2022/SCZ_Trubetskoy2022_hg19.tsv"
Ncase, Ncont = 67390, 94015
N = Ncase + Ncont
n_written = 0
with gzip.open(infile, "rt") as fin, open(outfile, "w") as fout:
    fout.write("snp\tchr\tpos\tref\talt\tbeta\tse\tp\tN\tNcase\tNcont\n")
    reader = csv.DictReader(fin, delimiter="\t")
    for row in reader:
        try:
            chrom = str(row.get("CHROM", row.get("CHR",""))).lstrip("chr")
            pos   = int(row.get("POS", row.get("BP","")))
            ref   = row.get("A1", row.get("REF","")).upper()
            alt   = row.get("A2", row.get("ALT","")).upper()
            beta  = row.get("BETA", row.get("OR",""))
            se    = row.get("SE","")
            pval  = row.get("P","")
            if not all([chrom, pos, ref, alt, beta, se, pval]):
                continue
            snp_id = f"chr{chrom}:{pos}:{ref}:{alt}"
            fout.write(f"{snp_id}\t{chrom}\t{pos}\t{ref}\t{alt}\t{beta}\t{se}\t{pval}\t{N}\t{Ncase}\t{Ncont}\n")
            n_written += 1
        except (ValueError, KeyError):
            continue
print(f"  Written {n_written:,} variants -> {outfile}")
print(f"SCZ Trubetskoy2022 written to {outfile}")
SCZPY
else
  echo "  SKIP: SCZ raw file not downloaded."
  echo "  Manual download: https://pgc.unc.edu/for-researchers/download-results/"
  echo "  Save as: $SCZ_RAW"
fi

# ===========================================================================
# 4a. BIPOLAR DISORDER — O'Connell et al. 2024 (PGC bip2024, EUR ex-23andMe)
#     Figshare: https://figshare.com/articles/dataset/bip2024/27216117
#     File: bip2024_eur_no23andMe.gz (daner format, hg19)
#     ~59,287 cases / ~781,022 controls (per-variant Nca/Nco in file)
# ===========================================================================
echo "--- BD: bip2024 EUR (ex-23andMe) ---"
BD24_DIR="$OUTDIR/BD_bip2024"
mkdir -p "$BD24_DIR"
BD24_RAW="$BD24_DIR/raw.gz"
BD24_OUT="$BD24_DIR/BD_bip2024_hg19.tsv"
BD24_URL="https://ndownloader.figshare.com/files/49760772"

if [ ! -s "$BD24_RAW" ]; then
  echo "Downloading bip2024 EUR (~332 MB, via figshare)..."
  curl -sL --max-time 600 --retry 2 -o "$BD24_RAW" "$BD24_URL" && \
    [ -s "$BD24_RAW" ] || { rm -f "$BD24_RAW"; echo "  WARN: figshare download failed."; }
fi
if [ -s "$BD24_RAW" ]; then
  python3 - <<'BD24PY'
import gzip, math, csv
infile  = "results/coloc/disease_gwas/BD_bip2024/raw.gz"
outfile = "results/coloc/disease_gwas/BD_bip2024/BD_bip2024_hg19.tsv"
n_written = 0
with gzip.open(infile, "rt") as fin, open(outfile, "w") as fout:
    fout.write("snp\tchr\tpos\tref\talt\tbeta\tse\tp\tN\tNcase\tNcont\n")
    header = fin.readline().strip().split()
    for line in fin:
        parts = line.strip().split()
        if len(parts) < 17:
            continue
        try:
            chrom = str(parts[1]).lstrip("chr")
            pos   = int(parts[2])
            ref   = parts[3].upper()
            alt   = parts[4].upper()
            or_   = float(parts[6])
            beta  = math.log(or_) if or_ > 0 else None
            se    = float(parts[7])
            pval  = float(parts[8])
            nca   = int(float(parts[-4]))
            nco   = int(float(parts[-3]))
            n_v   = nca + nco
            if beta is None or se <= 0 or pval <= 0:
                continue
            snp_id = f"chr{chrom}:{pos}:{ref}:{alt}"
            fout.write(f"{snp_id}\t{chrom}\t{pos}\t{ref}\t{alt}\t{beta}\t{se}\t{pval}\t{n_v}\t{nca}\t{nco}\n")
            n_written += 1
        except (ValueError, IndexError):
            continue
print(f"  Written {n_written:,} variants -> {outfile}")
BD24PY
  echo "BD bip2024 written to $BD24_OUT"
else
  echo "  SKIP: $BD24_RAW not downloaded."
  echo "  Manual: $BD24_URL"
fi

# ===========================================================================
# 4b. BIPOLAR DISORDER — Mullins et al. 2021 (legacy; kept for reference)
#    Figshare: https://figshare.com/ndownloader/files/26603681
#    N=413,466 (41,917 cases, 371,549 controls)
#    PMID: 34385711
# ===========================================================================
echo "--- BD: Mullins 2021 (legacy) ---"
BD_DIR="$OUTDIR/BD_Mullins2021"
mkdir -p "$BD_DIR"
BD_RAW="$BD_DIR/raw.tsv.gz"
BD_OUT="$BD_DIR/BD_Mullins2021_hg19.tsv"
BD_URL="https://figshare.com/ndownloader/files/26603681"

if [ ! -s "$BD_RAW" ]; then
  echo "Downloading Mullins 2021 BD (~140 MB, via figshare)..."
  curl -sL --max-time 300 --retry 2 -o "$BD_RAW" "$BD_URL" && \
    [ -s "$BD_RAW" ] || { rm -f "$BD_RAW"; echo "  WARN: figshare download failed."; }
fi
if [ -s "$BD_RAW" ]; then
  python3 - <<'BDPY'
import gzip, csv
infile  = "results/coloc/disease_gwas/BD_Mullins2021/raw.tsv.gz"
outfile = "results/coloc/disease_gwas/BD_Mullins2021/BD_Mullins2021_hg19.tsv"
Ncase, Ncont = 41917, 371549
N = Ncase + Ncont
n_written = 0
with gzip.open(infile, "rt") as fin, open(outfile, "w") as fout:
    fout.write("snp\tchr\tpos\tref\talt\tbeta\tse\tp\tN\tNcase\tNcont\n")
    reader = csv.DictReader(fin, delimiter="\t")
    for row in reader:
        try:
            chrom = str(row.get("CHR","")).lstrip("chr")
            pos   = int(row.get("BP",""))
            ref   = row.get("A1","").upper()
            alt   = row.get("A2","").upper()
            beta  = row.get("BETA", row.get("LOG_ODDS",""))
            se    = row.get("SE","")
            pval  = row.get("P","")
            if not all([chrom, pos, ref, alt, beta, se, pval]):
                continue
            snp_id = f"chr{chrom}:{pos}:{ref}:{alt}"
            fout.write(f"{snp_id}\t{chrom}\t{pos}\t{ref}\t{alt}\t{beta}\t{se}\t{pval}\t{N}\t{Ncase}\t{Ncont}\n")
            n_written += 1
        except (ValueError, KeyError):
            continue
print(f"  Written {n_written:,} variants -> {outfile}")
print(f"BD Mullins2021 written to {outfile}")
BDPY
else
  echo "  SKIP: BD raw file not downloaded."
  echo "  Manual download: https://pgc.unc.edu/for-researchers/download-results/"
  echo "  Save as: $BD_RAW"
fi

# ===========================================================================
# 5a. MAJOR DEPRESSIVE DISORDER — PGC MDD2025 EUR (ex-23andMe, GRCh37)
#     Figshare: https://figshare.com/articles/dataset/GWAS_summary_statistics_for_major_depression_PGC_MDD2025_/27061255
#     File: pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz
#     N=412,305 cases / 1,588,397 controls
# ===========================================================================
echo "--- MDD: PGC MDD2025 EUR (ex-23andMe) ---"
MDD25_DIR="$OUTDIR/MDD_MDD2025"
mkdir -p "$MDD25_DIR"
MDD25_RAW="$MDD25_DIR/raw.gz"
MDD25_OUT="$MDD25_DIR/MDD_MDD2025_hg19.tsv"
MDD25_URL="https://ndownloader.figshare.com/files/51487019"

if [ ! -s "$MDD25_RAW" ]; then
  echo "Downloading MDD2025 EUR (~223 MB, via figshare)..."
  curl -sL --max-time 600 --retry 2 -o "$MDD25_RAW" "$MDD25_URL" && \
    [ -s "$MDD25_RAW" ] || { rm -f "$MDD25_RAW"; echo "  WARN: figshare download failed."; }
fi
if [ -s "$MDD25_RAW" ]; then
  python3 - <<'MDD25PY'
import gzip, csv
infile  = "results/coloc/disease_gwas/MDD_MDD2025/raw.gz"
outfile = "results/coloc/disease_gwas/MDD_MDD2025/MDD_MDD2025_hg19.tsv"
Ncase, Ncont = 412305, 1588397
N = Ncase + Ncont
n_written = 0
with gzip.open(infile, "rt") as fin, open(outfile, "w") as fout:
    fout.write("snp\tchr\tpos\tref\talt\tbeta\tse\tp\tN\tNcase\tNcont\n")
    for line in fin:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 14:
            continue
        try:
            chrom = str(parts[0]).lstrip("chr")
            pos   = int(parts[1])
            ref   = parts[3].upper()  # EA
            alt   = parts[4].upper()  # NEA
            beta  = float(parts[5])
            se    = float(parts[6])
            pval  = float(parts[7])
            nca   = int(float(parts[12])) if parts[12] else Ncase
            nco   = int(float(parts[13])) if parts[13] else Ncont
            n_v   = nca + nco
            if se <= 0 or pval <= 0:
                continue
            snp_id = f"chr{chrom}:{pos}:{ref}:{alt}"
            fout.write(f"{snp_id}\t{chrom}\t{pos}\t{ref}\t{alt}\t{beta}\t{se}\t{pval}\t{n_v}\t{nca}\t{nco}\n")
            n_written += 1
        except (ValueError, IndexError):
            continue
print(f"  Written {n_written:,} variants -> {outfile}")
MDD25PY
  echo "MDD MDD2025 written to $MDD25_OUT"
else
  echo "  SKIP: $MDD25_RAW not downloaded."
  echo "  Manual: $MDD25_URL"
fi

# ===========================================================================
# 5b. MAJOR DEPRESSIVE DISORDER — Wray et al. 2018 (legacy)
#    GWAS Catalog: GCST005839  (PGC MDD 2018, excluding 23andMe)
#    N~173,005 (59,851 cases, 113,154 controls)
#    PMID: 29700475
# ===========================================================================
echo "--- MDD: Wray 2018 (legacy, ex-23andMe) ---"
MDD_DIR="$OUTDIR/MDD_Wray2018"
mkdir -p "$MDD_DIR"
MDD_RAW="$MDD_DIR/raw.gz"
MDD_OUT="$MDD_DIR/MDD_Wray2018_hg19.tsv"
MDD_URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005839/MDD2018_ex23andMe.gz"

if [ ! -s "$MDD_RAW" ]; then
  echo "Downloading Wray 2018 MDD (~469 MB)..."
  dl "$MDD_URL" "$MDD_RAW" || true
fi
if [ -s "$MDD_RAW" ]; then
  python3 - <<'MDDPY'
import gzip, csv, math
infile  = "results/coloc/disease_gwas/MDD_Wray2018/raw.gz"
outfile = "results/coloc/disease_gwas/MDD_Wray2018/MDD_Wray2018_hg19.tsv"
# Columns: CHR SNP BP A1 A2 FRQ_A_59851 FRQ_U_113154 INFO OR SE P ngt Direction HetISqt HetDf HetPVa Nca Nco Neff
n_written = 0
with gzip.open(infile, "rt") as fin, open(outfile, "w") as fout:
    fout.write("snp\tchr\tpos\tref\talt\tbeta\tse\tp\tN\tNcase\tNcont\n")
    reader = csv.DictReader(fin, delimiter="\t")
    for row in reader:
        try:
            chrom = str(row.get("CHR","")).lstrip("chr")
            pos   = int(row.get("BP",""))
            ref   = row.get("A1","").upper()
            alt   = row.get("A2","").upper()
            or_   = float(row.get("OR",""))
            beta  = str(math.log(or_)) if or_ > 0 else "NA"
            se    = row.get("SE","")
            pval  = row.get("P","")
            # Use per-variant Nca/Nco if present
            nca   = row.get("Nca", "59851")
            nco   = row.get("Nco", "113154")
            try:
                n_v = int(float(nca)) + int(float(nco))
            except ValueError:
                n_v = 172953
            if not all([chrom, pos, ref, alt, beta, se, pval]):
                continue
            snp_id = f"chr{chrom}:{pos}:{ref}:{alt}"
            fout.write(f"{snp_id}\t{chrom}\t{pos}\t{ref}\t{alt}\t{beta}\t{se}\t{pval}\t{n_v}\t{nca}\t{nco}\n")
            n_written += 1
        except (ValueError, KeyError):
            continue
print(f"  Written {n_written:,} variants -> {outfile}")
print(f"MDD Wray2018 written to {outfile}")
MDDPY
else
  echo "  SKIP: $MDD_RAW not downloaded. Get manually from:"
  echo "  $MDD_URL"
fi

# ===========================================================================
# 6. FTD/FTLD — van der Lee et al. 2019 meta-analysis (hg19)
#    *** MANUAL DOWNLOAD REQUIRED ***
#    Request from: https://www.niagads.org/ or contact FTD consortium
# ===========================================================================
echo "--- FTD: Manual download required ---"
FTD_DIR="$OUTDIR/FTD_vanderLee2019"
mkdir -p "$FTD_DIR"
cat > "$FTD_DIR/README.txt" <<'FTDREADME'
FTD/FTLD GWAS - Manual Download Required
==========================================
Recommended: van der Lee et al. 2019 (GCST009258) or newer FTD consortia data.
Request access from: https://www.niagads.org/ or contact the FTD consortium.

Once downloaded, standardize to:
  snp (chr:pos:REF:ALT), chr, pos, ref, alt, beta, se, p, N, Ncase, Ncont
and save as: FTD_vanderLee2019_hg19.tsv in this directory.
FTDREADME
echo "  See $FTD_DIR/README.txt for instructions."

# ===========================================================================
# 7. LEWY BODY DEMENTIA — Chia et al. 2021 (Nature Genetics)
#    GWAS Catalog: GCST90001390  (harmonised)
#    N=8,392 (2,591 cases, 5,801 controls)
#    PMID: 33589841
# ===========================================================================
echo "--- LBD: Chia 2021 ---"
LBD_DIR="$OUTDIR/LBD_Chia2021"
mkdir -p "$LBD_DIR"
LBD_RAW="$LBD_DIR/raw.tsv.gz"
LBD_OUT="$LBD_DIR/LBD_Chia2021_hg19.tsv"
LBD_URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90001001-GCST90002000/GCST90001390/harmonised/33589841-GCST90001390-EFO_0006792.h.tsv.gz"

if [ ! -s "$LBD_RAW" ]; then
  echo "Downloading Chia 2021 LBD (~391 MB)..."
  dl "$LBD_URL" "$LBD_RAW" || true
fi
if [ -s "$LBD_RAW" ]; then
  standardize_harmonised "$LBD_RAW" "$LBD_OUT" 8392 2591 5801
  echo "LBD Chia2021 written to $LBD_OUT"
else
  echo "  SKIP: $LBD_RAW not downloaded. Get manually from:"
  echo "  $LBD_URL"
fi

# ===========================================================================
# 8. ALS — van Rheenen et al. 2021 (Nature Genetics)
#    GWAS Catalog: GCST90012201  (harmonised, PMID 32915819)
#    N=138,086 (27,205 cases, 110,881 controls)
# ===========================================================================
echo "--- ALS: van Rheenen 2021 ---"
ALS_DIR="$OUTDIR/ALS_vanRheenen2021"
mkdir -p "$ALS_DIR"
ALS_RAW="$ALS_DIR/raw.tsv.gz"
ALS_OUT="$ALS_DIR/ALS_vanRheenen2021_hg19.tsv"
ALS_URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90012001-GCST90013000/GCST90012201/harmonised/32915819-GCST90012201-EFO_0005680.h.tsv.gz"

if [ ! -s "$ALS_RAW" ]; then
  echo "Downloading van Rheenen 2021 ALS (~12 MB)..."
  dl "$ALS_URL" "$ALS_RAW" || true
fi
if [ -s "$ALS_RAW" ]; then
  standardize_harmonised "$ALS_RAW" "$ALS_OUT" 138086 27205 110881
  echo "ALS vanRheenen2021 written to $ALS_OUT"
else
  echo "  SKIP: $ALS_RAW not downloaded. Get manually from:"
  echo "  $ALS_URL"
fi

echo ""
echo "=== Download complete. Check $OUTDIR for standardized files. ==="
echo "NOTE: FTD requires manual download — see $OUTDIR/FTD_vanderLee2019/README.txt"
echo "NOTE: SCZ/BD/MDD use figshare — if download failed, get from:"
echo "  SCZ:  https://pgc.unc.edu/for-researchers/download-results/ (PGC3 SCZ wave3)"
echo "  BD:   https://figshare.com/articles/dataset/bip2024/27216117 (bip2024_eur_no23andMe.gz)"
echo "  MDD:  https://figshare.com/articles/dataset/GWAS_summary_statistics_for_major_depression_PGC_MDD2025_/27061255"
echo "Active coloc diseases: use BD_bip2024 and MDD_MDD2025 (legacy Mullins/Wray kept in archive/)"
