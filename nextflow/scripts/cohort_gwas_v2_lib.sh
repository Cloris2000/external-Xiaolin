#!/bin/bash
# Shared helpers for the cohort_gwas_v2 re-run (preflight + driver + meta).
# Source after site_env.sh.
#
# Layout (all on scratch; nothing is written under /project):
#   ${RESULTS_ROOT}/cohort_gwas_v2/results/<cohort>/   per-cohort pipeline output
#   ${RESULTS_ROOT}/cohort_gwas_v2/compare/<cohort>/   comparison vs the SCC run
#   ${RESULTS_ROOT}/cohort_gwas_v2/preflight/          input manifest
# The extra "results/" level exists because META_COHORT_AUDIT globs
# {base_dir}/results/{cohort}/regenie_step2, so base_dir = .../cohort_gwas_v2.

V2_ROOT="${V2_ROOT:-${RESULTS_ROOT}/cohort_gwas_v2}"
V2_RESULTS="${V2_ROOT}/results"
V2_COMPARE="${V2_ROOT}/compare"
V2_PREFLIGHT="${V2_ROOT}/preflight"
V2_WORK="${WORK_ROOT}/cohort_gwas_v2"
V2_LOGS="${LOG_ROOT}/cohort_gwas_v2"
V2_MAP="${NF_DIR}/validation/cohort_trillium_paths.tsv"
SCC_RES="${SCC_DIR}/results"

ALL_COHORTS=(ROSMAP ROSMAP_array Mayo MSBB CMC_MSSM CMC_PENN CMC_PITT GTEx_v10 NABEC
             NIMH_HBCC_1M NIMH_HBCC_h650 NIMH_HBCC_Omni5M GVEX AMP_AD_Rush AMP_AD_Mayo)

# Read a single-quoted or double-quoted scalar param from a Nextflow config.
cfg_param() {  # cfg_param <config> <param>
    grep -oP "^\s*${2}\s*=\s*['\"]\K[^'\"]*" "$1" | head -1
}

# QC study name (drives <study>.QC.<chr>.normalized.vcf.gz under normalized_vcf_dir
# and the <study>.QC.final pgen prefix).  HBCC configs set genotyping_study.
qc_study_name() {  # qc_study_name <config>
    local g; g=$(cfg_param "$1" genotyping_study)
    [ -n "$g" ] && echo "$g" || cfg_param "$1" study
}

# GVEX: nextflow.config.combined.gvex reads the raw, uncompressed dosage VCFs
# (BrainGVEX_chr{N}.vcf, ~196 GB).  This run uses the bcftools `norm -m-` copies
# of the same files instead (same 400 samples and GT/DS/GP fields, multiallelics
# split, bgzipped, 9 GB) - the form every other cohort's input already has.
# Approved 2026-09-12.  Applied as an extra -c config (nextflow.config.gvex_normalized_input)
# loaded after the cohort config, so the cohort config is not edited.  Not as
# `--vcf_pattern ''`: Nextflow turns an empty CLI value into Boolean true.
GVEX_NORMALIZED_DIR="/project/rrg-shreejoy/GVEX/Genotype/BrainGVEX_combined_vcf_normalized"

# Append per-cohort genotype-input overrides (extra -c files) to the named array.
# Must be appended AFTER the cohort config's own -c so the override wins.
genotype_override_args() {  # genotype_override_args <cohort> <array-name>
    local -n _arr="$2"
    case "$1" in
        GVEX) _arr+=(-c "${NF_DIR}/nextflow.config.gvex_normalized_input") ;;
    esac
}

# Absolute path of the per-chromosome genotype input for chr $2, honouring the
# overrides above when the cohort is given.
resolve_vcf() {  # resolve_vcf <config> <chrom> [cohort]
    local cfg="$1" chrom="$2" cohort="${3:-}" pat nd study
    study=$(qc_study_name "$cfg")
    if [ "$cohort" = "GVEX" ]; then
        echo "${GVEX_NORMALIZED_DIR}/${study}.QC.${chrom}.normalized.vcf.gz"; return
    fi
    pat=$(cfg_param "$cfg" vcf_pattern)
    if [ -n "$pat" ]; then
        echo "${pat//\{chrom\}/$chrom}"; return
    fi
    nd=$(cfg_param "$cfg" normalized_vcf_dir)
    echo "${nd}/${study}.QC.${chrom}.normalized.vcf.gz"
}

# samples_to_keep: the configs reference ${projectDir}/results/<cohort>/samples_to_keep.txt,
# which this checkout does not have.  The SCC tree carries the files that were
# actually used (ROSMAP_array: the 170 array-only donors not in WGS; Mayo; MSBB).
#
# pheno_prep.R reads the file with read.table(sep="\t", header=TRUE) and takes the
# IID column.  The SCC Mayo file has a SPACE-separated header ("FID IID") over
# tab-separated rows, which R parses as one column named FID.IID -> no IID column
# -> zero samples kept.  So the file is copied to scratch with the header's
# whitespace normalised to tabs (rows untouched) and that copy is what gets passed.
# Prints the path to pass on the CLI, or nothing if the cohort has none.
resolve_samples_to_keep() {  # resolve_samples_to_keep <cohort> <config>
    local cohort="$1" cfg="$2" raw src
    raw=$(cfg_param "$cfg" samples_to_keep)
    [ -z "$raw" ] && return 0
    local local_path="${raw//\$\{projectDir\}/$NF_DIR}"
    if   [ -s "$local_path" ];                              then src="$local_path"
    elif [ -s "${SCC_RES}/${cohort}/samples_to_keep.txt" ]; then src="${SCC_RES}/${cohort}/samples_to_keep.txt"
    else return 0; fi   # missing -> pheno_prep.R warns and keeps all phenotyped samples
    local dst_dir="${V2_PREFLIGHT}/samples_to_keep" dst
    dst="${dst_dir}/${cohort}.tsv"
    mkdir -p "$dst_dir"
    awk 'NR==1 { gsub(/[ \t]+/, "\t") } { print }' "$src" > "$dst"
    echo "$dst"
}

# Cohort map row -> variables (empty fields become "-").
# Usage: while read_map_row; do ...; done < <(map_rows)
map_rows() {
    sed -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' -e 's/\t$/\t-/' "${V2_MAP}" | tail -n +2
}

# BGZF EOF marker check (truncated bgzip files lack it).
bgzf_eof_ok() {  # bgzf_eof_ok <file.gz>
    python3 - "$1" <<'EOF'
import sys, os
p = sys.argv[1]; eof = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
sz = os.path.getsize(p)
with open(p, "rb") as fh:
    fh.seek(max(0, sz - 28)); sys.exit(0 if fh.read(28) == eof else 1)
EOF
}
