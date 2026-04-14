/*
 * RUN_LDSC_H2
 *
 * Runs munge_sumstats.py followed by ldsc.py --h2 for one cell type.
 */

process RUN_LDSC_H2 {
    label 'medium_memory'
    tag "${cell_type}_ldsc"

    publishDir "${sumstats_output_dir}", mode: 'copy', overwrite: true, pattern: "*.sumstats.gz"
    publishDir "${sumstats_output_dir}", mode: 'copy', overwrite: true, pattern: "*.sumstats.gz.log"
    publishDir "${sumstats_output_dir}", mode: 'copy', overwrite: true, pattern: "*.progress.log"
    publishDir "${results_output_dir}", mode: 'copy', overwrite: true, pattern: "${cell_type}.h2*"

    input:
    tuple val(cell_type),
          path(ldsc_input_tsv),
          val(ldsc_conda_env),
          val(merge_alleles),
          val(ref_ld_chr),
          val(w_ld_chr),
          val(sumstats_output_dir),
          val(results_output_dir)

    output:
    path "${cell_type}.sumstats.gz", emit: munged_sumstats
    path "${cell_type}.sumstats.gz.log", emit: munge_log
    path "${cell_type}.progress.log", emit: progress_log
    path "${cell_type}.h2.log", emit: h2_log
    path "${cell_type}.h2.results", emit: h2_results

    script:
    def conda_init = '''
        export xml_catalog_files_libxml2="${xml_catalog_files_libxml2:-}"
        set +u
        if [ -f "$HOME/.anaconda3/etc/profile.d/conda.sh" ]; then
            source "$HOME/.anaconda3/etc/profile.d/conda.sh" 2>/dev/null || true
        elif [ -f "/nethome/kcni/xzhou/.anaconda3/etc/profile.d/conda.sh" ]; then
            source "/nethome/kcni/xzhou/.anaconda3/etc/profile.d/conda.sh" 2>/dev/null || true
        elif [ -f "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" ]; then
            source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
        fi
        set -u
    '''
    def conda_activate = ldsc_conda_env ? """
        set +u
        if ! conda activate ${ldsc_conda_env} 2>/dev/null; then
            echo "ERROR: Failed to activate conda environment '${ldsc_conda_env}'" >&2
            exit 1
        fi
        set -u
    """ : ""
    def merge_alleles_arg = merge_alleles ? "--merge-alleles \"${merge_alleles}\"" : ""
    def ref_ld_chr_arg = ref_ld_chr?.endsWith('/') ? ref_ld_chr : "${ref_ld_chr}/"
    def w_ld_chr_arg = w_ld_chr?.endsWith('/') ? w_ld_chr : "${w_ld_chr}/"
    """
    ${conda_init}
    ${conda_activate}

    if ! command -v munge_sumstats.py >/dev/null 2>&1; then
        echo "ERROR: munge_sumstats.py not found on PATH" >&2
        exit 1
    fi
    if ! command -v ldsc.py >/dev/null 2>&1; then
        echo "ERROR: ldsc.py not found on PATH" >&2
        exit 1
    fi

    set -o pipefail

    progress_log="${cell_type}.progress.log"
    : > "\${progress_log}"
    log_step() {
        printf '%s\t%s\n' "\$(date '+%Y-%m-%d %H:%M:%S')" "\$1" | tee -a "\${progress_log}"
    }

    if [ -s "${ldsc_input_tsv}" ]; then
        input_rows=\$(awk 'END{print NR-1}' "${ldsc_input_tsv}")
    else
        input_rows=0
    fi
    log_step "START cell_type=${cell_type} input=${ldsc_input_tsv} input_rows=\${input_rows}"
    log_step "MUNGE_SUMSTATS_START"

    munge_sumstats.py \\
        --sumstats "${ldsc_input_tsv}" \\
        --out "${cell_type}.sumstats" \\
        --snp SNP \\
        --a1 A1 \\
        --a2 A2 \\
        --p P \\
        --N-col N \\
        --signed-sumstats BETA,0 \\
        --chunksize 500000 \\
        ${merge_alleles_arg} 2>&1 | tee "${cell_type}.munge.stdout.log"

    if [ ! -s "${cell_type}.sumstats.sumstats.gz" ]; then
        log_step "MUNGE_SUMSTATS_FAILED missing=${cell_type}.sumstats.sumstats.gz"
        exit 1
    fi
    munged_bytes=\$(wc -c < "${cell_type}.sumstats.sumstats.gz")
    log_step "MUNGE_SUMSTATS_DONE munged_bytes=\${munged_bytes}"
    log_step "LDSC_H2_START ref_ld_chr=${ref_ld_chr} w_ld_chr=${w_ld_chr}"

    ldsc.py \\
        --h2 "${cell_type}.sumstats.sumstats.gz" \\
        --ref-ld-chr "${ref_ld_chr_arg}" \\
        --w-ld-chr "${w_ld_chr_arg}" \\
        --out "${cell_type}.h2" 2>&1 | tee "${cell_type}.ldsc.stdout.log"

    # ldsc.py --h2 writes ${cell_type}.h2.log, not ${cell_type}.h2.results.
    # Keep a .h2.results artifact for downstream compatibility and Nextflow outputs.
    cp "${cell_type}.h2.log" "${cell_type}.h2.results"

    mv "${cell_type}.sumstats.sumstats.gz" "${cell_type}.sumstats.gz"
    mv "${cell_type}.sumstats.log" "${cell_type}.sumstats.gz.log"

    h2_line=\$(awk '/Total Observed scale h2/ {line=\$0} END{print line}' "${cell_type}.h2.log")
    intercept_line=\$(awk '/Intercept:/ {line=\$0} END{print line}' "${cell_type}.h2.log")
    if [ -n "\${h2_line}" ]; then
        log_step "LDSC_H2_SUMMARY \${h2_line}"
    fi
    if [ -n "\${intercept_line}" ]; then
        log_step "LDSC_H2_SUMMARY \${intercept_line}"
    fi
    log_step "DONE results=${cell_type}.h2.results"
    """
}
