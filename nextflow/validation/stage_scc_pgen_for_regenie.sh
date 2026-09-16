#!/bin/bash
# Login-node only (/project is writable here). Points REGENIE step 2's
# ${projectDir}/results/<cohort>/*.QC.final at readable pgen files.
# Does not write under Xiaolin/SCC/.

set -euo pipefail
NF_DIR="${NF_DIR:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow}"
SCC_RES="${SCC_RES:-/project/rrg-shreejoy/zhoux156/Xiaolin/SCC/nextflow/results}"
LOCAL_OUT="${LOCAL_OUT:-/scratch/zhoux156/results/validation/gwas_Astrocyte}"

cfg_prefix() {
    case "$1" in
        NIMH_HBCC_*) echo "CMC_HBCC.QC.final" ;;
        *)           echo "${1}.QC.final" ;;
    esac
}

stage_one() {
    local cohort="$1"
    local cfgp src dest
    cfgp=$(cfg_prefix "${cohort}")
    src=""
    if [ -s "${SCC_RES}/${cohort}/${cfgp}.pgen" ]; then
        src="${SCC_RES}/${cohort}/${cfgp}"
    else
        local f
        for f in "${SCC_RES}/${cohort}"/*.QC.final.pgen; do
            if [ -s "${f}" ]; then
                src="${f%.pgen}"
                break
            fi
        done
    fi
    if [ -z "${src}" ] || [ ! -s "${src}.pgen" ]; then
        if [ -s "${LOCAL_OUT}/${cohort}/${cfgp}.pgen" ]; then
            src="${LOCAL_OUT}/${cohort}/${cfgp}"
            echo "${cohort}: LOCAL ${src}.pgen (SCC pgen missing/dangling)"
        else
            echo "${cohort}: MISSING" >&2
            return 1
        fi
    else
        echo "${cohort}: SCC  ${src}.pgen"
    fi
    dest="${NF_DIR}/results/${cohort}"
    mkdir -p "${dest}"
    ln -sfn "${src}.pgen" "${dest}/${cfgp}.pgen"
    ln -sfn "${src}.psam" "${dest}/${cfgp}.psam"
    ln -sfn "${src}.pvar" "${dest}/${cfgp}.pvar"
    # sanity
    [ -s "${dest}/${cfgp}.pgen" ] && [ -s "${dest}/${cfgp}.psam" ] && [ -s "${dest}/${cfgp}.pvar" ]
}

fail=0
for cohort in ROSMAP ROSMAP_array Mayo MSBB CMC_MSSM CMC_PENN CMC_PITT \
              GTEx_v10 NABEC NIMH_HBCC_1M NIMH_HBCC_h650 NIMH_HBCC_Omni5M \
              GVEX AMP_AD_Rush AMP_AD_Mayo; do
    stage_one "${cohort}" || fail=1
done
exit "${fail}"
