"""
Normalize AMP-AD Diverse Cohort WGS VCFs
=========================================
Splits multi-allelic records into biallelic using bcftools norm -m-any.
This mirrors the NABEC preprocessing step and is required before the
combined_pipeline_v2.nf genotyping QC, which expects biallelic input.

Run from the nextflow project root:
    python scripts/normalize_amp_ad_diverse_vcfs.py

Outputs to:
    /external/rprshnas01/netdata_kcni/stlab/AMP_AD_Diverse/Genotype/normalized/
    AMP_AD_Diverse.chr{N}.normalized.vcf.gz  (+  .tbi index)

After all jobs complete, update the three AMP-AD Nextflow configs:
    normalized_vcf_dir = ".../Genotype/normalized"
    vcf_pattern        = ".../Genotype/normalized/AMP_AD_Diverse.chr{chrom}.normalized.vcf.gz"
"""

import os
import subprocess

# ── paths ─────────────────────────────────────────────────────────────────────
BCFTOOLS       = "/nethome/kcni/xzhou/.anaconda3/envs/bcftools_env/bin/bcftools"
INPUT_VCF_DIR  = "/external/rprshnas01/netdata_kcni/stlab/AMP_AD_Diverse/Genotype"
INPUT_PATTERN  = "DivCo_GRCh38_743Samples_JointCalls_Recalibrated_Annotated_05-25-2025.chr{chrom}.vcf.gz"
OUTPUT_DIR     = os.path.join(INPUT_VCF_DIR, "normalized")
OUTPUT_PATTERN = "AMP_AD_Diverse.chr{chrom}.normalized.vcf.gz"
LOG_DIR        = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/normalize_amp_ad"

# ── SLURM settings ─────────────────────────────────────────────────────────────
PARTITION      = "long"       # 1-day limit; chr1/2 at 25GB can take 6-10h
MEM_PER_CPU    = "32000M"     # 32 GB – bcftools norm is single-threaded, I/O heavy
TIME_LIMIT     = "0-10:0"     # 10 hours; enough for any chromosome

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(LOG_DIR,    exist_ok=True)

# ── helper: run a SLURM job ────────────────────────────────────────────────────
def run_slurm(cmd, job_name, time=TIME_LIMIT, mem=MEM_PER_CPU,
              partition=PARTITION, log_file=None):
    log = log_file or f"{LOG_DIR}/{job_name}.log"
    slurm_cmd = (
        f"sbatch --job-name={job_name} "
        f"--partition={partition} "
        f"--mem-per-cpu={mem} "
        f"--time={time} "
        f"--output={log} "
        f"--error={log} "
        f"--wrap='{cmd}'"
    )
    print(f"Submitting: {job_name}")
    result = subprocess.run(slurm_cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr.strip()}")
    else:
        print(f"  {result.stdout.strip()}")


# ── submit one job per autosome ────────────────────────────────────────────────
submitted = 0
skipped   = 0

for chrom in range(1, 23):
    input_vcf  = os.path.join(INPUT_VCF_DIR, INPUT_PATTERN.format(chrom=chrom))
    output_vcf = os.path.join(OUTPUT_DIR,    OUTPUT_PATTERN.format(chrom=chrom))
    output_tbi = output_vcf + ".tbi"

    if not os.path.exists(input_vcf):
        print(f"chr{chrom}: input VCF not found, skipping ({input_vcf})")
        continue

    if os.path.exists(output_vcf) and os.path.exists(output_tbi):
        print(f"chr{chrom}: already done, skipping")
        skipped += 1
        continue

    # bcftools norm -m-any  → split multi-allelic into biallelic records
    # (no -f ref.fasta needed for split-only; left-alignment skipped intentionally
    #  to keep consistent with NABEC preprocessing approach)
    norm_cmd = (
        f"{BCFTOOLS} norm -m-any "
        f"{input_vcf} "
        f"-Oz -o {output_vcf} "
        f"&& {BCFTOOLS} index -t {output_vcf}"
    )

    # Larger chromosomes get a bit more time
    time_limit = "0-12:0" if chrom <= 5 else TIME_LIMIT

    run_slurm(
        norm_cmd,
        job_name=f"amp_ad_norm_chr{chrom}",
        time=time_limit,
        log_file=f"{LOG_DIR}/normalize_chr{chrom}.log",
    )
    submitted += 1

print(f"\nDone. {submitted} jobs submitted, {skipped} chromosomes already complete.")
print(f"\nMonitor with:  squeue -u {os.environ.get('USER','xzhou')} | grep amp_ad_norm")
print(f"Logs in:       {LOG_DIR}/")
print(f"\nOnce all jobs finish, update the three AMP-AD configs:")
print(f"  normalized_vcf_dir = \"{OUTPUT_DIR}\"")
print(f"  vcf_pattern        = \"{OUTPUT_DIR}/{OUTPUT_PATTERN.format(chrom='{chrom}')}\"")
print(f"  (and remove vcf_pattern_raw / keep vcf_pattern pointing at normalized files)")
