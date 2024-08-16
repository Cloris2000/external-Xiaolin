# This script is inherited from QC_ROSMAP_WGS.py, and we want to add sex check in this script

import pandas as pd
import os, subprocess, pandas as pd

#run
def run(cmd, pipefail=True, **kwargs):
    # Often advisable to manually run the command with stdbuf -i0 -o0 -e0
    return subprocess.run(f'set -eu{"o pipefail" if pipefail else ""}; {cmd}',
                          check=True, shell=True, executable='/bin/bash',
                          **kwargs)

#run_slurm
def run_slurm(cmd, job_name, time, log_file='/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/QC/ROSMAP_new_log/rosmap_joint_WGS.log', num_threads=6,
              mem_per_cpu='16000M', num_nodes=1):
    # 16000 MiB = 15.625 GiB, the memory per CPU on the 12-CPU SCC nodes:
    # $ scontrol show node node03 | grep CfgTRES
    # CfgTRES=cpu=12,mem=187.50G,billing=58
    assert ' ' not in job_name
    from tempfile import NamedTemporaryFile
    try:
        with NamedTemporaryFile('w', dir='', suffix='.sh',
                                delete=False) as temp_file:
            print(f'#!/bin/bash\n'
                  f'#SBATCH --job-name={job_name}\n'
                  f'#SBATCH --time={time}\n'
                  f'#SBATCH --cpus-per-task={num_threads}\n'
                  f'#SBATCH --mem-per-cpu={mem_per_cpu}\n'
                  f'#SBATCH --nodes={num_nodes}\n'
                  f'#SBATCH --output={log_file}\n'
                  f'export MKL_NUM_THREADS={num_threads}\n'
                  f'set -euo pipefail; {cmd}\n',
                  file=temp_file)
        sbatch_message = run(f'sbatch {temp_file.name}',
                             capture_output=True).stdout.decode().rstrip('\n')
        print(f'{sbatch_message} ("{job_name}")')
    finally:
        try:
            os.unlink(temp_file.name)
        except NameError:
            pass  #define run_slurm module



# Set to path where outputs should be saved
os.chdir('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/QC/ROSMAP_joint_WGS')
# Michael has the most updated plink
plink = '/external/rprshnas01/kcni/mwainberg/software/plink2'
plink_old = '/external/rprshnas01/kcni/mwainberg/software/plink'
bcftools = '/external/rprshnas01/netdata_kcni/stlab/Xiaolin/software/bcftools-1.12/bcftools'
# vdf file path
# mv *.normalized.vcf.gz /external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_WGS_vcf_normalized
normalized_vcf_path = '/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_WGS_vcf_normalized'
vcf_path = '"/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_joint_WGS"'
log_dir = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/QC/ROSMAP_new_log"
# path to rs id mapping files TODO: Make sure the rs id mapping file is correct for AMP-AD
#rs_path = '/external/rprshnas01/external_data/rosmap/genotype/TOPmed_imputed/vcf/merged'
# Study for output file naming
study = "ROSMAP"
# Want to filter for relatedness? Default is False as regenie can deal with relatedness
relatedness = True
# Indicate if sex check is needed
sex = True


# Error: Writing ROSMAP.QC.1.pgen ... 0%/var/spool/slurm/d/job2296689/slurm_script: line 8: 27273 Segmentation fault      (core dumped)
# bcftools view -H DEJ_11898_B01_GRM_WGS_2017-05-15_1.recalibrated_variants.vcf.gz | head
################################################################################
# Normalize the VCF File
# To split multiallelic variants into biallelic records, we can use bcftools norm:
# Do this step per chromosome
# If X chromosome is available use: for chrom in list(range(1, 23)) + ['X']]:
# Can be a bit slow, might consider using run_slurm for this step instead, esp on larger cohorts
################################################################################
# for chrom in range(1, 23):
#     if not os.listdir(normalized_vcf_path):
#         run_slurm(f'{bcftools} '
#             f'norm -m-any '
#             f'{vcf_path}/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_{chrom}.recalibrated_variants.Broad_Rush.vcf.gz '
#             f'-Oz -o {normalized_vcf_path}/{study}.QC.{chrom}.normalized.vcf.gz', job_name="ROSMAP_WGS_remove_multiallelic", time="0-4:0", log_file=f'{log_dir}/ROSMAP_WGS_remove_multiallelic.{chrom}.log')

################################################################################
# Create pgen file
# Do this step per chromosome
# If X chromosome is available use: for chrom in list(range(1, 23)) + ['X']]:
# Can be a bit slow, might consider using run_slurm for this step instead, esp on larger cohorts
################################################################################
study = "ROSMAP_new"
for chrom in range(1, 23):
    if not os.path.exists(f'{study}.QC.{chrom}.pgen'):
        run_slurm(f'{plink} ' 
            f'--vcf {normalized_vcf_path}/ROSMAP.QC.{chrom}.normalized.vcf.gz '
            # f'--vcf-min-gq 20 '
            # f'--vcf-min-dp 10 '
            # f'--vcf-min-qual 30 '
            f'--double-id '
            f'--make-pgen '
            f'--new-id-max-allele-len 2000 '
            f'--set-all-var-ids chr@:#:\$r:\$a '  # Set variant IDs
            f'--out {study}.QC.{chrom}', job_name="ROSMAP_WGS_QC", time="0-4:0", log_file=f'{log_dir}/{study}.QC.{chrom}.pgenlog') 

        
################################################################################
# Merge pgen files TODO: Check if we want to include X chromosome?
# If X chromosome is available: run(f'ls {study}.QC.{{{{1..22}},X}}.pgen | sed -e "s/.pgen//" > {study}.merge-list.txt')
################################################################################
if not os.path.exists(f'{study}.merge-list.txt'):
    run(f'ls {study}.QC.{{1..22}}.pgen | sed -e "s/.pgen//" > {study}.merge-list.txt')

if not os.path.exists(f'{study}.QC.merged.pgen'):
    run(f'{plink} ' 
        f'--pmerge-list {study}.merge-list.txt ' 
        f'--out {study}.QC.merged')
    
# ###  check pgen files format
# for chrom in range(1, 23):
#     if not os.path.exists(f'{study}.QC.validate_new'):
#         run(f'{plink} ' 
#             f'--pfile {study}.QC.{chrom} '  # Use --pfile to specify the .pgen, .pvar, .psam files without extension
#             f'--validate '  # Use --validate without any arguments
#             f'--out {study}.QC.validate_new')      


################################################################################
# Standard GWAS QC
# Threshold considerations:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/table/mpr1608-tbl-0001/
# Note for larger cohorts hwe should be set less strict (1e-9 or 15) UKB sometimes -50
# Too stringent hwe might remove interesting SNPs - we used 1*10^(-6) here as used in other related studies
# mach is an imputation filter the most common one
# minimach is the less common one
################################################################################
study = "ROSMAP_sex" #minor allele freq from 0.01 to 0.05 + include sex check    

if not os.path.exists(f'ROSMAP_new.QC.fam') or \
        not os.path.exists(f'{study}.QC.snplist'):
    run(f'{plink} ' 
        f'--pfile ROSMAP_new.QC.merged '
        #f'--mach-r2-filter 0.8 '  #imputation quality
        f'--maf 0.05 '             #minor allele freq from 0.01 to 0.05
        f'--hwe 1e-6 '
        f'--geno 0.02 '            # SNPs missingness - exclude SNPs with missing genotype rate > 2%
        #f'--mind 0.1 '            # excludes individuals with a missing genotype rate higher than 10% (0.1).
        f'--rm-dup error list '    # Check each group of duplicate-ID variants for equality. (Alleles are considered unequal even if the codes are the same, just in a different order; FILTER/INFO are considered unequal if the strings don't match exactly, even if they're semantically identical.) If any mismatches are found, this errors out, and writes a list of mismatching variant IDs to plink2.rmdup.mismatch
        f'--write-snplist '
        f'--make-just-fam '
        f'--out {study}.QC')
################################################################################
# Pruning
# Better to used kb thresholds in indep-pairwise (more robust, less arbitrary)
# Remove variants that are in linkage disequilibrium (LD) with each other to a high degree
################################################################################

if not os.path.exists(f'{study}.QC.prune.in'):
    run(f'{plink} '
        f'--pfile ROSMAP_new.QC.merged '
        f'--keep {study}.QC.fam '
        f'--extract {study}.QC.snplist '
        f'--indep-pairwise 500kb 1 0.2 ' #Broad GWAS: Initial scans in GWAS are often done with higher r2 thresholds (0.2-0.5) to reduce the computational burden and focus on stronger genetic signals. For 0.2 threshold: 9409752/10213410 variants removed; for 0.5: 8505886/10213410 variants removed. fpr 0/8:  7273297/10213410 variants removed.
        f'--out {study}.QC')


################################################################################
# Calculate heterozygosity rates
################################################################################

if not os.path.exists(f'{study}.QC.het'):
    run(f'{plink} '
        f'--pfile ROSMAP_new.QC.merged '
        f'--extract {study}.QC.prune.in '
        f'--keep {study}.QC.fam '
        f'--het '
        f'--out {study}.QC')
    
################################################################################
# Remove individuals with heterozogisity F coefficients that are more than 3 standard
# deviation (SD) units from the mean
################################################################################

if not os.path.exists(f'{study}.QC.valid'):
    data = pd.read_table(f'{study}.QC.het')
    data[["m"]] = data[["F"]].mean()
    data[["s"]] = data[["F"]].std()
    valid = data.query('F <= m+3*s and F >= m-3*s').iloc[:, [0, 1]]
    #valid = data.query('F <= 0.1 and F >= -0.1').iloc[:, [0, 1]]
    valid.to_csv(f'{study}.valid.sample', index=False) # use if performing sex check
    #valid.to_csv(f'{study}.QC.valid', sep="\t", header=False, index=False) # use if not performing sex check
################################################################################
################################################################################
# Sex check
# Can only perform sex check if X chromosome is included else skip this step
# TODO: need to exclude pseudoautosomal regions from X chromosome before doing sex check.
# TODO: Need to manually look this up in genome that you are using.
################################################################################
# #transfer pgen to bed files for plink1.9
# if not os.path.exists(f'{study}.QC.new') and sex:
#     run(f'{plink} '
#         f'--pfile ROSMAP_new.QC.merged '
#         f'--extract {study}.QC.prune.in '
#         f'--keep {study}.valid.sample '
#         f'--make-bed '
#         f'--out {study}.QC.new')
    
# #use plink1.9 for sex check
# plink1_9 = '/external/rprshnas01/kcni/mwainberg/software/plink'
if not os.path.exists(f'{study}.QC.valid') and sex:
    run(f'{plink} '
        f'--pfile ROSMAP_new.QC.merged '
        f'--extract {study}.QC.prune.in '
        f'--keep {study}.valid.sample '
        f'--keep-females '
        f'--out {study}.QC.female')
    

if not os.path.exists(f'{study}.QC.valid') and sex:
    run(f'{plink} '
        f'--pfile ROSMAP_new.QC.merged '
        f'--extract {study}.QC.prune.in '
        f'--keep {study}.valid.sample '
        f'--keep-males '
        f'--out {study}.QC.male')

    valid_FID = pd.read_table(f'{study}.valid.sample')["FID"] #valid_FID = pd.read_table(f'{study}.valid.sample')["#IID"] CHECK COLNAMES in valid file
    sex_data = pd.read_table(f'{study}.QC.sexcheck')
    valid_sex = sex_data.query('STATUS=="OK" and FID in @valid_DIF)').filter(["FID", "IID"]) # Check colnames in files
    valid_sex.to_csv(f'{study}.QC.valid', sep="\t", header=False, index=False)

################################################################################
# Generate final QC'ed target data file with relatedness filter
# Can take care of relatedness during GWAS, but we follow the method mentioned in the thesis here
# Exclude --keep step if sex check wasn't performed
################################################################################

if not os.path.exists(f'{study}.QC.final.pgen') and relatedness:
    run(f'{plink} '
        f'--pfile ROSMAP_new.QC.merged '
        f'--make-pgen '
        #f'--extract {study}.QC.prune.in ' # in the final pgen file, we want all vars, not only the pruned ones, so we exclude this line
        f'--keep {study}.QC.valid ' 
        f'--out {study}.QC.final '
        f'--extract {study}.QC.snplist')
    
################################################################################
# Calculate PCA
# Plink is slow but fine to use for small cohorts
# Technically need to reprune first, but likely won't change much
# Do PCs on {study}.QC
################################################################################
if not os.path.exists(f'{study}.QC.eigenvec'):
    run(f'{plink} '
        f'--pfile {study}.QC.final '
        f'--extract {study}.QC.prune.in '
        f'--pca '
        f'--out {study}.QC')

pca = pd.read_table(f'{study}.QC.eigenvec', header=None).tail(-1)
pca.columns = ["FID", "IID"] + [f'PC{str(e)}' for e in range(1, len(pca.columns) -1)]
pca.to_csv('pca_new.csv')