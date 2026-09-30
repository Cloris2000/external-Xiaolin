#!/usr/bin/env python3
"""Does the ROSMAP phenotype file join to the analyzed genotype .psam?"""
import sys

PSAM = "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/results/ROSMAP/ROSMAP.QC.final.psam"
PHENO = ("/scratch/zhoux156/nf_work/cohort_gwas_v2/ROSMAP/69/"
         "bd77bf9584bcb61c44ed1e2ab1b675/phenotypes_RINT.txt")

geno = set()
with open(PSAM) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        geno.add(line.rstrip("\n").split("\t")[1])

pheno = set()
with open(PHENO) as fh:
    header = next(fh)
    for line in fh:
        pheno.add(line.rstrip("\n").split("\t")[1])

print(f"PSAM  : {PSAM}\n        N={len(geno)}  sample={sorted(geno)[:3]}")
print(f"PHENO : {PHENO}\n        N={len(pheno)}  sample={sorted(pheno)[:3]}")
print(f"\nintersection = {len(geno & pheno)}")
print(f"psam-only    = {len(geno - pheno)}")
print(f"pheno-only   = {len(pheno - geno)}")
