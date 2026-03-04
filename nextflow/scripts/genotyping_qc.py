#!/usr/bin/env python3
"""
Parameterized Genotyping QC Script
Supports WGS and SNP array data for multiple cohorts
"""

import argparse
import os
import subprocess
import pandas as pd
import sys
from pathlib import Path

def run(cmd, pipefail=True, **kwargs):
    """Run a shell command"""
    return subprocess.run(f'set -eu{"o pipefail" if pipefail else ""}; {cmd}',
                          check=True, shell=True, executable='/bin/bash',
                          **kwargs)

def run_slurm(cmd, job_name, time, log_file, num_threads=6, mem_per_cpu='16000M', num_nodes=1):
    """Submit a SLURM job"""
    assert ' ' not in job_name
    from tempfile import NamedTemporaryFile
    try:
        with NamedTemporaryFile('w', dir='', suffix='.sh', delete=False) as temp_file:
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
            pass

def create_pgen(args):
    """Create pgen files from VCF"""
    print(f"Creating pgen files for {args.study}...")
    
    chroms = list(range(1, 23)) + ['X'] if args.include_x else list(range(1, 23))
    
    for chrom in chroms:
        # Determine input VCF
        if args.normalized_vcf_dir and os.path.exists(f"{args.normalized_vcf_dir}/{args.study}.QC.{chrom}.normalized.vcf.gz"):
            input_vcf = f"{args.normalized_vcf_dir}/{args.study}.QC.{chrom}.normalized.vcf.gz"
        elif args.vcf_pattern:
            input_vcf = args.vcf_pattern.format(chrom=chrom).strip('"')  # Remove quotes if present
        elif args.vcf_dir:
            input_vcf = f"{args.vcf_dir}/chr{chrom}.vcf.gz"
        else:
            print(f"  ERROR: Cannot determine input VCF for chromosome {chrom}")
            continue
        
        output_prefix = f"{args.study}.QC.{chrom}"
        
        # For CMC cohorts: subset VCF first if samples_to_keep is provided and we're not normalizing
        # (meaning we need to subset from already normalized VCFs)
        if args.samples_to_keep and not args.normalize_vcf:
            # Subset the VCF using bcftools
            subset_vcf = f"{output_prefix}.subset.vcf.gz"
            if not os.path.exists(subset_vcf):
                print(f"  Chromosome {chrom}: Subsetting VCF...")
                # bcftools view -S needs sample IDs only (no header), but samples_to_keep has "FID IID" header
                # Create temp file with just the IID column (2nd column), skipping header
                samples_no_header = f"{output_prefix}.samples_no_header.txt"
                run(f"tail -n +2 {args.samples_to_keep} | awk '{{print $2}}' > {samples_no_header}")
                # Use bcftools view -S to subset samples
                run(f'{args.bcftools_path} view -S {samples_no_header} -O z -o {subset_vcf} {input_vcf}')
            input_vcf = subset_vcf
        
        # Build plink command
        plink_cmd = f'{args.plink_path} --vcf {input_vcf}'
        
        # Add VCF quality filters if specified
        if args.vcf_min_gq:
            plink_cmd += f' --vcf-min-gq {args.vcf_min_gq}'
        if args.vcf_min_dp:
            plink_cmd += f' --vcf-min-dp {args.vcf_min_dp}'
        if args.vcf_min_qual:
            plink_cmd += f' --vcf-min-qual {args.vcf_min_qual}'
        
        # Add cohort-specific options
        if args.sort_vars:
            plink_cmd += ' --sort-vars'
        if args.autosome_only:
            plink_cmd += ' --autosome'
        
        plink_cmd += f' --double-id --make-pgen --new-id-max-allele-len 2000 --set-all-var-ids chr@:#:\\$r:\\$a --out {output_prefix}'
        
        run(plink_cmd)
        print(f"  Chromosome {chrom}: pgen created")

def create_pgen_single_chr(args):
    """Create pgen file for a single chromosome (for parallel processing)"""
    chrom = args.chromosome
    print(f"Creating pgen file for {args.study} chromosome {chrom}...")
    
    # Determine input VCF
    if args.normalized_vcf_dir and os.path.exists(f"{args.normalized_vcf_dir}/{args.study}.QC.{chrom}.normalized.vcf.gz"):
        input_vcf = f"{args.normalized_vcf_dir}/{args.study}.QC.{chrom}.normalized.vcf.gz"
    elif args.vcf_pattern:
        input_vcf = args.vcf_pattern.format(chrom=chrom).strip('"')  # Remove quotes if present
    elif args.vcf_dir:
        input_vcf = f"{args.vcf_dir}/chr{chrom}.vcf.gz"
    else:
        print(f"  ERROR: Cannot determine input VCF for chromosome {chrom}")
        return
    
    output_prefix = f"{args.study}.QC.{chrom}"
    use_plink_keep = False
    
    # When samples_to_keep is provided and we're not normalizing: try bcftools subset first.
    # If bcftools fails (e.g. NABEC VCF has contig/tag header issues), fall back to using
    # the original VCF with plink --keep instead.
    def _subset_vcf_valid(path):
        """Check if subset VCF exists and has valid header (not truncated from failed bcftools)"""
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            return False
        try:
            r = subprocess.run(
                f'{args.bcftools_path} view -h {path}',
                shell=True, capture_output=True, text=True, timeout=30
            )
            return r.returncode == 0 and '#CHROM' in (r.stdout or '')
        except Exception:
            return False

    if args.samples_to_keep and not args.normalize_vcf:
        subset_vcf = f"{output_prefix}.subset.vcf.gz"
        if os.path.exists(subset_vcf) and _subset_vcf_valid(subset_vcf):
            # Reuse existing valid subset
            input_vcf = subset_vcf
        else:
            if os.path.exists(subset_vcf):
                os.remove(subset_vcf)  # Remove corrupt/truncated file from previous run
            print(f"  Chromosome {chrom}: Subsetting VCF with bcftools...")
            samples_no_header = f"{output_prefix}.samples_no_header.txt"
            run(f"tail -n +2 {args.samples_to_keep} | awk '{{print $2}}' > {samples_no_header}")
            try:
                run(f'{args.bcftools_path} view -S {samples_no_header} -O z -o {subset_vcf} {input_vcf}')
                input_vcf = subset_vcf
            except subprocess.CalledProcessError as e:
                print(f"  Chromosome {chrom}: bcftools subset failed ({e}), using original VCF with plink --keep")
                if os.path.exists(subset_vcf):
                    os.remove(subset_vcf)  # Remove truncated output
                use_plink_keep = True
    
    # Build plink command
    plink_cmd = f'{args.plink_path} --vcf {input_vcf}'
    
    # Add VCF quality filters if specified
    if args.vcf_min_gq:
        plink_cmd += f' --vcf-min-gq {args.vcf_min_gq}'
    if args.vcf_min_dp:
        plink_cmd += f' --vcf-min-dp {args.vcf_min_dp}'
    if args.vcf_min_qual:
        plink_cmd += f' --vcf-min-qual {args.vcf_min_qual}'
    
    # Add cohort-specific options
    if args.sort_vars:
        plink_cmd += ' --sort-vars'
    if args.autosome_only:
        plink_cmd += ' --autosome'
    if use_plink_keep:
        plink_cmd += f' --keep {args.samples_to_keep}'
    
    plink_cmd += f' --double-id --make-pgen --new-id-max-allele-len 2000 --set-all-var-ids chr@:#:\\$r:\\$a --out {output_prefix}'
    
    run(plink_cmd)
    print(f"  Chromosome {chrom}: pgen created successfully")

def fix_psam_format(pgen_prefix):
    """
    Fix .psam file to ensure it has both #FID and IID columns.
    This is required for --keep to work properly in plink2.
    Backward compatible with all cohorts.
    """
    psam_file = f'{pgen_prefix}.psam'
    
    if not os.path.exists(psam_file):
        print(f"  Warning: {psam_file} not found, skipping psam fix")
        return
    
    # Read the current psam file
    psam = pd.read_table(psam_file)
    
    # Check if #FID column already exists (backward compatibility)
    if '#FID' in psam.columns:
        print(f"  {psam_file} has #FID column, verifying format...")
        
        # Only check that the essential #FID and IID columns have non-null values.
        # Other columns like PAT, MAT, SEX may legitimately be NA/0 and must NOT
        # trigger a false "corrupted" detection (which would restore a stale backup
        # with a different sample count than the pgen).
        has_iid = 'IID' in psam.columns
        first_row_ok = len(psam) > 0 and has_iid and pd.notna(psam['#FID'].iloc[0]) and pd.notna(psam['IID'].iloc[0])
        if has_iid and first_row_ok:
            print(f"  {psam_file} format verified, no fix needed")
            normalize_duplicated_ids_in_psam(pgen_prefix)
            return
        else:
            print(f"  {psam_file} appears corrupted (missing IID or null FID/IID in first row)")
            return
    
    # Check if we have #IID column
    if '#IID' not in psam.columns:
        print(f"  Warning: {psam_file} has unexpected format, skipping fix")
        return
    
    print(f"  Fixing {psam_file} to include #FID column...")
    
    # Backup original file BEFORE any modifications
    backup_file = f'{psam_file}.backup'
    if not os.path.exists(backup_file):
        import shutil
        shutil.copy2(psam_file, backup_file)
        print(f"  Original file backed up to {backup_file}")
    
    # Create proper format with #FID, IID, and preserve all other columns
    # Build column dictionary in correct order
    new_columns = {
        '#FID': psam['#IID'],
        'IID': psam['#IID']
    }
    
    # Add all other columns (like SEX) from original file
    for col in psam.columns:
        if col not in ['#IID', '#FID', 'IID']:
            new_columns[col] = psam[col]
    
    # Ensure SEX column exists (required by plink2)
    if 'SEX' not in new_columns:
        new_columns['SEX'] = 'NA'
    
    # Create DataFrame with explicit column order: #FID, IID, SEX, then any others
    column_order = ['#FID', 'IID', 'SEX'] + [c for c in new_columns.keys() if c not in ['#FID', 'IID', 'SEX']]
    psam_fixed = pd.DataFrame({col: new_columns[col] for col in column_order})
    
    # CRITICAL: Replace NaN values with "NA" string to prevent empty fields in TSV
    # When pandas writes NaN to TSV, it creates empty fields which breaks plink2
    for col in psam_fixed.columns:
        psam_fixed[col] = psam_fixed[col].fillna('NA')
    
    # Verify we have the right structure
    print(f"  New columns: {list(psam_fixed.columns)}")
    print(f"  Rows: {len(psam_fixed)}, Columns: {len(psam_fixed.columns)}")
    
    # Save the fixed file with proper formatting
    # Use na_rep='NA' as additional safety to ensure any remaining NaN values are written as "NA"
    psam_fixed.to_csv(psam_file, sep='\t', index=False, na_rep='NA')
    
    # Verify the saved file
    import subprocess
    try:
        # Check field counts in saved file
        result = subprocess.run(['awk', '{print NF}', psam_file], 
                              capture_output=True, text=True, check=True)
        field_counts = set(result.stdout.strip().split('\n'))
        if len(field_counts) == 1:
            print(f"  ✓ Fixed {psam_file}: {len(psam_fixed)} samples, {len(psam_fixed.columns)} columns per row")
        else:
            print(f"  ⚠ Warning: Field count mismatch detected: {field_counts}")
    except:
        print(f"  Fixed {psam_file}: {len(psam_fixed)} samples with {len(psam_fixed.columns)} columns")
        print(f"  (Could not verify field counts - awk not available)")

    # Normalize duplicated IDs (X_X -> X) for phenotype matching (e.g. NABEC)
    normalize_duplicated_ids_in_psam(pgen_prefix)


def normalize_duplicated_ids_in_psam(pgen_prefix):
    """
    Normalize sample IDs in psam when FID==IID and format is X_X (identical parts).
    E.g. KEN-1066_KEN-1066 -> KEN-1066. Required for phenotype/covariate matching
    when VCF/PLINK produces doubled IDs but biospecimen mapping uses single IDs.
    """
    psam_file = f'{pgen_prefix}.psam'
    if not os.path.exists(psam_file):
        return
    psam = pd.read_table(psam_file)
    fid_col = '#FID' if '#FID' in psam.columns else 'FID'
    iid_col = 'IID' if 'IID' in psam.columns else (psam.columns[1] if len(psam.columns) > 1 else None)
    if iid_col is None:
        return

    def _norm_id(s):
        s = str(s).strip()
        if '_' in s:
            parts = s.split('_', 1)
            if len(parts) == 2 and parts[0] == parts[1]:
                return parts[0]
        return s

    orig_fid = psam[fid_col].astype(str)
    orig_iid = psam[iid_col].astype(str)
    new_fid = orig_fid.apply(_norm_id)
    new_iid = orig_iid.apply(_norm_id)
    n_changed = ((new_fid != orig_fid) | (new_iid != orig_iid)).sum()
    if n_changed > 0:
        psam[fid_col] = new_fid
        psam[iid_col] = new_iid
        psam.to_csv(psam_file, sep='\t', index=False, na_rep='NA')
        print(f"  Normalized {n_changed} duplicated IDs in {psam_file} (X_X -> X)")


def _keep_file_matching_psam(samples_to_keep, psam_path, out_path):
    """Create a keep file with IDs matching the psam format (handles single vs double-ID mismatch)."""
    psam = pd.read_table(psam_path)
    fid_col = '#FID' if '#FID' in psam.columns else 'FID'
    iid_col = 'IID'

    def norm(x):
        s = str(x).strip()
        if '_' in s:
            parts = s.split('_', 1)
            if len(parts) == 2 and parts[0] == parts[1]:
                return parts[0]
        return s

    psam_map = {}
    for _, row in psam.iterrows():
        fid, iid = str(row[fid_col]).strip(), str(row[iid_col]).strip()
        psam_map[norm(iid)] = (fid, iid)
        psam_map[iid] = (fid, iid)

    keep_df = pd.read_table(samples_to_keep)
    fid_k = 'FID' if 'FID' in keep_df.columns else keep_df.columns[0]
    iid_k = 'IID' if 'IID' in keep_df.columns else keep_df.columns[1]
    matched = []
    for _, row in keep_df.iterrows():
        our_id = str(row[iid_k]).strip()
        key = norm(our_id)
        if key in psam_map:
            matched.append(psam_map[key])
        elif our_id in psam_map:
            matched.append(psam_map[our_id])
    if not matched:
        return False
    pd.DataFrame(matched, columns=['FID', 'IID']).to_csv(out_path, sep='\t', index=False)
    return True

def merge_pgen(args):
    """
    Merge pgen files across chromosomes.
    
    Strategy: Try PLINK2 --pmerge-list first (fast), fall back to bcftools if it fails.
    PLINK2 --pmerge-list has known issues with non-concatenating merges, but works for some cohorts.
    The bcftools approach is slower but more reliable for all cohorts.
    """
    print(f"Merging pgen files for {args.study}...")
    print(f"  Current directory: {os.getcwd()}")
    
    merged_prefix = f"{args.study}.QC.merged"
    
    # Determine chromosome range
    if args.include_x:
        chroms = list(range(1, 23)) + ['X']
    else:
        chroms = list(range(1, 23))
    
    # Find all existing pgen files
    existing_chroms = []
    for chrom in chroms:
        pgen_file = f"{args.study}.QC.{chrom}.pgen"
        if os.path.exists(pgen_file):
            existing_chroms.append(chrom)
        else:
            print(f"  WARNING: Missing pgen file for chromosome {chrom}")
    
    if len(existing_chroms) < 2:
        print(f"  ERROR: Found only {len(existing_chroms)} pgen file(s). Need at least 2.")
        sys.exit(1)
    
    print(f"  Found {len(existing_chroms)} pgen files to merge")
    
    # Try Method 1: PLINK2 --pmerge-list (fast but may fail)
    print("  Attempting Method 1: PLINK2 --pmerge-list...")
    merge_list = f"{args.study}.merge-list.txt"
    
    # Create merge list
    with open(merge_list, 'w') as f:
        for chrom in existing_chroms:
            f.write(f"{args.study}.QC.{chrom}\n")
    
    try:
        # Try PLINK2 merge
        run(f'{args.plink_path} --pmerge-list {merge_list} --out {merged_prefix}')
        print("  ✓ Method 1 succeeded: PLINK2 --pmerge-list")
        # Clean up merge list
        if os.path.exists(merge_list):
            os.remove(merge_list)
        # Fix psam format to ensure #FID column exists
        fix_psam_format(merged_prefix)
        return
    except subprocess.CalledProcessError as e:
        print(f"  ✗ Method 1 failed: PLINK2 --pmerge-list returned error {e.returncode}")
        print(f"  → Falling back to Method 2: bcftools concat (more reliable)")
    
    # Method 2: bcftools concat approach (slower but always works)
    # This is the proven method from the user's own NABEC QC script
    
    # Step 1: Convert each pgen to VCF
    print("  Step 1/5: Converting pgen files to VCF...")
    for chrom in existing_chroms:
        vcf_file = f"{args.study}.QC.{chrom}.vcf.gz"
        if not os.path.exists(vcf_file):
            print(f"    Converting chromosome {chrom}...")
            run(f'{args.plink_path} '
                f'--pfile {args.study}.QC.{chrom} '
                f'--export vcf bgz '
                f'--out {args.study}.QC.{chrom}')
    
    # Step 2: Index all VCF files (required for bcftools concat)
    print("  Step 2/5: Indexing VCF files...")
    for chrom in existing_chroms:
        vcf_file = f"{args.study}.QC.{chrom}.vcf.gz"
        csi_file = f"{vcf_file}.csi"
        if not os.path.exists(csi_file):
            print(f"    Indexing chromosome {chrom}...")
            run(f'{args.bcftools_path} index {vcf_file}')
    
    # Step 3: Create VCF list for concatenation in proper order
    vcf_list = f"{args.study}.vcf.list.txt"
    with open(vcf_list, 'w') as f:
        for chrom in existing_chroms:
            f.write(f"{args.study}.QC.{chrom}.vcf.gz\n")
    
    # Step 4: Concatenate VCFs with bcftools
    print("  Step 3/5: Concatenating VCFs with bcftools...")
    merged_vcf = f"{merged_prefix}.vcf.gz"
    run(f'{args.bcftools_path} concat '
        f'--file-list {vcf_list} '
        f'-a '  # Allow overlaps (important for some cohorts)
        f'-Oz -o {merged_vcf}')
    
    # Step 5: Convert merged VCF back to pgen (with --sort-vars so chromosomes are
    # in ascending numeric order, which is required by downstream tools like REGENIE)
    print("  Step 4/5: Converting merged VCF back to pgen...")
    run(f'{args.plink_path} '
        f'--vcf {merged_vcf} '
        f'--make-pgen '
        f'--sort-vars '
        f'--out {merged_prefix}')
    
    # Step 6: Clean up intermediate files to save disk space
    print("  Step 5/5: Cleaning up intermediate files...")
    for chrom in existing_chroms:
        vcf_file = f"{args.study}.QC.{chrom}.vcf.gz"
        csi_file = f"{vcf_file}.csi"
        if os.path.exists(vcf_file):
            os.remove(vcf_file)
        if os.path.exists(csi_file):
            os.remove(csi_file)
    if os.path.exists(vcf_list):
        os.remove(vcf_list)
    if os.path.exists(merged_vcf):
        os.remove(merged_vcf)
    if os.path.exists(merge_list):
        os.remove(merge_list)
    
    print("  ✓ Method 2 succeeded: Merge completed using bcftools concat approach")
    print(f"  → Final output: {merged_prefix}.pgen")
    
    # Fix psam format to ensure #FID column exists
    fix_psam_format(merged_prefix)

def standard_qc(args):
    """Apply standard GWAS QC filters"""
    print(f"Applying standard QC filters for {args.study}...")
    
    qc_prefix = f"{args.study}.QC"
    
    plink_cmd = f'{args.plink_path} --pfile {args.study}.QC.merged'

    # Add imputation quality filter if specified
    if args.mach_r2_filter:
        plink_cmd += f' --mach-r2-filter {args.mach_r2_filter}'

    # Add sample filtering if specified and file exists
    if args.samples_to_keep and os.path.exists(args.samples_to_keep):
        merged_psam = f'{args.study}.QC.merged.psam'
        keep_file = args.samples_to_keep
        if os.path.exists(merged_psam):
            keep_converted = f'{args.study}.QC.keep_converted.txt'
            if _keep_file_matching_psam(args.samples_to_keep, merged_psam, keep_converted):
                keep_file = keep_converted
        plink_cmd += f' --keep {keep_file}'
    elif args.samples_to_keep:
        print(f"  WARNING: samples_to_keep file not found: {args.samples_to_keep}")
        print(f"  Continuing without sample filtering...")
    
    # Add QC thresholds
    plink_cmd += f' --maf {args.maf_threshold}'
    plink_cmd += f' --hwe {args.hwe_threshold}'
    plink_cmd += f' --geno {args.geno_threshold}'
    plink_cmd += f' --mind {args.mind_threshold}'
    
    plink_cmd += f' --write-snplist --make-just-fam --out {qc_prefix}'
    
    run(plink_cmd)
    
    # Copy QC files to output directory if specified (useful for debugging)
    if args.output_dir:
        import shutil
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Copy key QC files for reference
        for qc_file in [f"{qc_prefix}.fam", f"{qc_prefix}.snplist"]:
            if os.path.exists(qc_file):
                output_file = f"{args.output_dir}/{qc_file}"
                shutil.copy2(qc_file, output_file)
                print(f"  Copied {qc_file} to output directory")
    
    print("  Standard QC completed")

def ld_pruning(args):
    """Perform LD pruning"""
    print(f"Performing LD pruning for {args.study}...")
    
    prune_prefix = f"{args.study}.QC"
    
    run(f'{args.plink_path} '
        f'--pfile {args.study}.QC.merged '
        f'--keep {prune_prefix}.fam '
        f'--extract {prune_prefix}.snplist '
        f'--rm-dup exclude-all '  # Remove duplicate variant IDs before LD pruning
        f'--indep-pairwise {args.prune_window_size} {args.prune_step_size} {args.prune_r2_threshold} '
        f'--out {prune_prefix}')
    
    # Copy prune.in file to output directory if specified (needed by GWAS pipeline)
    if args.output_dir:
        import shutil
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Copy prune.in (required by GWAS pipeline)
        if os.path.exists(f"{prune_prefix}.prune.in"):
            output_prune_in = f"{args.output_dir}/{prune_prefix}.prune.in"
            shutil.copy2(f"{prune_prefix}.prune.in", output_prune_in)
            print(f"  Copied prune.in to: {output_prune_in}")
        
        # Also copy prune.out for completeness (useful for debugging)
        if os.path.exists(f"{prune_prefix}.prune.out"):
            output_prune_out = f"{args.output_dir}/{prune_prefix}.prune.out"
            shutil.copy2(f"{prune_prefix}.prune.out", output_prune_out)
            print(f"  Copied prune.out to: {output_prune_out}")
    
    print("  LD pruning completed")

def heterozygosity(args):
    """Calculate heterozygosity and filter individuals"""
    print(f"Calculating heterozygosity for {args.study}...")
    
    het_file = f"{args.study}.QC.het"
    valid_file = f"{args.study}.QC.valid"
    
    # Calculate heterozygosity
    run(f'{args.plink_path} '
        f'--pfile {args.study}.QC.merged '
        f'--extract {args.study}.QC.prune.in '
        f'--keep {args.study}.QC.fam '
        f'--het '
        f'--out {args.study}.QC')
    
    # Filter individuals (within 3 SD of mean)
    data = pd.read_table(het_file)
    data["m"] = data["F"].mean()
    data["s"] = data["F"].std()
    filtered = data.query('F <= m+3*s and F >= m-3*s')
    
    # Check if het file has FID and IID or just IID
    # If first column is '#FID', we have both FID and IID
    # If first column is '#IID', we only have IID and need to duplicate it for FID
    if data.columns[0] == '#FID':
        # Has FID and IID - select both
        valid = filtered.iloc[:, [0, 1]]
    else:
        # Only has IID - duplicate it for FID (PLINK --keep expects FID IID)
        valid = filtered.iloc[:, [0]]
        valid = pd.concat([valid, valid.iloc[:, 0]], axis=1)
    
    valid.to_csv(valid_file, sep="\t", header=False, index=False)
    
    # Copy valid file to output directory if specified (useful for debugging)
    if args.output_dir:
        import shutil
        os.makedirs(args.output_dir, exist_ok=True)
        
        if os.path.exists(valid_file):
            output_valid = f"{args.output_dir}/{valid_file}"
            shutil.copy2(valid_file, output_valid)
            print(f"  Copied {valid_file} to output directory")
        
        if os.path.exists(het_file):
            output_het = f"{args.output_dir}/{het_file}"
            shutil.copy2(het_file, output_het)
            print(f"  Copied {het_file} to output directory")
    
    print(f"  Heterozygosity filtering completed: {len(valid)}/{len(data)} samples retained")

def final_qc(args):
    """Generate final QC'd pgen files"""
    print(f"Generating final QC'd files for {args.study}...")
    
    # Ensure merged psam file has correct format (safety check)
    merged_prefix = f"{args.study}.QC.merged"
    fix_psam_format(merged_prefix)
    
    # Determine output location
    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
        final_prefix = f"{args.output_dir}/{args.study}.QC.final"
    else:
        final_prefix = f"{args.study}.QC.final"
    
    run(f'{args.plink_path} '
        f'--pfile {args.study}.QC.merged '
        f'--make-pgen '
        f'--sort-vars '
        f'--keep {args.study}.QC.valid '
        f'--out {final_prefix} '
        f'--extract {args.study}.QC.snplist')
    print(f"  Final QC files created: {final_prefix}.pgen/psam/pvar")

def calculate_pca(args):
    """Calculate PCA"""
    print(f"Calculating PCA for {args.study}...")
    
    # Determine final pfile location
    if args.output_dir and os.path.exists(f"{args.output_dir}/{args.study}.QC.final.pgen"):
        final_pfile = f"{args.output_dir}/{args.study}.QC.final"
    else:
        final_pfile = f"{args.study}.QC.final"
    
    eigenvec_file = f"{args.study}.QC.eigenvec"
    pca_csv = f"{args.output_dir}/pca.csv" if args.output_dir else "pca.csv"
    
    # Calculate PCA
    run(f'{args.plink_path} '
        f'--pfile {final_pfile} '
        f'--extract {args.study}.QC.prune.in '
        f'--pca '
        f'--out {args.study}.QC')
    
    # Convert to CSV
    pca = pd.read_table(eigenvec_file, header=None).tail(-1)
    pca.columns = ["FID", "IID"] + [f'PC{str(e)}' for e in range(1, len(pca.columns) - 1)]
    pca.to_csv(pca_csv, index=False)
    
    print(f"  PCA completed: {pca_csv}")

def main():
    parser = argparse.ArgumentParser(description='Genotyping QC Pipeline')
    parser.add_argument('--step', required=True, 
                       choices=['normalize', 'create_pgen', 'create_pgen_single_chr', 'merge', 'standard_qc', 
                               'ld_pruning', 'heterozygosity', 'final_qc', 'pca'],
                       help='QC step to execute')
    parser.add_argument('--study', required=True, help='Study name')
    parser.add_argument('--chromosome', type=str, help='Chromosome number (for create_pgen_single_chr)')
    parser.add_argument('--plink_path', default='/external/rprshnas01/kcni/mwainberg/software/plink2',
                       help='Path to plink2')
    parser.add_argument('--bcftools_path', default='/external/rprshnas01/netdata_kcni/stlab/Xiaolin/software/bcftools-1.12/bcftools',
                       help='Path to bcftools')
    
    # VCF paths
    parser.add_argument('--vcf_dir', help='Directory containing VCF files')
    parser.add_argument('--vcf_pattern', help='VCF file pattern with {chrom} placeholder')
    parser.add_argument('--normalized_vcf_dir', help='Directory for normalized VCF files')
    
    # QC thresholds
    parser.add_argument('--maf_threshold', type=float, default=0.05, help='MAF threshold')
    parser.add_argument('--hwe_threshold', default='1e-15', help='HWE threshold')
    parser.add_argument('--geno_threshold', type=float, default=0.1, help='SNP missingness threshold')
    parser.add_argument('--mind_threshold', type=float, default=0.1, help='Individual missingness threshold')
    parser.add_argument('--mach_r2_filter', type=float, help='Imputation quality filter (R2)')
    
    # VCF quality filters
    parser.add_argument('--vcf_min_gq', type=int, help='Minimum genotype quality')
    parser.add_argument('--vcf_min_dp', type=int, help='Minimum depth')
    parser.add_argument('--vcf_min_qual', type=int, help='Minimum variant quality')
    
    # Pruning parameters
    parser.add_argument('--prune_window_size', default='500kb', help='LD pruning window size')
    parser.add_argument('--prune_step_size', type=int, default=1, help='LD pruning step size')
    parser.add_argument('--prune_r2_threshold', type=float, default=0.2, help='LD pruning R2 threshold')
    
    # Options
    parser.add_argument('--include_x', action='store_true', help='Include X chromosome')
    parser.add_argument('--sort_vars', action='store_true', help='Sort variants (for NABEC)')
    parser.add_argument('--autosome_only', action='store_true', help='Only include autosomes')
    parser.add_argument('--samples_to_keep', help='File with samples to keep (FID IID format)')
    parser.add_argument('--normalize_vcf', action='store_true', help='Deprecated: normalization step removed, use pre-normalized VCFs')
    parser.add_argument('--use_slurm', action='store_true', help='Use SLURM for normalization')
    parser.add_argument('--slurm_time', default='0-4:0', help='SLURM time limit')
    parser.add_argument('--log_dir', help='Directory for log files')
    parser.add_argument('--output_dir', help='Output directory')
    parser.add_argument('--work_dir', help='Working directory')
    
    args = parser.parse_args()
    
    # Create directories if they don't exist
    if args.work_dir:
        os.makedirs(args.work_dir, exist_ok=True)
        os.chdir(args.work_dir)
    
    if args.log_dir:
        os.makedirs(args.log_dir, exist_ok=True)
    
    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
    
    # Execute requested step (normalize step removed - use pre-normalized VCFs)
    if args.step == 'normalize':
        print("  Normalization step skipped. Use pre-normalized VCFs (vcf_pattern or normalized_vcf_dir).")
    elif args.step == 'create_pgen':
        create_pgen(args)
    elif args.step == 'create_pgen_single_chr':
        create_pgen_single_chr(args)
    elif args.step == 'merge':
        merge_pgen(args)
    elif args.step == 'standard_qc':
        standard_qc(args)
    elif args.step == 'ld_pruning':
        ld_pruning(args)
    elif args.step == 'heterozygosity':
        heterozygosity(args)
    elif args.step == 'final_qc':
        final_qc(args)
    elif args.step == 'pca':
        calculate_pca(args)

if __name__ == '__main__':
    main()

