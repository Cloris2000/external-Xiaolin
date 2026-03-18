/*
 * Module: Genotyping QC Step
 * Runs a single step of the genotyping QC pipeline using the parameterized Python script
 */

process GENOTYPING_QC_STEP {
    label 'high_memory'
    tag "${step}_${study}${chromosome ? '_chr' + chromosome : ''}"
    
    input:
    val step
    val study
    path script_file
    path prev_output  // Input from previous step to create dependency
    val qc_version  // Parameter version to invalidate cache when QC params change
    val work_dir
    val chromosome  // Optional: for parallel chromosome processing
    val vcf_dir
    val vcf_pattern
    val normalized_vcf_dir
    val plink_path
    val bcftools_path
    val maf_threshold
    val hwe_threshold
    val geno_threshold
    val mind_threshold
    val mach_r2_filter
    val vcf_min_gq
    val vcf_min_dp
    val vcf_min_qual
    val vcf_dosage
    val prune_window_size
    val prune_step_size
    val prune_r2_threshold
    val include_x
    val sort_vars
    val autosome_only
    val samples_to_keep
    val normalize_vcf
    val use_slurm
    val slurm_time
    val log_dir
    val output_dir
    
    output:
    path "*.done", emit: output_files, optional: true
    
    script:
    // Note: genotyping_qc.py uses system Python which has pandas
    // We don't activate conda_env here to use the system Python installation
    
    def vcf_dir_arg = (vcf_dir && vcf_dir != '') ? "--vcf_dir ${vcf_dir}" : ""
    def vcf_pattern_arg = (vcf_pattern && vcf_pattern != '') ? "--vcf_pattern \"${vcf_pattern}\"" : ""
    def normalized_vcf_dir_arg = (normalized_vcf_dir && normalized_vcf_dir != '') ? "--normalized_vcf_dir ${normalized_vcf_dir}" : ""
    def mach_r2_arg = (mach_r2_filter && mach_r2_filter != '') ? "--mach_r2_filter ${mach_r2_filter}" : ""
    def vcf_min_gq_arg = (vcf_min_gq && vcf_min_gq != '') ? "--vcf_min_gq ${vcf_min_gq}" : ""
    def vcf_min_dp_arg = (vcf_min_dp && vcf_min_dp != '') ? "--vcf_min_dp ${vcf_min_dp}" : ""
    def vcf_min_qual_arg = (vcf_min_qual && vcf_min_qual != '') ? "--vcf_min_qual ${vcf_min_qual}" : ""
    def vcf_dosage_arg = (vcf_dosage && vcf_dosage != '') ? "--vcf_dosage ${vcf_dosage}" : ""
    def include_x_arg = include_x ? "--include_x" : ""
    def sort_vars_arg = sort_vars ? "--sort_vars" : ""
    def autosome_only_arg = autosome_only ? "--autosome_only" : ""
    def samples_to_keep_arg = (samples_to_keep && samples_to_keep != '' && !samples_to_keep.toString().contains('DataflowVariable')) ? "--samples_to_keep ${samples_to_keep}" : ""
    def normalize_vcf_arg = normalize_vcf ? "--normalize_vcf" : ""
    def use_slurm_arg = use_slurm ? "--use_slurm" : ""
    def log_dir_arg = (log_dir && log_dir != '') ? "--log_dir ${log_dir}" : ""
    def output_dir_arg = (output_dir && output_dir != '') ? "--output_dir ${output_dir}" : ""
    def chromosome_arg = (chromosome && chromosome != '') ? "--chromosome ${chromosome}" : ""
    
    """
    # QC Version: ${qc_version} (increment when QC parameters change to force re-run)
    # Use the python_env conda environment which has pandas installed.
    PYTHON_BIN="/nethome/kcni/xzhou/.anaconda3/envs/python_env/bin/python"
    echo "Using Python: \$PYTHON_BIN"
    "\$PYTHON_BIN" "${script_file}" \\
        --step ${step} \\
        --study ${study} \\
        --work_dir ${work_dir} \\
        --plink_path ${plink_path} \\
        --bcftools_path ${bcftools_path} \\
        --maf_threshold ${maf_threshold} \\
        --hwe_threshold ${hwe_threshold} \\
        --geno_threshold ${geno_threshold} \\
        --mind_threshold ${mind_threshold} \\
        --prune_window_size ${prune_window_size} \\
        --prune_step_size ${prune_step_size} \\
        --prune_r2_threshold ${prune_r2_threshold} \\
        ${vcf_dir_arg} \\
        ${vcf_pattern_arg} \\
        ${normalized_vcf_dir_arg} \\
        ${mach_r2_arg} \\
        ${vcf_min_gq_arg} \\
        ${vcf_min_dp_arg} \\
        ${vcf_min_qual_arg} \\
        ${vcf_dosage_arg} \\
        ${include_x_arg} \\
        ${sort_vars_arg} \\
        ${autosome_only_arg} \\
        ${samples_to_keep_arg} \\
        ${normalize_vcf_arg} \\
        ${use_slurm_arg} \\
        ${log_dir_arg} \\
        ${output_dir_arg} \\
        ${chromosome_arg}
    
    # Create completion marker
    touch ${step}${chromosome ? '_chr' + chromosome : ''}.done
    """
}

