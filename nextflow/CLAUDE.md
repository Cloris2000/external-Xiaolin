# CLAUDE.md

## Project context

This repository contains a research Nextflow pipeline for genomic analyses, including
phenotype generation, genotype QC/PCA, REGENIE GWAS/QTL analyses, METAL
meta-analysis, and downstream summary/plotting steps.

The project is running on the Trillium HPC cluster.

The repository may contain ongoing and uncommitted development work. Preserve
existing behavior unless I explicitly ask for a refactor.

---

## HPC environment

This project is accessed through Trillium via SSH.

Important rules:

- Never run compute-intensive analyses on a login node (`tri-login*`).
- Do not launch full Nextflow pipelines, REGENIE analyses, large R/Python jobs,
  genotype processing, or other heavy workloads directly on a login node.
- Computational work must use the appropriate Slurm allocation/job.
- Before submitting a Slurm job, show me the command or job script and explain
  what it will run.
- Do not cancel, modify, or resubmit existing jobs unless I explicitly ask.
- It is fine to run lightweight commands on login nodes, such as:
  - `ls`, `find`, `grep`, `git status`, `git diff`
  - inspecting small text/config files
  - syntax checks or very lightweight validation
- If unsure whether a command is computationally expensive, ask before running it.

The Trillium configuration may use a whole-node allocation where Nextflow tasks run
locally within an allocated compute node rather than submitting every process as a
separate Slurm job. Preserve this design unless I explicitly request changes.

---

## Filesystem and data safety

Treat research data as valuable and potentially irreplaceable.

- Never delete raw data.
- Never modify source genotype, RNA-seq, metadata, VCF/BCF/PGEN, or other primary
  datasets in place.
- Do not use destructive commands such as `rm -rf` unless I explicitly approve the
  exact paths.
- Do not overwrite existing analysis results unless I explicitly ask.
- Prefer writing new or intermediate outputs to the locations already defined by
  the pipeline/configuration.
- Do not move large datasets between `/project`, `/scratch`, or other filesystems
  unless explicitly instructed.
- Before changing paths, verify whether they refer to input data, intermediate data,
  temporary storage, or final results.
- Do not expose credentials, tokens, private keys, or restricted data.

---

## Editing policy

Before making substantial changes:

1. Inspect the relevant workflow, config, and scripts.
2. Explain briefly what you think the current code is doing.
3. Describe the change you intend to make.
4. Identify which files will be modified.
5. Then make the smallest reasonable change.

For small obvious fixes, you may edit directly unless the change could affect
scientific results or pipeline behavior.

Prefer minimal, targeted edits over large rewrites.

Do not refactor unrelated code while fixing a specific problem.

Do not rename files, processes, channels, parameters, or output directories without
a clear reason.

Preserve backward compatibility with existing cohort configurations where practical.

---

## Scientific correctness

This is a research analysis pipeline. Scientific correctness is more important than
code elegance.

Never silently change:

- phenotype definitions
- sample inclusion/exclusion criteria
- covariates
- genotype QC thresholds
- allele handling
- variant filtering
- ancestry definitions
- REGENIE options
- meta-analysis models
- statistical transformations
- significance thresholds

If a requested code change could alter scientific interpretation or numerical
results, point this out before making the change.

Do not invent missing cohort metadata, sample mappings, parameters, file paths, or
analysis results.

If something is ambiguous, inspect the existing code/configuration first and then
ask me if necessary.

---

## Nextflow conventions

When modifying Nextflow:

- Follow the style and structure already used in this repository.
- Check both the `.nf` workflow and relevant `.config` files before changing process
  behavior.
- Preserve channel semantics and process dependencies.
- Avoid hard-coding cohort-specific paths into reusable workflow code.
- Prefer configuration parameters for cohort-specific paths and settings.
- Be careful with Nextflow variable scope, tuple structure, file staging, and
  process output declarations.
- Consider resume behavior and caching when changing process definitions.
- Do not remove `publishDir`, resource specifications, container/environment
  settings, or error-handling behavior without explaining why.
- When debugging, trace the problem to the earliest incorrect input/process rather
  than patching downstream output blindly.

---

## R, Python, shell, and genomics tools

Use the software/environment already established by the project whenever possible.

Do not install new packages or modify environments without asking first.

For shell scripts:

- quote paths and variables where appropriate
- use explicit error handling for important steps
- avoid destructive wildcard operations

For R/Python:

- prefer reproducible scripts over interactive one-off transformations
- preserve sample IDs and variant IDs exactly unless transformation is required
- explicitly check joins/merges for unexpected sample loss or duplication

For genomic files:

- pay particular attention to genome build, chromosome naming, REF/ALT orientation,
  effect allele, sample IDs, and variant IDs
- never assume two files use the same variant representation without checking
- report unexpected sample or variant count changes

---

## Git

The repository may contain uncommitted work.

Before significant modifications, inspect:

`git status`

and, when relevant:

`git diff`

Do not discard, reset, stash, or overwrite my existing changes unless I explicitly
ask.

Do not run:

`git reset --hard`
`git clean -fd`
`git checkout -- <file>`
`git restore <file>`

without explicit approval.

Do not commit or push changes unless I explicitly ask you to.

When finished with an editing task, summarize the files changed and the important
changes made.

---

## Validation

After editing code, perform the lightest useful validation that is safe on the
current machine.

Examples include:

- syntax/config checks
- inspecting generated commands
- checking Nextflow DAG/config interpretation when lightweight
- reviewing `git diff`
- testing small toy/example inputs if available

Do not launch a full dataset analysis merely to validate a code change.

If proper validation requires compute resources, tell me what should be submitted
and propose the appropriate test.

---

## Communication style

Be concise and technical.

When debugging, distinguish clearly between:

- what you observed in the code/logs
- what you infer
- what you are uncertain about
- what change you recommend

Do not claim a problem is fixed until there is evidence supporting that conclusion.

When reporting an error, include the relevant file/process and the likely root cause.

When multiple solutions exist, prefer the smallest and safest change first.

