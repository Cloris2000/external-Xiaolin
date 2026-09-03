#!/usr/bin/env python3
"""Rewrite SCC-only paths and SLURM settings in this repo for Trillium.

The pipeline was developed on CAMH's SCC and carries absolute paths to that
cluster's filesystem, Anaconda tree and partition names.  This rewrites them to
their Trillium equivalents and points the shell drivers at site_env.sh.

Run from the repo root:
    python3 scripts/port_scc_paths_to_trillium.py --dry-run
    python3 scripts/port_scc_paths_to_trillium.py --apply

Deliberately NOT handled here:
  * Raw per-cohort genotype/expression/metadata paths in nextflow.config.combined.*
    (VCFs, biospecimen CSVs, count matrices).  Those feed the per-cohort GWAS
    stages, which are not being re-run - the meta-analysis consumes the existing
    REGENIE outputs under $SCC_DIR/results.  Porting them needs a per-dataset
    mapping onto /project/rrg-shreejoy/{COHORT}/ and is a separate job.
  * SLURM job geometry (--mem, --array, --cpus-per-task).  Trillium allocates
    whole nodes, so the array drivers are repacked separately.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

NF_DIR = "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow"
REFS = "/project/rrg-shreejoy/pipeline_refs"
SCC_ROOT = "/project/rrg-shreejoy/zhoux156/Xiaolin/SCC"
# Spelled out rather than $HOME so the same substitution is valid in R, Python
# and Nextflow config files, where the shell would not expand it.
CONDA_ROOT = "/home/zhoux156/miniforge3"
SITE_ENV_LINE = f'source "${{SITE_ENV:-{NF_DIR}/site_env.sh}}"'

# Longest / most specific first: the bare ".../Xiaolin/nextflow" prefix would
# otherwise swallow paths that live outside the checkout.
PATH_MAP: list[tuple[str, str]] = [
    ("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/METAL/generic-metal/executables",
     f"{REFS}/tools/metal"),
    ("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/new_MTGnCgG_lfct2.5_Publication.csv",
     f"{REFS}/markers/new_MTGnCgG_lfct2.5_Publication.csv"),
    ("/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/hgnc_complete_set.txt",
     f"{REFS}/markers/hgnc_complete_set.txt"),
    ("/external/rprshnas01/kcni/mwainberg/software/regenie", f"{REFS}/tools/regenie"),
    ("/external/rprshnas01/kcni/mwainberg/software/plink2", f"{REFS}/tools/plink2"),
    ("/external/rprshnas01/kcni/mwainberg/ldsc/eur_w_ld_chr", f"{REFS}/ldsc/eur_w_ld_chr"),
    ("/external/rprshnas01/kcni/mwainberg/ldsc/w_hm3.snplist",
     f"{REFS}/ldsc/eur_w_ld_chr/w_hm3.snplist"),
    # ldsc.py / munge_sumstats.py are a git checkout, not conda entry points.
    ("/external/rprshnas01/kcni/mwainberg/ldsc", f"{REFS}/tools/ldsc"),
    ("/external/rprshnas01/netdata_kcni/stlab/cross_cohort_MGPs",
     "/project/rrg-shreejoy/cross_cohort_MGPs"),
    ("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow", NF_DIR),
    ("/nethome/kcni/xzhou/GWAS_tut/LDSC_Euro_ref", f"{SCC_ROOT}/ref_data/LDSC_Euro_ref"),
    ("/nethome/kcni/xzhou/.local/bin/nextflow", "nextflow"),
    ("/nethome/kcni/xzhou/.anaconda3/bin/Rscript", "Rscript"),
    ("/nethome/kcni/xzhou/.anaconda3/bin/python3", "python3"),
    ("/nethome/kcni/xzhou/.anaconda3/bin/python", "python"),
    ("/nethome/kcni/xzhou/.anaconda3", CONDA_ROOT),
    ("$HOME/.anaconda3", "$HOME/miniforge3"),
    ("~/.anaconda3", "~/miniforge3"),
    # R_LIBS pointing into a standalone SCC Anaconda tree; on Trillium R resolves
    # against the active conda env, so any survivor should point at `test`.
    ("/external/rprshnas01/kcni/xzhou/.anaconda3/lib/R/library",
     f"{CONDA_ROOT}/envs/test/lib/R/library"),
    ("/external/rprshnas01/kcni/xzhou/.anaconda3", CONDA_ROOT),
]

CODE_SUFFIXES = {".sbatch", ".sh", ".R", ".r", ".py", ".nf"}
CONFIG_PREFIX = "nextflow.config"
SKIP_DIRS = {".git", "work", "results", "run_dirs", ".nextflow", ".specstory", "_payload"}

SHELL_SUFFIXES = {".sbatch", ".sh"}

# Environment plumbing that site_env.sh now owns.
DROP_LINE_PATTERNS = [
    re.compile(r'^\s*export\s+R_LIBS=.*anaconda3.*$'),
    re.compile(r'^\s*export\s+PATH="?/nethome/kcni/xzhou/\.anaconda3/bin:\$PATH"?\s*$'),
    re.compile(r'^\s*source\s+.*[/~]\.?anaconda3/etc/profile\.d/conda\.sh\s*$'),
    re.compile(r'^\s*export\s+JAVA_HOME=.*anaconda3.*$'),
]

SCRATCH_ROOT = "/scratch/zhoux156"

# Compute nodes mount /project and $HOME read-only, so SLURM cannot create job
# log files anywhere except /scratch.  #SBATCH lines are parsed before the script
# runs, so these have to be literal paths rather than $LOG_ROOT.
LOG_RE = re.compile(r'^(#SBATCH\s+--(?:output|error)=)(\S+)\s*$')


def port_log_path(path: str) -> str:
    for prefix, repl in (
        (f"{NF_DIR}/logs/", f"{SCRATCH_ROOT}/logs/"),
        (f"{NF_DIR}/results/", f"{SCRATCH_ROOT}/results/"),
        (f"{NF_DIR}/manuscript_figure/", f"{SCRATCH_ROOT}/logs/manuscript_figure/"),
        ("logs/", f"{SCRATCH_ROOT}/logs/"),
    ):
        if path.startswith(prefix):
            return repl + path[len(prefix):]
    return path


# SCC partitions do not exist on Trillium.
PARTITION_RE = re.compile(r'^(#SBATCH\s+--partition=)(short|medium|mediumtmp|long)\s*$')
ACCOUNT_RE = re.compile(r'^#SBATCH\s+--account=')
SBATCH_RE = re.compile(r'^#SBATCH\s')

ASSIGN_RE = re.compile(
    r'^(?P<indent>\s*)(?P<var>NF_DIR|PROJECT_DIR|BASE_DIR|NEXTFLOW_DIR)='
    r'"?\'?' + re.escape(NF_DIR) + r'"?\'?\s*$'
)


def iter_files(root: Path):
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.is_symlink():
            continue
        if any(part in SKIP_DIRS for part in path.relative_to(root).parts):
            continue
        if path.suffix in CODE_SUFFIXES or path.name.startswith(CONFIG_PREFIX):
            yield path


def port_text(text: str) -> str:
    for old, new in PATH_MAP:
        text = text.replace(old, new)
    return text


def port_shell(lines: list[str], rel: str) -> list[str]:
    """Point a shell driver at site_env.sh and fix its SLURM directives."""
    out: list[str] = []
    sourced = False

    for line in lines:
        stripped = line.rstrip("\n")

        if any(p.match(stripped) for p in DROP_LINE_PATTERNS):
            continue

        m = PARTITION_RE.match(stripped)
        if m:
            out.append(f"{m.group(1)}compute\n")
            continue

        m = LOG_RE.match(stripped)
        if m:
            out.append(f"{m.group(1)}{port_log_path(m.group(2))}\n")
            continue

        m = ASSIGN_RE.match(stripped)
        if m:
            # The first such assignment becomes the site_env source; site_env.sh
            # already exports NF_DIR, so re-assigning it would defeat overrides.
            if not sourced:
                out.append(f"{m.group('indent')}{SITE_ENV_LINE}\n")
                sourced = True
            if m.group("var") != "NF_DIR":
                out.append(f'{m.group("indent")}{m.group("var")}="$NF_DIR"\n')
            continue

        out.append(line)

    # Alliance clusters require an accounting group on every job.
    if any(SBATCH_RE.match(l) for l in out) and not any(ACCOUNT_RE.match(l) for l in out):
        last = max(i for i, l in enumerate(out) if SBATCH_RE.match(l))
        out.insert(last + 1, "#SBATCH --account=rrg-shreejoy\n")

    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".", type=Path)
    ap.add_argument("--apply", action="store_true", help="write changes (default is dry-run)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    root = args.root.resolve()
    changed = 0

    for path in iter_files(root):
        rel = str(path.relative_to(root))
        if rel == "scripts/port_scc_paths_to_trillium.py" or rel == "site_env.sh":
            continue

        original = path.read_text(encoding="utf-8", errors="surrogateescape")
        text = port_text(original)

        if path.suffix in SHELL_SUFFIXES:
            text = "".join(port_shell(text.splitlines(keepends=True), rel))

        if text != original:
            changed += 1
            print(f"  {rel}")
            if args.apply:
                path.write_text(text, encoding="utf-8", errors="surrogateescape")

    verb = "rewrote" if args.apply else "would rewrite"
    print(f"\n{verb} {changed} files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
