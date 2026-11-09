"""Run Shreejoy's trainer x target CelMod matrix on our manifest, writing only to our scratch.

Engine: /scratch/shreejoy/composition/wt_gtex/stages/composition/retention_matrix.py
(git 134bb7b, which built results_2026-10-21c/retention), imported and run UNCHANGED
except for four module-level overrides made here:

  RM.OUT       -> $STRATEGY_OUT (our scratch); inputs still read from results_2026-10-21c
  RM.declared  -> the same function, source-identical, reading $STRATEGY_MANIFEST
                  instead of $CTCOMP_REPO/manifests/retention_pairs.tsv
  RM.MIN_FIT   -> $MIN_FIT (default 90 = Shreejoy's; 89 for the full run, decision
                  2026-11-04, admits NIMH_HBCC x Catching_2026_HBCC and no other pair)
  RM.Stage     -> a local stub. labcore's Stage deletes <LABCORE_RESULTS>/retention/
                  DONE.json on entry and writes a new one on exit; with Shreejoy's env.sh
                  that is HIS results root. The stub records notes/outputs to our OUT.

CelMod settings, CLR, consensus recipe, hold-out, gene space rule, light/raw bulk are
all Shreejoy's code paths. Nothing is written under /scratch/shreejoy.

Env (set by run_matrix.sbatch): CTCOMP_REPO, REGISTRY_RELEASE, STRATEGY_MANIFEST,
STRATEGY_OUT, MIN_FIT, MAX_TYPES (smoke), LEVELS.
"""
import inspect
import json
import os
import sys
import time
from pathlib import Path

OUT = Path(os.environ["STRATEGY_OUT"]).resolve()
MANIFEST = Path(os.environ["STRATEGY_MANIFEST"]).resolve()
SAFE = Path("/scratch/zhoux156")
if SAFE not in OUT.parents:
    raise SystemExit(f"refusing: STRATEGY_OUT {OUT} is not under {SAFE}")
os.environ["LABCORE_RESULTS"] = str(OUT / "_labcore")  # second guard; Stage is stubbed anyway
os.environ.pop("LABCORE_READS", None)

REPO = Path(os.environ["CTCOMP_REPO"])
sys.path.insert(0, str(REPO / "stages" / "composition"))
import retention_matrix as RM  # noqa: E402


class StubStage:
    """Stands in for labcore.contract.Stage: same calls, writes only to our OUT."""

    def __init__(self, name, description="", depends_on=(), **kw):
        self.name, self.description = name, description
        self.notes, self.outputs, self.t0 = [], {}, time.time()

    def __enter__(self):
        return self

    def note(self, text):
        self.notes.append(text)
        print(f"[note] {text}", flush=True)

    def output(self, path, rows=None):
        p = Path(path).resolve()
        if SAFE not in p.parents:
            raise RuntimeError(f"output outside {SAFE}: {p}")
        self.outputs[str(p)] = rows
        return path

    def log(self, *m):
        print(f"[{self.name}]", *m, flush=True)

    def input(self, path):  # labcore's Stage.input returns the path (rel_by_type reads it)
        return path

    def __getattr__(self, k):  # param/expect/observe etc.: accept, record nothing
        return lambda *a, **kw: a[0] if a else None

    def __exit__(self, et, e, tb):
        rec = dict(stage=self.name, status="ok" if et is None else f"failed: {et.__name__}: {e}",
                   wall_seconds=round(time.time() - self.t0, 1),
                   finished=time.strftime("%Y-%m-%dT%H:%M:%S"),
                   engine=str(REPO / "stages/composition/retention_matrix.py"),
                   manifest=str(MANIFEST), min_fit=RM.MIN_FIT, max_types=RM.MAX_TYPES,
                   level=RM.LEVEL, bulk_src=RM.BULK_SRC, bulk_tgt=RM.BULK_TGT,
                   counts=str(RM.COUNTS), notes=self.notes, outputs=self.outputs)
        (OUT / "RUN.json").write_text(json.dumps(rec, indent=2))
        return False


# declared(): Shreejoy's source with only the manifest path swapped
src = inspect.getsource(RM.declared)
OLD = 'REPO / "manifests" / "retention_pairs.tsv"'
assert src.count(OLD) == 1, "retention_matrix.declared() changed; re-check the override"
ns = dict(vars(RM))
ns["STRATEGY_MANIFEST_PATH"] = MANIFEST
exec(src.replace(OLD, "STRATEGY_MANIFEST_PATH"), ns)

RM.OUT = OUT
RM.declared = ns["declared"]
RM.MIN_FIT = int(os.environ.get("MIN_FIT", RM.MIN_FIT))
RM.Stage = StubStage

for name in ("OUT",):
    assert SAFE in Path(getattr(RM, name)).resolve().parents
print(f"[run_matrix] engine {RM.__file__}\n  manifest {MANIFEST}\n  out {OUT}\n"
      f"  MIN_FIT {RM.MIN_FIT}  MAX_TYPES {RM.MAX_TYPES}  LEVEL {RM.LEVEL}  "
      f"counts {RM.COUNTS}", flush=True)
if os.environ.get("DRY_RUN"):  # login-node check: manifest validation only, no fitting
    m = RM.declared()
    print(m.groupby("role").size().to_string())
    print(f"[dry] Stage is stub: {RM.Stage is StubStage}; OUT {RM.OUT}")
    sys.exit(0)
OUT.mkdir(parents=True, exist_ok=True)
RM.main()
