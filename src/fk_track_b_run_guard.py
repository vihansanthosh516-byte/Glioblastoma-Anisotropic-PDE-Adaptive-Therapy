"""Pre-run seal and real-data entry guard for the Track B LUMIERE protocol (draft; not adopted).

A seal is a JSON file written BEFORE any real-data run. It records the PDE provider (class and per-scenario
settings), the full FiLM/U-Net recipe, the bootstrap count, a fold-manifest hash, an eligibility/QC manifest hash
(patient id, qc_passed, interval weeks, spacing; no voxel data), SHA-256 of every execution code file and of the
protocol documents, git HEAD with a clean-tree flag, and python/numpy/scipy/torch versions.

``run_sealed_protocol`` is the only real-data entry. Real (synthetic=False) cases are refused by
``run_lumiere_cv`` and ``run_lumiere_protocol`` unless they carry the private ``_SEALED_TOKEN`` that only this
function passes. It re-derives everything from the current repository and refuses (``SealMismatch``) before any
case target is read if the seal is missing, its SHA-256 differs from the caller's ``expected_seal_sha256``, the git
state is unknown or dirty, versions differ, or anything else differs. It reads no patient data itself.

The injectable synthetic path (``fk_lumiere_protocol.run_lumiere_protocol`` with synthetic cases) is unchanged.
Limits: the seal does not hash patient input files; a caller who edits this module could still bypass it (the
token blocks accidents and casual misuse, not a determined editor; the code hash and git check expose the edit).
"""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

from src import fk_lumiere_protocol as lp
from src import fk_lumiere_result_writer as rw
from src import fk_lumiere_runner as lr
from src.fk_learned_baselines import UNET_RECIPE
from src.fk_lumiere_stats import N_BOOTSTRAP

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCHEMA_VERSION = 2
# Private token: only run_sealed_protocol passes it down; the runner and protocol refuse real cases without it.
_SEALED_TOKEN = object()


class SealMismatch(PermissionError):
    """The seal is missing, malformed, or differs from the settings about to be used."""


def _jsonable(x):
    if isinstance(x, Mapping):
        return {str(k): _jsonable(v) for k, v in sorted(x.items(), key=lambda kv: str(kv[0]))}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (str, int, float, bool)) or x is None:
        return x
    if callable(x):
        return f"{getattr(x, '__module__', '?')}.{getattr(x, '__qualname__', repr(x))}"
    return repr(x)


def _check_n_bootstrap(n_bootstrap) -> int:
    n = int(N_BOOTSTRAP if n_bootstrap is None else n_bootstrap)
    if n < 1:
        raise ValueError(f"n_bootstrap must be >= 1, got {n}")
    return n


def _git_probe(project_root, files) -> dict:
    """HEAD and clean-tree flag for the given files. 'unknown' if git is missing or fails."""
    def run(*args):
        return subprocess.run(["git", "-C", str(project_root), *args], capture_output=True, text=True,
                              timeout=60, check=True).stdout
    try:
        head = run("rev-parse", "HEAD").strip()
        dirty = run("status", "--porcelain", "--", *files).strip()
    except (OSError, subprocess.SubprocessError):
        return {"head": "unknown", "clean": "unknown"}
    return {"head": head or "unknown", "clean": not dirty}


def _env_versions() -> dict:
    out = {"python": platform.python_version()}
    for name in ("numpy", "scipy", "torch"):
        try:
            out[name] = __import__(name).__version__
        except Exception:  # not installed or broken: record it and compare it
            out[name] = "unavailable"
    return out


def eligibility_manifest_sha256(cases: Sequence[lr.LumiereCase]) -> str:
    """Hash of patient id and the non-voxel QC/interval fields per case. No image or target data is read."""
    lines = [f"{c.patient_id}|qc={bool(c.qc_passed)}|fw={c.forecast_weeks}|hw={c.history_weeks}|"
             f"sp={tuple(float(x) for x in c.spacing_mm)}" for c in sorted(cases, key=lambda c: c.patient_id)]
    return hashlib.sha256("\n".join(lines).encode()).hexdigest()


def describe_provider(provider) -> dict:
    """Class identity plus every public configuration attribute (not the call log)."""
    cls = type(provider)
    attrs = {k: v for k, v in vars(provider).items() if not k.startswith("_") and k != "calls"}
    return _jsonable({"class": f"{cls.__module__}.{cls.__qualname__}", "attributes": attrs})


def describe_provider_factory(provider_factory: Callable) -> dict:
    return {name: describe_provider(provider_factory(rho_b, d_b))
            for name, (rho_b, d_b) in lp.PRIOR_SCENARIOS.items()}


def fold_manifest_sha256(fold_by_patient: Mapping[str, int], patient_ids: Sequence[str]) -> str:
    lines = [f"{p}:{int(fold_by_patient[p])}" for p in sorted(patient_ids)]
    return hashlib.sha256("\n".join(lines).encode()).hexdigest()


def _current_state(*, provider_factory, recipe, n_bootstrap, fold_by_patient, patient_ids, project_root,
                   eligibility_sha256) -> dict:
    n_boot = _check_n_bootstrap(n_bootstrap)
    files = rw.execution_code_files(project_root)
    code = rw.file_hashes(files, project_root)
    protocol = rw.file_hashes(rw.PROTOCOL_FILES, project_root)
    return {
        "git": _git_probe(project_root, list(files) + list(rw.PROTOCOL_FILES)),
        "versions": _env_versions(),
        "eligibility_manifest_sha256": eligibility_sha256,
        "provider": describe_provider_factory(provider_factory),
        "recipe": _jsonable(asdict(recipe)),
        "n_bootstrap": n_boot,
        "fold_manifest_sha256": fold_manifest_sha256(fold_by_patient, patient_ids),
        "patient_ids": sorted(patient_ids),
        "code_sha256": code,
        "protocol_sha256": protocol,
    }


def build_seal(*, cases: Sequence[lr.LumiereCase], fold_by_patient: Mapping[str, int],
               provider_factory: Callable = lp.default_provider_factory, recipe=UNET_RECIPE,
               n_bootstrap: int | None = None, project_root=PROJECT_ROOT) -> dict:
    state = _current_state(provider_factory=provider_factory, recipe=recipe, n_bootstrap=n_bootstrap,
                           fold_by_patient=fold_by_patient, patient_ids=[c.patient_id for c in cases],
                           project_root=project_root, eligibility_sha256=eligibility_manifest_sha256(cases))
    return {"schema_version": SCHEMA_VERSION, "created_utc": datetime.now(timezone.utc).isoformat(),
            "record_type": "pre-run seal (code, protocol, settings, folds, eligibility manifest, git, versions); "
                           "patient input files NOT hashed",
            **state}


def write_seal(seal: Mapping, path) -> Path:
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite seal {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(seal, indent=1, sort_keys=True))
    return path


def file_sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_seal(path) -> dict:
    path = Path(path) if path is not None else None
    if path is None or not path.is_file():
        raise SealMismatch(f"no seal at {path}; a real-data run requires a pre-run seal")
    try:
        seal = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise SealMismatch(f"seal is not valid JSON: {exc}") from exc
    needed = ("provider", "recipe", "n_bootstrap", "fold_manifest_sha256", "patient_ids", "code_sha256",
              "protocol_sha256", "git", "versions", "eligibility_manifest_sha256")
    absent = [k for k in needed if k not in seal]
    if absent:
        raise SealMismatch(f"seal is missing fields: {absent}")
    return seal


def _diff_hashes(label, sealed: Mapping, now: Mapping) -> list[str]:
    out = [f"{label}: file changed {f}" for f in sorted(set(sealed) & set(now)) if sealed[f] != now[f]]
    out += [f"{label}: file not in seal {f}" for f in sorted(set(now) - set(sealed))]
    out += [f"{label}: sealed file now absent {f}" for f in sorted(set(sealed) - set(now))]
    return out


def verify_seal(seal: Mapping, current: Mapping) -> list[str]:
    problems = []
    if seal["provider"] != current["provider"]:
        problems.append("provider class or settings differ from the seal")
    if seal["recipe"] != current["recipe"]:
        problems.append("model recipe differs from the seal")
    if seal["n_bootstrap"] != current["n_bootstrap"]:
        problems.append(f"bootstrap count {current['n_bootstrap']} differs from sealed {seal['n_bootstrap']}")
    if seal["patient_ids"] != current["patient_ids"]:
        problems.append("patient set differs from the seal")
    if seal["fold_manifest_sha256"] != current["fold_manifest_sha256"]:
        problems.append("fold manifest differs from the seal")
    if seal["eligibility_manifest_sha256"] != current["eligibility_manifest_sha256"]:
        problems.append("eligibility/QC manifest differs from the seal")
    if seal["versions"] != current["versions"]:
        problems.append(f"package versions differ from the seal: sealed {seal['versions']}, now {current['versions']}")
    for label, git in (("seal", seal["git"]), ("current tree", current["git"])):
        if git.get("head") in (None, "unknown") or git.get("clean") is not True:
            problems.append(f"{label}: git HEAD unknown or sealed files not clean ({git})")
    problems += _diff_hashes("code", seal["code_sha256"], current["code_sha256"])
    problems += _diff_hashes("protocol", seal["protocol_sha256"], current["protocol_sha256"])
    return problems


def run_sealed_protocol(cases: Sequence[lr.LumiereCase], fold_by_patient: Mapping[str, int], *, seal_path,
                        expected_seal_sha256: str, gates: Mapping[str, bool] | None = None,
                        provider_factory: Callable | None = None, recipe=None, n_bootstrap: int | None = None,
                        project_root=PROJECT_ROOT, allow_synthetic: bool = False) -> dict:
    """Real-data entry. Verifies the seal first; no case target is read before it passes.

    Requires every case synthetic=False (``allow_synthetic=True`` exists only to test the guard itself), every
    LUMIERE gate True, the seal file's SHA-256 to equal ``expected_seal_sha256`` (recorded out of band), a clean
    git tree with known HEAD, and matching package versions. The optional overrides exist only so a caller can
    prove they equal the sealed values; any difference raises. There is no ``minimum`` argument.
    """
    if not cases:
        raise ValueError("no cases")
    if not allow_synthetic and any(c.synthetic for c in cases):
        raise PermissionError("refused: the real entry takes only synthetic=False cases")
    if not all(c.synthetic for c in cases):
        if any(c.synthetic for c in cases):
            raise PermissionError("refused: synthetic and real cases cannot be mixed")
        missing = [x for x in lr.LUMIERE_GATES if (gates or {}).get(x) is not True]
        if missing:
            raise PermissionError(f"scoring refused: gates not recorded as passed: {missing}")
    seal = load_seal(seal_path)
    actual = file_sha256(seal_path)
    if actual != expected_seal_sha256:
        raise SealMismatch(f"seal file SHA-256 {actual} differs from expected {expected_seal_sha256}")
    factory = provider_factory if provider_factory is not None else lp.default_provider_factory
    effective_recipe = recipe if recipe is not None else UNET_RECIPE
    try:
        current = _current_state(provider_factory=factory, recipe=effective_recipe, n_bootstrap=n_bootstrap,
                                 fold_by_patient=fold_by_patient, patient_ids=[c.patient_id for c in cases],
                                 project_root=project_root, eligibility_sha256=eligibility_manifest_sha256(cases))
    except ValueError as exc:
        raise SealMismatch(str(exc)) from exc
    except KeyError as exc:
        raise SealMismatch(f"fold manifest has no entry for patient {exc}") from exc
    except FileNotFoundError as exc:
        raise SealMismatch(str(exc)) from exc
    problems = verify_seal(seal, current)
    if problems:
        raise SealMismatch("seal check failed before any target read: " + "; ".join(problems))
    return lp.run_lumiere_protocol(cases, fold_by_patient, provider_factory=factory, gates=gates,
                                   recipe=effective_recipe, n_bootstrap=n_bootstrap, _sealed_token=_SEALED_TOKEN)
