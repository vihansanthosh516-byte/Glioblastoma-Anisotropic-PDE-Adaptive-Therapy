"""Small, data-safe provenance record for Track B exploratory runs."""
from __future__ import annotations

import hashlib
import platform
import subprocess
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Mapping


DEFAULT_DISTRIBUTIONS = {
    "numpy": "numpy",
    "scipy": "scipy",
    "scikit-learn": "scikit-learn",
    "nibabel": "nibabel",
}


def sha256_file(path: Path) -> str:
    """Hash one explicitly named metadata or source file in bounded memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_state(repo_root: Path) -> tuple[str | None, bool | None]:
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=repo_root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        return None, None
    return revision, bool(status.strip())


def _distribution_versions(distributions: Mapping[str, str]) -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for display_name, distribution_name in distributions.items():
        try:
            versions[display_name] = version(distribution_name)
        except PackageNotFoundError:
            versions[display_name] = None
    return versions


def build_run_provenance(
    repo_root: Path,
    *,
    manifests: Mapping[str, Path],
    source_files: Mapping[str, Path],
    distributions: Mapping[str, str] = DEFAULT_DISTRIBUTIONS,
) -> dict[str, object]:
    """Record revision, named file hashes, and runtime versions.

    This deliberately does not open MRI, segmentation, or other cohort files.
    """
    root = Path(repo_root).resolve()
    revision, dirty = _git_state(root)
    manifest_hashes = {
        name: sha256_file(path)
        for name, path in sorted(manifests.items())
    }
    source_hashes = {
        name: sha256_file(path)
        for name, path in sorted(source_files.items())
    }
    return {
        "schema_version": 1,
        "source_revision": revision,
        "working_tree_dirty": dirty,
        "manifest_sha256": manifest_hashes,
        "source_file_sha256": source_hashes,
        "python_version": platform.python_version(),
        "package_versions": _distribution_versions(distributions),
    }
