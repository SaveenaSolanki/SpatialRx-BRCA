"""Structured run logging.

Every CLI invocation records the git revision, package versions and outcome to
a JSON file under ``logs/remediation/`` so results have a machine-readable
provenance trail (the original pipeline only had free-text logs).
"""

from __future__ import annotations

import json
import platform
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

KEY_PACKAGES = [
    "numpy", "scipy", "pandas", "anndata", "scanpy", "squidpy", "decoupler",
    "scikit-learn", "statsmodels", "matplotlib", "seaborn", "harmonypy",
    "leidenalg", "igraph", "numba", "h5py", "torch", "zarr",
]


def _git(args: list[str], cwd: Path) -> str | None:
    try:
        out = subprocess.run(
            ["git", *args], cwd=str(cwd), capture_output=True, text=True, timeout=10
        )
        if out.returncode == 0:
            return out.stdout.strip()
    except Exception:
        pass
    return None


def git_info(root: Path) -> dict:
    return {
        "commit": _git(["rev-parse", "HEAD"], root),
        "short": _git(["rev-parse", "--short", "HEAD"], root),
        "branch": _git(["rev-parse", "--abbrev-ref", "HEAD"], root),
        "dirty": bool(_git(["status", "--porcelain"], root)),
    }


def package_versions() -> dict:
    from importlib.metadata import PackageNotFoundError, version

    out = {}
    for name in KEY_PACKAGES:
        try:
            out[name] = version(name)
        except PackageNotFoundError:
            out[name] = None
        except Exception:
            out[name] = "unknown"
    return out


@dataclass
class RunRecord:
    name: str
    command: str
    root: Path
    status: str = "running"
    result: dict = field(default_factory=dict)
    started: float = field(default_factory=time.time)
    finished: float | None = None

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "command": self.command,
            "status": self.status,
            "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(self.started)),
            "finished_utc": (
                time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(self.finished))
                if self.finished
                else None
            ),
            "duration_s": round((self.finished or time.time()) - self.started, 3),
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "git": git_info(self.root),
            "packages": package_versions(),
            "result": self.result,
        }


def write_record(log_dir: Path, record: RunRecord) -> Path:
    log_dir.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime(record.started))
    path = log_dir / f"{ts}_{record.name}.json"
    path.write_text(json.dumps(record.as_dict(), indent=2) + "\n")
    return path
