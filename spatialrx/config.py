"""Path and environment configuration.

All paths default to the repository layout but can be overridden with
environment variables so the pipeline can run against a relocated data tree:

    SPATIALRX_ROOT       repository root (default: two levels above this file)
    SPATIALRX_RAW        raw data dir        (default: <root>/data/raw)
    SPATIALRX_PROCESSED  processed data dir  (default: <root>/data/processed)
    SPATIALRX_OUTPUTS    outputs dir         (default: <root>/outputs)
    SPATIALRX_LOGS       logs dir            (default: <root>/logs)

Nothing in this module imports a scientific package, so it is safe to import at
test-collection time.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


def repo_root() -> Path:
    env = os.environ.get("SPATIALRX_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    return Path(__file__).resolve().parent.parent


@dataclass(frozen=True)
class Paths:
    root: Path
    raw: Path
    extracted: Path
    processed: Path
    outputs: Path
    figures: Path
    tables: Path
    logs: Path

    @classmethod
    def from_env(cls) -> "Paths":
        root = repo_root()
        raw = Path(os.environ.get("SPATIALRX_RAW", root / "data" / "raw")).expanduser()
        processed = Path(
            os.environ.get("SPATIALRX_PROCESSED", root / "data" / "processed")
        ).expanduser()
        outputs = Path(os.environ.get("SPATIALRX_OUTPUTS", root / "outputs")).expanduser()
        logs = Path(os.environ.get("SPATIALRX_LOGS", root / "logs")).expanduser()
        return cls(
            root=root,
            raw=raw,
            extracted=raw / "extracted",
            processed=processed,
            outputs=outputs,
            figures=outputs / "figures",
            tables=outputs / "tables",
            logs=logs,
        )

    def ingest_dir(self) -> Path:
        return self.processed / "00_ingest"

    def domains_dir(self) -> Path:
        return self.processed / "03_domains"

    def ecosystem_dir(self) -> Path:
        return self.processed / "05_ecosystem"


PATHS = Paths.from_env()
