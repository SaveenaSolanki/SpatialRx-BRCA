"""``spatialrx`` command-line entry point.

Implemented and tested commands:
    spatialrx validate-data   data/cohort integrity checks
    spatialrx reproduce       --smoke reproducibility check
    spatialrx env-info        versions + lock/git provenance

Commands are only exposed here once they have tests.  ``reproduce --full`` is
intentionally not implemented yet (Gates 1-3 scope).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .config import PATHS
from .runlog import RunRecord, write_record


def _emit(payload: dict, json_path: str | None) -> None:
    if json_path:
        p = Path(json_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(payload, indent=2, default=str) + "\n")
        print(f"  wrote {p}")


def cmd_validate_data(args) -> int:
    from .validate import print_report, run_validate_data

    rec = RunRecord("validate-data", " ".join(sys.argv), PATHS.root)
    report = run_validate_data(strict=args.strict, do_hash=args.hash)
    print("spatialrx validate-data")
    print_report(report)
    rec.status = "ok" if report["ok"] else "failed"
    rec.result = {
        "ok": report["ok"],
        "n_pass": report["n_pass"],
        "n_fail": report["n_fail"],
        "n_skip": report["n_skip"],
    }
    rec.finished = __import__("time").time()
    if not args.no_log:
        print(f"  log -> {write_record(PATHS.logs / 'remediation', rec)}")
    _emit(report, args.json)
    return 0 if report["ok"] else 1


def cmd_reproduce(args) -> int:
    if args.full:
        print(
            "spatialrx reproduce --full is not implemented in this remediation "
            "branch (Gate 1-3 scope). Use --smoke.",
            file=sys.stderr,
        )
        return 2
    from .smoke import run_smoke

    rec = RunRecord("reproduce-smoke", " ".join(sys.argv), PATHS.root)
    result = run_smoke()
    print("spatialrx reproduce --smoke")
    for name, c in result["checks"].items():
        mark = {"PASS": "✓", "FAIL": "✗", "SKIP": "–"}[c["status"]]
        print(f"  [{mark}] {name:16s} {c['detail']}")
    print(f"\n  smoke: {'PASS' if result['ok'] else 'FAIL'}")
    rec.status = "ok" if result["ok"] else "failed"
    rec.result = {"ok": result["ok"]}
    rec.finished = __import__("time").time()
    if not args.no_log:
        print(f"  log -> {write_record(PATHS.logs / 'remediation', rec)}")
    _emit(result, args.json)
    return 0 if result["ok"] else 1


def cmd_env_info(args) -> int:
    from .runlog import git_info, package_versions

    info = {
        "spatialrx_version": __version__,
        "python": sys.version.split()[0],
        "root": str(PATHS.root),
        "git": git_info(PATHS.root),
        "packages": package_versions(),
    }
    lock = PATHS.root / "requirements.lock"
    if lock.exists():
        info["requirements_lock"] = str(lock)
        info["requirements_lock_lines"] = len(lock.read_text().splitlines())
    print(json.dumps(info, indent=2))
    _emit(info, args.json)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="spatialrx", description=__doc__)
    p.add_argument("--version", action="version", version=f"spatialrx {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    v = sub.add_parser("validate-data", help="validate cohort/data integrity")
    v.add_argument("--strict", action="store_true", help="treat missing data as failure")
    v.add_argument("--hash", action="store_true", help="verify the 1.7 GB raw tar sha256")
    v.add_argument("--json", default=None, help="write the report as JSON")
    v.add_argument("--no-log", action="store_true", help="do not write a run record")
    v.set_defaults(func=cmd_validate_data)

    r = sub.add_parser("reproduce", help="reproducibility checks")
    g = r.add_mutually_exclusive_group(required=True)
    g.add_argument("--smoke", action="store_true", help="fast smoke test")
    g.add_argument("--full", action="store_true", help="full pipeline (not implemented)")
    r.add_argument("--json", default=None)
    r.add_argument("--no-log", action="store_true")
    r.set_defaults(func=cmd_reproduce)

    e = sub.add_parser("env-info", help="print versions and provenance")
    e.add_argument("--json", default=None)
    e.set_defaults(func=cmd_env_info)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
