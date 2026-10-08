"""Tests for the spatialrx CLI surface (only expose tested commands)."""

from __future__ import annotations

import json

from spatialrx.cli import build_parser, main


def test_parser_has_tested_commands():
    parser = build_parser()
    for cmd in ("validate-data", "reproduce", "env-info"):
        assert parser.parse_args([cmd] + (["--smoke"] if cmd == "reproduce" else []))


def test_env_info_writes_json(tmp_path):
    out = tmp_path / "env.json"
    rc = main(["env-info", "--json", str(out)])
    assert rc == 0
    payload = json.loads(out.read_text())
    assert "packages" in payload and "git" in payload
    assert payload["packages"]["scanpy"] is not None


def test_reproduce_full_is_explicitly_unimplemented():
    assert main(["reproduce", "--full", "--no-log"]) == 2


def test_validate_data_writes_json_without_logging(tmp_path):
    out = tmp_path / "validate.json"
    rc = main(["validate-data", "--no-log", "--json", str(out)])
    assert rc in (0, 1)
    payload = json.loads(out.read_text())
    assert "checks" in payload and "ok" in payload
