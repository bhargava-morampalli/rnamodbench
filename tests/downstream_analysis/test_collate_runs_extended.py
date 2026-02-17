"""
Tests for extended collate_runs outputs.
"""

import importlib.util
import json
from pathlib import Path
import sys

import pandas as pd


_mod_path = Path(__file__).resolve().parent.parent.parent / "bin" / "downstream_analysis" / "collate_runs.py"
_spec = importlib.util.spec_from_file_location("collate_runs", _mod_path)
_mod = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _mod
_spec.loader.exec_module(_mod)

collate = _mod.collate


def _write_csv(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(data).to_csv(path, index=False)


def test_collate_writes_all_expected_tables(tmp_path):
    run1 = tmp_path / "run1" / "collation"
    run2 = tmp_path / "run2" / "collation"
    run1.mkdir(parents=True)
    run2.mkdir(parents=True)

    pd.DataFrame({"run_id": ["r1"], "tool": ["tombo"]}).to_csv(run1 / "metrics_long.csv", index=False)
    pd.DataFrame({"run_id": ["r2"], "tool": ["xpore"]}).to_csv(run2 / "metrics_long.csv", index=False)

    pd.DataFrame({"run_id": ["r1"], "tool": ["tombo"], "window_size": [0]}).to_csv(
        run1 / "window_metrics_long.csv", index=False
    )
    pd.DataFrame({"run_id": ["r2"], "tool": ["xpore"], "delta": [1]}).to_csv(
        run2 / "lag_metrics_long.csv", index=False
    )

    out = tmp_path / "out"
    report = collate([run1.parent, run2.parent], out, quiet=True)

    expected = [
        "metrics_long.csv",
        "metrics_summary_long.csv",
        "window_metrics_long.csv",
        "window_metrics_summary_long.csv",
        "lag_metrics_long.csv",
        "lag_metrics_summary_long.csv",
    ]
    for name in expected:
        assert (out / name).exists(), f"missing {name}"

    merged = pd.read_csv(out / "metrics_long.csv")
    assert len(merged) == 2
    assert set(merged["tool"]) == {"tombo", "xpore"}
    assert report.runs_usable == 2


def test_empty_file_is_reported_not_crash(tmp_path):
    run1 = tmp_path / "run1" / "collation"
    run2 = tmp_path / "run2" / "collation"
    run1.mkdir(parents=True)
    run2.mkdir(parents=True)

    _write_csv(run1 / "metrics_long.csv", {"run_id": ["r1"], "tool": ["tombo"]})
    (run2 / "metrics_long.csv").write_text("", encoding="utf-8")

    out = tmp_path / "out"
    report = collate([run1.parent, run2.parent], out, tables=["metrics_long.csv"], quiet=True)

    merged = pd.read_csv(out / "metrics_long.csv")
    assert len(merged) == 1
    assert any(issue.code == "empty_file" for issue in report.issues)


def test_missing_table_is_warned_and_continues(tmp_path):
    run1 = tmp_path / "run1" / "collation"
    run2 = tmp_path / "run2" / "collation"
    run1.mkdir(parents=True)
    run2.mkdir(parents=True)

    _write_csv(run1 / "metrics_long.csv", {"run_id": ["r1"], "tool": ["tombo"]})

    out = tmp_path / "out"
    report = collate([run1.parent, run2.parent], out, tables=["metrics_long.csv"], quiet=True)

    merged = pd.read_csv(out / "metrics_long.csv")
    assert len(merged) == 1
    assert any(issue.code == "missing_table_file" for issue in report.issues)


def test_schema_drift_union_columns_and_reported(tmp_path):
    run1 = tmp_path / "run1" / "collation"
    run2 = tmp_path / "run2" / "collation"
    run1.mkdir(parents=True)
    run2.mkdir(parents=True)

    _write_csv(run1 / "metrics_long.csv", {"run_id": ["r1"], "tool": ["tombo"]})
    _write_csv(run2 / "metrics_long.csv", {"run_id": ["r2"], "tool": ["xpore"], "extra_col": [42]})

    out = tmp_path / "out"
    report = collate([run1.parent, run2.parent], out, tables=["metrics_long.csv"], quiet=True)

    merged = pd.read_csv(out / "metrics_long.csv")
    assert "extra_col" in merged.columns
    table_report = report.table_records[0]
    assert table_report.schema_drift is True
    assert any(issue.code == "schema_drift" for issue in report.issues)


def test_nonexistent_run_path_warns(tmp_path):
    run1 = tmp_path / "run1" / "collation"
    run1.mkdir(parents=True)
    _write_csv(run1 / "metrics_long.csv", {"run_id": ["r1"], "tool": ["tombo"]})

    out = tmp_path / "out"
    bad = tmp_path / "does_not_exist"
    report = collate([run1.parent, bad], out, tables=["metrics_long.csv"], quiet=True)

    merged = pd.read_csv(out / "metrics_long.csv")
    assert len(merged) == 1
    assert any(issue.code == "missing_input_path" for issue in report.issues)


def test_no_zero_byte_outputs_written(tmp_path):
    run1 = tmp_path / "run1" / "collation"
    run1.mkdir(parents=True)

    out = tmp_path / "out"
    collate([run1.parent], out, tables=["metrics_long.csv"], quiet=True)

    metrics = out / "metrics_long.csv"
    assert metrics.exists()
    assert metrics.stat().st_size > 0
    df = pd.read_csv(metrics)
    assert df.empty


def test_markdown_and_json_reports_written(tmp_path):
    run1 = tmp_path / "run1" / "collation"
    run1.mkdir(parents=True)
    _write_csv(run1 / "metrics_long.csv", {"run_id": ["r1"], "tool": ["tombo"]})

    out = tmp_path / "out"
    report = collate([run1.parent], out, tables=["metrics_long.csv"], quiet=True)

    report_json = out / "collate_report.json"
    report_md = out / "collate_report.md"
    assert report_json.exists()
    assert report_md.exists()

    payload = json.loads(report_json.read_text(encoding="utf-8"))
    assert payload["runs_usable"] == 1
    assert report.output_files["report_json"] == str(report_json)
    assert report.output_files["report_markdown"] == str(report_md)


def test_discovery_mode_finds_expected_runs(tmp_path):
    root = tmp_path / "covbench_results"
    run1 = root / "results_5x" / "downstream_analysis" / "collation"
    run2 = root / "results_10x" / "downstream_analysis" / "collation"
    run1.mkdir(parents=True)
    run2.mkdir(parents=True)

    _write_csv(run1 / "metrics_long.csv", {"run_id": ["r1"], "tool": ["tombo"]})
    _write_csv(run2 / "metrics_long.csv", {"run_id": ["r2"], "tool": ["xpore"]})

    out = tmp_path / "out"
    report = collate(
        [],
        out,
        tables=["metrics_long.csv"],
        runs_root=root,
        run_glob="results_*/downstream_analysis",
        quiet=True,
    )

    merged = pd.read_csv(out / "metrics_long.csv")
    assert len(merged) == 2
    assert report.runs_usable == 2
    assert report.runs_discovered == 2


def test_strict_mode_returns_nonzero_on_soft_issues(tmp_path):
    run1 = tmp_path / "run1" / "collation"
    run1.mkdir(parents=True)
    # No metrics file -> soft issue (missing table), still usable run.

    out = tmp_path / "out"
    rc = _mod.main(
        [
            "--runs",
            str(run1.parent),
            "--output",
            str(out),
            "--tables",
            "metrics_long.csv",
            "--strict",
            "--quiet",
        ]
    )
    assert rc == 1
