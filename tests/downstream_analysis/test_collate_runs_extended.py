"""
Tests for extended collate_runs outputs.
"""

import importlib.util
from pathlib import Path

import pandas as pd


_mod_path = Path(__file__).resolve().parent.parent.parent / "bin" / "downstream_analysis" / "collate_runs.py"
_spec = importlib.util.spec_from_file_location("collate_runs", _mod_path)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

collate = _mod.collate


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
    collate([run1.parent, run2.parent], out)

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
