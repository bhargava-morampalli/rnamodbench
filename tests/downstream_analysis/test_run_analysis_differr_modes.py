"""
Tests for DiffErr score-mode orchestration in run_analysis.py.
"""

import importlib.util
from pathlib import Path
import sys
import types

import pandas as pd


if "sklearn" not in sys.modules:
    sklearn_mod = types.ModuleType("sklearn")
    metrics_mod = types.ModuleType("sklearn.metrics")

    def _dummy_curve(*_args, **_kwargs):
        return [], [], []

    def _dummy_value(*_args, **_kwargs):
        return 0.0

    def _dummy_confusion(*_args, **_kwargs):
        return [[0, 0], [0, 0]]

    def _dummy_prfs(*_args, **_kwargs):
        return [0.0], [0.0], [0.0], [0]

    metrics_mod.precision_recall_curve = _dummy_curve
    metrics_mod.roc_curve = _dummy_curve
    metrics_mod.average_precision_score = _dummy_value
    metrics_mod.roc_auc_score = _dummy_value
    metrics_mod.confusion_matrix = _dummy_confusion
    metrics_mod.precision_recall_fscore_support = _dummy_prfs
    sklearn_mod.metrics = metrics_mod
    sys.modules["sklearn"] = sklearn_mod
    sys.modules["sklearn.metrics"] = metrics_mod


_mod_path = Path(__file__).resolve().parent.parent.parent / "bin" / "downstream_analysis" / "run_analysis.py"
_spec = importlib.util.spec_from_file_location("run_analysis", _mod_path)
_mod = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _mod
_spec.loader.exec_module(_mod)


def test_resolve_differr_runs_both_mode(tmp_path):
    output_dir = tmp_path / "downstream_analysis"
    planned = _mod._resolve_differr_runs(output_dir, "both")
    assert planned == [
        ("g_fdr_neglog10", output_dir),
        ("g_stat", output_dir / "differr_gstat"),
    ]


def test_main_dispatches_dual_mode_runs(tmp_path, monkeypatch):
    calls = []

    def _fake_load(_args, differr_score_field):
        assert differr_score_field in {"g_fdr_neglog10", "g_stat"}
        return {"tombo": pd.DataFrame({"score": [1.0]})}

    def _fake_run_analysis(**kwargs):
        calls.append((kwargs["differr_score_field"], Path(kwargs["output_dir"])))

    monkeypatch.setattr(_mod, "_load_tool_outputs_for_mode", _fake_load)
    monkeypatch.setattr(_mod, "run_analysis", _fake_run_analysis)

    out = tmp_path / "downstream_analysis"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_analysis.py",
            "--input-dir",
            str(tmp_path),
            "--output-dir",
            str(out),
            "--differr-score-field",
            "both",
        ],
    )

    _mod.main()

    assert calls == [
        ("g_fdr_neglog10", out),
        ("g_stat", out / "differr_gstat"),
    ]
