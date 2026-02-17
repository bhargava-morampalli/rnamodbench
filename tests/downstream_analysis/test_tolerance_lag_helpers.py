"""
Unit tests for tolerance/lag helper logic in run_analysis.py.
"""

import importlib.util
from pathlib import Path

import pandas as pd


_mod_path = Path(__file__).resolve().parent.parent.parent / "bin" / "downstream_analysis" / "run_analysis.py"
_spec = importlib.util.spec_from_file_location("run_analysis", _mod_path)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

_expand_positions_by_window = _mod._expand_positions_by_window
_shift_positions = _mod._shift_positions
_compute_labeled_scope_metrics = _mod._compute_labeled_scope_metrics


def test_expand_positions_window_respects_bounds_and_monotonicity():
    gt = {1, 5}
    ref_len = 5

    w1 = _expand_positions_by_window(gt, 1, ref_len)
    w2 = _expand_positions_by_window(gt, 2, ref_len)

    assert w1 == {1, 2, 4, 5}
    assert w2 == {1, 2, 3, 4, 5}
    assert w1.issubset(w2)


def test_shift_positions_respects_bounds():
    gt = {1, 3, 5}
    ref_len = 5

    plus_one = _shift_positions(gt, 1, ref_len)
    minus_two = _shift_positions(gt, -2, ref_len)

    assert plus_one == {2, 4}
    assert minus_two == {1, 3}


def test_compute_labeled_scope_metrics_single_class_returns_nan_metrics():
    df = pd.DataFrame(
        {
            "position": [1, 2, 3],
            "score_metric": [0.1, 0.2, 0.3],
        }
    )

    row = _compute_labeled_scope_metrics(df, positive_positions=set())
    assert row["metric_scope_note"] == "single_class_labels"
    assert pd.isna(row["auroc"])
    assert pd.isna(row["auprc"])
    assert row["n_positive"] == 0
    assert row["n_universe"] == 3


def test_compute_labeled_scope_metrics_binary_case_is_computable():
    df = pd.DataFrame(
        {
            "position": [1, 2, 3, 4],
            "score_metric": [0.1, 0.9, 0.2, 0.8],
        }
    )

    row = _compute_labeled_scope_metrics(df, positive_positions={2, 4})
    assert row["metric_scope_note"] == "ok"
    assert row["n_positive"] == 2
    assert row["n_universe"] == 4
    assert 0.0 <= row["auroc"] <= 1.0
    assert 0.0 <= row["auprc"] <= 1.0
