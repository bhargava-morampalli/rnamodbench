"""
Regression test for invalid-score gate behavior in run_analysis.py.

When invalid-score fraction is high for pvalue/fdr tools, QC flags should remain
set but computable metrics/curves should no longer be suppressed.
"""

import importlib.util
from pathlib import Path
import sys
import types

import numpy as np
import pandas as pd


try:
    import sklearn.metrics  # noqa: F401
except Exception:
    sklearn_mod = types.ModuleType("sklearn")
    metrics_mod = types.ModuleType("sklearn.metrics")

    def _dummy_pr_curve(*_args, **_kwargs):
        return [1.0, 0.5], [0.0, 1.0], [0.5]

    def _dummy_roc_curve(*_args, **_kwargs):
        return [0.0, 1.0], [0.0, 1.0], [1.0, 0.0]

    def _dummy_value(*_args, **_kwargs):
        return 0.0

    def _dummy_confusion(*_args, **_kwargs):
        return [[0, 0], [0, 0]]

    def _dummy_prfs(*_args, **_kwargs):
        return [0.0], [0.0], [0.0], [0]

    metrics_mod.precision_recall_curve = _dummy_pr_curve
    metrics_mod.roc_curve = _dummy_roc_curve
    metrics_mod.average_precision_score = _dummy_value
    metrics_mod.roc_auc_score = _dummy_value
    metrics_mod.confusion_matrix = _dummy_confusion
    metrics_mod.precision_recall_fscore_support = _dummy_prfs
    sklearn_mod.metrics = metrics_mod
    sys.modules["sklearn"] = sklearn_mod
    sys.modules["sklearn.metrics"] = metrics_mod


_mod_path = (
    Path(__file__).resolve().parent.parent.parent
    / "bin"
    / "downstream_analysis"
    / "run_analysis.py"
)
_spec = importlib.util.spec_from_file_location("run_analysis", _mod_path)
_mod = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _mod
_spec.loader.exec_module(_mod)


def _is_true(value) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def test_invalid_score_gate_keeps_computable_outputs(tmp_path, monkeypatch):
    # Keep this test focused on metrics behavior, not figure generation.
    monkeypatch.setattr(_mod, "plot_roc_tool_grid", lambda **_kwargs: None)
    monkeypatch.setattr(_mod, "plot_pr_tool_grid", lambda **_kwargs: None)
    monkeypatch.setattr(_mod, "plot_window_tolerance_panels", lambda **_kwargs: None)
    monkeypatch.setattr(_mod, "plot_lag_scan_panels", lambda **_kwargs: None)

    tool_df = pd.DataFrame(
        {
            "tool": ["nanocompore"] * 4,
            "reference": ["refA"] * 4,
            "position": [1, 2, 3, 4],
            "score": [np.nan, 1.2, -0.1, 0.4],  # 3/4 invalid raw pvalue scores
            "score_type": ["pvalue"] * 4,
            "pvalue": [np.nan, 1.2, -0.1, 0.4],
            "replicate": ["rep1"] * 4,
            "sample_id": ["sampleA"] * 4,
            "coverage": ["5x"] * 4,
            "_imputed": [False] * 4,
        }
    )
    tool_outputs = {"nanocompore": tool_df}
    ground_truth = pd.DataFrame({"reference": ["refA", "refA"], "position": [2, 4]})

    output_dir = tmp_path / "out_invalid_gate"
    _mod.run_analysis(
        tool_outputs=tool_outputs,
        ground_truth=ground_truth,
        output_dir=output_dir,
        references_csv=None,
        references_filter=None,
        min_replicates=2,
        score_threshold=None,
        expected_replicates=1,
        run_id="test_gate",
        coverage_label="5x",
        quality_label="unit",
        differr_score_field="g_fdr_neglog10",
    )

    metrics_path = output_dir / "by_reference" / "refA" / "04_metrics" / "metrics_per_replicate.csv"
    assert metrics_path.exists()
    metrics = pd.read_csv(metrics_path)
    row = metrics.loc[
        (metrics["tool"] == "nanocompore") & (metrics["replicate"].astype(str) == "rep1")
    ].iloc[0]

    assert not _is_true(row["metrics_valid"])
    assert row["metric_scope_note"] == "invalid_scores"
    assert int(row["invalid_score_n"]) == 3
    assert int(row["invalid_score_total"]) == 4
    assert float(row["invalid_score_fraction"]) > 0.5

    # Core assertion: metrics remain populated when computable, despite gate.
    for col in ["auprc", "auroc", "f1_optimal", "auprc_reported", "auroc_reported", "f1_reported"]:
        assert pd.notna(row[col]), f"{col} should remain computed for gate-triggered replicate"

    pr_path = output_dir / "by_reference" / "refA" / "05_curves" / "pr_curve_points.csv"
    roc_path = output_dir / "by_reference" / "refA" / "05_curves" / "roc_curve_points.csv"
    assert pr_path.exists()
    assert roc_path.exists()
    pr_df = pd.read_csv(pr_path)
    roc_df = pd.read_csv(roc_path)
    assert not pr_df.empty
    assert not roc_df.empty
    assert ((pr_df["tool"] == "nanocompore") & (pr_df["replicate"].astype(str) == "rep1")).any()
    assert ((roc_df["tool"] == "nanocompore") & (roc_df["replicate"].astype(str) == "rep1")).any()

    window_path = (
        output_dir
        / "by_reference"
        / "refA"
        / "09_tolerance_analysis"
        / "window_metrics_per_replicate.csv"
    )
    lag_path = (
        output_dir
        / "by_reference"
        / "refA"
        / "10_lag_diagnostics"
        / "lag_metrics_per_replicate.csv"
    )
    window_df = pd.read_csv(window_path)
    lag_df = pd.read_csv(lag_path)

    window_sub = window_df[
        (window_df["tool"] == "nanocompore") & (window_df["replicate"].astype(str) == "rep1")
    ].copy()
    lag_sub = lag_df[
        (lag_df["tool"] == "nanocompore") & (lag_df["replicate"].astype(str) == "rep1")
    ].copy()

    assert not window_sub.empty
    assert not lag_sub.empty
    assert (window_sub["metric_scope_note"] == "invalid_scores").any()
    assert (lag_sub["metric_scope_note"] == "invalid_scores").any()
    assert window_sub["auroc"].notna().any()
    assert lag_sub["auroc"].notna().any()
