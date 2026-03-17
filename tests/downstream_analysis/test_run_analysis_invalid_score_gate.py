"""
Regression test for invalid-score gate behavior in run_analysis.py.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


_mod_path = Path(__file__).resolve().parent.parent.parent / "bin" / "downstream_analysis" / "run_analysis.py"
_spec = importlib.util.spec_from_file_location("run_analysis", _mod_path)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)


class _DummyQualityReport:
    def __init__(self) -> None:
        self.warnings = []
        self.tool_statuses = []

    def save(self, path: Path) -> None:
        Path(path).write_text("ok\n", encoding="utf-8")

    def save_warnings(self, path: Path) -> None:
        Path(path).write_text("", encoding="utf-8")

    def summary(self) -> str:
        return "ok"


class _DummyMetrics:
    def __init__(self, tool: str, reference: str, replicate: str, n_predictions: int) -> None:
        self.tool = tool
        self.reference = reference
        self.replicate = replicate
        self.n_predictions = n_predictions

    def to_dict(self) -> dict:
        return {
            "tool": self.tool,
            "reference": self.reference,
            "replicate": self.replicate,
            "auprc": 0.61,
            "auroc": 0.57,
            "f1_optimal": 0.54,
            "precision_optimal": 0.55,
            "recall_optimal": 0.53,
            "optimal_threshold": 0.5,
            "n_true_positives": 2,
            "n_false_positives": 1,
            "n_false_negatives": 1,
            "n_true_negatives": 2,
            "n_ground_truth": 2,
            "n_predictions": self.n_predictions,
        }


class _DummyComparison:
    n_tools = 1
    n_union = 6
    n_intersection = 6
    n_positions_per_tool = {"nanocompore": 6}


def test_invalid_score_gate_keeps_metrics_and_curve_outputs(monkeypatch, tmp_path):
    tool_df = pd.DataFrame(
        {
            "tool": ["nanocompore"] * 6,
            "reference": ["r1"] * 6,
            "position": [1, 2, 3, 4, 5, 6],
            "replicate": ["rep1"] * 6,
            "score": [0.9, 0.1, 0.8, 0.2, 0.7, 0.3],
            "score_type": ["pvalue"] * 6,
        }
    )

    metric_ready = pd.DataFrame(
        {
            "tool": ["nanocompore"] * 6,
            "reference": ["r1"] * 6,
            "position": [1, 2, 3, 4, 5, 6],
            "replicate": ["rep1"] * 6,
            "score_type": ["pvalue"] * 6,
            "score_raw": [np.nan, np.inf, 1.2, 0.01, np.nan, np.nan],
            "score_metric": [0.1, 0.9, 0.2, 0.8, 0.3, 0.7],
            "is_reported": [True, True, True, True, False, False],
            "label": [0, 1, 0, 1, 0, 0],
        }
    )

    ground_truth = pd.DataFrame({"reference": ["r1", "r1"], "position": [2, 4]})

    def _fake_standardize(*_args, **_kwargs):
        return tool_df.copy()

    def _fake_build_metric_ready(_imputed, _gt_positions):
        return metric_ready.copy()

    def _fake_no_call(_metric_ready, gt_positions):
        n_reported = int(metric_ready["is_reported"].sum())
        n_universe = int(len(metric_ready))
        n_gt_total = int(len(gt_positions))
        n_gt_reported = int(metric_ready[metric_ready["position"].isin(gt_positions)]["is_reported"].sum())
        return pd.DataFrame(
            [
                {
                    "reference": "r1",
                    "tool": "nanocompore",
                    "replicate": "rep1",
                    "n_universe": n_universe,
                    "n_reported": n_reported,
                    "n_no_call": n_universe - n_reported,
                    "no_call_rate": float((n_universe - n_reported) / n_universe),
                    "n_gt_total": n_gt_total,
                    "n_gt_reported": n_gt_reported,
                    "gt_recall_raw": float(n_gt_reported / n_gt_total) if n_gt_total else np.nan,
                }
            ]
        )

    monkeypatch.setattr(_mod, "validate_all_outputs", lambda *_args, **_kwargs: _DummyQualityReport())
    monkeypatch.setattr(_mod, "generate_quality_summary", lambda *_args, **_kwargs: "ok")
    monkeypatch.setattr(_mod, "standardize_to_reference_universe", _fake_standardize)
    monkeypatch.setattr(_mod, "build_metric_ready_table", _fake_build_metric_ready)
    monkeypatch.setattr(_mod, "compute_no_call_diagnostics", _fake_no_call)
    monkeypatch.setattr(
        _mod,
        "calculate_metrics",
        lambda df, ground_truth=None, reference=None: _DummyMetrics(
            tool=str(df["tool"].iloc[0]),
            reference=str(reference or df["reference"].iloc[0]),
            replicate=str(df["replicate"].iloc[0]),
            n_predictions=len(df),
        ),
    )
    monkeypatch.setattr(_mod, "calculate_metrics_all_tools", lambda *_args, **_kwargs: pd.DataFrame())
    monkeypatch.setattr(_mod, "generate_upset_data", lambda *_args, **_kwargs: pd.DataFrame())
    monkeypatch.setattr(_mod, "calculate_pairwise_agreement", lambda *_args, **_kwargs: pd.DataFrame())
    monkeypatch.setattr(_mod, "get_consensus_positions", lambda *_args, **_kwargs: pd.DataFrame())
    monkeypatch.setattr(_mod, "compare_tools", lambda *_args, **_kwargs: _DummyComparison())
    monkeypatch.setattr(
        _mod,
        "_extract_curve_points",
        lambda *_args, **_kwargs: {
            "pr": pd.DataFrame(
                [
                    {
                        "reference": "r1",
                        "tool": "nanocompore",
                        "replicate": "rep1",
                        "scope": "universe",
                        "index": 0,
                        "precision": 0.5,
                        "recall": 0.5,
                        "threshold": 0.5,
                    }
                ]
            ),
            "roc": pd.DataFrame(
                [
                    {
                        "reference": "r1",
                        "tool": "nanocompore",
                        "replicate": "rep1",
                        "scope": "universe",
                        "index": 0,
                        "fpr": 0.0,
                        "tpr": 0.0,
                        "threshold": 1.0,
                    }
                ]
            ),
        },
    )
    monkeypatch.setattr(_mod, "plot_score_distributions", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(_mod, "generate_html_report", lambda *_args, **_kwargs: str(tmp_path / "report" / "index.html"))

    out_dir = tmp_path / "out"
    _mod.run_analysis(
        tool_outputs={"nanocompore": tool_df},
        ground_truth=ground_truth,
        output_dir=out_dir,
        references=["r1"],
        min_replicates=1,
        expected_replicates=1,
    )

    metrics = pd.read_csv(out_dir / "by_reference" / "r1" / "04_metrics" / "metrics_per_replicate.csv")
    assert len(metrics) == 1
    row = metrics.iloc[0]
    assert str(row["metrics_valid"]).lower() in {"false", "0"}
    assert float(row["invalid_score_fraction"]) > 0.5
    assert row["metric_scope_note"] == "invalid_scores"
    assert pd.notna(row["auprc"])
    assert pd.notna(row["auroc"])
    assert pd.notna(row["f1_optimal"])

    pr = pd.read_csv(out_dir / "by_reference" / "r1" / "05_curves" / "pr_curve_points.csv")
    roc = pd.read_csv(out_dir / "by_reference" / "r1" / "05_curves" / "roc_curve_points.csv")
    assert not pr.empty
    assert not roc.empty

    window = pd.read_csv(out_dir / "by_reference" / "r1" / "09_tolerance_analysis" / "window_metrics_per_replicate.csv")
    lag = pd.read_csv(out_dir / "by_reference" / "r1" / "10_lag_diagnostics" / "lag_metrics_per_replicate.csv")
    assert not window.empty
    assert not lag.empty
    assert (window["metric_scope_note"] == "invalid_scores").all()
    assert (lag["metric_scope_note"] == "invalid_scores").all()
    assert pd.to_numeric(window["auroc"], errors="coerce").notna().any()
    assert pd.to_numeric(lag["auroc"], errors="coerce").notna().any()
