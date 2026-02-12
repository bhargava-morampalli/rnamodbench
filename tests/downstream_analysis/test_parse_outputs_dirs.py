"""
Regression tests for directory-based tool output loading.

ELIGOS, EPINANO, and DRUMMER emit directories as their Nextflow output.
The loader must descend into these directories to find the actual result
files. Without the directory-handling fix, these tools are silently
omitted from downstream comparison (no crash, just missing data).
"""

import importlib.util
import sys
import tempfile
from pathlib import Path

import pytest
import pandas as pd

# Import parse_outputs directly from file to avoid pulling in the full
# downstream_analysis package (which has heavy dependencies like sklearn).
_mod_path = Path(__file__).resolve().parent.parent.parent / "bin" / "downstream_analysis" / "parse_outputs.py"
_spec = importlib.util.spec_from_file_location("parse_outputs", _mod_path)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

_load_standard_structure = _mod._load_standard_structure
_load_with_coverage_dirs = _mod._load_with_coverage_dirs
load_all_tool_outputs = _mod.load_all_tool_outputs


def _create_eligos_dir(parent: Path) -> None:
    """Create a minimal ELIGOS directory output."""
    d = parent / "sample1_eligos"
    d.mkdir(parents=True)
    (d / "test_paired_diff_mod_result.txt").write_text(
        "chrom\tstart_loc\tend_loc\tstrand\tpval\tadj_pval\toddR\tESB\n"
        "16S_rRNA\t100\t101\t+\t0.001\t0.01\t3.5\t0.8\n"
        "16S_rRNA\t200\t201\t+\t0.005\t0.05\t2.1\t0.6\n"
    )


def _create_epinano_dir(parent: Path) -> None:
    """Create a minimal EPINANO directory output."""
    d = parent / "sample1_epinano_error"
    d.mkdir(parents=True)
    (d / "sample_diff_err.csv").write_text(
        "Ref,pos,SumErr_z\n"
        "16S_rRNA,100,2.5\n"
        "16S_rRNA,200,3.1\n"
    )


def _create_drummer_dir(parent: Path) -> None:
    """Create a minimal DRUMMER directory output."""
    d = parent / "sample1_drummer"
    d.mkdir(parents=True)
    (d / "summary.txt").write_text(
        "chrom\tposition\tpval\todds_ratio\n"
        "16S_rRNA\t100\t0.002\t4.0\n"
        "16S_rRNA\t200\t0.01\t2.5\n"
    )


def _create_tombo_file(parent: Path) -> None:
    """Create a minimal TOMBO flat CSV file (file-based tool)."""
    (parent / "sample1_tombo.csv").write_text(
        "pos,stat\n"
        "100,0.003\n"
        "200,0.015\n"
    )


@pytest.fixture
def standard_tree(tmp_path):
    """Create a mock modifications/ tree with standard layout.

    Structure:
        tmp_path/
        ├── eligos/16S/sample1_eligos/test_paired_diff_mod_result.txt
        ├── epinano/16S/sample1_epinano_error/sample_diff_err.csv
        ├── drummer/16S/sample1_drummer/summary.txt
        └── tombo/16S/sample1_tombo.csv       (flat file, backward compat)
    """
    for tool in ["eligos", "epinano", "drummer", "tombo"]:
        (tmp_path / tool / "16S").mkdir(parents=True)

    _create_eligos_dir(tmp_path / "eligos" / "16S")
    _create_epinano_dir(tmp_path / "epinano" / "16S")
    _create_drummer_dir(tmp_path / "drummer" / "16S")
    _create_tombo_file(tmp_path / "tombo" / "16S")

    return tmp_path


@pytest.fixture
def coverage_tree(tmp_path):
    """Create a mock tree with coverage-structured layout.

    Structure:
        tmp_path/
        └── 100x/
            ├── eligos/16S/sample1_eligos/test_paired_diff_mod_result.txt
            ├── drummer/16S/sample1_drummer/summary.txt
            └── tombo/16S/sample1_tombo.csv
    """
    cov = tmp_path / "100x"
    for tool in ["eligos", "drummer", "tombo"]:
        (cov / tool / "16S").mkdir(parents=True)

    _create_eligos_dir(cov / "eligos" / "16S")
    _create_drummer_dir(cov / "drummer" / "16S")
    _create_tombo_file(cov / "tombo" / "16S")

    return tmp_path


class TestDirectoryOutputLoading:
    """Verify that directory-emitting tools are loaded by the standard loader."""

    def test_eligos_dir_loaded(self, standard_tree):
        results = _load_standard_structure(standard_tree, ["eligos"])
        assert "eligos" in results
        assert not results["eligos"].empty, "ELIGOS directory output was silently skipped"
        assert len(results["eligos"]) == 2

    def test_epinano_dir_loaded(self, standard_tree):
        results = _load_standard_structure(standard_tree, ["epinano"])
        assert "epinano" in results
        assert not results["epinano"].empty, "EPINANO directory output was silently skipped"
        assert len(results["epinano"]) == 2

    def test_drummer_dir_loaded(self, standard_tree):
        results = _load_standard_structure(standard_tree, ["drummer"])
        assert "drummer" in results
        assert not results["drummer"].empty, "DRUMMER directory output was silently skipped"
        assert len(results["drummer"]) == 2

    def test_tombo_flat_file_still_works(self, standard_tree):
        """Backward compatibility: file-based tools must still load."""
        results = _load_standard_structure(standard_tree, ["tombo"])
        assert "tombo" in results
        assert not results["tombo"].empty, "TOMBO flat file loading regressed"
        assert len(results["tombo"]) == 2

    def test_all_tools_loaded_together(self, standard_tree):
        """All four tools load when requested together."""
        tools = ["eligos", "epinano", "drummer", "tombo"]
        results = _load_standard_structure(standard_tree, tools)
        for tool in tools:
            assert tool in results, f"{tool} missing from results"
            assert not results[tool].empty, f"{tool} has empty DataFrame"

    def test_standardized_columns(self, standard_tree):
        """Directory-loaded data has the expected standardized columns."""
        results = _load_standard_structure(standard_tree, ["eligos"])
        df = results["eligos"]
        for col in ["tool", "reference", "position", "score", "score_type", "pvalue"]:
            assert col in df.columns, f"Missing column: {col}"


class TestCoverageDirectoryLoading:
    """Verify directory outputs load through the coverage-structured loader."""

    def test_eligos_dir_in_coverage_layout(self, coverage_tree):
        results = _load_with_coverage_dirs(coverage_tree, ["eligos"])
        assert "eligos" in results
        assert not results["eligos"].empty, "ELIGOS dir skipped in coverage layout"

    def test_drummer_dir_in_coverage_layout(self, coverage_tree):
        results = _load_with_coverage_dirs(coverage_tree, ["drummer"])
        assert "drummer" in results
        assert not results["drummer"].empty, "DRUMMER dir skipped in coverage layout"

    def test_tombo_flat_in_coverage_layout(self, coverage_tree):
        results = _load_with_coverage_dirs(coverage_tree, ["tombo"])
        assert "tombo" in results
        assert not results["tombo"].empty, "TOMBO flat file broken in coverage layout"


class TestLoadAllToolOutputs:
    """Integration test through the public entry point."""

    def test_standard_mode(self, standard_tree):
        results = load_all_tool_outputs(
            standard_tree, ["eligos", "drummer", "tombo"], coverage_dirs=False
        )
        for tool in ["eligos", "drummer", "tombo"]:
            assert tool in results
            assert not results[tool].empty, f"{tool} empty via load_all_tool_outputs"

    def test_coverage_mode(self, coverage_tree):
        results = load_all_tool_outputs(
            coverage_tree, ["eligos", "drummer", "tombo"], coverage_dirs=True
        )
        for tool in ["eligos", "drummer", "tombo"]:
            assert tool in results
            assert not results[tool].empty, f"{tool} empty via load_all_tool_outputs (coverage)"
