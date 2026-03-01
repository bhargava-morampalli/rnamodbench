#!/usr/bin/env python3
"""
Shared reporting helpers for per-run tool availability and historical backfill.
"""

from __future__ import annotations

import csv
import html
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

BIN_DIR = Path(__file__).resolve().parent
DOWNSTREAM_DIR = BIN_DIR / "downstream_analysis"
if str(DOWNSTREAM_DIR) not in sys.path:
    sys.path.insert(0, str(DOWNSTREAM_DIR))

from parse_outputs import _discover_tool_inputs, parse_tool_output  # type: ignore


TOOL_PROCESS_FAMILIES: Dict[str, List[str]] = {
    "tombo": ["TOMBO_DETECT_MODIFICATIONS", "TOMBO_TEXT_OUTPUT"],
    "yanocomp": ["YANOCOMP_PREPARE", "YANOCOMP_ANALYSIS"],
    "nanocompore": ["NANOCOMPORE_EVENTALIGN_COLLAPSE", "NANOCOMPORE_SAMPCOMP"],
    "xpore": ["XPORE_DATAPREP", "XPORE_DIFFMOD"],
    "eligos": ["ELIGOS_PAIR_DIFF_MOD"],
    "epinano": ["EPINANO_ERROR"],
    "differr": ["DIFFERR"],
    "drummer": ["DRUMMER"],
    "jacusa2": ["JACUSA2"],
}

TOOL_ORDER = [
    "tombo",
    "yanocomp",
    "nanocompore",
    "xpore",
    "eligos",
    "epinano",
    "differr",
    "drummer",
    "jacusa2",
]

PROCESS_STATUS_COLUMNS = [
    "run_name",
    "run_dir",
    "run_id",
    "coverage_label",
    "quality_label",
    "task_id",
    "native_id",
    "process_name",
    "process_family",
    "tool",
    "target",
    "replicate",
    "status",
    "exit_code",
    "error_action",
    "workdir",
    "trace_file",
]

LOG_EVENT_COLUMNS = [
    "run_name",
    "run_dir",
    "coverage_label",
    "tool",
    "target",
    "replicate",
    "log_file",
    "event_level",
    "event_code",
    "event_text",
    "excerpt",
    "reason_source",
]

TOOL_AVAILABILITY_COLUMNS = [
    "run_name",
    "run_dir",
    "run_id",
    "coverage_label",
    "quality_label",
    "tool",
    "target",
    "replicate",
    "nf_status",
    "exit_code_max",
    "required_processes",
    "observed_processes",
    "failed_processes",
    "output_state",
    "parsed_rows",
    "raw_output_paths",
    "log_file",
    "failure_reason_code",
    "failure_reason_text",
    "reason_source",
    "diagnostic_completeness",
]

RUN_SUMMARY_COLUMNS = [
    "run_name",
    "run_id",
    "coverage_label",
    "quality_label",
    "tool",
    "total_combinations",
    "parsed_nonempty",
    "parsed_empty",
    "raw_present_unparsed",
    "final_output_missing",
    "no_observed_artifacts",
    "failed_or_aborted",
    "partial_or_missing_trace",
    "warning_events",
    "error_events",
]

FAILURE_SUMMARY_COLUMNS = [
    "coverage_label",
    "tool",
    "total_combinations",
    "parsed_nonempty",
    "parsed_empty",
    "raw_present_unparsed",
    "final_output_missing",
    "no_observed_artifacts",
    "failed_or_aborted",
    "partial_or_missing_trace",
    "warning_events",
    "error_events",
    "success_rate",
]

LOG_PATTERN_SPECS: List[Tuple[str, str, re.Pattern[str]]] = [
    ("traceback", "ERROR", re.compile(r"traceback", re.IGNORECASE)),
    ("exception", "ERROR", re.compile(r"\bexception\b", re.IGNORECASE)),
    ("errno", "ERROR", re.compile(r"\berrno\b", re.IGNORECASE)),
    ("no_reads", "ERROR", re.compile(r"\bno reads?\b", re.IGNORECASE)),
    ("insufficient", "ERROR", re.compile(r"\binsufficient\b", re.IGNORECASE)),
    ("cannot", "ERROR", re.compile(r"\bcannot\b", re.IGNORECASE)),
    ("empty", "ERROR", re.compile(r"\bempty\b", re.IGNORECASE)),
    ("failed", "ERROR", re.compile(r"\bfailed\b", re.IGNORECASE)),
    ("error", "ERROR", re.compile(r"\berror\b", re.IGNORECASE)),
    ("warning", "WARNING", re.compile(r"\bwarning\b|\bwarn\b", re.IGNORECASE)),
]

EVENT_SEVERITY_RANK = {"ERROR": 0, "WARNING": 1}
TARGET_ORDER = {"16s": 0, "23s": 1, "5s": 2}
PROCESS_TO_TOOL = {
    process_family: tool
    for tool, families in TOOL_PROCESS_FAMILIES.items()
    for process_family in families
}


@dataclass(frozen=True)
class RunContext:
    run_name: str
    run_dir: Path
    run_id: str
    coverage_label: str
    quality_label: str
    trace_file: Optional[Path]
    extra_tools: Tuple[str, ...] = ()


def _empty_dataframe(columns: Sequence[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=list(columns))


def discover_run_dirs(runs_root: Path, run_glob: str) -> List[Path]:
    return sorted(path for path in runs_root.glob(run_glob) if path.is_dir())


def infer_coverage_label(run_name: str) -> str:
    match = re.search(r"results_(\d+)x", run_name, flags=re.IGNORECASE)
    if match:
        return f"{match.group(1)}x"
    match = re.search(r"(\d+)x", run_name, flags=re.IGNORECASE)
    if match:
        return f"{match.group(1)}x"
    return "unknown"


def _infer_target_from_text(text: str) -> Optional[str]:
    lowered = text.lower()
    for token in ("16s", "23s", "5s"):
        if token in lowered:
            return token
    return None


def _infer_target_from_reference(value: object) -> Optional[str]:
    if value is None:
        return None
    return _infer_target_from_text(str(value))


def _infer_replicate(text: str) -> Optional[str]:
    match = re.search(r"rep(?:licate)?[_-]?(\d+)", text, flags=re.IGNORECASE)
    if match:
        return f"rep{match.group(1)}"
    match = re.search(r"(?:^|[_-])r(\d+)(?:$|[_-])", text, flags=re.IGNORECASE)
    if match:
        return f"rep{match.group(1)}"
    return None


def _key_sort(target: str, replicate: str) -> Tuple[int, str, int, str]:
    rep_match = re.search(r"(\d+)", replicate or "")
    rep_num = int(rep_match.group(1)) if rep_match else 999
    return (TARGET_ORDER.get(target, 99), target, rep_num, replicate)


def _safe_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def _relative_to(path: Path, base: Path) -> str:
    try:
        return str(path.relative_to(base))
    except ValueError:
        return str(path)


def _latest_trace_file(run_dir: Path) -> Optional[Path]:
    traces = sorted((run_dir / "pipeline_info").glob("execution_trace_*.txt"))
    if not traces:
        return None
    return max(traces, key=lambda path: path.stat().st_mtime)


def _score_metadata_candidate(path: Path) -> Tuple[float, int]:
    name = path.parent.parent.name
    if name.startswith("downstream_") and "backup" not in name:
        priority = 2
    elif name == "downstream_analysis":
        priority = 1
    else:
        priority = 0
    return (path.stat().st_mtime, priority)


def _load_run_metadata(run_dir: Path) -> Dict[str, str]:
    candidates = list(run_dir.glob("downstream*/metadata/run_metadata.json"))
    if not candidates:
        return {}
    chosen = max(candidates, key=_score_metadata_candidate)
    try:
        payload = json.loads(_safe_text(chosen))
    except Exception:
        return {}
    if not isinstance(payload, dict):
        return {}
    return {str(key): str(value) for key, value in payload.items() if value is not None}


def build_run_context(run_dir: Path) -> RunContext:
    run_name = run_dir.name
    metadata = _load_run_metadata(run_dir)
    coverage_label = metadata.get("coverage_label") or infer_coverage_label(run_name)
    run_id = metadata.get("run_id") or run_name
    quality_label = metadata.get("quality_label") or "unknown"
    trace_file = _latest_trace_file(run_dir)

    modifications_dir = run_dir / "modifications"
    extra_tools: List[str] = []
    if modifications_dir.exists():
        for child in sorted(modifications_dir.iterdir()):
            if child.is_dir() and child.name not in TOOL_ORDER:
                extra_tools.append(child.name)

    return RunContext(
        run_name=run_name,
        run_dir=run_dir,
        run_id=run_id,
        coverage_label=coverage_label,
        quality_label=quality_label,
        trace_file=trace_file,
        extra_tools=tuple(extra_tools),
    )


def parse_trace_name(name: str) -> Dict[str, Optional[str]]:
    leaf = name.split(":")[-1].strip()
    match = re.match(r"(?P<family>.+?)\s*\((?P<label>[^()]*)\)\s*$", leaf)
    if match:
        process_family = match.group("family").strip()
        label = match.group("label").strip()
    else:
        process_family = leaf
        label = leaf

    return {
        "process_family": process_family,
        "label": label,
        "target": _infer_target_from_text(label or name),
        "replicate": _infer_replicate(label or name),
    }


def load_process_status(run_context: RunContext) -> pd.DataFrame:
    trace_file = run_context.trace_file
    if trace_file is None or not trace_file.exists():
        return _empty_dataframe(PROCESS_STATUS_COLUMNS)

    rows: List[Dict[str, object]] = []
    with open(trace_file, "r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for raw in reader:
            parsed = parse_trace_name(raw.get("name", ""))
            process_family = parsed["process_family"]
            tool = PROCESS_TO_TOOL.get(process_family or "")
            if tool is None:
                continue
            rows.append(
                {
                    "run_name": run_context.run_name,
                    "run_dir": str(run_context.run_dir),
                    "run_id": run_context.run_id,
                    "coverage_label": run_context.coverage_label,
                    "quality_label": run_context.quality_label,
                    "task_id": raw.get("task_id", ""),
                    "native_id": raw.get("native_id", ""),
                    "process_name": raw.get("name", ""),
                    "process_family": process_family,
                    "tool": tool,
                    "target": parsed["target"] or "",
                    "replicate": parsed["replicate"] or "",
                    "status": raw.get("status", ""),
                    "exit_code": raw.get("exit", ""),
                    "error_action": raw.get("error_action", ""),
                    "workdir": raw.get("workdir", ""),
                    "trace_file": str(trace_file),
                }
            )

    df = pd.DataFrame(rows)
    if df.empty:
        return _empty_dataframe(PROCESS_STATUS_COLUMNS)
    return df.reindex(columns=PROCESS_STATUS_COLUMNS)


def _infer_tool_from_log(rel_path: Path) -> Optional[str]:
    parts = rel_path.parts
    if len(parts) >= 2 and parts[0] == "logs":
        candidate = parts[1].lower()
        if candidate in TOOL_ORDER:
            return candidate
    lowered = str(rel_path).lower()
    for tool in TOOL_ORDER:
        if f"/{tool}/" in lowered or lowered.startswith(f"logs/{tool}/"):
            return tool
    return None


def _classify_log_event(line: str) -> Optional[Tuple[str, str]]:
    lowered = line.lower().strip()
    if not lowered:
        return None
    if lowered.startswith("===") and lowered.endswith("==="):
        return None
    if lowered.startswith("error threshold:"):
        return None
    if lowered == "warning message:":
        return None
    if lowered == "failed to locate timezone database":
        return None
    for event_code, event_level, pattern in LOG_PATTERN_SPECS:
        if pattern.search(line):
            return event_code, event_level
    return None


def _event_excerpt(lines: Sequence[str], index: int, line: str) -> str:
    cleaned = line.strip()
    if re.fullmatch(r"[*\s]*warning[*\s]*", cleaned, flags=re.IGNORECASE):
        for next_line in lines[index + 1 : index + 4]:
            next_cleaned = next_line.strip()
            if next_cleaned:
                return f"{cleaned} {next_cleaned}"[:500]
    return cleaned[:500]


def discover_logs(run_context: RunContext) -> Tuple[pd.DataFrame, pd.DataFrame]:
    logs_root = run_context.run_dir / "logs"
    if not logs_root.exists():
        log_records = pd.DataFrame(columns=["tool", "target", "replicate", "log_file"])
        return log_records, _empty_dataframe(LOG_EVENT_COLUMNS)

    log_rows: List[Dict[str, object]] = []
    event_rows: List[Dict[str, object]] = []

    for log_path in sorted(logs_root.rglob("*.log")):
        rel_path = log_path.relative_to(run_context.run_dir)
        rel_str = str(rel_path)
        tool = _infer_tool_from_log(rel_path)
        target = _infer_target_from_text(rel_str) or ""
        replicate = _infer_replicate(rel_str) or ""
        log_rows.append(
            {
                "tool": tool or "",
                "target": target,
                "replicate": replicate,
                "log_file": rel_str,
            }
        )

        seen_events = set()
        lines = _safe_text(log_path).splitlines()
        for index, raw_line in enumerate(lines):
            line = raw_line.strip()
            if not line:
                continue
            classified = _classify_log_event(line)
            if classified is None:
                continue
            event_code, event_level = classified
            dedupe_key = (event_code, line)
            if dedupe_key in seen_events:
                continue
            seen_events.add(dedupe_key)
            event_rows.append(
                {
                    "run_name": run_context.run_name,
                    "run_dir": str(run_context.run_dir),
                    "coverage_label": run_context.coverage_label,
                    "tool": tool or "",
                    "target": target,
                    "replicate": replicate,
                    "log_file": rel_str,
                    "event_level": event_level,
                    "event_code": event_code,
                    "event_text": _event_excerpt(lines, index, line),
                    "excerpt": _event_excerpt(lines, index, line),
                    "reason_source": "log",
                }
            )

    log_df = pd.DataFrame(log_rows)
    if log_df.empty:
        log_df = pd.DataFrame(columns=["tool", "target", "replicate", "log_file"])
    event_df = pd.DataFrame(event_rows)
    if event_df.empty:
        event_df = _empty_dataframe(LOG_EVENT_COLUMNS)
    else:
        event_df = event_df.reindex(columns=LOG_EVENT_COLUMNS)
    return log_df, event_df


def _infer_output_key_from_path(path: Path, run_dir: Path) -> Tuple[Optional[str], Optional[str]]:
    rel_str = _relative_to(path, run_dir)
    return _infer_target_from_text(rel_str), _infer_replicate(rel_str)


def discover_output_evidence(run_context: RunContext) -> pd.DataFrame:
    modifications_dir = run_context.run_dir / "modifications"
    rows: Dict[Tuple[str, str, str], Dict[str, object]] = {}

    if not modifications_dir.exists():
        return pd.DataFrame(
            columns=[
                "tool",
                "target",
                "replicate",
                "parsed_rows",
                "parse_failures",
                "raw_output_paths",
                "parse_attempts",
            ]
        )

    for tool in TOOL_ORDER:
        tool_dir = modifications_dir / tool
        if not tool_dir.exists():
            continue

        inputs = _discover_tool_inputs(tool, tool_dir)
        for input_path in inputs:
            target_hint, replicate_hint = _infer_output_key_from_path(input_path, run_context.run_dir)

            parse_failed = False
            parsed_groups: List[Tuple[str, str, int]] = []
            try:
                parsed = parse_tool_output(tool, input_path, coverage=run_context.coverage_label)
                if not parsed.empty:
                    temp = parsed.copy()
                    temp["target"] = temp["reference"].map(_infer_target_from_reference)
                    temp["replicate"] = temp["replicate"].astype(str)
                    grouped = (
                        temp[temp["target"].notna()]
                        .groupby(["target", "replicate"], dropna=False)
                        .size()
                        .reset_index(name="parsed_rows")
                    )
                    for _, row in grouped.iterrows():
                        parsed_groups.append(
                            (str(row["target"]), str(row["replicate"]), int(row["parsed_rows"]))
                        )
            except Exception:
                parse_failed = True

            rel_path = _relative_to(input_path, run_context.run_dir)
            if parsed_groups:
                for target, replicate, parsed_rows in parsed_groups:
                    key = (tool, target, replicate)
                    row = rows.setdefault(
                        key,
                        {
                            "tool": tool,
                            "target": target,
                            "replicate": replicate,
                            "parsed_rows": 0,
                            "parse_failures": 0,
                            "raw_output_paths": set(),
                            "parse_attempts": 0,
                        },
                    )
                    row["parsed_rows"] = int(row["parsed_rows"]) + parsed_rows
                    row["parse_attempts"] = int(row["parse_attempts"]) + 1
                    row["raw_output_paths"].add(rel_path)
                    if parse_failed:
                        row["parse_failures"] = int(row["parse_failures"]) + 1
            else:
                key = (tool, target_hint or "", replicate_hint or "")
                row = rows.setdefault(
                    key,
                    {
                        "tool": tool,
                        "target": target_hint or "",
                        "replicate": replicate_hint or "",
                        "parsed_rows": 0,
                        "parse_failures": 0,
                        "raw_output_paths": set(),
                        "parse_attempts": 0,
                    },
                )
                row["parse_attempts"] = int(row["parse_attempts"]) + 1
                row["raw_output_paths"].add(rel_path)
                if parse_failed:
                    row["parse_failures"] = int(row["parse_failures"]) + 1

    output_rows: List[Dict[str, object]] = []
    for row in rows.values():
        output_rows.append(
            {
                "tool": row["tool"],
                "target": row["target"],
                "replicate": row["replicate"],
                "parsed_rows": int(row["parsed_rows"]),
                "parse_failures": int(row["parse_failures"]),
                "raw_output_paths": "; ".join(sorted(row["raw_output_paths"])),
                "parse_attempts": int(row["parse_attempts"]),
            }
        )

    df = pd.DataFrame(output_rows)
    if df.empty:
        return pd.DataFrame(
            columns=[
                "tool",
                "target",
                "replicate",
                "parsed_rows",
                "parse_failures",
                "raw_output_paths",
                "parse_attempts",
            ]
        )
    return df.sort_values(["tool", "target", "replicate"]).reset_index(drop=True)


def _collect_expected_keys(
    process_df: pd.DataFrame,
    log_df: pd.DataFrame,
    output_df: pd.DataFrame,
) -> List[Tuple[str, str]]:
    keys = set()
    for df in (process_df, log_df, output_df):
        if df.empty:
            continue
        for _, row in df.iterrows():
            target = str(row.get("target", "") or "")
            replicate = str(row.get("replicate", "") or "")
            if not target or not replicate:
                continue
            keys.add((target, replicate))

    return sorted(keys, key=lambda item: _key_sort(item[0], item[1]))


def _collapse_statuses(statuses: Sequence[str], observed_count: int, required_count: int) -> str:
    if observed_count == 0:
        return "MISSING_TRACE"
    status_set = {status for status in statuses if status}
    if "FAILED" in status_set:
        return "FAILED"
    if "ABORTED" in status_set:
        return "ABORTED"
    if observed_count < required_count:
        return "PARTIAL"
    if status_set == {"COMPLETED"}:
        return "COMPLETED"
    return "PARTIAL"


def _exit_code_max(values: Iterable[object]) -> str:
    numeric: List[int] = []
    for value in values:
        try:
            numeric.append(int(value))
        except Exception:
            continue
    if not numeric:
        return ""
    return str(max(numeric))


def _best_event(events_df: pd.DataFrame) -> Optional[pd.Series]:
    if events_df.empty:
        return None
    ranked = events_df.copy()
    ranked["_rank"] = ranked["event_level"].map(EVENT_SEVERITY_RANK).fillna(99)
    ranked = ranked.sort_values(["_rank", "event_code", "log_file", "event_text"])
    return ranked.iloc[0]


def build_tool_availability(
    run_context: RunContext,
    process_df: pd.DataFrame,
    log_df: pd.DataFrame,
    events_df: pd.DataFrame,
    output_df: pd.DataFrame,
) -> pd.DataFrame:
    expected_keys = _collect_expected_keys(process_df, log_df, output_df)
    if not expected_keys:
        return _empty_dataframe(TOOL_AVAILABILITY_COLUMNS)

    rows: List[Dict[str, object]] = []

    for tool in TOOL_ORDER:
        required = TOOL_PROCESS_FAMILIES[tool]
        for target, replicate in expected_keys:
            proc_sub = process_df[
                (process_df["tool"] == tool)
                & (process_df["target"] == target)
                & (process_df["replicate"] == replicate)
            ].copy()
            log_sub = log_df[
                (log_df["tool"] == tool)
                & (log_df["target"] == target)
                & (log_df["replicate"] == replicate)
            ].copy()
            evt_sub = events_df[
                (events_df["tool"] == tool)
                & (events_df["target"] == target)
                & (events_df["replicate"] == replicate)
            ].copy()
            out_sub = output_df[
                (output_df["tool"] == tool)
                & (output_df["target"] == target)
                & (output_df["replicate"] == replicate)
            ].copy()

            observed_processes = sorted(proc_sub["process_family"].dropna().astype(str).unique().tolist())
            failed_processes = sorted(
                proc_sub[proc_sub["status"].isin(["FAILED", "ABORTED"])]["process_family"]
                .dropna()
                .astype(str)
                .unique()
                .tolist()
            )
            nf_status = _collapse_statuses(
                proc_sub["status"].dropna().astype(str).tolist(),
                observed_count=len(observed_processes),
                required_count=len(required),
            )

            parsed_rows = int(out_sub["parsed_rows"].sum()) if not out_sub.empty else 0
            parse_failures = int(out_sub["parse_failures"].sum()) if not out_sub.empty else 0
            raw_paths: List[str] = []
            if not out_sub.empty:
                for raw_value in out_sub["raw_output_paths"].dropna().astype(str):
                    raw_paths.extend([item.strip() for item in raw_value.split(";") if item.strip()])
            raw_paths = sorted(set(raw_paths))

            if parsed_rows > 0:
                output_state = "parsed_nonempty"
            elif raw_paths:
                output_state = "raw_present_unparsed" if parse_failures > 0 else "parsed_empty"
            elif not proc_sub.empty or not log_sub.empty:
                output_state = "final_output_missing"
            else:
                output_state = "no_observed_artifacts"

            best_event = _best_event(evt_sub)
            best_event_level = str(best_event["event_level"]) if best_event is not None else ""
            if best_event is not None and best_event_level == "ERROR":
                failure_reason_code = str(best_event["event_code"])
                failure_reason_text = str(best_event["event_text"])
                reason_source = "log"
            elif (
                best_event is not None
                and best_event_level == "WARNING"
                and nf_status == "COMPLETED"
                and output_state == "parsed_nonempty"
            ):
                failure_reason_code = str(best_event["event_code"])
                failure_reason_text = str(best_event["event_text"])
                reason_source = "log"
            elif nf_status in {"FAILED", "ABORTED"}:
                failure_reason_code = "process_failed"
                failure_reason_text = f"{tool} trace status {nf_status.lower()}"
                reason_source = "trace"
            elif nf_status == "PARTIAL":
                failure_reason_code = "partial_trace"
                failure_reason_text = "Not all required process families were observed in trace"
                reason_source = "trace"
            elif output_state == "final_output_missing":
                failure_reason_code = "final_output_missing"
                failure_reason_text = "No final published output was found for observed run artifacts"
                reason_source = "output_scan"
            elif output_state == "parsed_empty":
                failure_reason_code = "parsed_empty"
                failure_reason_text = "Final output exists but parsed to zero rows"
                reason_source = "output_scan"
            elif output_state == "raw_present_unparsed":
                failure_reason_code = "raw_present_unparsed"
                failure_reason_text = "Final output exists but parser could not recover usable rows"
                reason_source = "output_scan"
            elif output_state == "no_observed_artifacts":
                failure_reason_code = "no_observed_artifacts"
                failure_reason_text = "No trace, log, or output artifact was observed for this tool/key"
                reason_source = "unknown"
            else:
                failure_reason_code = ""
                failure_reason_text = ""
                reason_source = "unknown"

            has_trace = not proc_sub.empty
            has_log = not log_sub.empty
            has_output = parsed_rows > 0 or bool(raw_paths)
            if has_trace and has_log:
                diagnostic_completeness = "trace_and_log"
            elif has_trace:
                diagnostic_completeness = "trace_only"
            elif has_log:
                diagnostic_completeness = "log_only"
            elif has_output:
                diagnostic_completeness = "output_only"
            else:
                diagnostic_completeness = "limited"

            rows.append(
                {
                    "run_name": run_context.run_name,
                    "run_dir": str(run_context.run_dir),
                    "run_id": run_context.run_id,
                    "coverage_label": run_context.coverage_label,
                    "quality_label": run_context.quality_label,
                    "tool": tool,
                    "target": target,
                    "replicate": replicate,
                    "nf_status": nf_status,
                    "exit_code_max": _exit_code_max(proc_sub["exit_code"]) if not proc_sub.empty else "",
                    "required_processes": "; ".join(required),
                    "observed_processes": "; ".join(observed_processes),
                    "failed_processes": "; ".join(failed_processes),
                    "output_state": output_state,
                    "parsed_rows": parsed_rows,
                    "raw_output_paths": "; ".join(raw_paths),
                    "log_file": "; ".join(sorted(log_sub["log_file"].dropna().astype(str).unique().tolist())),
                    "failure_reason_code": failure_reason_code,
                    "failure_reason_text": failure_reason_text,
                    "reason_source": reason_source,
                    "diagnostic_completeness": diagnostic_completeness,
                }
            )

    availability_df = pd.DataFrame(rows)
    if availability_df.empty:
        return _empty_dataframe(TOOL_AVAILABILITY_COLUMNS)
    return availability_df.reindex(columns=TOOL_AVAILABILITY_COLUMNS)


def summarize_run_availability(
    availability_df: pd.DataFrame,
    events_df: pd.DataFrame,
    run_context: RunContext,
) -> pd.DataFrame:
    if availability_df.empty:
        return _empty_dataframe(RUN_SUMMARY_COLUMNS)

    event_counts = (
        events_df.groupby(["tool", "target", "replicate", "event_level"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
        if not events_df.empty
        else pd.DataFrame(columns=["tool", "target", "replicate", "ERROR", "WARNING"])
    )

    merged = availability_df.merge(
        event_counts,
        on=["tool", "target", "replicate"],
        how="left",
    )
    for column in ("ERROR", "WARNING"):
        if column not in merged.columns:
            merged[column] = 0
        merged[column] = merged[column].fillna(0).astype(int)

    rows: List[Dict[str, object]] = []
    for tool, group in merged.groupby("tool", dropna=False):
        rows.append(
            {
                "run_name": run_context.run_name,
                "run_id": run_context.run_id,
                "coverage_label": run_context.coverage_label,
                "quality_label": run_context.quality_label,
                "tool": tool,
                "total_combinations": int(len(group)),
                "parsed_nonempty": int((group["output_state"] == "parsed_nonempty").sum()),
                "parsed_empty": int((group["output_state"] == "parsed_empty").sum()),
                "raw_present_unparsed": int((group["output_state"] == "raw_present_unparsed").sum()),
                "final_output_missing": int((group["output_state"] == "final_output_missing").sum()),
                "no_observed_artifacts": int((group["output_state"] == "no_observed_artifacts").sum()),
                "failed_or_aborted": int(group["nf_status"].isin(["FAILED", "ABORTED"]).sum()),
                "partial_or_missing_trace": int(group["nf_status"].isin(["PARTIAL", "MISSING_TRACE"]).sum()),
                "warning_events": int(group["WARNING"].sum()),
                "error_events": int(group["ERROR"].sum()),
            }
        )

    return pd.DataFrame(rows).reindex(columns=RUN_SUMMARY_COLUMNS).sort_values("tool")


def build_run_report(run_dir: Path) -> Dict[str, object]:
    run_context = build_run_context(run_dir)
    process_df = load_process_status(run_context)
    log_df, events_df = discover_logs(run_context)
    output_df = discover_output_evidence(run_context)
    availability_df = build_tool_availability(run_context, process_df, log_df, events_df, output_df)
    summary_df = summarize_run_availability(availability_df, events_df, run_context)

    return {
        "context": run_context,
        "process_status": process_df,
        "log_index": log_df,
        "log_events": events_df,
        "output_evidence": output_df,
        "tool_availability": availability_df,
        "summary": summary_df,
    }


def _html_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "<p>No rows.</p>"
    headers = "".join(f"<th>{html.escape(str(col))}</th>" for col in df.columns)
    rows = []
    for _, row in df.iterrows():
        cells = "".join(f"<td>{html.escape(str(value))}</td>" for value in row.tolist())
        rows.append(f"<tr>{cells}</tr>")
    return f"<table><thead><tr>{headers}</tr></thead><tbody>{''.join(rows)}</tbody></table>"


def render_run_html(report: Dict[str, object], summary_df: pd.DataFrame) -> str:
    run_context: RunContext = report["context"]  # type: ignore[assignment]
    availability_df: pd.DataFrame = report["tool_availability"]  # type: ignore[assignment]
    events_df: pd.DataFrame = report["log_events"]  # type: ignore[assignment]

    issues_df = availability_df[
        (availability_df["output_state"] != "parsed_nonempty")
        | (availability_df["nf_status"] != "COMPLETED")
        | (availability_df["reason_source"] == "log")
    ].copy()
    issues_df = issues_df[
        [
            "tool",
            "target",
            "replicate",
            "nf_status",
            "output_state",
            "failure_reason_code",
            "failure_reason_text",
            "log_file",
        ]
    ]

    counts = {
        "rows": len(availability_df),
        "usable": int((availability_df["output_state"] == "parsed_nonempty").sum()),
        "problem_rows": len(issues_df),
        "error_events": int((events_df["event_level"] == "ERROR").sum()) if not events_df.empty else 0,
        "warning_events": int((events_df["event_level"] == "WARNING").sum()) if not events_df.empty else 0,
    }

    extra_tools_html = ""
    if run_context.extra_tools:
        extra_tools_html = (
            "<p><strong>Observed extra tool directories:</strong> "
            + html.escape(", ".join(run_context.extra_tools))
            + "</p>"
        )

    return """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Tool Availability Report</title>
  <style>
    body { font-family: sans-serif; margin: 24px; color: #1f2937; }
    h1, h2 { margin-bottom: 0.3rem; }
    .cards { display: flex; gap: 16px; margin: 16px 0 24px; flex-wrap: wrap; }
    .card { background: #f3f4f6; border-radius: 8px; padding: 12px 16px; min-width: 160px; }
    .card .value { font-size: 1.8rem; font-weight: 700; }
    table { border-collapse: collapse; width: 100%; margin: 16px 0 24px; }
    th, td { border: 1px solid #d1d5db; padding: 8px 10px; text-align: left; vertical-align: top; }
    th { background: #e5e7eb; }
    code { background: #f3f4f6; padding: 1px 4px; border-radius: 4px; }
  </style>
</head>
<body>
"""
    + f"<h1>Tool Availability Report</h1><p>Run: <code>{html.escape(run_context.run_name)}</code> | Coverage: <code>{html.escape(run_context.coverage_label)}</code> | Quality: <code>{html.escape(run_context.quality_label)}</code></p>"
    + extra_tools_html
    + f"""
<div class="cards">
  <div class="card"><div>Tool/key rows</div><div class="value">{counts['rows']}</div></div>
  <div class="card"><div>Usable outputs</div><div class="value">{counts['usable']}</div></div>
  <div class="card"><div>Problem rows</div><div class="value">{counts['problem_rows']}</div></div>
  <div class="card"><div>Error events</div><div class="value">{counts['error_events']}</div></div>
  <div class="card"><div>Warning events</div><div class="value">{counts['warning_events']}</div></div>
</div>
<h2>Per-Tool Summary</h2>
{_html_table(summary_df)}
<h2>Problem Rows</h2>
{_html_table(issues_df)}
<h2>Log Events</h2>
{_html_table(events_df)}
</body>
</html>
"""


def write_run_report(
    run_dir: Path,
    *,
    output_prefix: str = "error_summary",
) -> Dict[str, object]:
    report = build_run_report(run_dir)
    run_context: RunContext = report["context"]  # type: ignore[assignment]
    process_df: pd.DataFrame = report["process_status"]  # type: ignore[assignment]
    events_df: pd.DataFrame = report["log_events"]  # type: ignore[assignment]
    availability_df: pd.DataFrame = report["tool_availability"]  # type: ignore[assignment]
    summary_df: pd.DataFrame = report["summary"]  # type: ignore[assignment]

    pipeline_info = run_dir / "pipeline_info"
    pipeline_info.mkdir(parents=True, exist_ok=True)

    process_path = pipeline_info / "process_status.tsv"
    events_path = pipeline_info / "log_events.tsv"
    availability_path = pipeline_info / "tool_availability_per_run.tsv"
    summary_csv_path = pipeline_info / f"{output_prefix}.csv"
    summary_html_path = pipeline_info / f"{output_prefix}.html"

    process_df.to_csv(process_path, sep="\t", index=False)
    events_df.to_csv(events_path, sep="\t", index=False)
    availability_df.to_csv(availability_path, sep="\t", index=False)
    summary_df.to_csv(summary_csv_path, index=False)
    summary_html_path.write_text(render_run_html(report, summary_df), encoding="utf-8")

    report["paths"] = {
        "process_status": process_path,
        "log_events": events_path,
        "tool_availability": availability_path,
        "summary_csv": summary_csv_path,
        "summary_html": summary_html_path,
    }
    report["context"] = run_context
    return report


def load_run_reporting_artifact(run_root: Path, filename: str) -> pd.DataFrame:
    path = run_root / "pipeline_info" / filename
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.DataFrame()


def build_availability_matrix(run_roots: Sequence[Path]) -> pd.DataFrame:
    frames: List[pd.DataFrame] = []
    for run_root in run_roots:
        df = load_run_reporting_artifact(run_root, "tool_availability_per_run.tsv")
        if not df.empty:
            frames.append(df)
    if not frames:
        return _empty_dataframe(TOOL_AVAILABILITY_COLUMNS)
    merged = pd.concat(frames, ignore_index=True, sort=False)
    return merged.reindex(columns=TOOL_AVAILABILITY_COLUMNS)


def build_log_event_matrix(run_roots: Sequence[Path]) -> pd.DataFrame:
    frames: List[pd.DataFrame] = []
    for run_root in run_roots:
        df = load_run_reporting_artifact(run_root, "log_events.tsv")
        if not df.empty:
            frames.append(df)
    if not frames:
        return _empty_dataframe(LOG_EVENT_COLUMNS)
    merged = pd.concat(frames, ignore_index=True, sort=False)
    return merged.reindex(columns=LOG_EVENT_COLUMNS)


def build_failure_summary(
    availability_df: pd.DataFrame,
    log_events_df: pd.DataFrame,
) -> pd.DataFrame:
    if availability_df.empty:
        return _empty_dataframe(FAILURE_SUMMARY_COLUMNS)

    event_counts = (
        log_events_df.groupby(["coverage_label", "tool", "target", "replicate", "event_level"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
        if not log_events_df.empty
        else pd.DataFrame(columns=["coverage_label", "tool", "target", "replicate", "ERROR", "WARNING"])
    )

    merged = availability_df.merge(
        event_counts,
        on=["coverage_label", "tool", "target", "replicate"],
        how="left",
    )
    for column in ("ERROR", "WARNING"):
        if column not in merged.columns:
            merged[column] = 0
        merged[column] = merged[column].fillna(0).astype(int)

    rows: List[Dict[str, object]] = []
    for (coverage_label, tool), group in merged.groupby(["coverage_label", "tool"], dropna=False):
        total = int(len(group))
        parsed_nonempty = int((group["output_state"] == "parsed_nonempty").sum())
        rows.append(
            {
                "coverage_label": coverage_label,
                "tool": tool,
                "total_combinations": total,
                "parsed_nonempty": parsed_nonempty,
                "parsed_empty": int((group["output_state"] == "parsed_empty").sum()),
                "raw_present_unparsed": int((group["output_state"] == "raw_present_unparsed").sum()),
                "final_output_missing": int((group["output_state"] == "final_output_missing").sum()),
                "no_observed_artifacts": int((group["output_state"] == "no_observed_artifacts").sum()),
                "failed_or_aborted": int(group["nf_status"].isin(["FAILED", "ABORTED"]).sum()),
                "partial_or_missing_trace": int(group["nf_status"].isin(["PARTIAL", "MISSING_TRACE"]).sum()),
                "warning_events": int(group["WARNING"].sum()),
                "error_events": int(group["ERROR"].sum()),
                "success_rate": round(parsed_nonempty / total, 4) if total else 0.0,
            }
        )

    return pd.DataFrame(rows).reindex(columns=FAILURE_SUMMARY_COLUMNS).sort_values(
        ["coverage_label", "tool"]
    )


def render_availability_report_markdown(
    availability_df: pd.DataFrame,
    failure_summary_df: pd.DataFrame,
    log_events_df: pd.DataFrame,
) -> str:
    lines: List[str] = []
    lines.append("# Availability Report")
    lines.append("")
    lines.append("## Executive Summary")
    lines.append("")
    if availability_df.empty:
        lines.append("- No per-run reporting artifacts were available to collate.")
        return "\n".join(lines)

    runs = sorted(availability_df["run_name"].dropna().astype(str).unique().tolist())
    lines.append(f"- Runs represented: **{len(runs)}**")
    lines.append(f"- Tool/key rows: **{len(availability_df)}**")
    lines.append(
        f"- Usable outputs: **{int((availability_df['output_state'] == 'parsed_nonempty').sum())}**"
    )
    lines.append(
        f"- Failed or aborted rows: **{int(availability_df['nf_status'].isin(['FAILED', 'ABORTED']).sum())}**"
    )
    if not log_events_df.empty:
        lines.append(
            f"- Error log events: **{int((log_events_df['event_level'] == 'ERROR').sum())}**"
        )
        lines.append(
            f"- Warning log events: **{int((log_events_df['event_level'] == 'WARNING').sum())}**"
        )
    lines.append("")

    lines.append("## Per-Coverage Tool Summary")
    lines.append("")
    lines.append(
        "| Coverage | Tool | Total | Parsed Nonempty | Parsed Empty | Raw Present Unparsed | Final Output Missing | Failed/Aborted | Success Rate |"
    )
    lines.append("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for _, row in failure_summary_df.iterrows():
        lines.append(
            f"| {row['coverage_label']} | {row['tool']} | {int(row['total_combinations'])} | "
            f"{int(row['parsed_nonempty'])} | {int(row['parsed_empty'])} | "
            f"{int(row['raw_present_unparsed'])} | {int(row['final_output_missing'])} | "
            f"{int(row['failed_or_aborted'])} | {float(row['success_rate']):.2%} |"
        )
    lines.append("")

    problem_rows = availability_df[
        (availability_df["output_state"] != "parsed_nonempty")
        | (availability_df["nf_status"] != "COMPLETED")
    ].copy()
    lines.append("## Problem Rows")
    lines.append("")
    if problem_rows.empty:
        lines.append("- No problem rows detected.")
    else:
        lines.append("| Run | Coverage | Tool | Target | Replicate | Trace Status | Output State | Reason |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
        for _, row in problem_rows.sort_values(
            ["coverage_label", "tool", "target", "replicate", "run_name"]
        ).iterrows():
            reason = row["failure_reason_text"] or row["failure_reason_code"]
            lines.append(
                f"| {row['run_name']} | {row['coverage_label']} | {row['tool']} | "
                f"{row['target']} | {row['replicate']} | {row['nf_status']} | "
                f"{row['output_state']} | {reason} |"
            )
    lines.append("")
    return "\n".join(lines)
