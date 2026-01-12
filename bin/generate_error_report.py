#!/usr/bin/env python3
"""
Generate Error Summary Report for RNA Modification Pipeline

This script parses Nextflow trace files and log files to generate a centralized
error summary report in both HTML and CSV formats.

Usage:
    python generate_error_report.py --outdir <results_directory>
"""

import argparse
import csv
import os
import sys
from datetime import datetime
from pathlib import Path
import re
import glob


def parse_trace_file(trace_path):
    """Parse Nextflow trace file and extract process information."""
    processes = []

    if not os.path.exists(trace_path):
        return processes

    with open(trace_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            processes.append({
                'task_id': row.get('task_id', ''),
                'hash': row.get('hash', ''),
                'name': row.get('name', ''),
                'status': row.get('status', ''),
                'exit': row.get('exit', ''),
                'submit': row.get('submit', ''),
                'start': row.get('start', ''),
                'complete': row.get('complete', ''),
                'duration': row.get('duration', ''),
                'realtime': row.get('realtime', ''),
                'cpu_percent': row.get('%cpu', ''),
                'peak_rss': row.get('peak_rss', ''),
                'peak_vmem': row.get('peak_vmem', ''),
                'workdir': row.get('workdir', ''),
                'error_action': row.get('error_action', '')
            })

    return processes


def parse_log_file(log_path):
    """Parse a log file and extract error information."""
    errors = []
    warnings = []

    if not os.path.exists(log_path):
        return {'errors': errors, 'warnings': warnings, 'content': ''}

    with open(log_path, 'r') as f:
        content = f.read()

    # Look for common error patterns
    error_patterns = [
        r'(?i)error[:\s].*',
        r'(?i)exception[:\s].*',
        r'(?i)failed[:\s].*',
        r'(?i)traceback.*',
        r'(?i)errno.*',
    ]

    warning_patterns = [
        r'(?i)warning[:\s].*',
        r'(?i)warn[:\s].*',
    ]

    for pattern in error_patterns:
        matches = re.findall(pattern, content)
        errors.extend(matches[:5])  # Limit to first 5 matches per pattern

    for pattern in warning_patterns:
        matches = re.findall(pattern, content)
        warnings.extend(matches[:5])

    return {
        'errors': list(set(errors)),  # Remove duplicates
        'warnings': list(set(warnings)),
        'content': content[:5000] if len(content) > 5000 else content  # Truncate if too long
    }


def find_log_files(logs_dir):
    """Find all log files in the logs directory."""
    log_files = {}

    if not os.path.exists(logs_dir):
        return log_files

    for log_path in glob.glob(os.path.join(logs_dir, '**', '*.log'), recursive=True):
        rel_path = os.path.relpath(log_path, logs_dir)
        parts = rel_path.split(os.sep)

        tool = parts[0] if len(parts) > 0 else 'unknown'
        rrna = parts[1] if len(parts) > 1 else 'unknown'
        filename = os.path.basename(log_path)

        if tool not in log_files:
            log_files[tool] = {}
        if rrna not in log_files[tool]:
            log_files[tool][rrna] = []

        log_info = parse_log_file(log_path)
        log_info['path'] = log_path
        log_info['filename'] = filename
        log_files[tool][rrna].append(log_info)

    return log_files


def generate_html_report(processes, log_files, output_path):
    """Generate HTML error summary report."""

    # Filter failed processes
    failed_processes = [p for p in processes if p['status'] in ['FAILED', 'ABORTED']]
    ignored_processes = [p for p in processes if p['error_action'] == 'IGNORE']
    successful_processes = [p for p in processes if p['status'] == 'COMPLETED']

    # Count modification tools
    mod_tools = ['TOMBO', 'YANOCOMP', 'NANOCOMPORE', 'XPORE', 'ELIGOS', 'EPINANO', 'DIFFERR', 'DRUMMER', 'JACUSA2', 'NANORMS']
    mod_process_stats = {tool: {'total': 0, 'failed': 0, 'ignored': 0, 'success': 0} for tool in mod_tools}

    for p in processes:
        for tool in mod_tools:
            if tool in p['name'].upper():
                mod_process_stats[tool]['total'] += 1
                if p['status'] == 'COMPLETED':
                    mod_process_stats[tool]['success'] += 1
                elif p['status'] in ['FAILED', 'ABORTED']:
                    mod_process_stats[tool]['failed'] += 1
                if p['error_action'] == 'IGNORE':
                    mod_process_stats[tool]['ignored'] += 1
                break

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>RNA Modification Pipeline - Error Summary Report</title>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 20px; background-color: #f5f5f5; }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        h3 {{ color: #7f8c8d; }}
        .summary-cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin: 20px 0; }}
        .card {{ background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .card-title {{ font-size: 14px; color: #7f8c8d; margin-bottom: 5px; }}
        .card-value {{ font-size: 32px; font-weight: bold; }}
        .card-success {{ color: #27ae60; }}
        .card-warning {{ color: #f39c12; }}
        .card-error {{ color: #e74c3c; }}
        .card-info {{ color: #3498db; }}
        table {{ width: 100%; border-collapse: collapse; background: white; margin: 15px 0; border-radius: 8px; overflow: hidden; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        th {{ background-color: #3498db; color: white; padding: 12px; text-align: left; }}
        td {{ padding: 10px 12px; border-bottom: 1px solid #ecf0f1; }}
        tr:hover {{ background-color: #f8f9fa; }}
        .status-completed {{ color: #27ae60; font-weight: bold; }}
        .status-failed {{ color: #e74c3c; font-weight: bold; }}
        .status-ignored {{ color: #f39c12; font-weight: bold; }}
        .error-box {{ background: #fee; border-left: 4px solid #e74c3c; padding: 10px; margin: 10px 0; border-radius: 4px; }}
        .warning-box {{ background: #ffeaa7; border-left: 4px solid #f39c12; padding: 10px; margin: 10px 0; border-radius: 4px; }}
        .log-content {{ background: #2c3e50; color: #ecf0f1; padding: 15px; border-radius: 8px; font-family: monospace; font-size: 12px; overflow-x: auto; max-height: 300px; overflow-y: auto; }}
        .tool-section {{ background: white; padding: 20px; margin: 15px 0; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .collapsible {{ cursor: pointer; padding: 10px; background: #ecf0f1; border: none; width: 100%; text-align: left; font-size: 16px; border-radius: 4px; }}
        .collapsible:hover {{ background: #bdc3c7; }}
        .content {{ display: none; padding: 10px; }}
        .timestamp {{ color: #7f8c8d; font-size: 12px; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>RNA Modification Pipeline - Error Summary Report</h1>
        <p class="timestamp">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>

        <h2>Pipeline Summary</h2>
        <div class="summary-cards">
            <div class="card">
                <div class="card-title">Total Processes</div>
                <div class="card-value card-info">{len(processes)}</div>
            </div>
            <div class="card">
                <div class="card-title">Completed</div>
                <div class="card-value card-success">{len(successful_processes)}</div>
            </div>
            <div class="card">
                <div class="card-title">Failed</div>
                <div class="card-value card-error">{len(failed_processes)}</div>
            </div>
            <div class="card">
                <div class="card-title">Ignored Errors</div>
                <div class="card-value card-warning">{len(ignored_processes)}</div>
            </div>
        </div>

        <h2>Modification Tools Status</h2>
        <table>
            <tr>
                <th>Tool</th>
                <th>Total Runs</th>
                <th>Successful</th>
                <th>Failed</th>
                <th>Ignored</th>
                <th>Status</th>
            </tr>
"""

    for tool, stats in mod_process_stats.items():
        if stats['total'] > 0:
            status_class = 'status-completed' if stats['failed'] == 0 else 'status-failed'
            status_text = 'OK' if stats['failed'] == 0 else 'HAS ERRORS'
            html += f"""
            <tr>
                <td><strong>{tool}</strong></td>
                <td>{stats['total']}</td>
                <td class="status-completed">{stats['success']}</td>
                <td class="status-failed">{stats['failed']}</td>
                <td class="status-ignored">{stats['ignored']}</td>
                <td class="{status_class}">{status_text}</td>
            </tr>
"""

    html += """
        </table>
"""

    # Failed processes section
    if failed_processes:
        html += """
        <h2>Failed Processes</h2>
        <table>
            <tr>
                <th>Process Name</th>
                <th>Exit Code</th>
                <th>Duration</th>
                <th>Work Directory</th>
            </tr>
"""
        for p in failed_processes:
            html += f"""
            <tr>
                <td>{p['name']}</td>
                <td class="status-failed">{p['exit']}</td>
                <td>{p['duration']}</td>
                <td><code>{p['workdir']}</code></td>
            </tr>
"""
        html += "</table>"

    # Log files section
    html += """
        <h2>Log Files by Tool</h2>
"""

    for tool, rrna_logs in log_files.items():
        html += f"""
        <div class="tool-section">
            <h3>{tool.upper()}</h3>
"""
        for rrna, logs in rrna_logs.items():
            html += f"<h4>{rrna}</h4>"
            for log in logs:
                has_errors = len(log['errors']) > 0
                error_class = 'status-failed' if has_errors else 'status-completed'

                html += f"""
                <button class="collapsible">{log['filename']} - <span class="{error_class}">{len(log['errors'])} errors, {len(log['warnings'])} warnings</span></button>
                <div class="content">
"""
                if log['errors']:
                    html += '<div class="error-box"><strong>Errors:</strong><ul>'
                    for err in log['errors'][:10]:
                        html += f'<li>{err[:200]}</li>'
                    html += '</ul></div>'

                if log['warnings']:
                    html += '<div class="warning-box"><strong>Warnings:</strong><ul>'
                    for warn in log['warnings'][:10]:
                        html += f'<li>{warn[:200]}</li>'
                    html += '</ul></div>'

                html += f"""
                    <details>
                        <summary>View Full Log</summary>
                        <pre class="log-content">{log['content']}</pre>
                    </details>
                </div>
"""
        html += "</div>"

    html += """
        <h2>Resource Usage Summary</h2>
        <table>
            <tr>
                <th>Process</th>
                <th>Duration</th>
                <th>CPU %</th>
                <th>Peak Memory (RSS)</th>
                <th>Peak Virtual Memory</th>
            </tr>
"""

    # Sort by duration for resource summary
    sorted_processes = sorted(processes, key=lambda x: x.get('realtime', '0'), reverse=True)[:20]
    for p in sorted_processes:
        html += f"""
            <tr>
                <td>{p['name']}</td>
                <td>{p['duration']}</td>
                <td>{p['cpu_percent']}</td>
                <td>{p['peak_rss']}</td>
                <td>{p['peak_vmem']}</td>
            </tr>
"""

    html += """
        </table>
    </div>

    <script>
        var coll = document.getElementsByClassName("collapsible");
        for (var i = 0; i < coll.length; i++) {
            coll[i].addEventListener("click", function() {
                this.classList.toggle("active");
                var content = this.nextElementSibling;
                if (content.style.display === "block") {
                    content.style.display = "none";
                } else {
                    content.style.display = "block";
                }
            });
        }
    </script>
</body>
</html>
"""

    with open(output_path, 'w') as f:
        f.write(html)

    print(f"HTML report generated: {output_path}")


def generate_csv_report(processes, log_files, output_path):
    """Generate CSV error summary report."""

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Type', 'Tool', 'rRNA', 'Process Name', 'Status', 'Exit Code',
                        'Duration', 'Error Count', 'Warning Count', 'Work Directory', 'Log File'])

        # Write process information
        for p in processes:
            writer.writerow([
                'PROCESS',
                '',
                '',
                p['name'],
                p['status'],
                p['exit'],
                p['duration'],
                '',
                '',
                p['workdir'],
                ''
            ])

        # Write log file information
        for tool, rrna_logs in log_files.items():
            for rrna, logs in rrna_logs.items():
                for log in logs:
                    writer.writerow([
                        'LOG',
                        tool,
                        rrna,
                        '',
                        'HAS_ERRORS' if log['errors'] else 'OK',
                        '',
                        '',
                        len(log['errors']),
                        len(log['warnings']),
                        '',
                        log['path']
                    ])

    print(f"CSV report generated: {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Generate error summary report for RNA modification pipeline')
    parser.add_argument('--outdir', required=True, help='Results output directory')
    parser.add_argument('--output', default='error_summary', help='Output filename prefix')
    args = parser.parse_args()

    outdir = Path(args.outdir)

    # Find trace file
    trace_files = list(outdir.glob('pipeline_info/execution_trace_*.txt'))
    processes = []
    if trace_files:
        # Use the most recent trace file
        trace_file = max(trace_files, key=os.path.getmtime)
        print(f"Parsing trace file: {trace_file}")
        processes = parse_trace_file(str(trace_file))
    else:
        print("Warning: No trace file found")

    # Find log files
    logs_dir = outdir / 'logs'
    print(f"Searching for log files in: {logs_dir}")
    log_files = find_log_files(str(logs_dir))

    # Generate reports
    report_dir = outdir / 'pipeline_info'
    report_dir.mkdir(exist_ok=True)

    html_output = report_dir / f'{args.output}.html'
    csv_output = report_dir / f'{args.output}.csv'

    generate_html_report(processes, log_files, str(html_output))
    generate_csv_report(processes, log_files, str(csv_output))

    print(f"\nSummary:")
    print(f"  Total processes: {len(processes)}")
    print(f"  Failed: {len([p for p in processes if p['status'] in ['FAILED', 'ABORTED']])}")
    print(f"  Tools with logs: {len(log_files)}")


if __name__ == '__main__':
    main()
