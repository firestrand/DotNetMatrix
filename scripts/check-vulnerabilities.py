#!/usr/bin/env python3
"""Validate complete NuGet audit output; unknown/high/critical findings block."""
import json
import sys
from pathlib import Path


PROJECTS = {
    'DotNetMatrix/DotNetMatrix.csproj',
    'DotNetMatrix_Test/DotNetMatrix_Test.csproj',
    'tools/ApiBaseline/ApiBaseline.csproj',
    'samples/LeastSquares/LeastSquares.csproj',
    'benchmarks/DotNetMatrix.Benchmarks.csproj',
}


def check(report, root):
    if (report.get('version') != 1 or '--vulnerable' not in report.get('parameters', '')
            or '--include-transitive' not in report.get('parameters', '')
            or report.get('sources') != ['https://data.nuget.org/v3/index.json']
            or report.get('errors') or report.get('problems') or report.get('warnings')):
        raise ValueError('Missing or unsupported complete NuGet audit evidence.')
    projects = report.get('projects', [])
    paths = [Path(project['path']).resolve().relative_to(root.resolve()).as_posix() for project in projects]
    if len(paths) != len(set(paths)) or set(paths) != PROJECTS:
        raise ValueError('Audit must include every solution project exactly once.')
    findings = []
    for project in projects:
        if project.get('errors') or project.get('problems') or project.get('warnings'):
            raise ValueError('NuGet audit reported unavailable or incomplete evidence.')
        frameworks = project.get('frameworks')
        if (not isinstance(frameworks, list) or len(frameworks) != 1
                or frameworks[0].get('framework') != 'net10.0'
                or not isinstance(frameworks[0].get('topLevelPackages'), list)):
            raise ValueError('Audit evidence must explicitly include the net10.0 framework for every project.')
        for framework in frameworks:
            if framework.get('errors') or framework.get('problems') or framework.get('warnings'):
                raise ValueError('Framework audit reported unavailable or incomplete evidence.')
            for group in ('topLevelPackages', 'transitivePackages'):
                if not isinstance(framework.get(group, []), list):
                    raise ValueError('Malformed package audit group.')
                for package in framework.get(group, []):
                    vulnerabilities = package.get('vulnerabilities')
                    if not isinstance(vulnerabilities, list) or not vulnerabilities:
                        raise ValueError('A vulnerable-package result lacks vulnerability evidence.')
                    for vulnerability in vulnerabilities:
                        severity = str(vulnerability.get('severity', '')).lower()
                        if severity not in ('low', 'moderate', 'high', 'critical'):
                            raise ValueError('Unknown vulnerability severity; audit fails closed.')
                        findings.append((package['id'], severity))
                        if severity in ('high', 'critical'):
                            raise ValueError(f"{package['id']}: {severity} vulnerability blocks verification.")
    return findings


if __name__ == '__main__':
    try:
        if len(sys.argv) != 2:
            raise ValueError('Usage: python3 scripts/check-vulnerabilities.py REPORT_JSON')
        findings = check(json.loads(Path(sys.argv[1]).read_text(encoding='utf-8')), Path(__file__).resolve().parents[1])
        print(f'Fresh NuGet audit: all {len(PROJECTS)} projects; {len(findings)} low/moderate findings; no high/critical findings.')
    except (ValueError, OSError, KeyError, TypeError) as error:
        sys.exit(f'Vulnerability gate failed: {error}')
