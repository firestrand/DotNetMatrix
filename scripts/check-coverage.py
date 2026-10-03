#!/usr/bin/env python3
"""Reject missing/partial reports and enforce strict production coverage thresholds."""
import json
import sys
from pathlib import Path
from xml.etree import ElementTree


def check(directory):
    reports = list(directory.glob('coverage.cobertura*.xml'))
    results = list(directory.glob('*.trx'))
    if len(reports) != 1 or len(results) != 1:
        raise ValueError('Expected exactly one coverage report and one test result in this run.')
    root = ElementTree.parse(reports[0]).getroot()
    result = ElementTree.parse(results[0]).getroot()
    counters = result.find('.//{*}Counters')
    if counters is None:
        raise ValueError('Test report has no counters.')
    total = int(counters.attrib['total'])
    if total <= 0 or int(counters.attrib['passed']) != total:
        raise ValueError('Every discovered test must pass; empty or skipped suites are rejected.')
    expected = {
        'DotNetMatrix.GeneralMatrix', 'DotNetMatrix.CholeskyDecomposition',
        'DotNetMatrix.LUDecomposition', 'DotNetMatrix.QRDecomposition',
        'DotNetMatrix.SingularValueDecomposition', 'DotNetMatrix.EigenvalueDecomposition',
        'DotNetMatrix.Maths',
    }
    classes = root.findall('.//class')
    if len(classes) != len(expected) or {c.attrib['name'] for c in classes} != expected:
        raise ValueError('Coverage must include all seven production classes and no test assembly.')
    rates = {}
    for metric in ('lines', 'branches'):
        covered = int(root.attrib[f'{metric}-covered'])
        valid = int(root.attrib[f'{metric}-valid'])
        if valid <= 0 or covered < 0 or covered > valid or covered * 100 <= valid * 80:
            raise ValueError(f'{metric} coverage must strictly exceed 80%: {covered}/{valid}')
        rates[metric] = {'covered': covered, 'valid': valid, 'percent': covered * 100 / valid}
    modules = []
    for cls in classes:
        modules.append({
            'name': cls.attrib['name'],
            'line_percent': float(cls.attrib['line-rate']) * 100,
            'branch_percent': float(cls.attrib['branch-rate']) * 100,
        })
    summary = {'tests_passed': total, 'overall': rates, 'modules': modules}
    (directory / 'coverage-summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f"{total} tests passed; line coverage {rates['lines']['percent']:.2f}%; "
          f"branch coverage {rates['branches']['percent']:.2f}%.")
    for module in modules:
        print(f"  {module['name']}: lines {module['line_percent']:.2f}%, "
              f"branches {module['branch_percent']:.2f}%")
    print(f'Reports: {directory}')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit('Usage: python3 scripts/check-coverage.py RESULTS_DIRECTORY')
    try:
        check(Path(sys.argv[1]))
    except (ValueError, OSError, ElementTree.ParseError, KeyError) as error:
        sys.exit(f'Coverage gate failed: {error}')
