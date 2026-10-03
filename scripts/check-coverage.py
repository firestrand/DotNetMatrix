#!/usr/bin/env python3
"""Reject missing/partial reports and enforce strict production coverage thresholds."""
import json
import sys
from pathlib import Path
from xml.etree import ElementTree


def check(directory, expected_types, baseline=None):
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
    if (total <= 0 or int(counters.attrib['passed']) != total
            or any(int(counters.attrib.get(name, '0')) != 0
                   for name in ('failed', 'error', 'timeout', 'aborted', 'notExecuted', 'inconclusive'))):
        raise ValueError('Every discovered test must pass; empty or skipped suites are rejected.')
    if baseline is not None and total < baseline['minimumTests']:
        raise ValueError('Discovered test count is below the reviewed baseline.')
    expected = set(expected_types)
    if not expected or any(not name.startswith('DotNetMatrix.') for name in expected):
        raise ValueError('Production class inventory is empty or includes another assembly.')
    classes = root.findall('.//class')
    # Cobertura emits a separate class entry for each source file of a partial
    # type. Distinct files are valid; duplicated entries for one file are not.
    entries = [(c.attrib['name'], c.attrib.get('filename', '')) for c in classes]
    if len(entries) != len(set(entries)) or {c.attrib['name'] for c in classes} != expected:
        raise ValueError(f'Coverage class inventory mismatch. Missing: {sorted(expected - {c.attrib["name"] for c in classes})}; '
                         f'unexpected: {sorted({c.attrib["name"] for c in classes} - expected)}.')
    rates = {}
    for metric in ('lines', 'branches'):
        covered = int(root.attrib[f'{metric}-covered'])
        valid = int(root.attrib[f'{metric}-valid'])
        if valid <= 0 or covered < 0 or covered > valid or covered * 100 <= valid * 80:
            raise ValueError(f'{metric} coverage must strictly exceed 80%: {covered}/{valid}')
        rates[metric] = {'covered': covered, 'valid': valid, 'percent': covered * 100 / valid}
        if baseline is not None:
            previous = baseline[metric]
            if covered * previous['valid'] < previous['covered'] * valid:
                raise ValueError(f'{metric} coverage decreased below the reviewed baseline.')
    modules = []
    for cls in classes:
        modules.append({
            'name': cls.attrib['name'],
            'source_file': cls.attrib.get('filename', ''),
            'line_percent': float(cls.attrib['line-rate']) * 100,
            'branch_percent': float(cls.attrib['branch-rate']) * 100,
        })
    summary = {'tests_passed': total, 'overall': rates, 'modules': modules}
    (directory / 'coverage-summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(f"{total} tests passed; line coverage {rates['lines']['percent']:.2f}%; "
          f"branch coverage {rates['branches']['percent']:.2f}%.")
    for module in modules:
        print(f"  {module['name']}: lines {module['line_percent']:.2f}%, "
              f"branches {module['branch_percent']:.2f}%")
    print(f'Reports: {directory}')


if __name__ == '__main__':
    if len(sys.argv) not in (3, 4):
        sys.exit('Usage: python3 scripts/check-coverage.py RESULTS_DIRECTORY PRODUCTION_TYPES_JSON [BASELINE_JSON]')
    try:
        baseline = json.loads(Path(sys.argv[3]).read_text(encoding='utf-8')) if len(sys.argv) == 4 else None
        check(Path(sys.argv[1]), json.loads(Path(sys.argv[2]).read_text(encoding='utf-8')), baseline)
    except (ValueError, OSError, ElementTree.ParseError, KeyError) as error:
        sys.exit(f'Coverage gate failed: {error}')
