"""Exercise verifier failures independently of real coverage totals."""
import importlib.util
import tempfile
import unittest
from pathlib import Path
from xml.etree import ElementTree as ET


def load(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).parents[1] / f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


coverage = load('check-coverage')
api = load('check-api')
vulnerabilities = load('check-vulnerabilities')


class CoverageGateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.root = ET.Element('coverage', {
            'lines-covered': '81', 'lines-valid': '100',
            'branches-covered': '81', 'branches-valid': '100',
        })
        ET.SubElement(self.root, 'class', {'name': 'DotNetMatrix.Old', 'line-rate': '.81', 'branch-rate': '.81'})
        ET.SubElement(self.root, 'class', {'name': 'DotNetMatrix.New', 'line-rate': '.81', 'branch-rate': '.81'})
        self.counters = ET.Element('Counters', {'total': '2', 'passed': '2'})
        self.write()

    def write(self):
        ET.ElementTree(self.root).write(self.directory / 'coverage.cobertura.xml')
        result = ET.Element('TestRun')
        result.append(self.counters)
        ET.ElementTree(result).write(self.directory / 'result.trx')

    def check(self):
        coverage.check(self.directory, ['DotNetMatrix.Old', 'DotNetMatrix.New'])

    def test_above_threshold_with_new_class_passes(self):
        self.check()
        self.assertTrue((self.directory / 'coverage-summary.json').exists())

    def test_partial_class_in_distinct_source_files_passes(self):
        self.root.findall('class')[0].set('filename', 'Old.cs')
        ET.SubElement(self.root, 'class', {'name': 'DotNetMatrix.Old', 'filename': 'Old.Partial.cs',
                                         'line-rate': '.81', 'branch-rate': '.81'})
        self.write()
        self.check()
        ET.SubElement(self.root, 'class', {'name': 'DotNetMatrix.Old', 'filename': 'Old.Partial.cs'})
        self.write()
        with self.assertRaisesRegex(ValueError, 'inventory mismatch'):
            self.check()

    def test_exact_80_percent_rejected_for_both_metrics(self):
        for metric in ('lines', 'branches'):
            with self.subTest(metric=metric):
                self.root.set(f'{metric}-covered', '80')
                self.write()
                with self.assertRaisesRegex(ValueError, 'strictly exceed'):
                    self.check()
                self.root.set(f'{metric}-covered', '81')

    def test_missing_new_production_module_rejected(self):
        self.root.remove(self.root.findall('class')[1])
        self.write()
        with self.assertRaisesRegex(ValueError, 'Missing.*DotNetMatrix.New'):
            self.check()

    def test_test_assembly_and_duplicate_modules_rejected(self):
        for name in ('TestAssembly.Test', 'DotNetMatrix.New'):
            with self.subTest(name=name):
                extra = ET.SubElement(self.root, 'class', {'name': name})
                self.write()
                with self.assertRaisesRegex(ValueError, 'inventory mismatch'):
                    self.check()
                self.root.remove(extra)

    def test_missing_and_duplicate_reports_rejected(self):
        for filename in ('result.trx', 'coverage.cobertura.xml'):
            with self.subTest(filename=filename):
                path = self.directory / filename
                path.unlink()
                with self.assertRaisesRegex(ValueError, 'exactly one'):
                    self.check()
                self.write()
        (self.directory / 'second.trx').write_text('<TestRun/>')
        with self.assertRaisesRegex(ValueError, 'exactly one'):
            self.check()

    def test_failed_skipped_and_empty_suites_rejected(self):
        for counters in ({'total': '2', 'passed': '1', 'failed': '1'},
                         {'total': '2', 'passed': '1', 'notExecuted': '1'},
                         {'total': '0', 'passed': '0'},
                         {'total': '2', 'passed': '2', 'error': '1'}):
            with self.subTest(counters=counters):
                self.counters.attrib = counters
                self.write()
                with self.assertRaisesRegex(ValueError, 'Every discovered test'):
                    self.check()

    def test_invalid_counts_rejected(self):
        for covered, valid in (('-1', '100'), ('101', '100'), ('0', '0')):
            with self.subTest(covered=covered, valid=valid):
                self.root.set('lines-covered', covered)
                self.root.set('lines-valid', valid)
                self.write()
                with self.assertRaises(ValueError):
                    self.check()

    def test_reviewed_baseline_rejects_reduced_tests_or_coverage(self):
        baseline = {'minimumTests': 2, 'lines': {'covered': 81, 'valid': 100},
                    'branches': {'covered': 81, 'valid': 100}}
        coverage.check(self.directory, ['DotNetMatrix.Old', 'DotNetMatrix.New'], baseline)
        baseline['minimumTests'] = 3
        with self.assertRaisesRegex(ValueError, 'test count'):
            coverage.check(self.directory, ['DotNetMatrix.Old', 'DotNetMatrix.New'], baseline)
        baseline['minimumTests'] = 2
        for metric in ('lines', 'branches'):
            baseline[metric]['covered'] = 82
            with self.assertRaisesRegex(ValueError, 'decreased'):
                coverage.check(self.directory, ['DotNetMatrix.Old', 'DotNetMatrix.New'], baseline)
            baseline[metric]['covered'] = 81


class VulnerabilityGateTests(unittest.TestCase):
    def setUp(self):
        self.root = Path.cwd()
        self.report = {'version': 1, 'parameters': '--vulnerable --include-transitive',
                       'sources': ['https://data.nuget.org/v3/index.json'],
                       'projects': [{'path': str(self.root / name),
                                     'frameworks': [{'framework': 'net10.0', 'topLevelPackages': []}]}
                                    for name in sorted(vulnerabilities.PROJECTS)]}

    def test_complete_empty_result_passes(self):
        self.assertEqual([], vulnerabilities.check(self.report, self.root))

    def test_unknown_or_serious_transitive_vulnerability_blocks(self):
        finding = {'severity': 'low'}
        self.report['projects'][0]['frameworks'] = [{'framework': 'net10.0', 'topLevelPackages': [], 'transitivePackages': [
            {'id': 'Example.Package', 'vulnerabilities': [finding]}]}]
        self.assertEqual([('Example.Package', 'low')], vulnerabilities.check(self.report, self.root))
        for severity in ('high', 'critical', 'unknown', ''):
            finding['severity'] = severity
            with self.subTest(severity=severity), self.assertRaises(ValueError):
                vulnerabilities.check(self.report, self.root)

    def test_incomplete_failed_and_untrusted_evidence_rejected(self):
        for field, value in (('projects', []), ('errors', ['audit unavailable']),
                             ('sources', ['https://example.com/feed']), ('parameters', '--vulnerable')):
            with self.subTest(field=field), self.assertRaises(ValueError):
                vulnerabilities.check(dict(self.report, **{field: value}), self.root)
        self.report['projects'].append(self.report['projects'][0])
        with self.assertRaisesRegex(ValueError, 'exactly once'):
            vulnerabilities.check(self.report, self.root)

    def test_missing_framework_or_vulnerability_evidence_fails_closed(self):
        project = self.report['projects'][0]
        for frameworks in (None, [], [{'framework': 'net9.0', 'topLevelPackages': []}],
                           [{'framework': 'net10.0'}],
                           [{'framework': 'net10.0', 'topLevelPackages': [{'id': 'Unknown.Package'}]}],
                           [{'framework': 'net10.0', 'topLevelPackages': [], 'problems': ['audit failed']}]):
            with self.subTest(frameworks=frameworks), self.assertRaises(ValueError):
                project['frameworks'] = frameworks
                vulnerabilities.check(self.report, self.root)


class ApiGateTests(unittest.TestCase):
    def test_matching_api_passes_and_removed_member_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            baseline = Path(directory) / 'baseline'
            actual = Path(directory) / 'actual'
            baseline.write_text('type A\nmethod A.Solve()\n')
            actual.write_text(baseline.read_text())
            api.check(baseline, actual)
            actual.write_text('type A\n')
            with self.assertRaisesRegex(ValueError, 'Public API changed'):
                api.check(baseline, actual)


if __name__ == '__main__':
    unittest.main()
