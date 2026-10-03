"""Exercise actual MSBuild warning policy using isolated, package-free fixtures."""
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


class BuildPolicyTests(unittest.TestCase):
    def setUp(self):
        (ROOT / 'artifacts').mkdir(exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(prefix='warning-policy.', dir=ROOT / 'artifacts')
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        (self.directory / 'Probe.csproj').write_text('<Project Sdk="Microsoft.NET.Sdk" />\n', encoding='utf-8', newline='\n')
        result = self.command('restore', 'Probe.csproj', '-p:Configuration=Release')
        self.assertEqual(0, result.returncode, result.stdout)

    def command(self, *args):
        return subprocess.run(['dotnet', *args], cwd=self.directory, check=False,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True, timeout=90)

    def build(self, source, *properties):
        (self.directory / 'Probe.cs').write_text(source, encoding='utf-8', newline='\n')
        return self.command('build', 'Probe.csproj', '-c', 'Release', '--no-restore', *properties)

    def test_local_relaxation_preserves_nullable_and_platform_errors(self):
        advisory = '#warning ADVISORY\nnamespace BuildProbe;\ninternal static class Probe {}\n'
        strict = self.build(advisory)
        self.assertNotEqual(0, strict.returncode)
        self.assertIn('error CS1030', strict.stdout)
        relaxed = self.build(advisory, '-p:LocalWarningRelaxation=true', '-p:ContinuousIntegrationBuild=false')
        self.assertEqual(0, relaxed.returncode, relaxed.stdout)
        self.assertIn('warning CS1030', relaxed.stdout)
        for source, diagnostic in (
            ('namespace BuildProbe;\ninternal static class Probe { internal static string Read() => null; }\n', 'CS8603'),
            ('using System.Runtime.Versioning;\nnamespace BuildProbe;\ninternal static class Probe {\n'
             ' [SupportedOSPlatform("windows")] internal static void WindowsOnly() {}\n'
             ' internal static void Call() => WindowsOnly();\n}\n', 'CA1416'),
        ):
            with self.subTest(diagnostic=diagnostic):
                result = self.build(source, '-p:LocalWarningRelaxation=true', '-p:ContinuousIntegrationBuild=false')
                self.assertNotEqual(0, result.returncode, result.stdout)
                self.assertIn('error ' + diagnostic, result.stdout)

    def test_ci_rejects_relaxation_and_analyzer_disabling(self):
        for property_value in ('LocalWarningRelaxation=true', 'TreatWarningsAsErrors=false',
                               'RunAnalyzers=false', 'EnforceCodeStyleInBuild=false'):
            with self.subTest(property_value=property_value):
                result = self.build('namespace BuildProbe;\ninternal static class Probe {}\n',
                                    '-p:ContinuousIntegrationBuild=true', '-p:' + property_value)
                self.assertNotEqual(0, result.returncode, result.stdout)
                self.assertIn('error', result.stdout)
