"""Check canonical public API bytes under materially different process cultures."""
import os
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


class ApiCultureTests(unittest.TestCase):
    def test_public_api_is_invariant_and_uses_lf(self):
        expected = (ROOT / 'docs/public-api.txt').read_bytes()
        self.assertNotIn(b'\r', expected)
        for culture in ('en_US.UTF-8', 'fr_FR.UTF-8', 'tr_TR.UTF-8', 'ar_SA.UTF-8'):
            with self.subTest(culture=culture):
                environment = dict(os.environ, LANG=culture, LC_ALL=culture)
                result = subprocess.run([
                    'dotnet', 'tools/ApiBaseline/bin/Release/net10.0/ApiBaseline.dll',
                    'api', 'DotNetMatrix/bin/Release/net10.0/DotNetMatrix.dll',
                ], cwd=ROOT, env=environment, check=False, capture_output=True, timeout=30)
                self.assertEqual(0, result.returncode, result.stderr.decode(errors='replace'))
                self.assertEqual(expected, result.stdout)
