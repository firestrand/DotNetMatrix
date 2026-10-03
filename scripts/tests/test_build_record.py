"""Exercise provenance query limits and installation trust boundaries."""
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location(
    'build_record', Path(__file__).resolve().parents[1] / 'record-build.py')
assert SPEC is not None and SPEC.loader is not None
RECORD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RECORD)


class BuildRecordTests(unittest.TestCase):
    def test_success_drains_both_streams_and_preserves_argument_boundaries(self):
        arguments = ['spaces here', '"quotes"', '', '\\trailing\\', 'Unicode Ω']
        output = RECORD._run(sys.executable, ['-c',
            'import json,sys;sys.stderr.write("diagnostic");print(json.dumps(sys.argv[1:]))', *arguments])
        self.assertEqual(arguments, json.loads(output))

    def test_nonzero_exit_does_not_echo_tool_output(self):
        with self.assertRaisesRegex(RuntimeError, 'exit code 7') as error:
            RECORD._run(sys.executable, ['-c', 'import sys;print("private text");sys.exit(7)'])
        self.assertNotIn('private text', str(error.exception))

    def test_stdout_and_stderr_share_output_budget(self):
        with self.assertRaisesRegex(RuntimeError, 'output byte limit'):
            RECORD._run(sys.executable, ['-c',
                'import sys;sys.stdout.write("x"*20000);sys.stderr.write("y"*20000)'], output_limit=1000)

    def test_timeout_terminates_child(self):
        children = []
        original = subprocess.Popen

        def capture(*args, **kwargs):
            child = original(*args, **kwargs)
            children.append(child)
            return child

        with patch.object(RECORD.subprocess, 'Popen', side_effect=capture):
            with self.assertRaisesRegex(TimeoutError, 'deadline'):
                RECORD._run(sys.executable, ['-c', 'import time;time.sleep(20)'], timeout=0.1)
        self.assertIsNotNone(children[0].poll())

    def test_unknown_tool_or_query_rejected_without_execution(self):
        with patch.object(RECORD.subprocess, 'Popen') as start:
            with self.assertRaises(ValueError):
                RECORD.resolve_tool('python')
            with self.assertRaises(ValueError):
                RECORD.command('git', 'config', '--list')
            start.assert_not_called()

    def test_temporary_sdk_installation_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            executable = Path(directory) / ('dotnet.exe' if os.name == 'nt' else 'dotnet')
            executable.write_text('placeholder', encoding='utf-8')
            executable.chmod(0o755)
            with patch.dict(os.environ, {'DOTNET_ROOT': directory}, clear=False):
                with self.assertRaisesRegex(ValueError, 'temporary storage'):
                    RECORD.resolve_tool('dotnet')

    def test_checkout_installation_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            executable = Path(directory) / 'git'
            executable.write_text('placeholder', encoding='utf-8')
            executable.chmod(0o755)
            with patch.object(RECORD, 'root', Path(directory)):
                with self.assertRaisesRegex(ValueError, 'checkout'):
                    RECORD._validate_installation(executable)

    @unittest.skipIf(os.name == 'nt', 'Unix installation-mode qualification')
    def test_world_writable_installation_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            executable = Path(directory) / 'git'
            executable.write_text('placeholder', encoding='utf-8')
            executable.chmod(0o777)
            with patch.object(RECORD.tempfile, 'gettempdir', return_value='/unrelated-temp-root'):
                with self.assertRaisesRegex(ValueError, 'owned and protected'):
                    RECORD._validate_installation(executable)

    def test_relative_sdk_root_rejected(self):
        with patch.dict(os.environ, {'DOTNET_ROOT': 'relative/install'}, clear=False):
            with self.assertRaisesRegex(ValueError, 'absolute'):
                RECORD.resolve_tool('dotnet')

    def test_manifest_rejects_sdk_version_mismatch(self):
        with patch.object(RECORD, 'command', return_value='0.0.0'):
            with self.assertRaisesRegex(RuntimeError, 'does not match global.json'):
                RECORD.main()

    def test_limits_must_be_positive(self):
        for limits in ({'timeout': 0}, {'output_limit': 0}):
            with self.assertRaises(ValueError):
                RECORD._run(sys.executable, ['-c', 'pass'], **limits)


if __name__ == '__main__':
    unittest.main()
