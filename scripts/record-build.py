#!/usr/bin/env python3
"""Record non-secret build inputs; this is provenance, not a signed attestation.

Only fixed git/.NET queries run. Installations must be operator/runner managed:
Unix ownership/modes are checked (the operator primary group is trusted);
Windows installation ACLs remain an external
runner qualification responsibility. This does not attest installation integrity.
"""
import hashlib
import json
import os
import platform
import signal
import stat
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path


root = Path(__file__).resolve().parents[1]
QUERIES = {
    ('git', 'rev-parse', 'HEAD'),
    ('git', 'status', '--porcelain', '--untracked-files=normal'),
    ('dotnet', '--version'),
}


def _within(path, parent):
    return path == parent or parent in path.parents


def _validate_installation(path):
    """Reject workspace/temp executables; don't mistake an absolute path for trust."""
    path = path.resolve(strict=True)
    forbidden = [root.resolve(), Path(tempfile.gettempdir()).resolve()]
    for variable in ('GITHUB_WORKSPACE', 'RUNNER_TEMP', 'AGENT_TEMPDIRECTORY'):
        if os.environ.get(variable):
            forbidden.append(Path(os.environ[variable]).resolve())
    if any(_within(path, location) for location in forbidden):
        raise ValueError('Build tools must be installed outside checkout and temporary storage.')
    if not path.is_file() or not os.access(path, os.X_OK):
        raise ValueError('Build tool is not an executable file.')
    if os.name != 'nt':
        for location in (path, *path.parents):
            mode = location.stat()
            permissions = stat.S_IMODE(mode.st_mode)
            if (mode.st_uid not in (0, os.getuid()) or permissions & 0o002
                    or (permissions & 0o020 and mode.st_gid not in (0, os.getgid()))):
                raise ValueError('Build tool installation is not owned and protected by the operator or root.')
    return path


def resolve_tool(name):
    """Resolve an allowlisted executable only from explicit installation roots."""
    if name not in ('git', 'dotnet'):
        raise ValueError('Unsupported build tool.')
    if os.name == 'nt':
        program_files = Path(os.environ.get('ProgramFiles', 'C:/Program Files'))
        locations = ([program_files / 'Git/cmd', program_files / 'Git/bin'] if name == 'git'
                     else [program_files / 'dotnet', Path.home() / '.dotnet'])
        executable = name + '.exe'
    else:
        locations = ([Path('/usr/bin'), Path('/usr/local/bin'), Path('/opt/homebrew/bin')]
                     if name == 'git' else [Path.home() / '.dotnet', Path('/usr/share/dotnet'),
                                           Path('/usr/local/share/dotnet'), Path('/usr/local/bin'),
                                           Path('/opt/homebrew/bin')])
        executable = name
    # setup-dotnet owns these installations; the operator must control these
    # environment variables. They cannot authorize a workspace/temp executable.
    if name == 'dotnet':
        for variable in ('DOTNET_ROOT', 'DOTNET_ROOT_X64', 'DOTNET_ROOT_ARM64'):
            if os.environ.get(variable):
                location = Path(os.environ[variable])
                if not location.is_absolute():
                    raise ValueError('The .NET installation root must be absolute.')
                locations.insert(0, location)
    for location in locations:
        candidate = location / executable
        if candidate.exists():
            resolved = _validate_installation(candidate)
            # Symlinks are allowed only into other known installation locations,
            # including Homebrew's managed Cellar behind its bin links.
            approved = [item.resolve() for item in locations]
            if os.name != 'nt':
                approved += [Path('/opt/homebrew/Cellar'), Path('/usr/local/Cellar')]
            if not any(_within(resolved, item) for item in approved):
                raise ValueError('Build tool resolves outside approved installation roots.')
            return resolved
    raise FileNotFoundError(f'No approved {name} installation found.')


def _stop(process):
    try:
        if os.name == 'nt':
            # Fixed OS adapter: never look up taskkill using the checkout or PATH.
            system_root = Path(os.environ.get('SystemRoot', 'C:/Windows'))
            taskkill = _validate_installation(system_root / 'System32/taskkill.exe')
            subprocess.run([str(taskkill), '/PID', str(process.pid), '/T', '/F'],
                           stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL, timeout=5, check=False, shell=False)
        else:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
    finally:
        # A failed OS tree-cleanup adapter must still terminate the direct child.
        if process.poll() is None:
            process.kill()
        process.wait(timeout=5)


def _run(executable, arguments, timeout=30, output_limit=1_048_576):
    """Drain both pipes concurrently with a combined byte budget and deadline."""
    if timeout <= 0 or output_limit <= 0:
        raise ValueError('Command limits must be positive.')
    process = subprocess.Popen([str(executable), *arguments], cwd=root,
                               stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, shell=False,
                               start_new_session=os.name != 'nt')
    buffers = [bytearray(), bytearray()]
    exceeded = threading.Event()
    lock = threading.Lock()
    total = 0

    def drain(stream, destination):
        nonlocal total
        try:
            while chunk := stream.read(4096):
                with lock:
                    remaining = max(0, output_limit - total)
                    buffers[destination].extend(chunk[:remaining])
                    total += len(chunk)
                    if total > output_limit:
                        exceeded.set()
        finally:
            stream.close()

    readers = [threading.Thread(target=drain, args=(process.stdout, 0), daemon=True),
               threading.Thread(target=drain, args=(process.stderr, 1), daemon=True)]
    for reader in readers:
        reader.start()
    deadline = time.monotonic() + timeout
    try:
        while process.poll() is None or any(reader.is_alive() for reader in readers):
            if exceeded.is_set():
                raise RuntimeError('Build query exceeded its output byte limit.')
            if time.monotonic() >= deadline:
                raise TimeoutError('Build query exceeded its deadline.')
            time.sleep(min(0.01, max(0, deadline - time.monotonic())))
        if exceeded.is_set():
            raise RuntimeError('Build query exceeded its output byte limit.')
        if process.returncode != 0:
            # Do not echo arbitrary tool output into logs/errors.
            raise RuntimeError(f'Build query failed with exit code {process.returncode}.')
        return bytes(buffers[0]).decode('utf-8', errors='strict').strip()
    finally:
        if process.poll() is None or any(reader.is_alive() for reader in readers):
            _stop(process)
        for reader in readers:
            reader.join(timeout=5)
        if any(reader.is_alive() for reader in readers):
            raise RuntimeError('Build query output cleanup did not finish.')


def command(*args):
    if args not in QUERIES:
        raise ValueError('Unsupported provenance query.')
    return _run(resolve_tool(args[0]), args[1:])


def main():
    sdk_version = command('dotnet', '--version')
    expected_sdk = json.loads((root / 'global.json').read_text(encoding='utf-8'))['sdk']['version']
    if sdk_version != expected_sdk:
        raise RuntimeError('Provenance SDK version does not match global.json.')
    inputs = ['global.json', 'Directory.Build.props', 'Directory.Build.targets',
              'Directory.Packages.props', 'NuGet.Config', '.editorconfig']
    inputs += sorted(path.relative_to(root).as_posix() for path in root.glob('**/packages.lock.json')
                     if not any(part in ('artifacts', '.nuget', 'bin', 'obj') for part in path.relative_to(root).parts))
    record = {
        'schemaVersion': 1,
        'sourceCommit': command('git', 'rev-parse', 'HEAD'),
        'dirty': bool(command('git', 'status', '--porcelain', '--untracked-files=normal')),
        'sdkVersion': sdk_version,
        'buildProfile': {'configuration': 'Release', 'continuousIntegrationBuild': True, 'warnAsError': True},
        'host': {'system': platform.system(), 'release': platform.release(), 'machine': platform.machine()},
        'runner': {key: os.environ.get(key) for key in ('ImageOS', 'ImageVersion', 'EXPECTED_TEST_ARCH', 'GITHUB_RUN_ID')},
        'inputsSha256': {name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in inputs},
        'cachePolicy': 'repository-local approved-feed cache; no restored cross-job cache',
    }
    sys.stdout.buffer.write((json.dumps(record, indent=2) + '\n').encode('utf-8'))


if __name__ == '__main__':
    main()
