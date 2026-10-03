#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p artifacts
smoke_dir="$(mktemp -d artifacts/package-smoke.XXXXXX)"
smoke_dir="$(python3 -c 'import pathlib, sys; print(pathlib.Path(sys.argv[1]).resolve().as_posix())' "$smoke_dir")"
dotnet pack DotNetMatrix/DotNetMatrix.csproj -c Release --no-build --no-restore -p:ContinuousIntegrationBuild=true -warnaserror -o "$smoke_dir/feed"
mkdir -p "$smoke_dir/consumer"
cat > "$smoke_dir/consumer/Consumer.csproj" <<'PROJECT'
<Project Sdk="Microsoft.NET.Sdk">
  <PropertyGroup>
    <OutputType>Exe</OutputType>
    <TargetFramework>net10.0</TargetFramework>
    <ImplicitUsings>disable</ImplicitUsings>
    <Nullable>enable</Nullable>
    <TreatWarningsAsErrors>true</TreatWarningsAsErrors>
    <RestorePackagesWithLockFile>true</RestorePackagesWithLockFile>
  </PropertyGroup>
  <ItemGroup>
    <PackageReference Include="DotNetMatrix.LocalPreview" Version="2.0.0-preview.1" />
  </ItemGroup>
</Project>
PROJECT
cat > "$smoke_dir/consumer/Program.cs" <<'PROGRAM'
using System;
using System.Text.Json;
using DotNetMatrix;

GeneralMatrix matrix = new(new double[][] { new[] { 4.0, 2.0 }, new[] { 2.0, 10.0 } });
GeneralMatrix rhs = new(new double[][] { new[] { 10.0 }, new[] { 32.0 } });
GeneralMatrix solution = matrix.Chol().Solve(rhs);
if (Math.Abs(solution.Array[0][0] - 1) > 1e-12 || Math.Abs(solution.Array[1][0] - 3) > 1e-12)
    throw new InvalidOperationException("Packaged Cholesky contract failed.");
GeneralMatrix row = new(new double[][] { new[] { 1.0, 1.0 } });
GeneralMatrix minimumNorm = row.SolveMinimumNorm(new GeneralMatrix(new double[][] { new[] { 2.0 } }));
if (minimumNorm.RowDimension != 2 || Math.Abs(minimumNorm.Array[0][0] - 1) > 1e-12 || Math.Abs(minimumNorm.Array[1][0] - 1) > 1e-12)
    throw new InvalidOperationException("Packaged minimum-norm contract failed.");
JsonSerializerOptions jsonOptions = new();
jsonOptions.Converters.Add(new MatrixJsonConverter());
string payload = JsonSerializer.Serialize(matrix, jsonOptions);
GeneralMatrix? restored = JsonSerializer.Deserialize<GeneralMatrix>(payload, jsonOptions);
if (restored == null || !matrix.Equals(restored))
    throw new InvalidOperationException("Packaged JSON contract failed.");
Console.WriteLine("Independent local package consumer passed.");
PROGRAM
cat > "$smoke_dir/consumer/NuGet.Config" <<CONFIG
<configuration><packageSources><clear/><add key="local-preview" value="$smoke_dir/feed"/></packageSources></configuration>
CONFIG
# The generated consumer has no project reference and restores into a fresh cache.
# Run outside this repository so central properties/packages cannot influence it.
consumer_dir="$(python3 -c 'import pathlib, tempfile; print(pathlib.Path(tempfile.mkdtemp(prefix="dotnetmatrix-consumer.")).as_posix())')"
trap 'python3 -c "import shutil, sys; shutil.rmtree(sys.argv[1])" "$consumer_dir"' EXIT
cp "$smoke_dir/consumer/"* "$consumer_dir/"
cp global.json "$consumer_dir/"
dotnet restore "$consumer_dir/Consumer.csproj" --configfile "$consumer_dir/NuGet.Config" --packages "$smoke_dir/cache"
cp "$consumer_dir/packages.lock.json" "$smoke_dir/consumer/"
dotnet restore "$consumer_dir/Consumer.csproj" --configfile "$consumer_dir/NuGet.Config" --packages "$smoke_dir/cache" --locked-mode
dotnet run --project "$consumer_dir/Consumer.csproj" -c Release --no-restore
python3 - "$smoke_dir/feed/DotNetMatrix.LocalPreview.2.0.0-preview.1.nupkg" <<'PY'
import subprocess, sys, zipfile
from pathlib import Path
from xml.etree import ElementTree
with zipfile.ZipFile(sys.argv[1]) as package:
    names = set(package.namelist())
    for expected in ('lib/net10.0/DotNetMatrix.dll', 'lib/net10.0/DotNetMatrix.xml',
                     'README.md', 'LICENSE', 'docs/provenance-and-release.md'):
        if expected not in names:
            sys.exit(f'Package is missing {expected}')
    if package.read('LICENSE') != Path('LICENSE').read_bytes():
        sys.exit('Packaged license differs from the repository license.')
    nuspec = ElementTree.fromstring(package.read('DotNetMatrix.LocalPreview.nuspec'))
    license_element = nuspec.find('.//{*}license')
    if license_element is None or license_element.attrib.get('type') != 'expression' or license_element.text != 'MIT':
        sys.exit('Package must declare the MIT license expression.')
    repository = nuspec.find('.//{*}repository')
    commit = subprocess.run(['git', 'rev-parse', 'HEAD'], check=True, capture_output=True,
                            text=True, timeout=30).stdout.strip()
    if (repository is None or repository.attrib.get('commit') != commit
            or repository.attrib.get('url') != 'https://github.com/firestrand/DotNetMatrix'):
        sys.exit('Package repository metadata does not identify the exact source commit.')
with zipfile.ZipFile(Path(sys.argv[1]).with_suffix('.snupkg')) as symbols:
    if symbols.read('lib/net10.0/DotNetMatrix.pdb') != Path('DotNetMatrix/bin/Release/net10.0/DotNetMatrix.pdb').read_bytes():
        sys.exit('Symbol package does not contain the matching portable PDB.')
print('Package DLL, XML documentation, MIT license, commit metadata and matching symbols verified.')
PY
