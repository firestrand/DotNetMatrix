#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
# A changed feed/mapping policy gets an empty namespace rather than trusting an
# existing global cache. Never import arbitrary developer or cross-job caches.
export NUGET_PACKAGES="$(python3 -c 'import hashlib, pathlib; root = pathlib.Path.cwd(); fingerprint = hashlib.sha256((root / "NuGet.Config").read_bytes()).hexdigest(); print((root / ".nuget" / "packages" / fingerprint).as_posix())')"
mkdir -p artifacts
results_dir="$(mktemp -d artifacts/verification.XXXXXX)"
dotnet restore DotNetMatrix.sln --locked-mode -p:Configuration=Release -p:ContinuousIntegrationBuild=true -warnaserror
dotnet build DotNetMatrix.sln -c Release --no-restore -p:ContinuousIntegrationBuild=true -warnaserror
dotnet format DotNetMatrix.sln --verify-no-changes --no-restore
dotnet format style DotNetMatrix.sln --verify-no-changes --no-restore --severity info --diagnostics IDE0005
dotnet package list --project DotNetMatrix.sln --vulnerable --include-transitive --format json --no-restore > "$results_dir/vulnerabilities.json"
python3 scripts/check-vulnerabilities.py "$results_dir/vulnerabilities.json"
python3 -m unittest discover -s scripts/tests -v
dotnet restore tools/ApiBaseline/ApiBaseline.csproj --locked-mode -p:Configuration=Release -p:ContinuousIntegrationBuild=true -warnaserror
dotnet build tools/ApiBaseline/ApiBaseline.csproj -c Release --no-restore -p:ContinuousIntegrationBuild=true -warnaserror
dotnet tools/ApiBaseline/bin/Release/net10.0/ApiBaseline.dll api \
  DotNetMatrix/bin/Release/net10.0/DotNetMatrix.dll > "$results_dir/public-api.txt"
python3 scripts/check-api.py docs/public-api.txt "$results_dir/public-api.txt"
dotnet tools/ApiBaseline/bin/Release/net10.0/ApiBaseline.dll coverage-types \
  DotNetMatrix/bin/Release/net10.0/DotNetMatrix.dll > "$results_dir/production-types.json"
dotnet test --project DotNetMatrix_Test/DotNetMatrix_Test.csproj -c Release --no-build \
  -p:ContinuousIntegrationBuild=true \
  --results-directory "$results_dir" -- \
  --report-trx --coverlet --coverlet-include '[DotNetMatrix]*' \
  --coverlet-output-format cobertura --coverlet-threshold 80 \
  --coverlet-threshold-type line branch --coverlet-threshold-stat Total
python3 scripts/check-coverage.py "$results_dir" "$results_dir/production-types.json" docs/coverage-baseline.json
bash scripts/package-smoke.sh
dotnet run --project benchmarks/DotNetMatrix.Benchmarks.csproj -c Release --no-build --no-restore -- --validate
dotnet run --project samples/LeastSquares/LeastSquares.csproj -c Release --no-build --no-restore
# A clean hosted checkout must retrieve exact committed production sources.
if [[ "${CI:-}" == "true" ]]; then
  dotnet tools/ApiBaseline/bin/Release/net10.0/ApiBaseline.dll source-link \
    DotNetMatrix/bin/Release/net10.0/DotNetMatrix.dll > "$results_dir/source-link.json"
  python3 scripts/check-source-link.py "$results_dir/source-link.json" "$(git rev-parse HEAD)"
fi
# Unapproved exceptions deliberately fail the complete gate.
python3 scripts/check-standards.py
