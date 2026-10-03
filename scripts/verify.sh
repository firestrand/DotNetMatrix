#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p artifacts
results_dir="$(mktemp -d artifacts/verification.XXXXXX)"
dotnet restore DotNetMatrix.sln --locked-mode
dotnet build DotNetMatrix.sln -c Release --no-restore
dotnet format DotNetMatrix.sln --verify-no-changes --no-restore
python3 -m unittest discover -s scripts/tests -v
dotnet restore tools/ApiBaseline/ApiBaseline.csproj --locked-mode
dotnet build tools/ApiBaseline/ApiBaseline.csproj -c Release --no-restore
dotnet tools/ApiBaseline/bin/Release/net10.0/ApiBaseline.dll api \
  DotNetMatrix/bin/Release/net10.0/DotNetMatrix.dll > "$results_dir/public-api.txt"
python3 scripts/check-api.py docs/public-api.txt "$results_dir/public-api.txt"
dotnet tools/ApiBaseline/bin/Release/net10.0/ApiBaseline.dll coverage-types \
  DotNetMatrix/bin/Release/net10.0/DotNetMatrix.dll > "$results_dir/production-types.json"
dotnet test --project DotNetMatrix_Test/DotNetMatrix_Test.csproj -c Release --no-build \
  --results-directory "$results_dir" -- \
  --report-trx --coverlet --coverlet-include '[DotNetMatrix]*' \
  --coverlet-output-format cobertura --coverlet-threshold 80 \
  --coverlet-threshold-type line branch --coverlet-threshold-stat Total
python3 scripts/check-coverage.py "$results_dir" "$results_dir/production-types.json"
bash scripts/package-smoke.sh
dotnet run --project benchmarks/DotNetMatrix.Benchmarks.csproj -c Release --no-build --no-restore -- --validate
dotnet run --project samples/LeastSquares/LeastSquares.csproj -c Release --no-build --no-restore
