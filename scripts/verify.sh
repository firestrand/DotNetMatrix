#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p artifacts
results_dir="$(mktemp -d artifacts/verification.XXXXXX)"
dotnet restore DotNetMatrix.sln --locked-mode
dotnet build DotNetMatrix.sln -c Release --no-restore
dotnet format DotNetMatrix.sln --verify-no-changes --no-restore
dotnet test --project DotNetMatrix_Test/DotNetMatrix_Test.csproj -c Release --no-build \
  --results-directory "$results_dir" -- \
  --report-trx --coverlet --coverlet-include '[DotNetMatrix]*' \
  --coverlet-output-format cobertura --coverlet-threshold 80 \
  --coverlet-threshold-type line branch --coverlet-threshold-stat Total
python3 scripts/check-coverage.py "$results_dir"
