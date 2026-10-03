This code forked from:
http://www.codeproject.com/KB/recipes/psdotnetmatrix.aspx
A C# port of public domain Java Matrix library JAMA

Modern build and tests (2026-10-02)
----------------------------------
Prerequisites: .NET SDK 10.0.401 (runtime 10.0.12, C# 14), Bash, and Python 3
(standard library only, for checking coverage). The SDK is pinned in global.json.
Get the SDK from https://dotnet.microsoft.com/download/dotnet/10.0.

Run these commands from a fresh checkout:

  dotnet restore DotNetMatrix.sln --locked-mode
  dotnet build DotNetMatrix.sln -c Release --no-restore
  bash scripts/verify.sh

The verification command performs locked restore, Release build, formatting/
analyzer checks, all unit tests and production line/branch coverage. Both overall
coverage rates must strictly exceed 80%. It rejects absent/partial coverage and
empty or skipped test suites. Each run writes TRX, Cobertura and JSON summaries
under a new artifacts/verification.* directory, so stale reports cannot pass.
To execute the suite without coverage:

  dotnet test --project DotNetMatrix_Test/DotNetMatrix_Test.csproj -c Release

To collect coverage directly after building:

  dotnet test --project DotNetMatrix_Test/DotNetMatrix_Test.csproj -c Release --no-build --results-directory artifacts/manual -- --report-trx --coverlet --coverlet-include '[DotNetMatrix]*' --coverlet-output-format cobertura

Use a fresh results directory for each direct coverage run. The old VSTest
--collect syntax is replaced by the native Microsoft.Testing.Platform runner.
Dependencies are development-only MSTest 4.4.1, Coverlet MTP 10.1.0 and Microsoft
Testing Platform 2.4.1. Directory.Packages.props pins current stable direct and
transitive versions; both packages.lock.json files must be committed. There are
no third-party runtime dependencies in the library. The unused test telemetry
hook and assembly are excluded because the upstream extension expects
ApplicationInsights 2.x; its dependency graph is still restored and upgraded,
but its incompatible extension is never registered or loaded.

Verified result: 100/100 passing tests (including all 62 original unit tests and
the original numerical console harness), 99.62% line / 96.62% branch coverage.
Every production module exceeds 80% on both metrics. See README.md for current
verification commands, API documentation and compatibility guidance.

Compatibility and preserved legacy limitations
---------------------------------------------
The library now targets net10.0. Consumers must use .NET 10; it is no longer a
.NET Framework 4.0 binary. Assembly identity, public numerical API signatures,
virtual members, storage aliasing, exceptions and serialized field layout are
preserved. BinaryFormatter was removed by modern .NET; this migration does not
provide an unsafe compatibility formatter or claim old binary payload replay.
The legacy GeneralMatrix ISerializable callback already emitted no entries and
had no deserialization constructor; the empty callback remains unchanged.

Existing numerical/API defects remain explicitly characterized under the
requested no-behavior-change rule:
- Cholesky Solve is incorrect for some nonunit diagonal factors because its
  forward substitution divides too late. Reconstruction is unaffected.
- SolveTranspose returns the transpose of the documented solution.
- Equality throws when the left operand is null; hashes use backing-array
  identity despite value equality.
- Wide decompositions are not uniformly supported; LU may throw during
  construction and SVD has inconsistent wide factor dimensions.

These findings need a separate behavior-change decision. New tests cover both
correct supported equations and explicit legacy characterizations. Numerical
algorithm operation order was retained to avoid floating-point regressions.
The pointless matrix finalizer was removed, IDisposable retained, and random
factories use Random.Shared while preserving their range contracts.
