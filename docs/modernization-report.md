# DotNetMatrix modernization report

Verified 2026-10-02 on Ubuntu 24.04, Linux arm64. Verification completed on `modernization/dotnet10`. After a successful review,
the user authorized committing and pushing the reviewed changes. No merge or
release publication was requested.

## Outcome

The library builds on .NET 10 with zero warnings and zero errors. All 100 tests
pass, including all 62 original unit tests and the original numerical console
harness. Production coverage is **99.62% lines and 96.62% branches**. All seven
modules individually exceed 80% on both measures. The full verification gate
also passed in a fresh local clone with no build outputs and an empty NuGet
cache. Its library and test DLL/PDB files match the working build byte for byte.

## Baseline and upgrade strategy

Base revision: `683be7c9ad4ac33ef63a361485558681036cb14d`.
The only pre-existing untracked state was `.serena/`, preserved outside this
change. No staged or tracked edits existed.

| Area | Legacy baseline | Final configuration |
| --- | --- | --- |
| Runtime/language | .NET Framework 4.0; VS2010/MSBuild 4-era C# defaults | .NET 10.0.12; SDK 10.0.401; C# 14 |
| Build | non-SDK XML projects, manual file/assembly references | SDK projects, shared properties, deterministic path mapping |
| Dependencies | no runtime NuGet dependencies; VS2010 test DLL | no library runtime dependencies; current stable test graph, central pins and lockfiles |
| Tests | MSTest v1 DLL; 62 unit tests; standalone console harness | MSTest 4.4.1, native Microsoft.Testing.Platform 2.4.1; 100 cases |
| Coverage | none configured or measured | Coverlet MTP 10.1.0; line/branch Cobertura, TRX, JSON and strict gate |
| Compiler/lint/format | warning level 4; no linter/formatter config | nullable enabled, warnings as errors, current built-in analyzers, EditorConfig and dotnet format |
| CI/containers | none present | no existing CI/container migration required; local verification script supplied |

Original `dotnet build DotNetMatrix.sln` failed with exit 1 / MSB3644: .NET
Framework 4.0 reference assemblies are unavailable on this Linux host.
Original `dotnet test DotNetMatrix.sln` returned exit 0 but discovered/executed
no tests. That was **not** a passing baseline suite. Initial coverage was
unavailable, not 0%. Logs remain in `artifacts/baseline/`.

A direct jump was appropriate: the numerical library uses portable BCL arrays,
arithmetic and exceptions, without third-party application frameworks or
runtime packages requiring bridge migrations. Installed tooling and official
[.NET downloads](https://dotnet.microsoft.com/en-us/download/dotnet/10.0) were
verified; stable package versions were queried from NuGet's version index on
2026-10-02. `dotnet list ... package --outdated --include-transitive` reports no
updates for either project; the matching vulnerability audit reports none.
Versions are intentionally frozen, not floated to future releases.

Architecture: `GeneralMatrix` provides public storage/arithmetic/solve APIs and
constructs LU, QR, Cholesky, singular-value and eigenvalue decompositions.
`Maths.Hypot` supplies stable norm calculations. Tests and tooling are separate
from the numerical library.

## Changes and contract evidence

- Kept assembly names, version attributes, all public library signatures and
  virtual members, interfaces, exception behavior, alias/copy semantics,
  serialized field names/types, and numerical algorithm operation order.
- Annotated nullable equality inputs and constructor-path storage invariants.
  Random factories use `Random.Shared` with their existing bounds. Removed the
  empty finalizer while retaining `IDisposable`, usable managed storage after
  disposal, and suppression of derived finalizers.
- Migrated all 31 legacy exception attributes to exact-type exception assertions.
  Updated reference assertions to the current API and strengthened deep-copy
  checks. Corrected the inverted test-helper tolerance predicate to the same
  oracle already used by the original console harness. No original assertion
  was deleted, commented out or weakened.
- Added separate matrix, decomposition and eigen tests with analytic expected
  values, reconstruction equations, orthogonality, pivot isolation, multicolumn
  solves, dimension errors and overflow/underflow cases. Eigen tests include
  deterministic dense/scaled matrices; fixtures contain no production data.
- Replaced VSTest with the [native Microsoft.Testing.Platform
  runner](https://learn.microsoft.com/en-us/dotnet/core/testing/migrating-vstest-microsoft-testing-platform).
  Current direct and transitive packages are centrally pinned and locked.
  The unused telemetry extension is excluded completely: its generated hook
  and binary are absent from a clean build. This avoids invoking its old
  ApplicationInsights 2.x contract after upgrading that graph to 3.x.

An isolated comparison compiled **unchanged baseline library source** on .NET
10, then ran the same final suite against it: **100/100 passed**. Reflection
comparison found identical public signatures, virtual members, interfaces and
instance field schemas across **19 surfaces**. This is direct evidence of
source/domain compatibility on the modern runtime; it does not establish
compatibility with the old .NET Framework runtime or binary formatter.

## Measured coverage

All production classes are included, with no method/file/attribute exclusions
and no test assembly included in the denominator. Counts are taken from the
final Release Cobertura report, not estimated from source.

| Module | Line coverage | Branch coverage |
| --- | ---: | ---: |
| GeneralMatrix | 100.00% | 99.09% |
| CholeskyDecomposition | 100.00% | 100.00% |
| LUDecomposition | 100.00% | 100.00% |
| QRDecomposition | 100.00% | 100.00% |
| SingularValueDecomposition | 100.00% | 91.00% |
| EigenvalueDecomposition | 98.98% | 96.42% |
| Maths | 100.00% | 100.00% |
| **Overall** | **1318/1323 = 99.62%** | **828/857 = 96.62%** |

The checker uses integer counts to require **strictly greater than 80%** for
both overall metrics. It rejects absent/duplicate reports, a missing production
class, an empty suite, or any test that did not pass. Every run uses a newly
created output directory. Negative probes verified rejection of exactly 80%
branch coverage, omitted modules, skipped tests and absent reports.

## Reproducibility

Prerequisites: .NET SDK 10.0.401, Bash, Python 3 (standard library only), and
network access to NuGet for the first restore. From the repository root:

```bash
dotnet restore DotNetMatrix.sln --locked-mode
dotnet build DotNetMatrix.sln -c Release --no-restore
bash scripts/verify.sh
```

The script performs restore, Release build, format/analyzer checks, tests and
coverage. It stops on any failure. Latest reports are in the directory printed
by the command under `artifacts/verification.*`. See `ReadMe.txt` for direct test
and coverage commands. Both lockfiles must accompany the source.

Fresh clone verification used a separate empty `NUGET_PACKAGES` directory,
restored every package, and passed the same full gate with identical metrics.
SHA-256 matches across both checkout paths:

| Artifact | SHA-256 |
| --- | --- |
| Library DLL | `6c1ab75f4d9311362c28d8a310342219000a94bd947f0c1ea5784253c24489d9` |
| Library PDB | `2ca5a52a6a71b9410e048cc5049e184026b97cd2d355a9f575e11268c8d31d19` |
| Test DLL | `686e007e3d3a2256c64137e4fe6c9599b2c3fc587ed920db85f3d78995e1f246` |
| Test PDB | `aa92184df25d2622f16886f72a037355bd3b3317f5085b4b8cb546cc98318cc2` |

Compiler artifacts incorporate source-control metadata. A source archive
without Git metadata also passes every test and coverage check, but its PDB/DLL
hashes differ from a checkout. Artifact equality was verified using actual
clones at the same Git revision/remote and identical task source content.

## ATLAS hardness report

Edge cases tested: empty/packed/ragged matrices, invalid indices/dimensions,
array aliasing and mutation, IEEE zero division, singular and rank-deficient
solves, positive-definite detection, repeated/complex eigenvalues, scaled dense
systems, stable norms near overflow/underflow, random bounds, disposal and
explicit legacy behavior characterizations.

Tool/service failures: the original framework build failure is retained as
baseline evidence. Removed MSTest APIs were repaired with equivalent exact
exception checks. Registry checks skip project references. An initial archive
artifact mismatch led to a Git-metadata-preserving clone replay; no build check
was weakened. Package restore is the only external service required by the gate.

Concurrency/load: three workers owned separate new test files; the coordinator
alone changed shared projects, dependencies and numerical source, ran integrated
builds, and writes the Git index. An independent worker reviewed the final
integration. Tests are explicitly serialized to retain legacy harness behavior.
No application load benchmark or performance guarantee is claimed.

Security/data: no credentials or production data were introduced. Telemetry
registration is absent; no network/disk/database mocks are needed for domain
arithmetic. `NuGet.Config` restricts restore sources to the official registry.

Verification: current and fresh-clone `bash scripts/verify.sh` pass; all 100
baseline differential tests pass; API/schema comparison passes; dependency and
vulnerability audits are clean; independent final review passes. Detailed
workspace logs remain in `artifacts/`.

## Remaining limitations and standards

The target migration intentionally requires .NET 10 consumers. It does not
retain .NET Framework 4.0 binary compatibility. Modern .NET removed
[BinaryFormatter](https://learn.microsoft.com/en-us/dotnet/core/compatibility/serialization/9.0/binaryformatter-removal).
No roundtrip compatibility is claimed: the legacy `ISerializable` callback was
already empty and had no deserialization constructor, and remains unchanged.

Under the user's strict preservation rule, pre-existing defects were not
silently corrected:

- Nonunit-diagonal Cholesky solves can be numerically wrong; forward substitution
  divides too late. Separate tests establish factor reconstruction, correct
  unit-diagonal solves, and the independently derived legacy erroneous result.
- `SolveTranspose` returns the transpose of the documented solution.
- Null-left equality throws; equal-value matrices can have different hashes.
- Wide decomposition support is inconsistent; LU construction can throw and
  wide SVD factor dimensions are problematic.

These need separately scoped behavior-change decisions. Rare remaining eigen
and SVD branches are uncovered, as reflected above. Cross-platform execution
beyond Linux arm64 was not observed.

The supplied standards library has no C#/.NET guide. Official Microsoft guidance
was used without claiming compliance with an unavailable standard. The plan uses
Development Plan Creation Guide 2.7 and Autonomous Plan Guide reviewed 2026-10-01.
