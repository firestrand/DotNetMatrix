# DotNetMatrix

Dense `double` matrix arithmetic and Cholesky, LU, QR, singular-value and
eigenvalue decompositions for **.NET 10 / C# 14**. The library has no third-party
runtime dependencies. This repository contains a local `2.0.0-preview.1`
correctness release under the [MIT license](LICENSE). The package ID is for local
verification; public package publication remains a separate release action.
Read [provenance and release policy](docs/provenance-and-release.md).

## Build, test and try the package

Install the exact SDK in [global.json](global.json), Bash and Python 3. From a
fresh checkout:

```bash
dotnet restore DotNetMatrix.sln --locked-mode
dotnet build DotNetMatrix.sln -c Release --no-restore
bash scripts/verify.sh
```

The full gate checks formatting, all unit tests, strict **greater than 80%**
production line and branch coverage, verifier rejection cases, a reviewed public
API baseline and independent local NuGet consumption. Each run has fresh report
directories under `artifacts/`; no skipped tests or missing production classes
can pass. To run unit tests alone:

```bash
dotnet test --project DotNetMatrix_Test/DotNetMatrix_Test.csproj -c Release
```

To create and consume the local preview (after the Release build):

```bash
bash scripts/package-smoke.sh
```

Its generated consumer restores `DotNetMatrix.LocalPreview` from the new local
feed only, in a fresh cache, without a project reference. Inspect its generated
project under `artifacts/package-smoke.*/consumer` for a working install example.
There is no assumed nuget.org installation command for this unpublished ID.

## Solving systems

```csharp
using DotNetMatrix;

var a = new GeneralMatrix(new double[][] {
    new[] { 4.0, 2.0 }, new[] { 2.0, 10.0 }
});
var b = new GeneralMatrix(new double[][] { new[] { 10.0 }, new[] { 32.0 } });
var x = a.Chol().Solve(b); // x = [1, 3]^T; a*x = b

// Reuse the factorization for repeated right-hand sides.
var factor = a.Chol();
var secondSolution = factor.Solve(b);

// Underdetermined system: choose the solution with minimum Euclidean norm.
var row = new GeneralMatrix(new double[][] { new[] { 1.0, 1.0 } });
var minimumNorm = row.SolveMinimumNorm(
    new GeneralMatrix(new double[][] { new[] { 2.0 } })); // [1, 1]^T
```

`Solve` keeps its existing algorithm selection: LU for square systems, QR for
tall systems. It does not silently opt into a costly SVD or solve rank-deficient
systems. Use `PseudoInverse` / `SolveMinimumNorm` explicitly for square, tall,
wide or rank-deficient systems. Their relative cutoff retains singular values
**strictly greater** than `rcond * largestSingularValue`; the default is
`max(rows, columns) * 2^-52`. A deliberately discarded nonzero singular value
changes the represented system; test reconstruction against the truncated
matrix, not an impossible exact identity for the original input.

`SolveTranspose(b)` solves `X*A = b` and returns `X` in its documented
orientation. For a rectangular least-squares system, its residual need not be
zero: check the normal-equation optimality condition as well as dimensions.
`Inverse()` retains the legacy algorithm selection; it is not the new general
SVD pseudoinverse. Forming an inverse usually costs more and amplifies numerical
error compared with solving directly.

`SolveWithDiagnostics` and `SolveMinimumNormWithDiagnostics` return the solution,
residual Frobenius norm, numerical rank and convergence status. Condition number
computation is opt-in through `includeCondition: true`: ordinary solving needs a
separate SVD, while minimum-norm solving reuses its existing factors.
`TrySolve` returns false for expected singular/rank-deficient failures while
invalid arguments still throw. Convergence exhaustion is an exception, so a
returned diagnostic result is converged rather than a partial iterate.

## Shape and numerical contracts

| Operation | Supported domain |
| --- | --- |
| Matrix arithmetic / storage | Rectangular dimensions, including empty matrices where documented |
| Cholesky | Finite, nonempty square input; solve requires symmetric positive definite |
| LU | Finite, nonempty square or tall factorization; solve requires square nonsingular input |
| QR | Finite, nonempty square or tall input; solving requires full column rank |
| SVD | Finite, nonempty square, tall or wide input; economy factors for wide input |
| Eigenvalue decomposition | Finite, nonempty square input; symmetric and nonsymmetric paths |

Decompositions reject NaN/infinity and unsupported shapes at the boundary.
Iterative decompositions have a bounded convergence budget and report exhaustion
explicitly. Rank and conditioning depend on floating-point scale and tolerance;
near singularity is a numerical decision, not a universal boolean. Double
precision can overflow/underflow even for finite input. Validate residuals and
choose scale-aware tolerances for your problem. There is no promise that a
given factorization may be mutated or solved concurrently.

## Storage, equality and persistence

Matrices are mutable and store jagged arrays. The jagged-array constructor and
`Array` expose aliases. `ArrayCopy`, packed copies and `Copy()` produce independent
data. Legacy factor getters retain their alias behavior; use explicitly named
copy accessors when isolation matters. A read-only view does not guarantee an
immutable snapshot while another alias can mutate its storage.

| Factor access | Ownership |
| --- | --- |
| Cholesky `GetL()` | Aliases internal L; `GetLCopy()` is independent |
| LU `L`, `U`, `Pivot`, `DoublePivot` | Independent copies |
| QR `Q`, `R`, `H` | Independent copies |
| SVD `GetU()`, `GetV()`, `SingularValues` | Alias stored factors; `GetUCopy()`, `GetVCopy()`, `GetSingularValuesCopy()` isolate them |
| SVD `S` | Independent diagonal matrix |
| Eigen `GetV()`, `RealEigenvalues`, `ImagEigenvalues` | Alias stored factors; corresponding `Get*Copy()` methods isolate them |
| Eigen `D` | Independent real block-diagonal matrix |

Equality compares exact runtime type, dimensions and values with `Double.Equals`:
NaNs compare equal and signed zeros compare equal. Operators accept null safely;
hashes reflect values. **Do not mutate a matrix used as a dictionary/set key.**
Approximate comparison is a separate operation and validates finite nonnegative
absolute/relative tolerances. The convenience APIs add indexed access, independent
row/column extraction, diagonal/concatenation factories and caller-owned packed
copy buffers. Caller-provided `Random` supports repeatable sequences within the
same runtime; it is not a cross-runtime stream guarantee.

```csharp
GeneralMatrix firstRow = a.GetRow(0); // independent copy
GeneralMatrix diagonal = GeneralMatrix.Diagonal(1.0, 2.0);
double[] packed = new double[a.RowDimension * a.ColumnDimension];
a.CopyTo(packed); // row-major; columnMajor: true selects column-major
bool close = a.ApproximatelyEquals(a.Copy(), absoluteTolerance: 1e-12,
    relativeTolerance: 1e-12);
GeneralMatrix seeded = GeneralMatrix.Random(2, 3, new System.Random(42));
```

JSON persistence uses a versioned finite-only row-major schema with explicit
dimensions, storage order and flat values. Deserialization enforces an element
limit and rejects inconsistent dimensions and unknown versions. The historical
empty `ISerializable` callback remains preserved; .NET 10 cannot replay
BinaryFormatter payloads. The new JSON path does not claim binary compatibility.

```csharp
using System.Text.Json;

var jsonOptions = new JsonSerializerOptions();
jsonOptions.Converters.Add(new MatrixJsonConverter(maxElements: 1_000_000));
string json = JsonSerializer.Serialize(a, jsonOptions);
GeneralMatrix? restored = JsonSerializer.Deserialize<GeneralMatrix>(json, jsonOptions);
```

The field names are `formatVersion`, `rows`, `columns`, `storageOrder` and
`values`; version 1 uses `"row-major"`. Persisted matrices own independent storage.

See [CHANGELOG](CHANGELOG.md) for intentional breaking corrections and
[release policy](docs/provenance-and-release.md) for API baseline maintenance.

The [least-squares sample](samples/LeastSquares) demonstrates rank-aware fitting
with independent numerical checks and runs as part of the full gate. Performance
workloads and measured optimization decisions are documented in
[the benchmark project](benchmarks).
