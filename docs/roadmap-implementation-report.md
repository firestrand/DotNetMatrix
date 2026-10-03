# Roadmap implementation evidence

Date: 2026-10-02. Source: improvement-roadmap.md revision 1.2.
The implementation is a local correctness preview. Public distribution is not
claimed or authorized; public package identity remains a maintainer gate.
The maintainer subsequently selected MIT; LICENSE and package license metadata
now record that choice. Historical upstream attribution remains preserved.

## Requirement audit

| Item | Implemented result | Authoritative evidence |
| --- | --- | --- |
| R1 | Correct Cholesky normalization and transposed-solve orientation | DecompositionTests.CholeskyNonunitDiagonalSolveSatisfiesEquation; ModernMatrixTests.SolveTransposeReturnsTheDocumentedSolution; MatrixContractTests rectangular optimality oracle |
| R2 | Null-safe operators, exact runtime types, Double.Equals semantics, content hashing and IEquatable | MatrixContractTests null/subtype/NaN/signed-zero/dictionary tests; ModernMatrixTests hash/value assertions |
| R3 | Finite nonempty boundaries, LU/QR square/tall factors, square Cholesky/eigen, economy SVD for every nonempty rectangle | DecompositionContractTests validation, reconstruction, orthogonality, budget and ownership cases; DecompositionTests and EigenvalueTests |
| R4 | Linux x64/arm64 gates, Windows/macOS unit jobs, failure artifacts, API baseline, dynamic coverage inventory | .github/workflows/verify.yml; tools/ApiBaseline; docs/public-api.txt; nine Python gate tests; selected mutations below |
| R5 | Consumer guide, XML docs, changelog, provenance/release policy and local NuGet preview | README.md; CHANGELOG.md; docs/provenance-and-release.md; scripts/package-smoke.sh independent fresh-cache consumer |
| F1 | Pseudoinverse and minimum-norm solve with finite nonnegative relative cutoff and strict tie policy | MatrixCapabilityTests four Moore–Penrose identities for exact/truncated rank, analytic minimum norm/nullspace, zero/scaled/multiple-RHS systems and cutoff below/at/above tests |
| F2 | Residual/rank/optional original-system conditioning, successful convergence status, TrySolve and additive copies | MatrixCapabilityTests diagnostics and repeated same-instance LU/QR/Cholesky solves; fresh factor getters/input snapshots; borrowed-L edits affect subsequent solve, copy edits do not |
| F3 | Opt-in finite-only version-1 row-major JSON with strict fields and allocation bounds | MatrixJsonTests: 39 cases covering empty/rectangular roundtrips, duplicate/missing/unknown metadata, count/overflow/nonfinite and allocated-row limits |
| F4 | Indexer, copied row/column extraction, diagonal/concatenation, Span packed copy, approximate comparison and caller Random | MatrixCapabilityTests dimension/input preservation, buffer orders, seeded streams and special-value/tolerance/overflow checks |
| R6 | Development-only benchmark suite and measured blocking/reuse/SIMD experiments | 113 measured cases in five benchmarks/baselines CSV files; docs/performance.md; candidate correctness checks in verify.sh |
| Optional applied milestone | Least-squares sample with analytic full-rank and deficient designs, before adding a permanent helper | samples/LeastSquares runs coefficient/rank/residual oracles in verify.sh |

Generic numeric storage, sparse/GPU/native backends, automatic parallelism, new
target frameworks and additional applied features remain explicitly deferred
as specified by the roadmap. No production multiplication optimization was
promoted: this host's measurements did not establish a throughput benefit;
SIMD also changed a cancellation answer. Buffer reuse showed an allocation
benefit at a measured throughput cost. See the raw intervals and caveats in the
performance report.

## Verification and reproducibility

`bash scripts/verify.sh` passed locally with zero compiler warnings/errors,
formatting verification, **158 passing tests**, and all production source parts
included in the coverage inventory. Coverage: **1,590/1,603 lines (99.19%)** and
**1,046/1,088 branches (96.14%)**. Every production source part exceeds 80% for
both metrics. Nine Python verifier tests also pass. Local packaging confirms
DLL/XML/README/provenance contents and independently restores the package to
exercise Cholesky, minimum-norm solving and JSON.

Evidence: `artifacts/roadmap-final-verification.log` and
`artifacts/verification.BMCHE9/coverage-summary.json`. The last test refinement
re-fetches copied factor getters after each solve; its exact final state was
verified by the fresh-snapshot run below.

A local Git clone overlaid with the complete current source snapshot was tested
from an **initially empty NuGet cache**. Its full verification script passed with
the same 158 tests and coverage counts. Evidence:
`artifacts/roadmap-fresh-verification.log` (original reports in that temporary
snapshot's `artifacts/verification.035Mdr/`). No bin/obj or cached packages were
copied. At verification time the source was uncommitted; this was an exact
worktree snapshot rather than a test of the remote branch.

The local clone initially had a filesystem origin rather than the GitHub origin,
which changed Source Link metadata. After matching the origin and commit and
rebuilding, DLL and portable PDB hashes matched the workspace exactly:

```text
DLL bee8c7956f513c7d6d8b1e4b7c7c455e30ee7193abb85f9c84686172655de924
PDB 7de657f13b231f39b72856d7d62f3599b21673d06f466b232ac2a6cc222b72a5
```

Evidence: `artifacts/roadmap-reproducible-rebuild.log`. Reproduction requires
the pinned SDK, locked restore and matching source/control metadata; dependencies
are installed with `dotnet restore DotNetMatrix.sln --locked-mode`. Build with
`dotnet build DotNetMatrix.sln -c Release --no-restore`, then run
`bash scripts/verify.sh`. README and benchmarks/README give standalone commands.

## ATLAS hardness report

Edge cases tested: exact, least-squares and underdetermined equations; deficient
and deliberately truncated rank; multiple/zero-column right-hand sides; very
small scalar coefficients and near-maximum right-hand sides; null/subtype/NaN
equality; approximate infinity/NaN/zero/overflow; empty and malformed JSON;
unsupported decompositions and both eigen iteration paths. Minimum-norm solving
applies factors directly with scaled projections, avoiding a reciprocally
overflowing intermediate pseudoinverse and an overflowing U-transpose product.
Finite-input decompositions can still suffer floating-point overflow/underflow
or ill-conditioning; no universal arbitrary-scale accuracy guarantee is made.

Tool failures: the gate rejected Coverlet's initially empty report. Diagnostics
showed automatic dynamic exclusions included DotNetMatrix; explicit testconfig
now instruments the production assembly. Distinct source files of partial types
are accepted while duplicated class/file entries are rejected. A failed initial
benchmark child build exposed PrivateAssets propagation; BenchmarkDotNet remains
only in the nonpackable benchmark project and its child builds now run.

Selected mutation evidence: temporarily returning the transpose again killed
the SolveTranspose tests; changing cutoff `<=` to `<` killed the cutoff-tie test.
Both test runs exited 2 with assertion failures; original files were restored
before verification. Logs: `artifacts/transpose-result-mutation.log`,
`artifacts/cutoff-tie-mutation.log`, and `artifacts/selected-mutations.json`.
These two targeted mutations strengthen those contracts; they do not imply a
full mutation-coverage score.

Concurrency/load: legacy writable aliases remain; copy accessors are independent.
Sequential factor reuse is verified, concurrent mutation/reuse is not promised.
Benchmark timing is observational on a shared heterogeneous ARM host, not a CI
threshold. Library storage remains jagged; no runtime dependency was introduced.

Security: JSON validates dimensions/count/finite values before allocation and
bounds rows even with zero columns. Transport byte limits are caller-owned.
No credentials, production data, or external publication were used.

Independent implementation oracle reviewed the numerical APIs, tests and gates,
requested two numerical scaling corrections and fresh factor re-fetch checks,
and confirmed the final result: **No more changes needed.**

Passed: all locally executable roadmap gates and numerical acceptance checks.
Unobserved: hosted Linux x64/arm64, Windows and macOS workflow results; workflows
are configured but were not pushed/run remotely. Public release follow-ups:
upstream provenance review, registered package identity and assembly metadata
policy require maintainer decisions before public distribution. The local preview
uses the maintainer-selected MIT license and preserves historic ownership claims. Publication is
still a separate authorization.
