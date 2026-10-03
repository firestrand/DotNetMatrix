# DotNetMatrix improvement and feature roadmap

Date: 2026-10-02\
Reviewed revision: `e14019df1f4d63254fd31210ddca206f44b04709`\
Status: proposal for maintainer prioritization; implementation is not authorized by this document.
Revision: 1.2 — independent oracle review refinements.

## Recommendation

Make mathematical correctness and release reliability the next milestone. Follow
with reliable rectangular solvers and a small set of useful public APIs. Measure
performance before changing storage or arithmetic order.

Assumed direction: a small, understandable dense linear algebra library for .NET
applications, with no mandatory third-party runtime dependencies. Actual consumer
workloads and compatibility requirements have not been supplied. Sparse matrices,
GPU execution, and a general-purpose scientific computing framework should remain
outside the initial roadmap unless users establish a concrete need.

## Review baseline

The current library targets .NET 10 and C# 14 and implements dense `double`
matrices, arithmetic, norms, LU, QR, Cholesky, SVD, and real/complex eigensystems.
The modern build, locked dependency graph, and extensive tests are a good base.

For this review, `bash scripts/verify.sh` passed with zero build warnings/errors,
all **100 tests passing**, **99.62% line coverage**, and **96.62% branch coverage**.
The observed run is logged in `artifacts/roadmap-review-verification.log`; its
reports are in `artifacts/verification.9Jl7g4/`.

High coverage does not establish correctness of every public contract. Some tests
intentionally protect incorrect legacy behavior until a behavior-changing release
is approved. The earlier [modernization report](modernization-report.md) records
why those behaviors were retained.

| Finding | Evidence in the current repository | Consequence |
| --- | --- | --- |
| Cholesky forward substitution divides a row after using it to update subsequent rows. | `CholeskyDecomposition.Solve`; `CholeskyNonunitDiagonalSolvePreservesLegacyForwardSubstitutionDefect` in `DecompositionTests.cs` | Some valid positive-definite systems return incorrect solutions. |
| `SolveTranspose` returns the solution transpose. | `Matrix.cs`, `GeneralMatrix.SolveTranspose`; `SolveTransposePreservesLegacyTransposedResult` | Result shape and the documented equation `X*A=B` disagree. |
| Equality operators dereference a null left operand; value equality and hashing use different definitions. | `GeneralMatrix.operator ==`, `operator !=`, `Equals`, `GetHashCode` | Null comparisons can throw; equal matrices may not behave as equal dictionary/set keys. |
| Rectangular solver support is incomplete. | LU's `j < m & LU[j][j] != 0` accesses storage despite the guard; SVD allocates U with `min(m,n)` columns but `GetU` advertises `min(m+1,n)`; `S` indexes all n singular values although storage has `min(m+1,n)` entries | Wide inputs have unsafe indexing or inconsistent factor dimensions. These paths need targeted executable reproducers; they were inspected, not newly executed in this review. |
| `Inverse()` advertises a pseudoinverse but simply calls `Solve(Identity(...))`. | `GeneralMatrix.Inverse`, `Solve`; `QRDecomposition.Solve` rejects rank deficiency | A general rank-deficient or wide Moore–Penrose pseudoinverse is not implemented. Full-column-rank tall matrices are a narrower supported case. |
| Iterative decompositions count iterations without enforcing a limit. | Eigenvalue `tql2`/`hqr2`; SVD constructor's main iteration loop and placeholder comment | Nonconvergence can lack a bounded failure path. Hanging behavior was not reproduced; this is a code-level risk. |
| Decomposition accessors expose mutable internal storage. | Cholesky `GetL`; SVD `GetU`, `GetV`, `SingularValues`; eigenvalue `GetV`, `RealEigenvalues`, `ImagEigenvalues` | Editing returned factors/arrays can alter later solves or diagnostics; repeated-use guarantees need an explicit ownership contract. |
| Custom legacy serialization emits no entries. | `ISerializable.GetObjectData`; `LegacySerializationContractEmitsNoPayload` | There is no working documented persistence roundtrip for the matrix. |
| Automation and distribution documentation are incomplete. | No CI workflow, benchmark project, or LICENSE file found; `ReadMe.txt` is mostly migration guidance; assembly metadata still contains historic placeholders | Platform regressions and consumer adoption need better support. License/provenance and ownership were not legally evaluated. |

This is a code/product review, not a formal coding-standards compliance audit.
The supplied standards library has no C#/.NET standard; no substitute compliance
claim is made.

## Priorities and delivery order

Effort is relative: **S** is a focused API or infrastructure change, **M** touches
several APIs/tests, and **L** involves numerical algorithms or a new representation.
These are planning estimates, not calendar commitments.

| Order | Milestone | Main dependency | Effort | Release outcome |
| --- | --- | --- | --- | --- |
| 1 | Correctness and explicit contracts | Approve changes to documented legacy defects | M–L | Solvers and equality can be trusted within a declared input domain. |
| 2 | CI, documentation, and package readiness | Can begin alongside milestone 1; public release waits for it | M | Every change is checked and consumers can install/use a documented package. |
| 3 | Complete rectangular solving and numerical diagnostics | Define dimensions/convergence policy in milestone 1 | L | Rank-deficient and underdetermined problems have explicit, tested solutions. |
| 4 | Persistence and everyday matrix APIs | Stable shape/equality/ownership contracts | M | Matrices are easier to construct, compare, inspect, and exchange. |
| 5 | Performance measurements and targeted optimization | Correctness gate and representative workloads | M–L | Improvements have measured speed/allocation benefits and bounded numerical impact. |
| 6 | Optional applied features | Core solver/SVD stability and real demand | M–L | Small, demonstrated applications of the core library. |

### Milestone 1 — Correctness and explicit contracts

**R1. Correct Cholesky solving and transposed solving (highest priority, M).**
Move the Cholesky diagonal normalization ahead of dependent row updates. Return
`X` rather than `Xᵀ` from transposed solving. Before implementation, add failing
mathematical regression cases with nonunit diagonal factors, multiple right-hand
sides, and nonsymmetric coefficient matrices. For example, the current Cholesky
fixture `A=[[4,2],[2,10]], B=[[10,0],[32,36]]` should return
`X=[[1,-2],[3,4]]`.

Acceptance: for square consistent systems, verify `A*X≈B` for ordinary/Cholesky
solves and `X*A≈B` for transposed solves. Also test an inconsistent rectangular
transposed solve where `Aᵀ` is tall and full-column-rank: verify expected solution
and shape, and the optimality condition `(X*A−B)*Aᵀ≈0`, not a zero equation
residual. For example, `A=[[1,0,0],[0,1,0]], B=[[2,3,4]]` has least-squares
solution `X=[[2,3]]` and a nonzero residual. Other rectangular/deficient cases
belong to R3/F1. Verify input preservation and existing dimension/error contracts.
Record the approved behavior change and update the corresponding legacy
characterization oracle; retain all unrelated assertions. Use release notes and
a migration example rather than presenting these corrections as behavior-neutral.

**R2. Make equality coherent (high priority, M).**
Define null, subtype, NaN, and signed-zero semantics. Make `==`, `!=`, typed/object
`Equals`, and hashing agree. Consider implementing `IEquatable<GeneralMatrix>`.
Content-based hashing can agree with value equality, but mutable matrices still
must not be modified while used as hash keys. Provide a documented comparer or
an immutable snapshot if consumers need stable keys. Approximate comparison
belongs in a separate method and must not determine dictionary equality.

Acceptance: null comparisons are symmetric; independently constructed equal
matrices have equal hashes and work in dictionary/set lookup; unequal values or
shapes compare appropriately; special values/subtypes follow the approved policy.
The [Microsoft equality guidance](https://learn.microsoft.com/en-us/dotnet/csharp/programming-guide/statements-expressions-operators/how-to-define-value-equality-for-a-type)
and [mutable hash guidance](https://learn.microsoft.com/en-us/dotnet/api/system.object.gethashcode?view=net-10.0)
support these requirements.

**R3. Specify supported shapes and bounded failure (high priority, L).**
Document square/tall/wide and empty-input support per decomposition. Add executable
wide LU/SVD reproducers before choosing either safe rejection or full support.
Define whether decomposition entry points reject nonfinite values, and introduce
explicit convergence budgets/status for iterative algorithms. Avoid changing
all existing constructor exception types or IEEE arithmetic behavior accidentally.
Characterize copy-versus-alias behavior for every public factor/eigenvalue accessor
as well as `GeneralMatrix.Array`. Preserve existing aliases unless a behavior
change is approved; identify accessors whose output mutation changes subsequent
results.

Acceptance: unsupported inputs fail predictably before out-of-range indexing;
factor dimensions match actual storage; documented supported shapes reconstruct;
nonconvergence has a tested, bounded outcome; ordinary existing finite inputs
retain their numerical results. Budgets must be validated against difficult
converging cases, not chosen merely to make a timeout test pass. Ownership tests
must establish which returned arrays/matrices borrow factor storage and which
are independent copies.

### Milestone 2 — CI, documentation, and package readiness

**R4. Add continuous verification and durable regression gates (high priority, M).**
Add a proposed `.github/workflows/verify.yml` for pull requests and branch pushes.
Run the existing gate on Linux first; add Windows/macOS smoke and unit runs, then
architecture checks where runner access allows. Install both Bash and Python
where needed, or introduce a platform-neutral verification entry point. Upload
TRX, coverage, and useful failure logs; retain locked restore and zero-warning
build checks. Add tests for the coverage checker's rejection paths and a durable
public API compatibility baseline. Independent mathematical oracles and selected
mutation tests should supplement coverage rather than chasing 100% mechanically.

Acceptance: a failing test, formatter check, missing module, missing report, or
coverage of exactly 80% makes CI fail; intentional API changes require baseline
review. Evidence and failures are downloadable. Cross-platform results are
reported individually rather than inferred from a Linux pass.

**R5. Create a consumer-facing release package and guide (high priority, M).**
Add a concise README with installation, supported shapes, solve/decomposition
examples, numerical caveats, matrix and decomposition factor mutability/aliasing,
and the current defects until fixed. State that callers must not mutate borrowed
factors while reusing a decomposition; concurrent reuse is not promised without
separate verification. Generate XML API documentation and package metadata; define package ID,
versioning, changelog, and release process. Review provenance/attribution and
select an appropriate LICENSE with the maintainer before public distribution;
the existing origin URL alone is not a completed license review. Replace historic
metadata only with verified project information.

Acceptance: `dotnet pack -c Release` produces a package that a separate sample
project can restore and use; the package carries the agreed metadata, attribution,
and API docs. Test the public NuGet-style consumption path, not just a project
reference. Publishing remains a separately authorized action. Additional target
frameworks should follow actual consumer demand, not automatically reintroduce
unsupported runtimes.

### Milestone 3 — Complete rectangular solving and numerical diagnostics

**F1. Add an explicit SVD pseudoinverse and minimum-norm solver (high value, L).**
Proposed APIs: `PseudoInverse(...)` and `SolveMinimumNorm(...)`, with a documented
relative singular-value cutoff. Support square, tall, wide, and rank-deficient
matrices once SVD factor dimensions are corrected. Keep `Inverse()` behavior
separate until its migration policy is approved; do not silently route every
ordinary solve through a slower/different algorithm.

Acceptance: define the retained-rank matrix `Aτ` by zeroing singular values at the
approved cutoff. Verify the four Moore–Penrose identities against `Aτ` and its
computed pseudoinverse, using scale-aware roundoff tolerances. Against the original
`A`, separately bound truncation error by the discarded singular values; do not
require `A*A⁺*A≈A` at ordinary roundoff tolerance after deliberately discarding
nonzero singular values. For example, `A=diag(1,1e-4)` with relative cutoff `1e-3`
retains `Aτ=diag(1,0)` and has spectral truncation error `1e-4`.

Check analytic rank-deficient examples, zero matrices, scale changes, multiple
right-hand sides, and singular values just below/at/above the cutoff. Define
finite, nonnegative tolerance inputs and tie behavior; reported numerical rank
must use the same policy. Distinguish exact consistent solves from least-squares
and minimum-norm results at the retained numerical rank: deliberate truncation
can increase the original-system residual. The exact-rank cases should satisfy
the identities against the original matrix. This is a new capability, not a claim
that the existing `Inverse()` already meets these criteria.

Test `SolveMinimumNorm` itself independently of pseudoinverse identities. For the
untruncated underdetermined system `A=[[1,1]], B=[[2]]`, require the 2×1 result
`X=[[1],[1]]`, `A*X≈B`, and orthogonality to the nullspace direction `[1,-1]ᵀ`.
These checks distinguish the minimum-norm answer from other exact solutions.
Also test inconsistent systems against the retained-rank least-squares policy.

**F2. Add solution diagnostics and reusable solve workflows (medium value, M).**
Offer residual norms, numerical rank, conditioning information, and convergence
status through a small proposed result/options API. Document how callers can
reuse the existing LU/QR/Cholesky decomposition objects for multiple solves before
adding another factorization abstraction. Explain that current factor accessors
can expose mutable storage. If safe snapshots are needed, propose additive copy
or read-only accessors; do not silently change the existing alias contracts. A
read-only view is not an immutable snapshot if another alias can still mutate it.
Consider an additive `TrySolve` for expected singular/rank-deficient failures;
invalid arguments should remain explicit errors.

Acceptance: diagnostics are checked against independently calculated residuals;
repeated solves with unmodified factor storage produce consistent results and
leave inputs/factors unchanged. Test and document how edits through existing
borrowed accessors affect results; proposed copy accessors must be independent.
Unavailable estimates are clearly distinguished from exact or computed values.
Do not compute an expensive SVD automatically just to populate optional diagnostics.

### Milestone 4 — Persistence and everyday matrix APIs

**F3. Add explicit, versioned JSON persistence (medium value, M).**
Use a separate DTO/converter built on `System.Text.Json`; do not retrofit an
undocumented BinaryFormatter payload. Proposed schema: format version, row and
column counts, storage order, and flat numeric values. Agree on special-value
handling and input-size limits before freezing that new external contract.
[Microsoft's custom converter documentation](https://learn.microsoft.com/en-us/dotnet/standard/serialization/system-text-json/converters-how-to)
describes the supported mechanism.

Acceptance: roundtrip rectangular/empty matrices with dimensions and values
preserved; reject malformed shapes, length mismatches, overflow, and unknown
versions; nonfinite values follow an explicit policy. Existing matrix ownership
and legacy serialization behavior remain separate from the new format.

**F4. Add small matrix conveniences (medium value, S–M per API).**
Proposed additions: an indexer, explicitly named row/column extraction, diagonal
construction, row/column concatenation, caller-provided packed-buffer copying,
and `ApproximatelyEquals` with explicit absolute/relative tolerances. Add seeded
or caller-supplied random generation without promising identical streams across
runtime versions. Preserve the existing `Array` aliasing contract and identify
whether each new method copies or borrows storage.

Acceptance: examples are concise; operations check dimensions and preserve inputs
unless explicitly mutating; caller buffers are sized correctly; approximate
comparison has documented NaN/infinity/zero behavior; fixed random seeds reproduce
results within the documented environment. Select a few demonstrated use cases
first rather than adding the entire list at once.

### Milestone 5 — Performance measurements and targeted optimization

**R6. Establish benchmarks before optimizing (medium priority, M).**
Create a proposed benchmark project for multiplication, transpose, norms, copy,
factorization, and repeated solves across small/medium/large square and rectangular
inputs. Record time and allocated bytes; separate cold setup from repeated use.
[BenchmarkDotNet](https://benchmarkdotnet.org/articles/overview.html) provides
reproducible Release-mode measurement tools; add it as a development-only package
if this milestone is selected. Performance of the current code was not measured
in this review, so no speedup target or bottleneck is asserted.

Then evaluate cache blocking, reusable destination buffers, and SIMD individually.
Floating-point reassociation can change results: require numerical regression
checks and measured benefit. A contiguous-storage matrix could be a future new
type or backend, but replacing jagged storage in `GeneralMatrix` would disrupt
public `Array` aliases and field layout. Avoid that change as an incidental
optimization. Native BLAS is optional only after a deployment/maintenance decision.

Acceptance: checked-in benchmark instructions reproduce a baseline on specified
hardware/runtime; each optimization shows a useful workload-specific benefit,
its allocation impact, and acceptable numerical error. Routine CI should not fail
on noisy shared-runner timing unless stable hardware and budgets are established.

### Milestone 6 — Optional applied features

Choose one after confirming a real use case:

| Feature | User value | Dependency | Acceptance |
| --- | --- | --- | --- |
| Least-squares regression example/helper | Fits coefficients with residual/rank reporting | F1/F2 for deficient systems | Known coefficient fits, residual checks, and documented deficient-data behavior. |
| PCA and low-rank approximation example | Demonstrates SVD for dimension reduction | Correct SVD shapes, explicit centering policy | Centering/reconstruction verified; retained variance and low-rank error match analytic examples. |
| Immutable matrix snapshot | Safe shared data and value-based keys | R2 and ownership decisions | No exposed writable storage; equality/hash stability; copying cost documented. |

Start with samples before promoting these into permanent core APIs. Generic
numeric matrices, sparse storage, GPU kernels, and automatic parallel execution
are deferred: they require demand, benchmarks, and an architecture decision that
is not supported by the present repository evidence.

## Decisions before implementation

1. Approve a release strategy for changes to Cholesky, `SolveTranspose`, equality,
   input validation, and convergence failures. Recommend a clearly identified
   correctness release with a migration guide; decide major/additive compatibility
   treatment based on actual consumers. This roadmap does not override the prior
   modernization requirement to preserve behavior.
2. Confirm the dense-library direction and the most valuable consumer workloads.
3. Confirm package identity, provenance/license, and supported target frameworks.
4. For numerical features, agree on supported shapes, relative tolerance defaults,
   special-value policies, and failure reporting before documenting the APIs.
5. Expand each selected milestone into a focused implementation plan with tests
   before behavior changes, a GOTCHA spec where required, and ATLAS evidence at
   completion. Retain `bash scripts/verify.sh` and >80% overall line/branch coverage
   as minimum gates; add capability-specific numerical checks.

Recommended first deliverable: R1 regression tests and an approved correction
policy, R4 Linux CI, and R5 consumer documentation. Complete R2/R3 before claiming
a broadly reliable numerical release; then prioritize F1 over optional applied
features. No calendar dates are proposed without maintainer priorities and capacity.

## Independent review outcome

An independent oracle reviewed the original draft and revisions 1.1 and 1.2.
Final verdict for revision 1.2: **No more changes needed.**

A second independent oracle checked revision 1.2 against the current solver,
equality, decomposition, serialization, and coverage-checker code and confirmed
the same verdict. No further substantive revisions were recommended.

The revisions distinguish exact and least-squares solve acceptance, define
truncated-SVD pseudoinverse and numerical-rank checks, make decomposition ownership
explicit, and independently verify the minimum-norm solver. The approval covers
this roadmap proposal, not implementation authorization or new build results.
