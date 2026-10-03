# Dense matrix performance experiments

## Scope and reproduction

The development-only `benchmarks/DotNetMatrix.Benchmarks.csproj` uses
[BenchmarkDotNet 0.15.8](https://benchmarkdotnet.org/changelog/v0.15.8.html),
the latest stable version verified against NuGet on 2026-10-02. Its lockfile pins
the full dependency graph. The library has no BenchmarkDotNet runtime dependency.
Follow [the benchmark commands](../benchmarks/README.md) after a locked restore.

The workload suite measures time and managed allocations for multiplication,
transpose, copy, one/infinity/Frobenius norms, LU/QR/Cholesky/SVD/eigen
factorization, and repeated solves. Sizes are 16, 64, and 192; square dimensions
are size×size, tall matrices 2size×size, wide matrices size×2size. Rectangular
multiplication uses a compatible right matrix with size columns. LU/QR
factorizations include square/tall inputs; Cholesky/eigen require square inputs;
SVD includes all three shapes. Every solve uses four right-hand sides.

`GlobalSetup` creates matrices and reusable factors outside the measured operation.
Factorization benchmarks include factor creation and its allocation, while repeated
LU/QR/Cholesky solve benchmarks reuse a factor. The repeated minimum-norm workload
reuses a pseudoinverse created in setup; it measures its application, including the
result allocation. Each fresh-factorization benchmark therefore measures the setup
cost separately from repeated application. JIT startup/process cold start are not
represented by these steady-state means; use BenchmarkDotNet's ColdStart jobs
separately when application startup matters.

## Observed baseline

Measured on 2026-10-02 on Linux Ubuntu 24.04.5 LTS, Cortex-X925/Cortex-A725,
10 physical cores, Arm64 RyuJIT armv8.0-a; SDK 10.0.401 and runtime 10.0.12.
The short job used one launch, one warmup, three measured iterations targeting
100 ms each. This shared heterogeneous host was not isolated or frequency-pinned.
The full default job and repeated runs on stable hardware are required for a
production optimization decision. CSV snapshots retain the means, confidence
intervals, standard deviations and managed allocated bytes per operation.

[Multiplication candidate snapshot](../benchmarks/baselines/2026-10-02-multiplication.csv):

| Candidate | 32×32 mean | Bytes/op | 128×128 mean | Bytes/op |
| --- | ---: | ---: | ---: | ---: |
| Existing multiplication | 13.25 µs | 9,552 | 1,076.54 µs | 136,272 |
| Cache blocking | 27.65 µs | 9,272 | 1,754.75 µs | 135,224 |
| Reused destination and scratch | 20.56 µs | 0 | 1,241.17 µs | 0 |
| SIMD with packing included | 17.85 µs | 18,544 | 1,083.85 µs | 270,448 |

At size 128, the baseline's reported 99.9% confidence interval half-width was
1,607.003 µs, exceeding its mean; SIMD's was 612.259 µs. These sparse measurements
cannot establish a small throughput difference. The 32×32 blocking experiment
was about twice as slow in this run; it is not evidence against every possible
blocking arrangement or workload. The allocation measurements provide a clearer
tradeoff: destination/scratch reuse eliminated per-call allocations, while packing
for SIMD approximately doubled them.

The core suite also completed all **105 cases**, covering every configured size,
shape, operation, factorization, and repeated solve. Together with eight candidates,
113 benchmarks executed with measured results. The core run took 120.2 seconds;
the candidate run took 11.55 seconds including generated-project compilation.
The recorded snapshots are:

- [Matrix operations: 54 cases](../benchmarks/baselines/2026-10-02-operations.csv).
- [Square/tall LU and QR: 18 cases](../benchmarks/baselines/2026-10-02-lu-qr.csv).
- [Square Cholesky/eigen and repeated LU/Cholesky: 15 cases](../benchmarks/baselines/2026-10-02-square-factors.csv).
- [Square/tall/wide SVD and repeated minimum norm: 18 cases](../benchmarks/baselines/2026-10-02-svd.csv).

For example, size-192 square SVD factorization measured 49.42 ms and 886.69 KiB
allocated, while applying the setup pseudoinverse to four right-hand sides measured
71.71 µs and 13.58 KiB. Size-192 wide SVD measured 78.04 ms and 2,064.29 KiB;
its repeated application measured 142.47 µs and 25.58 KiB. These are different
operations whose setup distinction is deliberate; they are not interchangeable
solver performance claims. The snapshots retain all error intervals rather than
hiding uncertainty behind rounded headline numbers.

## Candidate correctness and decisions

The three candidate implementations remain inside the benchmark project. No
production arithmetic, jagged-array layout, alias, or public API changed for these
experiments. Cache blocking uses tiles of 32 and preserves each element's ascending
summation order. Reuse preserves the original column-copy/dot-product order and
writes each destination cell from a new sum on every call. SIMD includes the cost
of packing/transposing the right matrix and changes summation grouping.

`--validate` checks square/tall/wide products at sizes 1, 3, 16, 32, 64, 128,
and 192 against existing multiplication. It verifies unchanged inputs and repeated
use after deliberately poisoning an output cell with NaN. The finite normwise error
bound is `1e-12 * ||A||F * ||B||F`. This is a forward-error comparison to the existing
algorithm, not a promise of a small relative error when the result nearly cancels.

A separate cancellation fixture multiplies `[1e100,1,-1e100,2]` by a column of
ones. Existing multiplication returns 2 on this runtime; SIMD returns 3. The error
of 1 passes a normwise bound proportional to the huge operands, but it is material
relative to the result. This explicitly demonstrates why reassociation needs a
consumer-specific numerical policy and cannot silently replace the original method.

The initial decisions are:

- Keep the blocking prototype out of production: measured throughput decreased.
- Keep reusable buffers experimental: zero allocations is useful for repeated
  workloads, but this run did not show a throughput gain, and buffer ownership and
  alias rules need a separate API design before promotion.
- Keep SIMD experimental: packing increased allocations, a throughput benefit is
  unproven, and cancellation changed the observed result.

A future contiguous-storage type or native BLAS backend requires its own deployment
and compatibility decision. These experiments provide no basis to change
`GeneralMatrix.Array` aliases or promise shared-factor thread safety. Routine CI
runs correctness/build gates; it does not impose timing thresholds on shared runners.
