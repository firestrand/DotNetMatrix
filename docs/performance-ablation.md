# Performance optimization and ablation study

## GOTCHA spec

**Goals:** find and implement measured performance improvements without changing
numerical operations, public APIs, serialization, storage ownership or exceptions.

**Objectives:** retain only candidates with repeatable throughput or allocation
benefits; isolate each mechanism with an otherwise identical ablation; verify
unchanged behavior and non-regressing production coverage.

**Tasks:** inspect existing evidence and hot paths, benchmark isolated candidates,
reject inconclusive/regressing candidates, implement demonstrated improvements,
then validate the actual production implementation and document results.

**Capabilities:** local C# changes, pinned .NET/BenchmarkDotNet, tests, native CPU
affinity and statistical reports. No runtime/dependency changes, unsafe memory,
arithmetic reassociation, unapproved API or storage redesign, or remote publication.

**Health:** no user budget specified; use two process launches, warmup and repeated
measurements, with bounded build time and persisted logs. Re-measure only when new
changes, high uncertainty or contradictory evidence justify it.

**Attributes and constraints:** preserve summation order and floating-point bits;
zero shapes must not begin inspecting previously ignored storage. Malformed
borrowed arrays must retain their exception types. Every promotion needs its own
ablation, not a combined before/after explanation.

**Users:** maintainers and numerical library consumers. Public contracts remain
unchanged; timing claims apply only to the observed runtime/hardware/workloads.

**Runtime:** per-run benchmark processes and owned fixture arrays; persisted raw
reports and source variants. Process isolation and a fixed CPU reduce migration
noise but do not eliminate shared-host load or thermal effects.

**Beliefs and intentions:** previous multiplication candidates were slower or
inconclusive and must not be promoted on that evidence. Investigate data movement
and repeated jagged-array indexing first, then other measured hot paths.

## Study protocol

Benchmark the existing production entry point as a calibration alongside a frozen
scalar kernel and isolated candidates. Row caching versus the scalar kernel removes
only repeated jagged indexing; bulk copying versus row caching changes only the
copy primitive. Norm experiments retain each row's original ascending sum order.
Use square, tall, wide, single-row, single-column, tiny and empty shapes. Report
managed allocations, uncertainty and all negative/inconclusive results.

The study promoted a guarded bulk-copy implementation in `GeneralMatrix.Copy`.
The evidence and limitations are recorded below; other experiments remain
unpromoted.

## Discovery evidence

The [discovery reports](../benchmarks/baselines/data-movement/discovery/) retain
all measurements, allocations and the statistical comparison for 63 cases. On
Ubuntu 24.04.5 Arm64, .NET 10.0.12 / SDK 10.0.401 and BenchmarkDotNet 0.15.8,
row caching reduced the scalar copy kernel's mean time by 15–27% for medium,
large, tall and wide matrices. Adding bulk copying reduced those means by a
further 50–72%. These isolate indexing and copy primitive changes respectively.
The unmodified production copy calibrated closely to the scalar control.

Unconditional bulk copying regressed the 4096×1 case by 21%. Consequently it
will not be promoted unchanged. A guarded candidate retains the original loop
for zero/one-column matrices; additional 2/4/8-column workloads check that the
guard does not simply move the regression to another narrow shape.

The infinity-norm row-cache candidate benefited several larger shapes, but
regressed the single-column kernel by 9% and had inconclusive small results.
It remains an experiment; production norm arithmetic is unchanged. Earlier
multiplication experiments also remain unpromoted, as documented in
[the existing performance report](performance.md).

Each case used two launches, four warmup iterations per launch and twelve
100 ms measured iterations per launch, with outlier removal disabled and a 5%
Mann–Whitney comparison. Parent and child processes inherited Linux `taskset`
CPU 5 affinity (verified from `/proc`); BenchmarkDotNet's own affinity-setting
attempt failed in the discovery run, and high-priority scheduling was denied.
The guarded run uses inherited affinity without requesting that failing setting.
Shared-host load and thermal variation remain limitations; this is evidence for
the observed machine, not a portable guarantee.

## Guard selection and production confirmation

The [second study](../benchmarks/baselines/data-movement/guarded/) ran 60 cases
before changing production. Allowing bulk copying from two columns regressed
4096×2: 48.05 µs for the guarded candidate versus 39.92 µs for production.
Four/eight-column results were mixed. The final implementation therefore keeps
the original scalar loop below 16 columns. Sixteen is a conservative threshold
supported by these workloads, not a universal crossover for every processor.

The [production confirmation](../benchmarks/baselines/data-movement/production/)
ran 48 cases with the same job settings. It compares the frozen scalar control,
two-column guard, 16-column guard and actual production copy. Allocation counts
are identical for every corresponding copy variant. The two guard variants use
the same implementation with different thresholds, isolating the guard decision.

| Shape | Original production mean (µs) | Final production mean (µs) | Before/after speedup | Final control/final production |
| --- | ---: | ---: | ---: | ---: |
| 16×16 | 1.612 | 0.380 | 4.24× | 4.27× |
| 64×64 | 10.568 | 2.378 | 4.44× | 4.13× |
| 192×192 | 45.240 | 16.295 | 2.78× | 2.79× |
| 384×192 | 87.644 | 34.065 | 2.57× | 2.57× |
| 192×384 | 87.892 | 31.890 | 2.76× | 2.71× |
| 1×4096 | 9.967 | 0.973 | 10.24× | 9.21× |
| 4096×1 | 37.138 | 36.052 | 1.03× | 0.96× |
| 4096×2 | 39.919 | 40.768 | 0.98× | 1.00× |
| 4096×4 | 50.926 | 49.462 | 1.03× | 1.01× |
| 4096×8 | 70.534 | 68.111 | 1.04× | 1.00× |
| 1×1 | 0.0429 | 0.0442 | 0.97× | 0.93× |
| 32×0 | 0.491 | 0.424 | 1.16× | 0.99× |

Before and after production measurements are separate jobs; the final column
provides a same-job control to help expose drift. The mean reductions for the
wide-row cases agree with the isolated experiments, and BenchmarkDotNet classifies
all six as faster than the final scalar control at its 5% comparison threshold.
Raw reports include uncertainty and every sample; ratios here divide arithmetic
means rather than reproducing BenchmarkDotNet's distribution-based ratio column.

There is a small tradeoff: the 1×1 production copy increased from 42.90 to
44.16 ns (about 3%) across jobs and was classified slower than the final scalar
control. Narrow and empty production cases were classified the same as that
control. The empty before/after difference is not claimed as an optimization:
its unchanged scalar path and same-job control indicate run variation. The
implementation is retained for its substantial measured wide-row gains, with
the tiny-case cost explicitly disclosed. No universal speedup is claimed.

Only `Copy` changed. `Clone` and owned decomposition getters delegate to it,
but their end-to-end performance was not measured here. `ArrayCopy`, numerical
kernels, dependency versions and public APIs remain unchanged. The larger norm
candidate would need another guarded study before promotion; current mixed
evidence does not justify changing it.

## Reproduction and provenance

[Benchmark instructions](../benchmarks/README.md#data-movement-ablations) provide
the exact build and run commands. Use the final run's method filters to reproduce
the 48-case confirmation instead of timing every copy variant:

```bash
taskset -c 5 dotnet benchmarks/bin/Release/net10.0/DotNetMatrix.Benchmarks.dll --filter '*CopyAblation.ScalarControl*' '*CopyAblation.GuardedBulkCopy*' '*CopyAblation.WideRowsOnly*' '*CopyAblation.Production*' --warmupCount 4 --iterationCount 12 --launchCount 2 --iterationTime 100 --outliers DontRemove --statisticalTest 5% --exporters JSON CSV --artifacts artifacts/performance-ablation/production-reproduction
```

Each report directory retains CSV summaries, full JSON measurements and the
exact harness source as `DataMovementAblation.cs.txt`. To reproduce an earlier
phase, copy its snapshot over the benchmark source in a disposable checkout.
Original production is available at commit `4f83ec3`; use that commit for the
discovery and guarded phases. The [study manifest](../benchmarks/baselines/data-movement/study-manifest.json)
records source hashes, job settings and verification results. A fresh checkout
uses `dotnet restore DotNetMatrix.sln --locked-mode` before the documented build.

## ATLAS hardness report

**Edge cases tested:** exact signed-zero, infinity and custom NaN bits; independent
rows even when source rows alias; logical prefixes of longer rows; both scalar
and bulk-path null/short/missing rows; ignored storage for empty dimensions;
single-column copies. Bulk copying keeps the original allocation order and
falls back to the original element loop for truncated borrowed rows, retaining
`IndexOutOfRangeException` rather than `Array.Copy`'s argument exception.

**Failure handling:** benchmark setup rejects incorrect values/ownership before
timing. The failed BenchmarkDotNet affinity setting was replaced by verified
inherited Linux affinity. Scheduling priority remained unavailable and is
disclosed. No assertions, thresholds or approval checks were weakened.

**Concurrency/load:** benchmarks ran sequentially on CPU 5, with two independent
launches per case and no concurrent builds/tests. Shared-host and thermal effects
are not eliminated. Matrices retain borrowed mutable storage; this change does
not introduce an atomic snapshot or support concurrent source mutation.

**Security:** no unsafe memory, new dependency, network operation or credential
handling was added to production. Reports contain runtime/workload evidence.

**Verification:** `bash scripts/verify.sh` completed strict Release build with
zero warnings/errors, both formatting checks, vulnerability checks, all 61
Python verification tests, unchanged API comparison, all 179 C# tests, coverage
comparison, independent package consumer, numerical-candidate validation and the
least-squares sample. Production coverage is 1598/1611 lines (99.19%) and
1054/1096 branches (96.17%), passing the reviewed non-regression baselines.
The entire command exits 1 at the pre-existing final standards gate because
`EXC-2026-001` still lacks independent protected approval/issue evidence. This
study neither activates that pending exception nor claims the full gate is green.

**Disposition:** retain the measured wide-row copy improvement; retain the scalar
path below 16 columns; reject unconditional copying and the insufficient guard;
leave norm and multiplication experiments unpromoted. The reproducible harness
can detect recurrence of copy regressions. Further hardware-specific tuning is
optional and requires new evidence.
