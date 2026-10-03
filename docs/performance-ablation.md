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

Production changes are pending measurements. This document is not a performance
claim until raw measured results and correctness checks are recorded below.
