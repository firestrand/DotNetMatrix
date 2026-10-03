# Roadmap implementation plan

Version: 1.0.0; 2026-10-02. Source: improvement-roadmap.md revision 1.2.
Guide: Development Plan Creation Guide 2.7; autonomous guide reviewed 2026-10-01.
Mode: implementation; bounded parallel; rolling-wave task expansion.
Baseline: e14019d plus user-owned ignore and roadmap edits, preserved.
Owner: repository maintainer. Integration/review: root agent.

## Execution contract and GOTCHA spec

Goals: implement the roadmap's dense numerical capabilities and release gates.
Objectives: correct R1/R2; define and test R3 domains; implement F1–F4;
verify R4/R5 locally; establish R6 measurements and an applied sample.
Tasks: triggered by the implementation goal; terminate only with evidence for
each requirement or an explicitly unresolved external decision.
Capabilities: local tools, numerical tests, registry/documentation lookups;
publication and invented licensing are off limits.
Health: no fixed time/token budget; inspect every failed check before retry;
report phase progress and retain >80% line and branch coverage.
Attributes: preserve unrelated assertions, jagged storage aliases and legacy
serialization; release approved corrections as a documented major change.
Constraints: .NET 10, no mandatory third-party library runtime dependencies.
Users: dense double matrix consumers; no asserted GPU/sparse demand.
Runtime: persistent code/tests/docs; bounded workers with disjoint file ownership;
integrate shared API changes centrally, then verify the complete solution.
Beliefs: existing 100-test coverage protects both valid behavior and named defects.
Intentions: replace only approved defect characterizations with mathematical oracles.

Policies: equality uses exact runtime type and Double.Equals (NaNs equal, signed
zeros equal); content hashes require unmodified keys. Decompositions accept finite
nonempty supported shapes, use bounded iteration, and retain existing aliases.
SVD exposes economy factors for wide matrices. Cutoff retains singular values
strictly above rcond*smax; default rcond=max(rows,columns)*2^-52.
JSON is finite-only, row-major, version 1, with a configurable element limit.
Ordinary Solve/Inverse retain algorithm selection; new minimum-norm APIs use SVD.
Package/license identity and public publication remain maintainer decisions.

## Compliance and boundaries

| Requirement | Disposition |
| --- | --- |
| Language standard | No C# standard exists in supplied library; compliance unverified. |
| Behavioral tests before implementation | Required per packet; named legacy expectations updated for approved R1/R2/R3 policies. |
| Data fidelity | Analytic roadmap fixtures and existing tests; no production data. |
| Provider boundary | NuGet/build tooling only; JSON uses native System.Text.Json. |
| Publication/license | Local package smoke proceeds; no license invented or publication performed. |
| Full integration | Locked restore, build, format, all tests, coverage, package smoke, API and checker gates. |

## Fact ledger and evidence index

| ID | Required fact | Evidence to establish |
| --- | --- | --- |
| R1 | Cholesky and transpose solve satisfy equations/LS optimality | Numerical regression tests; verify.sh |
| R2 | Equality, hashing, null/special values/subtypes agree | Matrix contract tests; verify.sh |
| R3 | Supported shapes, ownership, finite input and convergence bounded | Decomposition contract tests; verify.sh |
| R4 | CI gates reject missing/invalid evidence and APIs are reviewable | Workflow, checker tests, API baseline comparison |
| R5 | Documented local package consumed independently | Pack and package-consumer smoke command |
| F1 | Cutoff pseudoinverse/minimum norm meet analytic and MP conditions | Numerical capability tests; verify.sh |
| F2 | Residual/rank/optional conditioning, reusable solves, independent copies | Diagnostics and ownership tests |
| F3 | Versioned JSON accepts valid shapes and rejects malformed input | Persistence boundary tests |
| F4 | Convenience APIs validate shapes/tolerances and preserve ownership | Matrix convenience tests |
| R6 | Benchmarks measure workloads before optimization | Benchmark instructions/results and numerical comparison |
| LOCAL-SAMPLE | One applied sample demonstrates rank-aware solving | Independently verified least-squares sample |

Real data manifest: none; pure numerical fixtures from roadmap and existing suite.
Provider manifest: no external boundary in numerical behavior; package smoke is
isolated from the library's project reference and uses a local package feed.

## Work packets and phase DAG

| Packet | Owner | Exclusive write scope | Dependency |
| --- | --- | --- | --- |
| P1 Matrix capabilities | root | Matrix.cs, new matrix partial/API files, matrix capability tests | baseline |
| P2 Decompositions | decomposition worker | decomposition classes, numerical guards/options, decomposition tests | policies above |
| P3 Release gates | tooling worker | scripts, CI, package configuration, README/CHANGELOG, API checker | baseline; final API integration |
| P4 Persistence | JSON worker | MatrixJson.cs and JSON tests | policies above |

Workers may run concurrently only in those scopes. A contract change pauses
dependent packets until integration agrees the new contract. Root owns this plan
and final evidence. Failed workers return findings; root resumes without resetting
unrelated work. Same DAG is executable sequentially. Max three workers; no
automatic repeated retries.

### Phase 0 — Baseline (support)
Protect all inherited modernization facts. Run `bash scripts/verify.sh`; observe
100 passing tests, locked clean build, coverage. Rollback: no behavior edits.

### Phase 1 — Mathematical contracts (capability)
Introduce R1/R2/R3/F1/F2. Tests precede implementation. Numerical equations,
cutoff ties, malformed/nonfinite input, ownership and iteration-budget cases are
observable acceptance. Run full verification after integrated worker changes.
Rollback: revert only this phase's scoped changes; never discard user edits.

### Phase 2 — Consumer APIs and persistence (capability)
Introduce F3/F4. Test rectangular/empty JSON, malformed schema and overflow;
test extraction/concatenation/buffers, special-value approximate equality and
seed reproducibility. Implement then verify the complete suite.

### Phase 3 — Release reliability (support)
Enable R4/R5 with platform CI, checker rejection tests, reviewed API snapshot,
consumer README/XML docs/changelog and local package smoke. Cross-platform CI
results cannot be claimed from this Linux machine; license/publication remain
separate human decisions. Local reproduction must be verified.

### Phase 4 — Measurements and applied sample (capability)
Establish R6 workloads and LOCAL-SAMPLE after numerical integration. Expand
benchmark tasks from observed results. Apply optimizations only with measured
benefit and numerical evidence; optional representations remain deferred.

### Phase 5 — Completion audit (support)
Inspect every ledger requirement against actual code/tests/results. Write ATLAS
report with edge cases, errors, concurrency assumptions, security and command
results. Never equate a coverage percentage with numerical correctness.

## Dispatch horizon

For P1/P2/P4: Test task depends on phase 0 and the policies above; uses only
named analytic inputs and current tests; produces failing regression evidence.
Implement task depends on that test and makes those assertions green; verify
with the complete suite after integration. P3 introduces no domain behavior,
protects existing gates, and verifies script negative paths plus package smoke.
P2/P4 independent GREEN and the integrated 156-test gate closed the initial
numerical dispatch horizon. R6 then dispatched numerical prototype validation
before benchmark experiments; it retained all prototypes outside production.
The final report binds the completed phases to observed commands and evidence.

## Changelog

1.0.0: initial full-scope plan and bounded work packet contracts.
