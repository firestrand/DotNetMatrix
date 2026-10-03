# DotNetMatrix modernization plan v1.0.1

Guide version: Development Plan Creation Guide 2.7 (2026-09-27), Autonomous Plan Guide reviewed 2026-10-01.
Mode: migration; plan type: brownfield; planning horizon: full; strategy: bounded parallel.
Baseline: 683be7c9ad4ac33ef63a361485558681036cb14d, master; execution branch modernization/dotnet10.
Authority: local implementation, package restore, generated mathematical unit fixtures, verification. Closeout: stage; no commits, publication, or deployment.
Ownership baseline: pre-existing untracked .serena/ remains user-owned. No staged/tracked changes existed.
Scope: seven numerical source files, two projects, existing tests, build/test/coverage documentation and automation.
Standards: no C#/.NET guide exists in the supplied library; official Microsoft guidance fallback used under autonomous best judgment; no compliance claim. No C# standards-compliance claim is possible. Graphify/Context7 unavailable; Serena available.

## Acceptance and compliance

| Requirement | Observable gate | Status at baseline |
| --- | --- | --- |
| R1: latest stable runtime/language/build | .NET 10 SDK pinned, C# 14, SDK projects | Pending |
| R2: current stable dependencies | NuGet registry queried; versions and transitive packages locked | Pending |
| R3: preserve domain/API/schema | original tests retained; additional reconstruction and contract tests | Pending |
| R4: reproducible clean build | locked restore, Release build with warnings as errors, clean clone replay | Pending |
| R5: >80% line AND branch coverage | production-only Cobertura, strict threshold script | Pending |
| R6: inspectable phases | execution ledger, evidence logs, independent review | Pending |
| Library C# standard | required source absent | Unavailable; official Microsoft guidance used; unavailable library standard remains disclosed |

## Fact ledger and evidence index

Each Tier-1 fact is owned by the operator in the product/platform capacity. Evidence will be bound to the final integrated snapshot.

| ID | Given / When / Then | Kind / lifecycle | Trace | Evidence / oracle |
| --- | --- | --- | --- | --- |
| DNM.BUILD.v1 | Given the documented SDK and a fresh checkout, when locked restore and Release build run, then compilation succeeds reproducibly. | Requirement / Active | R1,R2,R4 | scripts/verify.sh; SDK/config/lockfile |
| DNM.CONTRACT.v1 | Given supported matrix inputs, when public operations run, then dimensions, storage aliasing, arithmetic, error types and existing serialized contract remain compatible. | Requirement / Active | R3 | original MatrixTests plus ModernMatrixTests; baseline source and explicit expected values |
| DNM.NUMERICS.v1 | Given supported decomposition inputs, when decompositions run, then factors reconstruct inputs and valid legacy solves satisfy their equations within floating-point tolerance; the pre-existing nonunit Cholesky solve defect remains characterized separately. | Requirement / Active | R3 | DecompositionTests, EigenvalueTests; independent mathematical equations |
| DNM.COVERAGE.v1 | Given all unit tests, when coverage is collected, then overall production line and branch coverage each exceed 80%. | Requirement / Active | R5 | scripts/check-coverage.py; Coverlet Cobertura |

Unexpected legacy defects are characterization findings, not silent changes to Tier-1 contracts. Preserve behavior unless a correctness repair is authorized by the requested outcome and verified by a regression test; call out any such repair explicitly.

## Data and provider boundaries

No production/customer data or external service boundaries exist. Fixtures are existing test values and deterministic mathematical examples generated within the user-authorized coverage expansion. New test packets use these fixtures only. Package provider: NuGet, development-only MSTest/Test SDK/Coverlet; no library runtime package dependencies.
Registry verified 2026-10-02: MSTest.TestFramework and MSTest.TestAdapter 4.4.1, coverlet.MTP 10.1.0, Microsoft.Testing.Platform/TRX extensions 2.4.1. VSTest and its Test SDK/collector were replaced by native MTP. All transitive dependencies pinned to current stable in Directory.Packages.props; unused telemetry extension excluded to avoid the upstream old ApplicationInsights binary contract. SDK installed/latest stable official release 10.0.401, runtime 10.0.12. Pin all versions and lock transitives. Network allowed for restore; no secrets required.

## Execution packets

Coordinator/integrator: root (GPT-6); maximum three workers plus coordinator. Shared interface revision: baseline public API, net10.0, MSTest 4.4.1. Only root writes project/config/lockfiles, original tests, production source, ledgers, Git index and integration build outputs. Workers do not run contested build resources concurrently; coordinator runs aggregate checks. No worktree is needed for workers with disjoint new-file ownership. Execution branch isolates baseline branch; pre-existing .serena/ is preserved in place.

| Task | Dependencies | Owner | Owned files | Acceptance / handoff |
| --- | --- | --- | --- | --- |
| P0.1 baseline | none | root | artifacts/baseline, docs | record build/test outcomes and missing coverage |
| P1.1 toolchain | P0.1 | root | projects, global.json, common properties, configs | dotnet restore/build, lockfiles |
| P2.1 API migration | P1.1 | root | production, original tests | original assertions retained; test suite passes |
| P3.1 decomposition tests | baseline interface; P1.1 for execution | decomposition_tests | DecompositionTests.cs | meaningful factors/solve/error tests; handoff findings |
| P3.2 eigen tests | baseline interface; P1.1 for execution | eigen_tests | EigenvalueTests.cs | symmetric/nonsymmetric eigen equations; handoff findings |
| P3.3 matrix tests | baseline interface; P1.1 for execution | matrix_tests | ModernMatrixTests.cs | constructors, aliasing, arithmetic, errors, legacy contract |
| P3.4 integrate coverage | P2.1,P3.1,P3.2,P3.3 | root | coverage scripts/settings | scripts/verify.sh; strict >80% both |
| P4.1 independent review | P3.4 | reassigned finished worker | read-only integrated changes | numerical/API/build/test review |
| P4.2 fresh replay and closeout | P4.1 | root | docs, ledger, task-owned index | fresh build/test gate, hashes, staged review |

One packet retry after evidence-backed correction; larger uncertainty requires replan. Failed workers must stop before reassignment. Shared-interface changes pause dependent packets; source-only compatible fixes require updated tests and integrated rerun. Independent work continues when another packet is blocked. Same DAG is executable serially.

## Phases

All phases are support/migration phases; no new user-facing numerical capability is introduced. Each protects all four facts above. All commands run from repository root. No deployment, benchmarks, CLI/service observability or new numerical dependencies are in scope.

### P0: baseline and specification
Trigger: user modernization request. Terminate with manifest/API map, recorded legacy outcomes and this GOTCHA spec. Gate: dotnet build DotNetMatrix.sln and dotnet test DotNetMatrix.sln, capturing the legacy failure rather than pretending it passes. Outcome: old .NET Framework 4.0 cannot build on Linux (MSB3644); legacy test command discovers no tests and has no coverage tooling. Rollback: remove only task documentation if abandoned.

### P1: toolchain foundation
Depends on P0. Convert projects to SDK style; target .NET 10 / C#14, exact SDK pin, lockfiles, deterministic build, compiler nullability and built-in analyzers. Retain assembly identity and serialization attributes. Replace VS2010 test assembly with latest MSTest and coverage collector. Gate: dotnet restore and dotnet build DotNetMatrix.sln -c Release. Outcome: successful compile. Rollback: restore task-owned project diff only.

### P2: compatibility migration
Depends on P1. Repair removed MSTest exception attributes with explicit exact-type exception assertions; do not delete or weaken assertions. Initialize nullable decomposition storage according to actual constructor paths, annotate equality null inputs, use modern array cloning/copy patterns when safe. Preserve numerical operation order, virtual APIs, serialized field names, exception types, random ranges. Add regression evidence before any correctness change. Gate/demo: dotnet test DotNetMatrix.sln -c Release. Outcome: original suite and legacy console harness assertions pass. Rollback: task-owned source/test diff.

### P3: test coverage
Depends on P1,P2 and packets. Tests first, with numerical reconstruction, analytic fixtures, exact exceptions, aliasing and error coverage. Observe coverage, target remaining branches with meaningful cases. Gate/demo: bash scripts/verify.sh. Outcome: no failing tests and >80% line and branch coverage across all production code, without excluding difficult algorithms. Rollback: task-owned new tests/config only.

### P4: verification and handoff
Depends on P3. Independent source/API/test evidence review, ATLAS report, exact fresh-checkout commands, clean replay with empty outputs/cache, dependency status, formatted files, deterministic artifact comparison and recorded final hashes. Gate/demo: bash scripts/verify.sh in current and fresh workspaces. Outcome: verified staged changes; no commit/publish. Preserve user files/index. If gate fails, return to relevant packet and rerun affected evidence.

## GOTCHA spec

Goals: modern supported, reproducible numerical library without regressions.
Objectives: R1–R6 gates above. Tasks: trigger user request, terminate only after final verified gate and review. Capabilities: shell, official docs/NuGet, Serena, bounded agents; off-limits user work and external publication. Health: no user-set time/token budget; monitor tests/build/coverage, one unchanged-attempt retry maximum. Attributes: preserve API/schema/business behavior; never weaken assertions or hide coverage gaps. Constraints: Linux arm64 available, no legacy VS2010 framework runtime. Users: library consumers and maintainer. Runtime: persisted source/config/tests/evidence; decision loop follows dependency DAG and failures trigger evidence-based replanning. Beliefs: no runtime dependencies and seven source modules; intentions: direct supported-runtime migration because algorithms use portable BCL APIs.

## Changelog

1.0.0 — 2026-10-02: baseline and execution contract.

## Execution decisions and final gate

- Strict no-behavior-change rule takes precedence over correcting unrelated existing defects. Cholesky was not changed; its independent rational-output characterization documents the defect. No original assertions were removed or weakened.
- All 62 original unit tests were migrated; removed exception attributes became 31 exact-type assertions. Reference-identity assertions now use AreSame/AreNotSame; deep-copy evidence was strengthened. The inverted numerical test-helper predicate was corrected to match the original console harness tolerance oracle.
- Modern source edits: nullable equality parameters and constructor-path invariant annotations; Random.Shared; eliminate empty finalizer while retaining IDisposable and SuppressFinalize; consistent compiler formatting. Numerical operation order and fields retained.
- Native MTP package/runner configuration replaces legacy VSTest settings. Central transitive pins satisfy latest-stable requirement; dormant telemetry hook and assembly are excluded before registration because its old API cannot safely load ApplicationInsights 3.x.
- Verify and verify-final: bash scripts/verify.sh. Coverage checks reject absent/partial/empty/skipped reports and strictly enforce >80% overall line and branch coverage. Every one of seven production classes must be represented.
- Current result: 100/100 passing, 1318/1323 lines (99.62%), 828/857 branches (96.62%). All individual modules exceed both thresholds. Independent review passed. Fresh clone/cold-cache replay and deterministic DLL/PDB comparison are P4 completion gates.
- Fact surfaces: projects, global.json, Directory.Build.props, Directory.Packages.props, NuGet.Config, lockfiles, all test files, scripts/verify.sh and scripts/check-coverage.py. Generated fixture authorization covers P3.1–P3.3 only; no production data needed.
- New-source baseline differential: identical final test suite passes 100/100 against unchanged baseline source compiled on .NET10. Public/virtual/interface signatures and field schemas match 19 reflected surfaces. This is source compatibility evidence, not proof of old runtime or BinaryFormatter replay.

1.0.1 — 2026-10-02: lock final native MTP graph, guardrail disposition, evidence and verification commands.

Final P4 closure: current and cold-cache clone gates pass; all four DLL/PDB hashes match; final public API/schema comparison and baseline differential tests pass. Independent review passed. Execution ledger: AGENT_STATE.md.
