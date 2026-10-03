# Standards adoption audit

Audit date: 2026-10-02. Adopted source: ENG-STD-DOTNET-XPLAT
1.1.0-draft.1, supplied by the repository owner. Its corporate approval status
remains proposed; this audit does not supply that approval.

## GOTCHA: repository enforcement changes

**Goal:** make applicable code and repository controls conform to the supplied
standard while preserving numerical behavior, public signatures, JSON schema and
every existing test assertion.

**Outcome:** strict Release build, deterministic API/serialization checks,
non-regressing line and branch coverage, explicit native platform qualification,
controlled package sources and fail-closed exception validation.

**Tradeoffs:** keep the legacy solution and test identities because this is an
existing repository. Use the pinned .NET 10 SDK analyzer baseline plus explicit
style, security and platform rules; illustrative recommended profiles are not
additional mandatory rules. Preserve legacy exception types and serialization
assertions. The necessary obsolete metadata constructors need independent approval.

**Constraints:** no fabricated approver, issue, organization approval, external
branch policy or native-platform result. No publication or remote policy changes
are authorized by this audit. Embedded sources are disclosed, not assumed to be
retrievable from an uncommitted source revision.

**Hazards:** local warning relaxation must not affect CI, nullable or CA1416 checks;
package cache hits must not masquerade as fresh trusted restores; a successful test
process must not hide an empty or reduced suite; checkout-controlled approval
metadata must not authorize its own exception.

**Acceptance:** verify each applicable rule against observable evidence, document
non-applicable host requirements, and retain an explicit failed gate for any
unapproved exception or unavailable external control.

**Tasks and users:** the maintainer's adoption request triggers this work; local
verification and an explicit disposition for every rule terminate the audit.
Independent approvers own the unresolved exception and protected release decisions.

**Capabilities and health:** use local edits, canonical CLI checks and read-only
remote evidence. Do not publish, alter branch rules or create approval records.
No task cost budget was specified; verification uses bounded operations and a
30-minute CI job limit, with reports recording failures rather than blind retries.

**Runtime and intentions:** helpers use per-run reports and temporary directories;
the catalog, registry, source and locks persist in the repository. Unknown evidence
fails closed. Existing behavior is characterized by the unchanged API snapshot
and tests; preserve it while repairing mechanical and governance gaps.

## Findings requiring owner action

1. **FAIL — DOTNET-REPO-004 / DOTNET-GOV-004:** the serialization test still needs
   obsolete framework metadata constructors to preserve its existing assertions.
   [EXC-2026-001](serialization-exception.md) is a concrete, narrowly scoped
   **pending** proposal, not authorization. An independently verifiable approver,
   open HTTPS remediation issue and protected evidence are missing. The full gate
   fails on this finding. Replacing the constructors with reflection to conceal
   their warning, or deleting the assertions, would not resolve it.
2. **PARTIAL — DOTNET-GOV-001/002/003/004:** local catalog, exception validation,
   CODEOWNERS and daily CI checks are present. The current remote branch has no
   branch protection (read-only GitHub API returned 404). Required independent
   owner reviews, an authenticated protected evidence job, and owner notification
   for exceptions have not been established. A digest supplied by PR code is not
   authentication. No external settings, issues or approvals were created.
3. **PARTIAL — DOTNET-BASE-006 / DOTNET-CI-003:** native Ubuntu 24.04 Arm64 tests
   passed locally. The revised four-platform workflow has not run on hosted CI.
   Its profiles are qualification targets, not certified support claims.
4. **UNVERIFIED — new-revision source retrieval:** Source Link retrieval passed
   for committed `911ce18`, and fails as expected for modified sources mapped to
   that old SHA. Current changes must be committed and available remotely before
   their own mapped checksums can pass. Local symbol identity and reproducibility
   checks are separate evidence.
5. **UNVERIFIED — DOTNET-BASE-004 / DOTNET-GOV-006:** organization-wide patch SLA,
   responsible corporate approval roles and global effective date are outside
   this repository. The supplied standard explicitly remains a draft. The project
   adoption does not fabricate those decisions.

## Changes and preserved contracts

Production C# uses file-scoped namespaces, explicit accessibility, PascalCase
private methods and complete public XML documentation. Unused imports were removed.
The SDK analyzer baseline is pinned to `10.0`; style, security and platform
diagnostics are explicit. IDE0005 runs through a dedicated info-severity formatter
gate because its build analyzer requires generated XML documentation on every
project; test and tool projects do not need artificial public XML comments.
Unused private code is enforced by IDE0051/IDE0052. SDK sample
`10-recommended` is illustrative: selecting it would require changing preserved
legacy exception types or approving further exceptions, neither warranted here.

Release builds explicitly select CI mode and fail compiler/analyzer and MSBuild
warnings. Local relaxation remains optional and cannot weaken nullable or CA1416
errors. A CI target rejects relaxation and disabled analyzers. No public API,
numeric algorithm, exception type, JSON schema or existing assertion was changed.
An independent executable-token comparison confirmed the production changes are
limited to namespaces, imports, private method names and the internal Maths type.

NuGet uses a single explicit HTTPS feed, exact mappings for all locked package IDs,
an explicit vulnerability source and a repository-local cache namespace derived
from the feed/mapping configuration hash. Trust-policy changes use a new empty
namespace; arbitrary global or cross-job caches are not imported. All five lock
graphs remain unchanged. Weekly dependency/action updates are configured, and
updates still need ordinary review plus deliberate mapping/lock regeneration.

Packages include MIT licensing, XML documentation, repository commit metadata and
matching portable symbols. The Source Link checker validates approved URLs,
bounded downloads and PDB source checksums, and explicitly reports embedded source
not verified by that retrieval check. The only embedded production document in the
observed build is the generated framework assembly-attribute file; it contains no
credentials or business data. It remains disclosed rather than claimed remotely
retrievable.

## Platform and restore profiles

| Qualification target | Fixed runner | Required test process | Current evidence |
| --- | --- | --- | --- |
| Ubuntu 24.04 LTS x64 | `ubuntu-24.04` | X64 | Revised hosted job not run |
| Ubuntu 24.04 LTS Arm64 | `ubuntu-24.04-arm` | Arm64 | Local Ubuntu 24.04.5 / .NET 10.0.12 / Arm64 passed; revised hosted job not run |
| Windows Server 2025 x64 | `windows-2025` | X64 | Revised hosted job not run |
| macOS 15 Arm64 | `macos-15` | Arm64 | Revised hosted job not run |

No Windows client, Windows Arm64, macOS x64, Alpine, enterprise Linux or additional
TFM is newly certified by this audit. Runner labels were checked against
[GitHub's official runner reference](https://docs.github.com/en/actions/reference/runners/github-hosted-runners).
Fixed OS labels still receive image servicing; CI records the observed image,
OS, runtime, architecture, source SHA and configuration hashes.

All five projects use SDK 10.0.401, net10.0, C# 14, Release,
`ContinuousIntegrationBuild=true`, deterministic portable PDBs and strict warnings.
The library ships portable IL with no RID, trimming or AOT profile. Test, sample,
benchmark and API-tool executables are framework-dependent development tools;
there is no deployed host, container or native package dependency. Every matrix
job has isolated `obj`, `bin` and artifacts. Restore uses the same Release/CI
properties as consuming builds; lock drift fails without unlocked recovery.

## Rule dispositions

PASS is local evidence within the stated scope. PARTIAL/FAIL/UNVERIFIED must not be
interpreted as overall compliance. Every canonical ID appears below; the catalog
records the exact source ranges and enforcement contracts.

| Rule | Verdict | Evidence or boundary |
| --- | --- | --- |
| DOTNET-BASE-001 | PASS | Root stable SDK pin; maintenance ownership retained |
| DOTNET-BASE-002 | PASS | No preview SDK or silent major roll-forward |
| DOTNET-BASE-003 | PASS | SDK and net10.0 TFM configured separately |
| DOTNET-BASE-004 | PARTIAL | Project adoption explicit; corporate SLA/approval unverified |
| DOTNET-BASE-005 | PASS | net10.0 / C# 14 throughout |
| DOTNET-BASE-006 | PARTIAL | Exact targets above; new hosted qualification outstanding |
| DOTNET-REPO-001 | N/A | Existing legacy repository may retain .sln |
| DOTNET-REPO-002 | PASS | Production, tests, tools, scripts and docs separated |
| DOTNET-REPO-003 | PASS | Stable public namespace; acyclic references; purpose-specific test assembly |
| DOTNET-REPO-004 | FAIL | Automated style/warning controls pass; obsolete test exception pending; new-revision remote source evidence outstanding |
| DOTNET-DEP-001 | PASS | Central package versions; reviewed update automation configured |
| DOTNET-DEP-002 | PASS | Five locked graphs; matching Release profile; cold/warm restores verified locally |
| DOTNET-DEP-003 | PASS | Exact reviewed IDs, cleared feeds, controlled policy-specific cache namespace |
| DOTNET-DEP-004 | PASS | Fresh audit of all five projects; no reported vulnerabilities; unknown/high/critical findings rejected |
| DOTNET-API-001 | PASS | Package validation enabled; unchanged complete public API snapshot and consumer test |
| DOTNET-API-002 | N/A | No API compatibility suppression |
| DOTNET-API-003 | PASS | Existing local SemVer preview identity; no publication or replacement of deployed artifacts |
| DOTNET-CI-001 | PASS | Thin provider wrapper around scripts/verify.sh |
| DOTNET-CI-002 | PARTIAL | Strict build/artifact provenance prepared; current changes uncommitted and new hosted artifacts absent |
| DOTNET-CI-003 | PARTIAL | Fixed native profiles and architecture tripwire; hosted execution outstanding |
| DOTNET-TEST-001 | PASS | MTP, existing assertions preserved, nonempty successful TRX, production inventory, non-regressing line/branch baseline |
| DOTNET-PKG-001 | PASS | Portable library package and symbols; host deployment mode does not apply |
| DOTNET-CTR-001 | N/A | No container product or service |
| DOTNET-PLAT-001 | PASS | Native path APIs for tooling; canonical LF machine output; no runtime file I/O |
| DOTNET-PROC-001 | PARTIAL | No C# process boundary; provenance queries use verified absolute installations, bounded concurrent output and cleanup; Windows ACL/hosted cleanup qualification outstanding |
| DOTNET-FS-001 | N/A | No secret-bearing/private runtime file creation |
| DOTNET-FS-002 | PASS | Platform temp APIs, exclusive temporary directories and scoped cleanup |
| DOTNET-PLAT-002 | PASS | CA1416 enabled and remains an error under local relaxation; no platform-specific production calls |
| DOTNET-PLAT-003 | N/A | No P/Invoke/native library dependency |
| DOTNET-SEC-001 | PARTIAL | Runtime, audit, analyzer, input guards and source trust verified; protected release governance outstanding |
| DOTNET-CFG-001 | N/A | Pure numerical library; validated explicit options, no service IConfiguration host |
| DOTNET-SEC-002 | PASS | No credentials or production data added; build provenance records only explicit non-secret metadata |
| DOTNET-OBS-001 | N/A | No deployed host/logging boundary; explicit command/report output only |
| DOTNET-OBS-002 | N/A | No service/distributed trace boundary |
| DOTNET-PERF-001 | PASS | Existing measured benchmark record retained; no hot-path arithmetic change; candidate correctness checks pass |
| DOTNET-TIME-001 | PASS | No civil-time contract; governance uses explicit UTC, exclusive expiration, leap/date boundary tests |
| DOTNET-GLOB-001 | PASS | Ordinal/invariant/LF machine contracts; JSON tests across four cultures and API byte checks across locale environments |
| DOTNET-GOV-001 | PARTIAL | All 43 IDs/source digest synchronized; required remote controls explicitly unimplemented |
| DOTNET-GOV-002 | PARTIAL | Canonical CLI and repository owner identified; corporate/independent approver roles unresolved |
| DOTNET-GOV-003 | PARTIAL | Local automated checks and independent source reviews complete; protected approvals/hosted evidence incomplete |
| DOTNET-GOV-004 | FAIL | Pending exception rejected; no independent authenticated approval, open issue or protected branch/evidence job |
| DOTNET-GOV-005 | PASS | Broad legacy adoption explicitly requested; existing contracts/assertions preserved |
| DOTNET-GOV-006 | PARTIAL | Exact draft version pinned; no invented global RFC/effective date |

## Verification commands and results

`bash scripts/verify.sh` restores the policy-specific cache in locked mode, builds,
checks both formatter profiles, audits all projects, runs Python gate tests,
compares the API, runs MTP with line/branch coverage, consumes the package, and
checks numerical candidates and samples. In clean hosted CI it also retrieves
committed production sources. Its final standards gate is deliberately not skipped.

Final full local run: **168 unit tests and 61 Python gate tests passed**. The
Release build produced **zero warnings and zero errors**, and both formatting
profiles passed. Coverage was exactly 1590/1603 lines
(99.19%) and 1046/1088 branches (96.14%), matching committed 911ce18. No existing
coverage baseline was weakened; the minimum suite size is now 168.
Package DLL/docs/license/commit/symbol checks and
the independent fresh-cache consumer passed. NuGet reported no vulnerabilities in
any of the five projects. CI relaxation rejection and nullable/platform errors
under local relaxation were exercised through real MSBuild fixtures.

`bash scripts/verify.sh` exited **1 solely at the final standards gate**, which
listed the pending exception, missing independent approval, missing HTTPS issue
and missing protected evidence. This is an intentional unresolved failure, not a
green gate. Results: `artifacts/verification.Ci8Hq1/coverage-summary.json` and
`artifacts/verification.Ci8Hq1/*.trx`; full run log was captured locally at
`/tmp/dotnet-standards-final-verify.log`.

A separate clean build-output directory at a different absolute checkout path,
with the same uncommitted source overlay and HEAD metadata, reproduced the library
DLL and portable PDB byte for byte. The locked restore used the trusted warm cache
after the separately verified cold restore. Evidence:
`artifacts/standards-reproducibility.json`. DLL SHA-256:
`ac034382dc1f582e921e0102b6f1cf4bf4210529679970b9179782743315bb16`;
PDB SHA-256:
`74f933ea9edc7cf13be2c7d3741ae72669e66a48efa3df349e4f629c1651ebdf`.
These are local reproducibility results, not newly committed release provenance.

An initial parallel cold restore terminated MSBuild worker processes; the
evidence-based serial locked retry succeeded. A subsequent cold restore into the
new feed-policy cache namespace and its warm consuming build succeeded. No lock
file was regenerated to recover.

Source retrieval report: `artifacts/standards-source-link-committed.json` verifies
12 remote source checksums at committed 911ce18 and discloses one generated
embedded document. Current dirty-source retrieval demonstrably exits 1 on a
   checksum mismatch, as required. DLL/PDB mismatch, architecture mismatch, missing
CI architecture, stale/closed/unknown approval evidence and suppression/catalog
bypasses are tested as failures. Build inputs are recorded in
`artifacts/standards-build-provenance.json`; these reports are ignored local
evidence, not signed release attestations.

## ATLAS hardness report

**Edges:** exact coverage thresholds, regression ratios, missing modules/reports,
empty/skipped suites, serious/unknown vulnerability findings, four cultures,
native-architecture mismatches, malformed/duplicate/unknown clauses and waivers,
UTC/leap expiration, stale or closed issue evidence, scope expansion, multiline
suppression forms, source-map ambiguity/commit/checksum drift, redirect and byte
limits are covered by targeted tests.

**Failures:** unknown audit or approval evidence fails closed; bounded subprocess
and HTTPS operations report errors; the pending waiver prevents a green gate.
No retry switches to unlocked dependency resolution.

**Concurrency:** independent source/test/governance edits used separate ownership;
native CI jobs isolate build outputs. A serial restore resolved the observed local
worker failure. Numerical algorithms and concurrency contracts remain unchanged.

**Security:** exact package/source hosts, policy-specific caches, no shell
interpolation of untrusted values, no secret-bearing provenance and no fabricated
approval. Protected authentication must exist outside PR-controlled code before
active exception evidence is accepted.

The offline governance scanner recognizes ordinary explicitly accessible C# method
bodies and requires the matching restore inside the same method. Unsupported
scopes fail closed. It is not a Roslyn semantic parser, cannot audit arbitrary
outside-repository MSBuild imports, and cannot authenticate review identities.
The provenance tool resolves git/dotnet from explicit operator/runner installation
roots, verifies canonical locations outside checkout/temp storage, checks Unix
ownership/writeability, bounds retained stdout/stderr jointly to 1 MiB, enforces a
30-second deadline and checks SDK identity against global.json. Fixed queries
avoid shell parsing; failures do not echo arbitrary tool output. Eleven targeted
tests exercised command/argument bounds, failed exits, output limits, timeout
cleanup, untrusted locations and SDK mismatch. Windows installation ACLs and its
process-tree cleanup adapter still require hosted qualification. Test harnesses
and canonical Bash commands run in the operator/managed-runner toolchain context;
this audit does not attest the integrity of an arbitrary machine. Durable
ignored artifacts are retained intentionally for review; ephemeral package
consumers and reproduction clones use scoped platform cleanup.

**Disposition:** local code/configuration repairs are reviewable. Full standards
compliance is **not achieved** until the identified approval, protection and hosted
qualification requirements are resolved. Nothing in this report waives them.
