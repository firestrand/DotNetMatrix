# Adopted .NET standard and project exception policy

The repository uses the supplied **ENG-STD-DOTNET-XPLAT 1.1.0-draft.1** as its
technical standard. [engineering-standard.md](engineering-standard.md) remains
an exact copy of that draft; this repository does not claim organization-wide
approval or a corporate effective date.

[dotnet-rules.json](dotnet-rules.json) records all 43 rule/group IDs, exact source
ranges and enforcement references. The project policy below records the
maintainer-approved exception process.

## Maintainer approval

At the repository maintainer's explicit direction, this project accepts recorded
maintainer approval for exceptions. This replaces the draft's independent-review,
protected evidence, issue-status and mandatory-expiration process for this
repository. No GitHub review, remediation issue, branch-protection setting or
external approval artifact is required to activate a maintainer-approved record.

The registry is [overrides.yaml](overrides.yaml), using JSON syntax (a YAML 1.2
subset). Active entries require a recorded approver, rationale, owner,
security/data-integrity assessment, compensating controls, known rule and exact
clause, diagnostic, file and symbol scope. Pending and retired entries do not
approve a live suppression. Unregistered or expanded suppressions still fail.

`expires_on: null` records an accepted ongoing exception. If a maintainer specifies
an expiration date, the gate enforces it exclusively at 00:00 UTC on that date.
Malformed dates and references, duplicate IDs and source/catalog drift still fail.
The linter checks recorded approval and scope; it does not authenticate a human
identity or replace maintainers' responsibility for registry changes.

[EXC-2026-001](serialization-exception.md) is approved by the repository maintainer
for the existing legacy serialization regression test. Its original assertions
and narrow `SYSLIB0050` pragma remain intact.

Run the same offline gate locally and in CI:

```bash
python3 scripts/check-standards.py
python3 -m unittest discover -s scripts/tests -p 'test_standards.py' -v
bash scripts/verify.sh
```

This project policy changes exception administration. It does not disable
warnings-as-errors, analyzers, dependency audits, test assertions, coverage gates
or public API checks.
