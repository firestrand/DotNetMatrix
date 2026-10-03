# Adopted .NET standard and governance controls

The user adopted the exact supplied **ENG-STD-DOTNET-XPLAT 1.1.0-draft.1** for this
repository audit. [engineering-standard.md](engineering-standard.md) is an exact
copy, including its proposed/not-effective corporate status. This repository does
not declare a global effective date, organization approval or corporate SDK policy.

[dotnet-rules.json](dotnet-rules.json) records all 43 declared rule/group IDs, exact
source ranges, owner roles, clause selectors and required enforcement contracts.
A group owns its range until the next group or top-level section; individually
numbered rules retain precedence within that range. Enforcement names describe
requirements, not evidence of successful execution. [audit.md](audit.md) records
applicability, verification and external controls separately.

The registry uses **JSON syntax, a YAML 1.2 subset**, so Python's standard-library
JSON parser can reject duplicate fields without adding a dependency. General YAML
features are intentionally unsupported. Pending records are proposals; retired
records preserve history. Only independently verified active records can authorize
an exact matching suppression. The current serialization proposal is unapproved.

Run offline lint from any working directory:

```bash
python3 scripts/check-standards.py
python3 -m unittest discover -s scripts/tests -p 'test_standards.py' -v
```

The linter checks source/catalog consistency, known IDs/clauses, registry fields,
warning suppression markers, scoped diagnostics, independent authorization,
expiration and fresh issue evidence. The full gate deliberately fails on the
pending legacy serialization suppression. It must not be skipped to claim full
standards compliance.

For merge/release authorization, a **separately protected trusted job** must verify
required owner review and designated independent approver identity from real review
records, query the HTTPS issue using credentials unavailable to untrusted PR code,
and produce integrity-protected evidence bound to the registry revision. The job
passes the authenticated artifact outside the checkout plus its expected digest:

```bash
python3 scripts/check-standards.py --release \
  --trusted-evidence "$TRUSTED_EVIDENCE_PATH" \
  --evidence-sha256 "$TRUSTED_EVIDENCE_SHA256"
```

Those inputs are permitted only from independently protected workflow/policy;
**neither an outside-checkout path nor a SHA-256 alone proves authenticity**. This
CLI verifies binding/age/claims after that external authentication. It cannot prove
that GitHub branch protection, approver authority or workflow ownership exists.
A PR cannot supply its own expected digest or invoke its own modified gate as the
trusted authority. Missing trusted setup is an unmet external control.

Evidence shape (illustrative only, not real approval): top-level
`registry_sha256`, `owner_review_verified`, and `exceptions` keyed by exception ID.
Each exception requires `record_sha256` (SHA-256 of UTF-8 sorted-key JSON encoding
of the registry record), `approved_by`, `approval_verified`, `issue`, `issue_state`
(`open`) and `checked_at` (explicit UTC timestamp). The artifact must come from
protected authenticated processing. Missing/API-failed/closed/unknown issue state,
future or older-than-24-hour evidence, record scope drift and unverified reviews
fail closed. Unit tests simulate these claims; they are not trusted approvals.

Owners must protect the catalog, registry, analyzer configurations, gates and
trusted jobs using required independent review and branch/ruleset controls; add a
daily scheduled trusted check and owner notification for expirations/closed issues.
Repository configuration alone does not enforce those settings. No active exception
can be accepted while those prerequisites remain unavailable. Expiration is at
00:00 UTC on `expires_on`; review due is never automatic renewal. Existing
production artifacts follow a separate risk/remediation process.
