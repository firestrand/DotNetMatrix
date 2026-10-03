# Pending legacy serialization test exception

**Status: oracle technical approval and operator approval received; protected
governance activation pending.** `EXC-2026-001` does not yet authorize a suppression
under the protected gate. The gate must fail until the required independent review,
open remediation issue and authenticated evidence are available.

The adopted `DOTNET-REPO-004` warning-policy clause requires configured compiler
warnings to fail official builds unless an approved exception applies. The exact
source clause and test file/symbol are recorded in [overrides.yaml](overrides.yaml).
The existing `SYSLIB0050` pragma preserves the original direct regression assertion
that the legacy `ISerializable` callback emits zero members and retains the matrix
object type. The test constructs metadata; it invokes no binary formatter.

Proposed owner: repository maintainer. Proposed exclusive deadline: **2026-11-01
00:00 UTC**. The maintainer must identify an independent responsible approver and
an open HTTPS remediation/review issue. Security review is required if the
responsible owner identifies security implications. Neither identity nor approval
has been fabricated. Until approval is recorded and verified through protected
review records, the waiver remains invalid.

The assertions must remain while the legacy contract remains supported. Only a
separately authorized, compatibility-reviewed future contract migration may retire
this test. Changing or removing assertions merely to make the standards gate pass
is prohibited. The existing assertions and narrow pragma are currently retained.

The proposed exception permits no production changes, formatter or deserialization
invocation, runtime enablement, dependency additions or suppression scope expansion.

Approval requires a reviewed registry change, exact file/symbol/diagnostic scope,
security/data-integrity assessment and compensating controls. A trusted job must
verify review authorization and the issue's open state, bind evidence to the
registry and record hashes, and enforce the exclusive UTC expiration and maximum
24-hour issue-evidence age. An agent-written `approved_by` field is not approval.

## Oracle technical review

**Verdict: APPROVE the technical exception.** Reviewer:
`/root/serialization_exception_oracle`, independently spawned at the maintainer's
explicit request. Review recorded on **2026-10-03 UTC**. This is an agent technical
review record, not an authenticated protected governance approval.

The oracle verified that both original assertions match HEAD, the production
`ISerializable.GetObjectData` callback remains empty, and the pragma covers only
metadata construction and the direct callback. Assertions remain outside the
suppression. The focused test passed **1/1**. The test uses no formatter,
deserialization, external input, runtime enablement or added dependency.

The oracle approved only `SYSLIB0050` within
`ModernMatrixTests.LegacySerializationContractEmitsNoPayload`, through the exclusive
deadline **2026-11-01 00:00 UTC**. Reflection-based constructor replacement would
conceal the diagnostic; JSON-only assertions would change the verification intent.
Neither is an appropriate substitute for this narrow compatibility test.

The oracle explicitly recommended keeping the registry **pending**: its review does
not establish a responsible approval identity, protected review records, an open
remediation issue or fresh authenticated evidence. No `approved_by`, issue or
active status has been fabricated. Formal activation remains blocked by those
requirements under `DOTNET-GOV-004`.

## Operator approval and activation

The trusted operator subsequently replied **“Approved proceed”**, approving this
exception and continued activation work. Approval was received on **2026-10-03
UTC**. The approved scope and exclusive expiration remain exactly those reviewed
by the oracle; no production changes or broader suppression are authorized.

This approval is recorded faithfully as a session decision. It is not represented
as a GitHub review, a named independent review identity or an authenticated CI
artifact. The registry remains pending until those independently verifiable
requirements can be fulfilled.

A read-only repository check found no open remediation issue or pull request to
attach the exception to. Activation therefore still needs:

1. The independent reviewer's GitHub identity and an open HTTPS remediation issue.
2. Protected required owner review for the catalog, registry, analyzer configuration
   and gates, with real approval records covering the exact registry revision.
3. A protected authenticated job that verifies those records and the issue's open
   state, then supplies revision-bound evidence no older than 24 hours.

The approval removes the need to reconfirm the exception's technical scope. It
does not authorize substituting invented review/issue data or weakening the gate.
