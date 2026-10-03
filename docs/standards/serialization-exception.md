# Approved legacy serialization test exception

**Status: active. Approver: repository maintainer, by explicit operator approval.**
The maintainer accepted this exception and subsequently directed removal of the
additional independent-review, issue and protected-evidence requirements. That
recorded decision is sufficient under the [project exception policy](README.md).
No further approval step or artificial expiration applies.

## Accepted scope

`EXC-2026-001` permits only `SYSLIB0050` in
`DotNetMatrix_Test.ModernMatrixTests.LegacySerializationContractEmitsNoPayload`,
in `DotNetMatrix_Test/ModernMatrixTests.cs`. The exact rule, warning-policy clause,
compensating controls and active approval are recorded in
[overrides.yaml](overrides.yaml).

The test constructs `SerializationInfo` with `FormatterConverter` and invokes
`ISerializable.GetObjectData` directly. Both original assertions remain: zero
serialized members and the original matrix object type. This preserves the
legacy public contract without invoking a binary formatter or deserializer.

The existing oracle technical review approved this narrow test-only use and
confirmed that the production callback remains empty. The responsible maintainer
has accepted the security assessment: metadata construction only, with no external
input, formatter/deserialization operation, runtime enablement, added dependency
or production change.

## Ongoing controls

The exception remains valid while this compatibility test is supported. It does
not authorize broader suppressions, other files/symbols/diagnostics or changes to
its assertions. The gate continues to reject unregistered and out-of-scope
suppressions. A future contract migration must intentionally update the tests and
retire the unused exception record.

The original draft standard remains unchanged. Independent corporate approval
and protected governance records are not claimed or fabricated; the maintainer
explicitly chose the simpler project approval policy.
