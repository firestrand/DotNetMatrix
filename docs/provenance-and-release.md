# Provenance and release policy

## License and provenance

The maintainer selected the [MIT license](../LICENSE) on 2026-10-02. The package
declares the SPDX expression `MIT` and includes the license text. The canonical
text is published by the [Open Source Initiative](https://opensource.org/license/mit).

The historical [ReadMe.txt](../ReadMe.txt) states this fork came from
`http://www.codeproject.com/KB/recipes/psdotnetmatrix.aspx` and describes it as
“A C# port of public domain Java Matrix library JAMA.” Historical assembly
attributes name Microsoft and copyright 2010. Those are recorded source facts,
not proof of authorship or redistribution rights. A public-domain claim about JAMA does not establish the
license of a C# port, contributions or this fork.

Historical attribution and ownership attributes are retained. The maintainer's
MIT selection does not establish new facts about the original port's authorship
or permissions; upstream provenance remains documented for release review.

`DotNetMatrix.LocalPreview` version `2.0.0-preview.1` is a local verification ID.
The repository URL identifies this checkout's remote; the contributor label is
descriptive and makes no claim of a registered package namespace. Assembly
identity remains historical until the maintainer approves an identity migration.
The package target is net10.0; no untested additional frameworks are advertised.

## Local release gate

1. Run `bash scripts/verify.sh` from a clean checkout with the pinned SDK.
2. Review mathematical behavior changes and the changelog, not only coverage.
3. Review `docs/public-api.txt` against the intended contract. The reflection
   snapshot captures exported types, inheritance/interfaces, member signatures,
   accessibility of accessors, virtual/static flags and parameter defaults. It
   supplements tests; it is not a proof of numeric semantics or binary identity.
4. For an intentional API change, build the library, then regenerate the snapshot
   with the command below. Inspect and commit the diff with the corresponding
   tests/release notes; never regenerate merely to silence the gate.
5. Inspect the local nupkg, XML docs, README and independent consumer evidence.
6. Verify hosted CI on Linux x64/arm64, Windows and macOS before claiming those
   platforms passed. Local Linux success does not prove hosted results.
7. Review upstream provenance, public package identity, version/assembly policy and
   distribution approval before any publishing command. Publication is separate
   explicit authorization.

```bash
dotnet restore tools/ApiBaseline/ApiBaseline.csproj --locked-mode
dotnet build tools/ApiBaseline/ApiBaseline.csproj -c Release --no-restore
dotnet tools/ApiBaseline/bin/Release/net10.0/ApiBaseline.dll api \
  DotNetMatrix/bin/Release/net10.0/DotNetMatrix.dll > docs/public-api.txt
```

The production coverage inventory is derived from nonhidden portable-PDB
sequence points of the actual built assembly. New executable classes therefore
must appear in the report; a stale hardcoded class list cannot hide them.
Compiler-generated nested methods are attributed to their containing class.
Python verifier tests reject exact-80% rates, missing/duplicate modules or
reports, failed/skipped/empty suites and malformed coverage counts.

Generated XML documentation suppresses only CS1591 for historically undocumented
public members. Other warnings remain errors; filling the historical comment
gaps remains documentation work rather than a relaxed compiler gate.

The Linux arm64 job uses the documented standard `ubuntu-24.04-arm` runner;
see [GitHub's hosted runner reference](https://docs.github.com/en/actions/reference/runners/github-hosted-runners).
Action versions were checked against their official release tags and pinned by
commit SHA. Benchmark runs are observational and have no flaky timing threshold.

Packaging follows [Microsoft's package authoring guidance](https://learn.microsoft.com/en-us/nuget/create-packages/package-authoring-best-practices).
The project license is selected; package publication remains a separate release action.
