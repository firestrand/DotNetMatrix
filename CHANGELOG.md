# Changelog

## 2.0.0-preview.1 — local correctness preview

Intentional behavior corrections require a major consumer release:

- Cholesky forward substitution divides before eliminating following rows;
  nonunit diagonal factors now solve the stated equations.
- `SolveTranspose` returns the documented solution orientation.
- Equality handles null operands, exact runtime types, NaNs and signed zeros
  consistently; hashing uses contents. Mutable keys must remain unmodified.
- Decompositions validate finite nonempty inputs and supported shapes; wide SVD
  provides coherent economy factors. Iterative solvers fail on a bounded budget.

Additive capabilities include rank-cutoff pseudoinverse/minimum-norm solving,
optional numerical diagnostics, factor-copy accessors, matrix conveniences and
versioned JSON persistence. Legacy storage aliases and the empty serialization
callback remain. Ordinary `Solve`/`Inverse` keep their algorithm selection.

Release support adds locked reproducible builds, public API snapshots,
multi-platform CI, strict complete coverage gates and independent package smoke.
The repository is licensed under MIT, with the license included in the package
and declared in NuGet metadata. This package remains a local preview.

## Earlier modernization

The original .NET Framework 4.0 code was moved to .NET 10 / C# 14 with SDK-style
projects and modern MSTest/Microsoft.Testing.Platform tooling. The historical
modernization baseline had 100 passing tests.
