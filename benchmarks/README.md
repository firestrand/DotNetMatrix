# Dense matrix benchmarks

Run from the repository root with the pinned .NET SDK:

```bash
dotnet restore benchmarks/DotNetMatrix.Benchmarks.csproj --locked-mode
dotnet build benchmarks/DotNetMatrix.Benchmarks.csproj -c Release --no-restore
dotnet run --project benchmarks/DotNetMatrix.Benchmarks.csproj -c Release --no-build -- --validate
dotnet run --project benchmarks/DotNetMatrix.Benchmarks.csproj -c Release --no-build -- --filter '*' --artifacts artifacts/benchmarks
```

The default BenchmarkDotNet run takes longer and provides more measurements than
our checked-in exploratory short baseline. To reproduce that baseline's job:

```bash
dotnet run --project benchmarks/DotNetMatrix.Benchmarks.csproj -c Release --no-build -- --filter '*MultiplicationCandidates*' --job short --warmupCount 1 --iterationCount 3 --launchCount 1 --iterationTime 100 --artifacts artifacts/benchmarks
dotnet run --project benchmarks/DotNetMatrix.Benchmarks.csproj -c Release --no-build -- --filter '*MatrixOperations*' '*Factorization*' --job short --warmupCount 1 --iterationCount 3 --launchCount 1 --iterationTime 100 --artifacts artifacts/benchmarks-workloads
```

See [performance results and interpretation](../docs/performance.md), including
hardware, runtime, uncertainty, numerical checks, and allocation tradeoffs.
BenchmarkDotNet is confined to this nonpackable development project; the matrix
library does not reference it. Its generated benchmark projects need access to
its compile/runtime assets, so this project's package reference intentionally
has no `PrivateAssets="all"` restriction.
