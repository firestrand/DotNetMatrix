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

## Data movement ablations

See [the ablation study](../docs/performance-ablation.md) and its checked-in raw
reports under `baselines/data-movement/`. The harness retains the original scalar
loop, row caching alone, unconditional bulk copying, guarded bulk copying and the
actual production entry point. Setup checks copied bits, storage ownership and
norm equivalence before timing. Square, rectangular, narrow and empty fixtures
are included; the library has no benchmark dependency.

```bash
dotnet build benchmarks/DotNetMatrix.Benchmarks.csproj -c Release --no-restore -p:ContinuousIntegrationBuild=true -warnaserror
taskset -c 5 dotnet benchmarks/bin/Release/net10.0/DotNetMatrix.Benchmarks.dll --filter '*CopyAblation*' --warmupCount 4 --iterationCount 12 --launchCount 2 --iterationTime 100 --outliers DontRemove --statisticalTest 5% --exporters JSON CSV --artifacts artifacts/performance-ablation/reproduction
```

`taskset` is Linux-specific: choose an available CPU on your host, or omit it and
record the scheduling limitation. CPU 5 was used for the checked-in study. The
discovery run additionally requested `--affinity 32`; that BenchmarkDotNet setting
failed, while inherited `taskset` affinity was verified. For the separate norm
experiment, filter `'*InfinityNormAblation*'` with the same job settings. Original
production measurements were taken before changing `GeneralMatrix.Copy`; the
frozen scalar control remains runnable after the optimization.
