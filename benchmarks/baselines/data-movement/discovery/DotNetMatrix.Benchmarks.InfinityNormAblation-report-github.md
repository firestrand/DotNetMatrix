```

BenchmarkDotNet v0.15.8, Linux Ubuntu 24.04.5 LTS (Noble Numbat)
Cortex-X925, Cortex-A725, 10 physical cores
.NET SDK 10.0.401
  [Host]     : .NET 10.0.12 (10.0.12, 10.0.1226.42308), Arm64 RyuJIT armv8.0-a
  Job-RQTNKT : .NET 10.0.12 (10.0.12, 10.0.1226.42308), Arm64 RyuJIT armv8.0-a

OutlierMode=DontRemove  Affinity=100000  IterationCount=12
IterationTime=100ms  LaunchCount=2  WarmupCount=4

```
| Method         | Shape        | Mean          | Error       | StdDev        | Median        | Ratio | MannWhitney(5%) | RatioSD | Allocated | Alloc Ratio |
|--------------- |------------- |--------------:|------------:|--------------:|--------------:|------:|---------------- |--------:|----------:|------------:|
| **ScalarControl**  | **Tiny**         |     **12.258 ns** |   **1.5816 ns** |     **2.0565 ns** |     **10.911 ns** |  **1.02** | **Baseline**        |    **0.22** |         **-** |          **NA** |
| RowCachingOnly | Tiny         |     14.170 ns |   0.5067 ns |     0.6588 ns |     14.694 ns |  1.18 | Same            |    0.17 |         - |          NA |
| Production     | Tiny         |     13.210 ns |   0.4441 ns |     0.5774 ns |     12.922 ns |  1.10 | Same            |    0.16 |         - |          NA |
|                |              |               |             |               |               |       |                 |         |           |             |
| **ScalarControl**  | **Small**        |  **1,225.437 ns** | **121.6369 ns** |   **158.1623 ns** |  **1,115.940 ns** |  **1.02** | **Baseline**        |    **0.18** |         **-** |          **NA** |
| RowCachingOnly | Small        |  1,043.029 ns | 108.0108 ns |   140.4445 ns |    958.815 ns |  0.86 | Same            |    0.15 |         - |          NA |
| Production     | Small        |  1,170.025 ns | 111.8131 ns |   145.3886 ns |  1,096.794 ns |  0.97 | Same            |    0.17 |         - |          NA |
|                |              |               |             |               |               |       |                 |         |           |             |
| **ScalarControl**  | **Medium**       |  **7,207.053 ns** | **207.9213 ns** |   **270.3565 ns** |  **7,247.570 ns** |  **1.00** | **Baseline**        |    **0.06** |         **-** |          **NA** |
| RowCachingOnly | Medium       |  6,150.468 ns |  24.5795 ns |    31.9603 ns |  6,146.557 ns |  0.85 | Faster          |    0.04 |         - |          NA |
| Production     | Medium       |  7,370.059 ns |  82.7787 ns |   107.6356 ns |  7,395.968 ns |  1.02 | Same            |    0.05 |         - |          NA |
|                |              |               |             |               |               |       |                 |         |           |             |
| **ScalarControl**  | **Large**        | **25,235.351 ns** | **357.2300 ns** |   **464.5000 ns** | **25,325.416 ns** |  **1.00** | **Baseline**        |    **0.03** |         **-** |          **NA** |
| RowCachingOnly | Large        | 20,257.982 ns | 350.1148 ns |   455.2482 ns | 20,264.149 ns |  0.80 | Faster          |    0.02 |         - |          NA |
| Production     | Large        | 25,339.714 ns |  37.8575 ns |    49.2254 ns | 25,345.150 ns |  1.00 | Same            |    0.02 |         - |          NA |
|                |              |               |             |               |               |       |                 |         |           |             |
| **ScalarControl**  | **Tall**         | **45,296.755 ns** |  **44.9941 ns** |    **58.5051 ns** | **45,278.588 ns** |  **1.00** | **Baseline**        |    **0.00** |         **-** |          **NA** |
| RowCachingOnly | Tall         | 36,038.369 ns | 134.8797 ns |   175.3817 ns | 36,050.078 ns |  0.80 | Faster          |    0.00 |         - |          NA |
| Production     | Tall         | 45,259.973 ns |  55.1711 ns |    71.7380 ns | 45,278.352 ns |  1.00 | Same            |    0.00 |         - |          NA |
|                |              |               |             |               |               |       |                 |         |           |             |
| **ScalarControl**  | **Wide**         | **44,837.516 ns** |  **38.7475 ns** |    **50.3827 ns** | **44,834.907 ns** |  **1.00** | **Baseline**        |    **0.00** |         **-** |          **NA** |
| RowCachingOnly | Wide         | 39,281.358 ns |  47.2385 ns |    61.4234 ns | 39,279.251 ns |  0.88 | Faster          |    0.00 |         - |          NA |
| Production     | Wide         | 45,012.381 ns |  50.1235 ns |    65.1747 ns | 45,009.526 ns |  1.00 | Same            |    0.00 |         - |          NA |
|                |              |               |             |               |               |       |                 |         |           |             |
| **ScalarControl**  | **SingleRow**    |  **6,807.696 ns** | **847.8261 ns** | **1,102.4135 ns** |  **7,095.326 ns** |  **1.07** | **Baseline**        |    **0.48** |         **-** |          **NA** |
| RowCachingOnly | SingleRow    |  6,265.262 ns | 691.7537 ns |   899.4753 ns |  6,525.107 ns |  0.98 | Faster          |    0.44 |         - |          NA |
| Production     | SingleRow    |  6,992.835 ns | 909.4171 ns | 1,182.4992 ns |  7,399.225 ns |  1.10 | Same            |    0.50 |         - |          NA |
|                |              |               |             |               |               |       |                 |         |           |             |
| **ScalarControl**  | **SingleColumn** |  **6,111.122 ns** |  **83.0727 ns** |   **108.0180 ns** |  **6,071.651 ns** |  **1.00** | **Baseline**        |    **0.02** |         **-** |          **NA** |
| RowCachingOnly | SingleColumn |  6,688.526 ns |  33.9535 ns |    44.1491 ns |  6,678.867 ns |  1.09 | Slower          |    0.02 |         - |          NA |
| Production     | SingleColumn |  6,549.126 ns |  29.6050 ns |    38.4949 ns |  6,552.097 ns |  1.07 | Slower          |    0.02 |         - |          NA |
|                |              |               |             |               |               |       |                 |         |           |             |
| **ScalarControl**  | **Empty**        |    **153.862 ns** |  **13.7918 ns** |    **17.9332 ns** |    **142.292 ns** |  **1.01** | **Baseline**        |    **0.16** |         **-** |          **NA** |
| RowCachingOnly | Empty        |      2.343 ns |   0.0741 ns |     0.0964 ns |      2.351 ns |  0.02 | Faster          |    0.00 |         - |          NA |
| Production     | Empty        |    157.489 ns |  12.7793 ns |    16.6166 ns |    157.782 ns |  1.04 | Same            |    0.16 |         - |          NA |
