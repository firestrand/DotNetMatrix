```

BenchmarkDotNet v0.15.8, Linux Ubuntu 24.04.5 LTS (Noble Numbat)
Cortex-X925, Cortex-A725, 10 physical cores
.NET SDK 10.0.401
  [Host]     : .NET 10.0.12 (10.0.12, 10.0.1226.42308), Arm64 RyuJIT armv8.0-a
  Job-RRWQJM : .NET 10.0.12 (10.0.12, 10.0.1226.42308), Arm64 RyuJIT armv8.0-a

OutlierMode=DontRemove  IterationCount=12  IterationTime=100ms
LaunchCount=2  WarmupCount=4

```
| Method                | Shape        | Mean         | Error        | StdDev       | Median       | Ratio | MannWhitney(5%) | RatioSD | Gen0     | Gen1    | Allocated | Alloc Ratio |
|---------------------- |------------- |-------------:|-------------:|-------------:|-------------:|------:|---------------- |--------:|---------:|--------:|----------:|------------:|
| **ScalarControl**         | **Tiny**         |     **42.92 ns** |     **4.619 ns** |     **6.006 ns** |     **40.65 ns** |  **1.01** | **Baseline**        |    **0.17** |   **0.0229** |       **-** |      **96 B** |        **1.00** |
| RowCachingOnly        | Tiny         |     43.12 ns |     3.091 ns |     4.019 ns |     42.33 ns |  1.02 | Same            |    0.14 |   0.0226 |       - |      96 B |        1.00 |
| RowCachingAndBulkCopy | Tiny         |     42.15 ns |     1.631 ns |     2.120 ns |     41.71 ns |  1.00 | Same            |    0.11 |   0.0228 |       - |      96 B |        1.00 |
| GuardedBulkCopy       | Tiny         |     40.54 ns |     0.125 ns |     0.162 ns |     40.50 ns |  0.96 | Same            |    0.10 |   0.0227 |       - |      96 B |        1.00 |
| Production            | Tiny         |     42.90 ns |     0.178 ns |     0.231 ns |     42.83 ns |  1.01 | Same            |    0.10 |   0.0229 |       - |      96 B |        1.00 |
|                       |              |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**         | **Small**        |  **1,625.50 ns** |    **85.246 ns** |   **110.844 ns** |  **1,596.78 ns** |  **1.00** | **Baseline**        |    **0.09** |   **0.6221** |       **-** |    **2616 B** |        **1.00** |
| RowCachingOnly        | Small        |  1,329.23 ns |    19.496 ns |    25.351 ns |  1,324.51 ns |  0.82 | Faster          |    0.05 |   0.6213 |       - |    2616 B |        1.00 |
| RowCachingAndBulkCopy | Small        |    377.15 ns |     8.594 ns |    11.174 ns |    373.94 ns |  0.23 | Faster          |    0.01 |   0.6243 |       - |    2616 B |        1.00 |
| GuardedBulkCopy       | Small        |    371.14 ns |     2.066 ns |     2.687 ns |    370.59 ns |  0.23 | Faster          |    0.01 |   0.6232 |       - |    2616 B |        1.00 |
| Production            | Small        |  1,612.38 ns |    18.726 ns |    24.349 ns |  1,606.76 ns |  1.00 | Same            |    0.05 |   0.6103 |       - |    2616 B |        1.00 |
|                       |              |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**         | **Medium**       | **10,634.24 ns** |    **87.341 ns** |   **113.569 ns** | **10,616.91 ns** |  **1.00** | **Baseline**        |    **0.01** |   **8.2767** |       **-** |   **34872 B** |        **1.00** |
| RowCachingOnly        | Medium       |  7,808.75 ns |   130.572 ns |   169.780 ns |  7,771.48 ns |  0.73 | Faster          |    0.02 |   8.3281 |       - |   34872 B |        1.00 |
| RowCachingAndBulkCopy | Medium       |  2,383.53 ns |    26.239 ns |    34.118 ns |  2,369.49 ns |  0.22 | Faster          |    0.00 |   8.3238 |       - |   34872 B |        1.00 |
| GuardedBulkCopy       | Medium       |  2,319.75 ns |   163.880 ns |   213.091 ns |  2,371.47 ns |  0.22 | Faster          |    0.02 |   8.3381 |       - |   34872 B |        1.00 |
| Production            | Medium       | 10,568.24 ns |    82.018 ns |   106.647 ns | 10,613.28 ns |  0.99 | Same            |    0.01 |   8.2974 |       - |   34872 B |        1.00 |
|                       |              |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**         | **Large**        | **45,463.19 ns** |   **215.906 ns** |   **280.738 ns** | **45,366.89 ns** |  **1.00** | **Baseline**        |    **0.01** |  **64.3116** | **22.1920** |  **301112 B** |        **1.00** |
| RowCachingOnly        | Large        | 35,357.00 ns |   215.494 ns |   280.203 ns | 35,320.78 ns |  0.78 | Faster          |    0.01 |  64.2556 | 22.1208 |  301112 B |        1.00 |
| RowCachingAndBulkCopy | Large        | 16,217.99 ns |   515.198 ns |   669.903 ns | 16,324.42 ns |  0.36 | Faster          |    0.01 |  64.1234 | 22.5649 |  301112 B |        1.00 |
| GuardedBulkCopy       | Large        | 16,346.04 ns |   139.317 ns |   181.151 ns | 16,343.52 ns |  0.36 | Faster          |    0.00 |  64.1404 | 22.6378 |  301112 B |        1.00 |
| Production            | Large        | 45,239.58 ns |   291.822 ns |   379.451 ns | 45,385.25 ns |  1.00 | Same            |    0.01 |  64.3116 | 22.1920 |  301112 B |        1.00 |
|                       |              |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**         | **Tall**         | **87,213.79 ns** |   **583.610 ns** |   **758.858 ns** | **87,095.30 ns** |  **1.00** | **Baseline**        |    **0.01** | **118.0556** | **54.6875** |  **602168 B** |        **1.00** |
| RowCachingOnly        | Tall         | 67,721.14 ns |   298.535 ns |   388.180 ns | 67,634.72 ns |  0.78 | Faster          |    0.01 | 118.2065 | 54.3478 |  602168 B |        1.00 |
| RowCachingAndBulkCopy | Tall         | 33,722.02 ns |   162.248 ns |   210.968 ns | 33,754.42 ns |  0.39 | Faster          |    0.00 | 119.9597 | 53.0914 |  602168 B |        1.00 |
| GuardedBulkCopy       | Tall         | 33,947.73 ns |   321.867 ns |   418.519 ns | 34,092.12 ns |  0.39 | Faster          |    0.01 | 118.9840 | 54.1444 |  602168 B |        1.00 |
| Production            | Tall         | 87,644.18 ns |   258.269 ns |   335.822 ns | 87,653.34 ns |  1.01 | Same            |    0.01 | 117.9577 | 55.4577 |  602168 B |        1.00 |
|                       |              |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**         | **Wide**         | **86,324.97 ns** |   **769.401 ns** | **1,000.439 ns** | **86,507.96 ns** |  **1.00** | **Baseline**        |    **0.02** | **114.8649** | **75.1689** |  **596024 B** |        **1.00** |
| RowCachingOnly        | Wide         | 66,989.05 ns |   376.228 ns |   489.202 ns | 66,997.37 ns |  0.78 | Faster          |    0.01 | 115.0266 | 73.1383 |  596024 B |        1.00 |
| RowCachingAndBulkCopy | Wide         | 31,777.61 ns |   190.809 ns |   248.105 ns | 31,759.61 ns |  0.37 | Faster          |    0.01 | 115.2638 | 74.7487 |  596024 B |        1.00 |
| GuardedBulkCopy       | Wide         | 31,566.75 ns |   323.073 ns |   420.086 ns | 31,584.75 ns |  0.37 | Faster          |    0.01 | 115.3607 | 74.9378 |  596024 B |        1.00 |
| Production            | Wide         | 87,892.40 ns | 3,748.926 ns | 4,874.663 ns | 86,116.71 ns |  1.02 | Same            |    0.06 | 114.8649 | 75.1689 |  596024 B |        1.00 |
|                       |              |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**         | **SingleRow**    | **11,176.72 ns** |   **340.593 ns** |   **442.867 ns** | **10,999.25 ns** |  **1.00** | **Baseline**        |    **0.05** |   **7.7977** |       **-** |   **32856 B** |        **1.00** |
| RowCachingOnly        | SingleRow    |  7,696.20 ns |   328.733 ns |   427.446 ns |  7,740.16 ns |  0.69 | Faster          |    0.05 |   7.8028 |       - |   32856 B |        1.00 |
| RowCachingAndBulkCopy | SingleRow    |    981.43 ns |    10.267 ns |    13.350 ns |    980.63 ns |  0.09 | Faster          |    0.00 |   7.8040 |       - |   32856 B |        1.00 |
| GuardedBulkCopy       | SingleRow    |  1,025.03 ns |    59.339 ns |    77.157 ns |  1,001.77 ns |  0.09 | Faster          |    0.01 |   7.8112 |       - |   32856 B |        1.00 |
| Production            | SingleRow    |  9,966.74 ns |    69.085 ns |    89.830 ns |  9,949.19 ns |  0.89 | Faster          |    0.03 |   7.7751 |       - |   32856 B |        1.00 |
|                       |              |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**         | **SingleColumn** | **36,773.26 ns** |   **874.791 ns** | **1,137.475 ns** | **36,926.43 ns** |  **1.00** | **Baseline**        |    **0.04** |  **37.7095** | **12.2207** |  **163896 B** |        **1.00** |
| RowCachingOnly        | SingleColumn | 35,384.32 ns |   322.690 ns |   419.588 ns | 35,306.95 ns |  0.96 | Same            |    0.03 |  37.9213 | 12.2893 |  163896 B |        1.00 |
| RowCachingAndBulkCopy | SingleColumn | 47,310.02 ns | 3,037.712 ns | 3,949.884 ns | 46,995.53 ns |  1.29 | Slower          |    0.11 |  37.5000 | 12.0000 |  163896 B |        1.00 |
| GuardedBulkCopy       | SingleColumn | 37,416.18 ns | 2,422.256 ns | 3,149.617 ns | 36,941.16 ns |  1.02 | Same            |    0.09 |  37.7994 | 12.3503 |  163896 B |        1.00 |
| Production            | SingleColumn | 37,138.36 ns | 2,017.638 ns | 2,623.500 ns | 36,838.27 ns |  1.01 | Same            |    0.08 |  37.9464 | 12.6488 |  163896 B |        1.00 |
|                       |              |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**         | **Empty**        |    **441.96 ns** |    **27.029 ns** |    **35.145 ns** |    **426.04 ns** |  **1.01** | **Baseline**        |    **0.10** |   **0.2581** |       **-** |    **1080 B** |        **1.00** |
| RowCachingOnly        | Empty        |    297.00 ns |     2.646 ns |     3.441 ns |    296.05 ns |  0.68 | Faster          |    0.05 |   0.2576 |       - |    1080 B |        1.00 |
| RowCachingAndBulkCopy | Empty        |    297.17 ns |     2.498 ns |     3.249 ns |    297.39 ns |  0.68 | Faster          |    0.05 |   0.2582 |       - |    1080 B |        1.00 |
| GuardedBulkCopy       | Empty        |    490.68 ns |    21.059 ns |    27.383 ns |    484.60 ns |  1.12 | Same            |    0.10 |   0.2532 |       - |    1080 B |        1.00 |
| Production            | Empty        |    491.24 ns |    12.589 ns |    16.370 ns |    486.27 ns |  1.12 | Same            |    0.08 |   0.2573 |       - |    1080 B |        1.00 |
|                       |              |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**         | **TwoColumns**   | **40,316.74 ns** |   **293.256 ns** |   **381.316 ns** | **40,327.52 ns** |  **1.00** | **Baseline**        |    **0.01** |  **45.7803** | **11.9427** |  **196664 B** |        **1.00** |
| RowCachingOnly        | TwoColumns   | 43,456.20 ns | 3,271.085 ns | 4,253.335 ns | 43,292.10 ns |  1.08 | Same            |    0.10 |  45.8333 | 12.0833 |  196664 B |        1.00 |
| RowCachingAndBulkCopy | TwoColumns   | 46,444.62 ns |   489.156 ns |   636.041 ns | 46,338.26 ns |  1.15 | Slower          |    0.02 |  45.4545 | 11.8371 |  196664 B |        1.00 |
| GuardedBulkCopy       | TwoColumns   | 48,053.78 ns | 1,067.559 ns | 1,388.129 ns | 47,618.73 ns |  1.19 | Slower          |    0.04 |  45.6349 | 11.9048 |  196664 B |        1.00 |
| Production            | TwoColumns   | 39,918.79 ns |   444.791 ns |   578.354 ns | 39,743.66 ns |  0.99 | Same            |    0.02 |  45.5975 | 12.1855 |  196664 B |        1.00 |
|                       |              |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**         | **FourColumns**  | **54,738.90 ns** | **4,316.640 ns** | **5,612.851 ns** | **51,062.33 ns** |  **1.01** | **Baseline**        |    **0.14** |  **60.3070** | **19.7368** |  **262200 B** |        **1.00** |
| RowCachingOnly        | FourColumns  | 48,227.45 ns | 3,107.165 ns | 4,040.192 ns | 45,959.59 ns |  0.89 | Same            |    0.11 |  60.6061 | 19.8864 |  262200 B |        1.00 |
| RowCachingAndBulkCopy | FourColumns  | 48,909.11 ns |   734.470 ns |   955.019 ns | 48,959.54 ns |  0.90 | Same            |    0.09 |  60.5916 | 20.0382 |  262200 B |        1.00 |
| GuardedBulkCopy       | FourColumns  | 52,761.94 ns | 3,759.485 ns | 4,888.393 ns | 50,600.38 ns |  0.97 | Same            |    0.13 |  60.5620 | 19.8643 |  262200 B |        1.00 |
| Production            | FourColumns  | 50,925.90 ns | 1,019.079 ns | 1,325.091 ns | 50,848.08 ns |  0.94 | Same            |    0.09 |  60.4675 | 19.8171 |  262200 B |        1.00 |
|                       |              |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**         | **EightColumns** | **69,003.27 ns** |   **441.526 ns** |   **574.108 ns** | **69,019.30 ns** |  **1.00** | **Baseline**        |    **0.01** |  **75.6944** | **41.6667** |  **393272 B** |        **1.00** |
| RowCachingOnly        | EightColumns | 71,518.38 ns | 3,297.022 ns | 4,287.060 ns | 70,877.21 ns |  1.04 | Same            |    0.06 |  75.8333 | 44.1667 |  393272 B |        1.00 |
| RowCachingAndBulkCopy | EightColumns | 71,623.88 ns | 2,791.269 ns | 3,629.439 ns | 70,434.53 ns |  1.04 | Same            |    0.05 |  76.3889 | 40.9722 |  393272 B |        1.00 |
| GuardedBulkCopy       | EightColumns | 59,735.66 ns | 4,276.538 ns | 5,560.708 ns | 57,029.78 ns |  0.87 | Same            |    0.08 |  76.1905 | 42.2619 |  393272 B |        1.00 |
| Production            | EightColumns | 70,533.97 ns | 2,333.352 ns | 3,034.018 ns | 69,479.57 ns |  1.02 | Same            |    0.04 |  75.8427 | 41.4326 |  393272 B |        1.00 |
