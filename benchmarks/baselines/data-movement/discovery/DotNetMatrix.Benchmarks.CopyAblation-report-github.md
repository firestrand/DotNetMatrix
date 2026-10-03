```

BenchmarkDotNet v0.15.8, Linux Ubuntu 24.04.5 LTS (Noble Numbat)
Cortex-X925, Cortex-A725, 10 physical cores
.NET SDK 10.0.401
  [Host]     : .NET 10.0.12 (10.0.12, 10.0.1226.42308), Arm64 RyuJIT armv8.0-a
  Job-RQTNKT : .NET 10.0.12 (10.0.12, 10.0.1226.42308), Arm64 RyuJIT armv8.0-a

OutlierMode=DontRemove  Affinity=100000  IterationCount=12
IterationTime=100ms  LaunchCount=2  WarmupCount=4

```
| Method                | Shape        | Mean         | Error         | StdDev        | Ratio | MannWhitney(5%) | RatioSD | Gen0     | Gen1    | Allocated | Alloc Ratio |
|---------------------- |------------- |-------------:|--------------:|--------------:|------:|---------------- |--------:|---------:|--------:|----------:|------------:|
| **ScalarControl**         | **Tiny**         |     **43.99 ns** |      **5.154 ns** |      **6.702 ns** |  **1.02** | **Baseline**        |    **0.19** |   **0.0229** |       **-** |      **96 B** |        **1.00** |
| RowCachingOnly        | Tiny         |     43.76 ns |      2.937 ns |      3.818 ns |  1.01 | Same            |    0.14 |   0.0228 |       - |      96 B |        1.00 |
| RowCachingAndBulkCopy | Tiny         |     43.06 ns |      2.521 ns |      3.277 ns |  1.00 | Same            |    0.14 |   0.0229 |       - |      96 B |        1.00 |
| Production            | Tiny         |     45.04 ns |      2.741 ns |      3.563 ns |  1.04 | Same            |    0.14 |   0.0228 |       - |      96 B |        1.00 |
|                       |              |              |               |               |       |                 |         |          |         |           |             |
| **ScalarControl**         | **Small**        |  **1,726.49 ns** |    **133.947 ns** |    **174.169 ns** |  **1.01** | **Baseline**        |    **0.13** |   **0.6121** |       **-** |    **2616 B** |        **1.00** |
| RowCachingOnly        | Small        |  1,343.66 ns |     10.916 ns |     14.194 ns |  0.78 | Faster          |    0.06 |   0.6196 |       - |    2616 B |        1.00 |
| RowCachingAndBulkCopy | Small        |    376.32 ns |      1.616 ns |      2.101 ns |  0.22 | Faster          |    0.02 |   0.6250 |       - |    2616 B |        1.00 |
| Production            | Small        |  1,659.72 ns |     53.898 ns |     70.083 ns |  0.97 | Same            |    0.09 |   0.6121 |       - |    2616 B |        1.00 |
|                       |              |              |               |               |       |                 |         |          |         |           |             |
| **ScalarControl**         | **Medium**       |  **9,865.93 ns** |     **39.949 ns** |     **51.945 ns** |  **1.00** | **Baseline**        |    **0.01** |   **8.3070** |       **-** |   **34872 B** |        **1.00** |
| RowCachingOnly        | Medium       |  8,370.09 ns |     78.557 ns |    102.147 ns |  0.85 | Faster          |    0.01 |   8.3333 |       - |   34872 B |        1.00 |
| RowCachingAndBulkCopy | Medium       |  2,372.35 ns |     14.918 ns |     19.398 ns |  0.24 | Faster          |    0.00 |   8.3381 |       - |   34872 B |        1.00 |
| Production            | Medium       |  9,861.32 ns |     69.924 ns |     90.921 ns |  1.00 | Same            |    0.01 |   8.2418 |       - |   34872 B |        1.00 |
|                       |              |              |               |               |       |                 |         |          |         |           |             |
| **ScalarControl**         | **Large**        | **52,522.35 ns** | **10,970.924 ns** | **14,265.302 ns** |  **1.04** | **Baseline**        |    **0.33** |  **64.5161** | **24.1935** |  **301112 B** |        **1.00** |
| RowCachingOnly        | Large        | 38,192.76 ns |  1,705.249 ns |  2,217.306 ns |  0.76 | Faster          |    0.13 |  63.9451 | 22.0376 |  301112 B |        1.00 |
| RowCachingAndBulkCopy | Large        | 16,310.04 ns |    208.164 ns |    270.671 ns |  0.32 | Faster          |    0.05 |  64.1447 | 22.5329 |  301112 B |        1.00 |
| Production            | Large        | 47,591.34 ns |  4,421.394 ns |  5,749.062 ns |  0.94 | Same            |    0.19 |  64.3382 | 22.0588 |  301112 B |        1.00 |
|                       |              |              |               |               |       |                 |         |          |         |           |             |
| **ScalarControl**         | **Tall**         | **87,468.23 ns** |    **797.670 ns** |  **1,037.196 ns** |  **1.00** | **Baseline**        |    **0.02** | **118.0556** | **57.2917** |  **602168 B** |        **1.00** |
| RowCachingOnly        | Tall         | 68,647.25 ns |    326.341 ns |    424.336 ns |  0.78 | Faster          |    0.01 | 118.0556 | 54.8611 |  602168 B |        1.00 |
| RowCachingAndBulkCopy | Tall         | 34,572.26 ns |    217.127 ns |    282.326 ns |  0.40 | Faster          |    0.01 | 119.4444 | 53.1250 |  602168 B |        1.00 |
| Production            | Tall         | 88,345.09 ns |    589.676 ns |    766.745 ns |  1.01 | Same            |    0.01 | 118.0556 | 54.6875 |  602168 B |        1.00 |
|                       |              |              |               |               |       |                 |         |          |         |           |             |
| **ScalarControl**         | **Wide**         | **86,484.35 ns** |    **677.112 ns** |    **880.436 ns** |  **1.00** | **Baseline**        |    **0.01** | **114.5833** | **71.1806** |  **596024 B** |        **1.00** |
| RowCachingOnly        | Wide         | 67,301.26 ns |    422.773 ns |    549.724 ns |  0.78 | Faster          |    0.01 | 114.9194 | 72.5806 |  596024 B |        1.00 |
| RowCachingAndBulkCopy | Wide         | 31,358.54 ns |    341.005 ns |    443.403 ns |  0.36 | Faster          |    0.01 | 115.2146 | 75.1263 |  596024 B |        1.00 |
| Production            | Wide         | 85,860.99 ns |    520.744 ns |    677.114 ns |  0.99 | Same            |    0.01 | 114.7260 | 71.0616 |  596024 B |        1.00 |
|                       |              |              |               |               |       |                 |         |          |         |           |             |
| **ScalarControl**         | **SingleRow**    |  **8,605.50 ns** |    **899.130 ns** |  **1,169.123 ns** |  **1.04** | **Baseline**        |    **0.36** |   **7.7347** |       **-** |   **32856 B** |        **1.00** |
| RowCachingOnly        | SingleRow    |  7,466.86 ns |     69.677 ns |     90.600 ns |  0.91 | Faster          |    0.29 |   7.7842 |       - |   32856 B |        1.00 |
| RowCachingAndBulkCopy | SingleRow    |    970.85 ns |      4.114 ns |      5.349 ns |  0.12 | Faster          |    0.04 |   7.8089 |       - |   32856 B |        1.00 |
| Production            | SingleRow    |  8,587.46 ns |    910.919 ns |  1,184.452 ns |  1.04 | Same            |    0.36 |   7.7347 |       - |   32856 B |        1.00 |
|                       |              |              |               |               |       |                 |         |          |         |           |             |
| **ScalarControl**         | **SingleColumn** | **36,049.87 ns** |    **870.715 ns** |  **1,132.175 ns** |  **1.00** | **Baseline**        |    **0.04** |  **37.8571** | **12.5000** |  **163896 B** |        **1.00** |
| RowCachingOnly        | SingleColumn | 35,745.61 ns |  1,586.435 ns |  2,062.814 ns |  0.99 | Same            |    0.06 |  37.8571 | 12.5000 |  163896 B |        1.00 |
| RowCachingAndBulkCopy | SingleColumn | 43,412.38 ns |  2,109.190 ns |  2,742.543 ns |  1.21 | Slower          |    0.08 |  37.8521 | 12.3239 |  163896 B |        1.00 |
| Production            | SingleColumn | 35,034.86 ns |  2,104.044 ns |  2,735.852 ns |  0.97 | Same            |    0.08 |  37.8472 | 12.5000 |  163896 B |        1.00 |
|                       |              |              |               |               |       |                 |         |          |         |           |             |
| **ScalarControl**         | **Empty**        |    **423.00 ns** |      **3.146 ns** |      **4.090 ns** |  **1.00** | **Baseline**        |    **0.01** |   **0.2552** |       **-** |    **1080 B** |        **1.00** |
| RowCachingOnly        | Empty        |    292.69 ns |      1.544 ns |      2.008 ns |  0.69 | Faster          |    0.01 |   0.2566 |       - |    1080 B |        1.00 |
| RowCachingAndBulkCopy | Empty        |    293.37 ns |      1.538 ns |      1.999 ns |  0.69 | Faster          |    0.01 |   0.2555 |       - |    1080 B |        1.00 |
| Production            | Empty        |    429.06 ns |      3.417 ns |      4.443 ns |  1.01 | Same            |    0.01 |   0.2574 |       - |    1080 B |        1.00 |
