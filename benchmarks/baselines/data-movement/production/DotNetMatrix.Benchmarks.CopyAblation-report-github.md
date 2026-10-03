```

BenchmarkDotNet v0.15.8, Linux Ubuntu 24.04.5 LTS (Noble Numbat)
Cortex-X925, Cortex-A725, 10 physical cores
.NET SDK 10.0.401
  [Host]     : .NET 10.0.12 (10.0.12, 10.0.1226.42308), Arm64 RyuJIT armv8.0-a
  Job-RRWQJM : .NET 10.0.12 (10.0.12, 10.0.1226.42308), Arm64 RyuJIT armv8.0-a

OutlierMode=DontRemove  IterationCount=12  IterationTime=100ms
LaunchCount=2  WarmupCount=4

```
| Method          | Shape        | Mean         | Error        | StdDev       | Ratio | MannWhitney(5%) | RatioSD | Gen0     | Gen1    | Allocated | Alloc Ratio |
|---------------- |------------- |-------------:|-------------:|-------------:|------:|---------------- |--------:|---------:|--------:|----------:|------------:|
| **ScalarControl**   | **Tiny**         |     **41.21 ns** |     **3.480 ns** |     **4.525 ns** |  **1.01** | **Baseline**        |    **0.13** |   **0.0228** |       **-** |      **96 B** |        **1.00** |
| GuardedBulkCopy | Tiny         |     42.15 ns |     3.469 ns |     4.511 ns |  1.03 | Same            |    0.13 |   0.0229 |       - |      96 B |        1.00 |
| WideRowsOnly    | Tiny         |     40.95 ns |     2.105 ns |     2.737 ns |  1.00 | Same            |    0.10 |   0.0228 |       - |      96 B |        1.00 |
| Production      | Tiny         |     44.16 ns |     0.148 ns |     0.192 ns |  1.08 | Slower          |    0.08 |   0.0226 |       - |      96 B |        1.00 |
|                 |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**   | **Small**        |  **1,623.16 ns** |    **23.467 ns** |    **30.513 ns** |  **1.00** | **Baseline**        |    **0.03** |   **0.6148** |       **-** |    **2616 B** |        **1.00** |
| GuardedBulkCopy | Small        |    379.90 ns |    18.280 ns |    23.769 ns |  0.23 | Faster          |    0.01 |   0.6233 |       - |    2616 B |        1.00 |
| WideRowsOnly    | Small        |    374.79 ns |    10.467 ns |    13.610 ns |  0.23 | Faster          |    0.01 |   0.6252 |       - |    2616 B |        1.00 |
| Production      | Small        |    380.04 ns |     8.326 ns |    10.826 ns |  0.23 | Faster          |    0.01 |   0.6221 |       - |    2616 B |        1.00 |
|                 |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**   | **Medium**       |  **9,823.21 ns** |    **50.506 ns** |    **65.672 ns** |  **1.00** | **Baseline**        |    **0.01** |   **8.2547** |       **-** |   **34872 B** |        **1.00** |
| GuardedBulkCopy | Medium       |  2,394.72 ns |    21.554 ns |    28.026 ns |  0.24 | Faster          |    0.00 |   8.3365 |       - |   34872 B |        1.00 |
| WideRowsOnly    | Medium       |  2,352.06 ns |   244.962 ns |   318.519 ns |  0.24 | Faster          |    0.03 |   8.3365 |       - |   34872 B |        1.00 |
| Production      | Medium       |  2,377.62 ns |    14.722 ns |    19.143 ns |  0.24 | Faster          |    0.00 |   8.3159 |       - |   34872 B |        1.00 |
|                 |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**   | **Large**        | **45,453.34 ns** |   **511.172 ns** |   **664.668 ns** |  **1.00** | **Baseline**        |    **0.02** |  **64.3116** | **22.1920** |  **301112 B** |        **1.00** |
| GuardedBulkCopy | Large        | 16,133.40 ns |   422.751 ns |   549.696 ns |  0.36 | Faster          |    0.01 |  64.6104 | 23.5390 |  301112 B |        1.00 |
| WideRowsOnly    | Large        | 15,936.94 ns |   363.194 ns |   472.255 ns |  0.35 | Faster          |    0.01 |  64.2538 | 22.6403 |  301112 B |        1.00 |
| Production      | Large        | 16,294.92 ns |   153.120 ns |   199.099 ns |  0.36 | Faster          |    0.01 |  64.1319 | 22.5196 |  301112 B |        1.00 |
|                 |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**   | **Tall**         | **87,395.29 ns** | **2,935.218 ns** | **3,816.613 ns** |  **1.00** | **Baseline**        |    **0.06** | **118.0556** | **55.5556** |  **602168 B** |        **1.00** |
| GuardedBulkCopy | Tall         | 33,956.85 ns |   277.936 ns |   361.395 ns |  0.39 | Faster          |    0.01 | 119.2010 | 53.8015 |  602168 B |        1.00 |
| WideRowsOnly    | Tall         | 34,241.50 ns |   234.624 ns |   305.078 ns |  0.39 | Faster          |    0.01 | 119.2255 | 54.0082 |  602168 B |        1.00 |
| Production      | Tall         | 34,065.39 ns |   177.059 ns |   230.227 ns |  0.39 | Faster          |    0.01 | 119.1621 | 53.9148 |  602168 B |        1.00 |
|                 |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**   | **Wide**         | **86,314.74 ns** |   **469.089 ns** |   **609.948 ns** |  **1.00** | **Baseline**        |    **0.01** | **115.3169** | **75.7042** |  **596024 B** |        **1.00** |
| GuardedBulkCopy | Wide         | 31,308.96 ns |   229.978 ns |   299.036 ns |  0.36 | Faster          |    0.00 | 115.2638 | 74.7487 |  596024 B |        1.00 |
| WideRowsOnly    | Wide         | 31,452.50 ns |   182.745 ns |   237.620 ns |  0.36 | Faster          |    0.00 | 115.0641 | 74.6795 |  596024 B |        1.00 |
| Production      | Wide         | 31,890.02 ns |   214.025 ns |   278.293 ns |  0.37 | Faster          |    0.00 | 115.2638 | 74.7487 |  596024 B |        1.00 |
|                 |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**   | **SingleRow**    |  **8,962.33 ns** |    **49.949 ns** |    **64.948 ns** |  **1.00** | **Baseline**        |    **0.01** |   **7.8014** |       **-** |   **32856 B** |        **1.00** |
| GuardedBulkCopy | SingleRow    |    968.19 ns |     3.964 ns |     5.154 ns |  0.11 | Faster          |    0.00 |   7.8089 |       - |   32856 B |        1.00 |
| WideRowsOnly    | SingleRow    |    977.67 ns |     3.796 ns |     4.935 ns |  0.11 | Faster          |    0.00 |   7.8064 |       - |   32856 B |        1.00 |
| Production      | SingleRow    |    973.48 ns |     3.809 ns |     4.953 ns |  0.11 | Faster          |    0.00 |   7.8064 |       - |   32856 B |        1.00 |
|                 |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**   | **SingleColumn** | **34,590.43 ns** |   **170.901 ns** |   **222.219 ns** |  **1.00** | **Baseline**        |    **0.01** |  **37.9834** | **12.4309** |  **163896 B** |        **1.00** |
| GuardedBulkCopy | SingleColumn | 36,030.85 ns |   427.459 ns |   555.818 ns |  1.04 | Same            |    0.02 |  37.7825 | 12.3588 |  163896 B |        1.00 |
| WideRowsOnly    | SingleColumn | 35,476.72 ns |   252.362 ns |   328.142 ns |  1.03 | Same            |    0.01 |  37.8571 | 12.5000 |  163896 B |        1.00 |
| Production      | SingleColumn | 36,052.22 ns |   165.009 ns |   214.558 ns |  1.04 | Same            |    0.01 |  37.7155 | 12.2126 |  163896 B |        1.00 |
|                 |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**   | **Empty**        |    **422.00 ns** |     **3.558 ns** |     **4.626 ns** |  **1.00** | **Baseline**        |    **0.02** |   **0.2540** |       **-** |    **1080 B** |        **1.00** |
| GuardedBulkCopy | Empty        |    421.43 ns |     3.161 ns |     4.110 ns |  1.00 | Same            |    0.01 |   0.2564 |       - |    1080 B |        1.00 |
| WideRowsOnly    | Empty        |    419.91 ns |     2.792 ns |     3.630 ns |  1.00 | Same            |    0.01 |   0.2567 |       - |    1080 B |        1.00 |
| Production      | Empty        |    424.49 ns |     2.010 ns |     2.613 ns |  1.01 | Same            |    0.01 |   0.2580 |       - |    1080 B |        1.00 |
|                 |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**   | **TwoColumns**   | **40,785.37 ns** | **1,017.600 ns** | **1,323.168 ns** |  **1.00** | **Baseline**        |    **0.05** |  **45.6414** | **11.9243** |  **196664 B** |        **1.00** |
| GuardedBulkCopy | TwoColumns   | 46,701.79 ns |   638.883 ns |   830.728 ns |  1.15 | Slower          |    0.04 |  45.4545 | 11.8371 |  196664 B |        1.00 |
| WideRowsOnly    | TwoColumns   | 42,086.30 ns |   563.447 ns |   732.640 ns |  1.03 | Same            |    0.04 |  45.8333 | 12.0833 |  196664 B |        1.00 |
| Production      | TwoColumns   | 40,767.81 ns |   434.776 ns |   565.332 ns |  1.00 | Same            |    0.03 |  45.6414 | 11.9243 |  196664 B |        1.00 |
|                 |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**   | **FourColumns**  | **49,996.53 ns** |   **463.120 ns** |   **602.187 ns** |  **1.00** | **Baseline**        |    **0.02** |  **60.3814** | **20.1271** |  **262200 B** |        **1.00** |
| GuardedBulkCopy | FourColumns  | 48,771.61 ns |   427.309 ns |   555.622 ns |  0.98 | Same            |    0.02 |  60.4167 | 19.7917 |  262200 B |        1.00 |
| WideRowsOnly    | FourColumns  | 49,951.60 ns |   188.542 ns |   245.158 ns |  1.00 | Same            |    0.01 |  60.5315 | 20.1772 |  262200 B |        1.00 |
| Production      | FourColumns  | 49,462.18 ns |   644.808 ns |   838.433 ns |  0.99 | Same            |    0.02 |  60.5620 | 19.8643 |  262200 B |        1.00 |
|                 |              |              |              |              |       |                 |         |          |         |           |             |
| **ScalarControl**   | **EightColumns** | **67,930.42 ns** |   **352.294 ns** |   **458.082 ns** |  **1.00** | **Baseline**        |    **0.01** |  **76.0870** | **42.1196** |  **393272 B** |        **1.00** |
| GuardedBulkCopy | EightColumns | 56,265.64 ns |   314.540 ns |   408.991 ns |  0.83 | Faster          |    0.01 |  75.8929 | 43.5268 |  393272 B |        1.00 |
| WideRowsOnly    | EightColumns | 67,732.28 ns |   366.319 ns |   476.318 ns |  1.00 | Same            |    0.01 |  76.0870 | 42.1196 |  393272 B |        1.00 |
| Production      | EightColumns | 68,111.26 ns |   609.879 ns |   793.015 ns |  1.00 | Same            |    0.01 |  75.6944 | 41.6667 |  393272 B |        1.00 |
