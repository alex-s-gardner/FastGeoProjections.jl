[![Build Status](https://github.com/alex-s-gardner/FastGeoProjections.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/alex-s-gardner/FastGeoProjections.jl/actions/workflows/CI.yml?query=branch%3Amain)

!! UNDER DEVELOPMENT BY OVERCOMMITTED AND UNDER PAID DEVELOPERS !!

**FastGeoProjections** is intended to provide highly optimized native Julia geospatial coordinate transformations from one coordinate reference system (CRS) to another as defined by EPSG codes. It is not intended to to replace or to be as comprehensive as [Proj](https://github.com/JuliaGeo/Proj.jl). The package will natively support only the most common geospatial transformations and relies on **Proj.jl** for all others.

Benchmark of currently implemented EPSGs

**EPSG:4326 to EPSG:3413 [n = 1000]**





*Proj: single-thread*


BenchmarkTools.TrialEstimate: 


  time:             186.500 μs


  gctime:           0.000 ns (0.00%)


  memory:           64.44 KiB


  allocs:           36


  **max error:			0.0**





*Proj: multi-thread*


BenchmarkTools.TrialEstimate: 


  time:             5.358 ms


  gctime:           0.000 ns (0.00%)


  memory:           80.05 KiB


  allocs:           1658


  **max error:			0.0**





*FastGeoProjections: single-thread*


BenchmarkTools.TrialEstimate: 


  time:             41.333 μs


  gctime:           0.000 ns (0.00%)


  memory:           190.91 KiB


  allocs:           37


  **max error:			0.0001606**





*FastGeoProjections: multi-thread - Float64*


BenchmarkTools.TrialEstimate: 


  time:             19.042 μs


  gctime:           0.000 ns (0.00%)


  memory:           182.97 KiB


  allocs:           36


  **max error:			0.0001606**





*FastGeoProjections: multi-thread - Float32*


BenchmarkTools.TrialEstimate: 


  time:             7.667 μs


  gctime:           0.000 ns (0.00%)


  memory:           93.84 KiB


  allocs:           36


  **max error:			1.046**





*FastGeoProjections: M2 GPU - Float32 [including transfer time]*


BenchmarkTools.TrialEstimate: 


  time:             2.466 ms


  gctime:           0.000 ns (0.00%)


  memory:           146.69 KiB


  allocs:           4346


  **max error:			1.105**








**EPSG:3031 to EPSG:4326 [n = 1000]**





*Proj: single-thread*


BenchmarkTools.TrialEstimate: 


  time:             397.459 μs


  gctime:           0.000 ns (0.00%)


  memory:           64.44 KiB


  allocs:           36


  **max error:			0.0**





*Proj: multi-thread*


BenchmarkTools.TrialEstimate: 


  time:             5.245 ms


  gctime:           0.000 ns (0.00%)


  memory:           80.14 KiB


  allocs:           1661


  **max error:			0.0**





*FastGeoProjections: single-thread*


BenchmarkTools.TrialEstimate: 


  time:             44.917 μs


  gctime:           0.000 ns (0.00%)


  memory:           246.47 KiB


  allocs:           44


  **max error:			0.0**





*FastGeoProjections: multi-thread - Float64*


BenchmarkTools.TrialEstimate: 


  time:             20.250 μs


  gctime:           0.000 ns (0.00%)


  memory:           246.47 KiB


  allocs:           44


  **max error:			0.0**





*FastGeoProjections: multi-thread - Float32*


BenchmarkTools.TrialEstimate: 


  time:             8.500 μs


  gctime:           0.000 ns (0.00%)


  memory:           126.34 KiB


  allocs:           44


  **max error:			3.36e-5**





*FastGeoProjections: M2 GPU - Float32 [including transfer time]*


BenchmarkTools.TrialEstimate: 


  time:             3.191 ms


  gctime:           0.000 ns (0.00%)


  memory:           182.63 KiB


  allocs:           5466


  **max error:			3.654e-5**








**EPSG:4326 to EPSG:3413 [n = 1000000]**





*Proj: single-thread*


BenchmarkTools.TrialEstimate: 


  time:             107.239 ms


  gctime:           0.000 ns (0.00%)


  memory:           61.04 MiB


  allocs:           43


  **max error:			0.0**





*Proj: multi-thread*


BenchmarkTools.TrialEstimate: 


  time:             22.278 ms


  gctime:           0.000 ns (0.00%)


  memory:           76.30 MiB


  allocs:           1999666


  **max error:			0.0**





*FastGeoProjections: single-thread*


BenchmarkTools.TrialEstimate: 


  time:             63.616 ms


  gctime:           4.349 ms (6.84%)


  memory:           183.11 MiB


  allocs:           61


  **max error:			0.0001606**





*FastGeoProjections: multi-thread - Float64*


BenchmarkTools.TrialEstimate: 


  time:             24.591 ms


  gctime:           1.981 ms (8.06%)


  memory:           175.48 MiB


  allocs:           59


  **max error:			0.0001606**





*FastGeoProjections: multi-thread - Float32*


BenchmarkTools.TrialEstimate: 


  time:             9.596 ms


  gctime:           1.573 ms (16.39%)


  memory:           87.74 MiB


  allocs:           59


  **max error:			1.538**





*FastGeoProjections: M2 GPU - Float32 [including transfer time]*


BenchmarkTools.TrialEstimate: 


  time:             6.042 ms


  gctime:           0.000 ns (0.00%)


  memory:           7.77 MiB


  allocs:           4372


  **max error:			1.736**








**EPSG:3031 to EPSG:4326 [n = 1000000]**





*Proj: single-thread*


BenchmarkTools.TrialEstimate: 


  time:             349.651 ms


  gctime:           0.000 ns (0.00%)


  memory:           61.04 MiB


  allocs:           43


  **max error:			0.0**





*Proj: multi-thread*


BenchmarkTools.TrialEstimate: 


  time:             55.162 ms


  gctime:           0.000 ns (0.00%)


  memory:           76.30 MiB


  allocs:           1999666


  **max error:			0.0**





*FastGeoProjections: single-thread*


BenchmarkTools.TrialEstimate: 


  time:             73.833 ms


  gctime:           2.276 ms (3.08%)


  memory:           236.51 MiB


  allocs:           75


  **max error:			0.0**





*FastGeoProjections: multi-thread - Float64*


BenchmarkTools.TrialEstimate: 


  time:             26.440 ms


  gctime:           2.374 ms (8.98%)


  memory:           236.51 MiB


  allocs:           75


  **max error:			0.0**





*FastGeoProjections: multi-thread - Float32*


BenchmarkTools.TrialEstimate: 


  time:             10.800 ms


  gctime:           1.643 ms (15.21%)


  memory:           118.26 MiB


  allocs:           75


  **max error:			3.856e-5**





*FastGeoProjections: M2 GPU - Float32 [including transfer time]*


BenchmarkTools.TrialEstimate: 


  time:             6.938 ms


  gctime:           0.000 ns (0.00%)


  memory:           7.80 MiB


  allocs:           5499


  **max error:			4.1e-5**








**EPSG:4326 to EPSG:3413 [n = 10000000]**





*Proj: single-thread*


BenchmarkTools.TrialEstimate: 


  time:             1.111 s


  gctime:           12.125 ms (1.09%)


  memory:           610.35 MiB


  allocs:           43


  **max error:			0.0**





*Proj: multi-thread*


BenchmarkTools.TrialEstimate: 


  time:             187.946 ms


  gctime:           0.000 ns (0.00%)


  memory:           762.94 MiB


  allocs:           19999666


  **max error:			0.0**





*FastGeoProjections: single-thread*


BenchmarkTools.TrialEstimate: 


  time:             660.070 ms


  gctime:           20.549 ms (3.11%)


  memory:           1.79 GiB


  allocs:           61


  **max error:			0.0001606**





*FastGeoProjections: multi-thread - Float64*


BenchmarkTools.TrialEstimate: 


  time:             282.487 ms


  gctime:           17.656 ms (6.25%)


  memory:           1.71 GiB


  allocs:           59


  **max error:			0.0001606**





*FastGeoProjections: multi-thread - Float32*


BenchmarkTools.TrialEstimate: 


  time:             122.212 ms


  gctime:           9.555 ms (7.82%)


  memory:           877.38 MiB


  allocs:           59


  **max error:			1.602**





*FastGeoProjections: M2 GPU - Float32 [including transfer time]*


BenchmarkTools.TrialEstimate: 


  time:             62.645 ms


  gctime:           0.000 ns (0.00%)


  memory:           76.43 MiB


  allocs:           4372


  **max error:			1.776**








**EPSG:3031 to EPSG:4326 [n = 10000000]**





*Proj: single-thread*


BenchmarkTools.TrialEstimate: 


  time:             3.588 s


  gctime:           4.478 ms (0.12%)


  memory:           610.35 MiB


  allocs:           43


  **max error:			0.0**





*Proj: multi-thread*


BenchmarkTools.TrialEstimate: 


  time:             503.967 ms


  gctime:           0.000 ns (0.00%)


  memory:           762.94 MiB


  allocs:           19999667


  **max error:			0.0**





*FastGeoProjections: single-thread*


BenchmarkTools.TrialEstimate: 


  time:             800.905 ms


  gctime:           33.112 ms (4.13%)


  memory:           2.31 GiB


  allocs:           75


  **max error:			0.0**





*FastGeoProjections: multi-thread - Float64*


BenchmarkTools.TrialEstimate: 


  time:             329.724 ms


  gctime:           18.215 ms (5.52%)


  memory:           2.31 GiB


  allocs:           75


  **max error:			0.0**





*FastGeoProjections: multi-thread - Float32*


BenchmarkTools.TrialEstimate: 


  time:             146.504 ms


  gctime:           15.176 ms (10.36%)


  memory:           1.15 GiB


  allocs:           75


  **max error:			4.04e-5**





*FastGeoProjections: M2 GPU - Float32 [including transfer time]*


BenchmarkTools.TrialEstimate: 


  time:             122.944 ms


  gctime:           0.000 ns (0.00%)


  memory:           76.47 MiB


  allocs:           5499


  **max error:			4.184e-5**








