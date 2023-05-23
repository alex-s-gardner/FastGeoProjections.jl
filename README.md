[![Build Status](https://github.com/alex-s-gardner/FastGeoProjections.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/alex-s-gardner/FastGeoProjections.jl/actions/workflows/CI.yml?query=branch%3Amain)

!! UNDER DEVELOPMENT BY OVERCOMMITTED AND UNDER PAID DEVELOPERS !!

**FastGeoProjections** is intended to provide highly optimized native Julia geospatial coordinate transformations from one coordinate reference system (CRS) to another as defined by EPSG codes. It is not intended to to replace or to be as comprehensive as [Proj](https://github.com/JuliaGeo/Proj.jl). The package will natively support only the most common geospatial transformations and relies on **Proj.jl** for all others.

Benchmark of currently implemented EPSGs

**EPSG:4326 to EPSG:3413 [n = 1000]**

	

*Proj: single-thread*

	BenchmarkTools.TrialEstimate: 

	  time:             180.834 μs

	  gctime:           0.000 ns (0.00%)

	  memory:           64.44 KiB

	  allocs:           36

	  MAXIMUM ERROR:		0.0

	

*Proj: multi-thread*

	BenchmarkTools.TrialEstimate: 

	  time:             5.388 ms

	  gctime:           0.000 ns (0.00%)

	  memory:           80.05 KiB

	  allocs:           1658

	  MAXIMUM ERROR:		0.0

	

*FastGeoProjections: single-thread*

	BenchmarkTools.TrialEstimate: 

	  time:             40.958 μs

	  gctime:           0.000 ns (0.00%)

	  memory:           190.91 KiB

	  allocs:           37

	  MAXIMUM ERROR:		0.0001606

	

*FastGeoProjections: multi-thread - Float64*

	BenchmarkTools.TrialEstimate: 

	  time:             17.541 μs

	  gctime:           0.000 ns (0.00%)

	  memory:           182.97 KiB

	  allocs:           36

	  MAXIMUM ERROR:		0.0001606

	

*FastGeoProjections: multi-thread - Float32*

	BenchmarkTools.TrialEstimate: 

	  time:             7.792 μs

	  gctime:           0.000 ns (0.00%)

	  memory:           93.84 KiB

	  allocs:           36

	  MAXIMUM ERROR:		1.171

	

*FastGeoProjections: M2 GPU - Float32 [including transfer time]*

	BenchmarkTools.TrialEstimate: 

	  time:             2.431 ms

	  gctime:           0.000 ns (0.00%)

	  memory:           146.69 KiB

	  allocs:           4346

	  MAXIMUM ERROR:		1.228

	

	

**EPSG:3031 to EPSG:4326 [n = 1000]**

	

*Proj: single-thread*

	BenchmarkTools.TrialEstimate: 

	  time:             396.833 μs

	  gctime:           0.000 ns (0.00%)

	  memory:           64.44 KiB

	  allocs:           36

	  MAXIMUM ERROR:		0.0

	

*Proj: multi-thread*

	BenchmarkTools.TrialEstimate: 

	  time:             5.213 ms

	  gctime:           0.000 ns (0.00%)

	  memory:           80.08 KiB

	  allocs:           1659

	  MAXIMUM ERROR:		0.0

	

*FastGeoProjections: single-thread*

	BenchmarkTools.TrialEstimate: 

	  time:             44.333 μs

	  gctime:           0.000 ns (0.00%)

	  memory:           246.47 KiB

	  allocs:           44

	  MAXIMUM ERROR:		0.0

	

*FastGeoProjections: multi-thread - Float64*

	BenchmarkTools.TrialEstimate: 

	  time:             18.750 μs

	  gctime:           0.000 ns (0.00%)

	  memory:           246.47 KiB

	  allocs:           44

	  MAXIMUM ERROR:		0.0

	

*FastGeoProjections: multi-thread - Float32*

	BenchmarkTools.TrialEstimate: 

	  time:             8.583 μs

	  gctime:           0.000 ns (0.00%)

	  memory:           126.34 KiB

	  allocs:           44

	  MAXIMUM ERROR:		3.48e-5

	

*FastGeoProjections: M2 GPU - Float32 [including transfer time]*

	BenchmarkTools.TrialEstimate: 

	  time:             3.093 ms

	  gctime:           0.000 ns (0.00%)

	  memory:           182.63 KiB

	  allocs:           5466

	  MAXIMUM ERROR:		3.48e-5

	

	

**EPSG:4326 to EPSG:3413 [n = 1000000]**

	

*Proj: single-thread*

	BenchmarkTools.TrialEstimate: 

	  time:             108.872 ms

	  gctime:           0.000 ns (0.00%)

	  memory:           61.04 MiB

	  allocs:           43

	  MAXIMUM ERROR:		0.0

	

*Proj: multi-thread*

	BenchmarkTools.TrialEstimate: 

	  time:             22.950 ms

	  gctime:           0.000 ns (0.00%)

	  memory:           76.30 MiB

	  allocs:           1999666

	  MAXIMUM ERROR:		0.0

	

*FastGeoProjections: single-thread*

	BenchmarkTools.TrialEstimate: 

	  time:             59.284 ms

	  gctime:           1.810 ms (3.05%)

	  memory:           183.11 MiB

	  allocs:           61

	  MAXIMUM ERROR:		0.0001606

	

*FastGeoProjections: multi-thread - Float64*

	BenchmarkTools.TrialEstimate: 

	  time:             20.887 ms

	  gctime:           2.468 ms (11.81%)

	  memory:           175.48 MiB

	  allocs:           59

	  MAXIMUM ERROR:		0.0001606

	

*FastGeoProjections: multi-thread - Float32*

	BenchmarkTools.TrialEstimate: 

	  time:             7.766 ms

	  gctime:           0.000 ns (0.00%)

	  memory:           87.74 MiB

	  allocs:           59

	  MAXIMUM ERROR:		1.587

	

*FastGeoProjections: M2 GPU - Float32 [including transfer time]*

	BenchmarkTools.TrialEstimate: 

	  time:             5.861 ms

	  gctime:           0.000 ns (0.00%)

	  memory:           7.77 MiB

	  allocs:           4372

	  MAXIMUM ERROR:		1.732

	

	

**EPSG:3031 to EPSG:4326 [n = 1000000]**

	

*Proj: single-thread*

	BenchmarkTools.TrialEstimate: 

	  time:             350.848 ms

	  gctime:           0.000 ns (0.00%)

	  memory:           61.04 MiB

	  allocs:           43

	  MAXIMUM ERROR:		0.0

	

*Proj: multi-thread*

	BenchmarkTools.TrialEstimate: 

	  time:             54.614 ms

	  gctime:           0.000 ns (0.00%)

	  memory:           76.30 MiB

	  allocs:           1999666

	  MAXIMUM ERROR:		0.0

	

*FastGeoProjections: single-thread*

	BenchmarkTools.TrialEstimate: 

	  time:             73.675 ms

	  gctime:           3.215 ms (4.36%)

	  memory:           236.51 MiB

	  allocs:           75

	  MAXIMUM ERROR:		0.0

	

*FastGeoProjections: multi-thread - Float64*

	BenchmarkTools.TrialEstimate: 

	  time:             27.349 ms

	  gctime:           2.848 ms (10.42%)

	  memory:           236.51 MiB

	  allocs:           75

	  MAXIMUM ERROR:		0.0

	

*FastGeoProjections: multi-thread - Float32*

	BenchmarkTools.TrialEstimate: 

	  time:             10.781 ms

	  gctime:           1.947 ms (18.06%)

	  memory:           118.26 MiB

	  allocs:           75

	  MAXIMUM ERROR:		3.844e-5

	

*FastGeoProjections: M2 GPU - Float32 [including transfer time]*

	BenchmarkTools.TrialEstimate: 

	  time:             16.509 ms

	  gctime:           0.000 ns (0.00%)

	  memory:           7.80 MiB

	  allocs:           5499

	  MAXIMUM ERROR:		4.13e-5

	

	

**EPSG:4326 to EPSG:3413 [n = 10000000]**

	

*Proj: single-thread*

	BenchmarkTools.TrialEstimate: 

	  time:             1.087 s

	  gctime:           4.375 ms (0.40%)

	  memory:           610.35 MiB

	  allocs:           43

	  MAXIMUM ERROR:		0.0

	

*Proj: multi-thread*

	BenchmarkTools.TrialEstimate: 

	  time:             186.603 ms

	  gctime:           0.000 ns (0.00%)

	  memory:           762.94 MiB

	  allocs:           19999666

	  MAXIMUM ERROR:		0.0

	

*FastGeoProjections: single-thread*

	BenchmarkTools.TrialEstimate: 

	  time:             667.795 ms

	  gctime:           22.395 ms (3.35%)

	  memory:           1.79 GiB

	  allocs:           61

	  MAXIMUM ERROR:		0.0001606

	

*FastGeoProjections: multi-thread - Float64*

	BenchmarkTools.TrialEstimate: 

	  time:             286.764 ms

	  gctime:           18.025 ms (6.29%)

	  memory:           1.71 GiB

	  allocs:           59

	  MAXIMUM ERROR:		0.0001606

	

*FastGeoProjections: multi-thread - Float32*

	BenchmarkTools.TrialEstimate: 

	  time:             116.330 ms

	  gctime:           9.169 ms (7.88%)

	  memory:           877.38 MiB

	  allocs:           59

	  MAXIMUM ERROR:		1.598

	

*FastGeoProjections: M2 GPU - Float32 [including transfer time]*

	BenchmarkTools.TrialEstimate: 

	  time:             59.444 ms

	  gctime:           0.000 ns (0.00%)

	  memory:           76.43 MiB

	  allocs:           4372

	  MAXIMUM ERROR:		1.905

	

	

**EPSG:3031 to EPSG:4326 [n = 10000000]**

	

*Proj: single-thread*

	BenchmarkTools.TrialEstimate: 

	  time:             3.585 s

	  gctime:           0.000 ns (0.00%)

	  memory:           610.35 MiB

	  allocs:           43

	  MAXIMUM ERROR:		0.0

	

*Proj: multi-thread*

	BenchmarkTools.TrialEstimate: 

	  time:             544.866 ms

	  gctime:           0.000 ns (0.00%)

	  memory:           762.94 MiB

	  allocs:           19999667

	  MAXIMUM ERROR:		0.0

	

*FastGeoProjections: single-thread*

	BenchmarkTools.TrialEstimate: 

	  time:             805.097 ms

	  gctime:           29.917 ms (3.72%)

	  memory:           2.31 GiB

	  allocs:           75

	  MAXIMUM ERROR:		0.0

	

*FastGeoProjections: multi-thread - Float64*

	BenchmarkTools.TrialEstimate: 

	  time:             340.120 ms

	  gctime:           26.992 ms (7.94%)

	  memory:           2.31 GiB

	  allocs:           75

	  MAXIMUM ERROR:		0.0

	

*FastGeoProjections: multi-thread - Float32*

	BenchmarkTools.TrialEstimate: 

	  time:             148.571 ms

	  gctime:           17.788 ms (11.97%)

	  memory:           1.15 GiB

	  allocs:           75

	  MAXIMUM ERROR:		3.97e-5

	

*FastGeoProjections: M2 GPU - Float32 [including transfer time]*

	BenchmarkTools.TrialEstimate: 

	  time:             79.990 ms

	  gctime:           0.000 ns (0.00%)

	  memory:           76.47 MiB

	  allocs:           5499

	  MAXIMUM ERROR:		4.137e-5

	

	

