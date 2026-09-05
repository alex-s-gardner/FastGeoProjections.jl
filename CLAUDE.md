# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package is

Native-Julia coordinate transformations between EPSG-coded CRSs, faster than Proj for the
projections it implements and delegating to `Proj.jl` for everything else. Not a Proj replacement:
`fast_epsg_codes` in `src/epsg.jl` is the whole native set — 3031, 3413, 3857, 4326, 4978, 4979, plus
UTM 326XX/327XX.

Accuracy is defined against Proj. Every native projection is asserted against it in the test suite,
to 5.5e-9 m for UTM, ~2e-9 m for the geocentric conversions, and 2.4e-8 m for the polar
stereographic inverse — the least accurate of them, and the one place where the limit is a
truncated series rather than the last bits of a transcendental. See
[`PolarStereographicToLonLat`](@ref) on why it stops where it does.

## Commands

```bash
julia --project=. -e 'import Pkg; Pkg.test()'                    # full suite
julia --project=. --threads=8 -e 'import Pkg; Pkg.test()'        # threading paths need >1 thread

# single test file
julia --project=. -e 'using TestEnv; TestEnv.activate(); include("test/transformations.jl")'
```

CI runs 1.10 and nightly on ubuntu with `JULIA_NUM_THREADS=4`. No docs project, no formatter.

The committed `benchmark/` project pulls GLMakie, which will not precompile headless — for
benchmarking, make a scratch environment with `Pkg.develop(path=".")` plus BenchmarkTools rather
than using `--project=benchmark`. `benchmark/benchmark.jl` and the `benchmark.jpg` in the README
predate the point-operator restructure and no longer run.

## Architecture

A transformation is a **point operator**: a callable struct holding everything derived from the
projection parameters, mapping one point to one point. Threading and vectorization live in
`apply.jl`, not in the projections, so there is one implementation of each projection rather than
one per calling convention.

- `transformations.jl` — `abstract type GeoTransformation`, the traits every operator answers
  (`islanesafe`, `preservesz`, `adapt_eltype`), and `ComposedGeoTransformation`, a flat tuple that
  inlines into one pass.
- `apply.jl` — `transform`/`transform!` over collections. Picks a SIMD lane loop or a scalar loop
  per chunk, threads over chunks, and decides how a third coordinate is handled.
- `epsg.jl` — the EPSG registry. `pipeline()` composes `project_from_4326` with `project_to_4326`,
  since every projection is expressed relative to EPSG:4326 in (lon, lat) order.
- `coord.jl` — `Transformation`, the public entry point, wrapping an operator plus the EPSG pair
  and flags. Implements the CoordinateTransformations API.
- `proj.jl` — `ProjTransformation`, the fallback. Holds one cloned PROJ context per thread in a
  `Channel` pool, checked out via `borrow` once per chunk (a context is single-thread-only).
- `kernels.jl` — the `Math` submodule: `sin`, `cos`, `pow`, `cbrt`, … each taking a `MathKernel`
  first argument, dispatching to SLEEFPirates' `_fast` routines (`FastKernel`, the default), its
  fully-reduced ones (`SLEEFKernel`), or Base's libm (`BaseKernel`). These are plain Julia, so they
  inline into an operator and evaluate on `Vec` lanes as well as scalars — that is what lets a
  point operator reach array-kernel throughput.
- `projections/` — one file per projection family, each a `LonLatTo…`/`…ToLonLat` pair.

### Adding a projection

Add branches to both `project_to_4326` and `project_from_4326` in `epsg.jl` (marked by the
`## ⬇ ADD FAST PROJECTIONS HERE ⬇ ##` comment), append the code to `fast_epsg_codes`, and classify
it in `isgeographic`. The test suite walks `fast_epsg_codes` and checks each against Proj at both
axis orders, so an unclassified code fails rather than silently returning x and y reversed.

## Conventions specific to this code

- **EPSG codes hold a tuple**, not a scalar — read them as `first(epsg.val)`.
- **`always_xy` differs by layer.** `false` on `Transformation` (authority order, so EPSG:4326 is
  lat, lon); the operators underneath are *always* xy. Axis order is handled by composing `SwapXY`
  onto the geographic end of a pipeline, not by a per-point branch.
- **`islanesafe` and `preservesz` are independent traits.** Lane-safe means evaluable on `Vec`
  lanes; `preservesz` means a third coordinate passes through untouched. A map projection is both;
  the geocentric conversions are lane-safe but transform the height. Code that conflates them takes
  the wrong path silently.
- **A height is transformed, not carried, where the transformation changes one.** For a datum shift
  or a geocentric conversion, the two-coordinate call is the `h = 0` point — which moves x and y by
  metres — so `LonLatToGeocentric` and `GeocentricToLonLat` reject it rather than assume sea level.
- **Fusion.** Composing `GeocentricToLonLat` with a projection would form `atan(z, d)` and
  `atan(y, x)` and then immediately take their sines and cosines. `Direction` (`projections/direction.jl`)
  carries the unnormalized direction instead, a projection opts in with `project_direction`, and
  `pipeline` substitutes `FusedFromGeocentric`. This is the one place a pipeline is not literally
  `to ∘ from`.
- **`tranmerc.jl` is written for the compiler, not the reader.** Branches are `ifelse` and boolean
  arithmetic so the body vectorizes; the Clenshaw recurrences are spelled out rather than looped.
  `_angdiff`'s two-term compensated arithmetic is load-bearing at the ±180° branch cut and cannot
  be replaced by an angle-subtraction identity — which is why transverse Mercator fuses in latitude
  only.
- `@inbounds` in `apply.jl`'s lane and scalar loops is deliberate and the indices are derived from
  `eachindex`/lane masks. Don't add more without a benchmark.

## Performance notes

Measured on aarch64 (M-series), 1 thread, 100k points, out of place:

| pipeline | ns/point |
|---|---|
| 4326→3857 web Mercator | ~12 |
| 3857→4326 web Mercator, inverse | ~4 |
| 4326→3413 polar stereographic | ~10 |
| 3413→4326 polar stereographic, inverse | ~19 |
| 4326→32619 UTM | ~48 |
| 4978→3413 fused geocentric | ~34 |

The inverse costs twice the forward because the conformal→geodetic series is five
`Math.sin` calls against the forward's one `Math.conformal_ratio`.

Two things to know before optimizing:

- **`pick_vector_width(Float64)` is 2 on aarch64**, not 8. A `VecUnroll{U,2}` in a `@code_warntype`
  dump is an unroll factor, not a lane count. Per-lane speedups are 1.7–2.3×, not 8×.
- **Run-to-run noise on the array paths is about ±8%.** A single measurement showing a 5%
  improvement has measured nothing. Repeat a baseline three times before believing any variant.

`transform!` is allocation-free in place; out-of-place allocation is the destination array alone.
The scalar operators are type-stable with zero JET optimization reports.

**A general `pow` does not vectorize, and it caps a projection's throughput.** `Math.pow` costs
about five times any other transcendental in a loop, and `code_llvm` on that loop contains no
`<2 x double>` at all where the same loop over `Math.conformal_ratio` contains 25.

The cause is composition, not the function. `SLEEFPirates.pow_fast` is one line —
`exp2(y * log2_fast(x))` — and both halves vectorize alone: 1.83 and 1.46 ns/point separately, 3.34
in two passes over an array, but 6.55 fused into one loop and 10.28 through `pow_fast`. LLVM will
not vectorize a loop body holding both inlined. Nothing upstream can fix that, because a correct
general `pow` has to handle negative bases, integer exponents and the infinities, and is therefore
always too large to vectorize.

`Math.conformal_ratio` wins by not being general: `e < 0.09` and `|u| ≤ e` mean no range reduction
and no special cases, so twelve fused multiply-adds reach one ulp. Prefer a short series over
`Math.pow` wherever an exponent is fixed and its argument range is bounded — and check the IR for
`<2 x double>` rather than trusting a timing, since the array-path noise band swallows anything
under about 8%.

### Known gaps

- **A height-transforming operator never reaches the SIMD lane loop.** `_transform_pts!` sends it
  to the scalar path, because the interleaved loop writes rows 1 and 2 and preserves the rest. The
  cost is small — geocentric on the scalar path is ~12 ns/point, faster than any 2-D operator on
  the lane path — but a 3-coordinate lane loop would recover roughly 2× on that operator.
- **`_lane_range_aos!` is not worth tuning.** Running `Identity` through it measures 0.68–0.81
  ns/point, so the loop is ~4% of a projection's runtime and the strided AoS load costs 0.44
  ns/point over the contiguous SoA one. `UNROLL` (2/4/8/16) and `CHUNK` are all inside the noise
  band. Time is in the projection math.
