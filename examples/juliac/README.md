# `fastgeoproj` — a standalone CSV reprojection tool

A worked example of compiling FastGeoProjections ahead of time with
[`juliac`](https://github.com/JuliaLang/julia/tree/master/contrib/juliac) into
an executable with no Julia startup and no compilation at run time.

```console
$ ./fastgeoproj --input sample.csv --from 4326 --to 32619 --always-xy
500000.0,4.982950400226553e6
380244.57633207337,4.900734380124455e6
557548.832534876,5.149876586070551e6
426642.267622362,4.761207746467653e6
319952.79600128037,4.735400315222241e6
```

Input is strictly two columns of `Float64`, comma separated, no header. Blank
lines are skipped; `\r\n`, surrounding spaces and a missing final newline are
all accepted. Anything else is an error naming the line.

## Build

Needs **Julia 1.12 or later** and a C compiler on `PATH` for the final link.

```console
$ julia --project=examples/juliac examples/juliac/build.jl
```

That produces `examples/juliac/fastgeoproj` (about 17 MiB) in about half a
minute. `build.jl --trim=safe` makes the build fail on any unresolved call
instead of warning; see *What does not survive trimming* below for the two that
are expected.

## Use

```
Usage: fastgeoproj --input FILE --from EPSG --to EPSG [--output FILE]

  -i, --input FILE    input CSV
  -o, --output FILE   output CSV (default: stdout)
  -f, --from EPSG     source CRS, as 4326 or EPSG:4326
  -t, --to EPSG       target CRS
      --always-xy     read and write a geographic CRS as (lon, lat) rather
                      than in the authority order (lat, lon)
```

Supported codes are the ones FastGeoProjections implements natively — 4326,
3031, 3413, and the 120 WGS 84 UTM zones (32601–32660, 32701–32760). Any pair
of them works, composed through EPSG:4326: `--from 32619 --to 32620` is one
pass over the data, not two.

Without `--always-xy`, EPSG:4326 columns are `(lat, lon)`, which is the axis
order the authority defines and what `cs2cs` expects. Projected CRSs are always
`(easting, northing)`.

## How it is written

Two things in `fastgeoproj.jl` are what make it compile under `--trim`, and
both generalize to any app built on this package.

**The projection pipeline has to be a type, not a value.**
`Transformation(EPSG(from), EPSG(to))` picks a projection from a run-time
integer, so the type of what it returns is only known at run time — which the
trimmer cannot follow. `with_source`/`with_target` do the same job by passing a
*concrete* operator forward to a continuation instead of returning one:

```julia
@inline function with_source(k::K, code::Int, always_xy::Bool)::Nothing where {K}
    if code == 4326
        always_xy ? k(Identity()) : k(SwapXY())
    elseif code == 3031
        k(PolarStereographicToLonLat{F}(; lat_ts = -71.0, lon_0 = 0.0))
    ...
```

Each branch specializes `k` on one operator type, so the transform is compiled
with the projection fully inlined and nothing dispatches at run time. Three
source operators and three target ones means a handful of copies of the loop in
the binary, all statically reachable.

**The points are already in the layout the SIMD path wants.** A CSV of pairs
parses naturally into a `Vector{NTuple{2,Float64}}`, which is a dense
interleaved buffer of coordinates — exactly what `transform!` reads and writes
on SIMD lanes. No repacking, no struct-of-arrays shuffle.

Whether a point type really has that layout is settled from the type rather
than from a vector's contents, and the answer is a `Val{T}` rather than a
`DataType` value, so the element type of the buffer is known statically. That
is a package-side concern rather than something this app does, but it is what
keeps the whole array path reachable ahead of time: returning the float type as
an ordinary value leaves `_transform_interleaved!` unresolved for every
operator, which is 40 verifier errors and no trimmed binary.

## What does not survive trimming

The build defaults to `--trim=safe`, which fails on any unresolved call. Any
verifier error means the program did not trim, whatever `--trim=unsafe-warn`
may go on to link.

**A preference is needed to get there.** `HostCPUFeatures.__init__` reaches a
`dlopen` of LLVM to read the host's CPU feature string, which the verifier
cannot resolve; nothing here uses that package, it arrives through
VectorizationBase. `LocalPreferences.toml` sets its `freeze_cpu_target`, which
makes the path statically dead — `build.jl` explains why in full, including
why HostCPUFeatures has to be a direct dependency for the setting to apply.

**Threads.** `Threads.@threads` does not work in a trimmed binary: the task
bodies it creates are reached only through the scheduler, so they are not in
the image and the loop dies with a `MethodError` on the first chunk. This is a
`juliac` limitation, not a FastGeoProjections one — a six-line program that
fills an array in parallel fails the same way. So `reproject!` passes
`threaded = false`.

That costs less than it sounds like: see below. If you do want it, the pattern
that works is an explicit `Task` over a *named* callable struct plus
`Base.Experimental.entrypoint(Tuple{YourJobType})`, which gives the trimmer
something to anchor on. Keep the struct monomorphic (dispatch the projection
*inside* the task, not by parameterizing the job) so one declaration covers
every CRS pair.

**Relocatability.** The binary links only `libjulia`, but `libproj_jll`
`dlopen`s its library from an absolute path in the Julia depot at startup. It runs anywhere that depot is; it is not a copy-anywhere artifact.
Shipping one would mean bundling those with `--relative-rpath`.

## Numbers

Measured on an M4 Pro, 1,000,000 random points, Julia 1.12.7.

**Agreement with Proj**, comparing the binary's output against `Proj.jl` over
grids covering each projection's domain:

| conversion | max abs. difference |
| --- | --- |
| 4326 → 3413 | 2.8e-9 m |
| 4326 → 3031 | 1.9e-9 m |
| 4326 → all 120 UTM zones | 5.2e-9 m |
| all 120 UTM zones → 4326 | 5.7e-14 ° |
| 3413 → 4326 | 7.6e-12 ° |

**Speed**, 4326 → 3413, one thread, output to `/dev/null`:

| | wall time |
| --- | --- |
| `fastgeoproj` | 0.18 s |
| `cs2cs` (PROJ 9.8.1) | 1.19 s |

Same values to the digits `cs2cs -d 10` prints. It is not a perfectly matched
comparison — `cs2cs` also emits a third `z` column, so it writes 27% more text
(52.6 MB against 41.4 MB) — but the gap is not close enough for that to
explain it.

Startup, on a one-line file, is 30 ms. The same script run through `julia`
instead of compiled takes 14.4 s for the same work, essentially all of it
package loading and JIT.

**Where the time goes**, in process, 1M points:

| phase | ms |
| --- | --- |
| read the file | 8 |
| parse text → `Float64` | 76 |
| **transform** | **14** |
| format `Float64` → text and write | 44 |

The projection is 9% of the run. Turning text into floats and back is the
expensive half of a CSV tool, which is why the writer goes straight to a byte
buffer with `Ryu.writeshortest` (about 4× an `IOBuffer` a field at a time, same
shortest-round-trip digits) — and why threading the transform, if it were
available, would buy under 10% end to end.
