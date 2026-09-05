# `fastgeoproj` — a standalone CSV reprojection tool

A worked example of compiling FastGeoProjections ahead of time with
[`juliac`](https://github.com/JuliaLang/julia/tree/master/contrib/juliac) into
an executable with no Julia startup and no compilation at run time.

> **Julia 1.13 or later.** `juliac` now ships as the [JuliaC
> package](https://github.com/JuliaLang/JuliaC.jl) rather than as part of the
> distribution, so it is a dependency of this project. More to the point, 1.13
> is the first version that keeps task bodies in a trimmed image, which makes it
> the first version where the threaded transform survives. On 1.12 the build
> succeeds, reports zero errors, and the binary dies at run time — see [what
> trimming does and does not reach](#what-trimming-does-and-does-not-reach).

```console
$ ./fastgeoproj --input sample.csv --from 4326 --to 32619 --always-xy
500000.0,4.982950400226553e6
380244.57633207337,4.900734380124453e6
557548.832534876,5.149876586070551e6
426642.267622362,4.761207746467653e6
319952.79600128037,4.735400315222241e6
```

Input is strictly two columns of `Float64`, comma separated, no header. Blank
lines are skipped; `\r\n`, surrounding spaces and a missing final newline are
all accepted. Anything else is an error naming the line.

## Build

Needs **Julia 1.13 or later** and a C compiler on `PATH` for the final link.

```console
$ julia --project=examples/juliac examples/juliac/build.jl
```

That produces `examples/juliac/fastgeoproj` (about 15 MiB) in about half a
minute. The build fails on any unresolved call; `build.jl --trim=unsafe-warn`
downgrades those to warnings and links anyway, which is only useful for seeing
the whole list at once.

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

Supported codes are 4326, 3031, 3413, 3857, and the 120 WGS 84 UTM zones
(32601–32660, 32701–32760). Any pair of them works, composed through EPSG:4326:
`--from 32619 --to 32620` is one pass over the data, not two. Anything else is
an error naming the code.

That is the two-dimensional part of what the package implements natively. It
also has EPSG:4978 (geocentric) and 4979 (geographic 3D), which this tool does
not offer because its CSV is exactly two columns; a third would be a real
addition rather than another branch in `with_source`/`with_target`.

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
with the projection fully inlined and nothing dispatches at run time. Four
source operators and four target ones means a handful of copies of the loop in
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

## What trimming does and does not reach

The build defaults to `--trim=safe`, which fails on any unresolved call. Any
verifier error means the program did not trim, whatever `--trim=unsafe-warn`
may go on to link.

**Threads work here, and did not on 1.12.** A task's function is stored in the
task object by `jl_new_task` and invoked later by the scheduler, from C. No
Julia call site refers to it — the optimized IR of a function that spawns work
contains `jl_new_task`, `enq_work` and `_wait`, and no edge to the body at all
— so a reachability walk over the visible call graph never finds it. On 1.12
the body was left out of the image and the loop died with a `MethodError` on
the first chunk.

That is the one failure trimming will not warn you about. There is no
unresolved *call site*, because there is no call site, so `--trim=safe` reports
zero errors and links a binary that cannot run its own threaded loop. The
verifier sees Julia; it cannot see through C.

1.13 roots lowering-generated closures used as task bodies. That covers
`Threads.@threads`, `Threads.@spawn`, `StableTasks.@spawn` and a plain
`Task(() -> ...)` alike — they all hand `Task` a closure, and which macro built
it makes no difference. The exception is a *named* callable struct passed
straight to `Task`, which is rooted only if you declare
`Base.Experimental.entrypoint(Tuple{YourJobType})` yourself. That declaration is
also the only thing that worked on 1.12, where nothing was rooted automatically;
keep such a struct monomorphic (dispatch the projection *inside* the task rather
than parameterizing the job) so one declaration covers every CRS pair.

**One workaround this example no longer needs.** On 1.12 the build could not
resolve a `dlopen` of libLLVM inside `HostCPUFeatures.__init__`, reached through
VectorizationBase, and the fix was a `freeze_cpu_target` preference plus a
direct dependency to make that preference apply. On 1.13 the build is clean
without any of it. `build.jl` records what it was, since the shape of the
problem outlives this instance of it.

**Relocatability.** The binary links only `libjulia`, but `libproj_jll`
`dlopen`s its library from an absolute path in the Julia depot at startup. It runs anywhere that depot is; it is not a copy-anywhere artifact.
Shipping one would mean bundling those with `--relative-rpath`.

## Numbers

Measured on an M4 Pro, 1,000,000 random points, Julia 1.13.0-rc4.

**Agreement with Proj**, over grids covering each projection's domain — the
polar caps at 1° of latitude by 5° of longitude, each UTM zone at ±3° of its
central meridian by 3° of latitude:

| conversion | max abs. difference |
| --- | --- |
| 4326 → 3413 | 2.2e-9 m |
| 4326 → 3031 | 2.1e-9 m |
| 4326 → all 120 UTM zones | 5.5e-9 m |
| all 120 UTM zones → 4326 | 5.7e-14 ° |
| 3413 → 4326 | 1.3e-12 ° |
| 3031 → 4326 | 1.4e-12 ° |

**Speed**, 4326 → 3413, one thread, output to `/dev/null`:

| | wall time |
| --- | --- |
| `fastgeoproj` | 0.18 s |
| `cs2cs` (PROJ 9.8.1) | 1.22 s |

The two agree to 3.6e-9 m over the million points. It is not a perfectly
matched comparison — `cs2cs` also emits a third `z` column, so it writes 28%
more text (50.6 MB against 39.6 MB) — but the gap is not close enough for that
to explain it.

Startup, on a one-line file, is 30 ms. The same script run through `julia`
instead of compiled takes 14.4 s for the same work, essentially all of it
package loading and JIT.

**Threads.** A trimmed binary reads `JULIA_NUM_THREADS` as usual and defaults to
one, so the threading is there to be asked for. End to end it is close to
invisible — 0.20 s on one thread against 0.17 s on eight, most of that spread
being run-to-run noise — because only the transform is threaded and the
transform is 5% of the run. In isolation it scales properly:

| threads | `transform!`, 1e6 points | |
| --- | --- | --- |
| 1 | 8.0 ms | |
| 2 | 4.9 ms | 1.7× |
| 4 | 2.3 ms | 3.5× |
| 8 | 1.6 ms | 5.1× |

Output is byte-identical at every thread count, and identical to what the 1.12
single-threaded build produced.

**Where the time goes**, in process, 1M points:

| phase | ms |
| --- | --- |
| read the file | 8 |
| parse text → `Float64` | 76 |
| **transform** | **8** |
| format `Float64` → text and write | 44 |

The projection is 5% of the run, and was 9% before the polar stereographic
forward traded its `pow` for a series — 14.2 ms against 7.4 ms measured back to
back on 1.12, where the change landed. Turning text into floats and back is the
expensive half of a CSV tool, which is why the writer goes straight to a byte
buffer with `Ryu.writeshortest` (about 4× an `IOBuffer` a field at a time, same
shortest-round-trip digits) — and why threading the transform, now that it
survives trimming, still buys under 5% end to end.
