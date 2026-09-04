[![Build Status](https://github.com/alex-s-gardner/FastGeoProjections.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/alex-s-gardner/FastGeoProjections.jl/actions/workflows/CI.yml?query=branch%3Amain)

**FastGeoProjections** is intended to provide highly optimized native Julia geospatial coordinate transformations from one coordinate reference system (CRS) to another as defined by EPSG codes. It is not intended to replace, nor to be as comprehensive as, [Proj](https://github.com/JuliaGeo/Proj.jl). The package will natively support only the most common geospatial transformations and relies on **Proj.jl** for all others.

*Supported Projection EPSGs*
- 3031:     WGS 84 / Antarctic Polar Stereographic
- 3413:     WGS 84 / NSIDC Sea Ice Polar Stereographic North
- 4326:     WGS84 - World Geodetic System 1984
- 326XX:    WGS 84 / UTM zone XXN
- 327XX:    WGS 84 / UTM zone XXS

*Example*
```julia
julia> using Pkg; Pkg.add("FastGeoProjections")
julia> using FastGeoProjections
julia> lat = [84.0, 83.0]; lon = [50.0, 51.0];
julia> trans = FastGeoProjections.Transformation(EPSG(4326), EPSG(3413))
Transformation
    source_epsg:    EPSG:4326
    target_epsg:    EPSG:3413
    threaded:       true
    always_xy:      false
    proj_only:      false
julia> x, y = trans(lat, lon)
([648059.0510298966, 755038.7580833684], [56697.82026048413, 79357.7712642986])
```

*Transformations are point operators*

A `Transformation` is a callable struct that maps **one point to one point**, holding
everything the projection can precompute. The same operator is what gets applied to whole
collections, so there is one implementation of the projection rather than one per calling
convention:

```julia
julia> trans = FastGeoProjections.Transformation(EPSG(4326), EPSG(3413); always_xy=true)

julia> trans(-45.0, 70.0)                   # two coordinates
(0.0, -2.187927649279021e6)

julia> points = [(-45.0, 70.0), (-44.0, 71.0)];
julia> transform(trans, points)             # a vector of points
julia> transform!(trans, points)            # ...or in place
julia> transform(trans, lon, lat)           # ...or as two coordinate vectors
```

A single argument is read as a [GeoInterface](https://github.com/JuliaGeo/GeoInterface.jl)
point, so anything with `PointTrait` works -- including the plain tuples the array API
uses. Only `x` and `y` are used; `z` is ignored.

```julia
julia> trans((-45.0, 70.0)) == trans(GI.Point(-45.0, 70.0)) == trans(-45.0, 70.0)
true
```

Vectors of GeoInterface points go the same way, and at full speed. A vector of isbits
points is already a dense interleaved buffer of coordinates, so the SIMD path reads and
writes it in place, column by column -- `SVector{2}`, `GeometryBasics.Point2`, a two-field
struct of your own and `NTuple{2}` all run at the same ns/point:

```julia
julia> transform(trans, [Point2(-45.0, 70.0), Point2(-44.0, 71.0)])
julia> transform!(trans, [SVector(-45.0, 70.0), SVector(-44.0, 71.0)])
```

A third component costs nothing (the load deinterleaves in hardware) and is carried
through untouched, so `Point3` and `NTuple{3}` keep their `z`. Whether that layout really
holds is settled from the point *type*, by laying sentinel coordinates out in memory and
asking the resulting point where its x and y are. GeoInterface cannot answer it:
`getcoord` may read out of a C API, with no Julia-side memory to be ordered, and
`coordnames` names coordinates rather than describing storage. A point type that stores
its components in the other order is transformed one at a time instead -- correctly,
about half as fast. Such a type needs a `FastGeoProjections.rebuildpoint(::Type{P}, x, y)`
method, since the scalar path has to construct each result.

Pipelines are flat. `f ∘ g ∘ h` is one `ComposedGeoTransformation` holding a tuple of
stages in application order, not a tree of nested pairs, so a pipeline of any length
inlines into a single pass over the data:

```julia
julia> FastGeoProjections.Transformation(EPSG(4326), EPSG(32619)).f
LonLatToUTM{Float64}(zone = 19, north) ∘ FastGeoProjections.SwapXY()
```

Because the operator is generic over its scalar type, it is evaluated on SIMD lanes rather
than one point at a time, which is what makes the point-by-point form as fast as a hand
vectorized array kernel. Threading is a property of the *application*, not of the
transformation, and can be set per call:

```julia
julia> transform!(trans, points; threaded=false)
```

`transform`/`transform!` fall back to a plain scalar loop for any transformation that
cannot be vectorized -- one backed by Proj, or one built with `kernel=BaseKernel()`.

*Projections on their own*

The projections are usable directly, without going through an EPSG pair:

Each is named for what it does, source to target:

```julia
julia> tm = LonLatToTransverseMercator(; lon0 = -69.0, lat0 = 0.0)
julia> tm(-68.0, 45.0)
(78846.84165337226, 4.985430940725587e6)

julia> inv(tm)                                 # TransverseMercatorToLonLat
julia> convergence_scale(tm, -68.0, 45.0)      # meridian convergence [°], point scale
(0.7071430455192697, 1.0000764118961942)
```

UTM is its own operator, carrying the zone and hemisphere as fields and folding the
zone scale factor and false origin into its constants:

```julia
julia> u = LonLatToUTM(19, true)
LonLatToUTM{Float64}(zone = 19, north)

julia> u(-69.0, 45.0)
(500000.0, 4.982950400226553e6)

julia> convergence_scale(u, -69.0, 45.0)       # k includes the 0.9996 zone factor
(0.0, 0.9995999999999997)
```

*Math kernels*

The transcendental functions a projection uses live in a `Math` submodule and take the
back-end as their first argument, so one implementation serves all of them:

```julia
julia> Math.sin(FastKernel(), 0.5)
```

`FastKernel` (the default) is SLEEFPirates' `*_fast` family -- what `@turbo` silently
lowered to. `SLEEFKernel` is the same routines with full argument reduction: sub-ULP
anywhere, at about twice the cost, and identical in practice over the ranges a
projection produces. `BaseKernel` is libm, which does not vectorize.

*Ahead-of-time compilation*

The projections compile with [`juliac`](https://docs.julialang.org/en/v1/devdocs/build/juliac/)
into a standalone executable -- no Julia startup, no JIT.
[`examples/juliac`](examples/juliac) is a CSV reprojection tool built that way;
its README covers the two patterns that make an app on this package trimmable,
and what does not survive `--trim`.

```console
$ ./fastgeoproj --input points.csv --from 4326 --to 32619 --always-xy
```

*Benchmark*

ME = Maximum Error

![benchmark](benchmark/benchmark.jpg)

**Note**
If you have recommendations for additional projections to support feel free to submit a an issue