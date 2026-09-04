# Throughput harness for the optimization pass: one row per (operator, call shape),
# reported as ns/point so array and scalar paths are comparable.
#
# Run as `julia --project=. --threads=8 benchmark/optimize_bench.jl [label]`.
# Results and the numerical snapshot are serialized next to this file so a later
# run can be compared without re-measuring the baseline.

using FastGeoProjections
using BenchmarkTools
using Serialization
using Printf
using StaticArrays

const N = 100_000
const LABEL = isempty(ARGS) ? "run" : ARGS[1]
const OUT = joinpath(@__DIR__, "optimize_$(LABEL).jls")

lons = collect(range(-179.0, 179.0; length = N))
lats = collect(range(-84.0, 84.0; length = N))
hs = collect(range(-400.0, 8000.0; length = N))

# every case is (name, thunk, npoints); thunks close over interpolated data
function cases()
    utm = LonLatToUTM(19, true)
    ps = LonLatToPolarStereographic(; lat_ts = 70.0, lon_0 = -45.0)
    geo = LonLatToGeocentric()
    igeo = GeocentricToLonLat()
    fused = Transformation(EPSG(4978), EPSG(3413); always_xy = true).f

    xy2 = [(lo, la) for (lo, la) in zip(lons, lats)]
    xy3 = [(lo, la, h) for (lo, la, h) in zip(lons, lats, hs)]
    sv3 = [SVector(lo, la, h) for (lo, la, h) in zip(lons, lats, hs)]
    ecef = [geo(lo, la, h) for (lo, la, h) in zip(lons, lats, hs)]

    [
     # 2-D operators, array path: the SIMD reference point
     ("utm  array 2d", () -> transform(utm, xy2), N),
     ("ps   array 2d", () -> transform(ps, xy2), N),
     ("utm  soa 2d", () -> transform(utm, lons, lats), N),
     # height-transforming operators, array path
     ("geo  array 3d", () -> transform(geo, xy3), N),
     ("geo  array 3d SVector", () -> transform(geo, sv3), N),
     ("igeo array 3d", () -> transform(igeo, ecef), N),
     ("fused array 3d", () -> transform(fused, ecef), N),
     # scalar, for reference
     ("utm  scalar", () -> utm(-69.0, 45.0), 1),
     ("ps   scalar", () -> ps(-45.0, 70.0), 1),
     ("geo  scalar", () -> geo(5.39, 52.16, 100.0), 1),
     ("igeo scalar", () -> igeo(3.9036404612786868e6, 368315.27616670664, 5.013823349039822e6), 1),
     ("fused scalar", () -> fused(3.9036404612786868e6, 368315.27616670664, 5.013823349039822e6), 1),
    ]
end

# A digest of what each case computed, so a later run can be checked for
# equivalence rather than only for speed.
function snapshot()
    geo = LonLatToGeocentric()
    igeo = GeocentricToLonLat()
    fused = Transformation(EPSG(4978), EPSG(3413); always_xy = true).f
    utm = LonLatToUTM(19, true)
    ps = LonLatToPolarStereographic(; lat_ts = 70.0, lon_0 = -45.0)
    idx = 1:997:N
    Dict(
        "utm" => [utm(lons[i], lats[i]) for i in idx],
        "ps" => [ps(lons[i], lats[i]) for i in idx],
        "geo" => [geo(lons[i], lats[i], hs[i]) for i in idx],
        "igeo" => [igeo(geo(lons[i], lats[i], hs[i])...) for i in idx],
        "fused" => [fused(geo(lons[i], lats[i], hs[i])...) for i in idx],
    )
end

function main()
    results = Dict{String,NamedTuple{(:ns_per_point, :allocs, :bytes),Tuple{Float64,Int,Int}}}()
    for (name, f, np) in cases()
        f()                                     # warm up
        # `evals=1` on a sub-100 ns scalar call measures the timer, not the
        # call, so a scalar case is evaluated many times per sample.
        b = np == 1 ? (@benchmark $f() samples=200 evals=1000) :
                      (@benchmark $f() samples=200 evals=1)
        m = minimum(b)
        results[name] = (ns_per_point = m.time / np, allocs = m.allocs, bytes = m.memory)
        @printf("%-24s %9.2f ns/pt  %6d allocs  %9d B\n", name, m.time / np, m.allocs, m.memory)
    end
    serialize(OUT, (results = results, snapshot = snapshot(),
                    threads = Threads.nthreads()))
    println("\nwrote ", OUT, "  (", Threads.nthreads(), " threads)")
end

main()
