# FastGeoProjections to Proj speed comparison, one panel per CRS pair.
#
# Run as `julia --project=benchmark --threads=8 benchmark/benchmark.jl`. Writes
# `benchmark.png` next to this file, prints a ns/point table at the largest
# point count, and serializes the raw results so a later run can be compared
# without re-measuring.
#
# The transformation is built outside the timed region. Constructing one is a
# one-off cost -- large for Proj, which parses the CRS definitions and resolves
# a pipeline -- and including it measures setup rather than throughput at small
# point counts.

using FastGeoProjections
using FastGeoProjections: Transformation, transform   # DataFrames exports `transform` too
import Proj                                          # a weak dependency: `proj_only` needs it loaded
using BenchmarkTools
using DataFrames
using CairoMakie
using Printf
using Serialization

const OUTFILE = joinpath(@__DIR__, "benchmark")
const NS = [100, 1000, 10_000, 100_000, 1_000_000]
const SYSTEM = "Apple M2 Max"
const SOLUTIONS = ["Proj: single-thread", "Proj: multi-thread",
                   "FGP: single-thread", "FGP: multi-thread"]

# A case is a CRS pair plus the shape its coordinates are passed in: `:soa` for
# two coordinate vectors, `:aos` for a vector of points. A geocentric end needs
# three coordinates, which only the point-vector form carries.
struct Case
    source::EPSG
    target::EPSG
    shape::Symbol
    label::String
end

# Units of the target CRS, which is what an error is quoted in.
units(case::Case) = first(case.target.val) in (4326, 4979) ? "°" : "m"

const CASES = [
    Case(EPSG(4326), EPSG(3413), :soa, "polar stereographic"),
    Case(EPSG(3031), EPSG(4326), :soa, "polar stereographic, inverse"),
    Case(EPSG(4326), EPSG(3857), :soa, "web Mercator"),
    Case(EPSG(4326), EPSG(32636), :soa, "UTM zone 36N"),
    Case(EPSG(32735), EPSG(4326), :soa, "UTM zone 35S, inverse"),
    # A geocentric source is fused into the projection after it, so each of the
    # three projections that opt into `project_direction` is measured that way:
    # what the fusion is worth depends on how much of the projection's own
    # trigonometry it cancels, which differs per projection.
    Case(EPSG(4978), EPSG(3413), :aos, "geocentric, fused"),
    Case(EPSG(4978), EPSG(3857), :aos, "geocentric to web Mercator, fused"),
    Case(EPSG(4978), EPSG(32636), :aos, "geocentric to UTM 36N, fused"),
]

# Longitude and latitude samples in the hemisphere each case covers, then
# projected into the source CRS by Proj -- so the input to a benchmark never
# comes from the code being benchmarked.
function lonlat(case::Case, n::Int)
    code = first((case.source == EPSG(4326) || case.source == EPSG(4978) ?
                  case.target : case.source).val)
    u = range(0.0, 1.0; length = n)
    if code == 3413                      # north polar
        (collect(u .* 360 .- 180), collect(u .* 30 .+ 60))
    elseif code == 3031                  # south polar
        (collect(u .* 360 .- 180), collect(-(u .* 30 .+ 60)))
    elseif code == 3857                  # web Mercator, valid to ±85°
        (collect(u .* 360 .- 180), collect(u .* 170 .- 85))
    elseif code == 32636                 # UTM zone 36, north
        (collect(u .* 6 .+ 30), collect(u .* 80))
    elseif code == 32735                 # UTM zone 35, south
        (collect(u .* 6 .+ 24), collect(u .* -80))
    else
        error("no sampling defined for EPSG:$code")
    end
end

# A topographic range of ellipsoidal heights, which is what a geocentric
# conversion is asked for in practice.
heights(n::Int) = collect(range(-400.0, 8000.0; length = n))

"""
    inputs(case, n)

Coordinates in `case.source`, in the shape `case.shape` names, and a callable
`run(trans)` that applies a transformation to them.
"""
function inputs(case::Case, n::Int)
    lon, lat = lonlat(case, n)
    if case.shape === :soa
        to_source = Transformation(EPSG(4326), case.source;
                                   proj_only = true, always_xy = true)
        X, Y = case.source == EPSG(4326) ? (lon, lat) : to_source(lon, lat)
        return ((X, Y), trans -> transform(trans, X, Y; threaded = trans.threaded))
    else
        h = heights(n)
        to_source = Transformation(EPSG(4979), case.source;
                                   proj_only = true, always_xy = true)
        pts = [to_source(lo, la, hh) for (lo, la, hh) in zip(lon, lat, h)]
        return (pts, trans -> transform(trans, pts; threaded = trans.threaded))
    end
end

# Maximum absolute difference from the Proj result over the projected
# coordinates, in the units of the target CRS.
#
# Only the first two components, even where the pipeline returns three. A
# geocentric source hands back an ellipsoidal height, and PROJ's own geocentric
# inverse loses accuracy in that height with altitude -- 8.8e-7 m over a
# topographic range against 3.9e-9 m here -- so including it would report PROJ's
# error as this package's. The test suite pins the height against an exactly
# computed position instead.
#
# A difference in degrees is reduced to (-180, 180]: at the ±180° branch cut the
# two implementations may land on opposite representations of the same meridian,
# and an unreduced difference of 360° there swamps every real one.
function maxerror(case::Case, got, reference)
    wrap = units(case) == "°"
    m = 0.0
    for (a, b) in zip(got, reference)
        d = (a[1] - b[1], a[2] - b[2])
        wrap && (d = mod.(d .+ 180, 360) .- 180)
        m = max(m, maximum(abs.(d)))
    end
    m
end

flatten(r::Tuple) = collect(zip(r...))   # (X, Y) from the SoA path
flatten(r::AbstractVector) = r

function measure(run, reference, case::Case, solution::String)
    proj_only = startswith(solution, "Proj")
    threaded = endswith(solution, "multi-thread")
    trans = Transformation(case.source, case.target;
                           threaded, proj_only, always_xy = true)
    b = @benchmark $run($trans) seconds = 2
    (time = minimum(b).time, err = maxerror(case, flatten(run(trans)), reference))
end

function main()
    df = DataFrame()
    for solution in SOLUTIONS
        df[!, solution * "_time"] = zeros(length(CASES) * length(NS))
        df[!, solution * "_err"] = zeros(length(CASES) * length(NS))
    end
    df[!, :npoints] = zeros(Int, length(CASES) * length(NS))
    df[!, :case] = zeros(Int, length(CASES) * length(NS))

    for (k, case) in enumerate(CASES), (i, n) in enumerate(NS)
        r = i + (k - 1) * length(NS)
        df[r, :npoints] = n
        df[r, :case] = k
        printstyled("EPSG:$(first(case.source.val)) => EPSG:$(first(case.target.val))",
                    " [n = $n]\n"; color = :blue)
        _, run = inputs(case, n)
        proj = Transformation(case.source, case.target;
                              threaded = false, proj_only = true, always_xy = true)
        reference = flatten(run(proj))
        for solution in SOLUTIONS
            m = measure(run, reference, case, solution)
            df[r, solution * "_time"] = m.time
            df[r, solution * "_err"] = m.err
            @printf("  %-22s %10.1f µs  %8.2f ns/pt   ME %.2e\n",
                    solution, m.time / 1000, m.time / n, m.err)
        end
    end

    println("\nns/point at n = $(maximum(NS)):\n")
    @printf("%-38s %9s %9s %9s %9s  %8s\n", "pipeline", "Proj 1t", "Proj Nt",
            "FGP 1t", "FGP Nt", "ME")
    for (k, case) in enumerate(CASES)
        r = findfirst((df.case .== k) .& (df.npoints .== maximum(NS)))
        @printf("%-38s %9.1f %9.1f %9.1f %9.1f  %8.1e %s\n",
                "$(first(case.source.val))→$(first(case.target.val)) $(case.label)",
                (df[r, s * "_time"] / maximum(NS) for s in SOLUTIONS)...,
                df[r, SOLUTIONS[3] * "_err"], units(case))
    end

    serialize("$OUTFILE.jls", (df = df, threads = Threads.nthreads(), system = SYSTEM))
    figure(df)
end

function figure(df)
    nrows = ceil(Int, length(CASES) / 2)
    f = Figure(size = (1500, 500 * nrows), fontsize = 22)
    col = Makie.wong_colors()
    for (i, case) in enumerate(CASES)
        r = ceil(Int, i / 2)
        c = i - 2 * (r - 1)
        rs = df.case .== i
        err = df[findfirst(rs .& (df.npoints .== maximum(NS))), SOLUTIONS[3] * "_err"]
        ax = Axis(f[r, c];
            xscale = log10, yscale = log10,
            title = "EPSG:$(first(case.source.val)) => EPSG:$(first(case.target.val)) — " *
                    "$(case.label)\nME = $(@sprintf("%.1e", err)) $(units(case))",
            titlesize = 22,
            xlabel = "points converted",
            ylabel = "compute time [µs]",
            yminorticksvisible = true, yminorgridvisible = true,
            yminorticks = IntervalsBetween(5),
        )
        lins = [lines!(ax, df[rs, :npoints], df[rs, s * "_time"] ./ 1000;
                       linewidth = 6, color = col[j])
                for (j, s) in enumerate(SOLUTIONS)]
        i == 1 && axislegend(ax, lins, SOLUTIONS; position = :lt, labelsize = 18)
    end
    Label(f[0, :], "FastGeoProjections.jl, $SYSTEM, $(Threads.nthreads()) threads";
          fontsize = 28)
    save("$OUTFILE.png", f)
    println("\nwrote $OUTFILE.png")
end

main()
