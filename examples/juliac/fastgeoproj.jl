# fastgeoproj -- reproject a two-column CSV from one CRS to another.
#
# Compiled ahead of time with juliac into a standalone executable: no Julia
# process to start, no compilation at run time.
#
#     fastgeoproj --input points.csv --from 4326 --to 32619 --output utm.csv
#
# The input is strictly two columns of Float64, comma separated, no header.
#
# Two things make this compile under `--trim`, and both are worth copying into
# an app of your own:
#
#   * the CRS pair is resolved by `with_source`/`with_target`, which pass a
#     *concrete* operator to a continuation instead of returning one, so the
#     projection pipeline is a compile-time type rather than a run-time value;
#   * the points live in a `Vector{NTuple{2,Float64}}` -- the layout a CSV of
#     pairs already has -- which is exactly the interleaved buffer
#     FastGeoProjections transforms on SIMD lanes.
#
# See README.md for the build command and for what does not survive trimming.

module FastGeoProjApp

using FastGeoProjections
using FastGeoProjections: Identity, SwapXY,
                          LonLatToPolarStereographic, PolarStereographicToLonLat,
                          LonLatToUTM, UTMToLonLat

const F = Float64
const Point = NTuple{2,F}
const VERSION = "0.1.0"

# ---------------------------------------------------------------------------
# diagnostics
#
# `Core.stdout` and `Core.stderr` are the raw streams. `Base.stdout` is a
# global of abstract type, so writing to it would be a run-time dispatch --
# exactly what a trimmed binary has no method table for.
# ---------------------------------------------------------------------------

@noinline function die(msg::String)
    print(Core.stderr, "fastgeoproj: ", msg, "\n")
    exit(2)
end

@noinline function badinput(path::String, line::Int, msg::String)
    print(Core.stderr, "fastgeoproj: ", path, ":", string(line), ": ", msg, "\n")
    exit(1)
end

const USAGE = """
Usage: fastgeoproj --input FILE --from EPSG --to EPSG [--output FILE]

Reproject a headerless, two-column CSV of Float64 coordinates.

  -i, --input FILE    input CSV
  -o, --output FILE   output CSV (default: stdout)
  -f, --from EPSG     source CRS, as 4326 or EPSG:4326
  -t, --to EPSG       target CRS
      --always-xy     read and write a geographic CRS as (lon, lat) rather
                      than in the authority order (lat, lon)
  -h, --help          show this message
      --version       show the version

Supported EPSG codes:
  4326            WGS 84 geographic
  3031            WGS 84 / Antarctic Polar Stereographic
  3413            WGS 84 / NSIDC Sea Ice Polar Stereographic North
  32601 - 32660   WGS 84 / UTM zones 1N - 60N
  32701 - 32760   WGS 84 / UTM zones 1S - 60S
"""

# ---------------------------------------------------------------------------
# command line
# ---------------------------------------------------------------------------

struct Options
    input::String
    output::String
    from::Int
    to::Int
    always_xy::Bool
end

issupported(code::Int) =
    code == 4326 || code == 3031 || code == 3413 ||
    (32601 <= code <= 32660) || (32701 <= code <= 32760)

function parse_epsg(s::String, flag::String)
    off = (startswith(s, "EPSG:") || startswith(s, "epsg:")) ? 6 : 1
    code = tryparse(Int, SubString(s, off))
    code === nothing && die(string(flag, ": not an EPSG code: ", s))
    issupported(code) ||
        die(string(flag, ": EPSG:", string(code),
                   " is not one of the supported codes (--help lists them)"))
    code
end

function value_of(args::Vector{String}, i::Int, flag::String)
    i <= length(args) || die(string(flag, " needs a value"))
    args[i]
end

function parse_args(args::Vector{String})
    input = ""
    output = "-"
    from = 0
    to = 0
    always_xy = false
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--input" || a == "-i"
            i += 1; input = value_of(args, i, "--input")
        elseif a == "--output" || a == "-o"
            i += 1; output = value_of(args, i, "--output")
        elseif a == "--from" || a == "-f"
            i += 1; from = parse_epsg(value_of(args, i, "--from"), "--from")
        elseif a == "--to" || a == "-t"
            i += 1; to = parse_epsg(value_of(args, i, "--to"), "--to")
        elseif a == "--always-xy"
            always_xy = true
        elseif a == "--help" || a == "-h"
            print(Core.stdout, USAGE)
            exit(0)
        elseif a == "--version"
            print(Core.stdout, VERSION, "\n")
            exit(0)
        else
            die(string("unknown argument: ", a))
        end
        i += 1
    end
    isempty(input) && die("--input is required (--help for usage)")
    from == 0 && die("--from is required (--help for usage)")
    to == 0 && die("--to is required (--help for usage)")
    Options(input, output, from, to, always_xy)
end

# ---------------------------------------------------------------------------
# reading
#
# The file is slurped and scanned as bytes. Fields go to the same C routine
# `parse(Float64, ::String)` uses, so a coordinate that was written by this
# program reads back bit for bit, and nothing is allocated per row.
# ---------------------------------------------------------------------------

const NL = UInt8('\n')
const CR = UInt8('\r')
const COMMA = UInt8(',')
const SPACE = UInt8(' ')
const TAB = UInt8('\t')

@inline istrimmable(b::UInt8) = b == SPACE || b == TAB || b == CR

@inline function field(buf::Vector{UInt8}, lo::Int, hi::Int)
    hi < lo && return nothing
    ok, v = ccall(:jl_try_substrtod, Tuple{Bool,F},
                  (Ptr{UInt8}, Csize_t, Csize_t), buf, lo - 1, hi - lo + 1)
    ok ? v : nothing
end

function slurp(path::String)
    isfile(path) || die(string("no such file: ", path))
    io = open(path, "r")
    try
        read(io)
    finally
        close(io)
    end
end

function read_points(buf::Vector{UInt8}, path::String)
    n = length(buf)
    # one row per newline, plus a last row if the file does not end in one
    rows = 0
    @inbounds for i in 1:n
        buf[i] == NL && (rows += 1)
    end
    (n > 0 && @inbounds(buf[n]) != NL) && (rows += 1)

    pts = Vector{Point}(undef, rows)
    k = 0
    line = 0
    i = 1
    while i <= n
        j = i
        @inbounds while j <= n && buf[j] != NL
            j += 1
        end
        line += 1
        lo, hi = i, j - 1
        @inbounds while lo <= hi && istrimmable(buf[lo]); lo += 1; end
        @inbounds while hi >= lo && istrimmable(buf[hi]); hi -= 1; end
        if hi >= lo                       # blank lines are skipped
            k += 1
            @inbounds pts[k] = read_row(buf, lo, hi, line, path)
        end
        i = j + 1
    end
    resize!(pts, k)
end

function read_row(buf::Vector{UInt8}, lo::Int, hi::Int, line::Int, path::String)
    c = 0
    @inbounds for p in lo:hi
        if buf[p] == COMMA
            c == 0 || badinput(path, line, "more than two columns")
            c = p
        end
    end
    c == 0 && badinput(path, line, "expected two comma-separated columns")
    x = field(buf, lo, c - 1)
    x === nothing && badinput(path, line, "column 1 is not a number")
    y = field(buf, c + 1, hi)
    y === nothing && badinput(path, line, "column 2 is not a number")
    (x, y)
end

# ---------------------------------------------------------------------------
# the pipeline
#
# `FastGeoProjections.Transformation(EPSG(from), EPSG(to))` would build this in
# one call, but it picks the projection from a run-time EPSG code, so the type
# of what it returns is only known at run time -- which `--trim` cannot follow.
#
# Passing the operator *forward* to a continuation instead of returning it
# moves that choice into the type domain: every branch below specializes `k` on
# a concrete operator, so the transform is compiled with the projection fully
# inlined and nothing is dispatched at run time. Three source operators and
# three target ones (plus SwapXY for authority axis order) means the binary
# carries a handful of copies of the transform loop.
# ---------------------------------------------------------------------------

@inline function with_source(k::K, code::Int, always_xy::Bool)::Nothing where {K}
    if code == 4326
        # the native projections speak (lon, lat); EPSG:4326 is latitude first
        # unless the caller asks otherwise
        always_xy ? k(Identity()) : k(SwapXY())
    elseif code == 3031
        k(PolarStereographicToLonLat{F}(; lat_ts = -71.0, lon_0 = 0.0))
    elseif code == 3413
        k(PolarStereographicToLonLat{F}(; lat_ts = 70.0, lon_0 = -45.0))
    elseif 32601 <= code <= 32660
        k(UTMToLonLat{F}(code - 32600, true))
    else
        k(UTMToLonLat{F}(code - 32700, false))
    end
end

@inline function with_target(k::K, code::Int, always_xy::Bool)::Nothing where {K}
    if code == 4326
        always_xy ? k(Identity()) : k(SwapXY())
    elseif code == 3031
        k(LonLatToPolarStereographic{F}(; lat_ts = -71.0, lon_0 = 0.0))
    elseif code == 3413
        k(LonLatToPolarStereographic{F}(; lat_ts = 70.0, lon_0 = -45.0))
    elseif 32601 <= code <= 32660
        k(LonLatToUTM{F}(code - 32600, true))
    else
        k(LonLatToUTM{F}(code - 32700, false))
    end
end

# `threaded = true` works because 1.13 keeps task bodies in a trimmed image; on
# 1.12 this call compiled and then died at run time (see README.md). Threads
# come from `JULIA_NUM_THREADS`, which a trimmed binary reads as usual, so this
# is one thread unless the caller asks for more.
#
# It buys little here: the projection is 5% of the run, and turning text into
# floats and back is the expensive half of a CSV tool. The transform itself
# scales about 5x on eight threads.
function reproject!(pts::Vector{Point}, o::Options)
    with_source(o.from, o.always_xy) do to_lonlat
        with_target(o.to, o.always_xy) do from_lonlat
            transform!(from_lonlat ∘ to_lonlat, pts; threaded = true)
            nothing
        end
    end
end

# ---------------------------------------------------------------------------
# writing
#
# Straight into a byte buffer with Ryu, which is what `print(io, ::Float64)`
# does underneath: same shortest-round-trip digits, about four times the
# throughput of going through an IOBuffer a field at a time.
# ---------------------------------------------------------------------------

const CHUNK = 1 << 16
const MAXFLOAT = Base.Ryu.neededdigits(F)

function write_points(pts::Vector{Point}, path::String)
    if path == "-"
        emit(pts, Core.stdout)
    else
        io = open(path, "w")
        try
            emit(pts, io)
        finally
            close(io)
        end
    end
end

function emit(pts::Vector{Point}, io::O) where {O}
    # room for one more row past the flush point, so a row is never split
    b = Vector{UInt8}(undef, CHUNK + 2 * MAXFLOAT + 2)
    pos = 1
    @inbounds for k in eachindex(pts)
        x, y = pts[k]
        pos = Base.Ryu.writeshortest(b, pos, x)
        b[pos] = COMMA; pos += 1
        pos = Base.Ryu.writeshortest(b, pos, y)
        b[pos] = NL; pos += 1
        if pos > CHUNK
            flush_bytes(io, b, pos - 1)
            pos = 1
        end
    end
    pos > 1 && flush_bytes(io, b, pos - 1)
    nothing
end

@inline flush_bytes(io, b::Vector{UInt8}, n::Int) =
    GC.@preserve b unsafe_write(io, pointer(b), n % UInt)

# ---------------------------------------------------------------------------

function main(args::Vector{String})
    o = parse_args(args)
    pts = read_points(slurp(o.input), o.input)
    reproject!(pts, o)
    write_points(pts, o.output)
    0
end

end # module

function (@main)(args::Vector{String})
    FastGeoProjApp.main(args)
end
