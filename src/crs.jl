# CRSs described by their parameters rather than by an EPSG code.
#
# `pipeline` reads an EPSG code at run time and returns whichever operator it names, so the type of a
# `Transformation` built from two codes is known only once the codes are. That is fine in a session
# and is a run-time dispatch at every later call in a program compiled ahead of time. A `NativeCRS`
# carries the same information as a type plus its parameters, so `Transformation(source, target)`
# over two of them has a type that follows from theirs, whatever zone or hemisphere they hold.

"""
    NativeCRS

A CRS FastGeoProjections implements natively, given by its parameters: [`GeographicCRS`](@ref),
[`UTMCRS`](@ref) or [`PolarStereographicCRS`](@ref). [`native_crs`](@ref) gives the one an EPSG code
names.

`Transformation(source::NativeCRS, target::NativeCRS)` builds the same operator as
`Transformation(epsg(source), epsg(target); always_xy = true)`, with a type that depends only on the
types of `source` and `target`.
"""
abstract type NativeCRS end

"""
    GeographicCRS()

WGS 84 longitude and latitude in degrees, EPSG:4326, in `(lon, lat)` order.
"""
struct GeographicCRS <: NativeCRS end

"""
    UTMCRS(zone, isnorth)

A WGS 84 UTM zone, EPSG:326zz (north) or EPSG:327zz (south).
"""
struct UTMCRS <: NativeCRS
    zone::Int
    isnorth::Bool
    function UTMCRS(zone::Integer, isnorth::Bool)
        1 <= zone <= 60 || throw(ArgumentError("UTM zone must be in 1:60, got " * string(zone)))
        return new(Int(zone), isnorth)
    end
end

"""
    PolarStereographicCRS(code)

One of the two WGS 84 polar stereographic CRSs FastGeoProjections implements: NSIDC Sea Ice Polar
Stereographic North, EPSG:3413 (true scale at 70°N, central meridian 45°W), or Antarctic Polar
Stereographic, EPSG:3031 (true scale at 71°S, central meridian 0°).
"""
struct PolarStereographicCRS <: NativeCRS
    code::Int
    lat_ts::Float64
    lon_0::Float64
    function PolarStereographicCRS(code::Integer)
        code == 3413 && return new(3413, 70.0, -45.0)
        code == 3031 && return new(3031, -71.0, 0.0)
        throw(ArgumentError("EPSG:" * string(code) * " is not a polar stereographic CRS " *
                            "FastGeoProjections implements; 3413 and 3031 are"))
    end
end

"""
    epsg(c::NativeCRS) -> EPSG

The EPSG code `c` is.
"""
epsg(::GeographicCRS) = EPSG(4326)
epsg(c::UTMCRS) = utmzone2epsg(c.zone, c.isnorth)
epsg(c::PolarStereographicCRS) = EPSG(c.code)

"""
    native_crs(code::Integer) -> NativeCRS

The [`NativeCRS`](@ref) an EPSG code names: 4326, a WGS 84 UTM zone, 3413 or 3031. Any other code is an
error. The result is one of three concrete types; a caller that needs a concrete type downstream
branches on it.
"""
function native_crs(code::Integer)
    code == 4326 && return GeographicCRS()
    32601 <= code <= 32660 && return UTMCRS(code - 32600, true)
    32701 <= code <= 32760 && return UTMCRS(code - 32700, false)
    (code == 3413 || code == 3031) && return PolarStereographicCRS(code)
    throw(ArgumentError("EPSG:" * string(code) * " has no NativeCRS; native_crs knows 4326, the " *
                        "WGS 84 UTM zones, 3413 and 3031"))
end

# The operators to and from EPSG:4326 in (lon, lat) order: the same objects `project_to_4326` and
# `project_from_4326` build for the EPSG code of each.
_to_lonlat(::GeographicCRS, ::Type{T}, kernel) where {T} = Identity()
_to_lonlat(c::UTMCRS, ::Type{T}, kernel) where {T} = UTMToLonLat{T}(c.zone, c.isnorth; kernel)
_to_lonlat(c::PolarStereographicCRS, ::Type{T}, kernel) where {T} =
    PolarStereographicToLonLat{T}(; lat_ts = c.lat_ts, lon_0 = c.lon_0, kernel)

_from_lonlat(::GeographicCRS, ::Type{T}, kernel) where {T} = Identity()
_from_lonlat(c::UTMCRS, ::Type{T}, kernel) where {T} = LonLatToUTM{T}(c.zone, c.isnorth; kernel)
_from_lonlat(c::PolarStereographicCRS, ::Type{T}, kernel) where {T} =
    LonLatToPolarStereographic{T}(; lat_ts = c.lat_ts, lon_0 = c.lon_0, kernel)

"""
    Transformation(source::NativeCRS, target::NativeCRS, T = Float64; threaded = true, kernel = FastKernel())

The transformation from `source` to `target` in `(x, y)` / `(lon, lat)` order, as
`Transformation(epsg(source), epsg(target); always_xy = true, T, kernel)` builds it, and with a type
that depends only on the types of `source` and `target`. The working precision `T` is positional so
that it is part of the method's signature rather than a run-time value.
"""
function Transformation(source::NativeCRS, target::NativeCRS, ::Type{T} = Float64;
                        threaded::Bool = true, kernel::MathKernel = DEFAULT_KERNEL) where {T}
    f = _from_lonlat(target, T, kernel) ∘ _to_lonlat(source, T, kernel)
    return Transformation(f, epsg(source), epsg(target), true, threaded, false)
end
