# Universal Transverse Mercator.
#
# UTM is transverse Mercator on the WGS 84 ellipsoid with a fixed scale factor
# and false origin, and a central meridian set by the zone. Rather than compose
# those onto the projection, both directions are their own operator: the zone
# is then a field you can read back, and the scale and false origin are folded
# into constants at construction.

const UTM_K0 = 0.9996          # scale factor on the central meridian
const UTM_FE = 5e5             # false easting
const UTM_FN = 1e7             # false northing, southern hemisphere only

"""
    utm_lon0(zone)

Central meridian of a UTM zone, in decimal degrees.
"""
utm_lon0(zone::Integer) = -183 + 6 * zone

_utm_check(zone::Integer) =
    1 <= zone <= 60 || throw(ArgumentError("UTM zone must be in 1:60, got $zone"))

"""
    LonLatToUTM(zone, isnorth; T = Float64, kernel = FastKernel())
    LonLatToUTM{T}(zone, isnorth; kernel = FastKernel())

Point operator taking geodetic `(lon, lat)` in decimal degrees to UTM
`(easting, northing)` in metres, in the given `zone` (1 to 60) and hemisphere.

    julia> t = LonLatToUTM(19, true)
    LonLatToUTM{Float64}(zone = 19, north)

    julia> t(-69.0, 45.0)
    (500000.0, 4.982950400226553e6)

[`convergence_scale`](@ref) gives the grid convergence and point scale, the
latter including the `0.9996` zone scale factor.
"""
struct LonLatToUTM{T,K<:MathKernel} <: GeoTransformation
    zone::Int
    isnorth::Bool
    tm::LonLatToTransverseMercator{T,K}
    k0::T
    fe::T
    fn::T
end

"""
    UTMToLonLat(zone, isnorth; T = Float64, kernel = FastKernel())
    UTMToLonLat{T}(zone, isnorth; kernel = FastKernel())

Point operator taking UTM `(easting, northing)` in metres back to geodetic
`(lon, lat)` in decimal degrees. See [`LonLatToUTM`](@ref).
"""
struct UTMToLonLat{T,K<:MathKernel} <: GeoTransformation
    zone::Int
    isnorth::Bool
    tm::TransverseMercatorToLonLat{T,K}
    inv_k0::T          # the false origin is removed before the projection, so
    dx::T              # these carry -fe/k0 and -fn/k0 already divided through
    dy::T
end

# `{T}` is the primary form: with the working precision a type parameter rather
# than a `Type`-valued keyword, the constructed operator has a concrete type, so
# a caller that picks a zone at run time still gets a statically typed pipeline.
# The keyword form forwards to it and keeps the shorthand.
function LonLatToUTM{T}(zone::Integer, isnorth::Bool;
                        kernel::MathKernel = DEFAULT_KERNEL) where {T}
    _utm_check(zone)
    tm = LonLatToTransverseMercator{T}(; lon0 = utm_lon0(zone), lat0 = 0, kernel)
    LonLatToUTM{T,typeof(kernel)}(zone, isnorth, tm, T(UTM_K0), T(UTM_FE),
                                  isnorth ? zero(T) : T(UTM_FN))
end

function UTMToLonLat{T}(zone::Integer, isnorth::Bool;
                        kernel::MathKernel = DEFAULT_KERNEL) where {T}
    _utm_check(zone)
    tm = TransverseMercatorToLonLat{T}(; lon0 = utm_lon0(zone), lat0 = 0, kernel)
    fn = isnorth ? zero(T) : T(UTM_FN)
    UTMToLonLat{T,typeof(kernel)}(zone, isnorth, tm, T(inv(UTM_K0)),
                                  T(-UTM_FE / UTM_K0), T(-fn / UTM_K0))
end

LonLatToUTM(zone::Integer, isnorth::Bool; T::Type = Float64, kwargs...) =
    LonLatToUTM{T}(zone, isnorth; kwargs...)
UTMToLonLat(zone::Integer, isnorth::Bool; T::Type = Float64, kwargs...) =
    UTMToLonLat{T}(zone, isnorth; kwargs...)

LonLatToUTM(epsg::EPSG; kwargs...) =
    (z = epsg2utmzone(epsg); LonLatToUTM(z.zone, z.isnorth; kwargs...))
UTMToLonLat(epsg::EPSG; kwargs...) =
    (z = epsg2utmzone(epsg); UTMToLonLat(z.zone, z.isnorth; kwargs...))
LonLatToUTM{T}(epsg::EPSG; kwargs...) where {T} =
    (z = epsg2utmzone(epsg); LonLatToUTM{T}(z.zone, z.isnorth; kwargs...))
UTMToLonLat{T}(epsg::EPSG; kwargs...) where {T} =
    (z = epsg2utmzone(epsg); UTMToLonLat{T}(z.zone, z.isnorth; kwargs...))

@inline function (t::LonLatToUTM)(lon, lat)
    x, y = t.tm(lon, lat)
    (muladd(x, t.k0, t.fe), muladd(y, t.k0, t.fn))
end

@inline (t::UTMToLonLat)(x, y) =
    t.tm(muladd(x, t.inv_k0, t.dx), muladd(y, t.inv_k0, t.dy))

"""
    convergence_scale(t::LonLatToUTM, lon, lat) -> (γ, k)
    convergence_scale(t::UTMToLonLat, x, y) -> (γ, k)

Grid convergence and point scale of a UTM zone. `k` includes the zone's
`0.9996` scale factor, so it is `0.9996` on the central meridian and rises to
about `1.0010` at the zone edge.
"""
@inline function convergence_scale(t::LonLatToUTM, lon, lat)
    gam, k = convergence_scale(t.tm, lon, lat)
    (gam, k * t.k0)
end

@inline function convergence_scale(t::UTMToLonLat, x, y)
    gam, k = convergence_scale(t.tm, muladd(x, t.inv_k0, t.dx), muladd(y, t.inv_k0, t.dy))
    (gam, k * inv(t.inv_k0))
end

Base.inv(t::LonLatToUTM{T}) where {T} =
    UTMToLonLat(t.zone, t.isnorth; T, kernel = t.tm.kernel)
Base.inv(t::UTMToLonLat{T}) where {T} =
    LonLatToUTM(t.zone, t.isnorth; T, kernel = t.tm.kernel)

islanesafe(t::LonLatToUTM) = islanesafe(t.tm)
islanesafe(t::UTMToLonLat) = islanesafe(t.tm)

adapt_eltype(t::LonLatToUTM{T}, ::Type{S}) where {T,S} =
    T === S ? t : LonLatToUTM(t.zone, t.isnorth; T = S, kernel = t.tm.kernel)
adapt_eltype(t::UTMToLonLat{T}, ::Type{S}) where {T,S} =
    T === S ? t : UTMToLonLat(t.zone, t.isnorth; T = S, kernel = t.tm.kernel)

Base.show(io::IO, t::LonLatToUTM{T}) where {T} =
    print(io, "LonLatToUTM{$T}(zone = ", t.zone, ", ", t.isnorth ? "north" : "south", ")")
Base.show(io::IO, t::UTMToLonLat{T}) where {T} =
    print(io, "UTMToLonLat{$T}(zone = ", t.zone, ", ", t.isnorth ? "north" : "south", ")")

# ---------------------------------------------------------------------------
# zone <-> EPSG
# ---------------------------------------------------------------------------

"""
    utm_epsg(lon::Real, lat::Real, always_xy = true)

The EPSG code of the UTM zone containing `(lon, lat)` -- or of the relevant
polar stereographic projection outside the UTM latitude limits.

modified from: https://github.com/JuliaGeo/Geodesy.jl/blob/master/src/utm.jl
"""
function utm_epsg(lon::Real, lat::Real, always_xy = true)

    if !always_xy
        lat, lon = (lon, lat)
    end

    if lat > 84
        # WGS 84 / Arctic Polar Stereographic
        return EPSG(3995)
    elseif lat < -80
        # Antarctic Polar Stereographic
        return EPSG(19992)
    end

    # make sure lon is from -180 to 180
    lon = lon - floor((lon + 180) / (360)) * 360

    # int versions
    ilat = floor(Int64, lat)
    ilon = floor(Int64, lon)

    # get the latitude band
    band = max(-10, min(9, fld((ilat + 80), 8) - 10))

    # and check for weird ones
    zone = fld((ilon + 186), 6)
    if ((band == 7) && (zone == 31) && (ilon >= 3)) # Norway
        zone = 32
    elseif ((band == 9) && (ilon >= 0) && (ilon < 42)) # Svalbard
        zone = 2 * fld((ilon + 183), 12) + 1
    end

    EPSG(lat >= 0 ? 32600 + zone : 32700 + zone)
end

"""
    utmzone2epsg(zone = 0, isnorth = true)

EPSG code of a UTM zone; zone `0` selects the polar projection of the given
hemisphere.
"""
function utmzone2epsg(zone::Int = 0, isnorth::Bool = true)
    if zone == 0
        return isnorth ? EPSG(3995) : EPSG(19992)
    end
    EPSG(isnorth ? 32600 + zone : 32700 + zone)
end

"""
    epsg2utmzone(epsg)

`(zone, isnorth)` of a UTM EPSG code.
"""
function epsg2utmzone(epsg::EPSG)
    code = first(epsg.val)
    if code == 3995
        (zone = 0, isnorth = true)
    elseif code == 19992
        (zone = 0, isnorth = false)
    elseif 32601 <= code <= 32660
        (zone = code - 32600, isnorth = true)
    elseif 32701 <= code <= 32760
        (zone = code - 32700, isnorth = false)
    else
        error("EPSG:$code is not a UTM EPSG code")
    end
end

"""
    isutm(epsg)

Whether `epsg` is a WGS 84 UTM zone.
"""
function isutm(epsg::EPSG)
    code = first(epsg.val)
    (32601 <= code <= 32660) || (32701 <= code <= 32760)
end
