# Geodetic to geocentric Cartesian, as a pair of point operators.
#
# These are the only operators here that transform a height rather than carry it
# across: a geocentric z is a Cartesian coordinate, and x and y depend on the
# height going in. So both directions take and return three coordinates, and
# `preservesz` is false for them -- see `transformations.jl`.
#
# The inverse is Vermeille's 2002 closed form. Geodetic latitude from a
# Cartesian position has no closed form in the usual sense; this one is exact in
# the sense that matters, recovering an exactly-computed position to 2e-9 m
# without iterating, which is what keeps it branch-free and lane-safe. PROJ's own
# inverse loses accuracy with height -- 4.0e-3 m at 700 km -- so the round trip
# rather than PROJ is what pins this direction.

"""
    LonLatToGeocentric(; ellips = ellipsoid(EPSG(7030)), kernel = FastKernel())
    LonLatToGeocentric{T}(; ...)

Point operator taking geodetic `(lon, lat, height)` -- degrees and metres -- to
geocentric Cartesian `(x, y, z)` in metres. This is EPSG:4979 to EPSG:4978 on
WGS 84.

Unlike a map projection, this transforms the height rather than carrying it
across, so it takes three coordinates and returns three. A two-coordinate call
would mean a height of zero, which moves x and y by metres, so it is an error
rather than an implied sea-level point.

    julia> t = LonLatToGeocentric()
    LonLatToGeocentric{Float64}(WGS_84)

    julia> t(5.39, 52.16, 100.0)
    (3.9036404612786868e6, 368315.27616670664, 5.013823349039822e6)
"""
struct LonLatToGeocentric{T,K<:MathKernel} <: GeoTransformation
    a::T
    e2::T
    one_minus_e2::T   # the polar term of the prime vertical radius, folded once
    d2r::T
    name::Union{Nothing,Symbol}
    kernel::K
end

"""
    GeocentricToLonLat(; ellips = ellipsoid(EPSG(7030)), kernel = FastKernel())
    GeocentricToLonLat{T}(; ...)

Point operator taking geocentric Cartesian `(x, y, z)` in metres to geodetic
`(lon, lat, height)` in degrees and metres. The inverse of
[`LonLatToGeocentric`](@ref); see it for why both take three coordinates.
"""
struct GeocentricToLonLat{T,K<:MathKernel} <: GeoTransformation
    a::T
    e2::T
    e4::T             # e2^2, and a2 = a^2: both appear in every evaluation
    a2::T
    one_minus_e2::T
    r2d::T
    name::Union{Nothing,Symbol}
    kernel::K
end

function _lonlat_to_geocentric(::Type{T}, a, e2, name, kernel::K) where {T,K}
    LonLatToGeocentric{T,K}(a, e2, 1 - e2, pi / 180, name, kernel)
end

function _geocentric_to_lonlat(::Type{T}, a, e2, name, kernel::K) where {T,K}
    GeocentricToLonLat{T,K}(a, e2, e2^2, a^2, 1 - e2, 180 / pi, name, kernel)
end

function LonLatToGeocentric{T}(; ellips::Ellipsoid = ellipsoid(EPSG(7030)),
                                kernel::MathKernel = DEFAULT_KERNEL) where {T}
    _lonlat_to_geocentric(T, ellips.a, ellips.e2, ellips.name, kernel)
end
LonLatToGeocentric(; kwargs...) = LonLatToGeocentric{Float64}(; kwargs...)

function GeocentricToLonLat{T}(; ellips::Ellipsoid = ellipsoid(EPSG(7030)),
                                kernel::MathKernel = DEFAULT_KERNEL) where {T}
    _geocentric_to_lonlat(T, ellips.a, ellips.e2, ellips.name, kernel)
end
GeocentricToLonLat(; kwargs...) = GeocentricToLonLat{Float64}(; kwargs...)

@inline function (t::LonLatToGeocentric{T})(lon, lat, h) where {T}
    K = t.kernel
    slat, clat = Math.sincos(K, lat * t.d2r)
    slon, clon = Math.sincos(K, lon * t.d2r)
    # Radius of curvature in the prime vertical.
    re = t.a / sqrt(1 - t.e2 * slat * slat)
    ((re + h) * clat * clon,
     (re + h) * clat * slon,
     (re * t.one_minus_e2 + h) * slat)
end

@inline function (t::GeocentricToLonLat{T})(x, y, z) where {T}
    dir, h = geocentric_direction(t, x, y, z)
    K = t.kernel
    (Math.atan(K, dir.y, dir.x) * t.r2d,
     Math.atan(K, dir.z, dir.d) * t.r2d,
     h)
end

"""
    geocentric_direction(t::GeocentricToLonLat, x, y, z) -> (Direction, height)

Vermeille's solution stopped one step short of the angles: the geodetic
[`Direction`](@ref) at the point, and the ellipsoidal height.

The height is complete, but the latitude and longitude are still the quotients
`z/d` and `y/x` rather than the arctangents of them. A projection consumes those
quotients, so a fused pipeline never forms the angles -- see
[`project_direction`](@ref). `t` itself forms them, which is the only difference
between it and this.
"""
@inline function geocentric_direction(t::GeocentricToLonLat{T}, x, y, z) where {T}
    K = t.kernel
    e2, e4 = t.e2, t.e4
    lat2 = x * x + y * y
    # Lateral and polar distances, normalized by the semi-major axis.
    p = lat2 / t.a2
    q = (t.one_minus_e2 * z * z) / t.a2
    r = (p + q - e4) / 6
    s = (e4 * p * q) / (4 * r * r * r)
    tt = Math.cbrt(K, 1 + s + sqrt(s * (2 + s)))
    u = r * (1 + tt + 1 / tt)
    rv = sqrt(u * u + e4 * q)
    w = (e2 * (u + rv - q)) / (2 * rv)
    k = sqrt(u + rv + w * w) - w
    d = (k * sqrt(lat2)) / (k + e2)
    h = ((k + e2 - 1) * sqrt(d * d + z * z)) / k
    (Direction(d, z, x, y), h)
end

# A height is what these compute, so it cannot be carried across them, and all
# three coordinates travel together.
preservesz(::LonLatToGeocentric) = false
preservesz(::GeocentricToLonLat) = false
ncoords(::LonLatToGeocentric) = 3
ncoords(::GeocentricToLonLat) = 3

# Two coordinates would mean a height of zero. That is not a sea-level default:
# it moves x and y by metres, so it is a caller error rather than an assumption
# to make silently.
(t::LonLatToGeocentric)(x, y) = throw(ArgumentError(
    "LonLatToGeocentric transforms (lon, lat, height); a two-coordinate call would mean a height " *
    "of zero, which moves x and y by metres. Pass a height."))
(t::GeocentricToLonLat)(x, y) = throw(ArgumentError(
    "GeocentricToLonLat transforms (x, y, z); a two-coordinate call would mean z = 0, a point on " *
    "the equatorial plane. Pass all three coordinates."))

Base.inv(t::LonLatToGeocentric{T}) where {T} =
    _geocentric_to_lonlat(T, t.a, t.e2, t.name, t.kernel)
Base.inv(t::GeocentricToLonLat{T}) where {T} =
    _lonlat_to_geocentric(T, t.a, t.e2, t.name, t.kernel)

islanesafe(t::LonLatToGeocentric) = vectorizes(t.kernel)
islanesafe(t::GeocentricToLonLat) = vectorizes(t.kernel)

adapt_eltype(t::LonLatToGeocentric{T}, ::Type{S}) where {T,S} =
    T === S ? t : _lonlat_to_geocentric(S, t.a, t.e2, t.name, t.kernel)
adapt_eltype(t::GeocentricToLonLat{T}, ::Type{S}) where {T,S} =
    T === S ? t : _geocentric_to_lonlat(S, t.a, t.e2, t.name, t.kernel)

Base.show(io::IO, t::LonLatToGeocentric{T}) where {T} =
    print(io, "LonLatToGeocentric{$T}(", something(t.name, "a = $(t.a)"), ")")
Base.show(io::IO, t::GeocentricToLonLat{T}) where {T} =
    print(io, "GeocentricToLonLat{$T}(", something(t.name, "a = $(t.a)"), ")")
