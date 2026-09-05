# Web Mercator, as a pair of point operators.
#
# EPSG:3857 is the projection web map tiles are drawn in. It is spherical Mercator on a
# sphere of the WGS 84 semi-major axis, applied to *geodetic* latitude -- the latitude is
# used as if it were spherical rather than being converted first. That mismatch is the
# projection's definition rather than an approximation of it, which is why the authority
# calls it Pseudo-Mercator: reading the same numbers as true ellipsoidal Mercator puts a
# point up to 20 km out in northing. Matching the authority means reproducing it exactly.
#
# So the eccentricity does not appear here at all. Both directions are closed form -- no
# series, no iteration, no conformal latitude -- which makes this the cheapest projection in
# the package.

"""
    LonLatToWebMercator(; ellips = ellipsoid(EPSG(7030)), kernel = FastKernel())
    LonLatToWebMercator{T}(; ...)

Point operator taking geodetic `(lon, lat)` in decimal degrees to web Mercator
`(x, y)` in metres. This is EPSG:4326 to EPSG:3857.

Spherical Mercator on a sphere of radius `ellips.a`, evaluated at the geodetic
latitude. See [`WebMercatorToLonLat`](@ref) for the inverse.

    julia> t = LonLatToWebMercator()
    LonLatToWebMercator{Float64}(R = 6.378137e6)

    julia> t(-45.0, 70.0)
    (-5.009377085697311e6, 1.1068715659379125e7)
"""
struct LonLatToWebMercator{T,K<:MathKernel} <: GeoTransformation
    r::T          # sphere radius, the ellipsoid's semi-major axis
    d2r::T
    kernel::K
end

"""
    WebMercatorToLonLat(; ellips = ellipsoid(EPSG(7030)), kernel = FastKernel())
    WebMercatorToLonLat{T}(; ...)

Point operator taking web Mercator `(x, y)` in metres to geodetic `(lon, lat)` in
decimal degrees. The inverse of [`LonLatToWebMercator`](@ref).

`lat = 2·atan(exp(y/R)) − π/2`, the Gudermannian in exponential form. That is
`atan(sinh(y/R))` rewritten so one `exp` does the work of a `sinh`: 3.7 ns/point
against 8.0, for 2.1e-14° where the other reaches 1.4e-14 -- both about 1e-9 m,
and both within a nanometre of Proj.

The `π/2` is subtracted in radians rather than folded into the degree conversion
as `−90`. Both forms are algebraically the same, but the folded one returns
`−3.9e-16` at `y = 0` instead of zero: `2·atan(1)` is not exactly `π/2` in
`Float64`, and scaling before subtracting scales the discrepancy up rather than
cancelling it. Subtracting first makes the equator exactly `0.0`.
"""
struct WebMercatorToLonLat{T,K<:MathKernel} <: GeoTransformation
    inv_r::T      # 1 / R, folded so the operator multiplies rather than divides
    r2d::T
    kernel::K
end

function _lonlat_to_webmerc(::Type{T}, r, kernel::K) where {T,K}
    LonLatToWebMercator{T,K}(r, pi / 180, kernel)
end

function _webmerc_to_lonlat(::Type{T}, r, kernel::K) where {T,K}
    WebMercatorToLonLat{T,K}(1 / r, 180 / pi, kernel)
end

function LonLatToWebMercator{T}(; ellips::Ellipsoid = ellipsoid(EPSG(7030)),
                                 kernel::MathKernel = DEFAULT_KERNEL) where {T}
    _lonlat_to_webmerc(T, ellips.a, kernel)
end
LonLatToWebMercator(; kwargs...) = LonLatToWebMercator{Float64}(; kwargs...)

function WebMercatorToLonLat{T}(; ellips::Ellipsoid = ellipsoid(EPSG(7030)),
                                 kernel::MathKernel = DEFAULT_KERNEL) where {T}
    _webmerc_to_lonlat(T, ellips.a, kernel)
end
WebMercatorToLonLat(; kwargs...) = WebMercatorToLonLat{Float64}(; kwargs...)

@inline function (t::LonLatToWebMercator{T})(lon, lat) where {T}
    K = t.kernel
    # `asinh(tan(φ))`, not the equivalent `log(tan(π/4 + φ/2))`. The log form is
    # 8.0 ns/point against 12.2 -- `Math.asinh` has no cheap SLEEF variant and costs
    # roughly four times `Math.log` -- but it is the worse function at both ends of
    # the range it has to cover: 2.2e-8 m against 7.5e-9 across ±85°, and `-Inf` at
    # the south pole, where its argument underflows and Proj returns a finite value.
    # A projection whose whole specification is agreement with Proj cannot take that
    # trade for 4 ns.
    (t.r * (lon * t.d2r), t.r * Math.asinh(K, Math.tan(K, lat * t.d2r)))
end

@inline function (t::WebMercatorToLonLat{T})(x, y) where {T}
    K = t.kernel
    (x * t.inv_r * t.r2d,
     (2 * Math.atan(K, Math.exp(K, y * t.inv_r)) - T(pi / 2)) * t.r2d)
end

Base.inv(t::LonLatToWebMercator{T}) where {T} = _webmerc_to_lonlat(T, t.r, t.kernel)
Base.inv(t::WebMercatorToLonLat{T}) where {T} = _lonlat_to_webmerc(T, 1 / t.inv_r, t.kernel)

islanesafe(t::LonLatToWebMercator) = vectorizes(t.kernel)
islanesafe(t::WebMercatorToLonLat) = vectorizes(t.kernel)

adapt_eltype(t::LonLatToWebMercator{T}, ::Type{S}) where {T,S} =
    T === S ? t : _lonlat_to_webmerc(S, t.r, t.kernel)
adapt_eltype(t::WebMercatorToLonLat{T}, ::Type{S}) where {T,S} =
    T === S ? t : _webmerc_to_lonlat(S, 1 / t.inv_r, t.kernel)

Base.show(io::IO, t::LonLatToWebMercator{T}) where {T} =
    print(io, "LonLatToWebMercator{$T}(R = ", t.r, ")")
Base.show(io::IO, t::WebMercatorToLonLat{T}) where {T} =
    print(io, "WebMercatorToLonLat{$T}(R = ", 1 / t.inv_r, ")")
