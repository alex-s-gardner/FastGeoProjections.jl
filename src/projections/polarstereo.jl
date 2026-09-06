# Polar stereographic projection, as a pair of point operators.
#
# Julia port by Alex Gardner (JPL/NASA, 2023) of a MATLAB implementation by
# Andy Bliss (2011); see Snyder for the derivation. Restructured here so that
# each direction is a callable struct holding everything that depends only on
# the projection parameters.

"""
    LonLatToPolarStereographic(; lon_0, lat_ts, ellips, kernel)
    LonLatToPolarStereographic{T}(; ...)

Point operator taking geodetic `(lon, lat)` in decimal degrees to polar
stereographic `(x, y)` in metres.

`lon_0` is the meridian along the positive `y` axis and `lat_ts` the standard
parallel (the latitude of true scale); a negative `lat_ts` selects the southern
hemisphere. `T` is the working precision, `Float64` by default.
"""
struct LonLatToPolarStereographic{T,K<:MathKernel} <: GeoTransformation
    a::T          # semi-major axis
    e::T          # eccentricity
    lon_0::T      # central meridian [rad], sign-flipped in the south
    t_c::T        # Snyder's t and m at the standard parallel
    m_c::T
    pm::T         # +1 north, -1 south
    d2r::T
    s_lon_0::T    # sine and cosine of the central meridian, for the angle
    c_lon_0::T    # subtraction a fused pipeline does in place of `sincos`
    kernel::K
end

"""
    PolarStereographicToLonLat(; lon_0, lat_ts, ellips, kernel)
    PolarStereographicToLonLat{T}(; ...)

Point operator taking polar stereographic `(x, y)` in metres to geodetic
`(lon, lat)` in decimal degrees. See [`LonLatToPolarStereographic`](@ref) for the
parameters.

Latitude comes from the conformal latitude by the sine series `φ = χ + Σ cₖ
sin(2kχ)`, carried to `sin(10χ)` with every coefficient complete through `e¹⁰`
and summed by [`Math.sin2_series`](@ref), which costs one `sincos` for the whole
series.
The series is truncated in two independent ways -- the number of harmonics and the
`e²` order of each coefficient -- and the second is the one that binds: at four
harmonics with coefficients through `e⁸` every term carries a comparable error, so
dropping `sin(10χ)` alone would cost little while stopping the coefficients at `e⁸`
costs two orders of magnitude. Both together reach 2.4e-8 m against Proj, at every
ellipsoid in [`ellipsoid`](@ref); the term beyond is what now limits it.

A Newton solve on the exact relation, as [`TransverseMercatorToLonLat`](@ref)
runs, reaches 1.5e-9 m instead. It is not used here: it costs about five times
this series, and it forms `tan(χ)` as `(1 - t²)/(2t)`, which is `NaN` at `t = 0`
-- the pole, where polar stereographic data actually sits. The series has no such
point.
"""
struct PolarStereographicToLonLat{T,K<:MathKernel} <: GeoTransformation
    a::T
    e::T
    lon_0::T
    t_c::T
    m_c::T
    pm::T
    t_c_over_am_c::T   # t_c / (a * m_c), folded into one multiply
    r2d::T
    c::NTuple{5,T}     # conformal -> geodetic latitude series, through e^10
    kernel::K
end

# the hemisphere flip and Snyder's t_c / m_c, shared by both directions
function _ps_constants(lon_0::Real, lat_ts::Real, ellips::Ellipsoid)
    lat_ts = lat_ts * pi / 180
    lon_0 = lon_0 * pi / 180
    if lat_ts < 0                     # standard parallel in the south: flip
        pm = -1.0
        lat_ts = -lat_ts
        lon_0 = -lon_0
    else
        pm = 1.0
    end
    e = ellips.e
    t_c = tan(pi / 4 - lat_ts / 2) /
          ((1 - e * sin(lat_ts)) / (1 + e * sin(lat_ts)))^(e / 2)
    m_c = cos(lat_ts) / sqrt(1 - e^2 * sin(lat_ts)^2)
    (a = ellips.a, e = e, lon_0 = lon_0, t_c = t_c, m_c = m_c, pm = pm)
end

# build either direction from the already-reduced parameter set, so that `inv`
# never has to recover degrees from radians
function _lonlat_to_ps(::Type{T}, a, e, lon_0, t_c, m_c, pm, kernel::K) where {T,K}
    s0, c0 = sincos(lon_0)
    LonLatToPolarStereographic{T,K}(a, e, lon_0, t_c, m_c, pm, pi / 180, s0, c0, kernel)
end

function _ps_to_lonlat(::Type{T}, a, e, lon_0, t_c, m_c, pm, kernel::K) where {T,K}
    PolarStereographicToLonLat{T,K}(
        a, e, lon_0, t_c, m_c, pm, t_c / (a * m_c), 180 / pi,
        (e^2 / 2 + 5 * e^4 / 24 + e^6 / 12 + 13 * e^8 / 360 + 3 * e^10 / 160,
         7 * e^4 / 48 + 29 * e^6 / 240 + 811 * e^8 / 11520 + 81 * e^10 / 2240,
         7 * e^6 / 120 + 81 * e^8 / 1120 + 3029 * e^10 / 53760,
         4279 * e^8 / 161280 + 883 * e^10 / 20160,
         2087 * e^10 / 161280),
        kernel)
end

function LonLatToPolarStereographic{T}(; lon_0::Real, lat_ts::Real,
                                        ellips::Ellipsoid = ellipsoid(EPSG(7030)),
                                        kernel::MathKernel = DEFAULT_KERNEL) where {T}
    c = _ps_constants(lon_0, lat_ts, ellips)
    _lonlat_to_ps(T, c.a, c.e, c.lon_0, c.t_c, c.m_c, c.pm, kernel)
end
LonLatToPolarStereographic(; kwargs...) = LonLatToPolarStereographic{Float64}(; kwargs...)

function PolarStereographicToLonLat{T}(; lon_0::Real, lat_ts::Real,
                                        ellips::Ellipsoid = ellipsoid(EPSG(7030)),
                                        kernel::MathKernel = DEFAULT_KERNEL) where {T}
    c = _ps_constants(lon_0, lat_ts, ellips)
    _ps_to_lonlat(T, c.a, c.e, c.lon_0, c.t_c, c.m_c, c.pm, kernel)
end
PolarStereographicToLonLat(; kwargs...) = PolarStereographicToLonLat{Float64}(; kwargs...)

@inline function (p::LonLatToPolarStereographic{T})(lon, lat) where {T}
    K = p.kernel
    latr = lat * p.d2r * p.pm
    s = Math.sin(K, latr)
    t = Math.tan(K, T(pi / 4) - latr / 2) /
        Math.conformal_ratio(K, p.e, p.e * s)
    rho = p.a * p.m_c * t / p.t_c            # true scale at lat_ts
    dl = lon * p.d2r * p.pm - p.lon_0
    sdl, cdl = Math.sincos(K, dl)
    (p.pm * rho * sdl, -p.pm * rho * cdl)
end

@inline function (p::PolarStereographicToLonLat{T})(x, y) where {T}
    K = p.kernel
    rho = sqrt(x * x + y * y)
    chi = T(pi / 2) - 2 * Math.atan(K, rho * p.t_c_over_am_c)   # conformal latitude

    lat = chi + Math.sin2_series(K, chi, p.c)
    lon = p.lon_0 + Math.atan(K, x, -p.pm * y)   # pm folds the southern flip of y

    lat = lat * p.pm * p.r2d
    lon = p.pm * (mod(p.pm * lon + T(pi), T(2pi)) - T(pi)) * p.r2d
    (lon, lat)
end

Base.inv(p::LonLatToPolarStereographic{T}) where {T} =
    _ps_to_lonlat(T, p.a, p.e, p.lon_0, p.t_c, p.m_c, p.pm, p.kernel)
Base.inv(p::PolarStereographicToLonLat{T}) where {T} =
    _lonlat_to_ps(T, p.a, p.e, p.lon_0, p.t_c, p.m_c, p.pm, p.kernel)

# Everything this projection wants from the geodetic angles is their sines and
# cosines, so a point arriving as a `Direction` skips both `atan`s that formed
# them and both `sincos`es that took them apart. `tan(pi/4 - lat/2)` is the
# half-angle identity `cos/(1 + sin)`, which in Cartesian terms is `d/(R + z)`.
fuses_direction(::LonLatToPolarStereographic) = true

@inline function project_direction(p::LonLatToPolarStereographic{T},
                                   dir::Direction) where {T}
    K = p.kernel
    # `pm` flips the hemisphere, which in this form is a flip of z's sign.
    zp = dir.z * p.pm
    R = sqrt(dir.d * dir.d + zp * zp)
    s = zp / R
    t = (dir.d / (R + zp)) /
        Math.conformal_ratio(K, p.e, p.e * s)
    rho = p.a * p.m_c * t / p.t_c
    # sin and cos of (lon * pm - lon_0), by angle subtraction on the direction's
    # own sine and cosine against the central meridian's precomputed pair.
    slon, clon = sincos_lon(dir)
    slon = slon * p.pm
    sdl = slon * p.c_lon_0 - clon * p.s_lon_0
    cdl = clon * p.c_lon_0 + slon * p.s_lon_0
    (p.pm * rho * sdl, -p.pm * rho * cdl)
end

islanesafe(p::LonLatToPolarStereographic) = vectorizes(p.kernel)
islanesafe(p::PolarStereographicToLonLat) = vectorizes(p.kernel)

adapt_eltype(p::LonLatToPolarStereographic{T}, ::Type{S}) where {T,S} =
    T === S ? p : _lonlat_to_ps(S, p.a, p.e, p.lon_0, p.t_c, p.m_c, p.pm, p.kernel)
adapt_eltype(p::PolarStereographicToLonLat{T}, ::Type{S}) where {T,S} =
    T === S ? p : _ps_to_lonlat(S, p.a, p.e, p.lon_0, p.t_c, p.m_c, p.pm, p.kernel)

Base.show(io::IO, p::LonLatToPolarStereographic{T}) where {T} =
    print(io, "LonLatToPolarStereographic{$T}(lon_0 = ", p.pm * p.lon_0 * 180 / pi, "°)")
Base.show(io::IO, p::PolarStereographicToLonLat{T}) where {T} =
    print(io, "PolarStereographicToLonLat{$T}(lon_0 = ", p.pm * p.lon_0 * 180 / pi, "°)")
