# Transverse Mercator, as a pair of point operators.
#
# Ported from the MATLAB implementation in geographiclib_toolbox-2.0 by Alex
# Gardner (JPL/NASA); restructured here so that each direction is a callable
# struct holding everything that depends only on the projection parameters.
#
# The series method is described in
#
#   C. F. F. Karney, Transverse Mercator with an accuracy of a few nanometers,
#   J. Geodesy 85(8), 475-485 (Aug. 2011);
#   Addenda: https://geographiclib.sourceforge.io/tm-addenda.html
#
# extending Krueger (1912) to sixth order in the flattening: errors are under
# 5 nm within 3900 km of the central meridian and under 1 mm within 7600 km,
# and the mapping continues accurately over the poles to the opposite meridian.
#
# Copyright (c) Charles Karney (2012-2022) <charles@karney.com>.

# ---------------------------------------------------------------------------
# angle helpers -- branch-free so they run on SIMD lanes
# ---------------------------------------------------------------------------

# sign, with GeographicLib's convention that zero is positive (`1 - 2*(x < 0)`).
# Base's `sign` returns 0 at 0, and VectorizationBase's returns -1; neither is
# what the projection wants.
@inline _signnz(x) = ifelse(x < zero(x), -one(x), one(x))

# reduce an angle to [-180, 180]
@inline function _angwrap(::Type{T}, x) where {T}
    y = rem(x, T(360))
    y = ifelse(y < -T(180), y + T(360), y)
    ifelse(y > T(180), y - T(360), y)
end

# reduce an angle to (-180, 180]
@inline function _angnormalize(::Type{T}, x) where {T}
    y = _angwrap(T, x)
    ifelse(y == -T(180), T(180), y)
end

"""
    _angdiff(T, lon00, lon0, x)

`x - lon0` reduced to (-180, 180], computed by two-sum so that the difference
of two nearby meridians keeps full precision. `lon00` is the precomputed
`_angwrap(-lon0)`.
"""
@inline function _angdiff(::Type{T}, lon00, lon0, x) where {T}
    b = _angwrap(T, x)

    s = lon00 + b                 # exact sum: value in s, rounding error in err
    up = s - b
    vpp = s - up
    up -= lon00
    err = b - vpp - up

    u = _angwrap(T, s)

    s2 = u + err
    up2 = s2 - err
    vpp2 = s2 - up2
    up2 -= u
    err -= (vpp2 + up2)

    # at the branch cut the sign is set by the discarded low-order part
    cut = (s2 == zero(s2)) | (abs(s2) == T(180))
    z = ifelse(err == zero(err), x - lon0, -err)
    ifelse(cut, copysign(s2, z), s2)
end

# ---------------------------------------------------------------------------
# Krueger / Karney series, evaluated once when the operator is built
# ---------------------------------------------------------------------------

# ascending-order polynomial evaluation of one GeographicLib coefficient block
@inline function _series_term(coeff, o, m, x)
    poly = zero(float(x))
    for j in 0:m
        poly = poly + coeff[o+m-j] * x^j
    end
    poly / coeff[o+m+1]
end

"""
    alpf(n)

Coefficients α_j of the geodetic-to-projected Krueger series, to sixth order.
"""
function alpf(n)
    coeff = (
        31564, -66675, 34440, 47250, -100800, 75600, 151200,
        -1983433, 863232, 748608, -1161216, 524160, 1935360,
        670412, 406647, -533952, 184464, 725760,
        6601661, -7732800, 2230245, 7257600,
        -13675556, 3438171, 7983360,
        212378941, 319334400,
    )
    _krueger_series(coeff, n)
end

"""
    betf(n)

Coefficients β_j of the projected-to-geodetic Krueger series, to sixth order.
"""
function betf(n)
    coeff = (
        384796, -382725, -6720, 932400, -1612800, 1209600, 2419200,
        -1118711, 1695744, -1174656, 258048, 80640, 3870720,
        22276, -16929, -15984, 12852, 362880,
        -830251, -158400, 197865, 7257600,
        -435388, 453717, 15966720,
        20648693, 638668800,
    )
    _krueger_series(coeff, n)
end

function _krueger_series(coeff, n)
    maxpow = 6
    out = zeros(float(typeof(n)), maxpow)
    o = 1
    d = n
    for l in 1:maxpow
        m = maxpow - l
        out[l] = d * _series_term(coeff, o, m, n)
        o += m + 2
        d *= n
    end
    NTuple{6}(out)
end

"""
    A1m1f(epsi)

`A_1 - 1`, Eq. (17) of Karney (2011): the correction from the rectifying to
the equatorial radius.
"""
function A1m1f(epsi)
    eps2 = epsi^2
    t = (0 + 64 * eps2 + 4 * eps2^2 + eps2^3) / 256
    (t + epsi) / (1 - epsi)
end

"""
    C1f(epsi)

Coefficients C_{1,k} of the meridian-distance series, used to place the
northing origin at `lat0`.
"""
function C1f(epsi)
    coeff = (
        -1, 6, -16, 32,
        -9, 64, -128, 2048,
        9, -16, 768,
        3, -5, 512,
        -7, 1280,
        -7, 2048,
    )
    nC1 = 6
    out = zeros(float(typeof(epsi)), nC1)
    eps2 = epsi^2
    o = 1
    d = epsi
    for l in 1:nC1
        m = (nC1 - l) ÷ 2
        out[l] = d * _series_term(coeff, o, m, eps2)
        o += m + 2
        d *= epsi
    end
    out
end

"""
    SinCosSeries(sinp, sinx, cosx, c)

Clenshaw summation of `sum(c[j] * sin(2j*x))` (`sinp`) or
`sum(c[j] * cos((2j-1)*x))`.
"""
function SinCosSeries(sinp::Bool, sinx, cosx, c)
    isempty(c) && return zero(sinx)
    n = length(c)
    ar = 2 * (cosx - sinx) * (cosx + sinx)
    y1 = zero(sinx)
    y0 = isodd(n) ? (n -= 1; c[n+1]) : y1
    while n > 0
        y1 = ar * y0 - y1 + c[n]
        n -= 1
        y0 = ar * y1 - y0 + c[n]
        n -= 1
    end
    sinp ? 2 * sinx * cosx * y0 : cosx * (y0 - y1)
end

"""
    norm2(x, y)

`(x, y)` scaled onto the unit circle.
"""
function norm2(x, y)
    r = hypot(x, y)
    (x / r, y / r)
end

# every parameter both directions need, derived once from (a, e, lon0, lat0)
function _tm_params(::Type{T}, a::Real, e::Real, lon0::Real, lat0::Real) where {T}
    e2 = e^2
    f = e2 / (1 + sqrt(1 - e2))
    e2m = 1 - e2
    cc = sqrt(e2m) * exp(e * atanh(e))
    n = f / (2 - f)
    b1 = (1 - f) * (A1m1f(n) + 1)
    a1 = b1 * a

    # northing of the projection origin
    y0 = if lat0 == 0
        zero(a1)
    else
        sbet0, cbet0 = norm2((1 - f) * sind(lat0), cosd(lat0))
        a1 * (atan(sbet0, cbet0) + SinCosSeries(true, sbet0, cbet0, C1f(n)))
    end

    (a = T(a), e = T(e), lon0 = T(lon0), lat0 = T(lat0),
     a1 = T(a1), b1 = T(b1), e2 = T(e2), e2m = T(e2m), cc = T(cc), n = T(n),
     y0 = T(y0), fmin = sqrt(floatmin(T)),
     alp = NTuple{6,T}(alpf(n)), bet = NTuple{6,T}(betf(n)))
end

# ---------------------------------------------------------------------------
# the conformal tangent
# ---------------------------------------------------------------------------

"""
    _taupf(K, tau, tau1, e, ehalf)

`taup = cosh(σ)·tau − sinh(σ)·hypot(1, tau)` with `σ = e·atanh(e·sinφ)`: the
tangent of the conformal latitude given the tangent of the geodetic one.

Written through `w = exp(σ) = ((1+u)/(1−u))^(e/2)` rather than as
`sinh(e·atanh(u))`. SLEEF has no fast `atanh` or `sinh`, and those two together
cost about three times the one `pow` this needs -- which matters because the
inverse evaluates this five times per point inside a Newton solve.
"""
@inline function _taupf(K, tau, tau1, e, ehalf)
    u = e * tau / tau1
    w = Math.pow(K, (one(u) + u) / (one(u) - u), ehalf)
    wi = inv(w)
    ((w + wi) * tau - (w - wi) * tau1) / 2
end

# ---------------------------------------------------------------------------
# operators
# ---------------------------------------------------------------------------

"""
    LonLatToTransverseMercator(; lon0, lat0, ellips, kernel)
    LonLatToTransverseMercator{T}(; ...)

Point operator taking geodetic `(lon, lat)` in decimal degrees to transverse
Mercator `(x, y)` in metres, with no scale factor or false origin applied
(see [`UTM`](@ref) for those).

`lon0` is the central meridian and `lat0` the latitude of the projection
origin, both in decimal degrees. `T` is the working precision.

[`convergence_scale`](@ref) returns the meridian convergence and point scale
at the same point.
"""
struct LonLatToTransverseMercator{T,K<:MathKernel} <: GeoTransformation
    a::T             # ellipsoid, kept so `inv` can rebuild the series
    e::T
    lon0::T
    lat0::T
    a1::T            # rectifying radius
    b1::T
    e2::T
    e2m::T
    cc::T            # scale at the pole
    n::T             # third flattening
    y0::T            # northing of lat0
    lon00::T         # _angwrap(-lon0), for the accurate longitude difference
    fmin::T
    alp::NTuple{6,T}
    kernel::K
end

"""
    TransverseMercatorToLonLat(; lon0, lat0, ellips, kernel)
    TransverseMercatorToLonLat{T}(; ...)

Point operator taking transverse Mercator `(x, y)` in metres back to geodetic
`(lon, lat)` in decimal degrees. See [`LonLatToTransverseMercator`](@ref) for
the parameters.
"""
struct TransverseMercatorToLonLat{T,K<:MathKernel} <: GeoTransformation
    a::T
    e::T
    lon0::T
    lat0::T
    a1::T
    b1::T
    e2::T
    e2m::T
    cc::T
    n::T
    y0::T
    lon00::T         # _angnormalize(lon0)
    bet::NTuple{6,T}
    kernel::K
end

function _lonlat_to_tm_op(::Type{T}, p, kernel::K) where {T,K}
    LonLatToTransverseMercator{T,K}(p.a, p.e, p.lon0, p.lat0, p.a1, p.b1, p.e2,
                                   p.e2m, p.cc, p.n, p.y0,
                                   _angwrap(T, -p.lon0), p.fmin, p.alp, kernel)
end

function _tm_to_lonlat_op(::Type{T}, p, kernel::K) where {T,K}
    TransverseMercatorToLonLat{T,K}(p.a, p.e, p.lon0, p.lat0, p.a1, p.b1, p.e2,
                                    p.e2m, p.cc, p.n, p.y0,
                                    _angnormalize(T, p.lon0), p.bet, kernel)
end

function LonLatToTransverseMercator{T}(; lon0::Real = 0, lat0::Real = 0,
                                       ellips::Ellipsoid = ellipsoid(EPSG(7030)),
                                       kernel::MathKernel = DEFAULT_KERNEL) where {T}
    _lonlat_to_tm_op(T, _tm_params(T, ellips.a, ellips.e, lon0, lat0), kernel)
end
LonLatToTransverseMercator(; kwargs...) = LonLatToTransverseMercator{Float64}(; kwargs...)

function TransverseMercatorToLonLat{T}(; lon0::Real = 0, lat0::Real = 0,
                                       ellips::Ellipsoid = ellipsoid(EPSG(7030)),
                                       kernel::MathKernel = DEFAULT_KERNEL) where {T}
    _tm_to_lonlat_op(T, _tm_params(T, ellips.a, ellips.e, lon0, lat0), kernel)
end
TransverseMercatorToLonLat(; kwargs...) = TransverseMercatorToLonLat{Float64}(; kwargs...)

# `inv` rebuilds the series for the other direction; a setup-time cost only
_tm_reparams(t::Union{LonLatToTransverseMercator{T},TransverseMercatorToLonLat{T}}) where {T} =
    _tm_params(T, t.a, t.e, t.lon0, t.lat0)

Base.inv(t::LonLatToTransverseMercator{T}) where {T} =
    _tm_to_lonlat_op(T, _tm_reparams(t), t.kernel)
Base.inv(t::TransverseMercatorToLonLat{T}) where {T} =
    _lonlat_to_tm_op(T, _tm_reparams(t), t.kernel)

islanesafe(t::LonLatToTransverseMercator) = vectorizes(t.kernel)
islanesafe(t::TransverseMercatorToLonLat) = vectorizes(t.kernel)

adapt_eltype(t::LonLatToTransverseMercator{T,K}, ::Type{S}) where {T,K,S} =
    T === S ? t :
    LonLatToTransverseMercator{S,K}(S(t.a), S(t.e), S(t.lon0), S(t.lat0), S(t.a1),
                                   S(t.b1), S(t.e2), S(t.e2m), S(t.cc), S(t.n),
                                   S(t.y0), S(t.lon00), sqrt(floatmin(S)),
                                   NTuple{6,S}(t.alp), t.kernel)
adapt_eltype(t::TransverseMercatorToLonLat{T,K}, ::Type{S}) where {T,K,S} =
    T === S ? t :
    TransverseMercatorToLonLat{S,K}(S(t.a), S(t.e), S(t.lon0), S(t.lat0), S(t.a1),
                                   S(t.b1), S(t.e2), S(t.e2m), S(t.cc), S(t.n),
                                   S(t.y0), S(t.lon00),
                                   NTuple{6,S}(t.bet), t.kernel)

Base.show(io::IO, t::LonLatToTransverseMercator{T}) where {T} =
    print(io, "LonLatToTransverseMercator{$T}(lon0 = ", t.lon0, "°, lat0 = ", t.lat0, "°)")
Base.show(io::IO, t::TransverseMercatorToLonLat{T}) where {T} =
    print(io, "TransverseMercatorToLonLat{$T}(lon0 = ", t.lon0, "°, lat0 = ", t.lat0, "°)")

# ---------------------------------------------------------------------------
# the kernels
#
# `FULL` also returns the meridian convergence and point scale. It is a
# compile-time parameter, so the plain `(x, y)` call drops the second Clenshaw
# recurrence (the derivative series) entirely rather than computing and
# discarding it.
# ---------------------------------------------------------------------------

@inline function _lonlat_to_tm(t::LonLatToTransverseMercator{T, KT}, lon_in, lat, ::Val{FULL}) where {T, KT, FULL}
    K = t.kernel
    d2r = T(pi / 180)

    lon = _angdiff(T, t.lon00, t.lon0, lon_in)

    latsign = _signnz(lat)
    lonsign = _signnz(lon)
    lon = lon * lonsign
    lat = lat * latsign

    # the far side of the central meridian is the near side reflected
    backside = lon > T(90)
    latsign = ifelse(backside & (lat == zero(lat)), -one(latsign), latsign)
    lon = ifelse(backside, T(180) - lon, lon)

    slam, clam = Math.sincos(K, lon * d2r)
    slat, clat = Math.sincos(K, lat * d2r)

    tau = slat / max(t.fmin, clat)
    tau1 = sqrt(tau * tau + one(tau))
    taup = _taupf(K, tau, tau1, t.e, t.e / 2)
    htc = sqrt(taup * taup + clam * clam)

    atpole = lat == T(90)
    xip = ifelse(atpole, T(pi) / 2, Math.atan(K, taup, clam))
    etap = ifelse(atpole, zero(taup), Math.asinh(K, slam / htc))

    # Clenshaw summation of the alpha series, in the complex variable
    # 2*(xip + i*etap)
    s0, c0 = Math.sincos(K, 2 * xip)
    sh0, ch0 = Math.sinhcosh(K, 2 * etap)

    ar = 2 * c0 * ch0
    ai = -2 * s0 * sh0

    a1_, a2_, a3_, a4_, a5_, a6_ = t.alp

    y1r_6 = a6_
    y0r_6 = ar * y1r_6 + a5_
    y0i_6 = ai * y1r_6

    y1r_4 = ar * y0r_6 - ai * y0i_6 - y1r_6 + a4_
    y1i_4 = ar * y0i_6 + ai * y0r_6

    y0r_4 = ar * y1r_4 - ai * y1i_4 - y0r_6 + a3_
    y0i_4 = ar * y1i_4 + ai * y1r_4 - y0i_6

    y1r_2 = ar * y0r_4 - ai * y0i_4 - y1r_4 + a2_
    y1i_2 = ar * y0i_4 + ai * y0r_4 - y1i_4

    y0r = ar * y1r_2 - ai * y1i_2 - y0r_4 + a1_
    y0i = ar * y1i_2 + ai * y1r_2 - y0i_4

    xi = xip + y0r * s0 * ch0 - y0i * c0 * sh0
    eta = etap + y0r * c0 * sh0 + y0i * s0 * ch0

    xi = ifelse(backside, T(pi) - xi, xi)

    xout = t.a1 * eta * lonsign
    yout = t.a1 * xi * latsign - t.y0

    if FULL
        # the convergence is set by the *conformal* tangent, the scale by the
        # geodetic one
        taup1 = sqrt(one(taup) + taup * taup)
        gam = ifelse(atpole, lon, Math.atan(K, slam * taup, clam * taup1) / d2r)
        k = ifelse(atpole, t.cc, sqrt(t.e2m + t.e2 * clat * clat) * tau1 / htc)

        # the same recurrence differentiated: gives d(xi + i*eta)/d(xip + i*etap)
        z1r_6 = 12 * a6_
        z0r_6 = ar * z1r_6 + 10 * a5_
        z0i_6 = ai * z1r_6

        z1r_4 = ar * z0r_6 - ai * z0i_6 - z1r_6 + 8 * a4_
        z1i_4 = ar * z0i_6 + ai * z0r_6

        z0r_4 = ar * z1r_4 - ai * z1i_4 - z0r_6 + 6 * a3_
        z0i_4 = ar * z1i_4 + ai * z1r_4 - z0i_6

        z1r_2 = ar * z0r_4 - ai * z0i_4 - z1r_4 + 4 * a2_
        z1i_2 = ar * z0i_4 + ai * z0r_4 - z1i_4

        z0r_2 = ar * z1r_2 - ai * z1i_2 - z0r_4 + 2 * a1_
        z0i_2 = ar * z1i_2 + ai * z1r_2 - z0i_4

        zr = one(ar) - z1r_2 + z0r_2 * ar / 2 - z0i_2 * ai / 2
        zi = -z1i_2 + z0r_2 * ai / 2 + z0i_2 * ar / 2

        gam = gam - Math.atan(K, zi, zr) / d2r
        k = k * (t.b1 * sqrt(zr * zr + zi * zi))

        gam = ifelse(backside, T(180) - gam, gam)
        gam = _angnormalize(T, gam * latsign * lonsign)
        return (xout, yout, gam, k)
    end

    (xout, yout)
end

# one Newton step recovering tau from taup (Karney eq. 20)
@inline function _tau_step(K, tau, taup_target, e, ehalf, e2m)
    tau1 = sqrt(tau * tau + one(tau))
    taupa = _taupf(K, tau, tau1, e, ehalf)
    tau + (taup_target - taupa) * (one(tau) + e2m * tau * tau) /
          (e2m * tau1 * sqrt(taupa * taupa + one(taupa)))
end

@inline function _tm_to_lonlat(t::TransverseMercatorToLonLat{T}, x, y, ::Val{FULL}) where {T,FULL}
    K = t.kernel
    d2r = T(pi / 180)

    xiA = (y + t.y0) / t.a1
    etaA = x / t.a1
    xisign = _signnz(xiA)
    etasign = _signnz(etaA)
    xiB = xiA * xisign
    eta = etaA * etasign

    backside = xiB > T(pi) / 2
    xi = ifelse(backside, T(pi) - xiB, xiB)

    s0, c0 = Math.sincos(K, 2 * xi)
    sh0, ch0 = Math.sinhcosh(K, 2 * eta)

    ar = 2 * c0 * ch0
    ai = -2 * s0 * sh0

    b1_, b2_, b3_, b4_, b5_, b6_ = t.bet

    y1r_6 = -b6_
    y0r_6 = ar * y1r_6 - b5_
    y0i_6 = ai * y1r_6

    y1r_4 = ar * y0r_6 - ai * y0i_6 - y1r_6 - b4_
    y1i_4 = ar * y0i_6 + ai * y0r_6

    y0r_4 = ar * y1r_4 - ai * y1i_4 - y0r_6 - b3_
    y0i_4 = ar * y1i_4 + ai * y1r_4 - y0i_6

    y1r_2 = ar * y0r_4 - ai * y0i_4 - y1r_4 - b2_
    y1i_2 = ar * y0i_4 + ai * y0r_4 - y1i_4

    y0r_2 = ar * y1r_2 - ai * y1i_2 - y0r_4 - b1_
    y0i_2 = ar * y1i_2 + ai * y1r_2 - y0i_4

    xip = xi + y0r_2 * s0 * ch0 - y0i_2 * c0 * sh0
    etap = eta + y0r_2 * c0 * sh0 + y0i_2 * s0 * ch0

    # `s` carries the longitude directly, so it needs the accurate `sinh`: the
    # cancelling exp form would lose its leading digits near the central meridian
    s = Math.sinh(K, etap)
    sxip, cxip = Math.sincos(K, xip)
    c = max(zero(xip), cxip)
    r = sqrt(s * s + c * c)
    lon = Math.atan(K, s, c) / d2r
    sr = sxip / r

    # Newton on tau. The initial guess taup/e2m is already good to O(e²) and
    # Newton doubles the digits each step, so for an Earth-like ellipsoid the
    # last ulp is reached after two; Karney's own solver stops there too. The
    # old array kernel ran five unconditionally because a `@turbo` loop cannot
    # break out early.
    eh = t.e / 2
    tau = sr / t.e2m
    tau = _tau_step(K, tau, sr, t.e, eh, t.e2m)
    tau = _tau_step(K, tau, sr, t.e, eh, t.e2m)
    tau = _tau_step(K, tau, sr, t.e, eh, t.e2m)

    lat = Math.atan(K, tau) / d2r

    # r == 0 is the pole, where the longitude is indeterminate
    onsphere = r != zero(r)
    lat = ifelse(onsphere, lat, T(90))
    lon = ifelse(onsphere, lon, zero(lon))

    lat = lat * xisign
    lon = ifelse(backside, T(180) - lon, lon)
    lon = _angnormalize(T, lon * etasign + t.lon00)

    if FULL
        z1r_6 = -12 * b6_
        z0r_6 = ar * z1r_6 - 10 * b5_
        z0i_6 = ai * z1r_6

        z1r_4 = ar * z0r_6 - ai * z0i_6 - z1r_6 - 8 * b4_
        z1i_4 = ar * z0i_6 + ai * z0r_6

        z0r_4 = ar * z1r_4 - ai * z1i_4 - z0r_6 - 6 * b3_
        z0i_4 = ar * z1i_4 + ai * z1r_4 - z0i_6

        z1r_2 = ar * z0r_4 - ai * z0i_4 - z1r_4 - 4 * b2_
        z1i_2 = ar * z0i_4 + ai * z0r_4 - z1i_4

        z0r_2 = ar * z1r_2 - ai * z1i_2 - z0r_4 - 2 * b1_
        z0i_2 = ar * z1i_2 + ai * z1r_2 - z0i_4

        zr = one(ar) - z1r_2 + z0r_2 * ar / 2 - z0i_2 * ai / 2
        zi = -z1i_2 + z0r_2 * ai / 2 + z0i_2 * ar / 2

        gam = Math.atan(K, zi, zr) / d2r
        k = t.b1 / sqrt(zr * zr + zi * zi)
        gam = gam + Math.atan(K, sxip * Math.tanh(K, etap), c) / d2r

        st = sqrt(one(tau) + tau * tau)
        k = k * ifelse(onsphere, sqrt(t.e2m + t.e2 / (st * st)) * st * r, one(r))
        k = ifelse(onsphere, k, k * t.cc)

        gam = ifelse(backside, T(180) - gam, gam)
        gam = _angnormalize(T, gam * xisign * etasign)
        return (lon, lat, gam, k)
    end

    (lon, lat)
end

@inline (t::LonLatToTransverseMercator)(lon, lat) = _lonlat_to_tm(t, lon, lat, Val(false))
@inline (t::TransverseMercatorToLonLat)(x, y) = _tm_to_lonlat(t, x, y, Val(false))

"""
    convergence_scale(t, x, y) -> (γ, k)
    convergence_scale(t, point) -> (γ, k)

Meridian convergence `γ` in decimal degrees and point scale factor `k` of the
transverse Mercator projection `t`, at `(lon, lat)` for a `LonLatTo…` operator
and at `(x, y)` for a `…ToLonLat` one. The point may also be given as a single
GeoInterface point, exactly as for a transformation call.

These come out of the same series the projection itself evaluates, so asking
for them costs roughly twice a plain projection call rather than a second one.
"""
@inline convergence_scale(t::LonLatToTransverseMercator, lon, lat) =
    (r = _lonlat_to_tm(t, lon, lat, Val(true)); (r[3], r[4]))
@inline convergence_scale(t::TransverseMercatorToLonLat, x, y) =
    (r = _tm_to_lonlat(t, x, y, Val(true)); (r[3], r[4]))

# same GeoInterface point handling as a transformation call
@inline convergence_scale(t, p) = _apply_convergence_scale(t, GI.geomtrait(p), p)
@inline _apply_convergence_scale(t, ::GI.PointTrait, p) = convergence_scale(t, GI.x(p), GI.y(p))
