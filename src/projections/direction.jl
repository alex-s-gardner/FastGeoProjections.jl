# The geodetic direction, and the protocol a projection implements to accept one.
#
# Composing `GeocentricToLonLat` with a projection pays for trigonometry twice
# over: the geocentric inverse forms `lat = atan(z, d)` and `lon = atan(y, x)`,
# and the projection immediately takes the sine and cosine of both. Neither
# needs the angle, only its sine and cosine, and those are available directly
# from the Cartesian quantities:
#
#     sin(lat) = z / hypot(d, z)      cos(lat) = d / hypot(d, z)
#     sin(lon) = y / hypot(x, y)      cos(lon) = x / hypot(x, y)
#     tan(lat) = z / d
#
# so the pair of `atan`s and the pair of `sincos`es cancel. What passes between
# the stages instead is `Direction`: the unnormalized geodetic direction, which
# is what `GeocentricToLonLat` has before it forms an angle and what a
# projection wants after it has taken one apart.
#
# A projection opts in by implementing `project_direction`. The default forms
# the angles and calls the projection, so a projection that does not implement
# it still composes correctly -- just without the saving.

"""
    Direction(d, z, x, y)

The geodetic direction at a point, unnormalized: `atan(z, d)` is its latitude
and `atan(y, x)` its longitude, but the quotients that a projection actually
consumes are available without forming either angle.

`d` is the distance from the axis to the point's normal intercept, which is what
the geocentric inverse computes on the way to a latitude; `x` and `y` are the
equatorial Cartesian coordinates. The four are not independent -- `hypot(x, y)`
and `d` differ only by the ellipsoid's curvature correction -- but both are
needed, `d` for the latitude and `(x, y)` for the longitude.

The type parameter follows the components, which are whatever the arithmetic that
produced them returned -- a scalar on the scalar path, a `Vec` of lanes on the
SIMD one. It is not the operator's own precision: pinning it to that would make
every `Direction` a scalar and keep the geocentric conversions off the lane loop.
"""
struct Direction{T}
    d::T
    z::T
    x::T
    y::T
end

"""
    sincos_lat(dir)

`(sin, cos)` of the geodetic latitude, without forming the angle.
"""
@inline function sincos_lat(dir::Direction)
    r = sqrt(dir.d * dir.d + dir.z * dir.z)
    (dir.z / r, dir.d / r)
end

"""
    sincos_lon(dir)

`(sin, cos)` of the geodetic longitude, without forming the angle.
"""
@inline function sincos_lon(dir::Direction)
    r = sqrt(dir.x * dir.x + dir.y * dir.y)
    (dir.y / r, dir.x / r)
end

"""
    lat_degrees(dir)
    lon_degrees(dir)

The geodetic angles in degrees. What the default [`project_direction`](@ref)
falls back to, and what a projection reaches for where it needs the angle itself
rather than its sine and cosine -- the reduction against a central meridian in
[`LonLatToTransverseMercator`](@ref), for one.
"""
#
# `oftype`, not `T(...)`: on the lane path `T` is a `Vec`, which cannot be
# constructed from a scalar constant.
@inline lat_degrees(dir::Direction) = atan(dir.z, dir.d) * oftype(dir.z, 180 / pi)
@inline lon_degrees(dir::Direction) = atan(dir.y, dir.x) * oftype(dir.y, 180 / pi)

"""
    project_direction(t, dir::Direction)

Apply the projection `t` to a point given as a [`Direction`](@ref) rather than as
`(lon, lat)`, so that the trigonometry the geocentric stage before it would have
undone is never performed.

The default forms both angles and calls `t`, which is correct for any
projection; implementing it is an optimization, and one worth making only for a
projection whose own first act is to take the sine and cosine of what it was
given. [`fuses_direction`](@ref) reports whether a projection has.
"""
@inline project_direction(t::GeoTransformation, dir::Direction) =
    t(lon_degrees(dir), lat_degrees(dir))

"""
    fuses_direction(t)

Whether `t` implements [`project_direction`](@ref) natively, so that composing a
geocentric stage onto it is worth replacing with [`FusedFromGeocentric`](@ref).

False by default: a projection that has not opted in still composes correctly
through the fallback, and fusing it would only add a layer.
"""
fuses_direction(::GeoTransformation) = false
