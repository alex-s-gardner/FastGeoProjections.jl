# Replacing a geocentric stage plus a projection with the fused pair.
#
# `Direction` and the `project_direction` protocol are in `direction.jl`, which
# has to precede the projections that implement it; this is the composition
# built from them, and it needs `GeocentricToLonLat` to exist.

"""
    FusedFromGeocentric(geo, proj)

`proj ∘ geo` with the trigonometry between them cancelled, where `geo` is a
[`GeocentricToLonLat`](@ref) and `proj` a projection that implements
[`project_direction`](@ref).

Built by [`pipeline`](@ref) in place of the composition whenever
[`fuses_direction`](@ref) holds for `proj`, and behaves identically -- the same
function to within a few ulps, at roughly half the cost. Takes three
coordinates, since the geocentric stage consumes a height; returns the
projection's two plus the geodetic height it recovered.
"""
struct FusedFromGeocentric{G<:GeocentricToLonLat,P<:GeoTransformation} <: GeoTransformation
    geo::G
    proj::P
end

@inline function (t::FusedFromGeocentric)(x, y, z)
    dir, h = geocentric_direction(t.geo, x, y, z)
    px, py = project_direction(t.proj, dir)
    (px, py, h)
end

(t::FusedFromGeocentric)(x, y) = throw(ArgumentError(
    "$(typeof(t.geo)) transforms (x, y, z); a two-coordinate call would mean z = 0, a point on " *
    "the equatorial plane. Pass all three coordinates."))

# The height is consumed by the geocentric stage, so it cannot be carried, and the
# operator takes all three coordinates.
preservesz(::FusedFromGeocentric) = false
ncoords(::FusedFromGeocentric) = 3

islanesafe(t::FusedFromGeocentric) = islanesafe(t.geo) && islanesafe(t.proj)

adapt_eltype(t::FusedFromGeocentric, ::Type{S}) where {S} =
    FusedFromGeocentric(adapt_eltype(t.geo, S), adapt_eltype(t.proj, S))

# Inverting gives back the ordinary composition: the reverse pipeline is
# `LonLatToGeocentric ∘ inv(proj)`, and fusing that is the separate question of
# whether `inv(proj)` can hand a `Direction` forward.
Base.inv(t::FusedFromGeocentric) = inv(t.geo) ∘ inv(t.proj)

Base.show(io::IO, t::FusedFromGeocentric) =
    print(io, t.proj, " ∘ ", t.geo, " [fused]")
