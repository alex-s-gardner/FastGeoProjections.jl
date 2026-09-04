"""
    GeoTransformation

Supertype of every coordinate transformation in FastGeoProjections.

A `GeoTransformation` is an *atomic point operator*: calling it on a single
point returns a single point,

    t(x, y)  ->  (x′, y′)

with everything the projection can precompute stored in the struct. Nothing
about arrays, threading, or vectorization is baked into the transformation
itself -- that is the job of [`transform`](@ref) / [`transform!`](@ref), which
apply the same operator to whole collections.

Because the operator is generic over its scalar type, it evaluates unchanged on
SIMD lanes, which is how the point-by-point form matches array-kernel
throughput.

A single argument is taken to be a point in the GeoInterface sense, so anything
with `GeoInterface.PointTrait` works -- including the plain `(x, y)` tuples the
array API uses:

    t((x, y))
    t(GI.Point(x, y))

Only the first two coordinates are used; any `z` or `m` is ignored.

Transformations compose with `∘` and reverse with `inv`.
"""
abstract type GeoTransformation <: CoordinateTransformations.Transformation end

# --- calling conventions ---------------------------------------------------
# `t(x, y)` is the primitive; each projection defines that method and nothing
# else. A single argument goes through GeoInterface, whose `PointTrait` already
# covers `NTuple{2}` (and tuples of `Vec`), so this costs nothing at run time.
@inline (t::GeoTransformation)(p) = _apply_point(t, GI.geomtrait(p), p)
@inline _apply_point(t::GeoTransformation, ::GI.PointTrait, p) = t(GI.x(p), GI.y(p))
@inline _apply_point(t::GeoTransformation, trait, p) =
    throw(ArgumentError("a $(typeof(t)) transforms points; got $(trait === nothing ? typeof(p) : trait)"))

"""
    islanesafe(t)

Whether `t` may be evaluated on SIMD lanes (`VectorizationBase.Vec`) rather
than one scalar point at a time. True for transformations built only from
branch-free arithmetic and a vectorizing [`MathKernel`](@ref); false for
anything that calls into Proj or Base's libm.
"""
islanesafe(::GeoTransformation) = false

"""
    adapt_eltype(t, ::Type{T})

Rebuild `t` with its precomputed parameters stored as `T`, so a transformation
can be applied to `Float32` data at `Float32` cost. Called once per
`transform`, never per point.
"""
adapt_eltype(t::GeoTransformation, ::Type{T}) where {T} = t

# --- identity --------------------------------------------------------------
"""
    Identity()

The transformation that returns its argument unchanged.
"""
struct Identity <: GeoTransformation end

@inline (::Identity)(x, y) = (x, y)
Base.inv(t::Identity) = t
islanesafe(::Identity) = true

# --- axis order ------------------------------------------------------------
"""
    SwapXY()

Swap the two components of a point. Used to bridge authority axis order
(latitude first, for a geographic CRS) and visualization order (`always_xy`).
"""
struct SwapXY <: GeoTransformation end

@inline (::SwapXY)(x, y) = (y, x)
Base.inv(t::SwapXY) = t
islanesafe(::SwapXY) = true

# --- composition -----------------------------------------------------------
"""
    ComposedGeoTransformation(transformations::Tuple)

A pipeline of transformations, applied in tuple order: the first element sees
the input point, the last produces the result.

Composition is flat. `f ∘ g ∘ h` is one `ComposedGeoTransformation` holding
`(h, g, f)`, not a tree of nested pairs, so a whole EPSG-to-EPSG pipeline is a
single concrete callable that inlines into one pass over the data however many
stages it has.
"""
struct ComposedGeoTransformation{T<:Tuple} <: GeoTransformation
    transformations::T
end

ComposedGeoTransformation(ts::GeoTransformation...) = ComposedGeoTransformation(ts)

@inline (c::ComposedGeoTransformation)(x, y) = _applychain(c.transformations, x, y)

@inline _applychain(::Tuple{}, x, y) = (x, y)
@inline function _applychain(ts::Tuple, x, y)
    p = first(ts)(x, y)
    _applychain(Base.tail(ts), p[1], p[2])
end

# `∘` keeps Base's meaning -- `outer ∘ inner` applies `inner` first -- so the
# tuple is built inner-first, and nesting flattens
_chain(t::GeoTransformation) = (t,)
_chain(c::ComposedGeoTransformation) = c.transformations
_chain(::Identity) = ()

_compose(ts::Tuple{}) = Identity()
_compose(ts::Tuple{Any}) = only(ts)
_compose(ts::Tuple) = ComposedGeoTransformation(ts)

Base.:∘(outer::GeoTransformation, inner::GeoTransformation) =
    _compose((_chain(inner)..., _chain(outer)...))

Base.inv(c::ComposedGeoTransformation) =
    _compose(map(inv, reverse(c.transformations)))

islanesafe(c::ComposedGeoTransformation) = all(islanesafe, c.transformations)
adapt_eltype(c::ComposedGeoTransformation, ::Type{T}) where {T} =
    ComposedGeoTransformation(map(t -> adapt_eltype(t, T), c.transformations))

function Base.show(io::IO, c::ComposedGeoTransformation)
    join(io, reverse(c.transformations), " ∘ ")
end
