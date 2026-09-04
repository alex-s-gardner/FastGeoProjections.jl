"""
    Transformation(source_epsg, target_epsg; always_xy=false, threaded=true, proj_only=false, T=Float64, kernel=FastKernel())

A transformation pipeline between two coordinate reference systems.
`Transformation` implements the
[CoordinateTransformations.jl](https://github.com/JuliaGeometry/CoordinateTransformations.jl)
API: call an instance like a function.

It is an *atomic point operator* -- `trans((x, y))` transforms one point -- and
that same operator is what [`transform`](@ref) and [`transform!`](@ref) apply
to whole collections, on SIMD lanes and across threads. Vector arguments are
accepted directly for convenience:

    trans(x, y)         # two numbers  -> (x′, y′)
    trans((x, y))       # one point    -> (x′, y′)
    trans(X, Y)         # two vectors  -> (X′, Y′)
    trans(points)       # vector of points -> vector of points

`source_epsg` and `target_epsg` are EPSG authority codes (see
<https://epsg.io/>), given as `EPSG(3413)` or `"EPSG:3413"`.

`always_xy` fixes the axis order to x,y (lon,lat). By default it is `false`,
meaning the order is the one defined by the authority in charge of the CRS, as
explained in [this PROJ FAQ entry](https://proj.org/faq.html#why-is-the-axis-ordering-in-proj-not-consistent)
-- so EPSG:4326 is latitude first.

`threaded` sets the default for whole-array calls; it can be overridden per
call with the `threaded` keyword of `transform`/`transform!`.

`proj_only` forces the use of Proj.jl even where a native FastGeoProjections
implementation exists. By default Proj.jl is used only when there is none.

`T` is the working precision and `kernel` the transcendental back-end; see
[`FastKernel`](@ref), [`SLEEFKernel`](@ref) and [`BaseKernel`](@ref).

# Examples
```julia
julia> trans = FastGeoProjections.Transformation(EPSG(4326), EPSG(3413); always_xy=true)
Transformation
    source_epsg:    EPSG:4326
    target_epsg:    EPSG:3413
    threaded:       true
    always_xy:      true
    proj_only:      false

julia> trans(-45.0, 70.0)
(0.0, -2.187927649279021e6)
```
"""
struct Transformation{F<:GeoTransformation} <: GeoTransformation
    f::F
    source_epsg::EPSG
    target_epsg::EPSG
    always_xy::Bool
    threaded::Bool
    proj_only::Bool
end

function Transformation(source_epsg::EPSG, target_epsg::EPSG;
                        threaded::Bool = true,
                        always_xy::Bool = false,
                        proj_only::Bool = false,
                        T::Type = Float64,
                        kernel::MathKernel = DEFAULT_KERNEL)
    f = pipeline(source_epsg, target_epsg; always_xy, proj_only, T, kernel)
    Transformation(f, source_epsg, target_epsg, always_xy, threaded, proj_only)
end

Transformation(source_epsg::String, target_epsg::String; kwargs...) =
    Transformation(EPSG(source_epsg), EPSG(target_epsg); kwargs...)

# the pipeline is transparent: a Transformation behaves exactly as the operator
# it wraps
@inline (t::Transformation)(x, y) = t.f(x, y)
@inline (t::Transformation)(x, y, z) = t.f(x, y, z)
islanesafe(t::Transformation) = islanesafe(t.f)
preservesz(t::Transformation) = preservesz(t.f)
adapt_eltype(t::Transformation, ::Type{T}) where {T} = adapt_eltype(t.f, T)
_transform_pts!(dest, t::Transformation, src, threaded) =
    _transform_pts!(dest, t.f, src, threaded)
_transform_soa!(Xd, Yd, t::Transformation, Xs, Ys, threaded) =
    _transform_soa!(Xd, Yd, t.f, Xs, Ys, threaded)

# Invert the pipeline rather than rebuilding one from the EPSG pair: `T` and
# `kernel` are carried by the operator's own type and are not fields here, so
# rebuilding silently reverted both to their defaults -- a Float32
# transformation inverted to a Float64 one, and one built with `BaseKernel`
# inverted to a lane-safe `FastKernel` one. Inverting the operator gives the
# same pipeline the EPSG pair would, at the precision and kernel it was
# built with.
Base.inv(t::Transformation) =
    Transformation(inv(t.f), t.target_epsg, t.source_epsg,
                   t.always_xy, t.threaded, t.proj_only)

function Base.show(io::IO, t::Transformation)
    print(io,
        """Transformation
            source_epsg:    EPSG:$(first(t.source_epsg.val))
            target_epsg:    EPSG:$(first(t.target_epsg.val))
            threaded:       $(t.threaded)
            always_xy:      $(t.always_xy)
            proj_only:      $(t.proj_only)
        """)
end

# --- whole-collection calls ------------------------------------------------
(t::Transformation)(X::AbstractVector{<:Real}, Y::AbstractVector{<:Real}) =
    transform(t, X, Y; threaded = t.threaded)

# any vector of GeoInterface points. A `Vector{<:Real}` is itself a point, not
# a collection of them, so it stays with the point methods above.
(t::Transformation)(pts::AbstractVector) = transform(t, pts; threaded = t.threaded)
@inline (t::Transformation)(p::AbstractVector{<:Real}) = t(GI.x(p), GI.y(p))
