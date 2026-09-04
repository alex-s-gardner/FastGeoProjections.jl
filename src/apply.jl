# Applying a point operator to a collection.
#
# Threading and vectorization live here, not in the transformation. A
# transformation is a pure point operator; this file decides how many points
# to hand it at once (scalar or SIMD lanes) and on how many threads.

# unrolling the lane loop costs nothing and buys a little on short pipelines
const UNROLL = 8

_default_threaded() = Threads.nthreads() > 1

# Number of Float64/Float32 lanes the machine transforms at once.
@inline lanewidth(::Type{T}) where {T} = Int(VectorizationBase.pick_vector_width(T))

# ...as a `Val`, so the width survives being captured by the chunk closures in
# `_transform_interleaved!` and `_transform_soa!`. A closure captures a *value*:
# `W::Int` would reach the lane loop as an ordinary integer and force it to be
# dispatched at run time, where `Val{W}` carries the width in its type.
@inline lanewidth_val(::Type{T}) where {T} = Val(lanewidth(T))

# A dense strided array of floats is the only layout the SIMD path can address.
_strideable(::Vector{T}) where {T<:Union{Float32,Float64}} = true
_strideable(::Any) = false

"""
    isxymajor(P)

Whether a dense `Vector{P}` is a valid interleaved coordinate buffer: whether
the first float of each element is `GI.x` and the second is `GI.y`.

This cannot be asked of GeoInterface. `getcoord` is an arbitrary function of
the point -- ArchGDAL reads it out of the GDAL C API, where there is no Julia
side memory to be ordered at all -- and `coordnames` names coordinates rather
than describing storage, so it is neither necessary nor sufficient here: a
type stored `(y, x)` may report the default `(:X, :Y)`, and
`NamedTuple{(:Y,:X)}`, which GeoInterface itself accepts as a point, reports
`(:Y, :X)` while a `(x, y)` type could too. Both pass `GI.testgeometry`.

So the layout is established from `P` alone, by laying out sentinel
coordinates in memory and asking the resulting point where its x and y are.
That is only sound when every one of `P`'s fields bottoms out in `T`, which
is what `_isallfloat` establishes: then every bit pattern is a valid
instance, and the probe cannot hand an accessor something it should not
dereference. `isbitstype` alone would not do -- `struct H; p::Ptr{Cvoid};
q::Ptr{Cvoid}; end` is isbits and the size of two `Float64`s.

`NTuple{2}` and `NTuple{3}`, `SVector{2}`, `GeometryBasics.Point2` and
`Point3`, `GeoInterface.Wrappers.Point` and a plain `struct P; x; y; end` all
qualify. Define `isxymajor(::Type{P}) = false` for a point type of your own to
keep it off the fast path. Defining it as `true` does nothing on its own: a
type the probe cannot see does not store its coordinates as a dense run of
floats, which is the layout the SIMD path addresses.
"""
isxymajor(::Type{P}) where {P} = _pointlayout(P) !== nothing

# `(T, ncomponents)` for a point type stored as `ncomponents` `T`s with x
# first and y second, or `nothing` for anything else.
function _pointlayout(::Type{P}) where {P}
    isbitstype(P) || return nothing
    for T in (Float64, Float32)
        _isallfloat(P, T) || continue
        n, r = divrem(sizeof(P), sizeof(T))
        (r == 0 && n >= 2 && _probexy(P, T, Val(n))) || return nothing
        return (T, n)
    end
    nothing
end

# Whether every field of `P` bottoms out in `T`. A field with no storage (a
# `crs::Nothing`, say) cannot be laid out wrong, so it does not disqualify.
_isallfloat(::Type{T}, ::Type{T}) where {T<:Union{Float32,Float64}} = true
function _isallfloat(::Type{P}, ::Type{T}) where {P,T}
    sizeof(P) == 0 && return true
    isstructtype(P) || return false
    fs = fieldtypes(P)
    !isempty(fs) && all(F -> _isallfloat(F, T), fs)
end

# Sentinels rather than sampled values: a `(y, x)` type whose sampled points
# happen to lie on `x == y` passes a comparison against its own contents, and
# a freshly allocated destination has no contents to compare against at all.
@inline function _probexy(::Type{P}, ::Type{T}, ::Val{N}) where {P,T,N}
    p = reinterpret(P, ntuple(i -> i <= 2 ? T(i) : zero(T), Val(N)))
    GI.geomtrait(p) isa GI.PointTrait && GI.x(p) === T(1) && GI.y(p) === T(2)
end

"""
    _interleaved(v)

`v` seen as a dense `ncomponents × length(v)` matrix of floats whose first two
rows are x and y -- or `nothing` when `v` is not laid out that way and its
points have to be visited one at a time.

A stride of three or four costs no more than two, because the load
deinterleaves in hardware, so carrying a z around does not push a vector off
the fast path. Whether the layout holds is a property of the element type
alone; see [`isxymajor`](@ref).
"""
function _interleaved(v::Array)
    P = eltype(v)
    isxymajor(P) || return nothing
    l = _pointlayout(P)
    l === nothing && return nothing
    reinterpret(reshape, l[1], v)
end
_interleaved(::Any) = nothing

# ---------------------------------------------------------------------------
# scalar loops
# ---------------------------------------------------------------------------
"""
    rebuildpoint(P, x, y)
    rebuildpoint(P, x, y, z)

Build a point of type `P` out of a transformed coordinate pair. Defaults to
`convert(P, (x, y))`, which already covers `NTuple{2}`, `SVector{2}` and
`GeometryBasics.Point2`.

The three-argument form is used where the destination has a third component
to fill, and gets its z from the source point, so `Point3` and `NTuple{3}`
keep their z on this path as they do on the SIMD one. Which form is called is
decided once per chunk from the destination type, not per point.

Only the scalar path calls either; a vector whose points are a dense
interleaved buffer is written component by component and never reaches them.
So a point type of your own needs a method only if its vectors do *not* have
that layout -- the components stored in the other order, say, or the type not
isbits.
"""
@inline rebuildpoint(::Type{P}, x, y) where {P} = convert(P, (x, y))
@inline rebuildpoint(::Type{Any}, x, y) = (x, y)
@inline rebuildpoint(::Type{P}, x, y, z) where {P} = convert(P, (x, y, z))
@inline rebuildpoint(::Type{Any}, x, y, z) = (x, y, z)

# Whether a destination of type `P` has a third component to fill. Keyed on
# the destination rather than on the source, since a three-component source
# feeding a two-component destination is allowed: the SIMD path writes rows 1
# and 2 of each and leaves whatever else the destination has alone.
_carriesz(::Type{P}) where {P} = (l = _pointlayout(P); l !== nothing && l[2] >= 3)
_carriesz(::Type{Any}) = false

@inline _rebuild(::Type{P}, x, y, src, ::Val{false}) where {P} = rebuildpoint(P, x, y)
@inline _rebuild(::Type{P}, x, y, src, ::Val{true}) where {P} =
    rebuildpoint(P, x, y, GI.z(src))

# `src` holds GeoInterface points; for a tuple or an SVector `GI.x` is a
# `getindex` and inlines away, so this costs nothing over `p[1]`.
#
# ...and as with `lanewidth_val`, the `Val` is built outside the loop:
# `_carriesz` is a property of the destination type, but it does not infer as
# a constant, so leaving it inside would put a dispatch on every point.
@inline function _scalar_range!(dst, src, t, lo, hi)
    P = eltype(dst)
    zv = Val(_carriesz(P))
    @inbounds for i in lo:hi
        p = src[i]
        x, y = t(GI.x(p), GI.y(p))
        dst[i] = _rebuild(P, x, y, p, zv)
    end
end

@inline function _scalar_range!(dstx, dsty, srcx, srcy, t, lo, hi)
    @inbounds for i in lo:hi
        dstx[i], dsty[i] = t(srcx[i], srcy[i])
    end
end

# ---------------------------------------------------------------------------
# SIMD-lane loops: the same operator, handed `Vec`s instead of scalars
#
# The leftover elements at the end of a range go through a *masked* vector step
# rather than a scalar loop. Scalar arithmetic does not contract to FMA where
# the vector path does, so a scalar tail would make the result depend on how
# the range was chunked -- i.e. on the thread count.
# ---------------------------------------------------------------------------
@inline _lanemask(::Val{W}, n) where {W} = VectorizationBase.mask(Val(W), n)
@generated function _lane_range!(pdx, pdy, psx, psy, t, lo, hi, ::Val{W}, ::Val{U}) where {W,U}
    quote
        i = lo
        @inbounds while i + $(W * U) - 1 <= hi
            Base.Cartesian.@nexprs $U u -> begin
                j_u = i + (u - 1) * $W
                o_u = t(vload(psx, (MM{$W}(j_u),)), vload(psy, (MM{$W}(j_u),)))
            end
            Base.Cartesian.@nexprs $U u -> begin
                vstore!(pdx, o_u[1], (MM{$W}(j_u),))
                vstore!(pdy, o_u[2], (MM{$W}(j_u),))
            end
            i += $(W * U)
        end
        @inbounds while i + $W - 1 <= hi
            o = t(vload(psx, (MM{$W}(i),)), vload(psy, (MM{$W}(i),)))
            vstore!(pdx, o[1], (MM{$W}(i),))
            vstore!(pdy, o[2], (MM{$W}(i),))
            i += $W
        end
        @inbounds if i <= hi
            m = _lanemask(Val($W), hi - i + 1)
            o = t(vload(psx, (MM{$W}(i),), m), vload(psy, (MM{$W}(i),), m))
            vstore!(pdx, o[1], (MM{$W}(i),), m)
            vstore!(pdy, o[2], (MM{$W}(i),), m)
        end
        nothing
    end
end

# array-of-tuples: the two components are interleaved, addressed by row
@generated function _lane_range_aos!(pd, ps, t, lo, hi, ::Val{W}, ::Val{U}) where {W,U}
    quote
        i = lo
        @inbounds while i + $(W * U) - 1 <= hi
            Base.Cartesian.@nexprs $U u -> begin
                j_u = i + (u - 1) * $W
                o_u = t(vload(ps, (1, MM{$W}(j_u))), vload(ps, (2, MM{$W}(j_u))))
            end
            Base.Cartesian.@nexprs $U u -> begin
                vstore!(pd, o_u[1], (1, MM{$W}(j_u)))
                vstore!(pd, o_u[2], (2, MM{$W}(j_u)))
            end
            i += $(W * U)
        end
        @inbounds while i + $W - 1 <= hi
            o = t(vload(ps, (1, MM{$W}(i))), vload(ps, (2, MM{$W}(i))))
            vstore!(pd, o[1], (1, MM{$W}(i)))
            vstore!(pd, o[2], (2, MM{$W}(i)))
            i += $W
        end
        @inbounds if i <= hi
            m = _lanemask(Val($W), hi - i + 1)
            o = t(vload(ps, (1, MM{$W}(i)), m), vload(ps, (2, MM{$W}(i)), m))
            vstore!(pd, o[1], (1, MM{$W}(i)), m)
            vstore!(pd, o[2], (2, MM{$W}(i)), m)
        end
        nothing
    end
end

# ---------------------------------------------------------------------------
# chunking
# ---------------------------------------------------------------------------
@inline function _chunk(n, k, nchunks)
    c = cld(n, nchunks)
    ((k - 1) * c + 1, min(k * c, n))
end

function _run!(body!, n, threaded)
    if threaded && n >= 2 * Threads.nthreads()
        nchunks = Threads.nthreads()
        Threads.@threads :static for k in 1:nchunks
            lo, hi = _chunk(n, k, nchunks)
            lo <= hi && body!(lo, hi)
        end
    else
        body!(1, n)
    end
    nothing
end

# ---------------------------------------------------------------------------
# transform! / transform
# ---------------------------------------------------------------------------
"""
    transform!(t, points; threaded = Threads.nthreads() > 1)
    transform!(dest, t, points; threaded)

Apply the transformation `t` to every point of `points`, writing the results
back into `points` (or into `dest`). Returns the destination.

`points` is any vector of GeoInterface points -- `NTuple{2}`, `SVector{2}`,
`GeometryBasics.Point2`, your own two-field struct. Where the vector turns out
to be a dense interleaved buffer of floats (see [`_interleaved`](@ref)) and `t`
is lane-safe, the points are transformed on SIMD lanes; otherwise they go one
at a time, still threaded. Either way only x and y are touched, so a `Point3`
keeps its z.
"""
transform!(t::GeoTransformation, pts::AbstractVector; threaded = _default_threaded()) =
    transform!(pts, t, pts; threaded)

function transform!(dest::AbstractVector, t::GeoTransformation, src::AbstractVector;
                    threaded = _default_threaded())
    axes(dest) == axes(src) || throw(DimensionMismatch("destination and source must match"))
    isempty(src) && return dest
    GI.geomtrait(@inbounds src[begin]) isa GI.PointTrait ||
        throw(ArgumentError("transform! expects a vector of points, got elements of type $(eltype(src))"))
    _transform_pts!(dest, t, src, threaded)
    dest
end

function _transform_pts!(dest, t::GeoTransformation, src, threaded)
    n = length(src)
    if islanesafe(t)
        Ms = _interleaved(src)
        Md = dest === src ? Ms : _interleaved(dest)
        if Ms !== nothing && Md !== nothing && eltype(Md) === eltype(Ms)
            _transform_interleaved!(Md, Ms, t, n, threaded)
            return dest
        end
    end
    _run!(n, threaded) do lo, hi
        _scalar_range!(dest, src, t, lo, hi)
    end
    dest
end

# `Md` and `Ms` need not have the same number of rows: only rows 1 and 2 are
# read and written, so a Point3 source can feed a Point2 destination.
function _transform_interleaved!(Md, Ms, t, n, threaded)
    T = eltype(Ms)
    tt = adapt_eltype(t, T)
    W = lanewidth_val(T)
    GC.@preserve Md Ms begin
        pd = stridedpointer(Md)
        ps = stridedpointer(Ms)
        _run!(n, threaded) do lo, hi
            _lane_range_aos!(pd, ps, tt, lo, hi, W, Val(UNROLL))
        end
    end
    nothing
end

"""
    transform(t, points; threaded = Threads.nthreads() > 1)

Out-of-place [`transform!`](@ref): returns a new vector of transformed points,
of the same type as `points`. Components past x and y are carried over.

A point type that cannot be built from an `(x, y)` tuple only works here if its
vector has the interleaved layout the SIMD path writes through; otherwise pass
your own destination to [`transform!`](@ref).
"""
function transform(t::GeoTransformation, pts::AbstractVector; threaded = _default_threaded())
    dest = similar(pts)
    M = _interleaved(pts)
    # anything past x and y belongs to the caller, so carry it over
    (M !== nothing && size(M, 1) > 2) && copyto!(dest, pts)
    transform!(dest, t, pts; threaded)
end

"""
    transform!(t, X, Y; threaded = Threads.nthreads() > 1)
    transform!(Xdest, Ydest, t, X, Y; threaded)

Struct-of-arrays form: apply `t` to each `(X[i], Y[i])` pair in place, or into
`(Xdest, Ydest)`.
"""
transform!(t::GeoTransformation, X::AbstractVector, Y::AbstractVector;
           threaded = _default_threaded()) = transform!(X, Y, t, X, Y; threaded)

function transform!(Xd::AbstractVector, Yd::AbstractVector, t::GeoTransformation,
                    Xs::AbstractVector, Ys::AbstractVector; threaded = _default_threaded())
    axes(Xs) == axes(Ys) == axes(Xd) == axes(Yd) ||
        throw(DimensionMismatch("all coordinate vectors must match"))
    _transform_soa!(Xd, Yd, t, Xs, Ys, threaded)
    (Xd, Yd)
end

function _transform_soa!(Xd, Yd, t::GeoTransformation, Xs, Ys, threaded)
    n = length(Xs)
    T = promote_type(eltype(Xs), eltype(Ys))
    if islanesafe(t) && T <: Union{Float32,Float64} &&
            _strideable(Xs) && _strideable(Ys) && _strideable(Xd) && _strideable(Yd)
        tt = adapt_eltype(t, T)
        W = lanewidth_val(T)
        GC.@preserve Xd Yd Xs Ys begin
            pdx = stridedpointer(Xd); pdy = stridedpointer(Yd)
            psx = stridedpointer(Xs); psy = stridedpointer(Ys)
            _run!(n, threaded) do lo, hi
                _lane_range!(pdx, pdy, psx, psy, tt, lo, hi, W, Val(UNROLL))
            end
        end
    else
        _run!(n, threaded) do lo, hi
            _scalar_range!(Xd, Yd, Xs, Ys, t, lo, hi)
        end
    end
    (Xd, Yd)
end

"""
    transform(t, X, Y; threaded = Threads.nthreads() > 1)

Out-of-place struct-of-arrays [`transform!`](@ref); returns `(X′, Y′)`.
"""
function transform(t::GeoTransformation, X::AbstractVector, Y::AbstractVector;
                   threaded = _default_threaded())
    T = float(promote_type(eltype(X), eltype(Y)))
    transform!(similar(X, T), similar(Y, T), t, X, Y; threaded)
end
