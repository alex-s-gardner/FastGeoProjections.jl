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
    _interleaved(v)

`v` seen as a dense `ncomponents × length(v)` matrix of floats whose first two
rows are x and y -- or `nothing` when `v` is not laid out that way and its
points have to be visited one at a time.

Any dense vector of isbits points with two or more float components qualifies:
`NTuple{2}` and `NTuple{3}`, `SVector{2}`, `GeometryBasics.Point2` and
`Point3`, `GeoInterface.Wrappers.Point`, a plain `struct P; x; y; end`. A
stride of three or four costs no more than two, because the load deinterleaves
in hardware, so carrying a z around does not push a vector off the fast path.

The layout is *checked* against `GI.x` and `GI.y` rather than assumed: a point
type that happens to store its components in the other order falls back to the
scalar loop instead of silently transposing every coordinate.
"""
function _interleaved(v::Array)
    isempty(v) && return nothing
    P = eltype(v)
    isbitstype(P) || return nothing
    p = @inbounds v[begin]
    GI.geomtrait(p) isa GI.PointTrait || return nothing
    T = typeof(GI.x(p))
    (T <: Union{Float32,Float64} && typeof(GI.y(p)) === T) || return nothing
    ncomp, r = divrem(sizeof(P), sizeof(T))
    (r == 0 && ncomp >= 2) || return nothing
    M = reinterpret(reshape, T, v)
    for i in (firstindex(v), lastindex(v))
        q = @inbounds v[i]
        (isequal(@inbounds(M[1, i]), GI.x(q)) && isequal(@inbounds(M[2, i]), GI.y(q))) ||
            return nothing
    end
    M
end
_interleaved(::Any) = nothing

# ---------------------------------------------------------------------------
# scalar loops
# ---------------------------------------------------------------------------
"""
    rebuildpoint(P, x, y)

Build a point of type `P` out of a transformed coordinate pair. Defaults to
`convert(P, (x, y))`, which already covers `NTuple{2}`, `SVector{2}` and
`GeometryBasics.Point2`.

Only the scalar path calls this; a vector whose points are a dense interleaved
buffer is written component by component and never reaches it. So a point type
of your own needs this method only if its vectors do *not* have that layout --
the components stored in the other order, say, or the type not isbits.
"""
@inline rebuildpoint(::Type{P}, x, y) where {P} = convert(P, (x, y))
@inline rebuildpoint(::Type{Any}, x, y) = (x, y)

# `src` holds GeoInterface points; for a tuple or an SVector `GI.x` is a
# `getindex` and inlines away, so this costs nothing over `p[1]`.
@inline function _scalar_range!(dst, src, t, lo, hi)
    P = eltype(dst)
    @inbounds for i in lo:hi
        p = src[i]
        dst[i] = rebuildpoint(P, t(GI.x(p), GI.y(p))...)
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
