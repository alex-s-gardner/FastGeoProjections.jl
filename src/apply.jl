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
isxymajor(::Type{P}) where {P} = _componentval(P) !== nothing

# `Val{T}` for the float type `P` stores its coordinates as, or `nothing`.
#
# A `Val` rather than the type itself, and two spelled-out branches rather
# than a loop over `(Float64, Float32)`, so that the answer travels in the
# type domain. Returned as an ordinary value it is a `DataType`, which leaves
# the element type of the reinterpreted matrix unknown: a dispatch on every
# call, and nothing below it resolvable ahead of time.
@inline function _componentval(::Type{P}) where {P}
    isbitstype(P) || return nothing
    _qualifies(P, Float64) && return Val(Float64)
    _qualifies(P, Float32) && return Val(Float32)
    nothing
end

@inline function _qualifies(::Type{P}, ::Type{T}) where {P,T}
    _isallfloat(P, T) || return false
    n, r = divrem(sizeof(P), sizeof(T))
    r == 0 && n >= 2 && _probexy(P, T, Val(n))
end

# Number of `T`s a qualifying point type is stored as; 0 for anything else.
@inline _ncomponents(::Type{P}) where {P} = _ncomponents(P, _componentval(P))
@inline _ncomponents(::Type{P}, ::Nothing) where {P} = 0
@inline _ncomponents(::Type{P}, ::Val{T}) where {P,T} = sizeof(P) ÷ sizeof(T)

# Whether every field of `P` bottoms out in `T`. A field with no storage (a
# `crs::Nothing`, say) cannot be laid out wrong, so it does not disqualify.
#
# The recursion walks field *indices* carried in a `Val`, not `fieldtypes(P)`:
# a tuple of types is a tuple of `DataType` values, so walking it would leave
# every field type unknown until run time. Here each `fieldtype(P, N)` is a
# literal, which is what lets the predicate fold away entirely.
_isallfloat(::Type{T}, ::Type{T}) where {T<:Union{Float32,Float64}} = true
@inline function _isallfloat(::Type{P}, ::Type{T}) where {P,T}
    sizeof(P) == 0 && return true
    isstructtype(P) || return false
    n = fieldcount(P)
    n == 0 ? false : _fieldsallfloat(P, T, Val(n))
end
@inline _fieldsallfloat(::Type{P}, ::Type{T}, ::Val{0}) where {P,T} = true
@inline _fieldsallfloat(::Type{P}, ::Type{T}, ::Val{N}) where {P,T,N} =
    _isallfloat(fieldtype(P, N), T) && _fieldsallfloat(P, T, Val(N - 1))

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
@inline function _interleaved(v::Array)
    P = eltype(v)
    isxymajor(P) || return nothing
    _reinterp(v, _componentval(P))
end
@inline _reinterp(::Array, ::Nothing) = nothing
@inline _reinterp(v::Array, ::Val{T}) where {T} = reinterpret(reshape, T, v)
_interleaved(::Any) = nothing

# ---------------------------------------------------------------------------
# scalar loops
# ---------------------------------------------------------------------------
"""
    borrow(f, t)

Run `f(op)` with `op` a point operator the calling task has exclusive use of
for the duration. Checking out once per chunk rather than once per point is
the whole point of it, so this wraps a range of work and never a single call.

Defaults to `f(t)`: only a transformation holding something that cannot be
shared between tasks -- a PJ object, say -- needs to do anything here. A
composition holding such a stage is not built by `pipeline`, and would fall
back to that stage's own per-point locking.
"""
@inline borrow(f, t::GeoTransformation) = f(t)

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
_carriesz(::Type{P}) where {P} = _ncomponents(P) >= 3
_carriesz(::Type{Any}) = false

@inline _rebuild(::Type{P}, x, y, src, ::Val{false}) where {P} = rebuildpoint(P, x, y)
@inline _rebuild(::Type{P}, x, y, src, ::Val{true}) where {P} =
    rebuildpoint(P, x, y, GI.z(src))

# `src` holds GeoInterface points; for a tuple or an SVector `GI.x` is a
# `getindex` and inlines away, so this costs nothing over `p[1]`.
#
# Two coordinates go to the operator and a third, if the destination has one, is
# carried across from the source point.
@inline function _scalar_range!(dst, src, t, lo, hi, ::Val{2})
    P = eltype(dst)
    # as with `lanewidth_val`, built outside the loop: `_carriesz` is a property
    # of the destination type, but it does not infer as a constant, so leaving it
    # inside would put a dispatch on every point
    zv = Val(_carriesz(P))
    @inbounds for i in lo:hi
        p = src[i]
        x, y = t(GI.x(p), GI.y(p))
        dst[i] = _rebuild(P, x, y, p, zv)
    end
end

# ...and where the operator takes three, the transformed z comes back from it
# rather than from the source point, and x and y depend on the z going in.
# Selected once per chunk by `_transform_pts!`, so which form is called is not a
# branch per point.
@inline function _scalar_range!(dst, src, t, lo, hi, ::Val{3})
    P = eltype(dst)
    @inbounds for i in lo:hi
        p = src[i]
        x, y, z = t(GI.x(p), GI.y(p), GI.z(p))
        dst[i] = rebuildpoint(P, x, y, z)
    end
end

# How many coordinates travel together between a `S` source and a `D` destination:
# what the operator takes ([`ncoords`](@ref)), capped by what the two ends have
# room for. A three-coordinate operator over two-component points is *not* run at
# an implied zero height -- it falls to 2 here and the operator's own
# two-coordinate method refuses the call.
#
# A `Val`, so the count reaches the loops in the type domain: as an ordinary
# integer it would leave the number of loads per lane-loop iteration, and which
# scalar loop to run, unknown until run time.
@inline function _ncoords_val(t, ::Type{S}, ::Type{D}) where {S,D}
    (ncoords(t) >= 3 && _ncomponents(S) >= 3 && _carriesz(D)) ? Val(3) : Val(2)
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

# array-of-tuples: the components are interleaved, addressed by row.
#
# `C` is how many of them the operator takes and returns -- 2 for a map
# projection, 3 where the operator transforms the height rather than leaving it
# for the destination to keep. Rows past `C` are not touched either way.
@generated function _lane_range_aos!(pd, ps, t, lo, hi, ::Val{W}, ::Val{U},
                                     ::Val{C}) where {W,U,C}
    # one `vload`/`vstore!` per coordinate row, spliced in at each of the three
    # loop bodies below
    load(idx) = :(t($((:(vload(ps, ($c, MM{$W}($idx)))) for c in 1:C)...)))
    loadm(idx) = :(t($((:(vload(ps, ($c, MM{$W}($idx)), m)) for c in 1:C)...)))
    store(o, idx) = Expr(:block, (:(vstore!(pd, $o[$c], ($c, MM{$W}($idx)))) for c in 1:C)...)
    storem(o, idx) = Expr(:block, (:(vstore!(pd, $o[$c], ($c, MM{$W}($idx)), m)) for c in 1:C)...)
    quote
        i = lo
        @inbounds while i + $(W * U) - 1 <= hi
            Base.Cartesian.@nexprs $U u -> begin
                j_u = i + (u - 1) * $W
                o_u = $(load(:j_u))
            end
            Base.Cartesian.@nexprs $U u -> $(store(:o_u, :j_u))
            i += $(W * U)
        end
        @inbounds while i + $W - 1 <= hi
            o = $(load(:i))
            $(store(:o, :i))
            i += $W
        end
        @inbounds if i <= hi
            m = _lanemask(Val($W), hi - i + 1)
            o = $(loadm(:i))
            $(storem(:o, :i))
        end
        nothing
    end
end

# ---------------------------------------------------------------------------
# chunking
# ---------------------------------------------------------------------------
# Chunks are sized in points rather than cut from the thread count. Several
# chunks per thread let a dynamic scheduler even out a thread that draws a
# slow one -- on a machine with both performance and efficiency cores that is
# worth about 18% at 1e6 points on eight threads over the thread-count split
# -- and the partition stops depending on how many threads happen to be
# available.
#
# 2^12 by measurement, on the principle that what matters is having enough
# chunks to go round rather than the size of any one. 2^10 through 2^14 are
# within noise of each other at 1e6 and 1e7 points; 2^16 and above collapse on
# smaller inputs, nearly 4x worse at 1e5 points where they yield one or two
# chunks for eight threads to share.
const CHUNK = 1 << 12

# `:dynamic`, not `:static`: `:static` throws when it is nested or run
# concurrently, and `threaded` defaults to true whenever there is more than
# one thread, so `:static` made `transform` unusable from inside a threaded
# loop of the caller's own. Nothing here needed the pinning -- chunk results
# are independent, and the tail of each is masked rather than scalar.
function _run!(body!, r::AbstractUnitRange, threaded)
    n = length(r)
    n == 0 && return nothing
    if threaded && Threads.nthreads() > 1 && n > CHUNK
        nchunks = cld(n, CHUNK)
        Threads.@threads :dynamic for k in 1:nchunks
            lo = first(r) + (k - 1) * CHUNK
            body!(lo, min(lo + CHUNK - 1, last(r)))
        end
    else
        body!(first(r), last(r))
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
at a time, still threaded.

A `Point3` keeps its z where `t` is a map projection, which is a function of x
and y and leaves a height alone. Where `t` can change one -- a datum shift, a
compound CRS with a geoid model -- the height is transformed along with x and y
instead, and the points go one at a time. See [`preservesz`](@ref).
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
    C = _ncoords_val(t, eltype(src), eltype(dest))
    if islanesafe(t)
        Ms = _interleaved(src)
        Md = dest === src ? Ms : _interleaved(dest)
        if Ms !== nothing && Md !== nothing && eltype(Md) === eltype(Ms)
            _transform_interleaved!(Md, Ms, t, length(src), threaded, C)
            return dest
        end
    end
    # `eachindex`, not `1:n`: the scalar path is reached by anything that is
    # not a dense interleaved buffer, which includes a vector that does not
    # start at 1. Handing it `1:n` walked off the end with bounds checks off.
    _run!(eachindex(src), threaded) do lo, hi
        # once per chunk, not per point: a transformation backed by a
        # resource its task must have to itself checks one out here
        borrow(t) do tt
            _scalar_range!(dest, src, tt, lo, hi, C)
        end
    end
    dest
end

# `Md` and `Ms` need not have the same number of rows: only the first `C` are
# read and written, so a Point3 source can feed a Point2 destination.
function _transform_interleaved!(Md, Ms, t, n, threaded, ::Val{C}) where {C}
    T = eltype(Ms)
    tt = adapt_eltype(t, T)
    W = lanewidth_val(T)
    # the lane loop addresses the buffer by linear offset from 1
    Base.require_one_based_indexing(Md, Ms)
    GC.@preserve Md Ms begin
        pd = stridedpointer(Md)
        ps = stridedpointer(Ms)
        _run!(1:n, threaded) do lo, hi
            _lane_range_aos!(pd, ps, tt, lo, hi, W, Val(UNROLL), Val(C))
        end
    end
    nothing
end

"""
    transform(t, points; threaded = Threads.nthreads() > 1)

Out-of-place [`transform!`](@ref): returns a new vector of transformed points,
of the same type as `points`. Components past x and y are carried over, except
for a height that `t` itself transforms; see [`preservesz`](@ref).

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
        Base.require_one_based_indexing(Xd, Yd, Xs, Ys)
        GC.@preserve Xd Yd Xs Ys begin
            pdx = stridedpointer(Xd); pdy = stridedpointer(Yd)
            psx = stridedpointer(Xs); psy = stridedpointer(Ys)
            _run!(1:n, threaded) do lo, hi
                _lane_range!(pdx, pdy, psx, psy, tt, lo, hi, W, Val(UNROLL))
            end
        end
    else
        _run!(eachindex(Xs), threaded) do lo, hi
            borrow(t) do tt
                _scalar_range!(Xd, Yd, Xs, Ys, tt, lo, hi)
            end
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
