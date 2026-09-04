"""
    ProjTransformation(source_epsg, target_epsg; always_xy = false)

Fallback point operator backed by Proj, used for any CRS pair
FastGeoProjections has no native implementation for.

A PJ object may not be shared across threads, so this holds a small pool of
them -- each with its own cloned context -- which a task checks out for the
range it is working on and returns afterwards. A task that finds the pool
empty waits for one to come back.

The pool is checked out rather than indexed by `threadid()`: a thread id is
not a stable identity for a task, which may migrate between yield points, and
a thread adopted after construction (a `@ccallable` entry from a foreign
thread, say) has an id past the end of any array sized when the object was
built.
"""
struct ProjPJ
    ctx::Ptr{Cvoid}
    pj::Proj.Transformation
end

mutable struct ProjTransformation <: GeoTransformation
    source_epsg::EPSG
    target_epsg::EPSG
    always_xy::Bool
    pool::Channel{ProjPJ}
    all::Vector{ProjPJ}
end

function ProjTransformation(source_epsg::EPSG, target_epsg::EPSG; always_xy::Bool = false)
    src = "EPSG:$(first(source_epsg.val))"
    tgt = "EPSG:$(first(target_epsg.val))"
    n = max(Threads.nthreads(), 1)
    all = map(1:n) do _
        ctx = Proj.proj_context_clone()
        ProjPJ(ctx, Proj.Transformation(src, tgt; ctx, always_xy))
    end
    pool = Channel{ProjPJ}(n)
    foreach(p -> put!(pool, p), all)
    t = ProjTransformation(source_epsg, target_epsg, always_xy, pool, all)
    # the PJ objects reference their context, so they must be released first
    finalizer(t) do x
        close(x.pool)
        foreach(p -> finalize(p.pj), x.all)
        foreach(p -> Proj.proj_context_destroy(p.ctx), x.all)
        empty!(x.all)
    end
    t
end

function borrow(f, t::ProjTransformation)
    p = take!(t.pool)
    try
        f(p.pj)
    finally
        put!(t.pool, p)
    end
end

# Calling one of these directly, outside a `borrow`, still has to be safe --
# it is the documented point-operator interface -- so it checks out a PJ of
# its own. That costs roughly 60ns over the ~70ns the transform itself takes;
# anything transforming more than a handful of points should be going through
# `transform`/`transform!`, which check out once per chunk instead.
function (t::ProjTransformation)(x, y)
    p = take!(t.pool)
    try
        p.pj(x, y)
    finally
        put!(t.pool, p)
    end
end

Base.inv(t::ProjTransformation) =
    ProjTransformation(t.target_epsg, t.source_epsg; always_xy = t.always_xy)

islanesafe(::ProjTransformation) = false

Base.show(io::IO, t::ProjTransformation) =
    print(io, "ProjTransformation(EPSG:", first(t.source_epsg.val),
          " → EPSG:", first(t.target_epsg.val), ")")
