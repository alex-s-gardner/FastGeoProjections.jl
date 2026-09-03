"""
    ProjTransformation(source_epsg, target_epsg; always_xy = false)

Fallback point operator backed by Proj, used for any CRS pair
FastGeoProjections has no native implementation for.

One `Proj.Transformation` (and its own cloned context) is created per thread,
since a PJ object may not be shared across threads.
"""
mutable struct ProjTransformation <: GeoTransformation
    source_epsg::EPSG
    target_epsg::EPSG
    always_xy::Bool
    contexts::Vector{Ptr{Cvoid}}
    transformations::Vector{Proj.Transformation}
end

function ProjTransformation(source_epsg::EPSG, target_epsg::EPSG; always_xy::Bool = false)
    src = "EPSG:$(first(source_epsg.val))"
    tgt = "EPSG:$(first(target_epsg.val))"
    contexts = [Proj.proj_context_clone() for _ in 1:Threads.maxthreadid()]
    transformations = [Proj.Transformation(src, tgt; ctx, always_xy) for ctx in contexts]
    t = ProjTransformation(source_epsg, target_epsg, always_xy, contexts, transformations)
    # the PJ objects reference their context, so they must be released first
    finalizer(t) do x
        foreach(finalize, x.transformations)
        foreach(Proj.proj_context_destroy, x.contexts)
        empty!(x.transformations)
        empty!(x.contexts)
    end
    t
end

@inline function (t::ProjTransformation)(x, y)
    t.transformations[Threads.threadid()](x, y)
end

Base.inv(t::ProjTransformation) =
    ProjTransformation(t.target_epsg, t.source_epsg; always_xy = t.always_xy)

islanesafe(::ProjTransformation) = false

Base.show(io::IO, t::ProjTransformation) =
    print(io, "ProjTransformation(EPSG:", first(t.source_epsg.val),
          " → EPSG:", first(t.target_epsg.val), ")")
