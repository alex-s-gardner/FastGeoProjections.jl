# The Proj-backed fallback, for any CRS pair FastGeoProjections has no native
# implementation for.
#
# Proj is a weak dependency: the native projections are the reason the package
# exists and none of them needs it, so its load cost is paid only by a caller who
# asks for a CRS outside `fast_epsg_codes`. Until this extension loads,
# `proj_transformation` throws and directs the caller to `import Proj`.

module FastGeoProjectionsProjExt

import Proj
import FastGeoProjections as FGP
using FastGeoProjections: EPSG, GeoTransformation, ProjTransformation,
                          borrow, islanesafe, preservesz, ncoords

"""
    ProjPJ(ctx, pj)

One PROJ transformation and the context it was created against, pooled by
[`ProjTransformation`](@ref).
"""
struct ProjPJ
    ctx::Ptr{Cvoid}
    pj::Proj.Transformation
end

function FGP.proj_transformation(source_epsg::EPSG, target_epsg::EPSG, always_xy::Bool)
    src = "EPSG:$(first(source_epsg.val))"
    tgt = "EPSG:$(first(target_epsg.val))"
    n = max(Threads.nthreads(), 1)
    all = map(1:n) do _
        ctx = Proj.proj_context_clone()
        ProjPJ(ctx, Proj.Transformation(src, tgt; ctx, always_xy))
    end
    pool = Channel{ProjPJ}(n)
    foreach(p -> put!(pool, p), all)
    t = ProjTransformation{ProjPJ}(source_epsg, target_epsg, always_xy, pool, all)
    # the PJ objects reference their context, so they must be released first
    finalizer(t) do x
        close(x.pool)
        foreach(p -> finalize(p.pj), x.all)
        foreach(p -> Proj.proj_context_destroy(p.ctx), x.all)
        empty!(x.all)
    end
    t
end

end # module
