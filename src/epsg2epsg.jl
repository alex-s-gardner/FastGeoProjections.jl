# List of FastGeoProjections native projections
fast_epsgs = [EPSG(3031), EPSG(3413), EPSG(4326)]

"""
    epsg2epsg(XY, epsg_from::EPSG, epsg_to::EPSG; threaded=true, proj_only=false)

Transforms a point Tuple, or a vector of point Tuples, `XY` from one coordinate reference systems defined by `epsg_from` to another define by `epsg_to`. Multithreading can be turned on and off with `threaded`. Optimized Julia native code used when available. To force use of Proj 

"""
function epsg2epsg(XY, epsg_from::EPSG, epsg_to::EPSG; threaded=true, proj_only=false)
    if (any(fast_epsgs .== Ref(epsg_from)) || isutm(epsg_from)) && (any(fast_epsgs .== Ref(epsg_to)) || isutm(epsg_to)) && !proj_only
        # if both EPSG codes have been implimented in then use native transformation
        xy = project_to(epsg_to; threaded=threaded)(project_from(epsg_from; threaded=threaded)(XY))
       # xy = [(xy[1][i], xy[2][i]) for i in 1:length(xy[1])]
    else
        if threaded
            # This will work with threads (and you can add your own Proj context in ctxs), but not on the GPU - that would pretty much require a Julia implementation of Proj.
            ctxs = [Proj.proj_context_clone() for _ in 1:Threads.nthreads()]
            transforms = [Proj.Transformation("EPSG:$(epsg_from.val)", "EPSG:$(epsg_to.val)"; ctx, always_xy=true) for ctx in ctxs]
            xy = similar(XY)
            Threads.@threads for i in eachindex(XY)
                xy[i] = transforms[Threads.threadid()](XY[i][1], XY[i][2])
            end
        else
            transform = Proj.Transformation("EPSG:$(epsg_from.val)", "EPSG:$(epsg_to.val)", always_xy=true)
            xy = transform.(Array(getindex.(XY, 1)), Array(getindex.(XY, 2)))
        end
        xy = (getindex.(xy, 1), getindex.(xy, 2))
    end
    return xy
end

# project from an EPSG => EPSG(4326)
function project_from(epsg::EPSG; threaded=true)
    if epsg == EPSG(4326)
        f = (xy) -> identity((xy))
    elseif epsg == EPSG(3031)
        f = (xy) -> polarstereo_inv(getindex.(xy, 1), getindex.(xy, 2); lat_ts=-71.0, lon_0=0.0, ellips=ellipsoid(EPSG(7030)), threaded=threaded)
    elseif epsg == EPSG(3413)
        f = (xy) -> polarstereo_inv(getindex.(xy, 1), getindex.(xy, 2); lat_ts=70.0, lon_0=-45.0, ellips=ellipsoid(EPSG(7030)), threaded=threaded)
    elseif isutm(epsg)
        f = (xy) -> utm_inv(getindex.(xy, 1), getindex.(xy, 2), threaded=threaded, epsg = epsg)
    end
end

# project from EPSG(4326) => EPSG
function project_to(epsg::EPSG; threaded=true)
    if epsg == EPSG(4326)
        f = (xy) -> identity((xy))
    elseif epsg == EPSG(3031)
        f = (xy) -> polarstereo_fwd(getindex.(xy, 1), getindex.(xy, 2); lat_ts=-71.0, lon_0=0.0, ellips=ellipsoid(EPSG(7030)), threaded=threaded)
    elseif epsg == EPSG(3413)
        f = (xy) -> polarstereo_fwd(getindex.(xy, 1), getindex.(xy, 2); lat_ts=70.0, lon_0=-45.0, ellips=ellipsoid(EPSG(7030)), threaded=threaded)
    elseif isutm(epsg)
        f = (xy) -> utm_fwd(getindex.(xy, 1), getindex.(xy, 2); threaded=threaded, epsg = epsg)
    end
end

