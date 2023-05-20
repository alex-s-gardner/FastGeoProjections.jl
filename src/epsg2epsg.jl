# List of FastGeoProjections native projections
fast_epsgs = [EPSG(3031), EPSG(3413), EPSG(4326)]

function epsg2epsg(XY::Tuple{T,T}, epsg_from::EPSG, epsg_to::EPSG, proj_only=false) where {T<:Real}

    if any(fast_epsgs .== Ref(epsg_from)) && any(fast_epsgs .== Ref(epsg_to)) && !proj_only
        # if both EPSG codes have been implimented in then use native transformation
        @time begin
            XY = project_to(epsg_to)(project_from(epsg_from)(XY))
        end
    else
        transform = Proj.Transformation("EPSG:$(epsg_from.val)", "EPSG:$(epsg_to.val)")
        XY = transform(getindex(XY, 1), getindex(XY, 2))
    end

    return X, Y
end

function epsg2epsg(XY::Vector{Tuple{T,T}}, epsg_from::EPSG, epsg_to::EPSG; threaded=true, proj_only = false) where {T<:Real}
    if any(fast_epsgs .== Ref(epsg_from)) && any(fast_epsgs .== Ref(epsg_to)) && !proj_only
        # if both EPSG codes have been implimented in then use native transformation
        if threaded 
            xy = similar(XY)
            Threads.@threads for i = eachindex(XY)
               xy[i] = project_to(epsg_to)(project_from(epsg_from)(XY[i]))
            end
        else
            xy = project_to(epsg_to).(project_from(epsg_from).(XY))
        end
    else
        if threaded 
            # This will work with threads (and you can add your own Proj context in ctxs), but not on the GPU - that would pretty much require a Julia implementation of Proj.
            ctxs = [Proj.proj_context_clone() for _ in 1:Threads.nthreads()]
            transforms = [Proj.Transformation("EPSG:$(epsg_from.val)", "EPSG:$(epsg_to.val)"; ctx, always_xy=true) for ctx in ctxs]
            xy = similar(XY)
            Threads.@threads for i in eachindex(XY)
                xy[i] = transforms[Threads.threadid()](XY[1][1], XY[1][2])
            end
        else
            transform = Proj.Transformation("EPSG:$(epsg_from.val)", "EPSG:$(epsg_to.val)", always_xy=true)
            xy = transform.(getindex.(XY, 1), getindex.(XY, 2))
        end
    end
    return xy
end


# project from an EPSG => EPSG(4326)
function project_from(epsg::EPSG)
    if epsg == EPSG(4326)
        f = (xy) -> identity((xy))
    elseif epsg == EPSG(3031)
        f = (xy) -> polarstereo_inv(xy[1], xy[2]; a=6378137.0, e=0.08181919, lat_ts=-71.0, lon_0=0.0)
    elseif epsg == EPSG(3413)
        f = (xy) -> polarstereo_inv(xy[1], xy[2]; a=6378137.0, e=0.08181919, lat_ts=70.0, lon_0=-45)
    end
end

# project from EPSG(4326) => EPSG
function project_to(epsg::EPSG)
    if epsg == EPSG(4326)
        f = (xy) -> identity((xy))
    elseif epsg == EPSG(3031)
        f = (xy) -> polarstereo_fwd(xy[1], xy[2]; a=6378137.0, e=0.08181919, lat_ts=-71.0, lon_0=0.0)
    elseif epsg == EPSG(3413)
        f = (xy) -> polarstereo_fwd(xy[1], xy[2]; a=6378137.0, e=0.08181919, lat_ts=70.0, lon_0=-45)
    end
end

