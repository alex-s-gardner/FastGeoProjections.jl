"""
    epsg2epsg(source_epsg::EPSG, target_epsg::EPSG; threaded=true, proj_only=false)

Returns the Transform for points defined by `x` and `y`from one coordinate reference systems defined by `source_epsg` to another define by `target_epsg`. Coodinates `x` and `y` can be either a scalar or a vector. Multithreading can be turned on and off with `threaded`. Optimized Julia native code used when available. To force use of Proj set proj_only = true

"""
function epsg2epsg(source_epsg::EPSG, target_epsg::EPSG; threaded=true, proj_only=false, always_xy=true)
    if isfastepsg(source_epsg, target_epsg) && !proj_only
        # if both EPSG codes have been implimented in then use native transformation
        f = function (x::Union{Real,Vector{<:Real}}, y::Union{Real,Vector{<:Real}}) 
            x, y = project_from(source_epsg; threaded=threaded, always_xy=always_xy)(x,y)
            xx, yy = project_to(target_epsg; threaded=threaded, always_xy=always_xy)(x,y)
            return xx, yy
        end
    else
        if threaded
            # This will work with threads (and you can add your own Proj context in ctxs)
            ctxs = [Proj.proj_context_clone() for _ in 1:Threads.nthreads()]
            trans = [Proj.Transformation("EPSG:$(source_epsg.val)", "EPSG:$(target_epsg.val)"; ctx, always_xy=always_xy) for ctx in ctxs]

            f = function (x::Union{Real,Vector{<:Real}}, y::Union{Real,Vector{<:Real}})
                xx = zeros(size(x))
                yy = zeros(size(x))
                Threads.@threads for i in eachindex(x)
                    xx[i], yy[i] = trans[Threads.threadid()](x[i], y[i])
                end
                return xx, yy
            end
        else
            f = function (x::Union{Real,Vector{<:Real}}, y::Union{Real,Vector{<:Real}})
                xx = zeros(size(x))
                yy = zeros(size(x))
                trans = Proj.Transformation("EPSG:$(source_epsg.val)", "EPSG:$(target_epsg.val)", always_xy=always_xy)
                for i in eachindex(x)
                    xx[i], yy[i] = trans(x[i],y[i])
                end
                return xx, yy
            end

        end
    end
end

## ⬇ ADD FAST PROJECTIONS HERE ⬇ ##

# project from an EPSG => EPSG(4326)
function project_from(epsg::EPSG; threaded=true, always_xy=true)
    if epsg == EPSG(4326)
        f = (x,y) -> identity((x,y))
    elseif epsg == EPSG(3031)
        f = (x,y) -> polarstereo_inv(x, y; lat_ts=-71.0, lon_0=0.0, ellips=ellipsoid(EPSG(7030)), threaded=threaded, always_xy=always_xy)
    elseif epsg == EPSG(3413)
        f = (x,y) -> polarstereo_inv(x, y; lat_ts=70.0, lon_0=-45.0, ellips=ellipsoid(EPSG(7030)), threaded=threaded, always_xy=always_xy)
    elseif isutm(epsg)
        f = (x,y) -> utm_inv(x, y, threaded=threaded, epsg=epsg, always_xy=always_xy)
    end
end

# project from EPSG(4326) => EPSG
function project_to(epsg::EPSG; threaded=true, always_xy=true)
    if epsg == EPSG(4326)
        f = (x,y) -> identity((x,y))
    elseif epsg == EPSG(3031)
        f = (x,y) -> polarstereo_fwd(x, y; lat_ts=-71.0, lon_0=0.0, ellips=ellipsoid(EPSG(7030)), threaded=threaded, always_xy=always_xy)
    elseif epsg == EPSG(3413)
        f = (x,y) -> polarstereo_fwd(x, y; lat_ts=70.0, lon_0=-45.0, ellips=ellipsoid(EPSG(7030)), threaded=threaded, always_xy=always_xy)
    elseif isutm(epsg)
        f = (x,y) -> utm_fwd(x, y; threaded=threaded, epsg=epsg, always_xy=always_xy)
    end
end

# List of FastGeoProjections native projections
fast_epsgs = [EPSG(3031), EPSG(3413), EPSG(4326)]

function isfastepsg(source_epsg, target_epsg)
    tf = (any(fast_epsgs .== Ref(source_epsg)) || isutm(source_epsg)) && (any(fast_epsgs .== Ref(target_epsg)) || isutm(target_epsg))
    return tf
end


