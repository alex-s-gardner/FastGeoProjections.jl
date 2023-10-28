"""
    utm_inv(lon, lat; epsg::EPSG=EPSG(0), zone, isnorth, threaded, always_xy)

Returns geodetic coordinates (EPSG:4326) of longitude and latitude [decimal degrees] given x 
and y coordinates [meteres] in UTM projection. The UTM projection is defined by kwargs of 
EPSG *or* zone and isnorth. Also returnes meridian convergence (gam) and scale factor (k).
"""
function utm_inv(
    x::Union{AbstractVector{<:AbstractFloat},AbstractFloat}, 
    y::Union{AbstractVector{<:AbstractFloat},AbstractFloat}; 
    epsg::EPSG=EPSG(0), 
    threaded = true, 
    always_xy = true
    )

    if epsg !== EPSG(0)
        zone, isnorth = epsg2utmzone(epsg::EPSG)
    end

     T = eltype(x)

    #UTM_INV  Inverse UTM projection
    lon0 = convert(T, -183 + 6 * zone);
    lat0 = zero(T);
    fe = convert(T,5e5); 
    fn = convert(T,100e5 * !isnorth); 
    k0 = convert(T,0.9996);
    x = copy(x)
    y = copy(y)
    
    if threaded && Threads.nthreads() > 1
        @turbo thread = true for i = eachindex(x)
            x[i] = (x[i] - fe) / k0
            y[i] = (y[i] - fn) / k0
        end
    else
        @turbo thread = false for i = eachindex(x)
            x[i] = (x[i] - fe) / k0
            y[i] = (y[i] - fn) / k0
        end
    end

    lon, lat = tranmerc_inv(x, y; lon0=lon0, lat0=lat0, ellips=ellipsoid(EPSG(7030)), threaded=threaded, always_xy=always_xy)

    return lon, lat
end

function ups_inv(
    x::Union{AbstractVector{<:AbstractFloat},AbstractFloat}, 
    y::Union{AbstractVector{<:AbstractFloat},AbstractFloat}; 
    isnorth::Bool=true)
    #UPS_INV  Inverse UPS projection

    fe = 20e5;
    fn = 20e5;
    k0 = 0.994;
    x = x .- fe; 
    y = y .- fn;
    x = x / k0; 
    y = y / k0;

    if isnorth
        lat_ts = 90;
        lon_0 = 0;
    else
        lat_ts = -90;
        lon_0 = 0;
    end

    lat, lon = polarstereo_inv(x, y; lat_ts=lat_ts, lon_0=lon_0, ellips=ellipsoid(EPSG(7030)))
    #k = k * k0;

    return (lon, lat)
end


"""
    utm_fwd(lon, lat; epsg::EPSG=EPSG(0), zone, isnorth, threaded, always_xy)

Returns x and y coordinates [meteres] in UTM projection given geodetic 
coordinates (EPSG:4326) of longitude and latitude [decimal degrees]. The UTM projection is 
defined by kwargs of EPSG *or* zone and isnorth. Also returnes meridian convergence (gam) 
and scale factor (k).
"""

function utm_fwd(
    lon::Union{AbstractVector{<:AbstractFloat},AbstractFloat}, 
    lat::Union{AbstractVector{<:AbstractFloat},AbstractFloat}; 
    epsg::EPSG=EPSG(0), 
    threaded=true, 
    always_xy=true
    )

    if epsg !== EPSG(0)
        zone, isnorth = epsg2utmzone(epsg::EPSG)
    end

    T = eltype(lon)
    lon0 = convert(T, -183 + 6 * zone); 
    lat0 = zero(T);
    fe = convert(T, 5e5);
    fn = convert(T, 100e5 * !isnorth);
    k0 = convert(T, 0.9996);
    x, y = tranmerc_fwd(lon, lat; lon0=lon0, lat0=lat0, ellips=ellipsoid(EPSG(7030)), threaded=threaded, always_xy=always_xy)

    if !isa(x, Array)
        x = x * k0 + fe
        y = y * k0 + fn
        #k = k * k0

    elseif threaded && Threads.nthreads()>1
        @turbo thread = true for i = eachindex(x)
            x[i] = x[i] * k0 + fe
            y[i] = y[i] * k0 + fn
            #k[i] = k[i] * k0;
        end
    else
        @turbo thread = false for i = eachindex(x)
            x[i] = x[i] * k0 + fe
            y[i] = y[i] * k0 + fn
            #k[i] = k[i] * k0;
        end
    end

    return x, y
end

function ups_fwd(lon::Union{AbstractVector{<:AbstractFloat},AbstractFloat}, lat::Union{AbstractVector{<:AbstractFloat},AbstractFloat}; isnorth::Bool=true)
    #UPS_FWD  Forward UPS projection

    if isnorth
        lat_ts = 90
        lon_0 = 0
    else
        lat_ts = -90
        lon_0 = 0
    end

    fe = 20e5; 
    fn = 20e5; 
    k0 = 0.994;
    x, y, gam, k = polarstereo_fwd(lon::Union{AbstractVector{<:AbstractFloat},AbstractFloat}, lat::Union{AbstractVector{<:AbstractFloat},AbstractFloat}, lon_0=lon0, lat_ts=lat_ts, ellips=ellipsoid(EPSG(7030)))
    x = x * k0; 
    y = y * k0; 
    k = k * k0;
    x = x + fe; 
    y = y + fn;

    return x, y
end

"""
    utm_epsg(lon::Real, lat::Real, always_xy=true)

returns the EPSG code for the intersecting universal transverse Mercator (UTM) zone -OR- 
the relevant polar stereographic projection if outside of UTM limits.

modified from: https://github.com/JuliaGeo/Geodesy.jl/blob/master/src/utm.jl    
"""
function utm_epsg(lon::Real, lat::Real, always_xy=true)

    if !always_xy
        lat, lon = (lon, lat)
    end

    if lat > 84
        # NSIDC Sea Ice Polar Stereographic North
        return epsg = 3995
    elseif lat < -80
        # Antarctic Polar Stereographic
        return epsg = 19992
    end

    # make sure lon is from -180 to 180
    lon = lon - floor((lon + 180) / (360)) * 360

    # int versions
    ilat = floor(Int64, lat)
    ilon = floor(Int64, lon)

    # get the latitude band
    band = max(-10, min(9, fld((ilat + 80), 8) - 10))

    # and check for weird ones
    zone = fld((ilon + 186), 6)
    if ((band == 7) && (zone == 31) && (ilon >= 3)) # Norway
        zone = 32
    elseif ((band == 9) && (ilon >= 0) && (ilon < 42)) # Svalbard
        zone = 2 * fld((ilon + 183), 12) + 1
    end

    if lat >= 0
        epsg = 32600 + zone
    else
        epsg = 32700 + zone
    end

    # convert to proj string
    epsg = EPSG(epsg)

    return epsg
end

function utmzone2epsg(zone::Int = 0, isnorth::Bool = true)

    if zone == 0 
        if isnorth
            # NSIDC Sea Ice Polar Stereographic North
            return epsg = EPSG(3995)
        else
            # Antarctic Polar Stereographic
            return epsg = EPSG(19992)
        end
    end

    if isnorth
        epsg = 32600 + zone
    else
        epsg = 32700 + zone
    end

    # convert to EPSG type
    epsg = EPSG(epsg)

    return epsg
end

function epsg2utmzone(epsg::EPSG)
    if epsg.val[1] == 3995
        isnorth = true
        zone = 0
    elseif epsg.val[1] == 19992
        isnorth = false
        zone = 0
    elseif Int32(floor(epsg.val[1][1], digits = -2)) == 32600
        isnorth = true
        zone = epsg.val[1] - 32600
    elseif Int32(floor(epsg.val[1], digits = -2)) == 32700
        isnorth = false
        zone = epsg.val[1] - 32700
    else
        error("supplied epsg is not a UTM epsg")
    end

    return (zone = zone, isnorth = isnorth)
end

function isutm(epsg::EPSG)
    tf = Int32(floor(epsg.val[1], digits = -2)) == 32600 || Int32(floor(epsg.val[1], digits = -2)) == 32700
    return tf
end