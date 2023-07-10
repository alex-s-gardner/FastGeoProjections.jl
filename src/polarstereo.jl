"""
    polarstereo_fwd(lon, lat; lon_0, lat_ts, ellips, threaded, always_xy)

Returns x and y coordinates [meteres] in Polar Stereographic (PS) coordinates given geodetic 
coordinates (EPSG:4326) of longitude and latitude [decimal degrees]. The PS projection is 
defined kwargs of: lon_0: the meridian along positive Y axis, lat_ts: standard parallel, 
which is the latitude of true scale and an ellipsoid that is define by an equatorial radius 
in meters (a) and its eccentricity (e).  Also returnes scale factor (k).
"""
function polarstereo_fwd(
    lon::Union{AbstractVector{<:AbstractFloat},AbstractFloat}, 
    lat::Union{AbstractVector{<:AbstractFloat},AbstractFloat}; 
    lon_0::Real, 
    lat_ts::Real, 
    ellips::Ellipsoid=ellipsoid(EPSG(7030)), 
    threaded=true,
    always_xy=true
    )

    # This is a Julia implimnetation, written by Alex Gardner - JPL/NASA, 2023 of a Matlab 
    # version written by Andy Bliss, 9/12/2011

    if !always_xy
        (lat, lon) = (lon, lat)
    end

    T = eltype(lon);
    lat_ts = lat_ts * pi / 180
    lon_0 = lon_0 * pi / 180

    # If the standard parallel is in Southern Hemisphere, switch signs.
    if lat_ts < 0
        pm = -1 # plus or minus, north lat. or south
        lat_ts = -lat_ts
        lon_0 = -lon_0
    else
        pm = 1
    end

    t_c = convert(T, tan(pi / 4 - lat_ts / 2) / ((1 - ellips.e * sin(lat_ts)) / (1 + ellips.e * sin(lat_ts)))^(ellips.e / 2))
    m_c = convert(T, cos(lat_ts) / sqrt(1 - ellips.e^2 * (sin(lat_ts))^2))

    p = convert(T, pi)
    e = convert(T, ellips.e)
    a = convert(T, ellips.a)
    d2r = convert(T, pi / 180)
    lon_0 = convert(T, lon_0)

    x = Vector{T}(undef, length(lat))
    y = Vector{T}(undef, length(lat))
    #k = Vector{T}(undef, length(lat))

    if !isa(lat, Array)
        lat = [lat]
        lon = [lon]
    end

    if threaded && Threads.nthreads() > 1
        @turbo thread = true for i = eachindex(lat)
            t = tan((p / 4) - (lat[i] * d2r * pm / 2)) / ((1 - e * sin(lat[i] * d2r * pm)) / (1 + e * sin(lat[i] * d2r * pm)))^(e / 2)
            #m = cos(lat[i]) / sqrt(1 - e^2 * (sin(lat[i]))^2)
            rho = a * m_c * t / t_c   # True Scale at Lat lat_ts
            x[i] = pm * rho * sin(lon[i] * d2r * pm - lon_0)
            y[i] = -pm * rho * cos(lon[i] * d2r * pm - lon_0)
            #k[i] = rho / (a * m)
        end
    else
        @turbo thread = false for i = eachindex(lat)
            t = tan((p / 4) - (lat[i] * d2r * pm / 2)) / ((1 - e * sin(lat[i] * d2r * pm)) / (1 + e * sin(lat[i] * d2r * pm)))^(e / 2)
            #m = cos(lat[i]) / sqrt(1 - e^2 * (sin(lat[i]))^2)
            rho = a * m_c * t / t_c   # True Scale at Lat lat_ts
            x[i] = pm * rho * sin(lon[i] * d2r * pm - lon_0)
            y[i] = -pm * rho * cos(lon[i] * d2r * pm - lon_0)
            #k[i] = rho / (a * m)
        end
    end

    if length(x) == 1
        x = x[1]
        y = y[1]
        #k = k[1]
    end

    return x, y
end


"""
    polarstereo_inv(x, y; lon_0, lat_ts, ellips, threaded, always_xy)

Returns  geodetic coordinates (EPSG:4326) of longitude and latitude [decimal degrees] given 
x and y coordinates [meteres] in Polar Stereographic (PS) coordinates. The PS projection is 
defined kwargs of: lon_0: the meridian along positive Y axis, lat_ts: standard parallel, 
which is the latitude of true scale and an ellipsoid that is define by an equatorial radius 
in meters (a) and its eccentricity (e).
"""
function polarstereo_inv(
    x::Union{AbstractVector{<:AbstractFloat},AbstractFloat}, 
    y::Union{AbstractVector{<:AbstractFloat},AbstractFloat}; 
    lon_0::Real,
    lat_ts::Real, 
    ellips::Ellipsoid=ellipsoid(EPSG(7030)), 
    threaded=false,
    always_xy=true)

    # This is a Julia implimnetation, written by Alex Gardner - JPL/NASA, 2023 of a Matlab 
    # version written by Andy Bliss, 9/12/2011

    # set types
    T = eltype(x)

    # convert to radians
    lat_ts = lat_ts * pi / 180
    lon_0 = lon_0 * pi / 180
    
    # if the standard parallel is in S.Hemi., switch signs.
    if lat_ts < 0
        pm = -1
        lat_ts = -lat_ts
        lon_0 = -lon_0
        x = -x
        y = -y
    else
        pm = 1
    end
    
    # See Snyder for details.
    t_c = convert(T, tan(pi / 4 - lat_ts / 2) / ((1 - ellips.e * sin(lat_ts)) / (1 + ellips.e * sin(lat_ts)))^(ellips.e / 2))
    m_c = convert(T, cos(lat_ts) / sqrt(1 - ellips.e^2 * (sin(lat_ts))^2))

    a = convert(T, ellips.a)
    e = convert(T, ellips.e)
    lon_0 = convert(T,lon_0)

    r2d = convert(T, 180 / pi)
    p = convert(T, pi)
    e = convert(T, ellips.e)

    lat = Vector{T}(undef, length(x))
    lon = Vector{T}(undef, length(x))

    if !isa(x, Array)
        x = [x]
        y = [y]
    end

    # looping caused issues with @turbo so use broadcast
    if threaded && isa(x, Array) && Threads.nthreads() > 1
        @turbo thread = true for i = eachindex(x)
            rho = sqrt((pm * x[i]) ^ 2 + (pm * y[i]) ^ 2)
            t = rho * t_c / (a * m_c)
            chi = p / 2 - 2 * atan(t)

            # find lat with a series instead of iterating
            lat[i] = chi + (e^2 / 2 + 5 * e^4 / 24 + e^6 / 12 + 13 * e^8 / 360) * sin(2 * chi) +
                (7 * e^4 / 48 + 29 * e^6 / 240 + 811 * e^8 / 11520) * sin(4 * chi) +
                (7 * e^6 / 120 + 81 * e^8 / 1120) * sin(6 * chi) +
                (4279 * e^8 / 161280) * sin(8 * chi)
            lon[i] = lon_0 + atan(pm * x[i], -y[i])

            # correct the signs and phasing
            lat[i] = lat[i] * pm * r2d
            lon[i] = pm * (mod(pm * lon[i] + p, 2 * p) - p) * r2d
        end
    else
        @turbo thread = false for i = eachindex(x)
            rho = sqrt((pm * x[i])^2 + (pm * y[i])^2)
            t = rho * t_c / (a * m_c)
            chi = p / 2 - 2 * atan(t)

            # find lat with a series instead of iterating
            lat[i] = chi + (e^2 / 2 + 5 * e^4 / 24 + e^6 / 12 + 13 * e^8 / 360) * sin(2 * chi) +
                     (7 * e^4 / 48 + 29 * e^6 / 240 + 811 * e^8 / 11520) * sin(4 * chi) +
                     (7 * e^6 / 120 + 81 * e^8 / 1120) * sin(6 * chi) +
                     (4279 * e^8 / 161280) * sin(8 * chi)
            lon[i] = lon_0 + atan(pm * x[i], -y[i])

            # correct the signs and phasing
            lat[i] = lat[i] * pm * r2d
            lon[i] = pm * (mod(pm * lon[i] + p, 2 * p) - p) * r2d
        end
    end

    if length(lon) == 1
        lon = lon[1]
        lat = lat[1]
    end

    if always_xy
         return lon, lat
    else
        return lat, lon
    end
end