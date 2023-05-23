"""
This function that performs a forward transformation from geographic coordinates (latitude and longitude) to Polar Stereographic (PS) coordinates. The PS projection is commonly used in mapping polar regions, because it minimizes distortion near the poles.

The input arguments to this function are:

latitude: latitude of the point(s) to transform
lambda: longitude of the point(s) to transform
a: radius of the ellipsoid (default value is for WGS84)
e: eccentricity of the ellipsoid (default value is for WGS84)
lat_ts: standard parallel, which is the latitude of true scale (default is -70 for the South Pole)
lon_0: meridian along positive Y axis (default is 0)

The output arguments are x and y, which are the PS coordinates in meters, and k, which is the true scale at the standard parallel.

Note: This function assumes that the input latitude and longitude are in degrees, and the function outputs x and y coordinates are in meters. 


This is a Julia implimnetation a Matlab version written by Andy Bliss, 9/12/2011
"""
function polarstereo_fwd(longitude, latitude; a::Real=6378137.0, e::Real=0.08181919, lat_ts::Real=-71.0, lon_0::Real=0.0, threaded = false)

    d2r = deg2rad(one(eltype(longitude)))
    p = convert(eltype(longitude), pi)
    a = convert(eltype(longitude), a)
    e = convert(eltype(longitude), e)
    lat_ts = convert(eltype(longitude), lat_ts)
    lon_0 = convert(eltype(longitude), lon_0)

    lat_ts = lat_ts * d2r
    lon_0 = lon_0 * d2r

    # If the standard parallel is in Southern Hemisphere, switch signs.
    if lat_ts < 0
        pm = -1 # plus or minus, north lat. or south
        lat_ts = -lat_ts
        lon_0 = -lon_0
    else
        pm = 1
    end

    t_c = tan(p/4 - lat_ts/2)/((1-e*sin(lat_ts))/(1+e*sin(lat_ts)))^(e/2)
    m_c = cos(lat_ts)/sqrt(1-e^2*(sin(lat_ts))^2)

    if threaded
  
        t = @turbo tan.((p / 4) .- (latitude * d2r * pm ./ 2)) ./ ((1 .- e * sin.(latitude * d2r * pm)) ./ (1 .+ e * sin.(latitude * d2r * pm))) .^ (e / 2)

        rho = @turbo a * m_c * t / t_c   # True Scale at Lat lat_ts

        x = @turbo pm * rho .* sin.(longitude * d2r * pm .- lon_0)
        y = @turbo -pm * rho .* cos.(longitude * d2r * pm .- lon_0)
    else
        t = tan.((p/4) .- (latitude * d2r * pm / 2)) ./ ((1 .- e * sin.(latitude * d2r * pm)) ./ (1 .+ e * sin.(latitude * d2r * pm))) .^ (e / 2)
        rho = a * m_c * t / t_c   # True Scale at Lat lat_ts

        x = pm * rho .* sin.(longitude * d2r * pm .- lon_0)
        y = -pm * rho .* cos.(longitude * d2r * pm .- lon_0)
    end
    # k = rho / (a * m)
    return x, y
end


"""
This is a function for inverse polar stereographic projection. In this function, the input x and y are the coordinates in the polar-stereographic projection plane in meters. The other inputs 'a', 'e', 'phi_c' and 'lon_0' are optional and have default values set to WGS84 ellipsoid parameters and reference meridian for the south pole stereo projection.

x/y: coordinates of a point in the polar stereographic projection in meters
latitude: latitude of the point(s) to transform in degrees
longitude: longitude of the point(s) to transform in degrees
a: radius of the ellipsoid (default value is for WGS84)
e: eccentricity of the ellipsoid (default value is for WGS84)


Note: This function assumes that the input x and y coordinates are in meters, and the function outputs latitudes and longitudes in degrees. Ensure you provide the correct input units when using this function.

This is a Julia implimnetation a Matlab version written by Andy Bliss, 9/12/2011
"""

function polarstereo_inv(x, y; a::Real=6378137.0, e::Real=0.08181919, lat_ts::Real=-71.0, lon_0::Real=0.0, threaded=false)

    # set types
    d2r = deg2rad(one(eltype(x)))
    r2d = rad2deg(one(eltype(x)))
    p = convert(eltype(x), pi)
    a = convert(eltype(x), a)
    e = convert(eltype(x), e)
    lat_ts = convert(eltype(x), lat_ts)
    lon_0 = convert(eltype(x), lon_0)

    # convert to radians
    lat_ts = lat_ts * d2r
    lon_0 = lon_0 * d2r
    
    # if the standard parallel is in S.Hemi., switch signs.
    if lat_ts < 0
        pm = -1
        lat_ts = -lat_ts
        lon_0 = -lon_0
    else
        pm = 1
    end
    
    # this is not commented very well. See Snyder for details.
    t_c = tan(p/4 - lat_ts/2)/((1-e * sin(lat_ts))/(1 + e * sin(lat_ts)))^(e/2)
    m_c = cos(lat_ts)/sqrt(1-e^2*(sin(lat_ts))^2)

    if threaded
        rho = @turbo sqrt.((pm * x).^2 .+ (pm * y).^2)
        t = @turbo rho * (t_c / (a * m_c))
        
        # find latitude with a series instead of iterating.
        chi = @turbo p/2 .- 2 * atan.(t)

        latitude = @turbo chi .+ (e^2/2 + 5 * e^4/24 + e^6/12 + 13 * e^8/360) * sin.(2 * chi) .+
            (7*e^4/48 + 29*e^6/240 + 811*e^8/11520) * sin.(4 * chi) +
            (7*e^6/120+81*e^8/1120) * sin.(6 * chi) .+
            (4279*e^8/161280) * sin.(8 * chi)


        longitude = @turbo lon_0 .+ atan.(pm * x, y)

        # correct the signs and phasing
        latitude = @turbo pm * latitude * r2d
        longitude = @turbo (mod.(pm * longitude .+ p, Ref(2 * p)) .- p) * r2d

    else
        rho = sqrt.((pm * x).^2 .+ (pm * y).^2)
        t = rho * (t_c / (a * m_c))
        
        # find latitude with a series instead of iterating.
        chi = p/2 .- 2 * atan.(t)

        latitude = chi .+ (e^2 / 2 + 5 * e^4 / 24 + e^6 / 12 + 13 * e^8 / 360) * sin.(2 * chi) .+
            (7 * e^4 / 48 + 29 * e^6 / 240 + 811 * e^8 / 11520) * sin.(4 * chi) +
            (7 * e^6 / 120 + 81 * e^8 / 1120) * sin.(6 * chi) .+
            (4279 * e^8 / 161280) * sin.(8 * chi)

        longitude = lon_0 .+ atan.(pm * x, y)

        # correct the signs and phasing
        latitude = pm * latitude * r2d
        longitude = (mod.(pm * longitude .+ p, Ref(2 * p)) .- p) * r2d
    end
    return longitude, latitude
end