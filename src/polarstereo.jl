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
function polarstereo_fwd(longitude::Real, latitude::Real; a::Real=6378137.0, e::Real=0.08181919, lat_ts::Real=-70.0, lon_0::Real=0.0, always_xy=true)

    if !always_xy
       # this flips the inputs if input in latitude, longitude order
       latitude, longitude = (longitude, latitude)
    end
    
    # Convert to Radians
    latitude = deg2rad(latitude)
    lat_ts = deg2rad(lat_ts)
    longitude = deg2rad(longitude)
    lon_0 = deg2rad(lon_0)

    # If the standard parallel is in Southern Hemisphere, switch signs.
    if lat_ts < 0.0 
        pm = -1.0 # plus or minus, north lat. or south
        latitude = -latitude
        lat_ts = -lat_ts
        longitude = -longitude
        lon_0 = -lon_0
    else
        pm = 1.0
    end

    # This is not commented very well. See Snyder for details.
    t = tan(pi/4 - latitude/2)/((1-e*sin(latitude))/(1+e*sin(latitude)))^(e/2)
    t_c = tan(pi/4 - lat_ts/2)/((1-e*sin(lat_ts))/(1+e*sin(lat_ts)))^(e/2)
    m_c = cos(lat_ts)/sqrt(1-e^2*(sin(lat_ts))^2)
    rho = a * m_c * t / t_c # True Scale at Lat lat_ts

    #m = cos(latitude)/sqrt(1-e^2*(sin(latitude))^2)
    x = pm * rho * sin(longitude-lon_0)
    y = -pm * rho * cos(longitude-lon_0)
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

function polarstereo_inv(x::Real, y::Real; epsg = nothing, a::Real=nothing, e::Real=nothing, lat_ts::Real=nothing, lon_0::Real=nothing, always_xy=true)

    # convert to radians
    lat_ts = deg2rad(lat_ts)
    lon_0 = deg2rad(lon_0)
    
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
    
    # this is not commented very well. See Snyder for details.
    t_c = tan(pi/4 - lat_ts/2)/((1-e * sin(lat_ts))/(1 + e * sin(lat_ts)))^(e/2)
    m_c = cos(lat_ts)/sqrt(1-e^2*(sin(lat_ts))^2)

    rho = sqrt(x^2 + y^2)
    t = rho * t_c / (a * m_c)
    
    # find latitude with a series instead of iterating.
    chi = pi/2 - 2 * atan(t)

    latitude = chi + (e^2/2 + 5*e^4/24 + e^6/12 + 13*e^8/360)*sin(2*chi) +
        (7*e^4/48 + 29*e^6/240 + 811*e^8/11520)*sin(4*chi) +
        (7*e^6/120+81*e^8/1120)*sin(6*chi) +
        (4279*e^8/161280)*sin(8*chi)
    longitude = lon_0 + atan(x, -y)

    # correct the signs and phasing
    latitude = pm * latitude
    longitude = pm * longitude
    longitude = mod(longitude + pi, 2*pi) - pi

    # convert back to degrees
    latitude = rad2deg(latitude)
    longitude = rad2deg(longitude)
    
    if !always_xy
        return (latitude, longitude)
    else
        return (longitude, latitude)
    end
end
