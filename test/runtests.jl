using FastGeoProjections
using Test
using Proj

@testset "FastGeoProjections.jl" begin
    ## [1] Test WGS 84 / NSIDC Sea Ice Polar Stereographic North 
    latitude = 84.;
    longitude = 50.;

    x, y = FastGeoProjections.polarstereo_fwd(longitude, latitude; lat_ts=70, lon_0=-45, ellips = ellipsoid(EPSG(7030)))
    
    lon, lat = FastGeoProjections.polarstereo_inv(x, y; lat_ts=70, lon_0=-45, ellips=ellipsoid(EPSG(7030)))

    # check accuracy
    @test round(lon, digits=9) == round(longitude, digits=9)
    @test round(lat, digits=9) == round(latitude, digits=9)

    # now inverse 
    longitude0, latitude0 = FastGeoProjections.polarstereo_inv(x, y; lon_0=-45, lat_ts=70, ellips = ellipsoid(EPSG(7030)))
    @test round.(longitude0, digits=10) == longitude

    ## [1] Test WGS 84 / NSIDC Sea Ice Polar Stereographic North 
    latitude = 84.0
    longitude = 50.0

    x, y = FastGeoProjections.polarstereo_fwd(longitude, latitude; lat_ts=70, lon_0=-45, ellips=ellipsoid(EPSG(7030)))

    # build transformation 
    trans = Proj.Transformation("EPSG:4326", "EPSG:3413", always_xy=true)

    # project points
    data = trans.(longitude, latitude)

    # check accuracy within 1 mm
    @test round(x, digits=9) == round(data[1], digits=9)
    @test round(y, digits=9) == round(data[2], digits=9)

    # now inverse 
    longitude0, latitude0 = FastGeoProjections.polarstereo_inv(x, y; lat_ts=70, lon_0=-45, ellips=ellipsoid(EPSG(7030)))
    @test round.(longitude0, digits=9) == longitude
    @test round.(latitude0, digits=9) == latitude

    ## [2] Test WGS 84 / Antarctic Polar Stereographic
    latitude = -84.
    longitude = 50.

    x, y = FastGeoProjections.polarstereo_fwd(longitude, latitude; lon_0=0, lat_ts=-71, ellips = ellipsoid(EPSG(7030)), threaded = true)

    # build transformation 
    trans = Proj.Transformation("EPSG:4326", "EPSG:3031", always_xy=true)

    # project points
    data = trans.(longitude, latitude)

    # check accuracy within 1 mm
    @test round(x, digits=8) == round(data[1], digits=8)
    @test round(y, digits=8) == round(data[2], digits=8)

    # now inverse 
    longitude0, latitude0 = FastGeoProjections.polarstereo_inv(x, y; lon_0=0, lat_ts=-71, ellips = ellipsoid(EPSG(7030)))
    @test round.(longitude0, digits=9) == longitude
    @test round.(latitude0, digits=9) == latitude

    ## [3] Test Transverse Mercator projection [UTM]
    x = [10000., 20000.] 
    y = [10000., 20000.] 

    longitude0, latitude0 = FastGeoProjections.utm_inv(x, y; epsg = EPSG(32619))

    trans = Proj.Transformation("EPSG:32619", "EPSG:4326", always_xy=true)
    data = trans.(x, y)
    @test round.(longitude0, digits=7) == round.(getindex.(data, 1), digits=7)
    @test round.(latitude0, digits=7) == round.(getindex.(data, 2), digits=7)

    longitude0, latitude0 = FastGeoProjections.utm_inv(x, y; epsg=EPSG(32619), threaded = false)
    @test round.(longitude0, digits=7) == round.(getindex.(data, 1), digits=7)
    @test round.(latitude0, digits=7) == round.(getindex.(data, 2), digits=7)

    lat = -[80., 40., 1.]
    lon = [30., 31., 34.]

    foo = FastGeoProjections.utm_epsg.(lon, lat)
    x0, y0 = FastGeoProjections.utm_fwd(lon, lat, epsg=foo[2])

    trans = Proj.Transformation("EPSG:4326", "EPSG:32736", always_xy=true)
    data = trans.(lon, lat)

    @test round.(x0, digits=7) == round.(getindex.(data, 1), digits=7)
    @test round.(y0, digits=7) == round.(getindex.(data, 2), digits=7)

    x0, y0 = FastGeoProjections.utm_fwd(lon, lat, epsg=foo[2], threaded=false)
    @test round.(x0, digits=7) == round.(getindex.(data, 1), digits=7)
    @test round.(y0, digits=7) == round.(getindex.(data, 2), digits=7)

    lat = [80.0, 40.0, 1.0]
    lon = [30.0, 31.0, 34.0]

    foo = FastGeoProjections.utm_epsg.(lon, lat)
    x0, y0 = FastGeoProjections.utm_fwd(lon, lat, epsg=foo[2])

    trans = Proj.Transformation("EPSG:4326", "EPSG:32636", always_xy=true)
    data = trans.(lon, lat)

    @test round.(x0, digits=6) == round.(getindex.(data, 1), digits=6)
    @test round.(y0, digits=6) == round.(getindex.(data, 2), digits=6)

    x0, y0 = FastGeoProjections.utm_fwd(lon, lat, epsg=foo[2], threaded=false)
    @test round.(x0, digits=6) == round.(getindex.(data, 1), digits=6)
    @test round.(y0, digits=6) == round.(getindex.(data, 2), digits=6)

    ## [4] make sure EPSG Type is working
    @test typeof(EPSG(3031)) <: EPSG
end