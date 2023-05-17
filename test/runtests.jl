using FastGeoProjections
using Test
using Proj

@testset "FastGeoProjections.jl" begin
    # Write your tests here.


    # [1] Test WGS 84 / NSIDC Sea Ice Polar Stereographic North 
    latitude = 84;
    longitude = 50;

    x, y = polarstereo_fwd(longitude, latitude; a=6378137.0, e=0.08181919, lat_ts=70, lon_0=-45)
    
    # build transformation 
    trans = Proj.Transformation("EPSG:4326", "EPSG:3413", always_xy=true)

    # project points
    data = trans.(longitude, latitude)

    # check accuracy within 1 mm
    @test round(x, digits=3) == round(data[1], digits = 3)
    @test round(y, digits=3) == round(data[2], digits=3)

    # now inverse 
    longitude0, latitude0 = polarstereo_inv(x, y; a=6378137.0, e=0.08181919, lat_ts=70, lon_0=-45)
    @test round(longitude0, digits=6) == longitude
    @test round(latitude0, digits=6) == latitude

    # [1] Test WGS 84 / Antarctic Polar Stereographi
    latitude = -84
    longitude = 50

    x, y = polarstereo_fwd(longitude, latitude; a=6378137.0, e=0.08181919, lat_ts=-71, lon_0= 0)

    # build transformation 
    trans = Proj.Transformation("EPSG:4326", "EPSG:3031", always_xy=true)

    # project points
    data = trans.(longitude, latitude)

    # check accuracy within 1 mm
    @test round(x, digits=3) == round(data[1], digits=3)
    @test round(y, digits=3) == round(data[2], digits=3)

    # now inverse 
    longitude0, latitude0 = polarstereo_inv(x, y; a=6378137.0, e=0.08181919, lat_ts=-71, lon_0=0)
    @test round(longitude0, digits=6) == longitude
    @test round(latitude0, digits=6) == latitude
end
