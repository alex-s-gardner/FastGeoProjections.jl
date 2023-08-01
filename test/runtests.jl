using FastGeoProjections
using Test
using Proj

@testset "FastGeoProjections.jl" begin
    ## [1] Test WGS 84 / NSIDC Sea Ice Polar Stereographic North 
    lat0 = 84.0;
    lon0 = 50.0;

    trans = Proj.Transformation("EPSG:4326", "EPSG:3413")
    x0, y0 = trans(lat0,lon0)
    
    trans = FastGeoProjections.Transformation("EPSG:4326", "EPSG:3413")
    x1, y1 = trans(lat0, lon0)

    # check accuracy
    @test round(x0, digits=5) == round(x1, digits=5)
    @test round(y0, digits=5) == round(y1, digits=5)

    # now inverse 
    lat1, lon1 = inv(trans)(x1, y1)
    @test round(lat0, digits=8) == round(lat1, digits=8)
    @test round(lon0, digits=8) == round(lon1, digits=8)

    ## [2] Test WGS 84 / Antarctic Polar Stereographic
    lat0 = -84.0
    lon0 = 50.0

    trans = Proj.Transformation("EPSG:4326", "EPSG:3031")
    x0, y0 = trans(lat0, lon0)

    trans = FastGeoProjections.Transformation("EPSG:4326", "EPSG:3031")
    x1, y1 = trans(lat0, lon0)

    # check accuracy
    @test round(x0, digits=5) == round(x1, digits=5)
    @test round(y0, digits=5) == round(y1, digits=5)

    # now inverse 
    lat1, lon1 = inv(trans)(x1, y1)
    @test round(lat0, digits=8) == round(lat1, digits=8)
    @test round(lon0, digits=8) == round(lon1, digits=8)

    ## [3] Test Transverse Mercator projection [UTM] and vector input for North
    x0 = [10000., 20000.] 
    y0 = [10000., 20000.] 

    trans = Proj.Transformation("EPSG:32619", "EPSG:4326")
    ll0 = trans.(x0, y0)
    lat0 = [i[1] for i in ll0]
    lon0 = [i[2] for i in ll0]
    
    trans = FastGeoProjections.Transformation("EPSG:32619", "EPSG:4326")
    lat1, lon1 = trans(x0, y0)

    # check accuracy
    @test round.(lat0, digits=8) == round.(lat1, digits=8)
    @test round.(lon0, digits=8) == round.(lon1, digits=8)
    
    # now inverse 
    x1, y1 = inv(trans)(lat1, lon1)

    # check accuracy
    @test round.(x0, digits=5) == round.(x1, digits=5)
    @test round.(y0, digits=5) == round.(y1, digits=5)

    ## [4] Test Transverse Mercator projection [UTM] and vector input for South
    lat0 = -[80., 40., 1.]
    lon0 = [30., 31., 34.]

    trans = Proj.Transformation("EPSG:4326", "EPSG:32736")
    xy0 = trans.(lat0, lon0)
    x0 = [i[1] for i in xy0]
    y0 = [i[2] for i in xy0]

    trans = FastGeoProjections.Transformation(EPSG(4326), EPSG(32736))
    x1, y1 = trans(lat0, lon0)

    # check accuracy
    @test round.(x0, digits=5) == round.(x1, digits=5)
    @test round.(y0, digits=5) == round.(y1, digits=5)

    # now inverse 
    lat1, lon1 = inv(trans)(x1, y1)

    # check accuracy
    @test round.(lat0, digits=8) == round.(lat1, digits=8)
    @test round.(lon0, digits=8) == round.(lon1, digits=8)

    ## [4] make sure EPSG Type is working
    @test typeof(EPSG(3031)) <: EPSG
end