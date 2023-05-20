# FastGeoProjections to Proj speed comparison
using FastGeoProjections
using Proj
using BenchmarkTools

n = 1000000;
r = rand(n)

for k = 1:3
    if k == 1
        epsg_from = EPSG(4326)
        epsg_to = EPSG(3413)
        Y = r .* 30 .+ 60;
        X = r .* 360 .- 180;
        XY = [(X[i], Y[i]) for i in eachindex(X)]
    elseif k == 2
        epsg_from = EPSG(3031)
        epsg_to = EPSG(4326)

        Y = -(r .* 30 .+ 60);
        X = r .* 360 .- 180;
        XY = FastGeoProjections.polarstereo_fwd.(X, Y; a=6378137.0, e=0.08181919, lat_ts=-71.0, lon_0=0);

    elseif k == 3
        epsg_from = EPSG(4326)
        epsg_to = EPSG(32609)

        Y =  r .* 30 .+ 30
        X = -(r .* 14 .+ 139)
        XY = [(X[i], Y[i]) for i in eachindex(X)]
    end
    
    printstyled("EPSG:$(epsg_from.val) to EPSG:$(epsg_to.val)\n", color=:blue)
   
    printstyled("FastGeoProjections: single-thread\n", color=:lightgrey)
    @btime XY2 = epsg2epsg($XY, $epsg_from, $epsg_to; threaded=false)

    printstyled("FastGeoProjections: multi-thread\n", color=:lightgrey)
    @btime XY2 = epsg2epsg($XY, $epsg_from, $epsg_to; threaded=true)

    printstyled("Proj: single-thread\n", color=:lightgrey)
    @btime XY2 = epsg2epsg($XY, $epsg_from, $epsg_to; threaded=false, proj_only = true)

    printstyled("Proj: multi-thread\n", color=:lightgrey)
    @btime XY2 = epsg2epsg($XY, $epsg_from, $epsg_to; threaded=true, proj_only=true)
end