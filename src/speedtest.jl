using FastGeoProjections
using Test
using Proj

n = 10000000;
r = rand(n)

# scale to range
latitude = r .* 30 .+ 60;
longitude = r .* 360 .- 180;

@time polarstereo_fwd.(longitude, latitude; a=6378137.0, e=0.08181919, lat_ts=70.0, lon_0=-45.0)

x = ones(size(latitude))
y = ones(size(latitude))
@time Threads.@threads for i in eachindex(latitude)
    x[i], y[i] = polarstereo_fwd.(longitude[i], latitude[i]; a=6378137.0, e=0.08181919, lat_ts=70.0, lon_0=-45.0)
end

@time begin
    # build transformation 
    trans = Proj.Transformation("EPSG:4326", "EPSG:3413", always_xy=true)

    # project points
    data = trans.(longitude, latitude)
end
