# FastGeoProjections to Proj speed comparison
using FastGeoProjections
using Proj
using BenchmarkTools
using DataFrames
using GLMakie
using Printf

outfile = abspath("./assets/benchmark");
ns = [100, 1000, 10000, 100000, 1000000]
epsg_from_to = [
    (EPSG(4326), EPSG(3413)), 
    (EPSG(3031), EPSG(4326)),
    (EPSG(4326), EPSG(32636)), 
    (EPSG(32735), EPSG(4326))
]
system = "Apple M2 Max"
solutions = ["Proj: single-thread", "Proj:multi-thread", "FGP:single-thread", "FGP:multi-thread-64", "FGP:multi-thread-32"] #, "FGP:GPU-32"]

df = DataFrame();
for solution in solutions
    df[!, solution*"_time"] = zeros(length(epsg_from_to) * length(ns))
    df[!, solution*"_err"] =  zeros(length(epsg_from_to) * length(ns))
end

df[!, :npoints] =  zeros(length(epsg_from_to) * length(ns));
df[!, :epsg_from_to] .= [(EPSG(4326), EPSG(3413))];

ellips = ellipsoid(EPSG(7030));

for (i, n) in enumerate(ns)
    rando = rand(n);

    for k = eachindex(epsg_from_to)
         epsg_from = epsg_from_to[k][1]
         epsg_to = epsg_from_to[k][2]

        if k == 1
            Y = rando * 30 .+ 60;
            X = rando * 360 .- 180;
            XY = [(X[i], Y[i]) for i in eachindex(X)];
        elseif k == 2
            Y = -(rando * 30 .+ 60);
            X = rando * 360 .- 180;
            XY = FastGeoProjections.polarstereo_fwd(X, Y; ellips, lat_ts=-71.0, lon_0=0.0);
            XY = [(XY[1][i], XY[2][i]) for i in eachindex(XY[1])]
        elseif k == 3
            Y = rando * 80.
            X = rando * 9 .+ 28.5
            XY = [(X[i], Y[i]) for i in eachindex(X)];
        elseif k == 4
            Y = rando * -80.0
            X = rando * 9 .+ 22.5
            XY = FastGeoProjections.utm_fwd(X, Y; epsg=EPSG(epsg_from))
            XY = [(XY[1][i], XY[2][i]) for i in eachindex(XY[1])]
        end
        
        r = i+(k-1)*length(ns);
        df[r, :npoints] = n;
        df[r, :epsg_from_to] = epsg_from_to[k];

        printstyled("**EPSG:$(epsg_from.val) to EPSG:$(epsg_to.val) [n = $n]**\n", color=:blue)
        XY0 = epsg2epsg(XY, epsg_from, epsg_to; threaded=false, proj_only=true)

        printstyled("*Proj: single-thread*\n", color=:lightgrey)
        b = @benchmark XY1 = epsg2epsg($XY, $epsg_from, $epsg_to; threaded=false, proj_only=true)
        display(minimum(b))
        XY1 = epsg2epsg(XY, epsg_from, epsg_to; threaded=false, proj_only=true)
        err = maximum(abs.(vcat(XY1[1] - XY0[1], XY1[2] - XY0[2])))
        printstyled("\n  MAXIMUM ERROR:\t\t$err\n\n", color=:lightgrey)
        df[r, solutions[1]*"_time"] = minimum(b).time
        df[r, solutions[1]*"_err"] = err

        printstyled("*Proj: multi-thread*\n", color=:lightgrey)
        b = @benchmark XY1 = epsg2epsg($XY, $epsg_from, $epsg_to; threaded=true, proj_only=true)
        display(minimum(b))
        XY1 = epsg2epsg(XY, epsg_from, epsg_to; threaded=true, proj_only=true)
        err = maximum(abs.(vcat(XY1[1] - XY0[1], XY1[2] - XY0[2])))
        printstyled("\n  MAXIMUM ERROR:\t\t$err\n\n", color=:lightgrey)
        df[r, solutions[2]*"_time"] = minimum(b).time
        df[r, solutions[2]*"_err"] = err

        printstyled("*FastGeoProjections: single-thread*\n", color=:lightgrey)
        b = @benchmark XY1 = epsg2epsg($XY, $epsg_from, $epsg_to; threaded=false)
        display(minimum(b))
        XY1 = epsg2epsg(XY, epsg_from, epsg_to; threaded=false)
        err = maximum(abs.(vcat(XY1[1] - XY0[1], XY1[2] - XY0[2])))
        printstyled("\n  MAXIMUM ERROR:\t\t$err\n\n", color=:lightgrey)
        df[r, solutions[3]*"_time"] = minimum(b).time
        df[r, solutions[3]*"_err"] = err

        printstyled("*FastGeoProjections: multi-thread - Float64*\n", color=:lightgrey)
        b = @benchmark XY1 = epsg2epsg($XY, $epsg_from, $epsg_to)
        display(minimum(b))
        XY1 = epsg2epsg(XY, epsg_from, epsg_to; threaded=true)
        err = maximum(abs.(vcat(XY1[1] - XY0[1], XY1[2] - XY0[2])))
        printstyled("\n  MAXIMUM ERROR:\t\t$err\n\n", color=:lightgrey)
        df[r, solutions[4]*"_time"] = minimum(b).time
        df[r, solutions[4]*"_err"] = err

        printstyled("*FastGeoProjections: multi-thread - Float32*\n", color=:lightgrey)
        XY = [Float32.(foo) for foo in XY]
        b = @benchmark XY1 = epsg2epsg($XY, $epsg_from, $epsg_to; threaded=true)
        display(minimum(b))
        XY1 = epsg2epsg(XY, epsg_from, epsg_to; threaded=true)
        err = maximum(abs.(vcat(XY1[1] - XY0[1], XY1[2] - XY0[2])))
        printstyled("\n  MAXIMUM ERROR:\t\t$err\n\n", color=:lightgrey)
        df[r, solutions[5]*"_time"] = minimum(b).time
        df[r, solutions[5]*"_err"] = err
    end
end


Makie.inline!(false)
f = Figure(resolution=(1500, 750 * ceil(length(epsg_from_to)/2)), fontsize = 25)

for (i, epsg) in enumerate(epsg_from_to)
   
    r = ceil(Int64, i / 2)
    c = i - 2*(r-1)

    ax = Axis(
        f[c, r],
        yscale=log10,
        xscale=log10,
        title="EPSG:$(epsg[1].val) => EPSG:$(epsg[2].val)",
        yminorticksvisible=true,
        yminorgridvisible=true,
        xlabel="number of points converted",
        ylabel="compute time [µs]",
        yminorticks=IntervalsBetween(5),
    )

    rs = df.epsg_from_to .== [epsg]
    re = (df.epsg_from_to .== [epsg]) .& (df[:, :npoints] .== maximum(ns))
    for  solution in solutions

        err = df[re, solution*"_err"]
        err = round(err[1], sigdigits=3, base=10)
        lines!(df[rs, :npoints], df[rs, solution*"_time"] ./ 1000, label= "$solution [ME=$err]", linewidth = 3)
    end

    axislegend(ax, framevisible=false, position=:lt)
    f
end

save(abspath("$outfile.jpg"), f)

