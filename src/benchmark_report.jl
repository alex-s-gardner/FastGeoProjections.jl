# FastGeoProjections to Proj speed comparison
using FastGeoProjections
using Proj
using BenchmarkTools
using Metal
using DataFrames
using GLMakie
using Printf

outfile = abspath("./assets/benchmark");
ns = [1000, 10000, 100000, 1000000, 10000000]
epsg_from_to = [
    (EPSG(4326), EPSG(3413)), 
    (EPSG(3031), EPSG(4326))
]
system = "Apple M2 Max"
solutions = ["Proj: single-thread", "Proj:multi-thread", "FGP:single-thread", "FGP:multi-thread-64", "FGP:multi-thread-32", "FGP:GPU-32"]

df = DataFrame();
for solution in solutions
    df[!, solution*"_time"] = zeros(length(epsg_from_to) * length(ns))
    df[!, solution*"_err"] =  zeros(length(epsg_from_to) * length(ns))
end


df[!, :npoints] =  zeros(length(epsg_from_to) * length(ns));
df[!, :epsg_from_to] .= [(EPSG(4326), EPSG(3413))];

open(abspath("$outfile.txt"), "w") do io
    for (i, n) in enumerate(ns)
        rando = rand(n);

        for k = eachindex(epsg_from_to)
            if k == 1
                epsg_from = epsg_from_to[k][1]
                epsg_to = epsg_from_to[k][2]
                Y = rando * 30 .+ 60;
                X = rando * 360 .- 180;
                XY = [(X[i], Y[i]) for i in eachindex(X)];
            elseif k == 2
                epsg_from = epsg_from_to[k][1]
                epsg_to = epsg_from_to[k][2]

                Y = -(rando * 30 .+ 60);
                X = rando * 360 .- 180;
                XY = FastGeoProjections.polarstereo_fwd.(X, Y; a=6378137.0, e=0.08181919, lat_ts=-71.0, lon_0=0.);
            end
            
            r = i+(k-1)*length(ns)
            df[r, :npoints] = n
            df[r, :epsg_from_to] = epsg_from_to[k]

            printstyled(io, "**EPSG:$(epsg_from.val) to EPSG:$(epsg_to.val) [n = $n]**\n", color=:blue)
            XY0 = epsg2epsg(XY, epsg_from, epsg_to; threaded=false, proj_only=true)
            println(io, "")

            printstyled(io, "*Proj: single-thread*\n", color=:lightgrey)
            b = @benchmark XY1 = epsg2epsg($XY, $epsg_from, $epsg_to; threaded=false, proj_only=true)
            XY1 = epsg2epsg(XY, epsg_from, epsg_to; threaded=false, proj_only=true)
            err = maximum(abs.(vcat(XY1[1] - XY0[1], XY1[2] - XY0[2])))
            show(io, MIME"text/plain"(), minimum(b))
            printstyled(io, "\n  MAXIMUM ERROR:\t\t$err\n", color=:lightgrey)
            println(io, "")
            df[r, solutions[1]*"_time"] = minimum(b).time
            df[r, solutions[1]*"_err"] = err

            printstyled(io, "*Proj: multi-thread*\n", color=:lightgrey)
            b = @benchmark XY1 = epsg2epsg($XY, $epsg_from, $epsg_to; threaded=true, proj_only=true)
            XY1 = epsg2epsg(XY, epsg_from, epsg_to; threaded=true, proj_only=true)
            err = maximum(abs.(vcat(XY1[1] - XY0[1], XY1[2] - XY0[2])))
            show(io, MIME"text/plain"(), minimum(b))
            printstyled(io, "\n  MAXIMUM ERROR:\t\t$err\n", color=:lightgrey)
            println(io, "")
            df[r, solutions[2]*"_time"] = minimum(b).time
            df[r, solutions[2]*"_err"] = err

            printstyled(io, "*FastGeoProjections: single-thread*\n", color=:lightgrey)
            b = @benchmark XY1 = epsg2epsg($XY, $epsg_from, $epsg_to; threaded=false)
            XY1 = epsg2epsg(XY, epsg_from, epsg_to; threaded=false)
            err = maximum(abs.(vcat(XY1[1] - XY0[1], XY1[2] - XY0[2])))
            show(io, MIME"text/plain"(), minimum(b))
            printstyled(io, "\n  MAXIMUM ERROR:\t\t$err\n", color=:lightgrey)
            println(io, "")
            df[r, solutions[3]*"_time"] = minimum(b).time
            df[r, solutions[3]*"_err"] = err

            printstyled(io, "*FastGeoProjections: multi-thread - Float64*\n", color=:lightgrey)
            b = @benchmark XY1 = epsg2epsg($XY, $epsg_from, $epsg_to; threaded=true)
            XY1 = epsg2epsg(XY, epsg_from, epsg_to; threaded=true)
            err = maximum(abs.(vcat(XY1[1] - XY0[1], XY1[2] - XY0[2])))
            show(io, MIME"text/plain"(), minimum(b))
            printstyled(io, "\n  MAXIMUM ERROR:\t\t$err\n", color=:lightgrey)
            println(io, "")
            df[r, solutions[4]*"_time"] = minimum(b).time
            df[r, solutions[4]*"_err"] = err

            printstyled(io, "*FastGeoProjections: multi-thread - Float32*\n", color=:lightgrey)
            XY = [Float32.(foo) for foo in XY]
            b = @benchmark XY1 = epsg2epsg($XY, $epsg_from, $epsg_to; threaded=true)
            XY1 = epsg2epsg(XY, epsg_from, epsg_to; threaded=true)
            err = maximum(abs.(vcat(XY1[1] - XY0[1], XY1[2] - XY0[2])))
            show(io, MIME"text/plain"(), minimum(b))
            printstyled(io, "\n  MAXIMUM ERROR:\t\t$err\n", color=:lightgrey)
            println(io, "")
            df[r, solutions[5]*"_time"] = minimum(b).time
            df[r, solutions[5]*"_err"] = err

            printstyled(io, "*FastGeoProjections: M2 GPU - Float32 [including transfer time]*\n", color=:lightgrey)
            b = @benchmark begin
                XY1 = epsg2epsg(MtlArray($XY), $epsg_from, $epsg_to; threaded=false)
                XY1 = (Array(getindex(XY1,1)), Array(getindex(XY1,2)))
            end
            XY1 = epsg2epsg(MtlArray(XY), epsg_from, epsg_to; threaded=false)
            XY1 = (Array(getindex(XY1, 1)), Array(getindex(XY1, 2)))
            err = maximum(abs.(vcat(XY1[1] - XY0[1], XY1[2] - XY0[2])))
            show(io, MIME"text/plain"(), minimum(b))
            printstyled(io, "\n  MAXIMUM ERROR:\t\t$err\n", color=:lightgrey)
            println(io, "")
            println(io, "")
            df[r, solutions[6]*"_time"] = minimum(b).time
            df[r, solutions[6]*"_err"] = err
        end
    end
end


Makie.inline!(false)
f = Figure(resolution=(1000, 750 * length(epsg_from_to)), fontsize = 25)

for (i, epsg) in enumerate(epsg_from_to)
   

    ax = Axis(
        f[i, 1],
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

    axislegend(ax, framevisible=false, position=:rb)
    f
end

save(abspath("$outfile.png"), f)

