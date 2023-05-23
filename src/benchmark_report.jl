# FastGeoProjections to Proj speed comparison
using FastGeoProjections
using Proj
using BenchmarkTools
using Metal

outfile = abspath("./assets/benchmark.txt");

open(abspath("./assets/benchmark.txt"), "w") do io
    for n in [1000, 1000000, 10000000]
        r = rand(n);

        for k = 1:2
            if k == 1
                epsg_from = EPSG(4326)
                epsg_to = EPSG(3413)
                Y = r * 30 .+ 60;
                X = r * 360 .- 180;
                XY = [(X[i], Y[i]) for i in eachindex(X)];
            elseif k == 2
                epsg_from = EPSG(3031)
                epsg_to = EPSG(4326)

                Y = -(r * 30 .+ 60);
                X = r * 360 .- 180;
                XY = FastGeoProjections.polarstereo_fwd.(X, Y; a=6378137.0, e=0.08181919, lat_ts=-71.0, lon_0=0);

            elseif k == 3
                epsg_from = EPSG(4326)
                epsg_to = EPSG(32609)

                Y =  r * 30 .+ 30
                X = -(r * 14 .+ 139)
                XY = [(X[i], Y[i]) for i in eachindex(X)]
            end
            
            printstyled(io, "**EPSG:$(epsg_from.val) to EPSG:$(epsg_to.val) [n = $n]**\n", color=:blue)
            XY0 = epsg2epsg(XY, epsg_from, epsg_to; threaded=false, proj_only=true)
            println(io, "")

            printstyled(io, "*Proj: single-thread*\n", color=:lightgrey)
            b = @benchmark XY1 = epsg2epsg($XY, $epsg_from, $epsg_to; threaded=false, proj_only=true)
            XY1 = epsg2epsg(XY, epsg_from, epsg_to; threaded=false, proj_only=true)
            err = Float16(maximum(abs.(vcat(XY1[1] - XY0[1], XY1[2] - XY0[2]))))
            show(io, MIME"text/plain"(), minimum(b))
            printstyled(io, "\n  MAXIMUM ERROR:\t\t$err\n", color=:lightgrey)
            println(io, "")

            printstyled(io, "*Proj: multi-thread*\n", color=:lightgrey)
            b = @benchmark XY1 = epsg2epsg($XY, $epsg_from, $epsg_to; threaded=true, proj_only=true)
            XY1 = epsg2epsg(XY, epsg_from, epsg_to; threaded=true, proj_only=true)
            err = Float16(maximum(abs.(vcat(XY1[1] - XY0[1], XY1[2] - XY0[2]))))
            show(io, MIME"text/plain"(), minimum(b))
            printstyled(io, "\n  MAXIMUM ERROR:\t\t$err\n", color=:lightgrey)
            println(io, "")

            printstyled(io, "*FastGeoProjections: single-thread*\n", color=:lightgrey)
            b = @benchmark XY1 = epsg2epsg($XY, $epsg_from, $epsg_to; threaded=false)
            XY1 = epsg2epsg(XY, epsg_from, epsg_to; threaded=false)
            err = Float16(maximum(abs.(vcat(XY1[1] - XY0[1], XY1[2] - XY0[2]))))
            show(io, MIME"text/plain"(), minimum(b))
            printstyled(io, "\n  MAXIMUM ERROR:\t\t$err\n", color=:lightgrey)
            println(io, "")

            printstyled(io, "*FastGeoProjections: multi-thread - Float64*\n", color=:lightgrey)
            b = @benchmark XY1 = epsg2epsg($XY, $epsg_from, $epsg_to; threaded=true)
            XY1 = epsg2epsg(XY, epsg_from, epsg_to; threaded=true)
            err = Float16(maximum(abs.(vcat(XY1[1] - XY0[1], XY1[2] - XY0[2]))))
            show(io, MIME"text/plain"(), minimum(b))
            printstyled(io, "\n  MAXIMUM ERROR:\t\t$err\n", color=:lightgrey)
            println(io, "")

            printstyled(io, "*FastGeoProjections: multi-thread - Float32*\n", color=:lightgrey)
            XY = [Float32.(foo) for foo in XY]
            b = @benchmark XY1 = epsg2epsg($XY, $epsg_from, $epsg_to; threaded=true)
            XY1 = epsg2epsg(XY, epsg_from, epsg_to; threaded=true)
            err = Float16(maximum(abs.(vcat(XY1[1] - XY0[1], XY1[2] - XY0[2]))))
            show(io, MIME"text/plain"(), minimum(b))
            printstyled(io, "\n  MAXIMUM ERROR:\t\t$err\n", color=:lightgrey)
            println(io, "")

            printstyled(io, "*FastGeoProjections: M2 GPU - Float32 [including transfer time]*\n", color=:lightgrey)
            b = @benchmark begin
                XY1 = epsg2epsg(MtlArray($XY), $epsg_from, $epsg_to; threaded=false)
                XY1 = (Array(getindex(XY1,1)), Array(getindex(XY1,2)))
            end
            XY1 = epsg2epsg(MtlArray(XY), epsg_from, epsg_to; threaded=false)
            XY1 = (Array(getindex(XY1, 1)), Array(getindex(XY1, 2)))
            err = Float16(maximum(abs.(vcat(XY1[1] - XY0[1], XY1[2] - XY0[2]))))
            show(io, MIME"text/plain"(), minimum(b))
            printstyled(io, "\n  MAXIMUM ERROR:\t\t$err\n", color=:lightgrey)
            println(io, "")
            println(io, "")
        end
    end
end