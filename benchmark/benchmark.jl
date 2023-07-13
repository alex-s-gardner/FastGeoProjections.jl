# FastGeoProjections to Proj speed comparison
using FastGeoProjections
using BenchmarkTools
using DataFrames
using GLMakie

outfile = abspath("./benchmark/benchmark");
ns = [100, 1000, 10000, 100000, 1000000]
epsg_target_source = [
    (EPSG(4326), EPSG(3413)), 
    (EPSG(3031), EPSG(4326)),
    (EPSG(4326), EPSG(32636)), 
    (EPSG(32735), EPSG(4326))
]
system = "Apple M2 Max"
threads = Threads.nthreads()

solutions = ["Proj: single-thread", "Proj: multi-thread", "FGP: single-thread", "FGP: multi-thread"] 
# Float32 and GPU yielded little benefit 

df = DataFrame();
for solution in solutions
    df[!, solution*"_time"] = zeros(length(epsg_target_source) * length(ns))
    df[!, solution*"_err"] =  zeros(length(epsg_target_source) * length(ns))
end

df[!, :npoints] =  zeros(length(epsg_target_source) * length(ns));
df[!, :epsg_target_source] .= [(EPSG(4326), EPSG(3413))];

function trans_bench(X, Y, X0, Y0, df, r, solution, source_epsg, target_epsg, threaded, proj_only, always_xy)
    b = @benchmark begin
        trans = FastGeoProjections.Transformation($source_epsg, $target_epsg; threaded=$threaded, proj_only=$proj_only, always_xy=$always_xy)
        X1, Y1 = trans($X, $Y)
    end

    trans = FastGeoProjections.Transformation(source_epsg, target_epsg; threaded=threaded, proj_only=proj_only, always_xy=always_xy)
    
    display(minimum(b))
    X1, Y1 = trans(X, Y)
    err = maximum(abs.([X1 - X0; Y1 - Y0]))
    printstyled("\n  MAXIMUM ERROR:\t\t$err\n\n", color=:lightgrey)
    df[r, solution*"_time"] = minimum(b).time
    df[r, solution*"_err"] = err
    return df
end


for (i, n) in enumerate(ns)
    rando = rand(n);

    for k = eachindex(epsg_target_source)
         source_epsg = epsg_target_source[k][1]
         target_epsg = epsg_target_source[k][2]

        if k == 1
            Y = rando * 30 .+ 60;
            X = rando * 360 .- 180;
        elseif k == 2
            Y = -(rando * 30 .+ 60);
            X = rando * 360 .- 180;
            X, Y = FastGeoProjections.polarstereo_fwd(X, Y; lat_ts=-71.0, lon_0=0.0);
        elseif k == 3
            Y = rando * 80.
            X = rando * 9 .+ 28.5
        elseif k == 4
            Y = rando * -80.0
            X = rando * 9 .+ 22.5
            X, Y = FastGeoProjections.utm_fwd(X, Y; epsg=EPSG(source_epsg))
        end
        
        r = i+(k-1)*length(ns);
        df[r, :npoints] = n;
        df[r, :epsg_target_source] = epsg_target_source[k];

        printstyled("**EPSG:$(source_epsg.val) to EPSG:$(target_epsg.val) [n = $n]**\n", color=:blue)
        trans = FastGeoProjections.Transformation(source_epsg, target_epsg; threaded=false, proj_only=true, always_xy=true)
        X0, Y0 = trans(X, Y)

        printstyled("*Proj: single-thread*\n", color=:lightgrey)
        threaded = false; proj_only = true; always_xy = true
        df = trans_bench(X, Y, X0, Y0, df, r, solutions[1], source_epsg, target_epsg, threaded, proj_only, always_xy)
       
        printstyled("*Proj: multi-thread*\n", color=:lightgrey)
        threaded = true; proj_only = true; always_xy = true
        df = trans_bench(X, Y, X0, Y0, df, r, solutions[2], source_epsg, target_epsg, threaded, proj_only, always_xy)


        printstyled("*FastGeoProjections: single-thread*\n", color=:lightgrey)
        threaded = false; proj_only = false; always_xy = true
        df = trans_bench(X, Y, X0, Y0, df, r, solutions[3], source_epsg, target_epsg, threaded, proj_only, always_xy)


        printstyled("*FastGeoProjections: multi-thread - Float64*\n", color=:lightgrey)
        threaded = true; proj_only = false; always_xy = true
        df = trans_bench(X, Y, X0, Y0, df, r, solutions[4], source_epsg, target_epsg, threaded, proj_only, always_xy)        
    end
end


Makie.inline!(false)
f = Figure(resolution=(1500, 750 * ceil(length(epsg_target_source)/2)), fontsize = 25)
col = Makie.wong_colors();
for (i, epsg) in enumerate(epsg_target_source)
   
    r = ceil(Int64, i / 2)
    c = i - 2*(r-1)

    ax = Axis(
        f[r, c],
        yscale=log10,
        xscale=log10,
        title="EPSG:$(epsg[1].val) => EPSG:$(epsg[2].val)",
        yminorticksvisible=true,
        yminorgridvisible=true,
        xlabel="points converted",
        ylabel="compute time [µs]",
        yminorticks=IntervalsBetween(5),
    )

    rs = df.epsg_target_source .== [epsg]
    re = (df.epsg_target_source .== [epsg]) .& (df[:, :npoints] .== maximum(ns))

    lins = [lines!(df[rs, :npoints], df[rs, solution*"_time"] ./ 1000, label="$solution", linewidth=6, color=col[i]) for (i, solution) in enumerate(solutions)]

    if i == 1
        legends = axislegend(ax, lins, solutions, position=:lt)
    end

end

supertitle = Label(f[0, :], "FastGeoProjections.jl benchmarks,  $system using $threads threads", fontsize=30)


save(abspath("$outfile.jpg"), f)

