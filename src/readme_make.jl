
# create a ReadMe markdown 
open(abspath("./README.md"), "w") do io

    # FastGeoProjections
    println(io,"[![Build Status](https://github.com/alex-s-gardner/FastGeoProjections.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/alex-s-gardner/FastGeoProjections.jl/actions/workflows/CI.yml?query=branch%3Amain)\n")
    println(io,"!! UNDER DEVELOPMENT BY OVERCOMMITTED AND UNDER PAID DEVELOPERS !!\n")
    println(io,"**FastGeoProjections** is intended to provide highly optimized native Julia geospatial coordinate transformations from one coordinate reference system (CRS) to another as defined by EPSG codes. It is not intended to to replace or to be as comprehensive as [Proj](https://github.com/JuliaGeo/Proj.jl). The package will natively support only the most common geospatial transformations and relies on **Proj.jl** for all others.\n")

    println(io, "Benchmark of currently implemented EPSGs\n")

    open(abspath("./assets/benchmark.txt"), "r") do f
        # read till end of file
        while !eof(f)
            # read a new / next line for every iteration          
            s = readline(f)
            println(io, "$s")
        end
    end
end

