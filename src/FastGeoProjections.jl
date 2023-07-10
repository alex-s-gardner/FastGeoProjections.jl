module FastGeoProjections
    using Proj # Proj dependancy included untill package is more mature
    using GeoFormatTypes
    using LoopVectorization

    include("ellipsoids.jl")
    include("polarstereo.jl")
    include("tranmerc.jl")
    include("utm_ups.jl")
    include("epsg2epsg.jl")

    export epsg2epsg
    export EPSG
    export ellipsoid

    precompile(tranmerc_fwd, (Vector{Float64}, Vector{Float64},))
    precompile(tranmerc_fwd, (Vector{Float32}, Vector{Float32},))
    precompile(tranmerc_inv, (Vector{Float64}, Vector{Float64},))
    precompile(tranmerc_inv, (Vector{Float32}, Vector{Float32},))
    precompile(polarstereo_fwd, (Vector{Float64}, Vector{Float64},))
    precompile(polarstereo_fwd, (Vector{Float32}, Vector{Float32},))
    precompile(polarstereo_inv, (Vector{Float64}, Vector{Float64},))
    precompile(polarstereo_inv, (Vector{Float32}, Vector{Float32},))
end