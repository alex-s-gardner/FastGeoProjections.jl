module FastGeoProjections
    using Proj # Proj dependancy included untill package is more mature
    using GeoFormatTypes

    include("polarstereo.jl")
    include("epsg2epsg.jl")
    export epsg2epsg
    export EPSG
end
