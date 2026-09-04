module FastGeoProjections
    using Proj # Proj dependancy included untill package is more mature
    using GeoFormatTypes
    using CoordinateTransformations
    using VectorizationBase
    using VectorizationBase: vload, vstore!, stridedpointer, MM
    import GeoInterface as GI

    include("kernels.jl")
    using .Math: Math, MathKernel, FastKernel, SLEEFKernel, BaseKernel,
                 DEFAULT_KERNEL, vectorizes

    include("ellipsoids.jl")
    include("transformations.jl")
    include("apply.jl")
    include("projections/direction.jl")
    include("projections/polarstereo.jl")
    include("projections/tranmerc.jl")
    include("projections/utm.jl")
    include("projections/geocentric.jl")
    include("projections/fuse.jl")
    include("proj.jl")
    include("epsg.jl")
    include("coord.jl")

    export EPSG
    export Transformation, transform, transform!, inv
    export preservesz
    export GeoTransformation
    export LonLatToPolarStereographic, PolarStereographicToLonLat
    export LonLatToTransverseMercator, TransverseMercatorToLonLat
    export LonLatToUTM, UTMToLonLat, convergence_scale
    export LonLatToGeocentric, GeocentricToLonLat
    export Math, MathKernel, BaseKernel, SLEEFKernel, FastKernel
end
