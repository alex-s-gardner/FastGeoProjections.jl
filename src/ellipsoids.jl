

"""
An ellipsoidal representation of the Earth [modified from Geodesy.jl]
"""
struct Ellipsoid
    a::Float64        # Semi-major axis
    b::Float64        # Semi-minor axis
    f::Float64        # Flattening
    e::Float64        # Eccentricity 
    name::Union{Nothing,Symbol}      # Conventional name - for clarity, should match the name
    epsg::EPSG         # epsg code
    # of the const instance in the package!
end

function Ellipsoid(; 
    a::Union{Nothing,Float64}=nothing, 
    b::Union{Nothing,Float64}=nothing, 
    f_inv::Union{Nothing,Float64}=nothing, 
    name::Union{Nothing,Symbol} = nothing, 
    epsg::Union{Nothing,EPSG} = nothing
    )

    if isnothing(a) || (isnothing(b) == isnothing(f_inv))
        throw(ArgumentError("Specify parameter 'a' and either 'b' or 'f_inv'"))
    end

    if isnothing(b)
        _ellipsoid_af(a, f_inv, name, epsg)
    else
        _ellipsoid_ab(a, b, name, epsg)
    end
end

function _ellipsoid_ab(a::Float64, b::Float64, name, epsg)
    f = 1 - b / a
    e = sqrt(f * (2 - f))
    Ellipsoid(a, b, f, e, name, epsg)
end
function _ellipsoid_af(a::Float64, f_inv::Float64, name, epsg)
    b = a * (1 - inv(f_inv))

    _ellipsoid_ab(a, b, name, epsg)
end

function Base.show(io::IO, el::Ellipsoid)
    if !isnothing(el.name)
        # To clarify that these are Ellipsoids, we wrap the name in
        # 'Ellipsoid', even though the name itself should resolve to the
        # correct ellipsoid instance.
        print(io, "Ellipsoid(name = $(el.name), epsg = $(el.epsg))")
    else
        print(io, "Ellipsoid(a=$(el.a), b=$(el.b))")
    end
end


"""
    ellipsoid(epsg::EPSG)
define an ellipsoid given an EPSG
"""
function ellipsoid(epsg::EPSG)
    if epsg.val == 7030
        ellips = Ellipsoid(; a = 6378137., f_inv = 298.257223563, name = :WGS_84, epsg = EPSG(7030))
    elseif epsg.val == 7019
        ellips = Ellipsoid(; a = 6378137., f_inv = 298.257222101, name = :GRS_1980, epsg = EPSG(7019))
    else
        error("$(epsg.val) ellisoid is not defined, you may need to add it to ellipsoids.jl")
    end
    return ellips
end
