# The Proj-backed fallback: everything about it that does not mention Proj.
#
# Proj is a weak dependency, so the PJ objects this pools are opaque here -- the
# element type `P` is filled in by the extension in `ext/`, which is also where
# they are built. What lives here is the pooling, the calling convention and the
# traits, none of which need Proj to be loaded to be compiled.

"""
    ProjTransformation

Fallback point operator backed by Proj, used for any CRS pair FastGeoProjections
has no native implementation for. Built by [`pipeline`](@ref), and only once
`Proj` has been loaded -- see [`proj_transformation`](@ref).

A PJ object may not be shared across threads, so this holds a small pool of them
-- each with its own cloned context -- which a task checks out for the range it is
working on and returns afterwards. A task that finds the pool empty waits for one
to come back.

The pool is checked out rather than indexed by `threadid()`: a thread id is not a
stable identity for a task, which may migrate between yield points, and a thread
adopted after construction (a `@ccallable` entry from a foreign thread, say) has
an id past the end of any array sized when the object was built.
"""
mutable struct ProjTransformation{P} <: GeoTransformation
    source_epsg::EPSG
    target_epsg::EPSG
    always_xy::Bool
    pool::Channel{P}
    all::Vector{P}
end

"""
    proj_transformation(source_epsg, target_epsg, always_xy)

Build the [`ProjTransformation`](@ref) for a CRS pair. Implemented by the package
extension that `Proj` loads; without it there is no fallback, so this throws and
says what to import.

The error is raised where the transformation is built rather than where it is
first applied, so the message can name the CRS pair that has no native
implementation.

The fallback defined here is deliberately untyped in its arguments: the extension
adds the `(::EPSG, ::EPSG, ::Bool)` method, and a method it could *overwrite*
rather than take precedence over would make the extension fail to precompile.
"""
function proj_transformation(source_epsg, target_epsg, always_xy)
    throw(ArgumentError("""
        FastGeoProjections has no native transformation from \
        EPSG:$(first(source_epsg.val)) to EPSG:$(first(target_epsg.val)).
        Run `import Proj` to use the Proj.jl fallback for this CRS pair; the \
        natively implemented codes are in `FastGeoProjections.fast_epsg_codes`."""))
end

function borrow(f, t::ProjTransformation)
    p = take!(t.pool)
    try
        f(p.pj)
    finally
        put!(t.pool, p)
    end
end

# Calling one of these directly, outside a `borrow`, still has to be safe --
# it is the documented point-operator interface -- so it checks out a PJ of
# its own. That costs roughly 60ns over the ~70ns the transform itself takes;
# anything transforming more than a handful of points should be going through
# `transform`/`transform!`, which check out once per chunk instead.
function (t::ProjTransformation)(x, y)
    p = take!(t.pool)
    try
        p.pj(x, y)
    finally
        put!(t.pool, p)
    end
end

# the third coordinate goes to Proj, which applies the pipeline to it
function (t::ProjTransformation)(x, y, z)
    p = take!(t.pool)
    try
        p.pj(x, y, z)
    finally
        put!(t.pool, p)
    end
end

Base.inv(t::ProjTransformation) =
    proj_transformation(t.target_epsg, t.source_epsg, t.always_xy)

islanesafe(::ProjTransformation) = false

# The pipeline Proj resolves the EPSG pair into may be a datum shift, which is
# a rotation and translation in Cartesian space: the height participates in it,
# so a transformed x and y depend on the z given and the transformed z is not
# the one passed in. Both directions of that matter -- dropping the z silently
# transforms the point as though its height were zero, which moves x and y as
# well as losing the height.
preservesz(::ProjTransformation) = false

# `ncoords` stays at the default 2 all the same, which the native height-transforming
# operators do not: 2D in and 2D out is a thing Proj resolves, and is what such a
# pipeline is usually asked for, so it is not an error here the way it is for
# `LonLatToGeocentric`. This is why the two traits are separate.

Base.show(io::IO, t::ProjTransformation) =
    print(io, "ProjTransformation(EPSG:", first(t.source_epsg.val),
          " → EPSG:", first(t.target_epsg.val), ")")
