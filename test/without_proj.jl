# Proj is a weak dependency, so what the package does *before* it is loaded is
# part of the interface: every native projection has to work, and a CRS pair
# outside `fast_epsg_codes` has to say what to import.
#
# Run in a subprocess, because the rest of the suite loads Proj and a package
# extension cannot be unloaded once it is in.

@testset "without Proj loaded" begin
    script = """
    using FastGeoProjections
    using Test

    # Proj must not come in as a transitive dependency of anything else either:
    # if it does, the extension loads and this file tests nothing.
    @test !haskey(Base.loaded_modules, Base.PkgId(
        Base.UUID("c94c279d-25a6-4763-9509-64d165bea63e"), "Proj"))

    # every native pipeline works
    @test Transformation(EPSG(4326), EPSG(3413); always_xy = true)(-45.0, 70.0)[2] ≈
          -2.187927649279021e6
    @test Transformation(EPSG(4326), EPSG(32619); always_xy = true)(-69.0, 45.0)[1] > 0
    @test Transformation(EPSG(4979), EPSG(4978); always_xy = true)(5.39, 52.16, 100.0)[1] ≈
          3.9036404612786868e6
    # ...including over an array, which is where the operator meets `transform`
    @test length(transform(Transformation(EPSG(4326), EPSG(3857); always_xy = true),
                           [(1.0, 2.0), (3.0, 4.0)])) == 2

    # ...and a pair with no native implementation names both codes and says what
    # to import, rather than failing later or on something incidental
    err = try
        Transformation(EPSG(4326), EPSG(3395))
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    msg = sprint(showerror, err)
    @test occursin("EPSG:4326", msg) && occursin("EPSG:3395", msg)
    @test occursin("import Proj", msg)

    # `proj_only` has no fallback to force, so it says the same thing
    @test_throws "import Proj" Transformation(EPSG(4326), EPSG(3413); proj_only = true)
    """
    # `--startup-file=no`: a user's startup.jl may load Proj, which would make
    # the assertion above fail for a reason that is nothing to do with the package.
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(joinpath(@__DIR__, "..")) -e $script`
    @test success(pipeline(cmd; stdout, stderr))
end
