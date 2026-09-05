#!/usr/bin/env julia
#
# Build the `fastgeoproj` executable.
#
#     julia --project=examples/juliac examples/juliac/build.jl
#
# Options:
#     --trim=safe          fail the build on any unresolved call (the default)
#     --trim=unsafe-warn   report them and build anyway
#     --output-exe PATH    where to put the binary
#
# Needs Julia 1.12 or later, and a C compiler on PATH for the final link.
#
# ---------------------------------------------------------------------------
# Why LocalPreferences.toml sets HostCPUFeatures.freeze_cpu_target
#
# Without it this program does not trim. `--trim=safe` fails with two
# unresolved calls:
#
#     [1] feature_string()   HostCPUFeatures/src/cpu_info.jl:3   <- dlopen(libLLVM)
#     [2] reset_features!()  HostCPUFeatures/src/cpu_info.jl:48
#     [3] redefine()         HostCPUFeatures/src/HostCPUFeatures.jl:62
#     [4] __init__()         HostCPUFeatures/src/HostCPUFeatures.jl:95
#
# `HostCPUFeatures.__init__` calls `redefine()` when the CPU name it saw at
# precompile time differs from the one it sees at run time, and that reaches a
# `dlopen` of libLLVM to read `LLVMGetHostCPUFeatures`. The verifier cannot
# prove the branch is dead, so it follows it. `--trim=unsafe-warn` reports the
# two and links anyway, which is easy to mistake for success: any verifier
# error means the program did not trim.
#
# Nothing here uses HostCPUFeatures. It arrives through VectorizationBase --
# the only package in the manifest that depends on it -- which FastGeoProjections
# uses for `vload`/`vstore!`/`stridedpointer`/`MM`/`mask`, and which
# SLEEFPirates depends on in turn. Its one job in that chain is to work out the
# vector register width at load time.
#
# `__init__` returns early when `freeze_cpu_target` is set, and that is a
# `const` read from a Preference, so setting it makes the whole path
# statically dead:
#
#     [HostCPUFeatures]
#     freeze_cpu_target = true
#
# For an ahead-of-time build this is the honest setting rather than a dodge.
# juliac already compiles with `-C native`, so the vector width is baked into
# the emitted lane loops; a binary that discovered a wider register at startup
# could not use it. HostCPUFeatures is a direct dependency of this project for
# no other reason: Preferences are only applied to direct dependencies, and as
# a transitive one the setting is read and ignored.
#
# Measured on aarch64: `pick_vector_width` is unchanged at 2 x Float64 /
# 4 x Float32, output is byte-identical, and 1e6 points still take 0.17 s. On
# x86 the setting freezes an *under-approximation* of the CPU features, so
# check `VectorizationBase.pick_vector_width` there before trusting it.
#
# Removing the dependency outright would mean replacing VectorizationBase's
# lane primitives and SLEEFPirates' kernels both -- the whole vectorization
# layer -- which is a much larger question than this example.
# ---------------------------------------------------------------------------

using Pkg
Pkg.instantiate()

const HERE = @__DIR__
const JULIAC = joinpath(Sys.BINDIR, Base.DATAROOTDIR, "julia", "juliac", "juliac.jl")

isfile(JULIAC) ||
    error("""
          juliac.jl is not at $JULIAC.
          Ahead-of-time compilation needs Julia 1.12 or later; this is $VERSION.
          """)

trim = "safe"
out = joinpath(HERE, "fastgeoproj")
i = 1
while i <= length(ARGS)
    a = ARGS[i]
    if startswith(a, "--trim=")
        global trim = a[(length("--trim=") + 1):end]
    elseif a == "--output-exe"
        i += 1
        i <= length(ARGS) || error("--output-exe needs a value")
        global out = ARGS[i]
    else
        error("unknown argument: $a")
    end
    global i = i + 1
end

cmd = `$(Base.julia_cmd()) --startup-file=no --history-file=no --project=$HERE
       $JULIAC --output-exe $out --experimental --trim=$trim
       $(joinpath(HERE, "fastgeoproj.jl"))`

@info "building" out trim
run(cmd)
@info "built" out size = Base.format_bytes(filesize(out))
