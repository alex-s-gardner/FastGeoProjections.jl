#!/usr/bin/env julia
#
# Build the `fastgeoproj` executable.
#
#     julia --project=examples/juliac examples/juliac/build.jl
#
# Options:
#     --trim=safe          fail the build on any unresolved call (see README)
#     --trim=unsafe-warn   report them and build anyway (the default)
#     --output-exe PATH    where to put the binary
#
# Needs Julia 1.12 or later, and a C compiler on PATH for the final link.

using Pkg
Pkg.instantiate()

const HERE = @__DIR__
const JULIAC = joinpath(Sys.BINDIR, Base.DATAROOTDIR, "julia", "juliac", "juliac.jl")

isfile(JULIAC) ||
    error("""
          juliac.jl is not at $JULIAC.
          Ahead-of-time compilation needs Julia 1.12 or later; this is $VERSION.
          """)

trim = "unsafe-warn"
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
