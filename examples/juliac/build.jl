#!/usr/bin/env julia
#
# Build the `fastgeoproj` executable.
#
#     julia --project=examples/juliac examples/juliac/build.jl
#
# Options:
#     --trim=safe          fail the build on any unresolved call (the default)
#     --trim=unsafe-warn   report them and build anyway
#     --output-exe NAME    executable name, written next to this script
#
# Needs Julia 1.13 or later, and a C compiler on PATH for the final link.
#
# ---------------------------------------------------------------------------
# Why 1.13 and not 1.12
#
# `juliac` moved out of the Julia distribution and into the JuliaC package, so
# it is a dependency of this project rather than a script under `Sys.BINDIR`.
# That is the mechanical difference. The substantive one is that 1.13 compiles
# task bodies that 1.12 discarded.
#
# A task's function is stored in the task object by `jl_new_task` and invoked
# later by the scheduler, in C. There is no Julia call site to it, so a
# reachability walk over the visible call graph never reaches it. On 1.12 the
# body is left out of the image and the program dies at run time with a
# `MethodError` -- and because there is no unresolved *call site*, `--trim=safe`
# reports zero errors and links happily. It is the one failure mode trimming
# does not catch.
#
# 1.13 roots lowering-generated closures used as task bodies, which covers
# `Threads.@threads`, `Threads.@spawn` and a hand-written `Task(() -> ...)`
# alike. So `transform!(...; threaded = true)` works here, where on 1.12 it
# needed `Base.Experimental.entrypoint` on a named callable struct -- the only
# shape 1.13 does *not* root for you.
#
# Moving to 1.13 also retired a workaround this example used to need.
# `HostCPUFeatures.__init__` reaches a `dlopen` of libLLVM to read the host CPU
# feature string, which 1.12's verifier could not resolve; the fix was to set
# that package's `freeze_cpu_target` preference, which made the branch
# statically dead, and to depend on it directly since Preferences apply only to
# direct dependencies. On 1.13 the build is clean without any of it, verified
# from a cleared cache with the preference reading `false`.
# ---------------------------------------------------------------------------

# Checked before instantiating: on an older Julia the resolve fails first, with
# a precompile error that says nothing about the actual problem.
VERSION >= v"1.13.0-" ||
    error("""
          Ahead-of-time compilation here needs Julia 1.13 or later; this is $VERSION.
          On 1.12 the build runs but the threaded transform is left out of the
          image, and `--trim=safe` does not report it. See the comment at the
          top of this file.
          """)

using Pkg
Pkg.instantiate()

const HERE = @__DIR__

trim = "safe"
out = "fastgeoproj"
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

# `--output-exe` takes a bare name and writes it to the working directory, so
# the build runs in this one.
cmd = Cmd(`$(Base.julia_cmd()) --startup-file=no --history-file=no --project=$HERE
           -m JuliaC --project $HERE --output-exe $out
           --experimental --trim=$trim fastgeoproj.jl`; dir = HERE)

@info "building" out trim
run(cmd)
@info "built" out size = Base.format_bytes(filesize(joinpath(HERE, out)))
