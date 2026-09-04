"""
    FastGeoProjections.Math

Transcendental functions for the projections, each taking a
[`MathKernel`](@ref) as its first argument:

    Math.sin(kernel, x)
    Math.sincos(kernel, x)
    Math.pow(kernel, x, y)

The names shadow `Base`'s inside this module and are meant to be called
qualified (`Math.sin(...)`), so a projection can be written once and evaluated
through whichever back-end the caller picked -- including on SIMD lanes, which
is what lets a point operator reach array-kernel throughput.

`@turbo` used to do this rewriting invisibly: it substituted SLEEFPirates'
`*_fast` routines for `sin`, `cos`, `tan` and `^`
(`LoopVectorization/src/modeling/costs.jl`). Those routines are plain Julia, so
they inline into a point operator and evaluate on `Vec`s just as well as they
did inside a `@turbo` loop. Base's libm does neither: it costs about 1.3x more
per scalar and blocks vectorization entirely.
"""
module Math

import SLEEFPirates

export MathKernel, FastKernel, SLEEFKernel, BaseKernel

"""
    MathKernel

Back-end used to evaluate transcendental functions. See [`FastKernel`](@ref),
[`SLEEFKernel`](@ref) and [`BaseKernel`](@ref).
"""
abstract type MathKernel end

"""
    FastKernel()

SLEEFPirates' `*_fast` routines -- the ones `@turbo` lowered to. Vectorizes;
about twice the speed of [`SLEEFKernel`](@ref) for the same ~1 ULP accuracy
over the argument ranges a projection uses. The difference is argument
reduction: the fast routines use a cheap reduction that only starts to lose
digits for arguments far from zero, which a projection never produces. This is
the default.
"""
struct FastKernel <: MathKernel end

"""
    SLEEFKernel()

SLEEFPirates' fully range-reduced routines: sub-ULP for any finite argument,
at roughly twice the cost of [`FastKernel`](@ref). Also vectorizes. Worth
choosing only if you feed a projection arguments far outside its normal
domain.
"""
struct SLEEFKernel <: MathKernel end

"""
    BaseKernel()

Base's libm. Correctly rounded, but it does not vectorize -- a transformation
built with it falls back to a scalar loop. Provided as an accuracy reference.
"""
struct BaseKernel <: MathKernel end

const DEFAULT_KERNEL = FastKernel()

"""
    vectorizes(kernel)

Whether `kernel`'s routines are pure Julia, and so can be evaluated on
`VectorizationBase.Vec` lanes.
"""
vectorizes(::MathKernel) = true
vectorizes(::BaseKernel) = false

# ---------------------------------------------------------------------------
# Base's libm
# ---------------------------------------------------------------------------
@inline sin(::BaseKernel, x) = Base.sin(x)
@inline cos(::BaseKernel, x) = Base.cos(x)
@inline tan(::BaseKernel, x) = Base.tan(x)
@inline asin(::BaseKernel, x) = Base.asin(x)
@inline atan(::BaseKernel, x) = Base.atan(x)
@inline atan(::BaseKernel, y, x) = Base.atan(y, x)
@inline sinh(::BaseKernel, x) = Base.sinh(x)
@inline cosh(::BaseKernel, x) = Base.cosh(x)
@inline tanh(::BaseKernel, x) = Base.tanh(x)
@inline asinh(::BaseKernel, x) = Base.asinh(x)
@inline atanh(::BaseKernel, x) = Base.atanh(x)
@inline exp(::BaseKernel, x) = Base.exp(x)
@inline log(::BaseKernel, x) = Base.log(x)
@inline pow(::BaseKernel, x, y) = x^y
@inline cbrt(::BaseKernel, x) = Base.cbrt(x)
@inline sincos(::BaseKernel, x) = Base.sincos(x)
@inline sinhcosh(::BaseKernel, x) = (Base.sinh(x), Base.cosh(x))

# ---------------------------------------------------------------------------
# SLEEF, fully range-reduced
# ---------------------------------------------------------------------------
@inline sin(::SLEEFKernel, x) = SLEEFPirates.sin(x)
@inline cos(::SLEEFKernel, x) = SLEEFPirates.cos(x)
@inline tan(::SLEEFKernel, x) = SLEEFPirates.tan(x)
@inline asin(::SLEEFKernel, x) = SLEEFPirates.asin(x)
@inline atan(::SLEEFKernel, x) = SLEEFPirates.atan(x)
@inline atan(::SLEEFKernel, y, x) = SLEEFPirates.atan(y, x)
@inline sinh(::SLEEFKernel, x) = SLEEFPirates.sinh(x)
@inline cosh(::SLEEFKernel, x) = SLEEFPirates.cosh(x)
@inline tanh(::SLEEFKernel, x) = SLEEFPirates.tanh(x)
@inline asinh(::SLEEFKernel, x) = SLEEFPirates.asinh(x)
@inline atanh(::SLEEFKernel, x) = SLEEFPirates.atanh(x)
@inline exp(::SLEEFKernel, x) = SLEEFPirates.exp(x)
@inline log(::SLEEFKernel, x) = SLEEFPirates.log(x)
@inline pow(::SLEEFKernel, x, y) = SLEEFPirates.pow(x, y)
@inline cbrt(::SLEEFKernel, x) = SLEEFPirates.cbrt(x)
@inline sincos(::SLEEFKernel, x) = SLEEFPirates.sincos(x)
@inline sinhcosh(::SLEEFKernel, x) = (SLEEFPirates.sinh(x), SLEEFPirates.cosh(x))

# ---------------------------------------------------------------------------
# SLEEF, cheap argument reduction
# ---------------------------------------------------------------------------
@inline sin(::FastKernel, x) = SLEEFPirates.sin_fast(x)
@inline cos(::FastKernel, x) = SLEEFPirates.cos_fast(x)
@inline tan(::FastKernel, x) = SLEEFPirates.tan_fast(x)
@inline asin(::FastKernel, x) = SLEEFPirates.asin_fast(x)
@inline atan(::FastKernel, x) = SLEEFPirates.atan_fast(x)
@inline atan(::FastKernel, y, x) = SLEEFPirates.atan_fast(y, x)
@inline sinh(::FastKernel, x) = SLEEFPirates.sinh(x)
@inline cosh(::FastKernel, x) = SLEEFPirates.cosh(x)
@inline tanh(::FastKernel, x) = SLEEFPirates.tanh_fast(x)
@inline asinh(::FastKernel, x) = SLEEFPirates.asinh(x)
@inline atanh(::FastKernel, x) = SLEEFPirates.atanh(x)
@inline exp(::FastKernel, x) = SLEEFPirates.exp(x)
@inline log(::FastKernel, x) = SLEEFPirates.log_fast(x)
@inline pow(::FastKernel, x, y) = SLEEFPirates.pow_fast(x, y)
@inline cbrt(::FastKernel, x) = SLEEFPirates.cbrt_fast(x)
@inline sincos(::FastKernel, x) = SLEEFPirates.sincos_fast(x)

"""
    conformal_ratio(kernel, e, u)

`((1 − u) / (1 + u))^(e/2)`, for `u = e·sin(φ)` on an ellipsoid of eccentricity
`e` — the isometric-latitude factor both the polar stereographic and transverse
Mercator projections need, and `exp(−e·atanh(u))` written so that it vectorizes.

A general `pow` is the wrong instrument here. It costs about five times any
other transcendental in a loop, and that loop does not vectorize at all: LLVM
will not vectorize a body holding both halves of `pow_fast`'s
`exp2(y * log2_fast(x))` inlined, though each half vectorizes on its own. So a
`pow` caps a projection's throughput however well the loop around it is written,
and no general implementation can avoid it — handling negative bases, integer
exponents and the infinities is what makes one too large to vectorize.

Neither series needs many terms, because both arguments are small. Ellipsoids in
use have `e < 0.09`, so `|u| ≤ e` bounds `atanh`'s odd series in `u²`, and the
resulting `|e·atanh(u)| < 7e-3` bounds `exp`'s. Six terms each reach one ulp
across the whole range, in twelve fused multiply-adds with no branch, no table
and no exponent manipulation.

The `e < 0.09` bound is the contract. Every ellipsoid in [`ellipsoid`](@ref)
satisfies it by a wide margin — Earth's flattening is 1/298 — and so does any
plausible addition, but a genuinely eccentric body would need more terms and is
outside what this is derived for.
"""
@inline function conformal_ratio(K::MathKernel, e, u)
    v = u * u
    # atanh(u) / u, in v = u²
    b = muladd(v, oftype(v, 1 / 11), oftype(v, 1 / 9))
    b = muladd(v, b, oftype(v, 1 / 7))
    b = muladd(v, b, oftype(v, 1 / 5))
    b = muladd(v, b, oftype(v, 1 / 3))
    b = muladd(v, b, one(v))
    w = -e * u * b
    # exp(w)
    p = muladd(w, oftype(w, 1 / 120), oftype(w, 1 / 24))
    p = muladd(w, p, oftype(w, 1 / 6))
    p = muladd(w, p, oftype(w, 1 / 2))
    p = muladd(w, p, one(w))
    muladd(w, p, one(w))
end

"""
    sinhcosh(kernel, x)

Both hyperbolic functions of `x`.

Under [`FastKernel`](@ref) one `exp` supplies both, for a third of the cost of
calling `sinh` and `cosh` separately. The subtraction cancels as `x` goes to
zero, so that form is accurate in the absolute sense only -- right where the
pair multiplies a series coefficient, wrong for a `sinh` whose own magnitude
carries the result.
"""
@inline function sinhcosh(::FastKernel, x)
    ex = SLEEFPirates.exp(x)
    exi = inv(ex)
    ((ex - exi) / 2, (ex + exi) / 2)
end

end # module Math
