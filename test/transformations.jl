# Tests for the point-operator interface: one callable struct per
# transformation, applied to points, arrays of points, or coordinate vectors.

using FastGeoProjections
using FastGeoProjections: transform, transform!, islanesafe, GeoTransformation,
                          LonLatToPolarStereographic, PolarStereographicToLonLat
import GeoInterface as GI
using Test
using Proj
using StaticArrays, GeometryBasics

# Two point types for the array tests, defined here at top level: methods added
# inside a @testset land in its local scope and its own world age.
#
# MyPt happens to store its components in the buffer in x, y order...
struct MyPt; x::Float64; y::Float64; end
GI.isgeometry(::Type{MyPt}) = true
GI.geomtrait(::MyPt) = GI.PointTrait()
GI.ncoord(::GI.PointTrait, ::MyPt) = 2
GI.getcoord(::GI.PointTrait, p::MyPt, i) = i == 1 ? p.x : p.y

# ... and FlipPt the other way round, which no check on the type alone could
# tell apart from MyPt.
# isbits and exactly two Float64s wide, but its storage is not coordinates.
# Nothing may reinterpret coordinates into this and then call an accessor.
struct HandlePt; p::Ptr{Cvoid}; q::Ptr{Cvoid}; end

# a vector whose indices do not start at 1, without an OffsetArrays dependency
struct OffsetVec{T} <: AbstractVector{T}
    data::Vector{T}
    off::Int
end
Base.size(v::OffsetVec) = size(v.data)
Base.axes(v::OffsetVec) = (v.off .+ (1:length(v.data)),)
Base.IndexStyle(::Type{<:OffsetVec}) = IndexLinear()
Base.getindex(v::OffsetVec, i::Int) = v.data[i - v.off]
Base.setindex!(v::OffsetVec, x, i::Int) = (v.data[i - v.off] = x)
Base.similar(v::OffsetVec, ::Type{T}) where {T} = OffsetVec(similar(v.data, T), v.off)

struct FlipPt; y::Float64; x::Float64; end
GI.isgeometry(::Type{FlipPt}) = true
GI.geomtrait(::FlipPt) = GI.PointTrait()
GI.ncoord(::GI.PointTrait, ::FlipPt) = 2
GI.getcoord(::GI.PointTrait, p::FlipPt, i) = i == 1 ? p.x : p.y
FastGeoProjections.rebuildpoint(::Type{FlipPt}, x, y) = FlipPt(y, x)

@testset "point operators" begin
    @testset "EPSG:$(first(epsg.val))" for (epsg, latmin, latmax) in
            ((EPSG(3413), 60.0, 89.0), (EPSG(3031), -89.0, -60.0))
        n = 1000
        r = rand(n)
        lon = r .* 360 .- 180
        lat = latmin .+ r .* (latmax - latmin)

        fwd = Transformation(EPSG(4326), epsg; always_xy = true)
        pj = Proj.Transformation("EPSG:4326", "EPSG:$(first(epsg.val))"; always_xy = true)
        truth = [pj(lon[i], lat[i]) for i in 1:n]

        # a transformation is a point operator first
        @test fwd isa GeoTransformation
        @test all(abs.(fwd((lon[1], lat[1])) .- truth[1]) .< 1e-6)
        @test fwd(lon[1], lat[1]) == fwd((lon[1], lat[1]))

        # ...applied to a vector of points, in place or not
        pts = collect(zip(lon, lat))
        out = transform(fwd, pts)
        @test maximum(i -> maximum(abs, out[i] .- truth[i]), 1:n) < 1e-6
        v = copy(pts)
        @test transform!(fwd, v) === v
        @test v == out

        # ...or to a pair of coordinate vectors
        X, Y = transform(fwd, lon, lat)
        @test X == first.(out)
        @test Y == last.(out)
        @test fwd(lon, lat) == (X, Y)

        # threading must not change the result
        @test transform(fwd, pts; threaded = false) == transform(fwd, pts; threaded = true)

        # inverse round trip
        back = transform(inv(fwd), out)
        @test maximum(i -> abs(back[i][1] - lon[i]), 1:n) < 1e-8
        @test maximum(i -> abs(back[i][2] - lat[i]), 1:n) < 1e-8

        # Float32 data stays Float32, at Float32 accuracy
        X32, Y32 = transform(fwd, Float32.(lon), Float32.(lat))
        @test eltype(X32) === Float32
        @test maximum(abs, X32 .- Float32.(X)) < 1.0f-5 * maximum(abs, X)
    end

    @testset "composition and axis order" begin
        t = Transformation(EPSG(3413), EPSG(3031); always_xy = true)
        @test islanesafe(t)

        # always_xy swaps the geographic end only
        xy = Transformation(EPSG(4326), EPSG(3413); always_xy = true)
        authority = Transformation(EPSG(4326), EPSG(3413); always_xy = false)
        @test xy((-45.0, 70.0)) == authority((70.0, -45.0))

        # ∘ and inv compose
        f = LonLatToPolarStereographic(; lat_ts = 70.0, lon_0 = -45.0)
        @test (inv(f) ∘ f)((-45.0, 70.0))[2] ≈ 70.0
        @test inv(inv(f)) isa LonLatToPolarStereographic
    end

    @testset "math kernels agree" begin
        fast = Transformation(EPSG(4326), EPSG(3413); always_xy = true)
        base = Transformation(EPSG(4326), EPSG(3413); always_xy = true, kernel = BaseKernel())
        sleef = Transformation(EPSG(4326), EPSG(3413); always_xy = true, kernel = SLEEFKernel())
        @test islanesafe(fast)
        @test islanesafe(sleef)
        @test !islanesafe(base)          # Base's libm blocks vectorization
        p = (-45.0, 70.0)
        @test maximum(abs, collect(fast(p)) .- collect(base(p))) < 1e-6
        @test maximum(abs, collect(sleef(p)) .- collect(base(p))) < 1e-6
    end

    @testset "Proj fallback" begin
        t = Transformation(EPSG(4326), EPSG(28992); always_xy = true)
        @test !islanesafe(t)
        pj = Proj.Transformation("EPSG:4326", "EPSG:28992"; always_xy = true)
        @test all(abs.(t((5.39, 52.16)) .- pj(5.39, 52.16)) .< 1e-6)
        X, Y = transform(t, [5.39, 5.40], [52.16, 52.17])
        @test length(X) == 2
    end
end

@testset "transverse Mercator" begin
    using FastGeoProjections: LonLatToTransverseMercator, TransverseMercatorToLonLat,
                              LonLatToUTM, UTMToLonLat, convergence_scale, utm_epsg,
                              utmzone2epsg, epsg2utmzone, isutm, UTM_K0

    @testset "UTM zone $zone $(isnorth ? "N" : "S")" for (zone, isnorth) in
            ((1, true), (19, true), (31, true), (32, true), (60, true),
             (5, false), (36, false), (55, false))
        epsg = utmzone2epsg(zone, isnorth)
        code = first(epsg.val)
        lon0 = FastGeoProjections.utm_lon0(zone)

        n = 500
        r = rand(n)
        lon = lon0 .+ (r .* 6 .- 3)
        lat = isnorth ? r .* 84 : -(r .* 80)

        fwd = Transformation(EPSG(4326), epsg; always_xy = true)
        pj = Proj.Transformation("EPSG:4326", "EPSG:$code"; always_xy = true)
        truth = [pj(lon[i], lat[i]) for i in 1:n]

        @test islanesafe(fwd)
        @test all(abs.(fwd((lon[1], lat[1])) .- truth[1]) .< 1e-6)

        X, Y = transform(fwd, lon, lat)
        @test maximum(i -> max(abs(X[i] - truth[i][1]), abs(Y[i] - truth[i][2])), 1:n) < 1e-6

        pts = collect(zip(lon, lat))
        out = transform(fwd, pts)
        @test out == collect(zip(X, Y))
        @test transform(fwd, pts; threaded = false) == transform(fwd, pts; threaded = true)

        back = transform(inv(fwd), out)
        @test maximum(i -> max(abs(back[i][1] - lon[i]), abs(back[i][2] - lat[i])), 1:n) < 1e-9

        X32, Y32 = transform(fwd, Float32.(lon), Float32.(lat))
        @test eltype(X32) === Float32
        @test maximum(abs, X32 .- Float32.(X)) < 1.0f-4 * maximum(abs, X)
    end

    @testset "convergence and scale" begin
        lon0 = -69.0
        f = LonLatToTransverseMercator(; lon0, lat0 = 0.0)
        i = TransverseMercatorToLonLat(; lon0, lat0 = 0.0)

        # on the central meridian the grid is aligned with true north and
        # unstretched
        @test convergence_scale(f, (lon0, 45.0))[1] == 0.0
        @test convergence_scale(f, (lon0, 45.0))[2] ≈ 1.0

        # γ and k are properties of the point: both directions must agree
        for (lon, lat) in ((-68.0, 45.0), (-66.13, 78.39), (-90.0, 10.0), (-40.0, -60.0))
            p = f((lon, lat))
            gf, kf = convergence_scale(f, (lon, lat))
            gi, ki = convergence_scale(i, p)
            @test gf ≈ gi atol = 1e-9
            @test kf ≈ ki rtol = 1e-12
        end

        # ...and must match the projection's own derivative
        function finitediff(lon, lat; h = 1e-6)
            x1, y1 = f((lon, lat + h))
            x2, y2 = f((lon, lat - h))
            dx = (x1 - x2) / 2h
            dy = (y1 - y2) / 2h
            phi = deg2rad(lat)
            ds = f.a * (1 - f.e2) / (1 - f.e2 * sin(phi)^2)^1.5 * (pi / 180)
            (atand(-dx, dy), hypot(dx, dy) / ds)
        end
        for (lon, lat) in ((-68.0, 45.0), (-90.0, 10.0), (-40.0, -60.0))
            g0, k0 = finitediff(lon, lat)
            g1, k1 = convergence_scale(f, (lon, lat))
            @test g1 ≈ g0 atol = 1e-6
            @test k1 ≈ k0 rtol = 1e-7
        end

        # UTM reports the same convergence, with k scaled by the zone factor
        u = LonLatToUTM(19, true)
        @test all(convergence_scale(u, (lon0, 45.0)) .≈ (0.0, UTM_K0))
        for (lon, lat) in ((-68.0, 45.0), (-66.0, 60.0))
            gu, ku = convergence_scale(u, (lon, lat))
            gt, kt = convergence_scale(f, (lon, lat))
            @test gu == gt
            @test ku ≈ kt * UTM_K0
            # ...and the inverse operator agrees at the same physical point
            gi, ki = convergence_scale(inv(u), u(lon, lat))
            @test gi ≈ gu atol = 1e-9
            @test ki ≈ ku rtol = 1e-12
        end
    end

    @testset "off-equator origin" begin
        # lat0 != 0 shifts the northing origin; the old array kernels could
        # never run this branch (their series helpers were undefined)
        pj = Proj.Transformation("EPSG:4326",
            "+proj=tmerc +lat_0=30 +lon_0=-69 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m";
            always_xy = true)
        f = LonLatToTransverseMercator(; lon0 = -69.0, lat0 = 30.0)
        for (lon, lat) in ((-69.0, 30.0), (-66.0, 45.0), (-72.0, 10.0), (-69.0, -20.0))
            @test all(abs.(f((lon, lat)) .- pj(lon, lat)) .< 1e-6)
            @test all(abs.(inv(f)(f((lon, lat))) .- (lon, lat)) .< 1e-9)
        end
        @test f((-69.0, 30.0))[2] ≈ 0.0 atol = 1e-6
    end

    @testset "zone bookkeeping" begin
        @test isutm(EPSG(32619))
        @test isutm(EPSG(32736))
        @test !isutm(EPSG(4326))
        @test !isutm(EPSG(32600))       # zone 0 is not a UTM zone
        @test !isutm(EPSG(32661))       # nor is zone 61
        @test epsg2utmzone(EPSG(32619)) == (zone = 19, isnorth = true)
        @test epsg2utmzone(EPSG(32736)) == (zone = 36, isnorth = false)
        @test utmzone2epsg(19, true) == EPSG(32619)
        @test utm_epsg(-69.0, 45.0) == EPSG(32619)
        @test utm_epsg(45.0, -69.0, false) == EPSG(32619)   # authority order
        @test utm_epsg(30.0, -40.0) == EPSG(32736)
        @test_throws ErrorException epsg2utmzone(EPSG(4326))
    end

    @testset "UTM operator" begin
        t = LonLatToUTM(19, true)
        @test t.zone == 19
        @test t.isnorth
        @test islanesafe(t)
        @test inv(t) isa UTMToLonLat
        @test inv(t).zone == 19
        @test inv(inv(t)) isa LonLatToUTM
        @test t(-69.0, 45.0)[1] ≈ 5e5           # false easting on the meridian
        @test LonLatToUTM(19, false)(-69.0, -45.0)[2] ≈ 1e7 - t(-69.0, 45.0)[2]
        @test_throws ArgumentError LonLatToUTM(0, true)
        @test_throws ArgumentError LonLatToUTM(61, true)
        @test_throws ArgumentError UTMToLonLat(61, true)
        # constructing from an EPSG code picks the same zone
        @test LonLatToUTM(EPSG(32736)).zone == 36
        @test !LonLatToUTM(EPSG(32736)).isnorth
    end

    @testset "UTM precision is a type parameter" begin
        # `LonLatToUTM{T}(zone, isnorth)` is the form that keeps the operator's
        # type inferable when the zone is only known at run time -- what an
        # ahead-of-time compiled caller needs (see examples/juliac). The keyword
        # form is the same operator.
        @test LonLatToUTM{Float64}(19, true) === LonLatToUTM(19, true)
        @test UTMToLonLat{Float64}(19, true) === UTMToLonLat(19, true)
        @test LonLatToUTM{Float64}(EPSG(32736)).zone == 36

        # `Base.return_types` rather than `infer_return_type`: the latter is
        # 1.11+, and the package supports 1.10.
        @test only(Base.return_types(z -> LonLatToUTM{Float64}(z, true), (Int,))) ===
              typeof(LonLatToUTM(19, true))
        @test only(Base.return_types(z -> UTMToLonLat{Float32}(z, false), (Int,))) ===
              typeof(UTMToLonLat(19, false; T = Float32))

        t32 = LonLatToUTM{Float32}(19, true)
        @test t32.k0 isa Float32
        @test all(x -> x isa Float32, t32(-69.0f0, 45.0f0))
        @test_throws ArgumentError LonLatToUTM{Float64}(61, true)
    end

    @testset "UTM zone to UTM zone" begin
        # a composed pipeline stays a single lane-safe point operator
        t = Transformation(EPSG(32619), EPSG(32620); always_xy = true)
        @test islanesafe(t)
        pj = Proj.Transformation("EPSG:32619", "EPSG:32620"; always_xy = true)
        for p in ((500000.0, 5.0e6), (400000.0, 4.0e6))
            @test all(abs.(t(p) .- pj(p...)) .< 1e-6)
        end
        X, Y = transform(t, [500000.0, 400000.0], [5.0e6, 4.0e6])
        @test length(X) == 2
    end
end

@testset "point input" begin
    using FastGeoProjections: LonLatToUTM, Math, FastKernel, SLEEFKernel, BaseKernel,
                              ComposedGeoTransformation, SwapXY, Identity
    import GeoInterface as GI

    t = LonLatToUTM(19, true)
    ref = t(-69.0, 45.0)

    @testset "GeoInterface points" begin
        # the two-argument form is the primitive; one argument is a GI point
        @test t((-69.0, 45.0)) == ref
        @test t([-69.0, 45.0]) == ref
        @test t(GI.Point(-69.0, 45.0)) == ref
        # z is carried by the point but not by the transformation
        @test t((-69.0, 45.0, 100.0)) == ref
        @test t(GI.Point(-69.0, 45.0, 100.0)) == ref
        # a non-point geometry is a usage error, not a silent wrong answer
        @test_throws ArgumentError t(GI.LineString([(0.0, 0.0), (1.0, 1.0)]))
    end

    @testset "composition is a flat tuple" begin
        f = LonLatToUTM(19, true)
        g = inv(f)
        c = g ∘ f
        @test c isa ComposedGeoTransformation
        @test c.transformations == (f, g)          # applied first to last
        @test all(c((-69.0, 45.0)) .≈ (-69.0, 45.0))

        # nesting flattens rather than building a tree
        c3 = SwapXY() ∘ (g ∘ f)
        @test c3.transformations == (f, g, SwapXY())
        @test length((SwapXY() ∘ SwapXY() ∘ g ∘ f).transformations) == 4

        # Identity drops out, and a one-element chain is just the transformation
        @test (Identity() ∘ f) === f
        @test (f ∘ Identity()) === f
        # reversing a chain reverses the tuple and inverts each stage
        c2 = f ∘ SwapXY()
        @test c2.transformations == (SwapXY(), f)
        @test inv(c2).transformations == (g, SwapXY())
        @test islanesafe(c)
    end

    @testset "Math submodule" begin
        # Base's own functions are untouched
        @test Math.sin !== Base.sin
        @test Math.sin(BaseKernel(), 0.5) === sin(0.5)
        for K in (FastKernel(), SLEEFKernel(), BaseKernel())
            @test Math.sin(K, 0.5) ≈ sin(0.5) atol = 1e-15
            @test Math.cos(K, 0.5) ≈ cos(0.5) atol = 1e-15
            @test Math.pow(K, 1.3, 0.25) ≈ 1.3^0.25 atol = 1e-15
            @test all(Math.sincos(K, 0.5) .≈ sincos(0.5))
            # the fast sinh/cosh pair is accurate absolutely, not relatively
            @test all(abs.(Math.sinhcosh(K, 0.7) .- (sinh(0.7), cosh(0.7))) .< 1e-15)
        end
        @test FastGeoProjections.vectorizes(FastKernel())
        @test FastGeoProjections.vectorizes(SLEEFKernel())
        @test !FastGeoProjections.vectorizes(BaseKernel())

        @testset "conformal_ratio is a pow to one ulp" begin
            # Its contract is `e < 0.09` and `|u| <= e`; inside that the series
            # must be indistinguishable from the `pow` it replaces, since the
            # projections' agreement with PROJ rests on it.
            for K in (FastKernel(), SLEEFKernel(), BaseKernel())
                worst = 0.0
                for e in (0.0, 0.0033528, 0.081819190842621278, 0.0818191910428, 0.089)
                    for u in range(-e, e; length = 401)
                        want = ((1 - u) / (1 + u))^(e / 2)
                        got = Math.conformal_ratio(K, e, u)
                        worst = max(worst, abs(got - want) / abs(want))
                    end
                end
                @test worst <= 2eps(Float64)
            end
            # e = 0 is the sphere: the factor is exactly 1 and must not drift.
            @test Math.conformal_ratio(FastKernel(), 0.0, 0.0) === 1.0
            # Float32 carries through without widening.
            @test Math.conformal_ratio(FastKernel(), 0.0818f0, 0.05f0) isa Float32
        end
    end
end

@testset "arrays of GeoInterface points" begin
    using FastGeoProjections: LonLatToUTM, _interleaved

    t = LonLatToUTM(19, true)
    lons = collect(range(-72.0, -66.0; length = 37))
    lats = collect(range(40.0, 47.0; length = 37))

    mk(P) = [P(lo, la) for (lo, la) in zip(lons, lats)]
    mk(::Type{Tuple}) = [(lo, la) for (lo, la) in zip(lons, lats)]
    mk(::Type{FlipPt}) = [FlipPt(la, lo) for (lo, la) in zip(lons, lats)]  # fields are (y, x)
    xy(p) = (GI.x(p), GI.y(p))
    approx(v, r; rtol = 1e-12) = all(map((a, b) -> all(isapprox.(xy(a), b; rtol)), v, r))

    # The lane path contracts to FMA where the scalar operator does not, so the
    # reference for the array API is the array result, itself checked against
    # the point operator to a relative tolerance.
    ref = transform(t, mk(Tuple); threaded = false)
    @test approx(ref, [t(lo, la) for (lo, la) in zip(lons, lats)])

    @testset "layout detection" begin
        for v in (mk(Tuple), [(lo, la, 7.0) for (lo, la) in zip(lons, lats)],
                  mk(SVector{2,Float64}), mk(SVector{2,Float32}), mk(Point2{Float64}),
                  [Point3{Float64}(lo, la, 0.0) for (lo, la) in zip(lons, lats)],
                  mk(MyPt), [GI.Point((lo, la)) for (lo, la) in zip(lons, lats)])
            M = _interleaved(v)
            @test M !== nothing
            @test size(M, 1) == GI.ncoord(v[1])
            @test M[1, 3] == GI.x(v[3]) && M[2, 3] == GI.y(v[3])
        end
        # A reversed layout is caught, and caught because of what the type is
        # rather than because of what a particular vector happens to hold:
        # FlipPt is a conformant point (its `getcoord` order is x, y) that is
        # simply stored the other way round.
        @test GI.geomtrait(mk(FlipPt)[1]) isa GI.PointTrait
        @test isbitstype(FlipPt) && sizeof(FlipPt) == 2 * sizeof(Float64)
        @test _interleaved(mk(FlipPt)) === nothing
        @test _interleaved([[lo, la] for (lo, la) in zip(lons, lats)]) === nothing
        @test _interleaved(view(mk(Tuple), 1:10)) === nothing
    end

    @testset "layout is a property of the type, not of the values" begin
        # Sampling a vector's contents cannot answer this. Endpoints that lie
        # on x == y are indistinguishable from an x-major layout, however many
        # of them are checked, and the interior is then read transposed.
        diag = [FlipPt(v, v) for v in range(45.0, 46.0; length = 64)]
        diag[2] = FlipPt(45.5, -69.0)
        @test all(p -> GI.x(p) == getfield(p, :x), diag)
        @test _interleaved(diag) === nothing
        @test !FastGeoProjections.isxymajor(FlipPt)

        # GeoInterface accepts a NamedTuple as a point and lets it name its
        # coordinates in either order, so this is a y-major isbits point that
        # ships in GeoInterface itself.
        YX = NamedTuple{(:Y, :X),Tuple{Float64,Float64}}
        @test GI.geomtrait(YX((45.0, -69.0))) isa GI.PointTrait
        @test GI.coordnames(YX((45.0, -69.0))) == (:Y, :X)
        @test isbitstype(YX) && sizeof(YX) == 2 * sizeof(Float64)
        @test !FastGeoProjections.isxymajor(YX)
        @test _interleaved([YX((la, lo)) for (lo, la) in zip(lons, lats)]) === nothing

        # ...and the honest layouts are still recognised, from the type alone.
        for P in (NTuple{2,Float64}, NTuple{3,Float64}, SVector{2,Float64},
                  SVector{2,Float32}, Point2{Float64}, Point3{Float64}, MyPt,
                  typeof(GI.Point((1.0, 2.0))))
            @test FastGeoProjections.isxymajor(P)
        end

        # `isbitstype` on its own would not be enough to make the probe safe:
        # this is isbits and exactly two Float64s wide, and reinterpreting
        # coordinates into it would hand an accessor a wild pointer.
        @test isbitstype(HandlePt) && sizeof(HandlePt) == 2 * sizeof(Float64)
        @test !FastGeoProjections.isxymajor(HandlePt)
    end

    @testset "a fresh destination has no values to check" begin
        # `Vector{P}(undef, n)` of this size comes from zeroed pages, so a
        # transposed destination reads (0.0, 0.0) wherever it is sampled and
        # passes a comparison against its own contents trivially.
        n = 200_000
        src = [MyPt(-69.0 + 1e-7i, 45.0) for i in 1:n]
        dst = Vector{FlipPt}(undef, n)
        transform!(dst, t, src)
        # a transposed destination is off by orders of magnitude, not by ulp,
        # so the scalar path's couple-of-ulp tolerance is plenty to catch it
        @test approx(dst, [t(GI.x(p), GI.y(p)) for p in src])
    end

    @testset "every layout gives the same answer" begin
        for P in (Tuple, SVector{2,Float64}, Point2{Float64}, MyPt)
            for threaded in (false, true)
                @test xy.(transform(t, mk(P); threaded)) == ref
                @test xy.(transform!(t, mk(P); threaded)) == ref
            end
            @test eltype(transform(t, mk(P))) === eltype(mk(P))
        end
        # FlipPt has to go point by point, where the arithmetic does not
        # contract to FMA the way the lane path does: a couple of ulp
        for threaded in (false, true)
            @test approx(transform(t, mk(FlipPt); threaded), ref)
            @test approx(transform!(t, mk(FlipPt); threaded), ref)
        end
        @test eltype(transform(t, mk(FlipPt))) === FlipPt

        o32 = transform(t, mk(SVector{2,Float32}))    # Float32 points, Float32 operator
        @test eltype(eltype(o32)) === Float32
        @test approx(o32, ref; rtol = 1e-5)
    end

    @testset "z is carried through" begin
        zs = float.(1:37)
        v = [Point3{Float64}(lo, la, z) for (lo, la, z) in zip(lons, lats, zs)]
        for threaded in (false, true)
            for o in (transform(t, v; threaded), transform!(t, copy(v); threaded))
                @test xy.(o) == ref
                @test [p[3] for p in o] == zs
            end
        end
        v3 = [(lo, la, z) for (lo, la, z) in zip(lons, lats, zs)]
        @test transform(t, v3) == [(r[1], r[2], z) for (r, z) in zip(ref, zs)]

        # ...and on the scalar path too. Any transformation that is not lane
        # safe goes there whatever the point type: `BaseKernel` does not
        # vectorize, and neither does anything backed by Proj.
        tb = LonLatToUTM(19, true; kernel = BaseKernel())
        # EPSG:3395 is true ellipsoidal Mercator, which the package does not implement --
        # unlike EPSG:3857, which applies the spherical form to a geodetic latitude and is
        # native. So this one really is Proj-backed.
        tp = Transformation(EPSG(4326), EPSG(3395); always_xy = true)
        @test !FastGeoProjections.islanesafe(tb)
        @test !FastGeoProjections.islanesafe(tp)
        for ts in (tb, tp)
            expect = [ts(lo, la) for (lo, la) in zip(lons, lats)]
            for threaded in (false, true)
                o3 = transform(ts, v3; threaded)
                @test [(p[1], p[2]) for p in o3] == expect
                @test [p[3] for p in o3] == zs

                op = transform(ts, v; threaded)
                @test eltype(op) === Point3{Float64}
                @test xy.(op) == expect
                @test [p[3] for p in op] == zs
            end
        end
    end

    @testset "z is transformed where the transformation changes one" begin
        # Geographic 3D to geocentric: the height becomes a Cartesian
        # coordinate, so a carried-across z is not merely stale but of the
        # wrong quantity. Proj resolves the pipeline, and is the reference.
        tz = Transformation(EPSG(4979), EPSG(4978); always_xy = true)
        @test !FastGeoProjections.preservesz(tz)

        pj = Proj.Transformation("EPSG:4979", "EPSG:4978"; always_xy = true)
        lon, lat, h = 5.39, 52.16, 100.0
        want = pj(lon, lat, h)

        @test all(tz(lon, lat, h) .≈ want)                    # scalar, three args
        for threaded in (false, true)
            o = transform(tz, [(lon, lat, h)]; threaded)
            @test all(o[1] .≈ want)
            @test o[1][3] != h                                 # computed, not carried
            op = transform(tz, [Point3{Float64}(lon, lat, h)]; threaded)
            @test all(Tuple(op[1]) .≈ want)
        end

        # Dropping the height would move x and y as well: the two-argument call
        # is the h = 0 point, which is 61 m away here. So a three-component
        # source must not reach a path that transforms only x and y.
        @test !isapprox(pj(lon, lat)[1], want[1]; atol = 1.0)

        # A native operator whose whole subject is the height says so rather than
        # taking the h = 0 point: there is no height to transform in two
        # coordinates, and no useful default for one.
        @test_throws "would mean a height of zero" tz(lon, lat)
        @test_throws "would mean a height of zero" transform(tz, [(lon, lat)])

        # A Proj-backed pipeline that can change a height still takes x and y
        # alone, because 2D in and 2D out is a thing Proj resolves and is what
        # such a pipeline is usually asked for.
        tp = Transformation(EPSG(4326), EPSG(3395); always_xy = true)
        @test !FastGeoProjections.preservesz(tp)
        @test all(transform(tp, [(lon, lat)])[1] .≈
                  Proj.Transformation("EPSG:4326", "EPSG:3395"; always_xy = true)(lon, lat))

        # A map projection is a function of x and y, so it does preserve one.
        @test FastGeoProjections.preservesz(Transformation(EPSG(4326), EPSG(3413)))
        # ...and a chain is only as preserving as its least preserving stage.
        @test !FastGeoProjections.preservesz(
            Transformation(EPSG(4326), EPSG(3413)).f ∘ tz.f)
    end

    @testset "source and destination of different widths" begin
        src = [Point3{Float64}(lo, la, 0.0) for (lo, la) in zip(lons, lats)]
        dst = Vector{Point2{Float64}}(undef, length(src))
        @test xy.(transform!(dst, t, src)) == ref              # 3 rows into 2 rows
        dst32 = [SVector{2,Float32}(0, 0) for _ in src]        # eltype differs: scalar
        @test approx(transform!(dst32, t, src), ref; rtol = 1e-5)
    end

    @testset "the thread count does not change the answer" begin
        # sizes either side of a chunk boundary, so the threaded path is
        # actually taken rather than falling through to the single-chunk case
        C = FastGeoProjections.CHUNK
        for n in (1, 2, 3, 7, 8, 15, 16, 17, 1000, 1001, C, C + 1, 3C, 3C + 7)
            v = [(-72.0 + 6i / n, 40.0 + 7i / n) for i in 1:n]
            @test transform(t, v; threaded = false) == transform(t, v; threaded = true)
        end
    end

    @testset "transform is usable from inside a threaded loop" begin
        # `@threads :static` throws when nested or concurrent, and `threaded`
        # defaults to true whenever there is more than one thread, so this is
        # what a caller parallelising over tiles of their own hits. Only says
        # anything with -t2 or more; CI runs the suite threaded.
        n = 4 * FastGeoProjections.CHUNK
        v = [(-69.0 + 6i / n, 45.0) for i in 1:n]
        want = transform(t, v; threaded = false)
        out = [similar(v) for _ in 1:4]
        Threads.@threads for k in 1:4
            transform!(out[k], t, v)
        end
        @test all(o -> o == want, out)

        # ...and concurrently, which is the other thing :static refuses
        @test all(fetch.([Threads.@spawn transform(t, v) for _ in 1:4]) .== Ref(want))
    end

    @testset "a vector that does not start at 1" begin
        # The scalar path is where such a vector lands -- `_interleaved` takes
        # an `Array` -- and it used to be handed 1:n regardless, walking off
        # the end with bounds checking disabled.
        n = 1000
        raw = [(-72.0 + 6i / n, 40.0 + 7i / n) for i in 1:n]
        src = OffsetVec(copy(raw), -1)                 # indices 0:n-1
        dst = OffsetVec(similar(raw), -1)
        @test firstindex(src) == 0 && lastindex(src) == n - 1
        # ...and it lands on the scalar path, which does not contract to FMA
        # the way the lane path does: a couple of ulp, as for FlipPt above
        want = transform(t, raw; threaded = false)
        for threaded in (false, true)
            fill!(dst.data, (0.0, 0.0))
            transform!(dst, t, src; threaded)
            @test approx(dst.data, want)
        end
        # ...and the struct-of-arrays form the same way
        X = OffsetVec([p[1] for p in raw], -1)
        Y = OffsetVec([p[2] for p in raw], -1)
        Xd, Yd = transform(t, X, Y)
        @test all(isapprox.(collect(Xd), [p[1] for p in want]; rtol = 1e-12))
        @test all(isapprox.(collect(Yd), [p[2] for p in want]; rtol = 1e-12))
    end

    @testset "inv keeps the precision and the kernel" begin
        # neither is a field of `Transformation`; both live in the type of the
        # operator it wraps, so an inverse rebuilt from the EPSG pair would
        # quietly come back at Float64 with the default kernel
        t32 = FastGeoProjections.Transformation(EPSG(4326), EPSG(3413); T = Float32)
        @test typeof(inv(t32).f) === typeof(FastGeoProjections.pipeline(
            EPSG(3413), EPSG(4326); T = Float32))
        @test inv(inv(t32)).f === t32.f

        tb = FastGeoProjections.Transformation(EPSG(4326), EPSG(3413);
                                               kernel = BaseKernel())
        @test !FastGeoProjections.islanesafe(tb)
        @test !FastGeoProjections.islanesafe(inv(tb))

        # ...and it is still the inverse, for a native pair and a Proj-backed one
        for pair in ((4326, 32619), (4326, 3395))
            tt = FastGeoProjections.Transformation(EPSG(pair[1]), EPSG(pair[2]);
                                                   always_xy = true)
            @test all(isapprox.(inv(tt)(tt(-69.0, 45.0)...), (-69.0, 45.0); atol = 1e-9))
            @test inv(tt).source_epsg == tt.target_epsg
            @test inv(tt).target_epsg == tt.source_epsg
        end
    end

    @testset "authority axis order matches Proj across the native set" begin
        # `always_xy = false` means the axis order the authority defines, and
        # where the swap goes is decided by `isgeographic`, which is a list of
        # one code. A native CRS added without being classified there would
        # come back with x and y reversed rather than failing, so walk the
        # package's own native set -- adding to it puts the new code here.
        # a valid (lat, lon) for each code; the last branch is what a newly
        # added CRS gets, so that it is checked rather than erroring here
        pt(code) = code == 3031 ? (-75.0, 100.0) :
                   code == 3413 ? (75.0, -45.0) :
                   (32601 <= code <= 32660) ? (45.0, FastGeoProjections.utm_lon0(code - 32600)) :
                   (32701 <= code <= 32760) ? (-45.0, FastGeoProjections.utm_lon0(code - 32700)) :
                   (45.0, -69.0)
        native = vcat(FastGeoProjections.fast_epsgs,
                      [EPSG(c) for c in (32601, 32619, 32660, 32701, 32733, 32760)])
        @test EPSG(4326) in FastGeoProjections.fast_epsgs
        # the one list really is the one list: everything on it is claimed as
        # native, and everything on it has a projection behind it
        for c in FastGeoProjections.fast_epsg_codes
            @test FastGeoProjections.isfastepsg(EPSG(c))
            @test FastGeoProjections.project_to_4326(EPSG(c)) !== nothing
            @test FastGeoProjections.project_from_4326(EPSG(c)) !== nothing
        end
        for e in native
            code = first(e.val)
            code == 4326 && continue
            @test FastGeoProjections.isfastepsg(e)
            lat, lon = pt(code)
            for always_xy in (false, true)
                ours = Transformation(EPSG(4326), e; always_xy)
                pj = Proj.Transformation("EPSG:4326", "EPSG:" * string(code); always_xy)
                # `always_xy` decides which way round the input is read, so
                # feed both the same pair and let them disagree if it does
                a, b = always_xy ? (lon, lat) : (lat, lon)
                # A CRS that transforms a height has none to transform from a
                # two-coordinate call and rejects one, so it is fed three. Proj
                # reads (lat, lon, h) in authority order just as it reads
                # (lat, lon), so the height goes last either way.
                if FastGeoProjections.preservesz(ours)
                    @test all(isapprox.(ours(a, b), pj(a, b); rtol = 1e-6))
                else
                    @test all(isapprox.(ours(a, b, 100.0), pj(a, b, 100.0); rtol = 1e-6))
                end
            end
        end
    end

    @testset "the Proj pool is checked out, not indexed" begin
        # More concurrent tasks than the pool holds, so checkout has to block
        # and recycle rather than hand two tasks the same PJ. Needs a genuinely
        # Proj-backed pair, or there is no pool to exercise.
        tp = FastGeoProjections.Transformation(EPSG(4326), EPSG(3395); always_xy = true)
        @test tp.f isa FastGeoProjections.ProjTransformation
        n = 20_000
        v = [(-69.0 + 1e-4i, 45.0) for i in 1:n]
        want = transform(tp, v; threaded = false)
        outs = fetch.([Threads.@spawn transform(tp, v) for _ in 1:4*Threads.nthreads()])
        @test all(o -> o == want, outs)
        @test transform(tp, v; threaded = true) == want
        # ...and a bare point call, which checks one out for itself
        @test tp(v[1]...) == want[1]
    end

    @testset "through a Transformation" begin
        trans = FastGeoProjections.Transformation(EPSG(4326), EPSG(32619); always_xy = true)
        v = [SVector(lo, la) for (lo, la) in zip(lons, lats)]
        @test approx(trans(v), ref)                    # vector of points
        @test trans([lons[1], lats[1]]) == trans(lons[1], lats[1])   # ...vs one point
        @test eltype(trans(v)) === SVector{2,Float64}
    end

    @testset "edges" begin
        @test transform(t, NTuple{2,Float64}[]) == NTuple{2,Float64}[]
        @test_throws ArgumentError transform!(t, [1, 2, 3])
        @test_throws DimensionMismatch transform!(Vector{NTuple{2,Float64}}(undef, 2),
                                                  t, mk(Tuple))
    end
end

@testset "geocentric" begin
    fwd = LonLatToGeocentric()
    inverse = GeocentricToLonLat()
    pj = Proj.Transformation("EPSG:4979", "EPSG:4978"; always_xy = true)

    # Heights from below sea level to satellite altitude: the operator is used at
    # orbit as well as on the ground, and the inverse's accuracy is the thing
    # that varies with height.
    cases = [(lo, la, h) for la in -85.0:17.0:85.0, lo in -175.0:35.0:175.0,
                             h in (-500.0, 0.0, 1e3, 8e3, 7e5)]

    @testset "forward agrees with Proj" begin
        worst = 0.0
        for (lo, la, h) in cases
            worst = max(worst, maximum(abs.(fwd(lo, la, h) .- pj(lo, la, h))))
        end
        @test worst < 1e-8
    end

    @testset "the round trip closes at every height" begin
        # Proj's own inverse is 4.0e-3 m out at 700 km, so above the troposphere
        # the round trip rather than Proj is what pins the inverse.
        worst_h = 0.0
        worst_ground = 0.0
        for (lo, la, h) in cases
            xyz = fwd(lo, la, h)
            back = inverse(xyz...)
            worst_h = max(worst_h, abs(back[3] - h))
            worst_ground = max(worst_ground, maximum(abs.(fwd(back...) .- xyz)))
        end
        @test worst_h < 1e-8
        @test worst_ground < 1e-8
    end

    @testset "inverse agrees with Proj below satellite altitude" begin
        # The round trip above cannot catch an error the two directions share --
        # a wrong `e2` would close it perfectly -- so the inverse is compared to
        # an independent implementation as well. Confined to heights where Proj's
        # own inverse is trustworthy: at 700 km it is 4.0e-3 m out, which is why
        # the round trip is what covers the rest.
        pj_i = Proj.Transformation("EPSG:4978", "EPSG:4979"; always_xy = true)
        worst_horiz = 0.0
        worst_h = 0.0
        for la in -89.0:2.0:89.0, lo in -179.0:6.0:179.0, h in (-400.0, 0.0, 3e3, 1e4)
            xyz = pj(lo, la, h)
            got = inverse(xyz...)
            want = pj_i(xyz...)
            # Angles as ground distance, which is comparable across latitudes.
            worst_horiz = max(worst_horiz, abs(got[2] - want[2]) * 111320,
                              abs(got[1] - want[1]) * 111320 * cosd(la))
            worst_h = max(worst_h, abs(got[3] - want[3]))
        end
        @test worst_horiz < 1e-5
        @test worst_h < 1e-5
    end

    @testset "the poles and the antimeridian" begin
        # Where a projection breaks if it is going to. A pole has no defined
        # longitude, so only the latitude and the height are compared there.
        pj_i = Proj.Transformation("EPSG:4978", "EPSG:4979"; always_xy = true)
        for (lo, la, lbl) in ((0.0, 90.0, "north pole"), (0.0, -90.0, "south pole"),
                              (180.0, 45.0, "antimeridian"), (-180.0, 45.0, "antimeridian west"),
                              (179.9999, 0.0, "just inside the antimeridian"),
                              (0.0, 0.0, "origin"))
            want = pj(lo, la, 100.0)
            @test all(isapprox.(fwd(lo, la, 100.0), want; atol = 1e-6))

            got = inverse(want...)
            ref = pj_i(want...)
            @test got[2] ≈ ref[2] atol = 1e-9
            @test got[3] ≈ ref[3] atol = 1e-6
            abs(la) > 89.999 || @test got[1] ≈ ref[1] atol = 1e-9
        end
    end

    @testset "every kernel and precision against Proj" begin
        # Only the default is exercised elsewhere. `FastKernel` is the default and
        # costs about twice the error of the other two, which is the trade it
        # exists to make; all three are far inside any tolerance that matters.
        for (K, tol, lbl) in ((BaseKernel(), 1e-9, "Base"),
                              (SLEEFKernel(), 1e-9, "SLEEF"),
                              (FastKernel(), 1e-8, "Fast"))
            f = LonLatToGeocentric(kernel = K)
            worst = 0.0
            for la in -85.0:5.0:85.0, lo in -175.0:15.0:175.0
                worst = max(worst, maximum(abs.(f(lo, la, 100.0) .- pj(lo, la, 100.0))))
            end
            @test worst < tol
        end

        # Float32 is bounded by its own representable resolution at this
        # magnitude -- `eps(Float32) * 6.4e6` is about 0.76 m -- not by the
        # projection. Asserted so that the cost of asking for it is on record.
        f32 = LonLatToGeocentric{Float32}()
        worst32 = 0.0
        for la in -85.0:5.0:85.0, lo in -175.0:15.0:175.0
            worst32 = max(worst32,
                          maximum(abs.(Float64.(f32(Float32(lo), Float32(la), 100.0f0)) .-
                                       pj(lo, la, 100.0))))
        end
        @test 0.1 < worst32 < 5.0
        @test worst32 > eps(Float32) * 6.4e6 / 2
    end

    @testset "a height is transformed, not carried" begin
        @test !FastGeoProjections.preservesz(fwd)
        @test !FastGeoProjections.preservesz(inverse)
        for threaded in (false, true)
            o = transform(fwd, [(5.39, 52.16, 100.0)]; threaded)
            @test all(o[1] .≈ pj(5.39, 52.16, 100.0))
            @test o[1][3] != 100.0
        end
    end

    @testset "two coordinates are refused, not read as h = 0" begin
        # The h = 0 point is 61 m away in x here, so there is no defensible
        # default and the call is an error.
        @test_throws "would mean a height of zero" fwd(5.39, 52.16)
        @test_throws "z = 0" inverse(3.9e6, 3.7e5)
        @test !isapprox(pj(5.39, 52.16)[1], pj(5.39, 52.16, 100.0)[1]; atol = 1.0)
    end

    @testset "inv, precision and kernel" begin
        @test inv(fwd) isa GeocentricToLonLat
        @test inv(inverse) isa LonLatToGeocentric
        @test all(inv(fwd)(fwd(5.39, 52.16, 100.0)...) .≈ (5.39, 52.16, 100.0))

        f32 = LonLatToGeocentric{Float32}(kernel = BaseKernel())
        @test inv(f32) isa GeocentricToLonLat{Float32}
        @test inv(f32).kernel isa BaseKernel
        @test !FastGeoProjections.islanesafe(LonLatToGeocentric(kernel = BaseKernel()))
        @test FastGeoProjections.islanesafe(fwd)
    end

    @testset "through the EPSG registry" begin
        for (s, t) in ((4979, 4978), (4978, 4979))
            ours = Transformation(EPSG(s), EPSG(t); always_xy = true)
            theirs = Proj.Transformation("EPSG:$s", "EPSG:$t"; always_xy = true)
            @test FastGeoProjections.isfastepsg(EPSG(s))
            @test !FastGeoProjections.preservesz(ours)
            args = s == 4979 ? (5.39, 52.16, 100.0) :
                               (3.9036404612786868e6, 368315.27616670664, 5.013823349039822e6)
            @test all(isapprox.(ours(args...), theirs(args...); rtol = 1e-9))
        end

        # Geocentric to a 2-D projected CRS: the height is consumed by the
        # geocentric stage, and what reaches the projection is the geodetic
        # (lon, lat) it implies.
        ours = Transformation(EPSG(4978), EPSG(3413); always_xy = true)
        theirs = Proj.Transformation("EPSG:4978", "EPSG:3413"; always_xy = true)
        xyz = (3.9036404612786868e6, 368315.27616670664, 5.013823349039822e6)
        @test all(isapprox.(ours(xyz...)[1:2], theirs(xyz...)[1:2]; rtol = 1e-9))
    end

    @testset "fusing the geocentric stage into the projection" begin
        ll2xyz = Proj.Transformation("EPSG:4979", "EPSG:4978"; always_xy = true)

        # The registry returns the fused operator where the projection accepts a
        # Direction, and the plain composition where it does not.
        @test Transformation(EPSG(4978), EPSG(3413)).f isa
              FastGeoProjections.FusedFromGeocentric
        @test Transformation(EPSG(4978), EPSG(32619)).f isa
              FastGeoProjections.FusedFromGeocentric
        # ...and geographic targets still go through the ordinary path, since
        # there is no projection to fuse into.
        @test !(Transformation(EPSG(4978), EPSG(4979)).f isa
                FastGeoProjections.FusedFromGeocentric)

        @testset "$lbl" for (code, proj, lbl, lons) in (
                (3413, LonLatToPolarStereographic(; lat_ts = 70.0, lon_0 = -45.0),
                 "polar stereographic north", -175.0:25.0:175.0),
                (3031, LonLatToPolarStereographic(; lat_ts = -71.0, lon_0 = 0.0),
                 "polar stereographic south", -175.0:25.0:175.0),
                # A transverse Mercator zone is only meaningful near its own
                # meridian; far outside it the projection diverges, composed as
                # well as fused. So each zone is swept over its own ±18°.
                (32601, LonLatToUTM(1, true), "UTM zone 1N",
                 FastGeoProjections.utm_lon0(1) .+ (-18.0:6.0:18.0)),
                (32619, LonLatToUTM(19, true), "UTM zone 19N",
                 FastGeoProjections.utm_lon0(19) .+ (-18.0:6.0:18.0)),
                (32631, LonLatToUTM(31, true), "UTM zone 31N",
                 FastGeoProjections.utm_lon0(31) .+ (-18.0:6.0:18.0)),
                (32660, LonLatToUTM(60, true), "UTM zone 60N",
                 FastGeoProjections.utm_lon0(60) .+ (-18.0:6.0:18.0)),
                (32701, LonLatToUTM(1, false), "UTM zone 1S",
                 FastGeoProjections.utm_lon0(1) .+ (-18.0:6.0:18.0)),
                (32733, LonLatToUTM(33, false), "UTM zone 33S",
                 FastGeoProjections.utm_lon0(33) .+ (-18.0:6.0:18.0)),
                (32760, LonLatToUTM(60, false), "UTM zone 60S",
                 FastGeoProjections.utm_lon0(60) .+ (-18.0:6.0:18.0)))
            fused = Transformation(EPSG(4978), EPSG(code); always_xy = true)
            composed = proj ∘ GeocentricToLonLat()
            pj = Proj.Transformation("EPSG:4978", "EPSG:$code"; always_xy = true)

            worst_c = 0.0
            worst_p = 0.0
            for la in -80.0:10.0:80.0, lo in lons, h in (-200.0, 0.0, 5e3)
                xyz = ll2xyz(lo, la, h)
                f = fused(xyz...)
                # The fused operator is the composition with identities applied,
                # so it is the same function to a few ulps rather than an
                # approximation of it.
                worst_c = max(worst_c, maximum(abs.(f[1:2] .- composed(xyz...)[1:2])))
                worst_p = max(worst_p, maximum(abs.(f[1:2] .- pj(xyz...)[1:2])))
                # the height comes back from the geocentric stage either way
                @test f[3] ≈ h atol = 1e-6
            end
            @test worst_c < 1e-6
            @test worst_p < 1e-5
        end

        @testset "properties carry through the fusion" begin
            f = Transformation(EPSG(4978), EPSG(3413); always_xy = true).f
            @test !FastGeoProjections.preservesz(f)
            @test FastGeoProjections.islanesafe(f)
            @test !FastGeoProjections.islanesafe(
                FastGeoProjections.FusedFromGeocentric(
                    GeocentricToLonLat(kernel = BaseKernel()),
                    LonLatToPolarStereographic(; lat_ts = 70.0, lon_0 = -45.0)))
            # Two coordinates would mean a point on the equatorial plane.
            @test_throws "z = 0" f(3.9e6, 3.7e5)
            # Inverting gives the reverse pipeline, an ordinary composition:
            # `PolarStereographicToLonLat` produces angles and `LonLatToGeocentric`
            # consumes them, so there is no Direction to hand across. It closes
            # the round trip, and needs the height that the forward one returned.
            xyz = (3.9036404612786868e6, 368315.27616670664, 5.013823349039822e6)
            @test !FastGeoProjections.preservesz(inv(f))
            @test all(isapprox.(inv(f)(f(xyz...)...), xyz; rtol = 1e-9))
        end

        @testset "a projection that does not fuse still composes" begin
            # The fallback forms the angles and calls the projection, so a
            # projection that has not opted in gives the same answer.
            @test !FastGeoProjections.fuses_direction(GeocentricToLonLat())
            dir, h = FastGeoProjections.geocentric_direction(
                GeocentricToLonLat(), 3.9036404612786868e6, 368315.27616670664,
                5.013823349039822e6)
            p = LonLatToPolarStereographic(; lat_ts = 70.0, lon_0 = -45.0)
            @test all(FastGeoProjections.project_direction(p, dir) .≈
                      p(FastGeoProjections.lon_degrees(dir),
                        FastGeoProjections.lat_degrees(dir)))
            @test h ≈ 100.0 atol = 1e-6
        end
    end

    @testset "a non-WGS84 ellipsoid" begin
        grs80 = LonLatToGeocentric(ellips = FastGeoProjections.ellipsoid(EPSG(7019)))
        @test grs80.a == 6378137.0
        # GRS 1980 and WGS 84 differ in flattening only, by 1e-11 relative in e2,
        # so the positions differ by millimetres rather than not at all.
        d = maximum(abs.(grs80(5.39, 52.16, 100.0) .- fwd(5.39, 52.16, 100.0)))
        @test 0 < d < 1e-2
    end
end

@testset "polar stereographic inverse accuracy" begin
    # The conformal -> geodetic series is the least accurate part of any native
    # pipeline, so it gets its own bound rather than riding on the 1e-6 the
    # operator tests use -- that one is loose enough to pass a series two orders
    # of magnitude worse, which is how a truncated one went unnoticed.
    #
    # Bounds are absolute metres of ground distance. A degree of latitude is
    # ~111 km, so the comparison is scaled by that rather than left in degrees.
    @testset "EPSG:$code" for (code, lat_ts, lon_0, ymin, ymax) in
            ((3413, 70.0, -45.0, -2.6e6, -1.4e6), (3031, -71.0, 0.0, 1.4e6, 2.6e6))
        t = PolarStereographicToLonLat(; lat_ts, lon_0)
        pj = Proj.Transformation("EPSG:$code", "EPSG:4326"; always_xy = true)
        worst = 0.0
        for x in range(-3.0e5, 3.0e5; length = 40), y in range(ymin, ymax; length = 40)
            got, want = t(x, y), pj(x, y)
            worst = max(worst, abs(got[1] - want[1]), abs(got[2] - want[2]))
        end
        @test worst * 111320 < 1e-7
    end

    @testset "round trip closes" begin
        fwd = LonLatToPolarStereographic(; lat_ts = 70.0, lon_0 = -45.0)
        rev = PolarStereographicToLonLat(; lat_ts = 70.0, lon_0 = -45.0)
        worst = 0.0
        for lon in range(-179.0, 179.0; length = 40), lat in range(30.0, 89.5; length = 40)
            x, y = fwd(lon, lat)
            lo, la = rev(x, y)
            worst = max(worst, abs(lo - lon), abs(la - lat))
        end
        @test worst * 111320 < 1e-7
    end

    @testset "the pole is a point, not a NaN" begin
        # The origin is the pole. A Newton solve on the exact relation would form
        # `(1 - t^2)/(2t)` here and return NaN; the series is defined at t = 0,
        # which is why it is the one used. Polar stereographic data sits here.
        rev = PolarStereographicToLonLat(; lat_ts = 70.0, lon_0 = -45.0)
        @test rev(0.0, 0.0) == (-45.0, 90.0)
        @test all(isfinite, rev(1e-9, 1e-9))
        south = PolarStereographicToLonLat(; lat_ts = -71.0, lon_0 = 0.0)
        @test south(0.0, 0.0) == (0.0, -90.0)
    end

    @testset "every kernel agrees" begin
        # The series is plain arithmetic, so the back-end changes nothing here.
        # A difference would mean a transcendental is being asked for accuracy it
        # does not have, rather than the series being the limit.
        args = (; lat_ts = 70.0, lon_0 = -45.0)
        base = PolarStereographicToLonLat(; args..., kernel = BaseKernel())
        for K in (FastKernel(), SLEEFKernel())
            t = PolarStereographicToLonLat(; args..., kernel = K)
            worst = 0.0
            for x in range(-3.0e5, 3.0e5; length = 30), y in range(-2.6e6, -1.4e6; length = 30)
                worst = max(worst, maximum(abs, t(x, y) .- base(x, y)))
            end
            @test worst * 111320 < 1e-8
        end
    end
end

@testset "web Mercator" begin
    # EPSG:3857 applies spherical Mercator to the *geodetic* latitude, which is the
    # authority's definition rather than an approximation of it. Reading the same numbers
    # as true ellipsoidal Mercator would put a point kilometres out, so agreement with
    # Proj is the whole specification here.
    fwd = Transformation(EPSG(4326), EPSG(3857); always_xy = true)
    rev = Transformation(EPSG(3857), EPSG(4326); always_xy = true)
    pj_f = Proj.Transformation("EPSG:4326", "EPSG:3857"; always_xy = true)
    pj_r = Proj.Transformation("EPSG:3857", "EPSG:4326"; always_xy = true)

    @testset "against Proj" begin
        worst_f = 0.0
        worst_r = 0.0
        for lon in -180.0:5.0:180.0, lat in -85.0:2.5:85.0
            x, y = fwd(lon, lat)
            px, py = pj_f(lon, lat)
            worst_f = max(worst_f, abs(x - px), abs(y - py))
            lo, la = rev(px, py)
            qlo, qla = pj_r(px, py)
            # Longitude is circular, so ±180 is one meridian rather than two.
            dlon = abs(lo - qlo)
            worst_r = max(worst_r, min(dlon, 360 - dlon), abs(la - qla))
        end
        @test worst_f < 1e-7
        @test worst_r < 1e-12
    end

    @testset "the tile boundary is exact" begin
        # The latitude web maps clip at, where y is the half-circumference. A projection
        # that reached it by a series rather than in closed form would not land on it.
        x, y = fwd(180.0, 85.05112877980659)
        @test x ≈ 20037508.342789244 rtol = 1e-14
        @test y ≈ 20037508.342789244 rtol = 1e-14
        @test fwd(0.0, 0.0) === (0.0, 0.0)
    end

    @testset "round trip closes" begin
        worst = 0.0
        for lon in -179.0:7.0:179.0, lat in -85.0:5.0:85.0
            lo, la = rev(fwd(lon, lat)...)
            worst = max(worst, abs(lo - lon), abs(la - lat))
        end
        # Degrees; 1e-13° is about 1e-8 m of ground distance.
        @test worst < 1e-12
    end

    @testset "operators invert and adapt" begin
        f = LonLatToWebMercator()
        @test inv(f) isa WebMercatorToLonLat
        @test inv(inv(f)) isa LonLatToWebMercator
        # `inv` must recover the same sphere rather than reverting to a default.
        @test inv(inv(f))(-45.0, 70.0) === f(-45.0, 70.0)
        @test FastGeoProjections.islanesafe(f)
        @test FastGeoProjections.preservesz(f)
        g = FastGeoProjections.adapt_eltype(f, Float32)
        @test g isa LonLatToWebMercator{Float32}
        @test FastGeoProjections.adapt_eltype(f, Float64) === f
    end

    @testset "every kernel agrees" begin
        base = LonLatToWebMercator(kernel = BaseKernel())
        for K in (FastKernel(), SLEEFKernel())
            t = LonLatToWebMercator(kernel = K)
            worst = 0.0
            for lon in -180.0:15.0:180.0, lat in -85.0:5.0:85.0
                worst = max(worst, maximum(abs, t(lon, lat) .- base(lon, lat)))
            end
            @test worst < 1e-7
        end
    end

    @testset "composes with the other projections" begin
        # Web Mercator is projected, so it sits on the x,y side of a pipeline; a composition
        # through EPSG:4326 must still agree with Proj resolving the same pair directly.
        for (a, b) in ((3857, 3413), (3413, 3857), (3857, 32619), (32619, 3857))
            ours = Transformation(EPSG(a), EPSG(b); always_xy = true)
            pj = Proj.Transformation("EPSG:$a", "EPSG:$b"; always_xy = true)
            src = Transformation(EPSG(4326), EPSG(a); always_xy = true)
            for (lon, lat) in ((-45.0, 70.0), (-69.0, 45.0), (10.0, -30.0))
                x, y = src(lon, lat)
                @test all(isapprox.(ours(x, y), pj(x, y); rtol = 1e-9))
            end
        end
    end
end
