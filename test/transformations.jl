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
    end

    @testset "source and destination of different widths" begin
        src = [Point3{Float64}(lo, la, 0.0) for (lo, la) in zip(lons, lats)]
        dst = Vector{Point2{Float64}}(undef, length(src))
        @test xy.(transform!(dst, t, src)) == ref              # 3 rows into 2 rows
        dst32 = [SVector{2,Float32}(0, 0) for _ in src]        # eltype differs: scalar
        @test approx(transform!(dst32, t, src), ref; rtol = 1e-5)
    end

    @testset "the thread count does not change the answer" begin
        for n in (1, 2, 3, 7, 8, 15, 16, 17, 1000, 1001)
            v = [(-72.0 + 6i / n, 40.0 + 7i / n) for i in 1:n]
            @test transform(t, v; threaded = false) == transform(t, v; threaded = true)
        end
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
