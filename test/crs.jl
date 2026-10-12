# `NativeCRS` transformations: the same operator, bit for bit, as the EPSG-code constructor with
# `always_xy = true`, and a type that follows from the CRS types alone.

using FastGeoProjections: epsg

@testset "NativeCRS" begin
    codes = (4326, 32607, 32622, 32719, 3413, 3031)
    pts = Dict(4326 => [(-140.5, 61.2), (-45.3, 69.1), (-70.2, -45.0), (10.0, -78.0)],
               32607 => [(336885.0, 6.776715e6), (500000.0, 6.5e6)],
               32622 => [(531000.0, 7.66e6), (400000.0, 7.9e6)],
               32719 => [(500000.0, 5.0e6), (300000.0, 4.2e6)],
               3413 => [(-200000.0, -2.2e6), (500000.0, -1.0e6)],
               3031 => [(-1.5e6, 0.5e6), (2.0e5, -1.0e6)])

    @testset "EPSG:$s to EPSG:$t" for s in codes, t in codes
        s == t && continue
        typed = FastGeoProjections.Transformation(native_crs(s), native_crs(t))
        coded = FastGeoProjections.Transformation(EPSG(s), EPSG(t); always_xy = true)
        @test typeof(typed) === typeof(coded)
        for (x, y) in pts[s]
            @test typed(x, y) === coded(x, y)
            @test inv(typed)(typed(x, y)...) === inv(coded)(coded(x, y)...)
        end
        @test epsg(native_crs(s)) == EPSG(s)
    end

    # The type depends on the kinds of CRS, not on which zone or hemisphere.
    @test typeof(FastGeoProjections.Transformation(UTMCRS(7, true), PolarStereographicCRS(3413))) ===
          typeof(FastGeoProjections.Transformation(UTMCRS(58, false), PolarStereographicCRS(3031)))
    @test typeof(FastGeoProjections.Transformation(UTMCRS(7, true), UTMCRS(8, true), Float32)) ===
          typeof(FastGeoProjections.Transformation(EPSG(32607), EPSG(32608); always_xy = true, T = Float32))

    @test native_crs(32622) === UTMCRS(22, true)
    @test native_crs(32719) === UTMCRS(19, false)
    @test native_crs(4326) === GeographicCRS()
    @test_throws "EPSG:3857 has no NativeCRS" native_crs(3857)
    @test_throws "UTM zone must be in 1:60" UTMCRS(61, true)
    @test_throws "EPSG:3995 is not a polar stereographic CRS" PolarStereographicCRS(3995)
end
