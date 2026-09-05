# EPSG -> transformation registry.
#
# Every native projection is expressed relative to EPSG:4326 in (lon, lat)
# order, so an arbitrary pair of CRSs becomes `to_target ∘ from_source`, which
# composes into a single point operator.

## ⬇ ADD FAST PROJECTIONS HERE ⬇ ##

"""
    project_to_4326(epsg; T = Float64, kernel = DEFAULT_KERNEL)

Point operator from `epsg` to EPSG:4326 in `(lon, lat)` order, or `nothing`
when FastGeoProjections has no native implementation.
"""
function project_to_4326(epsg::EPSG; T::Type = Float64, kernel::MathKernel = DEFAULT_KERNEL)
    code = first(epsg.val)
    if code == 4326
        Identity()
    elseif code == 3031
        PolarStereographicToLonLat{T}(; lat_ts = -71.0, lon_0 = 0.0, kernel)
    elseif code == 3413
        PolarStereographicToLonLat{T}(; lat_ts = 70.0, lon_0 = -45.0, kernel)
    elseif code == 3857
        WebMercatorToLonLat{T}(; kernel)
    elseif code == 4978
        GeocentricToLonLat{T}(; kernel)
    elseif code == 4979
        Identity()
    elseif isutm(epsg)
        UTMToLonLat(epsg; T, kernel)
    else
        nothing
    end
end

"""
    project_from_4326(epsg; T = Float64, kernel = DEFAULT_KERNEL)

Point operator from EPSG:4326 in `(lon, lat)` order to `epsg`, or `nothing`
when FastGeoProjections has no native implementation.
"""
function project_from_4326(epsg::EPSG; T::Type = Float64, kernel::MathKernel = DEFAULT_KERNEL)
    code = first(epsg.val)
    if code == 4326
        Identity()
    elseif code == 3031
        LonLatToPolarStereographic{T}(; lat_ts = -71.0, lon_0 = 0.0, kernel)
    elseif code == 3413
        LonLatToPolarStereographic{T}(; lat_ts = 70.0, lon_0 = -45.0, kernel)
    elseif code == 3857
        LonLatToWebMercator{T}(; kernel)
    elseif code == 4978
        LonLatToGeocentric{T}(; kernel)
    elseif code == 4979
        Identity()
    elseif isutm(epsg)
        LonLatToUTM(epsg; T, kernel)
    else
        nothing
    end
end

"""
    fast_epsg_codes

The non-UTM CRSs FastGeoProjections implements natively, as authority codes.

One list, read by everything that needs it: [`isfastepsg`](@ref) and
[`fast_epsgs`](@ref) derive from it, and the test suite walks it to check each
code's axis order against Proj. That is the point of it being one list rather
than three -- a CRS added here but not classified in [`isgeographic`](@ref)
fails the suite instead of quietly returning x and y the wrong way round, and
one added here without a projection in [`project_to_4326`](@ref) fails too.
"""
const fast_epsg_codes = (3031, 3413, 3857, 4326, 4978, 4979)

"""
    fast_epsgs

[`fast_epsg_codes`](@ref) as `EPSG`s.
"""
const fast_epsgs = [EPSG(c) for c in fast_epsg_codes]

"""
    isfastepsg(epsg)

Whether FastGeoProjections implements `epsg` natively.
"""
isfastepsg(epsg::EPSG) = first(epsg.val) in fast_epsg_codes || isutm(epsg)
isfastepsg(source::EPSG, target::EPSG) = isfastepsg(source) && isfastepsg(target)

"""
    isgeographic(epsg)

Whether `epsg` is a geographic CRS, i.e. one whose authority axis order is
latitude first. These are the CRSs affected by `always_xy`.

Defined over the CRSs the package handles natively -- [`fast_epsg_codes`](@ref)
plus UTM -- and not meant as a general answer: EPSG:4269 and EPSG:4258 are
geographic too, and this says otherwise. Anything else goes to Proj, which
knows. A native CRS added without being classified here would have its axis
order silently reversed under `always_xy = false` rather than failing, so the
test suite walks `fast_epsg_codes` and checks each against Proj at both axis
orders.
"""
isgeographic(epsg::EPSG) = first(epsg.val) in (4326, 4979)

"""
    pipeline(source_epsg, target_epsg; always_xy, proj_only, T, kernel)

Build the point operator that takes a coordinate in `source_epsg` to
`target_epsg`. Axis order is handled by composing [`SwapXY`](@ref) onto the
geographic end of the pipeline, so no per-point branch is needed.
"""
function pipeline(source_epsg::EPSG, target_epsg::EPSG;
                  always_xy::Bool = false, proj_only::Bool = false,
                  T::Type = Float64, kernel::MathKernel = DEFAULT_KERNEL)
    if proj_only || !isfastepsg(source_epsg, target_epsg)
        return ProjTransformation(source_epsg, target_epsg; always_xy)
    end
    from = project_to_4326(source_epsg; T, kernel)
    to = project_from_4326(target_epsg; T, kernel)
    # A geocentric source hands the projection after it a `Direction` rather than
    # a pair of angles, so the trigonometry between the two stages cancels; see
    # `fuse.jl`. This is the one place the pipeline is not literally `to ∘ from`,
    # and it is the same function either way.
    if from isa GeocentricToLonLat && fuses_direction(to) && !isgeographic(target_epsg)
        return FusedFromGeocentric(from, to)
    end
    # native projections speak (x, y) / (lon, lat); swap at a geographic end
    # when the caller wants authority order
    if !always_xy
        isgeographic(source_epsg) && (from = from ∘ SwapXY())
        isgeographic(target_epsg) && (to = SwapXY() ∘ to)
    end
    to ∘ from
end
