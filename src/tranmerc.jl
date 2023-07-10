"""
    tranmerc_fwd(lon, lat; lon0, lat0, ellips, threaded, always_xy

Returns x and y coordinates [meteres] in transverse Mercator (TM) projection given geodetic 
coordinates (EPSG:4326) of longitude and latitude [decimal degrees]. The TM projection is 
defined by kwargs of longitude (lon0) and latitude (lat0), which specify the center of the 
projeciton, and an ellipsoid that is define by an equatorial radius in meters (a) and its 
eccentricity (e). Also returnes meridian convergence (gam) and scale factor (k).
"""

function tranmerc_fwd(
    lon::Union{AbstractVector{<:AbstractFloat},AbstractFloat}, 
    lat::Union{AbstractVector{<:AbstractFloat},AbstractFloat}; 
    lon0::Real=0, 
    lat0::Real=0, 
    ellips::Ellipsoid=ellipsoid(EPSG(7030)), 
    threaded = true, 
    always_xy = true
    )
# This code has been heavily modified for maximum preformance in Julia from the MATLAB 
# implimentation of geographiclib_toolbox-2.0 by Alex Gardner JPL/NASA.
#
# This implementation of the projection is based on the series method
#   described in
#
#     C. F. F. Karney, Transverse Mercator with an accuracy of a few
#     nanometers, J. Geodesy 85(8), 475-485 (Aug. 2011);
#     Addenda: https://geographiclib.sourceforge.io/tm-addenda.html
#
#   This extends the series given by Krueger (1912) to sixth order in the
#   flattening. In particular the errors in the projection
#   are less than 5 nanometers within 3900 km of the central meridian (and
#   less than 1 mm within 7600 km of the central meridian).  The mapping
#   can be continued accurately over the poles to the opposite meridian.
#
# Copyright (c) Charles Karney (2012-2022) <charles@karney.com>.

    if  !always_xy
        (lat, lon) = (lon, lat)
    end

    if isa(lon, AbstractFloat)
        lon = [lon]
        lat = [lat]
    end

    T = eltype(lon)
    oneT = one(T)
    zeroT = zero(T)

    lat = copy(lat)
    lon = copy(lon)

    if isa(lat, AbstractArray)
        lat = vec(lat)
        lon = vec(lon)
    end

    # parameters are wrapped in dispatch funciton for type stability
    d2r, p, lon0, lat0, e2, f, e2m, e2, cc, n, b1, a1 = 
        _tranmerc_parameters(lon[1], lat0, lon0, ellips.e, ellips.a)

    e = convert(T, ellips.e)

    if lat0 == 0
        y0 = zero(T)
    else
        sbet0, cbet0 = norm2((one(T) - f) * sin(lat0 * d2r), cose(lat0 * d2r))
        y0 = a1 * (atan(sbet0, cbet0) +
                   SinCosSeries(true, sbet0, cbet0, convert.(T, C1f(n))))
    end

    k = Vector{T}(undef, length(lon))
    gam = Vector{T}(undef, length(lon))
    x = Vector{T}(undef, length(lon))
    y = Vector{T}(undef, length(lon))
    latsign = Vector{T}(undef, length(lon))
    lonsign = Vector{T}(undef, length(lon))
    backside = Vector{Bool}(undef, length(lon))
    xip = Vector{T}(undef, length(lon))
    etap = Vector{T}(undef, length(lon))
    xi = Vector{T}(undef, length(lon))
    eta = Vector{T}(undef, length(lon))

    fmin = sqrt(floatmin(lon0))
    twoT = convert(T, 2)
    alp = convert.(T, alpf(n))
    T90 = convert(T, 90)
    T180 = convert(T, 180)
    T360 = convert(T, 360)

    lon00 = rem(-lon0, T360)
    ind = lon00 < (-T360 / twoT)
    lon00 = (lon00 + T360) * ind + lon00 * !ind
    ind = lon00 > (T360 / twoT)
    lon00 = (lon00 - T360) * ind + lon00 * !ind

    # the specific implementation of this code is very deliberate to maximize the 
    # performance provided by, and to work within the limits of, 
    # LoopVectorization.jl v0.12.159. 
    
    if threaded && Threads.nthreads()>1
        @turbo thread = true for i = eachindex(lon)
            b_1 = rem(lon[i], T360)
            ind_1 = b_1 < (-T360 / twoT)
            b_2 = (b_1 + T360) * ind_1 + b_1 * !ind_1
            ind_2 = b_2 > (T360 / twoT)
            b_3 = (b_2 - T360) * ind_2 + b_2 * !ind_2

            s_1 = lon00 + b_3
            up_1 = s_1 - b_3
            vpp_1 = s_1 - up_1
            up_1 -= lon00
            t_1 = b_3 - vpp_1 - up_1
            
            u_1 = rem(s_1, T360)
            ind_3 = u_1 < (-T360 / twoT)
            u_2 = (u_1 + T360) * ind_3 + u_1 * !ind_3
            ind_4 = u_2 > (T360 / twoT)
            u_3 = (u_2 - T360) * ind_4 + u_2 * !ind_4

            s_2 = u_3 + t_1
            up_2 = s_2 - t_1
            vpp_2 = s_2 - up_2
            up_2 -= u_3
            t_1 -= (vpp_2 + up_2)

            l = ((s_2 == 0) | (abs(s_2) == T180))
            z_1 = (lon[i] - lon0) * l
            ll = (t_1 * l) != 0
            z_2 = -t_1 * ll + z_1 * !ll

            lon[i] = (copysign(s_2, z_2)*l + (s_2 * !l))

            latsign[i] = sign(lat[i])
            lonsign[i] = sign(lon[i])

            lon[i] = lon[i] * lonsign[i] 
            lat[i] = lat[i] * latsign[i]

            backside[i] = lon[i] > T90
            bs = backside[i]>0.5
            b_4 = (bs & (lat[i] == 0))
            latsign[i] = (-oneT * b_4) + (latsign[i] * !b_4)
            lon[i] = (T180 - lon[i]) * bs + lon[i] * !bs

            slam = sin(lon[i] * d2r) 
            clam = cos(lon[i] * d2r) 

            tau = sin(lat[i] * d2r) / max(fmin, cos(lat[i] * d2r))
            
            tau1 = sqrt(tau^2 + oneT)
            sig = sinh(e * atanh(e * tau / tau1))
            taup = sqrt(sig^2 + oneT) * tau - sig * tau1

            h1t = sqrt(oneT + tau^2)
            htc = sqrt(taup^2 + clam^2) 

            xip[i] = atan(taup, clam) 
            etap[i] = asinh(slam / htc) 

            gam[i] = atan(slam * taup, clam * h1t) / d2r
            k[i] = sqrt(e2m + e2 * cos(lat[i] * d2r)^2) * h1t / htc 

            c_3 = lat[i] == T90
            xip[i] = p / twoT * c_3 + xip[i] * !c_3
            etap[i] = zeroT * c_3 + etap[i] * !c_3
            gam[i] = lon[i] * c_3 + gam[i] * !c_3
            k[i] = cc * c_3 + k[i] * !c_3
        end

        # splitting function in 2 greatly improves compile time 
        @turbo thread = true for i = eachindex(lon)
            c0 = cos(twoT * xip[i]) # turbo
            ch0 = cosh(twoT * etap[i]) # turbo
            s0 = sin(twoT * xip[i]) # turbo
            sh0 = sinh(twoT * etap[i]) # turbo

            ar = twoT * c0 * ch0
            ai = twoT * -s0 * sh0

            # --- j = 6 ------
            y1r_6 = alp[6]
            y1i_6 = zeroT

            z1r_6 = twoT * 6 * alp[6]
            z1i_6 = zeroT

            y0r_6 = ar * y1r_6 - ai * y1i_6 + alp[5]
            y0i_6 = ar * y1i_6 + ai * y1r_6

            z0r_6 = ar * z1r_6 - ai * z1i_6 + twoT * 5 * alp[5]
            z0i_6 = ar * z1i_6 + ai * z1r_6

            # --- j = 4 ------
            y1r_4 = ar * y0r_6 - ai * y0i_6 - y1r_6 + alp[4]
            y1i_4 = ar * y0i_6 + ai * y0r_6 - y1i_6

            z1r_4 = ar * z0r_6 - ai * z0i_6 - z1r_6 + twoT * 4 * alp[4]
            z1i_4 = ar * z0i_6 + ai * z0r_6 - z1i_6

            y0r_4 = ar * y1r_4 - ai * y1i_4 - y0r_6 + alp[3]
            y0i_4 = ar * y1i_4 + ai * y1r_4 - y0i_6

            z0r_4 = ar * z1r_4 - ai * z1i_4 - z0r_6 + twoT * 3 * alp[3]
            z0i_4 = ar * z1i_4 + ai * z1r_4 - z0i_6

            # --- j = 2 ------
            y1r_2 = ar * y0r_4 - ai * y0i_4 - y1r_4 + alp[2]
            y1i_2 = ar * y0i_4 + ai * y0r_4 - y1i_4

            z1r_2 = ar * z0r_4 - ai * z0i_4 - z1r_4 + twoT * 2 * alp[2]
            z1i_2 = ar * z0i_4 + ai * z0r_4 - z1i_4

            y0r = ar * y1r_2 - ai * y1i_2 - y0r_4 + alp[1]
            y0i = ar * y1i_2 + ai * y1r_2 - y0i_4

            z0r_2 = ar * z1r_2 - ai * z1i_2 - z0r_4 + twoT * alp[1]
            z0i_2 = ar * z1i_2 + ai * z1r_2 - z0i_4

            z1r = oneT - z1r_2 + z0r_2 * ar / twoT - z0i_2 * ai / twoT
            z1i = -z1i_2 + z0r_2 * ai / twoT + z0i_2 * ar / twoT

            xi[i] = xip[i] + y0r * s0 * ch0 - y0i * c0 * sh0
            eta[i] = etap[i] + y0r * c0 * sh0 + y0i * s0 * ch0

            gam[i] -= -atan(z1i, z1r) / d2r
            k[i] *= (b1 * sqrt(z1r^2 + z1i^2))


            bs = backside[i] > 0.5
            xi[i] = (p - xi[i]) * bs + (xi[i] * !bs)

            y[i] = a1 * xi[i] * latsign[i]
            x[i] = a1 * eta[i] * lonsign[i]

            gam[i] = ((T180 - gam[i]) * bs) + (gam[i] * !bs)

            yx_1 = rem(gam[i] * latsign[i] * lonsign[i], T360)
            ind_lt = yx_1 < (-T360 / twoT)
            yx_2 = (yx_1 + T360) * ind_lt + yx_1 * !ind_lt
            ind_gt = yx_2 > (T360 / twoT)
            yx = (yx_2 - T360) * ind_gt + (yx_2 * !ind_gt)
            
            ind_5 = yx == -T180
            gam[i] = (-yx * ind_5) + (yx * !ind_5)
            y[i] -= y0
        end
    else 
        @turbo thread = false for i = eachindex(lon)
            b_1 = rem(lon[i], T360)
            ind_1 = b_1 < (-T360 / twoT)
            b_2 = (b_1 + T360) * ind_1 + b_1 * !ind_1
            ind_2 = b_2 > (T360 / twoT)
            b_3 = (b_2 - T360) * ind_2 + b_2 * !ind_2

            s_1 = lon00 + b_3
            up_1 = s_1 - b_3
            vpp_1 = s_1 - up_1
            up_1 -= lon00
            t_1 = b_3 - vpp_1 - up_1
            
            u_1 = rem(s_1, T360)
            ind_3 = u_1 < (-T360 / twoT)
            u_2 = (u_1 + T360) * ind_3 + u_1 * !ind_3
            ind_4 = u_2 > (T360 / twoT)
            u_3 = (u_2 - T360) * ind_4 + u_2 * !ind_4

            s_2 = u_3 + t_1
            up_2 = s_2 - t_1
            vpp_2 = s_2 - up_2
            up_2 -= u_3
            t_1 -= (vpp_2 + up_2)

            l = ((s_2 == 0) | (abs(s_2) == T180))
            z_1 = (lon[i] - lon0) * l
            ll = (t_1 * l) != 0
            z_2 = -t_1 * ll + z_1 * !ll

            lon[i] = (copysign(s_2, z_2)*l + (s_2 * !l))

            latsign[i] = sign(lat[i])
            lonsign[i] = sign(lon[i])

            lon[i] = lon[i] * lonsign[i] 
            lat[i] = lat[i] * latsign[i]

            backside[i] = lon[i] > T90
            bs = backside[i]>0.5
            b_4 = (bs & (lat[i] == 0))
            latsign[i] = (-oneT * b_4) + (latsign[i] * !b_4)
            lon[i] = (T180 - lon[i]) * bs + lon[i] * !bs

            slam = sin(lon[i] * d2r) 
            clam = cos(lon[i] * d2r) 

            tau = sin(lat[i] * d2r) / max(fmin, cos(lat[i] * d2r))
            
            tau1 = sqrt(tau^2 + oneT)
            sig = sinh(e * atanh(e * tau / tau1))
            taup = sqrt(sig^2 + oneT) * tau - sig * tau1

            h1t = sqrt(oneT + tau^2)
            htc = sqrt(taup^2 + clam^2) 

            xip[i] = atan(taup, clam) 
            etap[i] = asinh(slam / htc) 

            gam[i] = atan(slam * taup, clam * h1t) / d2r
            k[i] = sqrt(e2m + e2 * cos(lat[i] * d2r)^2) * h1t / htc 

            c_3 = lat[i] == T90
            xip[i] = p / twoT * c_3 + xip[i] * !c_3
            etap[i] = zeroT * c_3 + etap[i] * !c_3
            gam[i] = lon[i] * c_3 + gam[i] * !c_3
            k[i] = cc * c_3 + k[i] * !c_3
        end

        # splitting function in 2 greatly improves compile time 
        @turbo thread = false for i = eachindex(lon)
            c0 = cos(twoT * xip[i]) # turbo
            ch0 = cosh(twoT * etap[i]) # turbo
            s0 = sin(twoT * xip[i]) # turbo
            sh0 = sinh(twoT * etap[i]) # turbo

            ar = twoT * c0 * ch0
            ai = twoT * -s0 * sh0

            # --- j = 6 ------
            y1r_6 = alp[6]
            y1i_6 = zeroT

            z1r_6 = twoT * 6 * alp[6]
            z1i_6 = zeroT

            y0r_6 = ar * y1r_6 - ai * y1i_6 + alp[5]
            y0i_6 = ar * y1i_6 + ai * y1r_6

            z0r_6 = ar * z1r_6 - ai * z1i_6 + twoT * 5 * alp[5]
            z0i_6 = ar * z1i_6 + ai * z1r_6

            # --- j = 4 ------
            y1r_4 = ar * y0r_6 - ai * y0i_6 - y1r_6 + alp[4]
            y1i_4 = ar * y0i_6 + ai * y0r_6 - y1i_6

            z1r_4 = ar * z0r_6 - ai * z0i_6 - z1r_6 + twoT * 4 * alp[4]
            z1i_4 = ar * z0i_6 + ai * z0r_6 - z1i_6

            y0r_4 = ar * y1r_4 - ai * y1i_4 - y0r_6 + alp[3]
            y0i_4 = ar * y1i_4 + ai * y1r_4 - y0i_6

            z0r_4 = ar * z1r_4 - ai * z1i_4 - z0r_6 + twoT * 3 * alp[3]
            z0i_4 = ar * z1i_4 + ai * z1r_4 - z0i_6

            # --- j = 2 ------
            y1r_2 = ar * y0r_4 - ai * y0i_4 - y1r_4 + alp[2]
            y1i_2 = ar * y0i_4 + ai * y0r_4 - y1i_4

            z1r_2 = ar * z0r_4 - ai * z0i_4 - z1r_4 + twoT * 2 * alp[2]
            z1i_2 = ar * z0i_4 + ai * z0r_4 - z1i_4

            y0r = ar * y1r_2 - ai * y1i_2 - y0r_4 + alp[1]
            y0i = ar * y1i_2 + ai * y1r_2 - y0i_4

            z0r_2 = ar * z1r_2 - ai * z1i_2 - z0r_4 + twoT * alp[1]
            z0i_2 = ar * z1i_2 + ai * z1r_2 - z0i_4

            z1r = oneT - z1r_2 + z0r_2 * ar / twoT - z0i_2 * ai / twoT
            z1i = -z1i_2 + z0r_2 * ai / twoT + z0i_2 * ar / twoT

            xi[i] = xip[i] + y0r * s0 * ch0 - y0i * c0 * sh0
            eta[i] = etap[i] + y0r * c0 * sh0 + y0i * s0 * ch0

            gam[i] -= -atan(z1i, z1r) / d2r
            k[i] *= (b1 * sqrt(z1r^2 + z1i^2))


            bs = backside[i] > 0.5
            xi[i] = (p - xi[i]) * bs + (xi[i] * !bs)

            y[i] = a1 * xi[i] * latsign[i]
            x[i] = a1 * eta[i] * lonsign[i]

            gam[i] = ((T180 - gam[i]) * bs) + (gam[i] * !bs)

            yx_1 = rem(gam[i] * latsign[i] * lonsign[i], T360)
            ind_lt = yx_1 < (-T360 / twoT)
            yx_2 = (yx_1 + T360) * ind_lt + yx_1 * !ind_lt
            ind_gt = yx_2 > (T360 / twoT)
            yx = (yx_2 - T360) * ind_gt + (yx_2 * !ind_gt)
            
            ind_5 = yx == -T180
            gam[i] = (-yx * ind_5) + (yx * !ind_5)
            y[i] -= y0
        end
    end
    
    if length(x) == 1;
        x = x[1]
        y = y[1]
        gam = gam[1]
        k = k[1]
    end

    return x, y, gam, k
end


"""
    tranmerc_inv(x, y; lon0, lat0, ellips, threaded, always_xy

Returns of longitude and latitude in geodetic coordinates (EPSG:4326) coordinates [decimal 
degrees] given x and y coodinates [meteres] in transverse Mercator (TM) projection. The 
TM projection is defined by kwargs of longitude (lon0) and latitude (lat0), which specify 
the center of the  projeciton, and an ellipsoid that is define by an equatorial radius 
[meters] (a) and its  eccentricity (e). Also returnes meridian convergence (gam) and scale 
factor (k).
"""

function tranmerc_inv(
    x::Union{AbstractVector{<:AbstractFloat},AbstractFloat}, 
    y::Union{AbstractVector{<:AbstractFloat},AbstractFloat}; 
    lon0::Real=0, 
    lat0::Real=0, 
    ellips::Ellipsoid=ellipsoid(EPSG(7030)), 
    threaded=true, 
    always_xy=true
    )

    # This code has been heavily modified for maximum preformance in Julia from the MATLAB 
    # implimentation of geographiclib_toolbox-2.0 by Alex Gardner JPL/NASA.
    #
    # This implementation of the projection is based on the series method
    #   described in
    #
    #     C. F. F. Karney, Transverse Mercator with an accuracy of a few
    #     nanometers, J. Geodesy 85(8), 475-485 (Aug. 2011);
    #     Addenda: https://geographiclib.sourceforge.io/tm-addenda.html
    #
    #   This extends the series given by Krueger (1912) to sixth order in the
    #   flattening. In particular the errors in the projection
    #   are less than 5 nanometers within 3900 km of the central meridian (and
    #   less than 1 mm within 7600 km of the central meridian).  The mapping
    #   can be continued accurately over the poles to the opposite meridian.
    #
    # Copyright (c) Charles Karney (2012-2022) <charles@karney.com>.

    T = eltype(x)

    if isa(x, AbstractFloat)
        x = [x]
        y = [y]
    end

    # parameters are wrapped in dispatch funciton for type stability
    d2r, p, lon0, lat0, e2, f, e2m, e2, cc, n, b1, a1 =
        _tranmerc_parameters(x[1], lat0, lon0, ellips.e, ellips.a)
    
    e = convert(T, ellips.e)
    if lat0 == 0
        y0 = 0
    else
        sbet0, cbet0 = norm2((one(T) .- f) .* sin(lat0 * d2r), cos(lat0 * d2r))
        y0 = a1 * (atan(sbet0, cbet0) + 
                    SinCosSeries(true, sbet0, cbet0,  convert(T,C1f(n))))
    end
    
    p2 = p / 2
    bet = convert.(T, betf(n))
    lon00 = AngNormalize(lon0)

    gam = Vector{T}(undef, length(x))
    k = Vector{T}(undef, length(x))
    lat = Vector{T}(undef, length(x))
    lon = Vector{T}(undef, length(x))

    xip = Vector{T}(undef, length(x))
    etap = Vector{T}(undef, length(x))
    gam = Vector{T}(undef, length(x))
    tau = Vector{T}(undef, length(x))
    k = Vector{T}(undef, length(x))
    xisign = Vector{T}(undef, length(x))
    etasign = Vector{T}(undef, length(x))
    backside = Vector{Bool}(undef, length(x))

    twoT = convert(T, 2)
    oneT = one(T)
    zeroT = zero(T)
    T90 = convert(T, 90)
    T180 = convert(T, 180)
    T360 = convert(T, 360)

    # the specific implementation of this code is very deliberate to maximize the 
    # performance provided by and work within the limits of LoopVectorization.jl v0.12.159. 

    if threaded && Threads.nthreads()>1
        @tturbo thread = true for i = eachindex(y)
            yA = y[i] + y0
            xiA = yA / a1
            etaA = x[i] / a1
            xisign[i] = sign(xiA)
            etasign[i] = sign(etaA)
            xiB = xiA * xisign[i]
            eta = etaA * etasign[i]
            backside[i] = xiB > p2

            bs = backside[i]>0.5
            xi = ((p - xiB) * bs) + (xiB * !bs)

            c0 = cos(twoT * xi)
            ch0 = cosh(twoT * eta)
            s0 = sin(twoT * xi)
            sh0 = sinh(twoT * eta)

            ar = twoT * c0 * ch0
            ai = twoT * -s0 * sh0

            # --- j = 6 ------
            y1r_6 = - bet[6] 
            y1i_6 = zeroT

            z1r_6 =  -twoT * 6 * bet[6] 
            z1i_6 = zeroT

            y0r_6 = ar * y1r_6 - ai * y1i_6 - bet[6-1]
            y0i_6 = ar * y1i_6 + ai * y1r_6

            z0r_6 = ar * z1r_6 - ai * z1i_6 - twoT * (6 - 1) * bet[6-1]
            z0i_6 = ar * z1i_6 + ai * z1r_6

            # --- j = 4 ------
            y1r_4 = ar * y0r_6 - ai * y0i_6 - y1r_6 - bet[4]
            y1i_4 = ar * y0i_6 + ai * y0r_6 - y1i_6

            z1r_4 = ar * z0r_6 - ai * z0i_6 - z1r_6 - twoT * 4 * bet[4]
            z1i_4 = ar * z0i_6 + ai * z0r_6 - z1i_6

            y0r_4 = ar * y1r_4 - ai * y1i_4 - y0r_6 - bet[4-1]
            y0i_4 = ar * y1i_4 + ai * y1r_4 - y0i_6

            z0r_4 = ar * z1r_4 - ai * z1i_4 - z0r_6 - twoT * (4 - 1) * bet[4-1]
            z0i_4 = ar * z1i_4 + ai * z1r_4 - z0i_6

            # --- j = 2 ------
            y1r_2 = ar * y0r_4 - ai * y0i_4 - y1r_4 - bet[2]
            y1i_2 = ar * y0i_4 + ai * y0r_4 - y1i_4

            z1r_2 = ar * z0r_4 - ai * z0i_4 - z1r_4 - twoT * 2 * bet[2]
            z1i_2 = ar * z0i_4 + ai * z0r_4 - z1i_4

            y0r_2 = ar * y1r_2 - ai * y1i_2 - y0r_4 - bet[2-1]
            y0i_2 = ar * y1i_2 + ai * y1r_2 - y0i_4

            z0r_2 = ar * z1r_2 - ai * z1i_2 - z0r_4 - twoT * (2 - 1) * bet[2-1]
            z0i_2 = ar * z1i_2 + ai * z1r_2 - z0i_4

            z1r = oneT - z1r_2 + z0r_2 * ar / twoT - z0i_2 * ai / twoT
            z1i = -z1i_2 + z0r_2 * ai / twoT + z0i_2 * ar / twoT
            
            xip[i] = xi + y0r_2 * s0 * ch0 - y0i_2 * c0 * sh0
            etap[i] = eta + y0r_2 * c0 * sh0 + y0i_2 * s0 * ch0

            gam[i] = atan(z1i, z1r) / d2r
            k[i] = b1 / sqrt(z1r^2 + z1i^2)

            s = sinh(etap[i])
            c = max(zeroT, cos(xip[i]))
            r = sqrt(s^2 + c^2)
            lon[i] = atan(s, c) / d2r

            sxip = sin(xip[i])
            sr = sxip / r

            gam[i] += atan(sxip * tanh(etap[i]), c) / d2r

            tau[i] = sr / e2m
            
            #numit = 5, iter = 1
            tau1_1 = sqrt(tau[i]^2 + oneT)
            sig_1 = sinh(e * atanh(e * tau[i]/ tau1_1))
            taupa_1 = sqrt(sig_1^2 + oneT) * tau[i] - sig_1 * tau1_1

            tau[i] +=
                (sr - taupa_1) * (oneT + e2m * tau[i]^2) /
                (e2m * tau1_1 * sqrt(taupa_1^2 + oneT))

            #numit = 5, iter = 2
            tau1_2 = sqrt(tau[i]^2 + oneT)
            sig_2 = sinh(e * atanh(e * tau[i] / tau1_2))
            taupa_2 = sqrt(sig_2^2 + oneT) * tau[i] - sig_2 * tau1_2
            tau[i] += (sr - taupa_2) * (oneT + e2m * tau[i]^2) / 
                (e2m * tau1_2 * sqrt(taupa_2^2 + oneT))
        
            #numit = 5, iter = 3
            tau1_3 = sqrt(tau[i]^2 + oneT)
            sig_3 = sinh(e * atanh(e * tau[i] / tau1_3))
            taupa_3 = sqrt(sig_3^2 + oneT) * tau[i] - sig_3 * tau1_3
            tau[i] += (sr - taupa_3) * (oneT + e2m * tau[i]^2) / 
                (e2m * tau1_3 * sqrt(taupa_3^2 + oneT))

            #numit = 5, iter = 4
            tau1_4 = sqrt(tau[i]^2 + oneT)
            sig_4 = sinh(e * atanh(e * tau[i] / tau1_4))
            taupa_4 = sqrt(sig_4^2 + oneT) * tau[i] - sig_4 * tau1_4
            tau[i] +=
                (sr - taupa_4) * (oneT + e2m * tau[i]^2) /
                (e2m * tau1_4 * sqrt(taupa_4^2 + oneT))

            #numit = 5, iter = 5
            tau1_5 = sqrt(tau[i]^2 + oneT)
            sig_5 = sinh(e * atanh(e * tau[i] / tau1_5))
            taupa_5 = sqrt(sig_5^2 + oneT) * tau[i] - sig_5 * tau1_5
            tau[i] +=
                (sr - taupa_5) * (oneT + e2m * tau[i]^2) /
                (e2m * tau1_5 * sqrt(taupa_5^2 + oneT))

            lat[i] = atan(tau[i], oneT) / d2r

            ca = r != zeroT
            k[i] *= ((e2m + e2 / sqrt(oneT + tau[i]^2)) * sqrt(oneT + tau[i]^2) * r * ca) + !ca

            cb = !ca
            lat[i] = T90 * cb + lat[i] * !cb
            lon[i] = zeroT * cb + lon[i] * !cb
            k[i] = k[i] * cc * cb + k[i] * !cb

            lat[i] *= xisign[i]

            bs = backside[i] > 0.5;
            lon[i] = ((T180 - lon[i]) * bs) + (lon[i] * !bs)

            yx_1 = rem(lon[i] * etasign[i] + lon00, T360)
            ind_lt = yx_1 < (-T360 / twoT)
            yx_2 = ((yx_1 + T360) * ind_lt) + (yx_1 * !ind_lt)
            ind_gt = yx_2 > (T360 / twoT)
            yx_3 = ((yx_2 - T360) * ind_gt) + (yx_2 * !ind_gt)

            ind_1 = yx_3 == -T180
            lon[i] = (-yx_3 * ind_1) + (yx_3 * !ind_1)

            bs2 = backside[i] > 0.5
            gam[i] = ((T180 - gam[i]) * bs2) + (gam[i] * !bs2)

            yx_4 = rem(gam[i] * xisign[i] * etasign[i], T360)
            ind_lt2 = yx_4 < (-T360 / twoT)
            yx_5 = ((yx_4 + T360) * ind_lt2) + (yx_4 * !ind_lt2)
            ind_gt2 = yx_5 > (T360 / twoT)
            yx_6 = ((yx_5 - T360) * ind_gt2) + (yx_5 * !ind_gt2)

            ind_2 = yx_6 == -T180
            gam[i] = -yx_6 * ind_2 + yx_6 * !ind_2
        end
    else
        @tturbo thread = false for i = eachindex(y)
            yA = y[i] + y0
            xiA = yA / a1
            etaA = x[i] / a1
            xisign[i] = sign(xiA)
            etasign[i] = sign(etaA)
            xiB = xiA * xisign[i]
            eta = etaA * etasign[i]
            backside[i] = xiB > p2

            bs = backside[i]>0.5
            xi = ((p - xiB) * bs) + (xiB * !bs)

            c0 = cos(twoT * xi)
            ch0 = cosh(twoT * eta)
            s0 = sin(twoT * xi)
            sh0 = sinh(twoT * eta)

            ar = twoT * c0 * ch0
            ai = twoT * -s0 * sh0

            # --- j = 6 ------
            y1r_6 = - bet[6] 
            y1i_6 = zeroT

            z1r_6 =  -twoT * 6 * bet[6] 
            z1i_6 = zeroT

            y0r_6 = ar * y1r_6 - ai * y1i_6 - bet[6-1]
            y0i_6 = ar * y1i_6 + ai * y1r_6

            z0r_6 = ar * z1r_6 - ai * z1i_6 - twoT * (6 - 1) * bet[6-1]
            z0i_6 = ar * z1i_6 + ai * z1r_6

            # --- j = 4 ------
            y1r_4 = ar * y0r_6 - ai * y0i_6 - y1r_6 - bet[4]
            y1i_4 = ar * y0i_6 + ai * y0r_6 - y1i_6

            z1r_4 = ar * z0r_6 - ai * z0i_6 - z1r_6 - twoT * 4 * bet[4]
            z1i_4 = ar * z0i_6 + ai * z0r_6 - z1i_6

            y0r_4 = ar * y1r_4 - ai * y1i_4 - y0r_6 - bet[4-1]
            y0i_4 = ar * y1i_4 + ai * y1r_4 - y0i_6

            z0r_4 = ar * z1r_4 - ai * z1i_4 - z0r_6 - twoT * (4 - 1) * bet[4-1]
            z0i_4 = ar * z1i_4 + ai * z1r_4 - z0i_6

            # --- j = 2 ------
            y1r_2 = ar * y0r_4 - ai * y0i_4 - y1r_4 - bet[2]
            y1i_2 = ar * y0i_4 + ai * y0r_4 - y1i_4

            z1r_2 = ar * z0r_4 - ai * z0i_4 - z1r_4 - twoT * 2 * bet[2]
            z1i_2 = ar * z0i_4 + ai * z0r_4 - z1i_4

            y0r_2 = ar * y1r_2 - ai * y1i_2 - y0r_4 - bet[2-1]
            y0i_2 = ar * y1i_2 + ai * y1r_2 - y0i_4

            z0r_2 = ar * z1r_2 - ai * z1i_2 - z0r_4 - twoT * (2 - 1) * bet[2-1]
            z0i_2 = ar * z1i_2 + ai * z1r_2 - z0i_4

            z1r = oneT - z1r_2 + z0r_2 * ar / twoT - z0i_2 * ai / twoT
            z1i = -z1i_2 + z0r_2 * ai / twoT + z0i_2 * ar / twoT
            
            xip[i] = xi + y0r_2 * s0 * ch0 - y0i_2 * c0 * sh0
            etap[i] = eta + y0r_2 * c0 * sh0 + y0i_2 * s0 * ch0

            gam[i] = atan(z1i, z1r) / d2r
            k[i] = b1 / sqrt(z1r^2 + z1i^2)

            s = sinh(etap[i])
            c = max(zeroT, cos(xip[i]))
            r = sqrt(s^2 + c^2)
            lon[i] = atan(s, c) / d2r

            sxip = sin(xip[i])
            sr = sxip / r

            gam[i] += atan(sxip * tanh(etap[i]), c) / d2r

            tau[i] = sr / e2m
            
            #numit = 5, iter = 1
            tau1_1 = sqrt(tau[i]^2 + oneT)
            sig_1 = sinh(e * atanh(e * tau[i]/ tau1_1))
            taupa_1 = sqrt(sig_1^2 + oneT) * tau[i] - sig_1 * tau1_1

            tau[i] +=
                (sr - taupa_1) * (oneT + e2m * tau[i]^2) /
                (e2m * tau1_1 * sqrt(taupa_1^2 + oneT))

            #numit = 5, iter = 2
            tau1_2 = sqrt(tau[i]^2 + oneT)
            sig_2 = sinh(e * atanh(e * tau[i] / tau1_2))
            taupa_2 = sqrt(sig_2^2 + oneT) * tau[i] - sig_2 * tau1_2
            tau[i] += (sr - taupa_2) * (oneT + e2m * tau[i]^2) / 
                (e2m * tau1_2 * sqrt(taupa_2^2 + oneT))
        
            #numit = 5, iter = 3
            tau1_3 = sqrt(tau[i]^2 + oneT)
            sig_3 = sinh(e * atanh(e * tau[i] / tau1_3))
            taupa_3 = sqrt(sig_3^2 + oneT) * tau[i] - sig_3 * tau1_3
            tau[i] += (sr - taupa_3) * (oneT + e2m * tau[i]^2) / 
                (e2m * tau1_3 * sqrt(taupa_3^2 + oneT))

            #numit = 5, iter = 4
            tau1_4 = sqrt(tau[i]^2 + oneT)
            sig_4 = sinh(e * atanh(e * tau[i] / tau1_4))
            taupa_4 = sqrt(sig_4^2 + oneT) * tau[i] - sig_4 * tau1_4
            tau[i] +=
                (sr - taupa_4) * (oneT + e2m * tau[i]^2) /
                (e2m * tau1_4 * sqrt(taupa_4^2 + oneT))

            #numit = 5, iter = 5
            tau1_5 = sqrt(tau[i]^2 + oneT)
            sig_5 = sinh(e * atanh(e * tau[i] / tau1_5))
            taupa_5 = sqrt(sig_5^2 + oneT) * tau[i] - sig_5 * tau1_5
            tau[i] +=
                (sr - taupa_5) * (oneT + e2m * tau[i]^2) /
                (e2m * tau1_5 * sqrt(taupa_5^2 + oneT))

            lat[i] = atan(tau[i], oneT) / d2r

            ca = r != zeroT
            k[i] *= ((e2m + e2 / sqrt(oneT + tau[i]^2)) * sqrt(oneT + tau[i]^2) * r * ca) + !ca

            cb = !ca
            lat[i] = T90 * cb + lat[i] * !cb
            lon[i] = zeroT * cb + lon[i] * !cb
            k[i] = k[i] * cc * cb + k[i] * !cb

            lat[i] *= xisign[i]

            bs = backside[i] > 0.5;
            lon[i] = ((T180 - lon[i]) * bs) + (lon[i] * !bs)

            yx_1 = rem(lon[i] * etasign[i] + lon00, T360)
            ind_lt = yx_1 < (-T360 / twoT)
            yx_2 = ((yx_1 + T360) * ind_lt) + (yx_1 * !ind_lt)
            ind_gt = yx_2 > (T360 / twoT)
            yx_3 = ((yx_2 - T360) * ind_gt) + (yx_2 * !ind_gt)

            ind_1 = yx_3 == -T180
            lon[i] = (-yx_3 * ind_1) + (yx_3 * !ind_1)

            bs2 = backside[i] > 0.5
            gam[i] = ((T180 - gam[i]) * bs2) + (gam[i] * !bs2)

            yx_4 = rem(gam[i] * xisign[i] * etasign[i], T360)
            ind_lt2 = yx_4 < (-T360 / twoT)
            yx_5 = ((yx_4 + T360) * ind_lt2) + (yx_4 * !ind_lt2)
            ind_gt2 = yx_5 > (T360 / twoT)
            yx_6 = ((yx_5 - T360) * ind_gt2) + (yx_5 * !ind_gt2)

            ind_2 = yx_6 == -T180
            gam[i] = -yx_6 * ind_2 + yx_6 * !ind_2
        end
    end

    if length(lon) == 1
        lon = lon[1]
        lat = lat[1]
        gam = gam[1]
        k = k[1]
    end

    if always_xy
        return lon, lat, gam, k
    else
        return lat, lon, gam, k
    end
end

function _tranmerc_parameters(x::Float32, lat0, lon0, e, a)
    d2r = Float32(pi / 180)
    p = Float32(pi)
    lat0 = Float32(lat0)
    lon0 = Float32(lon0)
    e2 = e^2
    f = e2 / (1 + sqrt(1 - e2))
    e2m = 1 - e2
    e2 = Float32(e2)
    cc = Float32(sqrt(e2m) * exp(e * atanh(e * 1)))
    e2m = Float32(e2m)
    n = f / (2 - f)
    b1 = Float32((1 - f) * (A1m1f(n) + 1))
    a1 = Float32(b1 * a)

    return d2r, p, lon0, lat0, e2, f, e2m, e2, cc, n, b1, a1
end

function _tranmerc_parameters(x::Float64, lat0, lon0, e, a)
    d2r = Float64(pi / 180)
    p = Float64(pi)
    lat0 = Float64(lat0)
    lon0 = Float64(lon0)
    e2 = e^2
    f = e2 / (1 + sqrt(1 - e2))
    e2m = 1 - e2
    e2 = Float64(e2)
    cc = Float64(sqrt(e2m) * exp(e * atanh(e * 1)))
    e2m = Float64(e2m)
    n = f / (2 - f)
    b1 = Float64((1 - f) * (A1m1f(n) + 1))
    a1 = Float64(b1 * a)

    return d2r, p, lon0, lat0, e2, f, e2m, e2, cc, n, b1, a1
end

function alpf(n)

    alpcoeff = [
        31564, -66675, 34440, 47250, -100800, 75600, 151200,
        -1983433, 863232, 748608, -1161216, 524160, 1935360,
        670412, 406647, -533952, 184464, 725760,
        6601661, -7732800, 2230245, 7257600,
        -13675556, 3438171, 7983360,
        212378941, 319334400,
    ]
    maxpow = 6
    alp = zeros(maxpow)
    o = 1
    d = n
    pwr = 0:length(alpcoeff)

    for l = 1:maxpow
        m = maxpow - l
        coeff = reverse((alpcoeff[o:(o+m)]));
        poly = sum(coeff .* n .^ pwr[1:length(coeff)])
        alp[l] = d * poly / alpcoeff[o+m+1]

        o = o + m + 2
        d = d * n
    end
    return alp
end

function betf(n)
    betcoeff = [
        384796, -382725, -6720, 932400, -1612800, 1209600, 2419200,
        -1118711, 1695744, -1174656, 258048, 80640, 3870720,
        22276, -16929, -15984, 12852, 362880,
        -830251, -158400, 197865, 7257600,
        -435388, 453717, 15966720,
        20648693, 638668800
    ]
    maxpow = 6
    bet = zeros(maxpow)
    o = 1
    d = n
    pwr = 0:length(betcoeff)

    for l = 1:maxpow
        m = maxpow - l

        coeff = reverse((betcoeff[o:(o+m)]))
        poly = sum(coeff .* n .^ pwr[1:length(coeff)])
        bet[l] = d * poly / betcoeff[o+m+1]

        o += m + 2
        d *= n
    end
    return bet
end

function A1m1f(epsi)
    # A1M1F  Evaluate A_1 - 1
    #
    #   A1m1 = A1M1F(epsi) evaluates A_1 - 1 using Eq. (17).  epsi and A1m1 are
    #   K x 1 arrays.

    eps2 = epsi^2
    coeff = [0, 64, 4, 1]
    pwr = 0:(length(coeff)-1)
    t = sum(coeff .* eps2 .^ pwr) / 256
    A1m1 = (t + epsi) / (1 - epsi)

    return A1m1
end

function AngNormalize(x)
#ANGNORMALIZE  Reduce angle to range (-180, 180]
#
#   x = ANGNORMALIZE(x) reduces angles to the range (-180, 180].  x can be
#   any shape.

    y = remx.(x, convert(eltype(x), 360))
    ind = y .== -180
    if any(ind)
        y[ind] .*= -one(eltype(x))
    end

    return y
end


function remx(x, y)
    #REMX   The remainder function
    #
    #   REMX(x, y) is the remainder of x on division by y.  Result is in [-y/2,
    #   y/2].  x can be compatible shapes.  y should be a  positive scalar.

    z = rem(x, y);
    if z < -y/2
        z += y
    elseif z > y/2
        z -= y
    end
  return z
end