# FastGeoProjections

[![Build Status](https://github.com/alex-s-gardner/FastGeoProjections.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/alex-s-gardner/FastGeoProjections.jl/actions/workflows/CI.yml?query=branch%3Amain)

!! UNDER DEVELOPMENT BY OVERCOMMITTED AND UNDER PAID DEVELOPERS !!

**FastGeoProjections** is intended to provide highly optimized native Julia geospatial coordinate transformations from one coordinate reference system (CRS) to another as defined by EPSG codes. It is not intended to to replace or to be as comprehensive as [Proj](https://github.com/JuliaGeo/Proj.jl). The package will natively support only the most common geospatial transformations and relies on **Proj.jl** for all others.