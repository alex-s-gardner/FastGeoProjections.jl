module FastGeoProjections
    
    using Proj # Proj dependancy included untill package is more mature
  
    include("polarstereo.jl")
    export polarstereo_fwd
    export polarstereo_inv

end
