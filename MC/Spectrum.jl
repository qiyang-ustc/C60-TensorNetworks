mutable struct Spectrum
    length::Int
    count ::Array{Int128,1}
    factor ::Array{Float64,1}
    level ::Array{Int,1}
    ground ::Int
    function Spectrum(length::Int)
        count = zeros(Int128,length)
        factor = zeros(Float64,length)
        level = zeros(Float64,length)
        ground = 0
        new(length,count,factor,level,ground)
    end
end

function set_ground!(spectrum::Spectrum,g::Int)
    spectrum.ground = g
    return 0
end
function add_count!(spectrum::Spectrum,energy::Int)
    spectrum.count[-spectrum.ground+energy]+=1
    return 0
end

function get_count(spectrum::Spectrum,energy::Int)
    return spectrum.count[-spectrum.ground+energy]
end
function add_factor!(spectrum::Spectrum,energy::Int,temp::Float64)
    return spectrum.factor[-spectrum.ground+energy]+=temp
end 

function get_factor(spectrum::Spectrum,energy::Int)
    return spectrum.factor[-spectrum.ground+energy]
end 

function normalize_factor(spectrum::Spectrum)
    max_factor = maximum(spectrum.factor)
    spectrum.factor = spectrum.factor .- max_factor
    return 1
end
