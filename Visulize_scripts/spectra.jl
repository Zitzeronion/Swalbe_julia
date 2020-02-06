using Plots, FFTW

"""
    omega(q, h0=1, q0=2.4)

Computes the spectrum omega of thin film with average height h0 and critical wavenumber q0.

# Agurments
- `q::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}` : Range for the wavenumbers to compute the specturm
- `h0::Float64` : Initial mean height value
- `q0::Float64` : Critical wave number based on interfacial potential

# Examples
```jldoctest
julia> A = [1 1 1 1]; 
julia> omega(A)             # call without optional arguments
3.50667  3.50667  3.50667  3.50667
```
"""
function omega(q, h0=1.0, q0=2.4)
    dispersion = h0^3 * q .*(2*q0^2 .- q.^2)/3.0
    return dispersion
end


"""
    initialstructure(L)

Computes the fourier transform of an initial height configuration.

# Agurments
- `L::Int64` : Length of the one dimensional height configuration

# Examples
```jldoctest
julia> L = 10
julia> initialstructure(L)
([3.50667  3.50667  3.50667  3.50667 ...], -L/2:1.0:(L/2-1)
```
"""
function initialstructure(L)
    h = ones(L)
    noise = rand(L)*0.001

    hint = h + noise
    ffth = fft(hint)
    ffth[1] = 0.001
    structure = fftshift(ffth)
    powerspectrum = (abs.(structure)).^2

    df = -L/2:1:(L/2-1)

    return powerspectrum, df
end
   
function spectrumdeterministic(L, N, t)
    powerspec, df = initialstructure(L)
    q = range(0.0001, stop = N*2.4, length = L)
    spectrum = powerspec .* exp.(2*omega(q)*t) 
    return spectrum, q
end

function justexponental(q)
    y = [exp.(2*omega(q)*2) exp.(2*omega(q)*2) + q.^2 ./omega(q).*(exp.(2*omega(q)*2).-1)]
    return y
end