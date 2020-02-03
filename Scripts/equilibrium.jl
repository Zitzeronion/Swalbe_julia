"""
Script for the implementation of the equilibrium distribution functions.
The two arrays for the weights and the lattice speeds are defined here 
and can be accessed by ever function following.
"""

# D2Q9 lattice weights 
Weightsd2q9 = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36]
# D2Q9 lattice velocity vector: (0,0) (1,0) (0,1) ....
Cd2q9 = [0 1 0 -1 0 1 -1 -1 1; 0 0 1 0 -1 1 1 -1 -1]
# D1Q3 lattice weights
Weightsd1q3 = [2/3 1/6 1/6]
# D1Q3 lattice velocity
Cd1q3 = [0 1 -1]

"""
    equiliriumnogd2q9(height, velocity_x, velocity_y)

Computes the `equilibrium distribution function` for D2Q9 shallow water lattice Boltzmann
based on `height`, `velocity_x` and `velocity_y` fields but without gravity.
See also: [`equiliriumd2q9`](@ref), [`equiliriumd1q3`](@ref).

# Arguments
- `height::Array{Float64,2}` : `height` at every lattice node
- `velocity_x::Array{Float64,2}` : `x`-component of the macroscopic velocity
- `velocity_y::Array{Float64,2}` : `y`-component of the macroscopic velocity

# Examples
```jldoctest
julia> A = [1 1; 1 1]; 
julia> B = [0 0; 0 0];                  # example without velocity
julia> C = equiliriumnogd2q9(A,B,B); 
julia> println(C[1, :, :])
[1.0 1.0; 1.0 1.0]
```
"""
function equiliriumnogd2q9(height, velocity_x, velocity_y)

    inputsize = size(height)
    
    Equi_dist = zeros(9, inputsize[1], inputsize[2])
    # The only special case is the one of the zero velocity component
    Equi_dist[1,:,:] = Weightsd2q9[1] .* height .*(9/4 * ones(inputsize) - 3/2 * (velocity_x.*velocity_x .+ velocity_y.*velocity_y))
    # All other 8 components are calculated with this formular
    for lat_vel = 2:9
        Equi_dist[lat_vel, :, :] = Weightsd2q9[lat_vel] .* height 
                                    .*(3 * (Cd2q9[1,lat_vel] .* velocity_x + Cd2q9[2,lat_vel] .* velocity_y) 
                                        + 9/2 * (Cd2q9[1,lat_vel] .* velocity_x .* Cd2q9[2,lat_vel] .* velocity_y)
                                        - 3/2 * (velocity_x.*velocity_x .+ velocity_y.*velocity_y))
    end

    return Equi_dist
end

"""
    equiliriumd2q9(height, velocity_x, velocity_y, gravity)

Computes the `equilibrium distribution function` for D2Q9 shallow water lattice Boltzmann
based on `height`, `velocity_x` and `velocity_y` fields with gravity.
See also: [`equiliriumnogd2q9`](@ref), [`equiliriumd1q3`](@ref), doi:https://doi.org/10.1357/002224099764805174.

# Arguments
- `height::Array{Float64,2}` : `height` at every lattice node
- `velocity_x::Array{Float64,2}` : `x`-component of the macroscopic velocity
- `velocity_y::Array{Float64,2}` : `y`-component of the macroscopic velocity
- `gravity::Float64` : Acceleration force due to `gravity`, needed if there is a `bed slope`

# Examples
```jldoctest
julia> A = [1 1; 1 1]; 
julia> B = [0 0; 0 0];                  # example without velocity
julia> C = equiliriumd2q9(A,B,B,0.1); 
julia> println(C[1, :, :])
[0.916667 0.916667; 0.916667 0.916667]
```
"""
function equiliriumd2q9(height, velocity_x, velocity_y, gravity)
    # Set up the output array
    inputsize = size(height)
    Equi_dist = Array{Float64,3}(undef, 9, inputsize[1], inputsize[2])

    # The only special case is the one of the zero velocity component
    Equi_dist[1,:,:] = Weightsd2q9[1] .* height .* (9/4 * ones(inputsize) - 15/8 * gravity * height - 
                                                    3/2 * (velocity_x.*velocity_x .+ velocity_y.*velocity_y))
    # All other 8 components are calculated with this formular
    for lat_vel  = 2:9
        Equi_dist[lat_vel, :, :] = Weightsd2q9[lat_vel] .* height .* (3/2 * gravity * height 
                                                                + 3 * (c[1,lat_vel] .* velocity_x + c[2,lat_vel] .* velocity_y) 
                                                                + 9/2 * (c[1,lat_vel] .* velocity_x .* c[2,lat_vel] .* velocity_y)
                                                                - 3/2 * (velocity_x.*velocity_x .+ velocity_y.*velocity_y))
    end

    return Equi_dist
end

"""
    equiliriumd1q3(height, velocity, gravity)

Computes the `equilibrium distribution function` for D1Q3 shallow water lattice Boltzmann
based on `height` and `velocity` fields with gravity.
See also: [`equiliriumnogd2q9`](@ref), [`equiliriumd2q9`](@ref), doi:10.1016/j.jcp.2010.06.022.

# Arguments
- `height::Array{Float64,1}` : `height` at every lattice node
- `velocity::Array{Float64,1}` : The macroscopic `velocity`
- `gravity::Float64` : Acceleration force due to `gravity`, needed if there is a `bed slope`

# Examples
```jldoctest
julia> A = [1 1 1 1]; 
julia> B = [0 0 0 0];                  # example without velocity
julia> C = equiliriumd1q3(A,B,0.1); 
julia> println(C[1, :])
[0.916667 0.916667 0.916667 0.916667]
```
"""
function equiliriumd1q3(height, velocity, gravity)
    # Set up the output array based on input
    inputsize = size(height)
    Equi_dist = Array{Float64,2}(undef, 3, inputsize[2])

    # Components of the equilibrium distribution
    Equi_dist[1,:] = height - 1/2 * gravity * height.^2 - height .* velocity .* velocity
    Equi_dist[2,:] = 1/4 * gravity * height.^2 + 1/2 * Cd1q3[2] * height .* velocity + 1/2 * height .* velocity .* velocity
    Equi_dist[3,:] = 1/4 * gravity * height.^2 + 1/2 * Cd1q3[3] * height .* velocity + 1/2 * height .* velocity .* velocity

    return Equi_dist
end