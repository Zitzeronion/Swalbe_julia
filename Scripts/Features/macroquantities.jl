"""
Computation of macroscopic quantities.
For the shallow water lattice Boltzmann method the sum of the distribution function is the height.
The sum of distribution function times the lattice velocity is the macroscopic velocity.
Theoretically this quantities are the moments of the distribution function.
"""


"""
    momentsd2q9(height, velocity_x, velocity_y, forces_x, forces_y, distribution, kindofforce="WFM")

Computes the macroscopic quantities `height` and `velocity` from a `distribution function` and if used `forces`.
See also: [`collideSRTd2q9`](@ref)

# Arguments
- `height::Array{Float64,2}` : `height`-value at every lattice node, similar to density in standard lattice Boltzmann
- `velocity_x::Array{Float64,2}` : The `x`-component of the macroscopic velocity
- `velocity_y::Array{Float64,2}` : The `y`-component of the macroscopic velocity
- `forces_x::Array{Float64,2}` : X-component of all applied `forces`
- `forces_y::Array{Float64,2}` : Y-component of all applied `forces`
- `distribution::Array{Float64,3}` : Internal distribution function used for the lattice Boltzmann method
- `kindofforce::String` : Determines which forcing model is used, higher order forcing are included in the moment calculation

# Examples
```jldoctest
julia> height = [1.0 1.0; 1.0 1.0]; 
julia> velo = [0.0 0.0; 0.0 0.0];                           # example without velocity
julia> distribution = fill(1.0, (9, 2, 2))
julia> h_out,velx_out,vely_out = momentsd2q9(height, velo, velo, distribution, "Guo"); 
julia> println(h_out)
[9.0 9.0; 9.0 9.0]
```

# Publication
- `The lattice Boltzmann method as a basis for ocean circulation modeling` : Rick Salmon, doi:10.1357/002224099764805174.
- `An evaluation of force terms in the lattice Boltzmann models in simulating shallow water flows over complex topography` : Weighting factor mode, doi:10.1002/fld.4726

"""
function momentsd2q9(height, velocity_x, velocity_y, forces_x, forces_y, distribution, kindofforce="WFM")
    # Get the size for the loops
    lat_vel, rows, cols = size(distribution)

    # Intermediate array to store the result
    sumresult = zeros(rows, cols)
    sumvelx = zeros(rows, cols)
    sumvely = zeros(rows, cols)
    
    # apparently Julia is fast when it comes to loops
    for node_x = 1:rows
        for node_y = 1:cols
            for q = 1:lat_vel
                sumresult[node_x, node_y] += distribution[q, node_x, node_y]                # Sum up for the height, zeroth moment
                sumvelx[node_x, node_y] += Cd2q9[1,q] * distribution[q, node_x, node_y]     # First moment 
                sumvely[node_x, node_y] += Cd2q9[2,q] * distribution[q, node_x, node_y]     # Second moment
            end
        end
    end
    # Compute the moments based on the forcing approach in collision kernel.
    if kindofforce == "WFM"
        # Normalize the values according to standard LBM
        height = sumresult
        velocity_x = sumvelx./sumresult
        velocity_y = sumvely./sumresult

    elseif kindofforce == "Buick" || kindofforce == "Guo"
        height = sumresult
        velocity_x = (sumvelx .+ (0.5 * forces_x))./sumresult
        velocity_y = (sumvely .+ (0.5 * forces_y))./sumresult
    
    else
        error("no active forcing scheme has been chosen, see collide.jl options are `WFM`, `Buick` and `Guo`")
    end
    
    # Return the computed values for the height velocity vector
    return height, velocity_x, velocity_y
end

"""
    momentsguod2q9(height, velocity_x, velocity_y, forces_x, forces_y, distribution)

Computes the macroscopic quantities `height` and `velocity` from a `distribution function` and the applied `forces`.
See also: [`momentsd2q9`](@ref), [`momentsd1q3`](@ref).

# Arguments
- `height::Array{Float64,2}` : `height`-value at every lattice node, similar to density in standard lattice Boltzmann
- `velocity_x::Array{Float64,2}` : The `x`-component of the macroscopic velocity
- `velocity_y::Array{Float64,2}` : The `y`-component of the macroscopic velocity
- `forces_x::Array{Float64,2}` : The `x`-component of summed up forces
- `forces_y::Array{Float64,2}` : The `y`-component of summed up forces
- `distribution::Array{Float64,3}` : Internal distribution function used for the lattice Boltzmann algorithm

# Examples
```jldoctest
julia> A = [1.0 1.0; 1.0 1.0]; 
julia> B = [0.0 0.0; 0.0 0.0];                  # example without velocity
julia> Fx = [1.0 1.0; 1.0 1.0];                 # And way to much forcing
julia> Fy = [0.0 0.0; 0.0 0.0];
julia> D = fill(1.0, (9, 2, 2))
julia> h_out,velx_out,vely_out = momentsd2q9(A,B,B,Fx,Fy,D); 
julia> println(h_out)  
[9.0 9.0; 9.0 9.0]
```

# Publication


"""
function momentsguod2q9(height, velocity_x, velocity_y, forces_x, forces_y, distribution)
    # Get the size for the loops
    lat_vel, rows, cols = size(distribution)

    # Intermediate array to store the result
    sumresult = zeros(rows, cols)
    sumvelx = zeros(rows, cols)
    sumvely = zeros(rows, cols)
    # apparently Julia is fast when it comes to loops
    for node_x = 1:rows
        for node_y = 1:cols
            for q = 1:lat_vel
                sumresult[node_x, node_y] += distribution[q, node_x, node_y]                # Sum up for the height, zeroth moment
                sumvelx[node_x, node_y] += Cd2q9[1,q] * distribution[q, node_x, node_y]     # First moment 
                sumvely[node_x, node_y] += Cd2q9[2,q] * distribution[q, node_x, node_y]     # Second moment
            end
        end
    end
    # Normalize the values according to standard LBM steps
    height = sumresult
    velocity_x = (sumvelx .+ (0.5 * forces_x))./sumresult
    velocity_y = (sumvely .+ (0.5 * forces_y))./sumresult
    # Return the computed values for the height velocity vector
    return height, velocity_x, velocity_y
end


"""
    momentsd1q3(height, velocity, distribution)

Computes the macroscopic quantities `height` and `velocity` from a `distribution function`.
See also: [`momentsd2q9`](@ref), [`momentsguod2q9`](@ref).

# Arguments
- `height::Array{Float64,1}` : `height`-value at every lattice node, similar to density in standard lattice Boltzmann
- `velocity::Array{Float64,1}` : The macroscopic `velocity` at every lattice node
- `distribution::Array{Float64,2}` : Internal distribution function used for the lattice Boltzmann algorithm

# Examples
```jldoctest
julia> A = [1.0 1.0 1.0 1.0]; 
julia> B = [0.0 0.0 0.0 0.0];                  # example without velocity
julia> D = fill(1.0, (3, 4))
julia> h_out,vel_out = momentsd2q9(A,B,D); 
julia> println(h_out)
[3.0 3.0 3.0 3.0]
```

# Publication
- `Study of the 1D Lattice Boltzmann Shallow Water Equation and Its Coupling to Build a Canal Network` : Guo-Chopard forcing, doi:10.1016/j.jcp.2010.06.022
"""
function momentsd1q3(height, velocity, distribution)
    # Get the size for the loops
    lat_vel, rows = size(distribution)

    # Intermediate array to store the result
    sumresult = zeros(rows)
    sumvel = zeros(rows)

    # apparently Julia is fast when it comes to loops
    for node_x = 1:rows
        for q = 1:lat_vel
            sumresult[node_x] += distribution[q, node_x]                # Sum up for the height, zeroth moment
            sumvel[node_x] += Cd1q3[q] * distribution[q, node_x]        # First moment 
        end
    end
    # Normalize the values according to standard LBM steps
    height = sumresult
    velocity = sumvel./sumresult
    
    # Return the computed values for height and velocity
    return height, velocity
end