"""
Julia script for various collision operators
The actual heart of the lattice Boltzmann method
"""

"""
    collideSRTd2q9(inputdistribution, equilibriumdistribution, forces_x, forces_y)

Computes an output `distribution function` from an input & equilibrium `distribution` and the `forces` acting.   
See also: [`collideSRTguod2q9`](@ref), [`collideTRTguod2q9`](@ref), [`collideSRTd1q3`](@ref).

# Arguments
- `inputdistribution::Array{Float64,3}` : `Input distribution` at time t, the first input is usually the equilibrium. 
- `equilibriumdistribution::Array{Float64,3}` : The `equilibrium distribution` see [`equiliriumd2q9`](@ref). 
- `forces_x::Array{Float64,2}` : x-component of all `forces` applied, contains only rows × cols elements. 
- `forces_y::Array{Float64,2}` : y-component of all `forces` applied.

# Examples
```jldoctest
julia> A = zeros(9,2,2);                    # Does not matter in case of τ=1
julia> B = ones(9,2,2);                     # Output is just the equilibrium
julia> C = zeros(2,2)                       # Forces x-component
julia> D = zeros(2,2)                       # Forces y-component
julia> f_out = collideSRTd2q9(A,B,C,D); 
julia> f_out == B
true
```

"""
function collideSRTd2q9(inputdistribution, equilibriumdistribution, forces_x, forces_y)
    # Get the size of the distributrion functions
    lat_vel, rows, cols = size(equilibriumdistribution)
    outputdistribution = Array{Float64,3}(undef, lat_vel, rows, cols)
    # Loop over all lattice velocities and grid points, ω = 1/τ
    for q = 1:lat_vel
        for node_x = 1:rows
            for node_y = 1:cols
                outputdistribution[q, node_x, node_y] = (1 - ω) * inputdistribution[q, node_x, node_y] .+                   # Contribution from previous time step
                ω * equilibriumdistribution[q, node_x, node_y] .+                                                           # Contribution from equilibrium
                Forceweight[q] .* (Cd2q9[1, q] .* forces_x[node_x, node_y] .+ Cd2q9[2, q] .* forces_y[node_x, node_y])      # Contribution from forces
            end
        end
    end
    # Return the newly calculated distribution function
    return outputdistribution
end
