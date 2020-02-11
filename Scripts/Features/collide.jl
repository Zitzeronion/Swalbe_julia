"""
Julia script for various collision operators
The actual heart of the lattice Boltzmann method
"""
# include("focringshemes.jl")
"""
    collideSRTd2q9(inputdistribution, equilibriumdistribution, velocity_x, velocity_y, forces_x, forces_y, kindofforcing = "WFM")

Computes an output `distribution function` from an input & equilibrium `distribution` and the `forces` acting.   
See also: [`collideTRTguod2q9`](@ref), [`collideSRTd1q3`](@ref).

# Arguments
- `inputdistribution::Array{Float64,3}` : `Input distribution` at time t, the first input is usually the equilibrium. 
- `equilibriumdistribution::Array{Float64,3}` : The `equilibrium distribution` see [`equiliriumd2q9`](@ref). 
- `velocity_x::Array{Float64,2}` : x-component of `velocity` calculated with [`momentsd2q9`](@ref) or [`momentsguod2q9`](@ref). 
- `velocity_y::Array{Float64,2}` : y-component of `velocity`. 
- `forces_x::Array{Float64,2}` : x-component of all `forces` applied, contains only rows × cols elements. 
- `forces_y::Array{Float64,2}` : y-component of all `forces` applied.
- `kindofforcing::String` : Determines which kind of forcing model is used, options are "WFM", "Buick" and "Guo"

# Examples
```jldoctest
julia> distfromprevioust = zeros(9,2,2);                    # Does not matter in case of τ=1
julia> distequilibrium = ones(9,2,2);                       # Output is just the equilibrium
julia> velo_x = zeros(2,2)                                  # Velocity x-component
julia> velo_y = zeros(2,2)                                  # Velocity y-component
julia> forc_x = zeros(2,2)                                  # Forces x-component
julia> forc_y = zeros(2,2)                                  # Forces y-component
julia> use_force = "WFM"
julia> f_out = collideSRTd2q9(distfromprevioust, distequilibrium, velo_x, velo_y, forc_x, forc_y, use_force); 
julia> f_out == B
true
```

# Publication
- `An evaluation of force terms in the lattice Boltzmann models in simulating shallow water flows over complex topography` : Weighting factor mode, doi:10.1002/fld.4726
- `Gravity in a lattice Boltzmann model` : Buick, J. M. and Greated, C. A., doi:10.1103/PhysRevE.61.5307
- `Discrete lattice effects on the forcing term in the lattice Boltzmann method` : Guo, Zheng, Shi, doi:10.1103/PhysRevE.65.046308

"""
function collideSRTd2q9(inputdistribution, equilibriumdistribution, velocity_x, velocity_y, forces_x, forces_y, kindofforcing = "WFM")
    # Get the size of the distributrion functions
    lat_vel, rows, cols = size(equilibriumdistribution)
    outputdistribution = Array{Float64,3}(undef, lat_vel, rows, cols)
    
    # Determine which kind of forcing should be used, default is weighting factor method, F = 3*w/(c²)*cᵢ.Fᵢ
    if kindofforcing == "WFM"
        # Loop over all lattice velocities and grid points, for clarity ω = 1/τ (relaxation rate)
        for q = 1:lat_vel
            for node_x = 1:rows
                for node_y = 1:cols
                    outputdistribution[q, node_x, node_y] = (1 - ω) * inputdistribution[q, node_x, node_y] .+                   # Contribution from previous time step
                    ω * equilibriumdistribution[q, node_x, node_y] .+                                                           # Contribution from equilibrium
                    Forceweight[q] .* (Cd2q9[1, q] .* forces_x[node_x, node_y] .+ Cd2q9[2, q] .* forces_y[node_x, node_y])      # Contribution from forces, WFM 
                end
            end
        end

    # Another kind is the Buick-Greated Scheme, F = (1-1/2τ)*3w/(c^2)cᵢFᵢ
    # Can only be use with [`equiliriumguod2q9`](@ref)
    elseif kindofforcing == "Buick"
        for q = 1:lat_vel
            for node_x = 1:rows
                for node_y = 1:cols
                    outputdistribution[q, node_x, node_y] = (1 - ω) * inputdistribution[q, node_x, node_y] .+                                   # Contribution from previous time step
                    ω * equilibriumdistribution[q, node_x, node_y] .+                                                                           # Contribution from equilibrium
                    (1-1/(2*τ)) * Forceweight[q] .* (Cd2q9[1, q] .* forces_x[node_x, node_y] .+ Cd2q9[2, q] .* forces_y[node_x, node_y])        # Contribution from forces
                end
            end
        end

    # The last for now is the Guo forcing term, F = (1-1/2τ)*3w*[(cᵢ - uᵢ)/(c^2) + 3*(cᵢ.uᵢ)/(c^2)^2 cᵢ]Fᵢ
    # Can only be use with [`equiliriumguod2q9`](@ref)
    elseif kindofforcing == "Guo"
        for q = 1:lat_vel
            for node_x = 1:rows
                for node_y = 1:cols
                    outputdistribution[q, node_x, node_y] = (1 - ω) * inputdistribution[q, node_x, node_y] .+                                   # Contribution from previous time step
                    ω * equilibriumdistribution[q, node_x, node_y] .+                                                                           # Contribution from equilibrium
                    (1-1/(2*τ)) * Forceweight[q] .* ((Cd2q9[1, q] .- velocity_x[node_x,node_y]) .* forces_x[node_x, node_y] .+                  # Term (c_(alpha i) - u_i) F_i
                                                     (Cd2q9[2, q] .- velocity_y[node_x,node_y]) .* forces_y[node_x, node_y] .+ 
                                                     3 * (Cd2q9[1, q] .* velocity_x[node_x,node_y] .+ Cd2q9[2, q] .* velocity_y[node_x,node_y])/C²[q] .*    # Term 3(e.u)/e^2*e.F 
                                                     (Cd2q9[1, q] .* forces_x[node_x, node_y] .+ Cd2q9[2, q] .* forces_y[node_x, node_y])) 
                end
            end
        end
    else
        error("no valid forcing scheme was chosen, see collide.jl for options")
    end

    # Return the newly calculated distribution function
    return outputdistribution
end


