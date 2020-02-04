"""
Computation of macroscopic quantities.
For the shallow water lattice Boltzmann method the sum of the distribution function is the height.
The sum of distribution function times the lattice velocity is the macroscopic velocity.
Theoretically this quantities are the moments of the distribution function.
"""

# D2Q9 lattice velocity vector: (0,0) (1,0) (0,1) ....
Cd2q9 = [0 1 0 -1 0 1 -1 -1 1; 0 0 1 0 -1 1 1 -1 -1]


"""
    momentsd2q9!(height, velocity_x, velocity_y, distribution)

Computes the macroscopic quantities `height` and `velocity` from a `distribution function`.
See also: [`equiliriumd2q9`](@ref), [`equiliriumd1q3`](@ref).

# Arguments
- `height::Array{Float64,2}` : `height`-value at every lattice node, similar to density in standard lattice Boltzmann
- `velocity_x::Array{Float64,2}` : The `x`-component of the macroscopic velocity
- `velocity_y::Array{Float64,2}` : The `y`-component of the macroscopic velocity
- `distribution::Array{Float64,3}` : Internal distribution function used for the lattice Boltzmann algorithm

# Examples
```jldoctest
julia> A = [1 1; 1 1]; 
julia> B = [0 0; 0 0];                  # example without velocity
julia> D = fill(1.0, (9, 2, 2))
julia> C = momentsd2q9!(A,B,B,D); 
julia> println(C[1, :, :])
[1.0 1.0; 1.0 1.0]
```

"""
function momentsd2q9!(height, velocity_x, velocity_y, distribution)
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
                sumresult[node_x, node_y] += distribution[q, node_x, node_y]                       # Sum up for the height, zeroth moment
                sumvelx[node_x, node_y] += Cd2q9[1,q] .* distribution[q, node_x, node_y]     # First moment 
                sumvely[node_x, node_y] += Cd2q9[2,q] .* distribution[q, node_x, node_y]     # Second moment
            end
        end
    end

    height = sumresult
    velocity_x = sumvelx./sumresult
    velocity_y = sumvely./sumresult
end