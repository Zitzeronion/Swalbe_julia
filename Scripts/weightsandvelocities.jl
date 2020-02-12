"""
Some defined lattice Boltzmann constants.
Of course the D2Q9 weights but also the lattice velocity vectors.

"""
# D2Q9:
# D2Q9 lattice weights 
Weightsd2q9 = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36]
# D2Q9 lattice velocity vector: (0,0) (1,0) (0,1) ....
Cd2q9 = [0 1 0 -1 0 1 -1 -1 1; 0 0 1 0 -1 1 1 -1 -1]

# D1Q3:
# D1Q3 lattice weights
Weightsd1q3 = [2/3 1/6 1/6]
# D1Q3 lattice velocity
Cd1q3 = [0 1 -1]

# Calculate force weighting based on doi:10.1002/fld.4726
which = "Fancy"
Forceweight = zeros(9)
C² = zeros(9)
if which == "Simple"
    for q in 2:9
        C²[q] = 1.0
        Forceweight[q] = 3 * Weightsd2q9[q] ./ C²[q]
    end
elseif which == "Fancy"
    for q in 2:9
        C²[q] = Cd2q9[1,q]*Cd2q9[1,q] + Cd2q9[2,q]*Cd2q9[2,q]
        Forceweight[q] = 3 * Weightsd2q9[q] ./ C²[q]
    end
end
# This due to Zhous "centered based forcing", to switch back to the usual lattice Boltzmann, divide by 1 instead of Cd2q9[1,q]*Cd2q9[1,q] + Cd2q9[2,q]*Cd2q9[2,q]

# Speed of sound and speed of sound squared
cₛ = 1/sqrt(3)
cₛ² = 1/3.0

"""
Used in [`momentsd2q9`](@ref), [`equilibriumd2q9`](@ref).
"""