#=
Core julia script that links together all the parts.
This script will be used to run lattice Boltzmann simulations.
The lattice Boltzmann methods is based on the shallow water lattice Boltzmann,
first derived by Rick Salmon doi:10.1357/002224099764805174
=#

println("Hello world!")
A = [1 1; 1 1]
include("weightsandvelocities.jl")
include("Features/equilibrium.jl")
include("Features/macroquantities.jl")
C = zeros(2,2)
B = equiliriumnogd2q9(A, C, C)

println(B[1,:,:])