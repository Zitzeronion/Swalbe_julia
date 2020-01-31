#= 
Julia script to compute the D2Q9 equilibria of the shallow water
lattice Boltzmann model. Several papers can be found that show 
how to implement them here is one of them doi:10.1357/002224099764805174
=#

function equilirium_d2q9(height, velocity)
    """
    Calculating the equilibrium distribution function
    based on the height field and the velocity field.
    """
    # D2Q9 lattice weights
    weights = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36]
    # D2Q9 lattice velocities
    c = [0 1 0 -1 0 1 -1 -1 1; 0 0 1 0 -1 1 1 -1 -1]

    println(c[1, :])
end
