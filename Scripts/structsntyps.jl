struct Latticepoint
    """
    Structure to address a lattice point on the grid.
    Contains a x-component and a y-component.
    """
    x :: Int64
    y :: Int64
end

struct Distributionfunction
    """
    Structure to declare the shape of the distribution functions.
    Generally all of them should have a name, e.g. equilibrium, temp...
    The number of components, for D2Q9 = 9, for D1Q3 = 3.
    Finally it should contain an address to a single lattice point
    """
    name :: Char
    components :: Int
    adr :: Latticepoint
    values :: Float64
end

struct Observable
    """
    Physical quantities which can be defined at ever Latticepoint
    """
    name :: Char
    adr :: Latticepoint
    values :: Float64
end
