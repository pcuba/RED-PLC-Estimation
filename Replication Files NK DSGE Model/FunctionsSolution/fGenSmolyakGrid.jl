# Calls Smolyak_Elem_Isotrop.jl and Smolyak_Grid.jl to produce a Smolyak GRID

function fGenSmolyakGrid(d,mu,sd::Float64=2)

    #Generate Smolyak element

    Smol_elem=fSmolyak_Elem_Isotrop(d,mu)

    Grid=fSmolyak_Grid(d,mu,Smol_elem)

    return Grid

end
