export cellaverage, cellcenters, coordinates 
export visualize

"""
function cellaverage(Q; M = nothing)

# Description
Compute the cell-average of Q given the mass matrix M.
Assumes that Q and M are the same size

# Arguments
`Q`: MPIStateArrays (array)

# Keyword Arguments
`M`: Mass Matrix (array)

# Return
The cell-average of Q
"""
function cellaverage(Q; M = nothing)
    if M==nothing
        return nothing
    end
    return (sum(M .* Q, dims = 1) ./ sum(M , dims = 1))[:]
end

"""
function coordinates(grid::DiscontinuousSpectralElementGrid)

# Description
Gets the (x,y,z) coordinates corresponding to the grid

# Arguments
- `grid`: DiscontinuousSpectralElementGrid

# Return
- `x, y, z`: views of x, y, z coordinates
"""
function coordinates(grid::DiscontinuousSpectralElementGrid)
    x = view(grid.vgeo, :, grid.x1id, :)   # x-direction	
    y = view(grid.vgeo, :, grid.x2id, :)   # y-direction	
    z = view(grid.vgeo, :, grid.x3id, :)   # z-direction
    return x, y, z
end

"""
function cellcenters(Q; M = nothing)

# Description
Get the cell-centers of every element in the grid

# Arguments
- `grid`: DiscontinuousSpectralElementGrid

# Return
- Tuple of cell-centers
"""
function cellcenters(grid::DiscontinuousSpectralElementGrid)
    x, y, z = coordinates(grid)
    M = view(grid.vgeo, :, grid.Mid, :)  # mass matrix
    xC = cell_average(x, M = M)
    yC = cell_average(y, M = M)
    zC = cell_average(z, M = M)
    return xC[:], yC[:], zC[:]
end

"""
function visualize(g::DiscontinuousSpectralElementGrid)

# Description
- Use Makie to visualize a 3D grid

# Arguments
- `g`: A DiscontinuousSpectralElementGrid
"""
function visualize(g::DiscontinuousSpectralElementGrid)
    x, y, z = coordinates(g)

    e = collect(1:size(x)[2])
    nfaces = 6
    faces = collect(1:nfaces)
    opacities = collect(range(0,1, length = 10))

    scene, layout = layoutscene()

    element = Node{Any}(e[1])
    face = Node{Any}(faces[6])
    opacity = Node{Any}(opacities[1])

    tmpx = @lift(x[:, $element])
    tmpy = @lift(y[:, $element])
    tmpz = @lift(z[:, $element])

    tmpxf = @lift(x[g.vmap⁻[:, $face, $element]])
    tmpyf = @lift(y[g.vmap⁻[:, $face, $element]])
    tmpzf = @lift(z[g.vmap⁻[:, $face, $element]])

    total_color = @lift(RGBA(0,0,0, $opacity))

    lscene = layout[1:4, 2:4] = LScene(scene)
    Makie.scatter!(lscene, x[:], y[:], z[:], color = total_color, markersize = 100.0, strokewidth = 0)
    Makie.scatter!(lscene, tmpx, tmpy, tmpz, color = RGBA(0,0,0,0.5), markersize = 100.0, strokewidth = 0, camera = cam3d!)
    Makie.scatter!(lscene, tmpxf, tmpyf, tmpzf, color = RGBA(1,0,0,1.0), markersize = 100.0, strokewidth = 0, camera = cam3d!)
    supertitle = layout[1,2] = LText(scene, " "^10 * " Gauss-Lobatto Points " * " "^10, textsize = 50, color = :black)

    menu  = LMenu(scene, options = zip(e,e))
    menu2 = LMenu(scene, options = zip(faces,faces))
    menu3 = LMenu(scene, options = zip(opacities,opacities))
    layout[1, 1] = vgrid!(
        LText(scene, "Element", width = nothing),
        menu,
        LText(scene, "Face", width = nothing),
        menu2,
        LText(scene, "Opacity", width = nothing),
        menu3,
    )
    on(menu.selection) do s
        element[] = s
    end
    on(menu2.selection) do s
        face[] = s
    end
    on(menu3.selection) do s
        opacity[] = s
    end
    display(scene)
end