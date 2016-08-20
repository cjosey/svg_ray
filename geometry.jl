include("cell.jl")
include("physics.jl")

# =============================================================================
# Types
# =============================================================================

"""
A geometry type.  Includes an array of cells ``cell_array`` and a count.
"""
type Geometry
    cell_array # Array of cells
    count      # Number of cells
end

# =============================================================================
# Functions
# =============================================================================

"""
Constructor for Geometry
"""
function Geometry(cell_array)
    count = length(cell_array)
    Geometry(cell_array, count)
end

"""
Calculates the next collision for the photon as it traverses the geometry.
Returns the distance to collision, the point of collision (reduces error),
the normal vector of the surface, and the cell id it collided with
"""
function calculate_collision(geo, photon)
    # The geometry is set up so that cell [i+1] overlaps [i]
    # As such, if a photon is in cell [i], it can only collide with [j], j>= i
    # If it leaves a cell, it can only enter cell [j], j < i, or the void.

    # Calculate nearest collision
    length = typemax(Float64)
    point = []
    normal = []
    cell = 0

    if photon.last == 0
        start_at = 1
    else
        start_at = photon.last
    end

    for i = start_at:geo.count
        inside, cell_point, cell_length, cell_normal = contains(geo.cell_array[i], photon.xy, photon.xy + photon.uv)

        if cell_length < length
            length = cell_length
            point = cell_point
            normal = cell_normal
            cell = i
        end
    end

    length, point, normal, cell
end

"""
Transports a photon through the geometry one step.  If it leaves the geometry,
it is transported a distance equivalent to maxdim.  Returns moved photon.
"""
function transport(geo, photon_last, maxdim)
    tol = 1.0e-6 # A bit large, but there may be a lot of rounding errors here.

    # I don't know if I need to do this.
    photon = deepcopy(photon_last)

    if photon.last == -1
        # Uninitialized
        photon.last = locate(geo, photon)
    end

    # Get collision point
    length, point, normal, cell_ind = calculate_collision(geo, photon)

    # If not collided, move anyways maxdim along uv and return.
    if cell_ind == 0
        photon.xy = photon.xy + maxdim * photon.uv
        return photon
    end

    # Move photon slightly into the next cell
    photon.xy = point - tol * photon.uv

    # Locate what next cell means
    if cell_ind == photon.last
        # Moving out of cell
        new_i = locate(geo, photon, photon.last - 1)
    else
        new_i = cell_ind
    end

    # Get index of refractions
    if photon.last == 0
        n_old = 1
    else
        n_old = calculate_n(geo.cell_array[photon.last].material, photon.wl)
    end

    if new_i == 0
        n_new = 1
    else
        n_new = calculate_n(geo.cell_array[new_i].material, photon.wl)
    end

    # Calculate new angle
    photon.uv, reflect = fresnel(n_old, n_new, photon, normal)

    # Move particle
    photon.xy = point + tol * photon.uv

    # Locate particle
    if !reflect
        photon.last = new_i
    end

    photon
end

"""
Calculates the index of the cell the photon is currently in.
If ``start_at`` is provided, we assume that the index of the new cell is at
least as high as ``start_at``, for performance reasons.
"""
function locate(geo, photon, start_at = -1)
    # Convert start_at to something useful
    if start_at == -1
        start_at = geo.count
    end

    # Default return "void"
    cell = 0

    for i = start_at:-1:1
        inside, new_point, min_length, normal = contains(geo.cell_array[i], photon.xy, photon.xy + photon.uv)
        if inside == 1
            cell = i
            break
        end
    end

    cell
end

"""
Plots all cells in the geometry and saves to filename
"""
function plot_geometry(geo::Geometry, xmax, ymax, filename)
    fig = figure()
    ax = fig[:add_subplot](111)
    for cell in geo.cell_array
        plot_cell(cell, ax)
    end
    ax[:set_aspect]("equal")
    xlim((0,xmax))
    ylim((0,ymax))
    ax[:invert_yaxis]()
    savefig(filename)
end
