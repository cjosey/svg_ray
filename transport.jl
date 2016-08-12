include("svg_convert.jl")
include("cell.jl")
include("physics.jl")
include("geometry.jl")
include("tally.jl")

"""
Performs photon transport on geo using the photon generator provided.
"""
function transport(geo, photon_gen, n, maxdimx, maxdimy)
    t = Tally(maxdimx, maxdimy)
    maxdim = sqrt(maxdimx^2 + maxdimy^2)*2
    for i = 1:n
        p = pgen()
        p.last = locate(geo, p)
        path = p.xy
        # Tally a path for at most 30 reflections
        # Note, does not count going from cell a to b just cell a to a.
        j = 1
        scaling = []
        while j < 30
            last_cell = p.last
            p = transport(geo, p, maxdim)
            path = vcat(path, p.xy)

            # Append luminance information
            if last_cell > 0
                scale = geo.cell_array[last_cell].material.scale
            else
                scale = 0.25
            end
            push!(scaling, scale)

            if p.xy[1] > maxdimx || p.xy[1] < 0 || p.xy[2] > maxdimy || p.xy[2] < 0
                break
            end
            if last_cell == p.last
                # Reflected
                j += 1
            end

            if size(path)[1] > 1000
                println(size(path), " ",j," ", p.xy, " ", last_cell, " ", p.last, " ", p.uv)
            end
        end
        tally_path(t, path, p.wl, scaling)
        if i % 100 == 0
            println("Sampled photon ", i)
        end
    end
    tally_to_xyz(t)
end