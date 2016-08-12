using HDF5

# This is an example script to demonstrate the capabilities of this program.
# It is recommended to run this script in parallel for performance reasons.

# Important configuration
# Total number of particles to run
total = 10000
# Resolution of the output file
# Test with small values first.  In parallel, this takes a ton of RAM
# With 4 threads, the below values will use around 6.5 GB of RAM.
x_res = 1600 * 16/9
y_res = 1600
# If you do not want to tweak the source and geometry, adjust the scaling
# value below.  Lower values use less RAM, but the resulting image is lower
# resolution
res_scaling = 1.0
# Ellipses aren't treated exactly.  The below number tells the code how 
# many linesegments to discretize an ellipse into. 20 usually works, but
# check the geometry file.
res = 20
# Filename for the ".h5" results
result_file = "result"
# Filename of the .svg to render
# Here is an example one with a lens and a prism in a meaningless config.
filename = "bg4.svg"
# Filename of geometry plot
# This will show you how the geometry will appear in the simulation.
# useful for debugging purposes
geo_filename = "geometry.pdf"
# It will also show you exactly how this program mangles the coordinate
# system, and if there are any spurious elements in the geometry.


# Calculate the number of threads and particles per thread
np = nworkers()
frac = round(Int, total/np)

# Parallel loop
test = @parallel (+) for i = 1:np
    # All direct includes must be inside a parallel loop
    # using statements can be outside.

    # Load the transport operator
    include("transport.jl")
    # Load in the example source distribution
    include("example_gen.jl")
    # Load in the svg capabilities
    include("example_gen.jl")
    # Load up definition of materials
    include("physics.jl")

    # Here are some example materials
    # Source
    # https://en.wikipedia.org/wiki/Sellmeier_equation
    # and 
    # http://refractiveindex.info/?shelf=glass&book=LASF9&page=SCHOTT
    # Note the * 1.0e6.  It converts um^2 to nm^2 required by the simulation.
    m1 = Material([1.03961212, 0.231792344, 1.01046945], [6.00069867e-3, 2.00179144e-2, 1.03560653e2] * 1.0e6, 1) # BK7 1.521 at 500 nm
    m2 = Material([1.5851495, 0.143559385, 1.08521269], [0.00926681282, 0.0424489805, 105.613573] * 1.0e6, 1) # BAF10 1.678 at 500 nm
    m3 = Material([1.43134930, 0.65054713, 5.3414021], [5.2799261e-3, 1.42382647e-2, 3.25017834e2] * 1.0e6, 1) # Sapphire (ordinary wave) 1.774 at 500 nm
    m4 = Material([1.5039759, 0.55069141, 6.5927379], [5.48041129e-3, 1.47994281e-2, 4.0289514e2] * 1.0e6, 1) # Sapphire (extraordinary wave) 1.766 at 500 nm
    m5 = Material([0.696166300, 0.407942600, 0.897479400], [4.67914826e-3, 1.35120631e-2, 97.9340025] * 1.0e6, 1) # Fused Silica 1.462 at 500 nm
    m6 = Material([2.00029547, 0.298926886, 1.80691843], [0.0121426017, 0.0538736236, 156.530829] * 1.0e6, 1) # LASF9 avg 1.8656 at 500 nm

    # There are two ways you can define materials.  The first is to give all
    # cells the same material.  
    materials = m1
    # or you can define materials for each cell.  A cell corresponds to a path
    # in the SVG file, and they're in the same order as defined in the file.
    # mat = [m2, m5, m3, m6, m3, m1, m6]
    # for example for a geometry with 7 paths.

    # You may also want to adjust the brightness.  By default, the void has a 
    # brightness of 0.25, and the above materials have a brightness of 1. 
    # For "exact" results, if you have the LAB coordinate of a color, 
    # the scale value should be sqrt(L / 100).
    materials.scale = 0.7

    # Our geometry is too small for the resolution of interest.  It is only
    # 608.08 pixels tall.  We can apply a transformation matrix during load.
    # This one scales everything so that the picture is 1600 pixels tall.
    # This is optional.  If you don't want it, don't pass it to 
    # load_geometry
    matr = [[2.631 0 0];[0 2.631 0];[0 0 1]]

    # Load geometry
    geo = load_geometry(filename, geo_filename, materials, res, matr)

    # Perform particle transport
    t = transport(geo, pgen, frac, x_res, y_res, res_scaling)

    # Save the data to hdf5
    h5open(result_file * string(i) * ".h5", "w") do file
        write(file, "data", t.data)
        write(file, "x", t.x)
        write(file, "y", t.y)
        write(file, "nw", t.nw)
        write(file, "scaling_factor", t.scaling_factor)
        write(file, "n_low", t.n_low)
        write(file, "n_high", t.n_high)
    end
    1
end