using HDF5

# This is an example script to demonstrate the capabilities of this program.
# It is recommended to run this script in parallel for performance reasons.

# Important configuration
# Total number of particles to run
total = 10000
# Resolution of the output file
x_res = 1600 * 16/9
y_res = 1600
# Ellipses aren't treated exactly.  The below number tells the code how 
# many linesegments to discretize an ellipse into. 20 usually works, but
# check the geometry file.
res = 20
# The program has the ability to read the colors of some files, namely, 
# if they're not gradients and are in the "style" command.  These can 
# be transformed into a spectrographic filter, or a luminant multiplier.
# The below variable controls which
color_mode = 1 # Assume all cells are pure white
# color_mode = 2 # Calculate greyscale information
# color_mode = 3 # Use color filters
# Additionally, the void has a scaling value. This is how bright the void is.
void_scale = 0.5
# Realism is color_mode = 1, void_scale = 1.0
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
# You will have to uncomment the line:
# plot_geometry(geo, geo_filename)
# in the loop below. It's commeented out as it occasionally crashes the program
# due to weird python-parallel julia interactions.


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

    # Our geometry is too small for the resolution of interest.  It is only
    # 608.08 pixels tall.  We can apply a transformation matrix during load.
    # This one scales everything so that the picture is 1600 pixels tall.
    # This is optional.  If you don't want it, don't pass it to 
    # load_geometry
    matr = [[2.631 0 0];[0 2.631 0];[0 0 1]]

    # Load geometry
    geo = load_geometry(filename, materials, res, color_mode, matr)

    # Optionally, plot geometry.
    # Python sometimes crashes if you do this in a parallel loop.
    # Useful for debugging your g
    # plot_geometry(geo, geo_filename)

    # Perform particle transport
    t = transport(geo, pgen, frac, x_res, y_res, void_scale)

    # Save the data to hdf5
    h5open(result_file * string(i) * ".h5", "w") do file
        write(file, "data", t.data)
        write(file, "x", t.x)
        write(file, "y", t.y)
        write(file, "scaling_factor", t.scaling_factor)
    end
    1
end