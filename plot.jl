using HDF5
include("tally.jl")

# This example script will load up all .h5 files in the directory and plot them
# Here are the values worth adjusting:
scale = 1.5
filename = "example.png" # Adjust file type using the filename
# The resulting image may need rotating and or flipping

t = []

for fn in readdir()
    if contains(fn, ".h5")
        if t == []
            data = []
            x = []
            y = []
            nw = []
            scaling_factor = []
            n_low = []
            n_high = []
            h5open(fn, "r") do file
                data = read(file, "data")
                x = read(file, "x")
                y = read(file, "y")
                nw = read(file, "nw")
                scaling_factor = read(file, "scaling_factor")
                n_low = read(file, "n_low")
                n_high = read(file, "n_high")
            end
            t = Tally(data, x, y, nw, scaling_factor, n_low, n_high)
        else
            h5open(fn, "r") do file
                t.data += read(file, "data")
            end
        end
    end
end

# Convert from XYZ to RGB
t = tally_to_rgb(t)

# Save as filename
save_tally_rgb(t, filename, scale)