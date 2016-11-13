push!(LOAD_PATH, ".")
using SvgRay

# This example script will load up all .h5 files in the directory and plot them
# Here are the values worth adjusting:
scale = 1.5
filename = "example.png" # Adjust file type using the filename
# The resulting image may need rotating and or flipping

t = SvgRay.load_tally_hdf5(".")

# Convert from XYZ to RGB
t = SvgRay.tally_to_rgb(t)

# Save as filename
SvgRay.save_tally_rgb(t, filename, scale)
