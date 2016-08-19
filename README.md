# svg_ray

This program performs linear 2D visible light transport through assorted glass objects as defined by an SVG file.  Each path in the SVG file forms a cell, which is then filled with some optically transparent material (glass, sapphire, etc.). Then, photons are traced throughout the geometry to illuminate components.  The result can then be written to an image.  For example:

![example ray traced image](https://raw.githubusercontent.com/cjosey/svg_ray/master/example.jpg)

This program was inspired by this (much faster, but less flexible) program: https://github.com/tunabrain/tantalum

Dependencies
------------
This program has a fairly long list of dependencies:
* Julia v0.4.x (does not yet work in v0.5-rc1)
  * PyPlot.jl
  * PyCall.jl
  * HDF5.jl (for saving the simulation data)
* Python
   * NumPy
   * Matplotlib
   * Scipy
   * Pillow (for ``scipy.imsave``)

Usage
-----
An example simulation is presented in ``example.jl``, which can be run via ``julia example.jl`` or ``julia -p <number of threads> example.jl`` for a parallel simulation.  It will generate a number of result files (one per thread) as well as a plot of how the code sees the geometry.  The results can be merged and converted to an image using ``julia plot.jl``.  The file itself has extensive comments on how it all works.

Limitations
-----------
There are quite a number of SVG features not yet supported.  Below is an incomplete list:
 * Clipping paths
 * Shapes
 * Gradients
 * Exact representations of elliptical paths (approximated by line segments for now)
 
It is hoped that with time some of these can be remedied, but they are not particularly high priority.  Additionally, the SVG reader has only been tested with a narrow range of SVG files.  There are likely many bugs in the SVG reader, so please report if a file does not work and it is not due to the above.

Another limitation is that color does not affect the transport operator.  Experiments were done to allow attenuation by color.  They looked terrible.  For now, the geometry is mathematically equivalent to a 2D clear glass model in which an infinitesmal quantity of light is reflected towards z.  This light then shines through a user defined filter (which can be a function of wavelength).

Polarization is not tracked.   While the Fresnel equations are used in the simulation (and thus reflective probability could be a function of s- and p-polarization), I was not clear how polarization should be tracked, especially as the trajectory changes.

Finally, in the simulation, photons are treated as particles. They do not interfere with each other or themselves, and so certain quantum effects cannot be demonstrated.  Converting the program to support wave interference would require a completely different solver model, and as such will probably never be added.
