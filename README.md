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
An example simulation is presented in ``example.jl``, which can be run via ``julia example.jl`` or ``julia -p <number of threads> example.jl``.  It will generate a number of result files (one per thread) as well as a plot of how the code sees the geometry.  The results can be merged and converted to an image using ``julia plot.jl``.  The file itself has extensive comments on how it all works.

It is worth mentioning that due to the size of the tallies, the example simulation will take 1.6 GB of RAM per thread.  This can be reduced by adjusting the value of ``res_scaling`` in line 17 of ``example.jl``.

Limitations
-----------
Although Bézier curves are treated exactly with this simulation, elliptic curves are not.  They are discretized into linear segments.  This is a fairly coarse (but adjustable) approximation, so it is recommended to avoid elliptic paths in the SVG file.  Hopefully, this limitation can be removed.

In the simulation, photons are treated as particles. They do not interfere with each other or themselves, and so certain quantum effects cannot be demonstrated.

The SVG reader has not been extensively tested.  If you have a file that cannot be read, bring it up as an issue.  Shapes, however, probably will not be added.  If possible, convert them to paths.

