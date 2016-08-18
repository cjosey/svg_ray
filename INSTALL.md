Windows Installation
====================

You will need two programs, Python and Julia, as well as a few support libraries.  This document will guide you on setting these up for Windows.

Python
------
Acquire Python from https://www.continuum.io/downloads .  You want the Python 3.5 Windows 64-bit version.  Anaconda saves a bit of time, as it already has matplotlib, pillow, and SciPy installed.

Run the installer.  None of the settings matter except for "Add Anaconda to my PATH environment variable".  Ensure this is checked.

Julia
-----
Download Julia from here: http://julialang.org/downloads/ . You want the 64-bit Windows Self-Extracting Archive.

Install it.  None of the options matter, but it is very important to copy down where Julia is installed for the next step.

Now, you will want to add Julia to your path.  Search through the start menu for "edit the system environment variables".  Click on the "Environment Variables..." button on the bottom.

Click "New...".  For "Variable Name" put ``JULIA_HOME``.  For "Variable Value" put the folder Julia was installed in, with "\bin" appended on the end.  So, for example, if your install directory was:
```
C:\Users\<user name>\AppData\Local\Julia-0.4.6
```
You want to put
```
C:\Users\<user name>\AppData\Local\Julia-0.4.6\bin
```

After that, edit the "Path" variable and (on Windows 10) add a new entry with ``%JULIA_HOME%`` in it.  On earlier versions, append ``;%JULIA_HOME%`` to the end.  Do not delete any of the text that was there orginally.

Open a new terminal, either ``cmd`` or preferrably PowerShell if it is installed.  Run ``julia``.  At the prompt, run the following:

```
Pkg.add("PyPlot")
Pkg.add("HDF5")
Pkg.add("LightXML")
```
It may get stuck on ``LightXML`` at ``Building WinRPM``.  If it gets stuck here too long (a few minutes) pressing Ctrl-C seemed to get it going again.

svg_ray
-------
Download the source code from github.

Open a terminal, and change directory to where the source was downloaded.

Run 
```
julia example.jl
```
Or if you have a large amount of RAM:
```
julia -p <number of CPU cores> example.jl
```

If it runs without errors, run
```
julia plot.jl
```
If that runs without errors, a new ``example.png`` should be created with the result of the simulation.

Now examine the ``example.jl`` file and figure out how it works.  It is heavily documented.

Theory
======
The basic idea behind making an image is that you start your photons somewhere, they hit a glass surface and go in another direction.  As the light flies along, it draws a line on an image.  After enough lines are drawn, the lines merge into a smooth distribution, making an aesthetically pleasing rendering.

There are two main inputs from the user.  The first is the ``.svg`` file.  The second is the particle source.  In the example file, ``example.jl``, this is the ``pgen`` function (in ``example_gen.jl``).

The particle source is a function that creates a photon object every time it is called.  A photon needs 4 data points:
 * A position vector, [x y].
 * A direction vector, [u v].
 * A cell ID (leave as -1 and the code will figure it out)
 * A wavelength.

It is recommended that all but the ID are randomly sampled in some way, as randomness is somewhat aesthetically pleasing.

Let's look at the example generator.  It first generates a wavelength using the ``reject_wl`` function.  What this function does is, given a temperature of light, randomly samples a wavelength using the blackbody spectrum at that temperature via rejection sampling of the blackbody PDF.  Since the rest of the code is written assuming 6500K as "white", anything above 6500 will yield a blue hue, and anything below will yield a red hue.

Then, it generates a position.  In this case, the position is at x = -1000 pixels, and y is anywhere between 600 and 1000 pixels, randomly distributed.  

Finally, since the light source is outside of the geometry, it saves time to make sure that only light pointing in the direction of the "world" is created.  This is what ``sample_valid_angle`` does.  It picks a random angle between 0 and 2pi, calculates the direction vector, and then checks if such a photon will ever intersect the geometry.  If it does, it returns the vector.  If it doesn't, it samples another angle and repeats.
