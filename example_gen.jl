using SvgRay

# This file contains an example set of particle generation parameters.  Tune
# to your artistic preference.  It is also worth noting, as this is a linear
# operation (until gamma correction in the plotting), you can always run the
# simulation twice with different sources, and add the results together.  

"""
This function uses rejection sampling to sample a wavelength from a blackbody
radiation source.  This simulation is rendered into D65 color space.
As such, if you use T = 6500, you'll get "white" light.
"""
function reject_wl(T::Real)
    # Samples a wavelength using black body radiation
    c1 = 3.74183e-16 
    c2 = 1.4388e-2

    maxp = c1 * T^5 / c2^5 * 21.2014356605 # Derived with mathematica
    p_lambda = 0
    rand1 = 0
    rand2 = 0
    rand1m = 0
    for i = 1:1000
        rand1 = rand() * 400 + 380
        rand2 = rand()
        rand1m = rand1 * 1.0e-9 # Convert to m
        p_lambda = c1 * rand1m^(-5)/(exp(c2/(rand1m * T)) - 1)
        if p_lambda > rand2 * maxp
            return rand1
        end
    end
    println(p_lambda, " ",rand1, " ",rand2, " ", rand1m)
end

"""
Given an xy coordinate and the size of the geometry, randomly samples an
isotropic angle that is guaranteed to intersect the image.  Improves 
performance for sources outside the geometry.
"""
function sample_valid_angle(xy)
    # Repeatedly samples an angle isotropically until one intersects the
    # geometry on at least one side.

    # Bounding box
    xmin = 0
    ymin = 0
    xmax = 1600 * 16/9
    ymax = 1600
    x1 = [xmin xmin xmax xmax]
    x2 = [xmin xmax xmax xmin]
    y1 = [ymin ymax ymax ymin]
    y2 = [ymax ymax ymin ymin]
    x3 = xy[1]
    y3 = xy[2]
    while true
        angle = rand() * 2 * pi
        u = cos(angle)
        v = sin(angle)

        x4 = xy[1] + u
        y4 = xy[2] + v

        # Calculate all intersections for the full lines
        denom = (x1 - x2) .* (y3 - y4) - (y1 - y2) .* (x3 - x4)

        p_x = ((x1 .* y2 - y1 .* x2).*(x3 - x4) - (x1 - x2).*(x3 .* y4 - y3 .* x4)) ./ denom
        p_y = ((x1 .* y2 - y1 .* x2).*(y3 - y4) - (y1 - y2).*(x3 .* y4 - y3 .* x4)) ./ denom

        int1 = (p_y[1] > ymin && p_y[1] < ymax && ((x3 < xmin && u > 0) || (x3 > xmin && u < 0)))
        int2 = (p_x[2] > xmin && p_x[2] < xmax && ((y3 < ymax && v > 0) || (y3 > ymax && v < 0)))
        int3 = (p_y[3] > ymin && p_y[3] < ymax && ((x3 < xmax && u > 0) || (x3 > xmax && u < 0)))
        int4 = (p_x[4] > xmin && p_x[4] < xmax && ((y3 < ymin && v > 0) || (y3 > ymin && v < 0)))

        if int1 || int2 || int3 || int4
            return [u v]
        end
    end
end

"""
The particle generator itself.  Samples a random particle from a blackbody
source of 6500K starting at x = -1000 and y = [600, 1200].  The angle is 
isotropic.
"""
function pgen()
    # Randomly generates a photon
    wl = reject_wl(6500.0)
    xy = [-1000 600 + 400*rand()]
    uv = sample_valid_angle(xy)

    p = SvgRay.Photon(xy, uv, -1, wl)
end