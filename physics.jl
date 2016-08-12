# =============================================================================
# Types
# =============================================================================

"""
Photon type.  Contains position, direction, current cell, and wavelength
in nanometers.
"""
type Photon
    xy::Array{Number,2} # Position
    uv::Array{Number,2} # Normalized direction
    last::Integer       # Last cell the photon was inside
    wl::Number          # Wavelength in nanometers
end

"""
Material type.  Contains Sellmeier coefficients and a scaling factor for 
the luminant intensity.
"""
type Material
    # Sellmeier coefficients for index of refraction
    b::Array{Number,1}
    c::Array{Number,1}
    
    scale::Number # Scaling factor for light for artistic reasons
end

# =============================================================================
# Functions
# =============================================================================

"""
Calculate the index of refraction for wavelength ``wl`` using the Sellmeier
equations.
"""
function calculate_n(mat, wl)
    # Using the sellmeier equations, calculate index of refraction
    n_sq = 1 + sum(mat.b .* wl^2 ./ (wl^2 - mat.c))

    return sqrt(n_sq)
end

"""
Performs a surface interaction using the Fresnel equations.
"""
function fresnel(n1, n2, photon, normal)
    # Fresnel formulation of reflection.

    # Find the angle between the two line segments.
    v1 = reshape(photon.uv,2)
    v1 = -v1 / norm(v1)
    
    v2 = reshape(normal,2)
    v2 = v2 / norm(v2) # should be normalized anyways
    
    # Calculate dot product 
    dot_p = dot(v1, v2)
    # For some reason, julia requires 3-long vectors for the builtin cross product.
    cross_p = v1[1] * v2[2] - v1[2] * v2[1]
    if cross_p < 0
        in_angle = acos(dot_p)
    else
        in_angle = -acos(dot_p)
    end

    # If in_angle is outside of [-pi/2, pi/2] we need to recalculate with a negative normal.
    # Concequence of cell intersecting itself.  Can't be caught ahead of time.
    if in_angle >= pi/2 || in_angle <= -pi/2
        v2 = -v2

        # Calculate dot product 
        dot_p = dot(v1, v2)
        # For some reason, julia requires 3-long vectors for the builtin cross product.
        cross_p = v1[1] * v2[2] - v1[2] * v2[1]
        if cross_p < 0
            in_angle = acos(dot_p)
        else
            in_angle = -acos(dot_p)
        end
    end

    # We then calculate what happens.  
    # A) sin_angle > 1 or < -1  - Total internal reflection
    # B) Fresnel equations indicate reflection (randomly sampled)
    # C) Transmit
    
    sin_angle = n1/n2*sin(in_angle)

    reflect = false

    tol = 1.0e-8

    if abs(n1 - n2) < tol
        # Short circuit continuous materials
        reflect = false
        out_angle = in_angle
    elseif sin_angle >= 1.0 || sin_angle <= -1.0
        # Total internal reflectance
        reflect = true
    else
        out_angle = asin(n1 / n2 * sin(in_angle))

        # Assuming non-magnetic material with random polarization
        # We could probably keep track of polarization but it would be a pain
        # to track the orientation.
        p_ref = (((n1 * cos(in_angle) - n2 * cos(out_angle))/(n1 * cos(in_angle) + n2 * cos(out_angle)))^2 +
                ((n1 * cos(out_angle) - n2 * cos(in_angle))/(n1 * cos(out_angle) + n2 * cos(in_angle)))^2) / 2

        if rand() < p_ref
            reflect = true
        end
    end

    if reflect
        out_angle = -in_angle

        # Calculate new unit vector uv with in_angle = out_angle from wall
        # Normal towards inside
        u_new = cos(out_angle) * v2[1] - sin(out_angle) * v2[2]
        v_new = sin(out_angle) * v2[1] + cos(out_angle) * v2[2]
    else
        # Calculate new unit vector uv with angle = out_angle from wall
        # Normal towards outside
        u_new = -cos(out_angle) * v2[1] + sin(out_angle) * v2[2]
        v_new = -sin(out_angle) * v2[1] - cos(out_angle) * v2[2]
    end

    uv = [u_new v_new]
    
    uv = uv / norm(uv)
    
    uv,  reflect
end