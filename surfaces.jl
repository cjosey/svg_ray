# =============================================================================
# Types
# =============================================================================

"""
LineArray contains a series of line segments
"""
type LineArray
    points::Array{Number,2} # Line segments stored as [x1 y1 x2 y2] columns

    bbox::Array{Number,2} # Bounding boxes, stored as [xmin, ymin, xmax, ymax] columns
    # We can't reduce the memory footprint by using a ring, unfortunately,
    # as it is possible these lines join other curve types.
end

"""
QBezier contains an array of quadratic bezier curves.
"""
type QBezier
    points::Array{Number,2} # Stored as [x1; y1; xc1; yc1; x2; y2] rows
end

"""
CBezier contains an array of cubic bezier curves.
"""
type CBezier
    points::Array{Number,2} # Stored as [x1; y1; xc1; yc1; xc2; yc2; x2; y2] rows
end

# Note, although ellipses are probably possible, SVG specification makes them
# convoluted messes, and it is recommended to avoid them.  Regardless, the code
# supports discretization of ellipses into a LineArray to avoid most issues.

# =============================================================================
# Functions
# =============================================================================

function LineArray(points)
    bbox = zeros(size(points))

    for i = 1:size(points)[1]
        # Calculate min, max for bounding box
        x1 = points[i,1]
        y1 = points[i,2]
        x2 = points[i,3]
        y2 = points[i,4]

        bbox[i,1] = min(x1, x2)
        bbox[i,2] = max(x1, x2)
        bbox[i,3] = min(y1, y2)
        bbox[i,4] = max(y1, y2)
    end
    LineArray(points, bbox)
end

"""
Applies transformation matrix mat to the array
"""
function apply_transform!(line::LineArray, mat)
    for i = 1:size(line.points)[1]
        # Apply to lines

        v1 = [line.points[i,1:2] 1]'
        v2 = [line.points[i,3:4] 1]'
        v1t = mat*v1
        v2t = mat*v2
        line.points[i,1:2] = v1t[1:2]
        line.points[i,3:4] = v2t[1:2]

        # Recalculate bounding box
        x1 = v1t[1]
        y1 = v1t[2]
        x2 = v2t[1]
        y2 = v2t[2]

        line.bbox[i,1] = min(x1, x2)
        line.bbox[i,2] = max(x1, x2)
        line.bbox[i,3] = min(y1, y2)
        line.bbox[i,4] = max(y1, y2)
    end
end

"""
Applies transformation matrix mat to the array
"""
function apply_transform!(qb::QBezier, mat)
    for i = 1:size(qb.points)[2]
        # Apply to each point
        # Fun fact, any linear transformation of the points of a
        # Bezier curve carries over to the path.

        v1 = [qb.points[1:2, i]; 1]
        v2 = [qb.points[3:4, i]; 1]
        v3 = [qb.points[5:6, i]; 1]

        v1t = mat*v1
        v2t = mat*v2
        v3t = mat*v3

        qb.points[1:2, i] = v1t[1:2]
        qb.points[3:4, i] = v2t[1:2]
        qb.points[5:6, i] = v3t[1:2]

    end
end

"""
Applies transformation matrix mat to the array
"""
function apply_transform!(cb::CBezier, mat)
    for i = 1:size(cb.points)[2]
        # Apply to each point
        # Fun fact, any linear transformation of the points of a
        # Bezier curve carries over to the path.

        v1 = [cb.points[1:2, i]; 1]
        v2 = [cb.points[3:4, i]; 1]
        v3 = [cb.points[5:6, i]; 1]
        v4 = [cb.points[7:8, i]; 1]

        v1t = mat*v1
        v2t = mat*v2
        v3t = mat*v3
        v4t = mat*v4

        cb.points[1:2, i] = v1t[1:2]
        cb.points[3:4, i] = v2t[1:2]
        cb.points[5:6, i] = v3t[1:2]
        cb.points[7:8, i] = v4t[1:2]
    end
end

"""
Calculates the number of, and nearest, intersections, as well as the normal line.
"""
function intersect(line::LineArray, p1, p2)
    # For some reason, when this was inside cell.jl, the vectored version 
    # worked well.  I have had to un-vector it to get it working now.
    tol = 1.0e-10

    x3 = p1[1]
    y3 = p1[2]

    # Pick any point as long as it is not coincident with first one.
    x4 = p2[1]
    y4 = p2[2]

    l_count = size(line.points)[1]

    new_point = []
    min_length = typemax(Float64)
    normal = []

    # We assume that if it's on a parallel line, then it is outside.  This also
    # should never happen in practice, but we don't care too much, really. 
    # This algorithm shouldn't be used for engineering applications.

    # If any of the points above are also in the bounding box of the segment, 
    # then a collision is recorded.
    intersects = zeros(Int, l_count)
    for i = 1:l_count
        # Calculate intersection point
        x1 = line.points[i,1]
        y1 = line.points[i,2]
        x2 = line.points[i,3]
        y2 = line.points[i,4]

        denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
        p_x = ((x1 * y2 - y1 * x2)*(x3 - x4) - (x1 - x2)*(x3 * y4 - y3 * x4)) / denom
        p_y = ((x1 * y2 - y1 * x2)*(y3 - y4) - (y1 - y2)*(x3 * y4 - y3 * x4)) / denom

        # Check if in bounding box (with some wiggle room)
        iv = (p_x > line.bbox[i,1] - tol) &&
             (p_x < line.bbox[i,2] + tol) &&
             (p_y > line.bbox[i,3] - tol) &&
             (p_y < line.bbox[i,4] + tol)
        # Check if in right direction
        if iv == 1
            if x3 < x4
                gx0 = p_x > x3
            elseif x4 < x3
                gx0 = p_x < x3
            else
                gx0 = true
            end
            iv = iv && gx0
        end
        if iv == 1
            if y3 < y4
                gy0 = (p_y > y3)
            elseif y4 < y3
                gy0 = (p_y < y3)
            else
                gy0 = true
            end
            iv = iv && gy0
        end
        intersects[i] = iv

        # If intersects, save some time, compute distance and normal
        if iv == 1
            path_length = sqrt((x3 - p_x)^2 + (y3 - p_y)^2)
            if path_length < min_length
                min_length = path_length
                new_point = [p_x p_y]

                # Calculate normal
                dx = x2 - x1
                dy = y2 - y1
                normal = [-dy dx]
                normal = normal / norm(normal)
            end
        end
    end

    # Count number of intersections
    int_count = sum(intersects)
    # Return count, point, distance, and normal
    int_count, new_point, min_length, normal
end

"""
Calculates the number of, and nearest, intersections, as well as the normal line.
"""
function intersect(cb::CBezier, p1, p2)
    # I apologize for the complexity of this monstrosity
    tol = 1.0e-10

    x11 = p1[1]
    y11 = p1[2]

    # Pick any point as long as it is not coincident with first one.
    x12 = p2[1]
    y12 = p2[2]

    new_point = []
    min_length = typemax(Float64)
    normal = []

    # I got the idea for this from
    # https://www.particleincell.com/2013/cubic-line-intersection/
    # Though, I didn't use the code.

    int_count = 0
    for i = 1:size(cb.points)[2]
        x01 = cb.points[1,i]
        y01 = cb.points[2,i]
        x02 = cb.points[3,i]
        y02 = cb.points[4,i]
        x03 = cb.points[5,i]
        y03 = cb.points[6,i]
        x04 = cb.points[7,i]
        y04 = cb.points[8,i]

        # Calculate coefficients of characteristic polynomial
        t3c = (x01*y11 - x01*y12 - 3*x02*y11 + 3*x02*y12 + 3*x03*y11 - 3*x03*y12 - x04*y11 + x04*y12 - x11*y01 + 3*x11*y02 - 3*x11*y03 + x11*y04 + x12*y01 - 3*x12*y02 + 3*x12*y03 - x12*y04)
        t2c = (-3*x01*y11 + 3*x01*y12 + 6*x02*y11 - 6*x02*y12 - 3*x03*y11 + 3*x03*y12 + 3*x11*y01 - 6*x11*y02 + 3*x11*y03 - 3*x12*y01 + 6*x12*y02 - 3*x12*y03)
        t1c = (3*x01*y11 - 3*x01*y12 - 3*x02*y11 + 3*x02*y12 - 3*x11*y01 + 3*x11*y02 + 3*x12*y01 - 3*x12*y02)
        t0c = -x01*y11 + x01*y12 + x11*y01 - x11*y12 - x12*y01 + x12*y11

        a2 = complex(t2c/t3c)
        a1 = complex(t1c/t3c)
        a0 = complex(t0c/t3c)

        # https://en.wikipedia.org/wiki/Cubic_function
        delta0 = a2^2 - 3*a1
        delta1 = 2*a2^3 - 9*a2*a1 + 27*a0

        C = ((delta1 + sqrt(delta1^2 - 4*delta0^3))/2)^(1/3)

        u1C = C
        u2C = C*(-1 + 1im*sqrt(3))/2
        u3C = C*(-1 - 1im*sqrt(3))/2

        z1 = -1/3*(a2 + u1C + delta0/u1C)
        z2 = -1/3*(a2 + u2C + delta0/u2C)
        z3 = -1/3*(a2 + u3C + delta0/u3C)

        roots = [z1, z2, z3]

        for root in roots
            # Check if their imaginary component is sufficiently small
            if abs(imag(root)) < tol
                t = real(root)
                # Check bounding box on Bezier
                if t <= 1.0 + tol && t >= -tol
                    # Get px, py
                    p_x = (1 - t)^3 * x01 + 3*(1 - t)^2 * t * x02 + 3*(1 - t)*t^2 * x03 + t^3 * x04
                    p_y = (1 - t)^3 * y01 + 3*(1 - t)^2 * t * y02 + 3*(1 - t)*t^2 * y03 + t^3 * y04

                    # Check that collision is forward from p1
                    if x11 < x12
                        gx0 = p_x > x11
                    elseif x12 < x11
                        gx0 = p_x < x11
                    else
                        gx0 = true
                    end
                    gy0 = false
                    if gx0
                        if y11 < y11
                            gy0 = (p_y > y11)
                        elseif y11 < y11
                            gy0 = (p_y < y11)
                        else
                            gy0 = true
                        end
                    end

                    if gy0
                        # Collided, in both bounding boxes
                        int_count += 1

                        path_length = sqrt((x11 - p_x)^2 + (y11 - p_y)^2)
                        if path_length < min_length
                            min_length = path_length
                            new_point = [p_x p_y]

                            # Calculate normal
                            dx = 3*(1-t)^2 * (x02 - x01) + 6*(1-t)*t*(x03 - x02) + 3*t^2*(x04 - x03)
                            dy = 3*(1-t)^2 * (y02 - y01) + 6*(1-t)*t*(y03 - y02) + 3*t^2*(y04 - y03)
                            normal = [-dy dx]
                            normal = normal / norm(normal)
                        end
                    end
                end
            end
        end
    end

    int_count, new_point, min_length, normal
end

"""
Calculates the number of, and nearest, intersections, as well as the normal line.
"""
function intersect(qb::QBezier, p1, p2)
    # I apologize for the complexity of this monstrosity
    tol = 1.0e-10

    x11 = p1[1]
    y11 = p1[2]

    # Pick any point as long as it is not coincident with first one.
    x12 = p2[1]
    y12 = p2[2]

    new_point = []
    min_length = typemax(Float64)
    normal = []

    # I got the idea for this from
    # https://www.particleincell.com/2013/cubic-line-intersection/
    # Though, I didn't use the code.

    int_count = 0
    for i = 1:size(qb.points)[2]
        x01 = qb.points[1,i]
        y01 = qb.points[2,i]
        x02 = qb.points[3,i]
        y02 = qb.points[4,i]
        x03 = qb.points[5,i]
        y03 = qb.points[6,i]

        # Calculate coefficients of characteristic polynomial
        t2c = (-x01*y11 + x01*y12 + 2*x02*y11 - 2*x02*y12 - x03*y11 + x03*y12 + x11*y01 - 2*x11*y02 + x11*y03 - x12*y01 + 2*x12*y02 - x12*y03)
        t1c = (2*x01*y11 - 2*x01*y12 - 2*x02*y11 + 2*x02*y12 - 2*x11*y01 + 2*x11*y02 + 2*x12*y01 - 2*x12*y02)
        t0c = - x01*y11 + x01*y12 + x11*y01 - x11*y12 - x12*y01 + x12*y11

        a = complex(t2c)
        b = complex(t1c)
        c = complex(t0c)

        # Quadratic function
        discr = b^2 - 4*a*c
        z1 = (-b + sqrt(discr))/(2*a)
        z2 = (-b - sqrt(discr))/(2*a)

        roots = [z1, z2]

        for root in roots
            # Check if their imaginary component is sufficiently small
            if abs(imag(root)) < tol
                t = real(root)
                # Check bounding box on Bezier
                if t <= 1.0 + tol && t >= -tol
                    # Get px, py
                    p_x = (1 - t)^2 * x01 + 2*(1-t)*t * x02 + t^2 * x03
                    p_y = (1 - t)^2 * y01 + 2*(1-t)*t * y02 + t^2 * y03

                    # Check that collision is forward from p1
                    if x11 < x12
                        gx0 = p_x > x11
                    elseif x12 < x11
                        gx0 = p_x < x11
                    else
                        gx0 = true
                    end
                    gy0 = false
                    if gx0
                        if y11 < y11
                            gy0 = (p_y > y11)
                        elseif y11 < y11
                            gy0 = (p_y < y11)
                        else
                            gy0 = true
                        end
                    end

                    if gy0
                        # Collided, in both bounding boxes
                        int_count += 1

                        path_length = sqrt((x11 - p_x)^2 + (y11 - p_y)^2)
                        if path_length < min_length
                            min_length = path_length
                            new_point = [p_x p_y]

                            # Calculate normal
                            dx = 2*(1-t)*(x02 - x01) + 2*t*(x03 - x02)
                            dy = 2*(1-t)*(y02 - y01) + 2*t*(y03 - y02)
                            normal = [-dy dx]
                            normal = normal / norm(normal)
                        end
                    end
                end
            end
        end
    end

    int_count, new_point, min_length, normal
end
