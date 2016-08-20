using LightXML
include("surfaces.jl")
include("geometry.jl")
include("surfaces.jl")
include("spectral_weight.jl")

# =============================================================================
# Approximation functions
# =============================================================================

"""
Calculate angle between u and v (at origin, for use with ellipse)
"""
function angle(u, v)
    sign_val = sign(u[1] * v[2] - u[2] * v[1])

    dotp = dot(u, v)
    denom = norm(u) * norm(v)

    sign_val * acos(dotp/denom)
end

"""
Approximate an ellipse with linear segments of quantity res
"""
function ellipse_approximate(x1, y1, rx, ry, phi_deg, fa, fs, x2, y2, res)
    
    # Why.

    phi = phi_deg / (2*pi)

    # https://www.w3.org/TR/SVG/implnote.html#ArcConversionEndpointToCenter
    # Step 1
    x1p = cos(phi) * ((x1 - x2)/2) + sin(phi) * ((y1 - y2)/2)
    y1p = -sin(phi) * ((x1 - x2)/2) + cos(phi) * ((y1 - y2)/2)

    # Adjust rx, ry to ensure reality
    rx = abs(rx)
    ry = abs(ry)

    lamb = x1p^2/rx^2 + y1p^2/ ry^2
    if lamb > 1
        rx *= sqrt(lamb)
        ry *= sqrt(lamb)
    end

    # Step 2
    coeff = sqrt((rx^2 * ry^2 - rx^2 * y1p^2 - ry^2 * x1p^2)/(rx^2 *y1p^2 + ry^2 * x1p^2))
    if fa == fs
        coeff = -coeff
    end

    cxp = coeff * rx * y1p / ry
    cyp = -coeff * ry * x1p / rx

    # Step 4 (deliberately out of order)
    th1 = angle([1, 0], [(x1p - cxp)/rx, (y1p - cyp)/ry])
    dth = angle([(x1p - cxp)/rx, (y1p - cyp)/ry], [(-x1p - cxp)/rx, (-y1p - cyp)/ry])

    # Demunge it
    if fs == 0 && dth >= 0
        dth -= 2*pi
    end
    if fs == 1 && dth <= 0
        dth += 2*pi
    end

    # Approximate with splines
    ellipse_approximate_spline((x1 + x2)/2, (y1 + y2)/2, rx, ry, phi, th1, th1+dth, res)
end

function ellipse_approximate_spline(cx, cy, rx, ry, phi, t1, t2, res_orig)
    # This uses quite a bit of theory from this excellent article
    # https://pomax.github.io/bezierinfo/#circles_cubic
    # Truthfully, I'd rather not do this, but it makes things so much easier,
    # since the code is highly tuned for Beziers.

    # First, we form the approximation from theta = 0 to theta = (t2 - t1)
    # for a circle
    csplines = []

    # Calculate stepwise resolution
    res = ceil(Int, (t2 - t1) / (2*pi) * res_orig) 
    # Calculate original spline
    dth = (t2 - t1) / res
    x0 = 1
    y0 = 0

    f = 4/3 * tan(dth/4)

    xc1 = 1
    yc1 = f

    xc2 = cos(dth) + f * sin(dth)
    yc2 = sin(dth) - f * cos(dth)

    x1 = cos(dth)
    y1 = sin(dth)

    for i = 1:res
        if csplines == []
            csplines = reshape([x0, y0, xc1, yc1, xc2, yc2, x1, y1], (8,1))
        else
            csplines = hcat(csplines, [x0, y0, xc1, yc1, xc2, yc2, x1, y1])
        end

        # Rotate by dth
        x0 = x1
        y0 = y1

        xc11 = cos(dth)*xc1 - sin(dth)*yc1
        yc1 = sin(dth)*xc1 + cos(dth)*yc1
        xc1 = xc11

        xc21 = cos(dth)*xc2 - sin(dth)*yc2
        yc2 = sin(dth)*xc2 + cos(dth)*yc2
        xc2 = xc21

        x11 = cos(dth)*x1 - sin(dth)*y1
        y1 = sin(dth)*x1 + cos(dth)*y1
        x1 = x11
    end

    # Now, convert into a surface for more transforms
    cb = CBezier(csplines)

    # mat1, rotate from t1 = 0 to t1 = t1
    mat1 = [[cos(t1) -sin(t1) 0];[sin(t1) cos(t1) 0];[0 0 1]]

    # mat2, stretch to match rx
    mat2 = [[rx 0 0];[0 1 0];[0 0 1]]

    # mat3, stretch to match ry
    mat3 = [[1 0 0];[0 ry 0];[0 0 1]]

    # mat4, shift to [cx, cy]
    mat4 = [[1 0 cx];[0 1 cy];[0 0 1]]

    # mat5, rotate last time to final phi
    mat5 = [[cos(phi) -sin(phi) 0];[sin(phi) cos(phi) 0];[0 0 1]]

    mat = mat5*mat4*mat3*mat2*mat1

    apply_transform!(cb, mat)
    cb.points
end

# =============================================================================
# SVG converter components
# =============================================================================

"""
Splits apart the components of a path "d" 
May need adjusting as I haven't tested all files.
"""
function split_svg(command)
    # First, we split by all commands.
    points_str = split(command, r"[a-df-zA-DF-Z]")

    # Convert strings into float lists
    points = []
    for svg_str in points_str[2:end] # Skip empty (pre-M) beginning component
        svg_str = strip(svg_str)
        if svg_str == ""
            # Just append a dummy value
            push!(points, [0,0])
        else
            # Required for Adobe Illustrator output
            # substitutes [float]-[float] to [float],-[float]
            # but not [float] -[float] as required for Inkscape strings
            svg_str = replace(svg_str, r"([^,e ])-", s"\1,-")
            # Required for Inkscape output
            svg_str = replace(svg_str, " ", ",")

            # Split on commas, convert to float array
            point_flt = [float(x) for x in split(svg_str, ",")]
            push!(points, point_flt)
        end
    end

    # Then we extract commands
    commands = []
    for x in matchall(r"[a-df-zA-DF-Z]", command)
        push!(commands, x)
    end

    points, commands
end

"""
Takes an SVG path command and converts to computational geometry components
"""
function path_to_rings(svg_str, bez_res)
    segments = []
    qsplines = []
    csplines = []

    tol = 1.0e-10

    points, commands = split_svg(svg_str)

    origin = [0 0]
    current_point = [0 0]
    last_control_point = [0 0]
    beginning_of_ring = [0 0]
    appended = false

    if svg_str == "z"
        return []
    end

    for i = 1:length(commands)
        appended = false
        com = commands[i]

        if com == "M" || com == "L" || com == "m" || com == "l"
            # move / line with possible lines afterwards
            n_p = Int(length(points[i]) / 2)

            for j = 0:n_p-1
                x = points[i][2*j + 1]
                y = points[i][2*j + 2]

                if com == "m" || com == "l"
                    x += current_point[1]
                    y += current_point[2]
                else
                    origin = [x y]
                end

                # Check for new ring
                if j == 0 && (i == 1 || commands[i-1] == "z" || commands[i-1] == "Z")
                    # New beginning of ring
                    beginning_of_ring = [x y]
                end

                if com == "L" || com == "l" || j != 0
                    # Draw line
                    if segments == []
                        segments = reshape([current_point x y], (1,4))
                    else
                        segments = vcat(segments, [current_point x y])
                    end
                end

                current_point = [x y]
                last_control_point = [x y]
            end
        elseif com == "h" || com == "H"
            # horizontal line
            n_p = Int(length(points[i]))
            for j = 0:n_p-1
                x = points[i][j + 1]
                y = current_point[2]

                if com == "h"
                    # Relative horizontal line
                    x += current_point[1]
                else
                    origin = [x y]
                end

                if segments == []
                    segments = reshape([current_point x y], (1,4))
                else
                    segments = vcat(segments, [current_point x y])
                end

                current_point = [x y]
                last_control_point = [x y]
            end
        elseif com == "v" || com == "V"
            # vertical line
            n_p = Int(length(points[i]))
            for j = 0:n_p-1
                x = current_point[1]
                y = points[i][j + 1]

                if com == "v"
                    # Relative vertical line
                    y += current_point[2]
                else
                    origin = [x y]
                end

                if segments == []
                    segments = reshape([current_point x y], (1,4))
                else
                    segments = vcat(segments, [current_point x y])
                end

                current_point = [x y]
                last_control_point = [x y]
            end
        elseif com == "c" || com == "C"
            # Cubic splines
            n_p = Int(length(points[i]) / 6)

            for j = 0:n_p-1
                x0 = current_point[1]
                y0 = current_point[2]

                xc1 = points[i][6*j+1]
                yc1 = points[i][6*j+2]
                xc2 = points[i][6*j+3]
                yc2 = points[i][6*j+4]
                x1 = points[i][6*j+5]
                y1 = points[i][6*j+6]

                if com == "c"
                    # Relative cubic splines
                    xc1 += x0
                    xc2 += x0
                    x1 += x0

                    yc1 += y0
                    yc2 += y0
                    y1 += y0
                else
                    origin = [x1 y1]
                end

                if csplines == []
                    csplines = reshape([x0, y0, xc1, yc1, xc2, yc2, x1, y1], (8,1))
                else
                    csplines = hcat(csplines, [x0, y0, xc1, yc1, xc2, yc2, x1, y1])
                end

                current_point = [x1 y1]
                last_control_point = [xc2 yc2]
            end
        elseif com == "s" || com == "S"
            # Cubic smooth splines
            n_p = Int(length(points[i]) / 4)

            for j = 0:n_p-1
                x0 = current_point[1]
                y0 = current_point[2]

                xc1 = 2*x0 - last_control_point[1]
                yc1 = 2*y0 - last_control_point[2]
                xc2 = points[i][4*j+1]
                yc2 = points[i][4*j+2]
                x1 = points[i][4*j+3]
                y1 = points[i][4*j+4]

                if com == "s"
                    # Relative cubic splines
                    xc2 += x0
                    x1 += x0

                    yc2 += y0
                    y1 += y0
                else
                    origin = [x1 y1]
                end

                if csplines == []
                    csplines = reshape([x0, y0, xc1, yc1, xc2, yc2, x1, y1], (8,1))
                else
                    csplines = hcat(csplines, [x0, y0, xc1, yc1, xc2, yc2, x1, y1])
                end

                current_point = [x1 y1]
                last_control_point = [xc2 yc2]
            end
        elseif com == "q" || com == "Q"
            # Quadratic splines
            n_p = Int(length(points[i]) / 4)

            for j = 0:n_p-1
                x0 = current_point[1]
                y0 = current_point[2]

                xc1 = points[i][4*j+1]
                yc1 = points[i][4*j+2]
                x1 = points[i][4*j+3]
                y1 = points[i][4*j+4]

                if com == "q"
                    # Relative quadratic splines
                    xc1 += x0
                    x1 += x0

                    yc1 += y0
                    y1 += y0
                else
                    origin = [x1 y1]
                end

                if qsplines == []
                    qsplines = reshape([x0, y0, xc1, yc1, x1, y1], (6,1))
                else
                    qsplines = hcat(qsplines, [x0, y0, xc1, yc1, x1, y1])
                end

                current_point = [x1 y1]
                last_control_point = [xc1 yc1]
            end
        elseif com == "t" || com == "T"
            # Smooth quadratic splines
            n_p = Int(length(points[i]) / 2)

            for j = 0:n_p-1
                x0 = current_point[1]
                y0 = current_point[2]

                xc1 = 2*x0 - last_control_point[1]
                yc1 = 2*y0 - last_control_point[2]
                x1 = points[i][2*j+1]
                y1 = points[i][2*j+2]

                if com == "t"
                    # Relative quadratic splines
                    x1 += x0

                    y1 += y0
                else
                    origin = [x1 y1]
                end

                if qsplines == []
                    qsplines = reshape([x0, y0, xc1, yc1, x1, y1], (6,1))
                else
                    qsplines = hcat(qsplines, [x0, y0, xc1, yc1, x1, y1])
                end

                current_point = [x1 y1]
                last_control_point = [xc1 yc1]
            end
        elseif com == "a" || com == "A"
            # Ellipse
            n_p = Int(length(points[i]) / 7)

            for j = 0:n_p-1
                x0 = current_point[1]
                y0 = current_point[2]

                rx = points[i][7*j + 1]
                ry = points[i][7*j + 2]
                x_rot = points[i][7*j + 3]
                large_arc = Int(points[i][7*j + 4])
                sweep = Int(points[i][7*j + 5])
                x1 = points[i][7*j + 6]
                y1 = points[i][7*j + 7]

                if com == "a"
                    # Relative ellipse
                    x1 += x0

                    y1 += y0
                else
                    origin = [x1 y1]
                end

                # Ellipses are hard, so we approximate them.
                points_ellipse = ellipse_approximate(x0, y0, rx, ry, x_rot, large_arc, sweep, x1, y1, bez_res)

                if csplines == []
                    csplines = points_ellipse
                else
                    csplines = hcat(csplines, points_ellipse)
                end

                current_point = [x1 y1]
                last_control_point = [x1 y1]
            end
        elseif com == "z" || com == "Z"
            # Terminate ring, append to ring list.
            appended = true
            if norm(current_point - beginning_of_ring) > tol
                if segments == []
                    segments = reshape([current_point beginning_of_ring], (1,4))
                else
                    segments = vcat(segments, [current_point beginning_of_ring])
                end
            end

            # Frankly, I just don't understand this part.
            # A capital command moves the origin from [0 0]
            # a lower case m at the beginning of a new path uses this origin
            # unless it's never set, at which point it's relative to the end of the last path.
            # ???
            if origin != [0 0]
                current_point = origin
            else
                current_point = beginning_of_ring
            end
        end
    end

    # Close ring if z is never called for some reason
    if !appended
        if norm(current_point - beginning_of_ring) > tol
            if segments == []
                segments = reshape([current_point beginning_of_ring], (1,4))
            else
                segments = vcat(segments, [current_point beginning_of_ring])
            end
        end
    end

    # If empty, replace with dummy 2d array
    if segments == []
        segments = reshape(zeros(0), (0,0))
    end

    if qsplines == []
        qsplines = reshape(zeros(0), (0,0))
    end

    if csplines == []
        csplines = reshape(zeros(0), (0,0))
    end
    la = LineArray(segments)
    qs = QBezier(qsplines)
    cs = CBezier(csplines)

    la, qs, cs
end

"""
Reads the string of a transformation matrix and forms the matrix

Supports only one transform per string currently
"""
function transform_matrix(mat_str)
    # Split off type, and numbers
    t1 = split(mat_str, "(")
    mat_type = t1[1]
    mat_data = split(split(t1[2], ")")[1], ",")
    if length(mat_data) == 1
        # Try again
        mat_data = split(split(t1[2], ")")[1], " ")
    end
    mat_data = [float(s) for s in mat_data]
    if mat_type == "matrix"
        a = mat_data[1]
        b = mat_data[2]
        c = mat_data[3]
        d = mat_data[4]
        e = mat_data[5]
        f = mat_data[6]
    elseif mat_type == "translate"
        x = mat_data[1]
        if length(mat_data) == 2
            y = mat_data[2]
        else
            y = 0
        end
        a = 1
        b = 0
        c = 0
        d = 1
        e = x
        f = y
    elseif mat_type == "scale"
        x = mat_data[1]
        if length(mat_data) == 2
            y = mat_data[2]
        else
            y = 0
        end
        a = x
        b = 0
        c = 0
        d = y
        e = 0
        f = 0
    elseif mat_type == "rotate"
        a = mat_data[1] * pi / 180
        if length(mat_data) == 2
            x = mat_data[2]
        else
            x = 0
        end

        if length(mat_data) == 3
            y = mat_data[3]
        else
            y = 0
        end

        # Compound transform possible
        m1 = [[1 0 x];[0 1 y];[0 0 1]]
        m2 = [[cos(a) -sin(a) 0];[sin(a) cos(a) 0];[0 0 1]]
        m3 = [[1 0 -x];[0 1 -y];[0 0 1]]
        return m1 * m2 * m3
    elseif mat_type == "skewX"
        a = mat_data[1] * pi / 180

        a = 1
        b = 0
        c = tan(a)
        d = 1
        e = 0
        f = 0
    elseif mat_type == "skewY"
        a = mat_data[1] * pi / 180

        a = 1
        b = tan(a)
        c = 0
        d = 1
        e = 0
        f = 0
    end

    [[a c e];[b d f];[0 0 1]]
end

"""
Given an XML path, loads x from either XML or from "style"
"""
function get_term(e, term)
    has_term = false
    # Check style attribute
    if has_attribute(e, "style")
        style_str = attributes_dict(e)["style"]
        style_parts = split(style_str, ";")
        for parts in style_parts
            part = split(parts, ":")
            if strip(part[1]) == term
                result = strip(part[2])
                has_term = true
                break
            end
        end
    elseif has_attribute(e, term)
        result = strip(attributes_dict(e)[term])
        has_term = true
    end

    if !has_term
        result = "none"
    end
    result
end

"""
Converts hex string to RGB coordinates
"""
function hex_to_rgb(str)
    try
        r = parse(Int, str[2:3],16)
        g = parse(Int, str[4:5],16)
        b = parse(Int, str[6:7],16)
        return [r, g, b]
    catch
        r = 255
        g = 255
        b = 255
        return [255, 255, 255]
    end
end

"""
Creates a rectangle path
"""
function form_rect(x, y, width, height, rx, ry, res)
    
    # Fix rx, ry components
    if rx > 0 && ry == 0
        ry = rx
    end

    if ry > 0 && rx == 0
        rx = ry
    end

    if rx > width / 2
        rx = width / 2
    end

    if ry > height / 2
        ry = height / 2
    end

    if rx < 0
        rx = 0
    end

    if ry < 0
        ry = 0
    end

    if rx != 0 || ry != 0
        curved_corners = true
    else
        curved_corners = false
    end

    csplines = []
    segments = []

    # Process if curved
    if curved_corners
        # Process counter clockwise from curve at origin
        # C1
        cx = x + rx
        cy = y + ry
        th1 = pi
        th2 = 3/2 * pi
        csplines = ellipse_approximate_spline(cx, cy, rx, ry, 0, th1, th2, res)

        # C2
        cx = x + width - rx
        cy = y + ry
        th1 = 3/2 * pi
        th2 = 2 * pi
        csplines = hcat(csplines, ellipse_approximate_spline(cx, cy, rx, ry, 0, th1, th2, res))

        # C3
        cx = x + width - rx
        cy = y + height - ry
        th1 = 0
        th2 = pi/2
        csplines = hcat(csplines, ellipse_approximate_spline(cx, cy, rx, ry, 0, th1, th2, res))

        # C4
        cx = x + rx
        cy = y + height - ry
        th1 = pi/2
        th2 = pi
        csplines = hcat(csplines, ellipse_approximate_spline(cx, cy, rx, ry, 0, th1, th2, res))
    else
        csplines = reshape(zeros(0), (0,0))
    end
    # Get lines
    # L1
    x0 = x + rx
    x1 = x + width - rx
    yl = y
    segments = reshape([x0 yl x1 yl], (1,4))

    # L2
    y0 = y + ry
    y1 = y + height - ry
    xl = x + width
    segments = vcat(segments, [xl y0 xl y1])

    # L3
    x0 = x + width - rx
    x1 = x + rx
    yl = y + height
    segments = vcat(segments, [x0 yl x1 yl])

    # L4
    y0 = y + height - ry
    y1 = y + ry
    xl = x
    segments = vcat(segments, [xl y0 xl y1])

    # Generate quadratic objects
    qsplines = reshape(zeros(0), (0,0))

    ls = LineArray(segments)
    qs = QBezier(qsplines)
    cs = CBezier(csplines)

    ls, qs, cs
end

"""
Loads a rectangle, converts to cells
"""
function parse_rect(e, res)
    x = float(attributes_dict(e)["x"])
    y = float(attributes_dict(e)["y"])
    width = float(attributes_dict(e)["width"])
    height = float(attributes_dict(e)["height"])

    curved_corners = false

    if has_attribute(e, "rx")
        # Oh boy
        rx = float(attributes_dict(e)["rx"])
        curved_corners = true
    else
        rx = 0
    end

    if has_attribute(e, "ry")
        ry = float(attributes_dict(e)["ry"])
        curved_corners = true
    else
        ry = 0
    end

    ls = []
    qs = []
    cs = []
    colors = []

    # Check if path has stroke width
    sw = get_term(e, "stroke-width")
    fill = get_term(e, "fill")
    stroke = get_term(e, "stroke")
    if sw != "none"
        sw = float(sw)
        # Oh boy
        if rx == 0 && ry == 0
            rx_low = 0
            ry_low = 0
            rx_high = 0
            ry_high = 0
        elseif rx == 0
            rx_low = ry - sw/2
            ry_low = ry - sw/2
            rx_high = ry + sw/2
            ry_high = ry + sw/2
        elseif ry == 0
            rx_low = rx - sw/2
            ry_low = rx - sw/2
            rx_high = rx + sw/2
            ry_high = rx + sw/2
        else
            rx_low = rx - sw/2
            ry_low = ry - sw/2
            rx_high = rx + sw/2
            ry_high = ry + sw/2
        end

        # Ensure reality
        if width - sw < 0 || height - sw < 0
            fill = "none"
            ls2 = LineArray(reshape(zeros(0), (0,0)))
            qs2 = QBezier(reshape(zeros(0), (0,0)))
            cs2 = CBezier(reshape(zeros(0), (0,0)))
        else
            # Inner path
            ls2, qs2, cs2 = form_rect(x + sw/2, y + sw/2, width - sw, height - sw, rx_low, ry_low, res)
        end

        # Outer path
        ls1, qs1, cs1 = form_rect(x - sw/2, y - sw/2, width + sw, height + sw, rx_high, ry_high, res)

        # Check if fill is none
        if fill == "none"
            # Merge inner, outer paths
            ls_merged = LineArray( vcat(ls1.points, ls2.points))
            qs_merged = QBezier( hcat(qs1.points, qs2.points))
            cs_merged = CBezier( hcat(cs1.points, cs2.points))

            # Add inner path
            push!(ls, ls_merged)
            push!(qs, qs_merged)
            push!(cs, cs_merged)
            push!(colors, hex_to_rgb(stroke))
        else
            # Layer fill on top of path around
            push!(ls, ls1)
            push!(qs, qs1)
            push!(cs, cs1)
            push!(colors, hex_to_rgb(stroke))

            push!(ls, ls2)
            push!(qs, qs2)
            push!(cs, cs2)
            push!(colors, hex_to_rgb(fill))
        end
    else
        fill = get_term(e, "fill")
        ls1, qs1, cs1 = form_rect(x, y, width, height, rx, ry, res)
        push!(ls, ls1)
        push!(qs, qs1)
        push!(cs, cs1)
        push!(colors, hex_to_rgb(fill))
    end

    # Get transform matrix
    if has_attribute(e, "transform")
        # Apply transform
        transform_str = attributes_dict(e)["transform"]
        transform_mat = transform_matrix(transform_str)
        [apply_transform!(ls_p, transform_mat) for ls_p in ls]
        [apply_transform!(qs_p, transform_mat) for qs_p in qs]
        [apply_transform!(cs_p, transform_mat) for cs_p in cs]
    end

    ls, qs, cs, colors
end

"""
Loads an ellipse, converts to cells
"""
function parse_ellipse(e, res)
    cx = float(attributes_dict(e)["cx"])
    cy = float(attributes_dict(e)["cy"])
    rx = float(attributes_dict(e)["rx"])
    ry = float(attributes_dict(e)["ry"])

    sw = get_term(e, "stroke-width")
    fill = get_term(e, "fill")
    stroke = get_term(e, "stroke")

    ls = []
    qs = []
    cs = []
    colors = []

    # Check if path has stroke width
    sw = get_term(e, "stroke-width")
    if sw != "none"
        sw = float(sw)

        # Outer path
        csplines = ellipse_approximate_spline(cx, cy, rx + sw/2, ry + sw/2, 0, 0, 2*pi, res)
        cs1 = CBezier(csplines)
        qs1 = QBezier(reshape(zeros(0), (0,0)))
        ls1 = LineArray(reshape(zeros(0), (0,0)))

        # Ensure reality
        if rx - sw/2 < 0 || ry - sw/2 < 0
            fill = "none"
            ls2 = LineArray(reshape(zeros(0), (0,0)))
            qs2 = QBezier(reshape(zeros(0), (0,0)))
            cs2 = CBezier(reshape(zeros(0), (0,0)))
        else
            csplines2 = ellipse_approximate_spline(cx, cy, rx - sw/2, ry - sw/2, 0, 0, 2*pi, res)
            cs2 = CBezier(csplines2)
            qs2 = QBezier(reshape(zeros(0), (0,0)))
            ls2 = LineArray(reshape(zeros(0), (0,0)))
        end

        # Check if fill is none
        if fill == "none"
            # Merge inner, outer paths
            ls_merged = LineArray( vcat(ls1.points, ls2.points))
            qs_merged = QBezier( hcat(qs1.points, qs2.points))
            cs_merged = CBezier( hcat(cs1.points, cs2.points))

            # Add inner path
            push!(ls, ls_merged)
            push!(qs, qs_merged)
            push!(cs, cs_merged)
            push!(colors, hex_to_rgb(stroke))
        else
            # Layer fill on top of path around
            push!(ls, ls1)
            push!(qs, qs1)
            push!(cs, cs1)
            push!(colors, hex_to_rgb(stroke))

            push!(ls, ls2)
            push!(qs, qs2)
            push!(cs, cs2)
            push!(colors, hex_to_rgb(fill))
        end
    else
        csplines = ellipse_approximate_spline(cx, cy, rx, ry, 0, 0, 2*pi, res)
        cs1 = CBezier(csplines)
        qs1 = QBezier(reshape(zeros(0), (0,0)))
        ls1 = LineArray(reshape(zeros(0), (0,0)))

        push!(ls, ls1)
        push!(qs, qs1)
        push!(cs, cs1)
        push!(colors, hex_to_rgb(fill))
    end

    # Get transform matrix
    if has_attribute(e, "transform")
        # Apply transform
        transform_str = attributes_dict(e)["transform"]
        transform_mat = transform_matrix(transform_str)
        [apply_transform!(ls_p, transform_mat) for ls_p in ls]
        [apply_transform!(qs_p, transform_mat) for qs_p in qs]
        [apply_transform!(cs_p, transform_mat) for cs_p in cs]
    end

    ls, qs, cs, colors
end

"""
Loads a circle, converts to cells
"""
function parse_circle(e, res)
    cx = float(attributes_dict(e)["cx"])
    cy = float(attributes_dict(e)["cy"])
    r = float(attributes_dict(e)["r"])

    sw = get_term(e, "stroke-width")
    fill = get_term(e, "fill")
    stroke = get_term(e, "stroke")

    ls = []
    qs = []
    cs = []
    colors = []

    # Check if path has stroke width
    sw = get_term(e, "stroke-width")
    if sw != "none"
        sw = float(sw)

        # Outer path
        csplines = ellipse_approximate_spline(cx, cy, r + sw/2, r + sw/2, 0, 0, 2*pi, res)
        cs1 = CBezier(csplines)
        qs1 = QBezier(reshape(zeros(0), (0,0)))
        ls1 = LineArray(reshape(zeros(0), (0,0)))

        # Ensure reality
        if r - sw/2 < 0
            fill = "none"
            ls2 = LineArray(reshape(zeros(0), (0,0)))
            qs2 = QBezier(reshape(zeros(0), (0,0)))
            cs2 = CBezier(reshape(zeros(0), (0,0)))
        else
            csplines2 = ellipse_approximate_spline(cx, cy, r - sw/2, r - sw/2, 0, 0, 2*pi, res)
            cs2 = CBezier(csplines2)
            qs2 = QBezier(reshape(zeros(0), (0,0)))
            ls2 = LineArray(reshape(zeros(0), (0,0)))
        end

        # Check if fill is none
        if fill == "none"
            # Merge inner, outer paths
            ls_merged = LineArray( vcat(ls1.points, ls2.points))
            qs_merged = QBezier( hcat(qs1.points, qs2.points))
            cs_merged = CBezier( hcat(cs1.points, cs2.points))

            # Add inner path
            push!(ls, ls_merged)
            push!(qs, qs_merged)
            push!(cs, cs_merged)
            push!(colors, hex_to_rgb(stroke))
        else
            # Layer fill on top of path around
            push!(ls, ls1)
            push!(qs, qs1)
            push!(cs, cs1)
            push!(colors, hex_to_rgb(stroke))

            push!(ls, ls2)
            push!(qs, qs2)
            push!(cs, cs2)
            push!(colors, hex_to_rgb(fill))
        end
    else
        csplines = ellipse_approximate_spline(cx, cy, r, r, 0, 0, 2*pi, res)
        cs1 = CBezier(csplines)
        qs1 = QBezier(reshape(zeros(0), (0,0)))
        ls1 = LineArray(reshape(zeros(0), (0,0)))

        push!(ls, ls1)
        push!(qs, qs1)
        push!(cs, cs1)
        push!(colors, hex_to_rgb(fill))
    end

    # Get transform matrix
    if has_attribute(e, "transform")
        # Apply transform
        transform_str = attributes_dict(e)["transform"]
        transform_mat = transform_matrix(transform_str)
        [apply_transform!(ls_p, transform_mat) for ls_p in ls]
        [apply_transform!(qs_p, transform_mat) for qs_p in qs]
        [apply_transform!(cs_p, transform_mat) for cs_p in cs]
    end

    ls, qs, cs, colors
end

"""
Find all paths and subsidiary g in the XML heirarchy
"""
function parse_g(root, bez_res)
    ls = []
    qs = []
    cs = []
    c_array = []
    for c in child_nodes(root)
        if is_elementnode(c)
            e = XMLElement(c)
            needs_fill = false

            if LightXML.name(c) == "g"
                # Nested g
                ls_g, qs_g, cs_g, colors_g = parse_g(e, bez_res)
                [push!(ls, ls_i) for ls_i in ls_g]
                [push!(qs, qs_i) for qs_i in qs_g]
                [push!(cs, cs_i) for cs_i in cs_g]
                [push!(c_array, colors_i) for colors_i in colors_g]
            elseif LightXML.name(c) == "path"
                # Actual path. Get d

                path_string = attributes_dict(e)["d"]
                ls_p, qs_p, cs_p = path_to_rings(path_string, bez_res)
                if has_attribute(e, "transform")
                    # Apply transform
                    transform_str = attributes_dict(e)["transform"]
                    transform_mat = transform_matrix(transform_str)
                    apply_transform!(ls_p, transform_mat)
                    apply_transform!(qs_p, transform_mat)
                    apply_transform!(cs_p, transform_mat)
                end
                push!(ls, ls_p)
                push!(qs, qs_p)
                push!(cs, cs_p)
                
                needs_fill = true
            elseif LightXML.name(c) == "ellipse"
                ls_e, qs_e, cs_e, colors_e = parse_ellipse(e, bez_res)
                [push!(ls, ls_i) for ls_i in ls_e]
                [push!(qs, qs_i) for qs_i in qs_e]
                [push!(cs, cs_i) for cs_i in cs_e]
                [push!(c_array, colors_i) for colors_i in colors_e]
            elseif LightXML.name(c) == "circle"
                ls_e, qs_e, cs_e, colors_e = parse_circle(e, bez_res)
                [push!(ls, ls_i) for ls_i in ls_e]
                [push!(qs, qs_i) for qs_i in qs_e]
                [push!(cs, cs_i) for cs_i in cs_e]
                [push!(c_array, colors_i) for colors_i in colors_e]
            elseif LightXML.name(c) == "rect"
                ls_e, qs_e, cs_e, colors_e = parse_rect(e, bez_res)
                [push!(ls, ls_i) for ls_i in ls_e]
                [push!(qs, qs_i) for qs_i in qs_e]
                [push!(cs, cs_i) for cs_i in cs_e]
                [push!(c_array, colors_i) for colors_i in colors_e]
            end

            if needs_fill
                # Check for a fill
                has_fill = false
                if has_attribute(e, "style")
                    style_str = attributes_dict(e)["style"]
                    style_parts = split(style_str, ";")
                    for parts in style_parts
                        if contains(parts, "fill:")
                            c_str = split(parts, ":")[2]
                            rgb = hex_to_rgb(c_str)
                            push!(c_array, rgb)
                            has_fill = true
                            break
                        end
                    end
                elseif has_attribute(e, "fill")
                    fill_str = attributes_dict(e)["fill"]
                    rgb = hex_to_rgb(fill_str)
                    has_fill
                end

                if !has_fill
                    push!(c_array, [255, 255, 255])
                end
            end
        end
    end

    if has_attribute(root, "transform")
        # Apply transform
        transform_str = attributes_dict(root)["transform"]
        transform_mat = transform_matrix(transform_str)
        for i = 1:length(ls)
            apply_transform!(ls[i], transform_mat)
            apply_transform!(qs[i], transform_mat)
            apply_transform!(cs[i], transform_mat)
        end
    end
    ls, qs, cs, c_array
end

"""
Load an XML file and convert it into the components for a geometry
"""
function svg_to_geo(filename, bez_res)
    xdoc = parse_file(filename)
    xroot = LightXML.root(xdoc)

    ls = []
    qs = []
    cs = []
    c_array = []

    for c in child_nodes(xroot)
        if is_elementnode(c)
            e = XMLElement(c)
            if LightXML.name(c) == "g"
                ls_g, qs_g, cs_g, colors_g = parse_g(e, bez_res)
                [push!(ls, ls_i) for ls_i in ls_g]
                [push!(qs, qs_i) for qs_i in qs_g]
                [push!(cs, cs_i) for cs_i in cs_g]
                [push!(c_array, colors_i) for colors_i in colors_g]
            end
        end
    end

    ls, qs, cs, c_array
end

"""
Helper function to load geometry
"""
function load_geometry(filename, materials, res, color_mode, mat = [[1 0 0];[0 1 0];[0 0 1]])
    ls, qs, cs, color_array = svg_to_geo(filename, res)

    # Apply transformation matrix
    for i = 1:length(ls)
        apply_transform!(ls[i], mat)
        apply_transform!(qs[i], mat)
        apply_transform!(cs[i], mat)
    end

    # Form geometry
    c_array = []
    c_count = 0
    for i = 1:length(ls)
        # Deal with colors
        if color_mode == 1
            scale = 1.0
        elseif color_mode == 2
            ~, scale = spec_weight(color_array[i])
        else
            scale, ~ = spec_weight(color_array[i])
        end
        if isa(materials, Array)
            mat1 = deepcopy(materials[i])
            mat1.scale = scale
            push!(c_array, Cell(ls[i], qs[i], cs[i], mat1))
        else
            mat1 = deepcopy(materials)
            mat1.scale = scale
            push!(c_array, Cell(ls[i], qs[i], cs[i], mat1))
        end
    end

    geo = Geometry(c_array)
    geo
end