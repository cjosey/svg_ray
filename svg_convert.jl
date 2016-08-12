using LightXML
include("surfaces.jl")
include("geometry.jl")

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

    # Traverse path
    t = linspace(th1, th1 + dth, res + 1)

    # Calculate xp, yp
    xp = rx * cos(t)
    yp = ry * sin(t)

    # Transform it back to real coordinates (step 3)
    x = cos(phi) * xp - sin(phi) * yp + (x1 + x2)/2
    y = sin(phi) * xp + cos(phi) * yp + (y1 + y2)/2

    hcat(x[1:end-1], y[1:end-1], x[2:end], y[2:end])
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

                if segments == []
                    segments = points_ellipse
                else
                    segments = vcat(segments, points_ellipse)
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
Find all paths and subsidiary g in the XML heirarchy
"""
function parse_g(root, bez_res)
    ls = []
    qs = []
    cs = []
    for c in child_nodes(root)
        if is_elementnode(c)
            e = XMLElement(c)
            if LightXML.name(c) == "g"
                # Nested g
                ls_g, qs_g, cs_g = parse_g(e, bez_res)
                [push!(ls, ls_i) for ls_i in ls_g]
                [push!(qs, qs_i) for qs_i in qs_g]
                [push!(cs, cs_i) for cs_i in cs_g]
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
    ls, qs, cs
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

    for c in child_nodes(xroot)
        if is_elementnode(c)
            e = XMLElement(c)
            if LightXML.name(c) == "g"
                ls_g, qs_g, cs_g = parse_g(e, bez_res)
                [push!(ls, ls_i) for ls_i in ls_g]
                [push!(qs, qs_i) for qs_i in qs_g]
                [push!(cs, cs_i) for cs_i in cs_g]
            end
        end
    end

    ls, qs, cs
end

"""
Helper function to load geometry
"""
function load_geometry(filename, geo_filename, materials, res, mat = [[1 0 0];[0 1 0];[0 0 1]])
    ls, qs, cs = svg_to_geo(filename, res)

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
        if isa(materials, Array)
            push!(c_array, Cell(ls[i], qs[i], cs[i], materials[i]))
        else
            push!(c_array, Cell(ls[i], qs[i], cs[i], materials))
        end
    end

    geo = Geometry(c_array)
    plot_geometry(geo, geo_filename)
    geo
end