import Base.+

# =============================================================================
# Types
# =============================================================================

"""
A tally object.  Contains the data, as well as information on dimensions and
wavelengths stored.
"""
type Tally
    data::Array{Float32, 3} # 3D array of dimension [x y color]
    x::Int                  # Size in x
    y::Int                  # Size in y
    scaling_factor::Number  # Ratio of actual dimensions to tally dimensions.  If > 1, increases tally size.
end

"""
A constructor for the tally with some defaults
"""
function Tally(x::Number, y::Number, scaling_factor=1)
    x = round(Int, x * scaling_factor)
    y = round(Int, y * scaling_factor)

    data = zeros(Float32, x, y, 3)

    Tally(data, x, y, scaling_factor)
end

"""
Add two tallies together
"""
function +(a::Tally, b::Tally)
    # Generate merged
    temp = a.data + b.data
    Tally(temp, a.x, a.y, a.scaling_factor)
end

"""
Tally a photon's path through the geometry
"""
function tally_path(t, path, wl, scale)
    # Calculate index in CIE table
    n_low = 380
    n_high = 780
    n_w = 80
    wl_ind = floor(Int, (wl - n_low) / (n_high - n_low) * (n_w-1)) + 1

    # Extract WL to CIE XYZ conversion
    cie_color_match = [[0.0014 0.0000 0.0065]; [0.0022 0.0001 0.0105]; [0.0042 0.0001 0.0201];
              [0.0076 0.0002 0.0362]; [0.0143 0.0004 0.0679]; [0.0232 0.0006 0.1102];
              [0.0435 0.0012 0.2074]; [0.0776 0.0022 0.3713]; [0.1344 0.0040 0.6456]; 
              [0.2148 0.0073 1.0391]; [0.2839 0.0116 1.3856]; [0.3285 0.0168 1.6230]; 
              [0.3483 0.0230 1.7471]; [0.3481 0.0298 1.7826]; [0.3362 0.0380 1.7721]; 
              [0.3187 0.0480 1.7441]; [0.2908 0.0600 1.6692]; [0.2511 0.0739 1.5281]; 
              [0.1954 0.0910 1.2876]; [0.1421 0.1126 1.0419]; [0.0956 0.1390 0.8130]; 
              [0.0580 0.1693 0.6162]; [0.0320 0.2080 0.4652]; [0.0147 0.2586 0.3533]; 
              [0.0049 0.3230 0.2720]; [0.0024 0.4073 0.2123]; [0.0093 0.5030 0.1582]; 
              [0.0291 0.6082 0.1117]; [0.0633 0.7100 0.0782]; [0.1096 0.7932 0.0573]; 
              [0.1655 0.8620 0.0422]; [0.2257 0.9149 0.0298]; [0.2904 0.9540 0.0203]; 
              [0.3597 0.9803 0.0134]; [0.4334 0.9950 0.0087]; [0.5121 1.0000 0.0057]; 
              [0.5945 0.9950 0.0039]; [0.6784 0.9786 0.0027]; [0.7621 0.9520 0.0021]; 
              [0.8425 0.9154 0.0018]; [0.9163 0.8700 0.0017]; [0.9786 0.8163 0.0014]; 
              [1.0263 0.7570 0.0011]; [1.0567 0.6949 0.0010]; [1.0622 0.6310 0.0008]; 
              [1.0456 0.5668 0.0006]; [1.0026 0.5030 0.0003]; [0.9384 0.4412 0.0002]; 
              [0.8544 0.3810 0.0002]; [0.7514 0.3210 0.0001]; [0.6424 0.2650 0.0000]; 
              [0.5419 0.2170 0.0000]; [0.4479 0.1750 0.0000]; [0.3608 0.1382 0.0000]; 
              [0.2835 0.1070 0.0000]; [0.2187 0.0816 0.0000]; [0.1649 0.0610 0.0000]; 
              [0.1212 0.0446 0.0000]; [0.0874 0.0320 0.0000]; [0.0636 0.0232 0.0000]; 
              [0.0468 0.0170 0.0000]; [0.0329 0.0119 0.0000]; [0.0227 0.0082 0.0000]; 
              [0.0158 0.0057 0.0000]; [0.0114 0.0041 0.0000]; [0.0081 0.0029 0.0000]; 
              [0.0058 0.0021 0.0000]; [0.0041 0.0015 0.0000]; [0.0029 0.0010 0.0000]; 
              [0.0020 0.0007 0.0000]; [0.0014 0.0005 0.0000]; [0.0010 0.0004 0.0000]; 
              [0.0007 0.0002 0.0000]; [0.0005 0.0002 0.0000]; [0.0003 0.0001 0.0000]; 
              [0.0002 0.0001 0.0000]; [0.0002 0.0001 0.0000]; [0.0001 0.0000 0.0000]; 
              [0.0001 0.0000 0.0000]; [0.0001 0.0000 0.0000]]'

    xyz_val = cie_color_match[:, wl_ind]


    for i = 1:size(path)[1] - 1
        p1 = path[i,:]
        p2 = path[i+1,:]

        tally_line_hq(t, p1, p2, xyz_val, scale[i])
    end
end

"""
Tally line using high quality track length estimation on a square grid

There is no low quality version implemented yet.
"""
function tally_line_hq(t, p1, p2, xyz_val, scale)
    # Tally a line from point p1 to p2.
    # Uses "exact" tallying (track length estimators on square grid)

    x1 = p1 * t.scaling_factor
    x2 = p2 * t.scaling_factor
    uv = x2 - x1

    uv = uv / norm(uv)

    if uv[1] > 0
        # Moving forward on x
        x_dir = 1
    else
        x_dir = -1
    end

    if uv[2] > 0
        # Moving forward on y
        y_dir = 1
    else
        y_dir = -1
    end

    # Move point to surface of tally if needed
    if x1[1] < 1 && x_dir == 1
        # Solve for distance
        distance = (1 - x1[1]) / uv[1]
        x1[1] = 1
        x1[2] += uv[2] * distance
    end 

    if x1[1] > t.x && x_dir == -1
        # Solve for distance
        distance = (t.x - x1[1]) / uv[1]
        x1[1] = t.x
        x1[2] += uv[2] * distance
    end 

    if x1[2] < 1 && y_dir == 1
        # Solve for distance
        distance = (1 - x1[2]) / uv[2]
        x1[1] += uv[1] * distance
        x1[2] = 1
    end 

    if x1[2] > t.y && y_dir == -1
        # Solve for distance
        distance = (t.y - x1[2]) / uv[2]
        x1[1] += uv[1] * distance
        x1[2] = t.y
    end 

    # Are we still outside?  If so, give up.
    if x1[1] < 1 || x1[1] > t.x || x1[2] < 1 || x1[2] > t.y
        return
    end

    x_ind = floor(Int, x1[1])
    y_ind = floor(Int, x1[2])

    # Tally ray
    done = false
    while !done
        # Find length in current pixel
        if x_dir == 1
            dx = (1 + x_ind - x1[1]) / uv[1]
        else
            dx = (x_ind - x1[1]) / uv[1]
        end
        if y_dir == 1
            dy = (1 + y_ind - x1[2]) / uv[2]
        else
            dy = (y_ind - x1[2]) / uv[2]
        end

        if dy == -Inf
            dy = Inf
        end

        if dx == -Inf
            dx = Inf
        end

        remain = sqrt((x1[1] - x2[1])^2 + (x1[2] - x2[2])^2)

        dt = min(dx, dy, remain)

        if dt == remain
            done = true
        end

        dt *= scale * xyz_val

        if isnan(dt[1]) || isnan(dt[2]) || isnan(dt[3])
            dt = [0, 0, 0]
        end

        # Tally
        t.data[x_ind, y_ind, :] += dt

        # println(x_ind," ", y_ind, " ", dx, " ", dy, " ", remain, " ", x1, " ", uv, )

        # Move
        if dx > dy
            y_ind += y_dir
            x1 += dy * uv
        else
            x_ind += x_dir
            x1 += dx * uv
        end

        # If hit wall, delete
        if x_ind > t.x || y_ind > t.y || x_ind < 1 || y_ind < 1
            done = true
        end
    end
end

"""
Converts HDR XYZ tally to HDR RGB tally
"""
function tally_to_rgb(t)
    x = t.x
    y = t.y

    # Convert to RGB
    t2 = Tally(zeros(Float32, x, y, 3), x, y, t.scaling_factor)
    
    # Source http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
    # sRGB D65 matrix
    invmat = [[3.2404542 -1.5371385 -0.4985314];[-0.9692660  1.8760108  0.0415560];[0.0556434 -0.2040259  1.0572252]]

    for i = 1:x
        for j = 1:y
            t2.data[i,j,:] = invmat * reshape(t.data[i,j,:],3)
        end
    end

    t2
end

"""
Plots using imshow an RGB tally
"""
function plot_tally_rgb(t, filename, scale, gamma)
    data = t.data

    # data -= minimum(data)

    data /= maximum(data) / scale

    data[data .< 0] = 0
    data[data .> 1] = 1

    data .^= (1/gamma)

    imshow(data, interpolation="none")
    savefig(filename, dpi=300)
end

"""
Exports an image using scipy.  Supports at least png, pdf.
"""
function save_tally_rgb(t, filename, scale)

    @pyimport scipy.misc as sm

    data = t.data

    data /= maximum(data) / scale

    data[data .< 0] = 0
    data[data .> 1] = 1

    data_new = zeros(t.y, t.x, 3)

    # https://en.wikipedia.org/wiki/SRGB
    for wl = 1:3, y = 1:t.y, x = 1:t.x
        p = data[x, y, wl]
        if p <= 0.0031308
            data_new[y, x, wl] = 12.92 * p
        else
            data_new[y, x, wl] = (1+0.055) * p^(1/2.4) - 0.055
        end
        if isnan(data_new[y, x, wl])
            data_new[y, x, wl] = 0.0
        end
    end

    sm.imsave(filename, data_new)
end

function save_tally_hdf5(t, filename)
    h5open(filename, "w") do file
        write(file, "data", t.data)
        write(file, "x", t.x)
        write(file, "y", t.y)
        write(file, "scaling_factor", t.scaling_factor)
    end
end

function load_tally_hdf5(dir)
    t = []

    for fn in readdir(dir)
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
                    scaling_factor = read(file, "scaling_factor")
                end
                t = SvgRay.Tally(data, x, y, scaling_factor)
            else
                h5open(fn, "r") do file
                    t.data += read(file, "data")
                end
            end
        end
    end
    t
end