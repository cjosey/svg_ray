import Base.contains

# =============================================================================
# Types
# =============================================================================

"""
Cell type.  Contains all the surfaces as well as a description of the material
within the cell.
"""
type Cell
    lines::LineArray   # Surface linear segments
    q_bez::QBezier     # Surface quadratic bezier curves
    c_bez::CBezier     # Surface cubic bezier curves
    material::Material # The material, of type Material
end

# =============================================================================
# Functions
# =============================================================================

"""
Check if x-y coordinate pair ``point`` exists inside of the cell, where it 
first collides, and the normal vector at said point.
"""
function contains(self::Cell, point::Array, point2=[-1,-1])
    int_count1, new_point1, min_length1, normal1 = intersect(self.lines, point, point2)
    int_count2, new_point2, min_length2, normal2 = intersect(self.c_bez, point, point2)
    int_count3, new_point3, min_length3, normal3 = intersect(self.q_bez, point, point2)

    int_count = int_count1 + int_count2 + int_count3

    inside = int_count % 2

    min_length = min(min_length1, min_length2, min_length3)

    # Find which possible surface it collides with
    if min_length == min_length1
        new_point = new_point1
        normal = normal1
    elseif min_length == min_length2
        new_point = new_point2
        normal = normal2
    else
        new_point = new_point3
        normal = normal3
    end


    inside, new_point, min_length, normal
end

"""
Plots cell using PyPlot.  Assumes somewhere else handles making and saving the figure.
"""
function plot_cell(cell, ax)
    @pyimport matplotlib.patches as patches
    @pyimport matplotlib.path as path
    # lines
    for i = 1:size(cell.lines.points)[1]
        x1 = cell.lines.points[i,1]
        y1 = cell.lines.points[i,2]
        x2 = cell.lines.points[i,3]
        y2 = cell.lines.points[i,4]
        ax[:plot]([x1, x2], [y1, y2], "k", linewidth=0.25)
    end

    # cubic bezier
    for i = 1:size(cell.c_bez.points)[2]
        x1 = cell.c_bez.points[1,i]
        y1 = cell.c_bez.points[2,i]
        x2 = cell.c_bez.points[3,i]
        y2 = cell.c_bez.points[4,i]
        x3 = cell.c_bez.points[5,i]
        y3 = cell.c_bez.points[6,i]
        x4 = cell.c_bez.points[7,i]
        y4 = cell.c_bez.points[8,i]

        verts = [(x1, y1), (x2, y2), (x3, y3), (x4, y4)]
        codes = [path.Path[:MOVETO], path.Path[:CURVE4], path.Path[:CURVE4], path.Path[:CURVE4]]

        pat = path.Path(verts, codes)
        pa = patches.PathPatch(pat, facecolor="none", linewidth=0.25)

        ax[:add_patch](pa)
    end

    # quad bezier
    for i = 1:size(cell.q_bez.points)[2]
        x1 = cell.q_bez.points[1,i]
        y1 = cell.q_bez.points[2,i]
        x2 = cell.q_bez.points[3,i]
        y2 = cell.q_bez.points[4,i]
        x3 = cell.q_bez.points[5,i]
        y3 = cell.q_bez.points[6,i]

        verts = [(x1, y1), (x2, y2), (x3, y3)]
        codes = [path.Path[:MOVETO], path.Path[:CURVE3], path.Path[:CURVE3]]

        pat = path.Path(verts, codes)
        pa = patches.PathPatch(pat, facecolor="none", linewidth=0.25)

        ax[:add_patch](pa)
    end
end