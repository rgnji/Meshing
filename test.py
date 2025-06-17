import numpy as np
"""
x=np.linspace(1, 10, 10)
y=np.linspace(2, 11, 10)
z=np.linspace(x,y,10)
z[0] = x
print(z)"""

import matplotlib.pyplot as plt
from numpy import pi, sin, cos, array
from lib.plot3dout import binary_fort12
from lib.laplace import solve_grid_laplace_sor

def o_grid(qua, r, r_square, arc_angle, axi_pos, axial, d1, d2, d3, expr=1):
    """
    Generate o-grid (no inner h-grid).
    Parameters:
        qua: 
            quadrant
        r: 
            radius of cylinder
        r_square: 
            half diagonal of inner square or radius of inner cylinder
        arc_angle: 
            angle of inner arc in degree
        axi_pos: 
            position in axial direction
        axial: 
            axial length
        d1: 
            number of cells in arc direction
        d2: 
            number of cells in radius direction
        d3: 
            number of cells in axial direction 
        tp:
            1 for log scale distribution, 0 for linear distribution
    Returns:
        X, Y, Z:
            d3+1 * d2+1 * d1+1
    """

    # find radius
    rot_qua = (qua - 1) * 0.5 * pi
    rotation = array([ [cos(rot_qua), sin(rot_qua)], [-sin(rot_qua), cos(rot_qua)] ])
    start = array([r_square, 0]) @ rotation
    end = array([0, r_square]) @ rotation
    angle = np.deg2rad(arc_angle)
    radius = np.linalg.norm(end - start) / sin(angle) * sin(0.5 * pi - 0.5 * angle)

    # find center
    dir_vec = end - start # start -> end
    dir_vec = dir_vec / np.linalg.norm(dir_vec)
    perp_vec = array([ -dir_vec[1], dir_vec[0] ]) # counter-clockwise
    midpoint = (end + start) / 2
    center = midpoint + perp_vec * (0.5 * np.linalg.norm(end - start) / sin(0.5 * angle) * sin(0.5 * (pi - angle)))

    # create boundary
    arc_div = np.linspace(0.25*pi + rot_qua - 0.5*angle, 0.25*pi + rot_qua + 0.5*angle, d1+1)
    arc_x = radius * cos(arc_div) + center[0]
    arc_y = radius * sin(arc_div) + center[1]
    arc_x[0] = start[0]
    arc_x[-1] = end[0]
    arc_y[0] = start[1]
    arc_y[-1] = end[1]
    
    # create outer boundary
    out_div = np.linspace(rot_qua, rot_qua + 0.5*pi, d1+1)
    out_x = r * cos(out_div)
    out_y = r * sin(out_div)
    
    # initial grid between boundaries
    grid_x = np.zeros((d2+1, d1+1))
    grid_x[0] = out_x
    grid_y = np.zeros((d2+1, d1+1))
    grid_y[0] = out_y
    
    Lx = arc_x - out_x
    Ly = arc_y - out_y
    ratios = 1 + (expr - 1) * np.linspace(0, 1, d2)
    total_ratio = np.sum(ratios)
    for i in range(d1+1):
        segmentx = Lx[i] * ratios / total_ratio
        segmenty = Ly[i] * ratios / total_ratio
        for j in range(1, d2+1):
            grid_x[j, i] = grid_x[j-1, i] + segmentx[j-1]
            grid_y[j, i] = grid_y[j-1, i] + segmenty[j-1]
    grid_x[-1] = arc_x
    grid_y[-1] = arc_y
    
    """
    grid_x = np.linspace(out_x, arc_x, d2+1)
    grid_y = np.linspace(out_y, arc_y, d2+1)
    """
    
    # solve laplace
    grid_x = np.round(grid_x, 14)
    grid_y = np.round(grid_y, 14)
    if arc_angle != 90:
        x, y, iterations = solve_grid_laplace_sor(grid_x, grid_y)
        print(f"SOR of O-grid in quadrant {qua} converged in {iterations} iterations.")
    else:
        x, y = grid_x, grid_y

    # create 3D mesh
    z = np.linspace(axi_pos - axial, axi_pos, d3+1)
    X = np.full( (d3+1, x.shape[0], x.shape[1]), x)  # shape(0, 1, 2) -> k, j, i 
    Y = np.full( (d3+1, y.shape[0], y.shape[1]), y)
    jj, Z, ii = np.meshgrid(np.arange(x.shape[0]),   # j, k, i
                          z,
                          np.arange(x.shape[1]))

    return X, Y, Z

# gas inlet
RGSOUT = 14.11
RGSHGD = 0.5*RGSOUT
ZGSSTR = 0
LENGS  = 6
DARCGS = 10
DRADGS = 7
DAXLGS = 7

mesh = ['o', RGSOUT, RGSHGD, 90, ZGSSTR, LENGS,  DARCGS, DRADGS, DAXLGS, 5]
XT = []
YT = []
ZT = []
x, y, z = o_grid(1, *mesh[1:])
XT.append(x)
YT.append(y)
ZT.append(z)

plt.scatter(x[0], y[0])
plt.show()