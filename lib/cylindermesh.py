import numpy as np
from numpy import array, cos, sin, pi
import sys

from lib.laplace import solve_grid_laplace_sor

#================= =================
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

    # input check
    if qua not in [1,2,3,4]:
        sys.exit("InputType Error: quadrant must be in 1,2,3,4.")
    if r <= 0 or r_square <= 0:
        sys.exit("InputType Error: radius must greater than 0.")
    if r < r_square:
        sys.exit("Input Error: the block is inside-out when r < r_square.")
    if arc_angle < 0:
        sys.exit("InputType Error: arc_angle must be positive or 0.")
    if type(d1) != int or type(d2) != int or type(d3) != int:
        sys.exit("InputType Error: number of cells must be integer.")
    if d1<=0 or d2<=0 or d3<=0:
        sys.exit("Input Error: number of cells must be positive integer.")

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
        x, y = grid_x, grid_y
        #x, y, iterations = solve_grid_laplace_sor(grid_x, grid_y)
        #print(f"SOR of O-grid in quadrant {qua} converged in {iterations} iterations.")
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


#================= =================
def h_grid(r_square, arc_angle, axi_pos, axial, d1, d3):
    """
    Generate h-grid (square).
    Parameters:
        r_square: 
            half diagonal of inner square or radius of inner cylinder
        arc_angle: 
            angle of arc in degree
        axi_pos: 
            position in axial direction
        axial: 
            axial length
        d1: 
            number of cells in arc direction
        d3: 
            number of cells in axial direction 
    Returns:
        X, Y, Z:
            d3+1 * d2+1 * d1+1
    """

    # input check
    if r_square <= 0:
        sys.exit("InputType Error: radius must greater than 0.")
    if arc_angle < 0:
        sys.exit("InputType Error: arc_angle must be positive or 0.")
    if type(d1) != int or type(d3) != int:
        sys.exit("InputType Error: number of cells must be integer.")
    if d1<=0 or d3<=0:
        sys.exit("Input Error: number of cells must be positive integer.")

    # boudary mesh
    arc_x = np.zeros((4, d1+1), float)
    arc_y = np.zeros((4, d1+1), float)
    for qua in range(1,5):
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
        arc_x[qua-1] = radius * cos(arc_div) + center[0]
        arc_y[qua-1] = radius * sin(arc_div) + center[1]
        arc_x[qua-1][0] = start[0]
        arc_x[qua-1][-1] = end[0]
        arc_y[qua-1][0] = start[1]
        arc_y[qua-1][-1] = end[1]
    
    # inner mesh
    grid_x = np.linspace(arc_x[0][1:-1], arc_x[2][1:-1][::-1], d1+1) # d1+1 * d1-1
    grid_y = np.linspace(arc_y[0][1:-1], arc_y[2][1:-1][::-1], d1+1)

    # combine mesh
    x = np.zeros((d1+1, d1+1), float)
    y = np.zeros((d1+1, d1+1), float)
    # 1
    x[0] = arc_x[0]
    y[0] = arc_y[0]
    # 3
    x[-1] = arc_x[2][::-1]
    y[-1] = arc_y[2][::-1]
    # 2
    x[:, -1] = arc_x[1]
    y[:, -1] = arc_y[1]
    # 4
    x[:, 0] = arc_x[3][::-1]
    y[:, 0] = arc_y[3][::-1]
    # inner
    x[1:-1, 1:-1] = grid_x[1:-1]
    y[1:-1, 1:-1] = grid_y[1:-1]

    # solve laplace
    x = np.round(x, 14)
    y = np.round(y, 14)
    x, y, iterations = solve_grid_laplace_sor(x, y)
    print(f"SOR of H-grid converged in {iterations} iterations.")

    # create 3D mesh
    z = np.linspace(axi_pos - axial, axi_pos, d3+1)
    X = np.full( (d3+1, x.shape[0], x.shape[1]), x)  # shape(0, 1, 2) -> k, j, i 
    Y = np.full( (d3+1, y.shape[0], y.shape[1]), y)
    jj, Z, ii = np.meshgrid(np.arange(x.shape[0]),   # j, k, i
                          z,
                          np.arange(x.shape[1]))

    return X, Y, Z