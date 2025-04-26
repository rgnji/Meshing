import numpy as np
import sys
import math

#========== create initial mesh of o-grid ==========
def initial_o(quadrant, ro, ri, arc_number, radius_number):
    """
    Input:
        ro: radius of cylinder
        ri: radius of the arc of inner square
        arc_number: number of division on arc direction
        radius_number: number of division on radius direction
        quadrant: one of the four segment of o-grid
    Output:
        x: meshgrid of x coordinates
        y: meshgrid of y coordinates
    """
    # type check
    if type(arc_number) != int or type(radius_number) != int or quadrant not in [1,2,3,4]:
        sys.exit("arc_number, radius_number and quadrant must be integer.\nprogram stopped.")

    # inner block boundary
    inner_square_angle = 20  # inner square angle (degree)
    theta = inner_square_angle * np.pi / 180 # radian
    phi = 0.5 * (0.5 * np.pi - theta) # radian   

    start = phi + 0.5 * np.pi * (quadrant-1)
    end   = phi + 0.5 * np.pi * (quadrant-1) + theta
    inner_range = np.linspace(start, end, arc_number+1)

    inner_x = np.cos(inner_range) * ri
    inner_y = np.sin(inner_range) * ri
    if quadrant % 2 == 1: # 1 or 3
        inner_x -= inner_x[-1]
        inner_y -= inner_y[0]
    elif quadrant % 2 == 0: # 2 or 4
        inner_x -= inner_x[0]
        inner_y -= inner_y[-1]

    # outer block boundary
    outer_angle = np.linspace(0.5 * np.pi * (quadrant-1), 0.5 * np.pi * quadrant, arc_number+1)
    outer_x = np.cos(outer_angle) * ro
    outer_y = np.sin(outer_angle) * ro

    # initial grid between boundaries
    ratio = np.array(range(1,radius_number)) / radius_number
    inner_grid_x, inner_ratio_x = np.meshgrid(inner_x[:], ratio)
    inner_grid_y, inner_ratio_y = np.meshgrid(inner_y[:], ratio)
    outer_grid_x, outer_ratio_x = np.meshgrid(outer_x[:], 1 - ratio)
    outer_grid_y, outer_ratio_y = np.meshgrid(outer_y[:], 1 - ratio)

    grid_x = inner_grid_x * inner_ratio_x + outer_grid_x * outer_ratio_x
    grid_y = inner_grid_y * inner_ratio_y + outer_grid_y * outer_ratio_y

    # initial grid
    x = np.concatenate(([outer_x], grid_x, [inner_x]), axis=0)
    y = np.concatenate(([outer_y], grid_y, [inner_y]), axis=0)

    return x, y


#========== create initial mesh of inner square ==========
"""

  start  /\ start
        /  \ 
       /    \ 
      /      \ 
 end /        \ end
     \        /
      \      /
       \    /
        \  /
         \/
"""

def initial_square(ri, arc_number):
    """
    Input:
        ri: radius of the arc of inner square
        arc_number: number of division on arc direction
    Output:
        x: meshgrid of x coordinates
        y: meshgrid of y coordinates
    """

    x = np.zeros((arc_number+1, arc_number+1), float)
    y = np.zeros((arc_number+1, arc_number+1), float)
    boundary = np.zeros((4,2), object)

    # create boundary mesh
    for i in range(1,5):
        inner_square_angle = 20  # inner square angle (degree)
        theta = inner_square_angle * np.pi / 180 # radian
        phi = 0.5 * (0.5 * np.pi - theta) # radian   

        start = phi + 0.5 * np.pi * (i-1)
        end   = phi + 0.5 * np.pi * (i-1) + theta
        inner_range = np.linspace(start, end, arc_number+1)

        boundary[i-1][0] = np.cos(inner_range) * ri # x
        boundary[i-1][1] = np.sin(inner_range) * ri # y

        if i % 2 == 1: # 1 or 3
            boundary[i-1][0] -= boundary[i-1][0][-1]
            boundary[i-1][1] -= boundary[i-1][1][0]
        elif i % 2 == 0: # 2 or 4
            boundary[i-1][0] -= boundary[i-1][0][0]
            boundary[i-1][1] -= boundary[i-1][1][-1]

    # reverse quadrant 1 and 4  
    boundary[0][0] = boundary[0][0][::-1]
    boundary[0][1] = boundary[0][1][::-1]

    boundary[3][0] = boundary[3][0][::-1]
    boundary[3][1] = boundary[3][1][::-1]

    # initial grid between boundaries
    ratio = np.array(range(1,arc_number)) / arc_number
    inner_grid_x, inner_ratio_x = np.meshgrid(boundary[0,0][1:-1], 1 - ratio)     # quadrant 1 -> up
    inner_grid_y, inner_ratio_y = np.meshgrid(boundary[0,1][1:-1], 1 - ratio)
    outer_grid_x, outer_ratio_x = np.meshgrid(boundary[2,0][1:-1], ratio) # quadrant 3 -> down
    outer_grid_y, outer_ratio_y = np.meshgrid(boundary[2,1][1:-1], ratio)
    grid_x = inner_grid_x * inner_ratio_x + outer_grid_x * outer_ratio_x
    grid_y = inner_grid_y * inner_ratio_y + outer_grid_y * outer_ratio_y

    # combine boudaries and inner grid
    x[0] = boundary[0, 0]     # quadrant 1
    y[0] = boundary[0, 1]   
     
    x[-1] = boundary[2, 0]    # quadrant 3
    y[-1] = boundary[2, 1]   
    
    x[:, 0] = boundary[1, 0]  # quadrant 2
    y[:, 0] = boundary[1, 1]    

    x[:, -1] = boundary[3, 0] # quadrant 4
    y[:, -1] = boundary[3, 1]
    
    x[1:-1, 1:-1] = grid_x    # inner grid
    y[1:-1, 1:-1] = grid_y
    
    return x, y


# parameter: quadrant, radius of cylinder, vertices of inner square, arc angle, number of division in three direction

def o_grid(qua, r, end, start, arc_angle, axial, d1, d2, d3):
    """
    Input:
        qua: quadrant
        r: radius of cylinder
        r_inner: ratio of half diagonal of inner square to radius of cylinder, between 0 and 1
        start: [x,y]
        end: [x,y]
        arc_angle: angle of inner arc
        axial: axial length
        d1: number of division in arc direction
        d2: number of division in radius direction
        d3: number of division in axial direction 
    """

    # find radius
    angle = np.deg2rad(arc_angle)
    radius = np.linalg.norm(end - start) / np.sin(angle) * np.sin(0.5 * np.pi - 0.5 * angle)

    dir_vec = end - start # start -> end
    dir_vec = dir_vec / np.linalg.norm(dir_vec)

    perp_vec = np.array([ -dir_vec[1], dir_vec[0] ]) # counter-clockwise
    midpoint = (end + start) / 2
    center = midpoint + perp_vec * (0.5 * np.linalg.norm(dir_vec) / np.sin(0.5 * angle) * np.sin(0.5 * (np.pi - angle)))

    vec_start = (start - center) / np.linalg.norm(start - center)
    vec_end = (end - center) / np.linalg.norm(end - center)

    perp_start = np.array([ -vec_start[1], vec_start[0] ])
    perp_end = np.array([ vec_end[1], -vec_end[0]])

    n = 2 * np.pi / angle # segment
    distance = 4/3 * radius *np.tan(np.pi / (2  * n))
    
    control_start = start + perp_start * distance
    control_end = end + perp_end * distance

    
