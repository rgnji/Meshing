import numpy as np

# projection of a quarter + outer region
def generate_block(center, r, R, theta_center, angle_start, angle_end, 
                   square_start, square_end, num_s, num_r, num_rl, num_arc, num_pipe, outer_z_start, outer_z_end, outer_theta_fixed):
    blocks = []

    #=================== thickness ===================
    # arc points
    angle = np.linspace(angle_start, angle_end, num_s+1)
    x = r * np.cos(angle)
    z = r * np.sin(angle)

    # square points
    square_x = np.linspace(square_start[0], square_end[0], num_s+1)
    square_z = np.linspace(square_start[1], square_end[1], num_s+1)

    # radial interpolation
    grid2d_x = np.linspace(x, square_x, num_r+1) + center[0]
    grid2d_z = np.linspace(z, square_z, num_r+1) + center[2]

    # project to cylinder
    theta = -np.arccos(grid2d_x / R) + theta_center

    # 3D mesh thickness
    R_grid = np.linspace(R, center[0], num_rl+1)
    jj, Rg, ii = np.meshgrid(np.arange(theta.shape[0]), R_grid, np.arange(theta.shape[1]))
    theta_grid = np.broadcast_to(theta, (num_rl+1, *theta.shape))

    grid_x = Rg * np.cos(theta_grid)
    grid_y = Rg * np.sin(theta_grid)
    grid_z = np.broadcast_to(grid2d_z, (num_rl+1, *grid2d_z.shape))

    # i,j,k to k,i,j
    grid_x = grid_x.transpose(1, 2, 0)
    grid_y = grid_y.transpose(1, 2, 0)
    grid_z = grid_z.transpose(1, 2, 0)

    # round
    grid_x = np.round(grid_x, 14)
    grid_y = np.round(grid_y, 14)
    grid_z = np.round(grid_z, 14)

    blocks.append((grid_x, grid_y, grid_z))

    #=================== another thickness ===================
    theta_grid_2 = theta_grid + np.pi
    grid_x_2 = Rg * np.cos(theta_grid_2)
    grid_y_2 = Rg * np.sin(theta_grid_2)

    # i,j,k to k,i,j
    grid_x_2 = grid_x_2.transpose(1, 2, 0)
    grid_y_2 = grid_y_2.transpose(1, 2, 0)

    # round
    grid_x_2 = np.round(grid_x_2, 14)
    grid_y_2 = np.round(grid_y_2, 14)

    blocks.append((grid_x_2, grid_y_2, grid_z))

    #=================== inlet pipe ===================
    # y interpolation
    grid2d_y = R * np.sin(-np.arccos(grid2d_x / R))
    inlet2d_y = np.full_like(grid2d_y, center[1])
    inlet_y = np.linspace(inlet2d_y, grid2d_y, num_pipe+1)
    inlet_x = np.broadcast_to(grid2d_x, (num_pipe+1, *grid2d_x.shape))
    inlet_z = np.broadcast_to(grid2d_z, (num_pipe+1, *grid2d_z.shape))

    # rotate
    theta_inlet = np.arctan(inlet_y / inlet_x)
    r_inlet = np.sqrt(inlet_y**2 + inlet_x**2)
    theta_inlet += theta_center
    inlet_x = r_inlet * np.cos(theta_inlet)
    inlet_y = r_inlet * np.sin(theta_inlet)

    # i,j,k to k,i,j
    inlet_x = inlet_x.transpose(1, 2, 0)
    inlet_y = inlet_y.transpose(1, 2, 0)
    inlet_z = inlet_z.transpose(1, 2, 0)

    # round
    inlet_x = np.round(inlet_x, 14)
    inlet_y = np.round(inlet_y, 14)
    inlet_z = np.round(inlet_z, 14)

    blocks.append((inlet_x, inlet_y, inlet_z))

    #=================== another pipe ===================
    theta_inlet_2 = theta_inlet + np.pi
    inlet_x_2 = r_inlet * np.cos(theta_inlet_2)
    inlet_y_2 = r_inlet * np.sin(theta_inlet_2)

    # i,j,k to k,i,j
    inlet_x_2 = inlet_x_2.transpose(1, 2, 0)
    inlet_y_2 = inlet_y_2.transpose(1, 2, 0)

    # round
    inlet_x_2 = np.round(inlet_x_2, 14)
    inlet_y_2 = np.round(inlet_y_2, 14)

    blocks.append((inlet_x_2, inlet_y_2, inlet_z))

    #=================== outer extension ===================
    bound_theta = np.full(num_s+1, outer_theta_fixed) # axial direction
    bound_z = np.linspace(outer_z_start, outer_z_end, num_s+1) # axial direction
    out2d_theta = np.linspace(bound_theta, theta[0], num_arc+1) # arc direction
    out2d_z = np.linspace(bound_z, grid2d_z[0], num_arc+1) # arc direction

    jj, Rg, ii = np.meshgrid(np.arange(out2d_theta.shape[0]), R_grid, np.arange(out2d_theta.shape[1]))
    out_theta = np.broadcast_to(out2d_theta, (num_rl+1, *out2d_theta.shape))

    out_x = Rg * np.cos(out_theta)
    out_y = Rg * np.sin(out_theta)
    out_z = np.broadcast_to(out2d_z, (num_rl+1, *out2d_z.shape))

    # i,j,k to k,i,j
    out_x = out_x.transpose(1, 2, 0)
    out_y = out_y.transpose(1, 2, 0)
    out_z = out_z.transpose(1, 2, 0)

    # round
    out_x = np.round(out_x, 14)
    out_y = np.round(out_y, 14)
    out_z = np.round(out_z, 14)

    blocks.append((out_x, out_y, out_z))

    #=================== another outer ===================
    out_theta_2 = out_theta + np.pi
    out_x_2 = Rg * np.cos(out_theta_2)
    out_y_2 = Rg * np.sin(out_theta_2)

    # i,j,k to k,i,j
    out_x_2 = out_x_2.transpose(1, 2, 0)
    out_y_2 = out_y_2.transpose(1, 2, 0)

    # round
    out_x_2 = np.round(out_x_2, 14)
    out_y_2 = np.round(out_y_2, 14)

    blocks.append((out_x_2, out_y_2, out_z))

    return blocks




def liquid_cylinder(r, R, R_inner, num_s, num_r, num_rl, num_arc, num_pipe, inlet_pipe_len):
    """
    Meshing for tangetial inlet.
    Parameters:
        r: inlet orifice radius
        R: liquid domain outer radius
        R_inner: liquid domain inner radius
        num_s: cells on orifice arc
        num_r: cells on orifice radius
        num_rl: cells on liquid domain thickness
        num_arc: cells on liquid domain arc
        num_pipe: inlet pipe cells
        inlet_pipe_len: distance from inlet face to y-axis
    """
    XT, YT, ZT = [], [], []

    center = np.array([R_inner, -inlet_pipe_len, -0.5]) # orifice center 1
    theta_center = np.arccos(center[0] / R) # angle after projection

    blocks_def = [
        (0, 0.5*np.pi, [0.5*r, 0], [0, 0.5*r], -0.5, 0, 0.5*np.pi),
        (0.5*np.pi, np.pi, [0, 0.5*r], [-0.5*r, 0], 0, -0.5, -0.5*np.pi),
        (np.pi, 1.5*np.pi, [-0.5*r, 0], [0, -0.5*r], -0.5, -1, -0.5*np.pi),
        (1.5*np.pi, 2*np.pi, [0, -0.5*r], [0.5*r, 0], -1, -0.5, 0.5*np.pi),
    ]

    #=================== four quarter ===================
    for angle_start, angle_end, square_start, square_end, outer_z_start, outer_z_end, outer_theta_fixed in blocks_def:
        blocks = generate_block(center, r, R, theta_center, 
                                 angle_start, angle_end, 
                                 square_start, square_end,
                                 num_s, num_r, num_rl, num_arc, num_pipe,
                                 outer_z_start, outer_z_end, outer_theta_fixed)
        for grid_x, grid_y, grid_z in blocks:
            XT.append(grid_x)
            YT.append(grid_y)
            ZT.append(grid_z)
    
    #=================== inner square thickness ===================
    square_start1 = blocks_def[0][2]
    square_end1 = blocks_def[0][3]
    square_start3 = blocks_def[2][2]
    square_end3 = blocks_def[2][3]

    square_up_x = np.linspace(square_start1[0], square_end1[0], num_s+1)
    square_up_z = np.linspace(square_start1[1], square_end1[1], num_s+1)
    square_down_x = np.linspace(square_start3[0], square_end3[0], num_s+1)
    square_down_z = np.linspace(square_start3[1], square_end3[1], num_s+1)

    square_x = np.linspace(square_up_x, square_down_x[::-1], num_s+1) + center[0]
    square_z = np.linspace(square_up_z, square_down_z[::-1], num_s+1) + center[2]

    # project
    square_theta = -np.arccos(square_x / R) + theta_center

    # 3D mesh
    R_grid = np.linspace(R, center[0], num_rl+1)
    jj, Rg, ii = np.meshgrid(np.arange(square_theta.shape[0]), R_grid, np.arange(square_theta.shape[1]))
    theta_grid = np.broadcast_to(square_theta, (num_rl+1, *square_theta.shape))

    grid_x = Rg * np.cos(theta_grid)
    grid_y = Rg * np.sin(theta_grid)
    grid_z = np.broadcast_to(square_z, (num_rl+1, *square_z.shape))

    # round
    grid_x = np.round(grid_x, 14)
    grid_y = np.round(grid_y, 14)
    grid_z = np.round(grid_z, 14)

    XT.append(grid_x)
    YT.append(grid_y)
    ZT.append(grid_z)

    #=================== another inner square thickness ===================
    theta_grid_2 = theta_grid + np.pi
    grid_x_2 = Rg * np.cos(theta_grid_2)
    grid_y_2 = Rg * np.sin(theta_grid_2)

    # round
    grid_x_2 = np.round(grid_x_2, 14)
    grid_y_2 = np.round(grid_y_2, 14)

    XT.append(grid_x_2)
    YT.append(grid_y_2)
    ZT.append(grid_z)

    #=================== pipe inner square ===================
    # y interpolation
    square_y = R * np.sin(-np.arccos(square_x / R))
    square_inlet_y = np.full_like(square_y, center[1])
    pipe_square_y = np.linspace(square_inlet_y, square_y, num_pipe+1)
    pipe_square_x = np.broadcast_to(square_x, (num_pipe+1, *square_x.shape))
    pipe_square_z = np.broadcast_to(square_z, (num_pipe+1, *square_z.shape))

    # rotate
    square_theta_inlet = np.arctan(pipe_square_y / pipe_square_x)
    square_r_inlet = np.sqrt(pipe_square_y**2 + pipe_square_x**2)
    square_theta_inlet += theta_center
    pipe_square_x = square_r_inlet * np.cos(square_theta_inlet)
    pipe_square_y = square_r_inlet * np.sin(square_theta_inlet)

    # round
    pipe_square_x = np.round(pipe_square_x, 14)
    pipe_square_y = np.round(pipe_square_y, 14)
    pipe_square_z = np.round(pipe_square_z, 14)

    XT.append(pipe_square_x)
    YT.append(pipe_square_y)
    ZT.append(pipe_square_z)

    #=================== another pipe inner square ===================
    square_theta_inlet_2 = square_theta_inlet + np.pi
    pipe_square_x_2 = square_r_inlet * np.cos(square_theta_inlet_2)
    pipe_square_y_2 = square_r_inlet * np.sin(square_theta_inlet_2)

    # round
    pipe_square_x_2 = np.round(pipe_square_x_2, 14)
    pipe_square_y_2 = np.round(pipe_square_y_2, 14)

    XT.append(pipe_square_x_2)
    YT.append(pipe_square_y_2)
    ZT.append(pipe_square_z)

    return XT, YT, ZT