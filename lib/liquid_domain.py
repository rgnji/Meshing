import numpy as np

# projection of a quarter + outer region
def generate_block(center, r, R, fix, theta_center, angle_start, angle_end, 
                   square_start, square_end, num_s, num_r, num_rl, num_arc, outer_z_start, outer_z_end, outer_theta_fixed):
    blocks = []

    # arc points
    angle = np.linspace(angle_start, angle_end, num_s+1)
    x = r * np.cos(angle) + center[0]
    z = r * np.sin(angle) + center[2]

    # square points
    square_x = np.linspace(square_start[0], square_end[0], num_s+1) + center[0]
    square_z = np.linspace(square_start[1], square_end[1], num_s+1) + center[2]

    # radial interpolation
    grid2d_x = np.linspace(x, square_x, num_r+1)
    grid2d_z = np.linspace(z, square_z, num_r+1)

    # project to cylinder
    d = fix[0] - grid2d_x
    theta = -np.arccos((R - d) / R) + theta_center

    # 3D mesh internal
    R_grid = np.linspace(R, center[0], num_rl+1)
    jj, Rg, ii = np.meshgrid(np.arange(theta.shape[0]), R_grid, np.arange(theta.shape[1]))
    theta_grid = np.broadcast_to(theta, (num_rl+1, *theta.shape))

    grid_x = Rg * np.cos(theta_grid)
    grid_y = Rg * np.sin(theta_grid)
    grid_z = np.broadcast_to(grid2d_z, (num_rl+1, *grid2d_z.shape))
    blocks.append((grid_x, grid_y, grid_z))

    # outer extension
    bound_theta = np.full(num_s+1, outer_theta_fixed)
    bound_z = np.linspace(outer_z_start, outer_z_end, num_s+1)
    out2d_theta = np.linspace(bound_theta, theta[0], num_arc)
    out2d_z = np.linspace(bound_z, grid2d_z[0], num_arc)

    jj, Rg, ii = np.meshgrid(np.arange(out2d_theta.shape[0]), R_grid, np.arange(out2d_theta.shape[1]))
    out_theta = np.broadcast_to(out2d_theta, (num_rl+1, *out2d_theta.shape))

    out_x = Rg * np.cos(out_theta)
    out_y = Rg * np.sin(out_theta)
    out_z = np.broadcast_to(out2d_z, (num_rl+1, *out2d_z.shape))
    blocks.append((out_x, out_y, out_z))

    return blocks


def liquid_cylinder():
    XT, YT, ZT = [], [], []

    # Parameters
    r = 0.5
    R = 15.61
    num_s = 3
    num_r = 3
    num_rl = 2
    num_arc = 10
    center = np.array([15.11, 0, -0.5])
    fix = np.array([15.61, 0, -0.5])
    d_center = fix[0] - center[0]
    theta_center = np.arccos((R - d_center) / R)

    blocks_def = [
        (0, 0.5*np.pi, [0.5*r, 0], [0, 0.5*r], -0.5, 0, 0.5*np.pi),
        (0.5*np.pi, np.pi, [0, 0.5*r], [-0.5*r, 0], 0, -0.5, -0.5*np.pi),
        (np.pi, 1.5*np.pi, [-0.5*r, 0], [0, -0.5*r], -0.5, -1, -0.5*np.pi),
        (1.5*np.pi, 2*np.pi, [0, -0.5*r], [0.5*r, 0], -1, -0.5, 0.5*np.pi),
    ]

    for angle_start, angle_end, square_start, square_end, outer_z_start, outer_z_end, outer_theta_fixed in blocks_def:
        blocks = generate_block(center, r, R, fix, theta_center, 
                                 angle_start, angle_end, 
                                 square_start, square_end,
                                 num_s, num_r, num_rl, num_arc,
                                 outer_z_start, outer_z_end, outer_theta_fixed)
        for grid_x, grid_y, grid_z in blocks:
            XT.append(grid_x)
            YT.append(grid_y)
            ZT.append(grid_z)

    return XT, YT, ZT