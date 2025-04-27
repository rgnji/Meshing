import numpy as np

def liquid_cylinder():
    XT = []
    YT = []
    ZT = []

    num_s = 4
    num_r = 4
    num_R = 3
    num_S1 = 30 # up
    num_S2 = 20 # down

    # 1
    # arc
    angle = np.linspace(0, 0.5*np.pi, num_s)
    y = 0
    center = np.array([15.11, 0, -0.5])
    r = 0.5
    x = r * np.cos(angle) + center[0]
    z = r * np.sin(angle) + center[2]
    # square
    start = np.array([0.5*r, 0])
    end = np.array([0, 0.5*r])
    square_x = np.linspace(start[0], end[0], num_s) + center[0]
    square_z = np.linspace(start[1], end[1], num_s) + center[2]
    # radial
    grid2d_x = np.linspace(x, square_x, num_r)
    grid2d_z = np.linspace(z, square_z, num_r)

    # project
    R = 15.61
    fix = np.array([15.61, 0, -0.5])
    d = fix[0] - grid2d_x
    theta = -np.arccos((R - d) / R) # R, theta, z

    # 3D mesh
    R_grid = np.linspace(15.61, 15.11, num_R)
    jj, R_grid, ii = np.meshgrid(np.arange(np.shape(theta)[0]),
                                 R_grid,
                                 np.arange(np.shape(theta)[1]))
    theta_grid = np.full((num_R, np.shape(theta)[0], np.shape(theta)[1]), theta)
    grid_x = R_grid * np.cos(theta_grid)
    grid_y = R_grid * np.sin(theta_grid)
    grid_z = np.full((num_R, np.shape(grid2d_z)[0], np.shape(grid2d_z)[1]), grid2d_z)

    # update
    XT.append(grid_x)
    YT.append(grid_y)
    ZT.append(grid_z)

    # outer
    bound_theta = np.full(num_s, 0.5*np.pi)
    bound_z = np.linspace(-0.5, 0, num_s)
    out2d_theta = np.linspace(bound_theta, theta[0], num_S1)
    out2d_z = np.linspace(bound_z, grid2d_z[0], num_S1)

    # 3D mesh
    R_grid = np.linspace(15.61, 15.11, num_R)
    jj, R_grid, ii = np.meshgrid(np.arange(np.shape(out2d_theta)[0]),
                                 R_grid,
                                 np.arange(np.shape(out2d_theta)[1]))
    out_theta = np.full((num_R, np.shape(out2d_theta)[0], np.shape(out2d_theta)[1]), out2d_theta)
    out_x = R_grid * np.cos(out_theta)
    out_y = R_grid * np.sin(out_theta)
    out_z = np.full((num_R, np.shape(out2d_z)[0], np.shape(out2d_z)[1]), out2d_z)

    # update
    XT.append(out_x)
    YT.append(out_y)
    ZT.append(out_z)


    # 2
    # arc
    angle = np.linspace(0.5*np.pi, np.pi, num_s)
    y = 0
    center = np.array([15.11, 0, -0.5])
    r = 0.5
    x = r * np.cos(angle) + center[0]
    z = r * np.sin(angle) + center[2]
    # square
    start = np.array([0, 0.5*r])
    end = np.array([-0.5*r, 0])
    square_x = np.linspace(start[0], end[0], num_s) + center[0]
    square_z = np.linspace(start[1], end[1], num_s) + center[2]
    # radial
    grid2d_x = np.linspace(x, square_x, num_r)
    grid2d_z = np.linspace(z, square_z, num_r)

    # project
    R = 15.61
    fix = np.array([15.61, 0, -0.5])
    d = fix[0] - grid2d_x
    theta = -np.arccos((R - d) / R) # R, theta, z

    # 3D mesh
    R_grid = np.linspace(15.61, 15.11, num_R)
    jj, R_grid, ii = np.meshgrid(np.arange(np.shape(theta)[0]),
                                 R_grid,
                                 np.arange(np.shape(theta)[1]))
    theta_grid = np.full((num_R, np.shape(theta)[0], np.shape(theta)[1]), theta)
    grid_x = R_grid * np.cos(theta_grid)
    grid_y = R_grid * np.sin(theta_grid)
    grid_z = np.full((num_R, np.shape(grid2d_z)[0], np.shape(grid2d_z)[1]), grid2d_z)

    # update
    XT.append(grid_x)
    YT.append(grid_y)
    ZT.append(grid_z)

    # outer
    bound_theta = np.full(num_s, -0.5*np.pi)
    bound_z = np.linspace(0, -0.5, num_s)
    out2d_theta = np.linspace(bound_theta, theta[0], num_S1)
    out2d_z = np.linspace(bound_z, grid2d_z[0], num_S1)

    # 3D mesh
    R_grid = np.linspace(15.61, 15.11, num_R)
    jj, R_grid, ii = np.meshgrid(np.arange(np.shape(out2d_theta)[0]),
                                 R_grid,
                                 np.arange(np.shape(out2d_theta)[1]))
    out_theta = np.full((num_R, np.shape(out2d_theta)[0], np.shape(out2d_theta)[1]), out2d_theta)
    out_x = R_grid * np.cos(out_theta)
    out_y = R_grid * np.sin(out_theta)
    out_z = np.full((num_R, np.shape(out2d_z)[0], np.shape(out2d_z)[1]), out2d_z)

    # update
    XT.append(out_x)
    YT.append(out_y)
    ZT.append(out_z)


    # 3
    # arc
    angle = np.linspace(np.pi, 1.5*np.pi, num_s)
    y = 0
    center = np.array([15.11, 0, -0.5])
    r = 0.5
    x = r * np.cos(angle) + center[0]
    z = r * np.sin(angle) + center[2]
    # square
    start = np.array([-0.5*r, 0])
    end = np.array([0, -0.5*r])
    square_x = np.linspace(start[0], end[0], num_s) + center[0]
    square_z = np.linspace(start[1], end[1], num_s) + center[2]
    # radial
    grid2d_x = np.linspace(x, square_x, num_r)
    grid2d_z = np.linspace(z, square_z, num_r)

    # project
    R = 15.61
    fix = np.array([15.61, 0, -0.5])
    d = fix[0] - grid2d_x
    theta = -np.arccos((R - d) / R) # R, theta, z

    # 3D mesh
    R_grid = np.linspace(15.61, 15.11, num_R)
    jj, R_grid, ii = np.meshgrid(np.arange(np.shape(theta)[0]),
                                 R_grid,
                                 np.arange(np.shape(theta)[1]))
    theta_grid = np.full((num_R, np.shape(theta)[0], np.shape(theta)[1]), theta)
    grid_x = R_grid * np.cos(theta_grid)
    grid_y = R_grid * np.sin(theta_grid)
    grid_z = np.full((num_R, np.shape(grid2d_z)[0], np.shape(grid2d_z)[1]), grid2d_z)

    # update
    XT.append(grid_x)
    YT.append(grid_y)
    ZT.append(grid_z)

    # outer
    bound_theta = np.full(num_s, -0.5*np.pi)
    bound_z = np.linspace(-0.5, -1, num_s)
    out2d_theta = np.linspace(bound_theta, theta[0], num_S1)
    out2d_z = np.linspace(bound_z, grid2d_z[0], num_S1)

    # 3D mesh
    R_grid = np.linspace(15.61, 15.11, num_R)
    jj, R_grid, ii = np.meshgrid(np.arange(np.shape(out2d_theta)[0]),
                                 R_grid,
                                 np.arange(np.shape(out2d_theta)[1]))
    out_theta = np.full((num_R, np.shape(out2d_theta)[0], np.shape(out2d_theta)[1]), out2d_theta)
    out_x = R_grid * np.cos(out_theta)
    out_y = R_grid * np.sin(out_theta)
    out_z = np.full((num_R, np.shape(out2d_z)[0], np.shape(out2d_z)[1]), out2d_z)

    # update
    XT.append(out_x)
    YT.append(out_y)
    ZT.append(out_z)


    # 4
    # arc
    angle = np.linspace(1.5*np.pi, 2*np.pi, num_s)
    y = 0
    center = np.array([15.11, 0, -0.5])
    r = 0.5
    x = r * np.cos(angle) + center[0]
    z = r * np.sin(angle) + center[2]
    # square
    start = np.array([0, -0.5*r])
    end = np.array([0.5*r, 0])
    square_x = np.linspace(start[0], end[0], num_s) + center[0]
    square_z = np.linspace(start[1], end[1], num_s) + center[2]
    # radial
    grid2d_x = np.linspace(x, square_x, num_r)
    grid2d_z = np.linspace(z, square_z, num_r)

    # project
    R = 15.61
    fix = np.array([15.61, 0, -0.5])
    d = fix[0] - grid2d_x
    theta = -np.arccos((R - d) / R) # R, theta, z

    # 3D mesh
    R_grid = np.linspace(15.61, 15.11, num_R)
    jj, R_grid, ii = np.meshgrid(np.arange(np.shape(theta)[0]),
                                 R_grid,
                                 np.arange(np.shape(theta)[1]))
    theta_grid = np.full((num_R, np.shape(theta)[0], np.shape(theta)[1]), theta)
    grid_x = R_grid * np.cos(theta_grid)
    grid_y = R_grid * np.sin(theta_grid)
    grid_z = np.full((num_R, np.shape(grid2d_z)[0], np.shape(grid2d_z)[1]), grid2d_z)

    # update
    XT.append(grid_x)
    YT.append(grid_y)
    ZT.append(grid_z)

    # outer
    bound_theta = np.full(num_s, 0.5*np.pi)
    bound_z = np.linspace(-1, -0.5, num_s)
    out2d_theta = np.linspace(bound_theta, theta[0], num_S1)
    out2d_z = np.linspace(bound_z, grid2d_z[0], num_S1)

    # 3D mesh
    R_grid = np.linspace(15.61, 15.11, num_R)
    jj, R_grid, ii = np.meshgrid(np.arange(np.shape(out2d_theta)[0]),
                                 R_grid,
                                 np.arange(np.shape(out2d_theta)[1]))
    out_theta = np.full((num_R, np.shape(out2d_theta)[0], np.shape(out2d_theta)[1]), out2d_theta)
    out_x = R_grid * np.cos(out_theta)
    out_y = R_grid * np.sin(out_theta)
    out_z = np.full((num_R, np.shape(out2d_z)[0], np.shape(out2d_z)[1]), out2d_z)

    # update
    XT.append(out_x)
    YT.append(out_y)
    ZT.append(out_z)

    return XT, YT, ZT
