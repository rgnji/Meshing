import numpy as np
import matplotlib.pyplot as plt

from lib.cylindermesh import o_grid, h_grid
from lib.liquid_domain import liquid_cylinder
from lib.plot3dout import unformatted_fort12, binary_fort12

# arc_number == num_arc

#========== parameters for cylinder ==========


#========== parameters for liquid domain ==========
r_orifice = 0.5 # inlet orifice radius
R_liquidout = 15.61 # liquid domain outer radius
R_liquidin = 15.11  # liquid domain inner radius
num_s = 3    # cells on orifice arc
num_r = 3    # cells on orifice radius
num_rl = 2   # cells on liquid domain thickness
num_arc = 15 # cells on liquid domain arc
num_pipe = 20 # inlet pipe cells
inlet_pipe_len = 15.11 # inlet pipe length


#================= grid generation =================
XT = []
YT = []
ZT = []

# meshing parameters
mesh = [['o', 14.11, 0.5*14.11, 20, 0, 6, 15, 28, 10], # gas inlet
        ['h', 0.5*14.11, 20, 0, 6, 15, 10], # gas inlet
        ['o', 15.61, 15.11, 90, -1, 5, 15, 2, 10], # liquid cylinder
        ['o', 14.11, 0.5*14.11, 20, -6, 28.22, 15, 28, 50], # gas recess
        ['h', 0.5*14.11, 20, -6, 28.22, 15, 50], # gas recess
        ['o', 15.11, 14.11, 90, -6, 28.22, 15, 4, 50], # wall recess
        ['o', 15.61, 15.11, 90, -6, 28.22, 15, 2, 50], # liquid recess
        ['o', 14.11, 0.5*14.11, 20, -34.22, 20, 15, 28, 85], # gas out
        ['h', 0.5*14.11, 20, -34.22, 20, 15, 85], # gas out
        ['o', 15.11, 14.11, 90, -34.22, 20, 15, 4, 85], # wall out
        ['o', 15.61, 15.11, 90, -34.22, 20, 15, 2, 85], # liquid out
        ['o', 25, 15.61, 90, -34.22, 20, 15, 40, 85], # out
        ['l', r_orifice, R_liquidout, R_liquidin, num_s, num_r, num_rl, num_arc, num_pipe, inlet_pipe_len] # liquid inlet
        ]

# generate mesh
for grid in mesh:
    if grid[0] == 'o':
        for qua in range(1,5): 
            x, y, z = o_grid(qua, *grid[1:])
            XT.append(x)
            YT.append(y)
            ZT.append(z)
    
    if grid[0] == 'h':
        x, y, z = h_grid(*grid[1:])
        XT.append(x)
        YT.append(y)
        ZT.append(z)

    if grid[0] == 'l':
        xx, yy, zz = liquid_cylinder(*grid[1:])
        XT += xx
        YT += yy
        ZT += zz

# establish fort12
txt = binary_fort12(XT, YT, ZT)
print(txt)

txt = unformatted_fort12(XT, YT, ZT)
print(txt)