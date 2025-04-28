import numpy as np
import matplotlib.pyplot as plt

from lib.cylindermesh import o_grid, h_grid
from lib.plot3dout import unformatted_fort12, binary_fort12
from lib.liquid_domain import liquid_cylinder

# arc_number == num_arc

#========== parameters for cylinder ==========
r = 5 # cylinder radius
r_square = 0.5 * r # diagonal of inner square (or radius of inner cylinder)
len_z = 6 # axial length of cylinder
pos_z = 0 # axial position of cylinder

# number of cells in arc, radius and axial direction
arc_number = 10
radius_number = 5
axis_number = 10

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


#================= main program =================
XT = []
YT = []
ZT = []



# gas inlet
for qua in range(1,5): 
    x, y, z = o_grid(qua, 14.11, 0.5*14.11, 20, 0, 6, 15, 28, 10)
    XT.append(x)
    YT.append(y)
    ZT.append(z)

x, y, z = h_grid(0.5*14.11, 20, 0, 6, 15, 10)
XT.append(x)
YT.append(y)
ZT.append(z)

# liquid inlet
xx, yy, zz = liquid_cylinder(r_orifice, R_liquidout, R_liquidin, num_s, num_r, num_rl, num_arc, num_pipe, inlet_pipe_len)

XT += xx
YT += yy
ZT += zz

# liquid cylinder
for qua in range(1,5): 
    x, y, z = o_grid(qua, 15.61, 15.11, 90, -1, 5, 15, 2, 10)
    XT.append(x)
    YT.append(y)
    ZT.append(z)

# gas recess
for qua in range(1,5): 
    x, y, z = o_grid(qua, 14.11, 0.5*14.11, 20, -6, 28.22, 15, 28, 50)
    XT.append(x)
    YT.append(y)
    ZT.append(z)

x, y, z = h_grid(0.5*14.11, 20, -6, 28.22, 15, 50)
XT.append(x)
YT.append(y)
ZT.append(z)

# wall recess
for qua in range(1,5): 
    x, y, z = o_grid(qua, 15.11, 14.11, 90, -6, 28.22, 15, 4, 50)
    XT.append(x)
    YT.append(y)
    ZT.append(z)

# liquid recess
for qua in range(1,5): 
    x, y, z = o_grid(qua, 15.61, 15.11, 90, -6, 28.22, 15, 2, 50)
    XT.append(x)
    YT.append(y)
    ZT.append(z)

# gas out
for qua in range(1,5): 
    x, y, z = o_grid(qua, 14.11, 0.5*14.11, 20, -34.22, 50, 15, 28, 85)
    XT.append(x)
    YT.append(y)
    ZT.append(z)

x, y, z = h_grid(0.5*14.11, 20, -34.22, 50, 15, 85)
XT.append(x)
YT.append(y)
ZT.append(z)

# wall out
for qua in range(1,5): 
    x, y, z = o_grid(qua, 15.11, 14.11, 90, -34.22, 50, 15, 4, 85)
    XT.append(x)
    YT.append(y)
    ZT.append(z)

# liquid out
for qua in range(1,5): 
    x, y, z = o_grid(qua, 15.61, 15.11, 90, -34.22, 50, 15, 2, 85)
    XT.append(x)
    YT.append(y)
    ZT.append(z)

# out
for qua in range(1,5): 
    x, y, z = o_grid(qua, 25, 15.61, 90, -34.22, 50, 15, 40, 85)
    XT.append(x)
    YT.append(y)
    ZT.append(z)

# establish fort12
txt = binary_fort12(XT, YT, ZT)
print(txt)



