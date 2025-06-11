from lib.cylindermesh import o_grid, h_grid
from lib.liquid_domain import liquid_cylinder
from lib.plot3dout import unformatted_fort12, binary_fort12, ascii_fort12

# arc_number == num_arc

#========== parameters for cylinder ==========


#========== parameters for liquid domain ==========
r_orifice = 0.5E-3 # inlet orifice radius
R_liquidout = 15.61E-3 # liquid domain outer radius
R_liquidin = 15.11E-3  # liquid domain inner radius
num_s = 3    # cells on orifice arc
num_r = 3    # cells on orifice radius
num_rl = 2   # cells on liquid domain thickness
num_arc = 10 # cells on liquid domain arc
num_pipe = 10 # inlet pipe cells
inlet_pipe_len = 15.11E-3 # inlet pipe length


#================= grid generation =================
XT = []
YT = []
ZT = []

# meshing parameters
mesh = [['o', 14.11E-3,     0.5*14.11E-3, 20,        0,         6E-3,     10, 7, 7], # gas inlet
        ['h', 0.5*14.11E-3, 20,           0,         6E-3,      10,       7], # gas inlet
        ['o', 15.61E-3,     15.11E-3,     90,        -1E-3,     5E-3,     10, 2,  5], # liquid cylinder
        ['o', 14.11E-3,     0.5*14.11E-3, 20,        -6E-3,     28.22E-3, 10, 7, 12], # gas recess
        ['h', 0.5*14.11E-3, 20,           -6E-3,     28.22E-3,  10,       12], # gas recess
        ['o', 15.11E-3,     14.11E-3,     90,        -6E-3,     28.22E-3, 10, 3,  12], # wall recess
        ['o', 15.61E-3,     15.11E-3,     90,        -6E-3,     28.22E-3, 10, 2,  12], # liquid recess
        ['o', 14.11E-3,     0.5*14.11E-3, 20,        -34.22E-3, 20E-3,    10, 7, 10], # gas out
        ['h', 0.5*14.11E-3, 20,           -34.22E-3, 20E-3,     10,       10], # gas out
        ['o', 15.11E-3,     14.11E-3,     90,        -34.22E-3, 20E-3,    10, 3,  10], # wall out
        ['o', 15.61E-3,     15.11E-3,     90,        -34.22E-3, 20E-3,    10, 2,  10], # liquid out
        ['o', 25E-3,        15.61E-3,     90,        -34.22E-3, 20E-3,    10, 10, 10], # out
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
txt = ascii_fort12(XT, YT, ZT)
print(txt)