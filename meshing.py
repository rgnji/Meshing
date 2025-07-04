from lib.cylindermesh import o_grid, h_grid
from lib.liquid_domain import liquid_cylinder
from lib.plot3dout import unformatted_fort12, binary_fort12, ascii_fort12

# arc_number == num_arc

#========== parameters for cylinder ==========
# gas inlet
RGSOUT = 14.11
RGSHGD = 0.5*RGSOUT
ZGSSTR = 0
LENGS  = 6
DARCGS = 20
DRADGS = 10
DAXLGS = 14

# liquid swirl
RLQOUT = 15.61
RLQINN = 15.11
ZLQSTR = -1
LENLQ  = 5
DARCLQ = DARCGS
DRADLQ = 3
DAXLLQ = 15

# gas recess
ZRCSTR = -6
LENRC  = 28.22
DRADWL = 3
DAXLRC = 60

#========== parameters for liquid domain ==========
r_orifice = 0.5 # inlet orifice radius
R_liquidout = 15.61 # liquid domain outer radius
R_liquidin = 15.11  # liquid domain inner radius
inlet_pipe_len = 15.11 # inlet pipe length
num_s = 4    # cells on orifice arc
num_r = 5    # cells on orifice radius
num_rl = DRADLQ   # cells on liquid domain thickness
num_arc = DARCLQ  # cells on liquid domain arc
num_pipe = 20 # inlet pipe cells



#================= grid generation =================
XT = []
YT = []
ZT = []

# meshing parameters
"""
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

mesh = [['o', 14.11,     0.5*14.11, 20,     0,      6,     10, 7, 7], # gas inlet
        ['h', 0.5*14.11,            20,     0,      6,     10,    7], # gas inlet
        ['o', 15.61,     15.11,     90,     -1,     5,     10, 2, 5], # liquid swirl
        ['o', 14.11,     0.5*14.11, 20,     -6,     28.22, 10, 7, 12], # gas recess
        ['h', 0.5*14.11,            20,     -6,     28.22, 10,    12], # gas recess
        ['o', 15.11,     14.11,     90,     -6,     28.22, 10, 3, 12], # wall recess
        ['o', 15.61,     15.11,     90,     -6,     28.22, 10, 2, 12], # liquid recess
        ['l', r_orifice, R_liquidout, R_liquidin, num_s, num_r, num_rl, num_arc, num_pipe, inlet_pipe_len] # liquid inlet
        ]
"""
mesh = [['o', RGSOUT, RGSHGD, 20, ZGSSTR, LENGS,  DARCGS, DRADGS, DAXLGS, 4, 1], # gas inlet
        ['h', RGSHGD,         20, ZGSSTR, LENGS,  DARCGS,         DAXLGS,    1], # gas inlet
        ['o', RLQOUT, RLQINN, 90, ZLQSTR, LENLQ,  DARCLQ, DRADLQ, DAXLLQ, 1, 4], # liquid swirl
        ['o', RGSOUT, RGSHGD, 20, ZRCSTR, LENRC,  DARCGS, DRADGS, DAXLRC, 4, 1], # gas recess
        ['h', RGSHGD,         20, ZRCSTR, LENRC,  DARCGS,         DAXLRC,    1], # gas recess
        ['o', RLQINN, RGSOUT, 90, ZRCSTR, LENRC,  DARCGS, DRADWL, DAXLRC, 1, 1], # wall recess
        ['o', RLQOUT, RLQINN, 90, ZRCSTR, LENRC,  DARCLQ, DRADLQ, DAXLRC, 1, 1], # liquid recess
        ['l', r_orifice, R_liquidout, R_liquidin, num_s, num_r, num_rl, num_arc, num_pipe, inlet_pipe_len, 1] # liquid inlet
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
"""
txt = ascii_fort12(XT, YT, ZT)
print(txt)
"""