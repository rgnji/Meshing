from lib.cylindermesh import h_grid, simplified_o_grid
from lib.liquid2 import liquid_inlet, pipe, wall
from lib.plot3dout import unformatted_fort12, binary_fort12
import numpy as np

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
DARCLQ = DARCGS * 2
DRADLQ = 4
DAXLLQ = 15

# gas recess
ZRCSTR = -6
LENRC  = 28.22
DRADWL = 4
DAXLRC = 60

#========== parameters for liquid domain ==========
PIPRIN = 0.5 # inlet orifice radius
PIPLEN = 5 # inlet pipe length
LQROUT = 15.61 # liquid domain outer radius
LQRINN = 15.11  # liquid domain inner radius
CNTINT = np.array([LQRINN, -PIPLEN, -PIPRIN])
NARC = 4    # cells on orifice arc
NRAD = 5    # cells on orifice radius
NPIP = 30 # inlet pipe cells

#================= grid generation =================
XT = []
YT = []
ZT = []

mesh = [["H", RGSOUT, 90, ZGSSTR, LENGS, DARCGS, DAXLGS],
        ["O", RLQOUT, RLQINN, ZLQSTR, LENLQ, DARCLQ, DRADLQ, DAXLLQ],
        ["H", RGSOUT, 90, ZRCSTR, LENRC, DARCGS, DAXLRC],
        ["O", RLQOUT, RLQINN, ZRCSTR, LENRC, DARCLQ, DRADLQ, DAXLRC],
        ["O", RLQINN, RGSOUT, ZRCSTR, LENRC, DARCLQ, DRADWL, DAXLRC],
        ["L", PIPRIN, PIPLEN, CNTINT, LQROUT, LQRINN, NARC, DRADLQ, NPIP, DARCLQ]]

for item in mesh:
    if item[0] == "H":
        x, y, z = h_grid(*item[1:])
        XT.append(x)
        YT.append(y)
        ZT.append(z)
    if item[0] == "O":
        for i in range(2):
            x, y, z = simplified_o_grid(i+1, *item[1:])
            XT.append(x)
            YT.append(y)
            ZT.append(z)
    if item[0] == "L":
        x, y, z = liquid_inlet(*item[1:])
        XT += x
        YT += y
        ZT += z
    
txt = binary_fort12(XT, YT, ZT)
print(txt)
txt = unformatted_fort12(XT, YT, ZT)
print(txt)