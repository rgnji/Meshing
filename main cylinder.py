import numpy as np
import matplotlib.pyplot as plt

from lib.cylindermesh import o_grid, h_grid
from lib.plot3dout import unformatted_fort12, binary_fort12
from lib.liquid_domain import liquid_cylinder

#========== parameters ==========
# cylinder radius
r = 5
# diagonal of inner square (or radius of inner cylinder)
r_square = 0.5 * r
# axial length of cylinder
len_z = 6
pos_z = 0

# number of cells in arc, radius and axial direction
arc_number = 10
radius_number = 5
axis_number = 10

#================= main program =================
XT = []
YT = []
ZT = []

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
"""
for qua in range(1,5): 
    x, y, z = o_grid(qua, r, r_square, 20, pos_z, len_z, arc_number, radius_number, axis_number)
    XT.append(x)
    YT.append(y)
    ZT.append(z)

x, y, z = h_grid(r_square, 20, pos_z, len_z, arc_number, axis_number)
XT.append(x)
YT.append(y)
ZT.append(z)

# establish fort.12
txt = binary_fort12(XT, YT, ZT)
print(txt)
"""

XT, YT, ZT = liquid_cylinder()



# plot
for q in range(len(XT)):
    ax.scatter(XT[q], YT[q], ZT[q])

plt.show()
