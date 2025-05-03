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


#================= grid generation =================
XT = []
YT = []
ZT = []

# group 2
izon_count = 0
ibnd_count = (5+5+5)+(8+4+4+5)
id_count = (4+(4+10)*2+8)+(4+4)+(4)

# group 3
IZT = []
JZT = []
KZT = []

# meshing parameters
mesh = [['o', 14.11, 0.5*14.11, 20, 0, 6, 15, 28, 10], # gas inlet
        ['h', 0.5*14.11, 20, 0, 6, 15, 10], # gas inlet
        ['l', r_orifice, R_liquidout, R_liquidin, num_s, num_r, num_rl, num_arc, num_pipe, inlet_pipe_len], # liquid inlet
        ['o', 15.61, 15.11, 90, -1, 5, 15, 2, 10], # liquid cylinder
        ['o', 14.11, 0.5*14.11, 20, -6, 28.22, 15, 28, 50], # gas recess
        ['h', 0.5*14.11, 20, -6, 28.22, 15, 50], # gas recess
        ['o', 15.11, 14.11, 90, -6, 28.22, 15, 4, 50], # wall recess
        ['o', 15.61, 15.11, 90, -6, 28.22, 15, 2, 50], # liquid recess
        ['o', 14.11, 0.5*14.11, 20, -34.22, 50, 15, 28, 85], # gas out
        ['h', 0.5*14.11, 20, -34.22, 50, 15, 85], # gas out
        ['o', 15.11, 14.11, 90, -34.22, 50, 15, 4, 85], # wall out
        ['o', 15.61, 15.11, 90, -34.22, 50, 15, 2, 85], # liquid out
        ['o', 25, 15.61, 90, -34.22, 50, 15, 40, 85], # out
        ]


# generate mesh
for grid in mesh:
    if grid[0] == 'o':
        for qua in range(1,5): 
            x, y, z = o_grid(qua, *grid[1:])
            XT.append(x)
            YT.append(y)
            ZT.append(z)
            izon_count += 1
            IZT.append(x.shape[2])
            JZT.append(x.shape[1])
            KZT.append(x.shape[0])
    
    if grid[0] == 'h':
        x, y, z = h_grid(*grid[1:])
        XT.append(x)
        YT.append(y)
        ZT.append(z)
        izon_count += 1
        IZT.append(x.shape[2])
        JZT.append(x.shape[1])
        KZT.append(x.shape[0])

    if grid[0] == 'l':
        xx, yy, zz = liquid_cylinder(*grid[1:])
        XT += xx
        YT += yy
        ZT += zz
        for i in range(len(xx)):
            izon_count += 1
            IZT.append(xx[i].shape[2])
            JZT.append(xx[i].shape[1])
            KZT.append(xx[i].shape[0])

# establish fort12
txt = binary_fort12(XT, YT, ZT)
print(txt)


def o_patch(startblock, endblock, direction):
    IFCYC = 1

    if direction == 1: # I
        f = 1
        first = 1
        second = 0
    elif direction == 2: # J
        f = 3
        first = 2
        second = 0
    elif direction == 3: # K
        f = 5
        first = 2
        second = 1

    IZB1 = startblock
    IZF1 = f
    IZB2 = endblock
    IZF2 = f + 1
    
    IJZ11 = 1
    IJZ12 = XT[startblock].shape[first] 
    JKZ11 = 1
    JKZ12 = XT[startblock].shape[second] 
    IJZ21 = 1
    IJZ22 = XT[startblock].shape[first] 
    JKZ21 = 1
    JKZ22 = XT[startblock].shape[second] 

    INONUF = 0

    return [[IFCYC, IZB1, IZF1, IJZ11, IJZ12, JKZ11, JKZ12, INONUF],
            [IZB2, IZF2, IJZ21, IJZ22, JKZ21, JKZ22]]


#================= fort.11 =================
with open("fort11.txt", 'w', encoding='UTF-8') as f:
    group1 = ['TITLE', 'IDIM']
    group2 = ['IZON', 'IZFACE', 'IBND', 'ID', 'ISNGL']
    group3 = ['IZT', 'JZT', 'KZT', 'LPROC', 'CBG1', 'CBG2', 'CBG3', 'CBV1', 'CBV2', 'CBV3']
    group4_1 = ['IFCYC', 'IZB1', 'IZF1', 'IJZ11', 'IJZ12', 'JKZ11', 'JKZ12', 'INONUF']
    group4_2 = ['', 'IZB2', 'IZF2', 'IJZ21', 'IJZ22', 'JKZ21', 'JKZ22', '']

    # group 1
    f.write(f'{group1[0]}:')
    f.write(' GCSC Injector A\n')

    f.write(f'{group1[1]},\n')
    f.write(f'{3:>4},\n')

    # group 2
    for g2 in group2:
        f.write(f'{g2:>6},')
    f.write('\n')

    f.write(f'{izon_count:>6},')

    izface_count = (izon_count * 6 - (id_count + ibnd_count)) * 0.5
    f.write(f'{izface_count:>6},')

    f.write(f'{ibnd_count:>6},')

    f.write(f'{id_count:>6},')

    f.write(f'{0:>6},')

    f.write('\n')

    # group 3
    for g3 in group3:
        f.write(f'{g3:>6},')
    f.write('\n')

    dim = [IZT, JZT, KZT]

    for i in range(len(IZT)):
        for j in range(3):
            f.write(f'{dim[j][i]:>6},')

        if i < 48:
            f.write(f'{1:>6},')
        else:
            f.write(f'{2:>6},')
        
        for k in range(6):
            f.write(f'{0.:>6},')
        
        f.write('\n')
    
    # group 4
    for g4 in group4_1:
        f.write(f'{g4:>6},')
    f.write('\n')
    for g4 in group4_2:
        f.write(f'{g4:>6},')
    f.write('\n')

    
    
    print("fort11.txt established.")
