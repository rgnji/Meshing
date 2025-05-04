import numpy as np
import matplotlib.pyplot as plt

from lib.cylindermesh import o_grid, h_grid
from lib.plot3dout import unformatted_fort12, binary_fort12
from lib.liquid_domain import liquid_cylinder
from lib.controlcards import gridz_size, patch, o_patch, h_patch, o_inner_patch

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
        ['o', 14.11, 0.5*14.11, 20, -34.22, 50, 15, 28, 85], # gas out
        ['h', 0.5*14.11, 20, -34.22, 50, 15, 85], # gas out
        ['o', 15.11, 14.11, 90, -34.22, 50, 15, 4, 85], # wall out
        ['o', 15.61, 15.11, 90, -34.22, 50, 15, 2, 85], # liquid out
        ['o', 25, 15.61, 90, -34.22, 50, 15, 40, 85], # out
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




#================= fort.11 =================
group4_cards_o = [1, 6, 10, 15, 19, 23, 28, 32, 36]
group4_cards_h = [1, 10, 23]
group4_cards_oin = [[15, 10], [19, 15], [28, 23], [32, 28], [36, 32]]
group4_cards_z = [[1, 10, 3], [2, 11, 3], [3, 12, 3], [4, 13, 3], [5, 14, 3],
                  [6, 19, 3], [7, 20, 3], [8, 21, 3], [9, 22, 3],
                  [10, 23, 3], [11, 24, 3], [12, 25, 3], [13, 26, 3], [14, 27, 3],
                  [15, 28, 3], [16, 29, 3], [17, 30, 3], [18, 31, 3],
                  [19, 32, 3], [20, 33, 3], [21, 34, 3], [22, 35, 3],
                  ]
group4_cards_l = [[40, 42, 3], [40, 44, 2], [40, 64, 2], [40, 46, 1], # pipe, outer, h-grid, o-grid
                  [46, 48, 3], [46, 50, 2], [46, 64, 2], [46, 52, 1],
                  [52, 54, 3], [52, 56, 2], [52, 64, 2], [52, 58, 1],
                  [58, 60, 3], [58, 62, 2], [58, 64, 2], [58, 40, 1],
                  [41, 43, 3], [41, 45, 2], [41, 65, 2], [41, 47, 1],
                  [47, 49, 3], [47, 51, 2], [47, 65, 2], [47, 53, 1],
                  [53, 55, 3], [53, 57, 2], [53, 65, 2], [53, 59, 1],
                  [59, 61, 3], [59, 63, 2], [59, 65, 2], [59, 41, 1],
                  [42, 48, 1], [48, 54, 1], [54, 60, 1], [60, 42, 1], # pipe 1 o-grid
                  [42, 66, 2], [48, 66, 2], [54, 66, 2], [60, 66, 2], # pipe 1 h-grid
                  [43, 49, 1], [49, 55, 1], [55, 61, 1], [61, 43, 1], # pipe 2 o-grid
                  [43, 67, 2], [49, 67, 2], [55, 67, 2], [61, 67, 2], # pipe 2 h-grid
                  [42, 48, ]]


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
    ibnd_count = (5+5+5)+(8+4+4+5)
    id_count = (4+(4+10)*2+8)+(4+4)+(4)

    for g2 in group2:
        f.write(f'{g2:>6},')
    f.write('\n')

    izface_count = (len(XT) * 6 - (id_count + ibnd_count)) * 0.5
    f.write(f'{len(XT):>6},')
    f.write(f'{izface_count:>6},')
    f.write(f'{ibnd_count:>6},')
    f.write(f'{id_count:>6},')
    f.write(f'{0:>6},')
    f.write('\n')

    # group 3
    for g3 in group3:
        f.write(f'{g3:>6},')
    f.write('\n')

    cards_3 = gridz_size(XT)

    for card in cards_3:
        for e in card:
            f.write(f'{e:>6},')
        f.write('\n')
            
    # group 4
    for g4 in group4_1:
        f.write(f'{g4:>6},')
    f.write('\n')
    for g4 in group4_2:
        f.write(f'{g4:>6},')
    f.write('\n')

    
    
    print("fort11.txt established.")
