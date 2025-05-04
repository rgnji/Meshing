import numpy as np
import matplotlib.pyplot as plt

from lib.cylindermesh import o_grid, h_grid
from lib.plot3dout import unformatted_fort12, binary_fort12
from lib.liquid_domain import liquid_cylinder
from lib.controlcards import gridz_size, patch, o_patch, h_patch, o_inner_patch, h_patch_liquid

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
# o_patch
group4_cards_o = [1, 6, 10, 15, 19, 23, 28, 32, 36]
# h_patch
group4_cards_h = [1, 10, 23]
# o_inner_patch
group4_cards_o_in = [[15, 10], [19, 15], [28, 23], [32, 28], [36, 32]]
# patch
group4_cards_z = [[1, 10, 3], [2, 11, 3], [3, 12, 3], [4, 13, 3], [5, 14, 3],
                  [6, 19, 3], [7, 20, 3], [8, 21, 3], [9, 22, 3],
                  [10, 23, 3], [11, 24, 3], [12, 25, 3], [13, 26, 3], [14, 27, 3],
                  [15, 28, 3], [16, 29, 3], [17, 30, 3], [18, 31, 3],
                  [19, 32, 3], [20, 33, 3], [21, 34, 3], [22, 35, 3],
                  ]
group4_cards_l = [[42, 40, 2], [44, 40, 1], [40, 46, 3], # pipe, outer, o-grid
                  [48, 46, 2], [50, 46, 1], [46, 52, 3],
                  [54, 52, 2], [56, 52, 1], [52, 58, 3],
                  [60, 58, 2], [62, 58, 1], [58, 40, 3],
                  [43, 41, 2], [45, 41, 1], [41, 47, 3],
                  [49, 47, 2], [51, 47, 1], [47, 53, 3],
                  [55, 53, 2], [57, 53, 1], [53, 59, 3],
                  [61, 59, 2], [63, 59, 1], [59, 41, 3],
                  [42, 48, 3], [48, 54, 3], [54, 60, 3], [60, 42, 3], # pipe 1 o-grid                 
                  [43, 49, 3], [49, 55, 3], [55, 61, 3], [61, 43, 3], # pipe 2 o-grid
                  [62, 44, 3], [50, 56, 3], [63, 45, 3], [51, 57, 3] # outer
                  ]
# h_patch_liquid
group4_liquid_h = [[40, 64], [41, 65], [42, 66], [43, 67]]
# patch (izf1 and 2 = 2)
group4_outer = [[44, 51], [62, 57], [50, 45], [56, 63]]


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

    group4_output = []

    for i in group4_cards_o:
        group4_output.extend(o_patch(XT, i))
    for i in group4_cards_h:
        group4_output.extend(h_patch(XT, i))
    for i in group4_cards_o_in:
        group4_output.extend(o_inner_patch(XT, *i))
    for i in group4_cards_z:
        if i[2] == 1:
            izf1 = 1
            izf2 = 2
        elif i[2] == 2:
            izf1 = 3
            izf2 = 4
        elif i[2] == 3:
            izf1 = 5
            izf2 = 6
        group4_output.extend(patch(XT, i[0], i[1], izf1, izf2, 1, 1, 1, 1))
    for i in group4_cards_l:
        if i[2] == 1:
            izf1 = 1
            izf2 = 2
        elif i[2] == 2:
            izf1 = 3
            izf2 = 4
        elif i[2] == 3:
            izf1 = 5
            izf2 = 6
        group4_output.extend(patch(XT, i[0], i[1], izf1, izf2, 1, 1, 1, 1))
    for i in group4_liquid_h:
        group4_output.extend(h_patch_liquid(XT, *i))
    for i in group4_outer:
        group4_output.extend(patch(XT, i[0], i[1], 2, 2, 1, 1, 1, 1))
    
    print(len(group4_output))
    print(group4_output[0][1])
    
    for item in group4_output:
        """
        print(type(group4_output[item][0])) # list ???
        break
        """
        for item_1 in item[0]:
            f.write(f'{item_1:>6},')
        f.write('\n')
        f.write('       ')
        for item_2 in item[1]:
            f.write(f'{item_2:>6},')
        f.write('\n')
        
    print("fort11.txt established.")
