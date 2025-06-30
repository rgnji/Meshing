import numpy as np
from struct import unpack

# =============================================
# ========== global variables =================
# =============================================
processor_num = 5
# ===== group 5 =====
# zonal index of flow boundary (1-based)
# information can be obtained from paraview
# inlet, outlet
IBCZON = [1,2,3,4,5,
          25,26,31,32,37,38,43,44,49,50,
          10,11,12,13,14,
          15,16,17,18,
          19,20,21,22]
IDBC = [5,5,5,5,5,
        5,5,5,5,5,5,5,5,5,5,
        6,6,6,6,6,
        6,6,6,6,
        6,6,6,6]
ITYBC = [-1,-1,-1,-1,-1,
         -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
         2,2,2,2,2,
         2,2,2,2,
         2,2,2,2]
# ===== group 6 =====
IWTM = 1
HQDOX = 0
IWALL = 0
DENNX = 0
VISWX = 0
# ===== group 8 =====
IDATA = 1
IGEO = 2    # any number > 0 but 1 and 9, which are for special cases
ITT = 1E4 
ITPNT = 50
ICOUP = 1   # for unstable flow, ICOUP > 1
NLIMT = 1
IAX = 1
ICYC = 0
# ===== group 9 =====
# reference time = XREF/UREF = 1.4E-5
DTT = 1E-5
IREC = 3
REC = 0.1
THETA = 1.0 #
BETAP = 1   # compressible for real fluid model
IEXX = 1 #
PRAT = 0 #
# ===== group 10 =====
IPC = 5793 # (I,J,K,IZ)=(18, 3,14, 1)
JPC = 1
IPEX = 199 # (I,J,K,IZ)=(10,10, 1,14)
JPEX = 14
IMN = 5793
JMN = 1
# ===== group 11 =====
VISC = 18.37e-6
IG = 2
ITURB = 1   # use low-Re model if calculating viscous boundary layer
AMC = 1     # compressible for real fluid model
GAMA = 1.455
CBE = 0
CBH = 0 #
EREXT = 5.0e-6
# ===== group 12 =====
ISWU = 93
ISWP = 97
ISWK = 93
ISKEW = 0
# ===== group 13 =====
U = 1
V = 1
W = 1
TM = 1
DK = 1
DE = 1
FL = 0
EQ = 0
VS = 1
FM = 1
SP = 1
# ===== group 14 =====
NGAS = 2 # h2o and air
NREACT = 0
IUNIT = 1
DENREF = 1.1550255578176354
UREF = 69.593
TREF = 300
XREF = 1E-3
PREF = 1E5
# ===== group 17 =====
IGDINN = 1
IOFINN = 1
IOFOUT = 1
IOFP3D = 1
# =============================================
# ================== read in ==================
# =============================================
with open("fort.12", 'rb') as f:
    buff = f.read(4*3)
    blocks = unpack('<3i', buff)[1]
    
    dim = []
    for i in range(blocks):
        buff = f.read(4*5)
        dim.append(unpack('<5i', buff)[1:4])
    
    XT, YT, ZT =[], [], []
    for i in range(blocks):
        size = dim[i][0]*dim[i][1]*dim[i][2]
        
        f.read(4)
        buff = f.read(4*size)
        xflat = np.array(unpack(f'<{size}f', buff))
        XT.append(xflat.reshape((dim[i][2], dim[i][1], dim[i][0])))
        f.read(4)
        
        f.read(4)
        buff = f.read(4*size)
        yflat = np.array(unpack(f'<{size}f', buff))
        YT.append(yflat.reshape((dim[i][2], dim[i][1], dim[i][0])))
        f.read(4)
        
        f.read(4)
        buff = f.read(4*size)
        zflat = np.array(unpack(f'<{size}f', buff))
        ZT.append(zflat.reshape((dim[i][2], dim[i][1], dim[i][0])))
        f.read(4)
# =============================================
# =============================================
# =============================================
# ===== group 2 =====
# eight vertices of each block
verticesX = []
verticesY = []
verticesZ = []
loc = [ [ 0, 0, 0],
        [-1, 0, 0],
        [-1,-1, 0],
        [ 0,-1, 0],
        [ 0, 0,-1],
        [-1, 0,-1],
        [-1,-1,-1],
        [ 0,-1,-1]]
for block in range(blocks):
    x = []
    y = []
    z = []
    for i, j, k in loc:
        x.append(XT[block][k][j][i])
        y.append(YT[block][k][j][i])
        z.append(ZT[block][k][j][i])
    verticesX.append(x)
    verticesY.append(y)
    verticesZ.append(z)

# neighboring blocks in the six direction
patched_interface = []
direction = [[3,6,7],
             [1,4,5],
             [3,7,8],
             [1,2,5],
             [6,7,8],
             [1,2,4]]
decimals = 5

for target in range(blocks):
    patch = [] # shared faces, 6 values

    # search for patched blocks in six direction
    for d in range(6):
        nb = [] # shared vertex, 3 sets

        # find the blocks sharing the same vertices of direction[d]
        for v in direction[d]:
            target_x = round(verticesX[target][v-1], decimals)
            target_y = round(verticesY[target][v-1], decimals)
            target_z = round(verticesZ[target][v-1], decimals)
            neighbor = set()

            blks = list(range(blocks))
            blks.remove(target)
            for block in blks:
                for i in range(8):
                    if (
                        target_x == round(verticesX[block][i], decimals) and
                        target_y == round(verticesY[block][i], decimals) and
                        target_z == round(verticesZ[block][i], decimals)
                    ):
                        neighbor.add(block)
                        break
            nb.append(neighbor)

        # find the block sharing the same three vertices
        if nb[0] & nb[1] & nb[2]:
            patch.append(list(nb[0] & nb[1] & nb[2])[0])
        else:
            patch.append(-1)
    
    patched_interface.append(patch)


IZON = blocks
IZFACE = 0
for blk4 in range(IZON):
        for d4 in range(6):
            if patched_interface[blk4][d4] and blk4 < patched_interface[blk4][d4]:
                IZFACE += 1
IBND = len(IBCZON)
ID = IZON * 6 - IZFACE * 2 - IBND
ISNGL = 0

# ===== group 4 =====
def gp4_re(block, face):
    
    def start(direction):
        if direction in [2,4,6]:
            vertex = 1
        elif direction == 1:
            vertex = 2
        elif direction == 3:
            vertex = 4
        elif direction == 5:
            vertex = 5
        return vertex
    
    # block 1 starting point
    target_start = start(face)
    target_vx = round(verticesX[block][target_start-1], decimals)
    target_vy = round(verticesY[block][target_start-1], decimals)
    target_vz = round(verticesZ[block][target_start-1], decimals)
    
    # block 2 face
    neighbor = patched_interface[block][face-1]
    face_neighbor = patched_interface[neighbor].index(block) + 1
    
    # block 2 shared point
    for neighbor_vertex in range(8): # 0-based
        neighbor_vx = round(verticesX[neighbor][neighbor_vertex], decimals)
        neighbor_vy = round(verticesY[neighbor][neighbor_vertex], decimals)
        neighbor_vz = round(verticesZ[neighbor][neighbor_vertex], decimals)
        if neighbor_vx == target_vx and neighbor_vy == target_vy and neighbor_vz == target_vz:
            break
    
    if face_neighbor in [1,2]:
        d = [1, 2] # j, k (dim)
        jj = loc[neighbor_vertex][1]
        kk = loc[neighbor_vertex][2]
    elif face_neighbor in [3,4]:
        d = [0, 2] # i, k
        jj = loc[neighbor_vertex][0]
        kk = loc[neighbor_vertex][2]
    elif face_neighbor in [5,6]:
        d = [0, 1] # i, j
        jj = loc[neighbor_vertex][0]
        kk = loc[neighbor_vertex][1]
    
    if jj == -1:
        ijz21 = dim[neighbor][d[0]]
        ijz22 = 1
    else:
        ijz21 = 1
        ijz22 = dim[neighbor][d[0]]
    if kk == -1:
        jkz21 = dim[neighbor][d[1]]
        jkz22 = 1
    else:
        jkz21 = 1
        jkz22 = dim[neighbor][d[1]]
    
    return ijz21, ijz22, jkz21, jkz22
    

# =============================================
# ================== fort.11 ==================
# =============================================
with open("fort.11", "w", encoding="UTF-8") as f:
    title = "GCSC INJECTOR A"

    # ================== group 1 ==================
    f.write(f"TITLE: {title}\n")
    f.write(f'IDIM,\n{3:>4},\n')

    # ================== group 2 ==================
    group2 = ['IZON', 'IZFACE', 'IBND', 'ID', 'ISNGL']
    for g2 in group2:
        f.write(f'{g2:>6},')
    f.write('\n')
    f.write(f'{IZON:>6},{IZFACE:>6},{IBND:>6},{ID:>6},{ISNGL:>6},\n')

    # ================== group 3 ==================
    group3 = ['IZT', 'JZT', 'KZT', 'LPROC', 'CBG1', 'CBG2', 'CBG3', 'CBV1', 'CBV2', 'CBV3']

    for g3 in group3:
        f.write(f'{g3:>6},')
    f.write('\n')

    processor = 0
    for blk3 in range(blocks):
        f.write(f'{dim[blk3][0]:>6},{dim[blk3][1]:>6},{dim[blk3][2]:>6},')

        f.write(f'{(processor % processor_num) + 1:>6},')
        processor += 1

        for gg3 in range(6):
            f.write(f'{0.:>6},')

        f.write('\n')

    # ================== group 4 ==================
    group4_1 = ['IFCYC', 'IZB1', 'IZF1', 'IJZ11', 'IJZ12', 'JKZ11', 'JKZ12', 'INONUF']
    group4_2 = ['IZB2', 'IZF2', 'IJZ21', 'IJZ22', 'JKZ21', 'JKZ22']

    for g4 in group4_1:
        f.write(f'{g4:>6},')
    f.write('\n')
    f.write(f'{"":>6} ')
    for g4 in group4_2:
        f.write(f'{g4:>6},')
    f.write('\n')

    IFCYC = 1
    INONUF = 0
    for blk4 in range(blocks):
        for d4 in range(6):
            if patched_interface[blk4][d4] and blk4 < patched_interface[blk4][d4]:
                IZB1 = blk4 + 1
                IZF1 = d4 + 1
                IZB2 = patched_interface[blk4][d4] + 1
                IZF2 = patched_interface[IZB2-1].index(blk4) + 1
                
                IJZ11 = 1
                JKZ11 = 1
                if (d4+1) in [1,2]:
                    IJZ12 = dim[blk4][1]
                    JKZ12 = dim[blk4][2]
                elif (d4+1) in [3,4]:
                    IJZ12 = dim[blk4][0]
                    JKZ12 = dim[blk4][2]
                elif (d4+1) in [5,6]:
                    IJZ12 = dim[blk4][0]
                    JKZ12 = dim[blk4][1]
                
                IJZ21, IJZ22, JKZ21, JKZ22 = gp4_re(blk4, d4+1)

                f.write(f'{IFCYC:>6},{IZB1:>6},{IZF1:>6},')
                f.write(f'{IJZ11:>6},{IJZ12:>6},{JKZ11:>6},{JKZ12:>6},')
                f.write(f'{INONUF:>6},\n')
                f.write(f'{"":>6} {IZB2:>6},{IZF2:>6},{IJZ21:>6},{IJZ22:>6},{JKZ21:>6},{JKZ22:>6},')
                f.write('\n')

    # ================== group 5 ==================
    group5 = ['IBCZON','IDBC','ITYBC','IJBB','IJBS','IJBT','JKBS','JKBT']

    for g5 in group5:
        f.write(f'{g5:>7},')
    f.write('\n')

    for g5 in range(len(IBCZON)):
        f.write(f'{IBCZON[g5]:>7},{IDBC[g5]:>7},{ITYBC[g5]:>7},')
        
        if IDBC[g5] in [2,4,6]:
            IJBB = 1
        elif IDBC[g5] == 1:
            IJBB = dim[ IBCZON[g5]-1 ][0]
        elif IDBC[g5] == 3:
            IJBB = dim[ IBCZON[g5]-1 ][1]
        elif IDBC[g5] == 5:
            IJBB = dim[ IBCZON[g5]-1 ][2]

        IJBS = JKBS = 1
        if IDBC[g5] in [1,2]:
            IJBT = dim[ IBCZON[g5]-1 ][1]
            JKBT = dim[ IBCZON[g5]-1 ][2]
        elif IDBC[g5] in [3,4]:
            IJBT = dim[ IBCZON[g5]-1 ][0]
            JKBT = dim[ IBCZON[g5]-1 ][2]
        elif IDBC[g5] in [5,6]:
            IJBT = dim[ IBCZON[g5]-1 ][0]
            JKBT = dim[ IBCZON[g5]-1 ][1]

        f.write(f'{IJBB:>7},{IJBS:>7},{IJBT:>7},{JKBS:>7},{JKBT:>7},\n')

    # ================== group 6 ==================
    group6 = ['IWBZON','L1','L2','M1','M2','N1','N2','IWTM','HQDOX','IWALL','DENNX','VISWX']
    for g6 in group6:
        f.write(f'{g6:>7},')
    f.write('\n')

    def gp6(blk, face):
        if face == 1:
            istart= iend = dim[blk][0]
            jstart = kstart = 1
            jend = dim[blk][1]
            kend = dim[blk][2]
        elif face == 2:
            istart = iend = 1
            jstart = kstart = 1
            jend = dim[blk][1]
            kend = dim[blk][2]
        elif face == 3:
            jstart = jend = dim[blk][1]
            istart = kstart = 1
            iend = dim[blk][0]
            kend = dim[blk][2]
        elif face == 4:
            jstart = jend = 1
            istart = kstart = 1
            iend = dim[blk][0]
            kend = dim[blk][2]
        elif face == 5:
            kstart = kend = dim[blk][2]
            istart = jstart = 1
            iend = dim[blk][0]
            jend = dim[blk][1]
        elif face == 6:
            kstart = kend = 1
            istart = jstart = 1
            iend = dim[blk][0]
            jend = dim[blk][1]
        return istart, iend, jstart, jend, kstart, kend

    for blk6 in range(blocks):
        for face6 in range(6):
            if patched_interface[blk6][face6] == -1:
                if (blk6 + 1) not in IBCZON:
                    IWBZON = blk6 + 1
                    L1, L2, M1, M2, N1, N2 = gp6(blk6, face6 + 1)
                    f.write(f'{IWBZON:>7},')
                    f.write(f'{L1:>7},{L2:>7},{M1:>7},{M2:>7},{N1:>7},{N2:>7},')
                    f.write(f'{IWTM:>7},{HQDOX:>7},{IWALL:>7},{DENNX:>7},{VISWX:>7},\n')

                else:
                    index6 = [i for i, x in enumerate(IBCZON) if x == (blk6 + 1)]
                    if (face6 + 1) not in [x for i, x in enumerate(IDBC) if i in index6]:
                        IWBZON = blk6 + 1
                        L1, L2, M1, M2, N1, N2 = gp6(blk6, face6 + 1)
                        f.write(f'{IWBZON:>7},')
                        f.write(f'{L1:>7},{L2:>7},{M1:>7},{M2:>7},{N1:>7},{N2:>7},')
                        f.write(f'{IWTM:>7},{HQDOX:>7},{IWALL:>7},{DENNX:>7},{VISWX:>7},\n')
    
    # ================== group 7 ==================
    group7 = ['ISNZON', 'ISNBC', 'ISNAX', 'ISNBS', 'ISNBT']
    for g7 in group7:
        f.write(f'{g7:>6},')
    f.write('\n')

    # ================== group 8 ==================
    group8 = ['IDATA', 'IGEO', 'ITT', 'ITPNT', 'ICOUP', 'NLIMT', 'IAX', 'ICYC']
    for g8 in group8:
        f.write(f'{g8:>6},')
    f.write('\n')
    f.write(f'{IDATA:>6},{IGEO:>6},{ITT:>.3E},{ITPNT:>6},{ICOUP:>6},{NLIMT:>6},{IAX:>6},{ICYC:>6},\n')

    # ================== group 9 ==================
    group9 = ['DTT', 'IREC', 'REC', 'THETA', 'BETAP', 'IEXX', 'PRAT']
    for g9 in group9:
        f.write(f'{g9:>6},')
    f.write('\n')
    f.write(f'{DTT:>.3E},{IREC:>6},{REC:>6},{THETA:>6},{BETAP:>6},{IEXX:>6},{PRAT:>6},\n')

    # ================== group 10 ==================
    group10 = ['IPC', 'JPC', 'IPEX', 'JPEX', 'IMN', 'JMN']
    for g10 in group10:
        f.write(f'{g10:>6},')
    f.write('\n')
    f.write(f'{IPC:>6},{JPC:>6},{IPEX:>6},{JPEX:>6},{IMN:>6},{JMN:>6},\n')

    # ================== group 11 ==================
    group11 = ['VISC', 'IG', 'ITURB', 'AMC', 'GAMA', 'CBE', 'CBH', 'EREXT']
    for g11 in group11:
        f.write(f'{g11:>9},')
    f.write('\n')
    f.write(f'{VISC:>9},{IG:>9},{ITURB:>9},{AMC:>9},{GAMA:>9},{CBE:>9},{CBH:>9},{EREXT:>9},\n')

    # ================== group 12 ==================
    group12 = ['ISWU','ISWP','ISWK','ISKEW']
    for g12 in group12:
        f.write(f'{g12:>6},')
    f.write('\n')
    f.write(f'{ISWU:>6},{ISWP:>6},{ISWK:>6},{ISKEW:>6},\n')

    # ================== group 13 ==================
    group13 = ['U','V','W','TM','DK','DE','FL','8','EQ','VS','FM','SP']
    f.write('INSO(IEQ):\n')
    for g13 in group13:
        f.write(f'{g13:>3},')
    f.write('\n')
    f.write(f'{U:>3},{V:>3},{W:>3},{TM:>3},{DK:>3},{DE:>3},{FL:>3},{"0":>3},{EQ:>3},{VS:>3},{FM:>3},{SP:>3},\n')

    # ================== group 14 ==================
    group14 = ['NGAS','NREACT','IUNIT','DENREF','UREF','TREF','XREF','PREF']
    for g14 in group14:
        f.write(f'{g14:>6},')
    f.write('\n')
    f.write(f'{NGAS:>6},{NREACT:>6},{IUNIT:>6},{DENREF:>.4E},{UREF:>6},{TREF:>6},{XREF:>6},{PREF:>6},\n')

    # ================== group 15 ==================
    group15 = ['ISPARK','ISPKMIN','ISPMAX']
    for g15 in group15:
        f.write(f'{g15:>6},')
    f.write('\n')

    # ================== group 16 ==================
    group16 = ['ISPKON','ISPKZN','ISPKI1','ISPKIM','ISPKJ1','ISPKJM','ISPKK1','ISPKKM','TMSPK','ISPKDBG']
    for g16 in group16:
        f.write(f'{g16:>6},')
    f.write('\n')

    # ================== group 17 ==================
    group17 = ['IGDINN','IOFINN','IOFOUT','IOFP3D']
    for g17 in group17:
        f.write(f'{g17:>6},')
    f.write('\n')
    f.write(f'{IGDINN:>6},{IOFINN:>6},{IOFOUT:>6},{IOFP3D:>6},\n')

    # ================== group 18 ==================
    # polynomial coefficients are from dbase.dat
    name = ['H2O', 'AIR']
    mtmole = [18, 0.42*16+1.58*14]
    coef = [[.26340650E+01,  .31121900E-02, -.90278460E-06,  .12673050E-09,
            -.69164730E-14, -.29876260E+05,  .70823870E+01,  .41675560E+01,
            -.18106870E-02,  .59450870E-05, -.48670870E-08,  .15284140E-11,
            -.30289550E+05, -.73088000E+00,],
            [.30126640E+01,  .14229400E-02, -.53882480E-06,  .97471380E-10,
            -.66685940E-14, -.95470450E+03,  .58239980E+01,  .37210700E+01,
            -.17581880E-02,  .43543090E-05, -.30122530E-08,  .68264410E-12,
            -.10640230E+04,  .25293090E+01,]]
    for g18 in range(len(name)):
        f.write(f'{name[g18]:<66}')
        f.write(f'{mtmole[g18]:12.5f}\n')

        g18_count = 0
        for i in range(14):
            f.write(f'{coef[g18][i]:>15.8E}')
            g18_count += 1
            if g18_count % 5 == 0:
                f.write('\n')
        
        f.write('\n')

    # ================== group 19 ==================
    
    # ================== group 20 ==================

    # ================== group 21 ==================

    # ================== fluid entry ==================
    """
    f.write('FLUID\n')
    f.write('c Species(a20) , ideal gas(=0), real fluid(=1)\n')
    f.write(f"{'H2O':<20}{'1':>5}\n")
    f.write(f"{'AIR':<20}{'0':>5}\n")
    f.write(f"{'DONE':<20}{'0':>5}\n")
    """
    # ================== end ==================
    print('fort.11 established.')