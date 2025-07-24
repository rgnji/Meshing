import pandas as pd
import numpy as np
from struct import unpack

processor_num = 8
# ===== group 1 =====
title = 'INJECTOR A'
IDIM = 3
# ===== group 5 =====
# zonal index of flow boundary (1-based)
# information can be obtained from paraview
# inlet, outlet
"""
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
ITYBC = [-3,-3,-3,-3,-3,
         -3,-3,-3,-3,-3,-3,-3,-3,-3,-3,
         2,2,2,2,2,
         2,2,2,2,
         2,2,2,2]
"""
IBCZON = [8, 9,
          4, 5, 6, 7, 8]
IDBC = [5,5,
        6,6,6,6,6]
ITYBC = [-1,-1,
         2,2,2,2,2]

# ===== group 6 =====
IWTM = 1
HQDOX = 0
IWALL = 0
DENNX = 0
VISWX = 0
# ===== group 8 =====
IDATA = 1
IGEO = 2    # any number > 0 but 1 and 9, which are for special cases
ITT = 2E4
ITPNT = 50
ICOUP = 1   # for unstable flow, ICOUP > 1
NLIMT = 1
IAX = 1
ICYC = 0
# ===== group 9 =====
# reference time = XREF/UREF = 1.4E-5
DTT = 1E-10
IREC = 3
REC = 0.9
THETA = 1   # time marching scheme
BETAP = 1   # compressible for real fluid model
IEXX = 1 #
PRAT = 0 #
# ===== group 10 =====
"""
IPC = 5793 # (I,J,K,IZ)=(18, 3,14, 1)
JPC = 1
IPEX = 199 # (I,J,K,IZ)=(10,10, 1,14)
JPEX = 14
IMN = 5793
JMN = 1
"""
IPC = 9788
JPC = 1
IPEX = 324 # (I,J,K,IZ)=(10,10, 1,14)
JPEX = 4
IMN = 324
JMN = 1
# ===== group 11 =====
VISC = 18.37e-6
IG = 2      # laminar or turbulent
ITURB = 2   # use low-Re model if calculating viscous boundary layer
AMC = 1     # compressible for real fluid model
GAMA = 1.455
CBE = 0
CBH = -2 #
EREXT = 2E-8
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
FL = 1
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
with open("input files/injector A.bin", 'rb') as f:
    buff = f.read(4*3)
    IZON = unpack('<3i', buff)[1]
    
    dim = []
    buff = f.read(4)
    for i in range(IZON):
        buff = f.read(4*3)
        dim.append(unpack('<3i', buff))
    buff = f.read(4)
    
    XT, YT, ZT =[], [], []
    for i in range(IZON):
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
with open('input files/injector A.inp', 'r') as f:
    lines = f.readlines()

# =============================================
# =============================================
# =============================================
count_wall = 0
for line in lines:
    if line == 'I0:   GRID   SEGMENT    BCTYPE      JSTA      JEND      KSTA      KEND     NDATA\n':
        sta_i0 = count_wall + 1
        sta_idim = sta_i0 + 1 + IZON
        sta_j0 = sta_idim + 1 + IZON
        sta_jdim = sta_j0 + 1 + IZON
        sta_k0 = sta_jdim + 1 + IZON
        sta_kdim = sta_k0 + 1 + IZON
        
        bnd_1 = lines[sta_idim: sta_idim + IZON]
        bnd_2 = lines[sta_i0: sta_i0 + IZON]
        bnd_3 = lines[sta_jdim: sta_jdim + IZON]
        bnd_4 = lines[sta_j0: sta_j0 + IZON]
        bnd_5 = lines[sta_kdim: sta_kdim + IZON]
        bnd_6 = lines[sta_k0: sta_k0 + IZON]
        
    count_wall += 1

bnd_1 = [line.strip().split() for line in bnd_1]
bnd_2 = [line.strip().split() for line in bnd_2]
bnd_3 = [line.strip().split() for line in bnd_3]
bnd_4 = [line.strip().split() for line in bnd_4]
bnd_5 = [line.strip().split() for line in bnd_5]
bnd_6 = [line.strip().split() for line in bnd_6]

columns_wall_i = ['GRID', 'SEGMENT', 'BCTYPE', 'JSTA', 'JEND', 'KSTA', 'KEND', 'NDATA']
columns_wall_j = ['GRID', 'SEGMENT', 'BCTYPE', 'ISTA', 'IEND', 'KSTA', 'KEND', 'NDATA']
columns_wall_k = ['GRID', 'SEGMENT', 'BCTYPE', 'ISTA', 'IEND', 'JSTA', 'JEND', 'NDATA']

df_bnd_1 = pd.DataFrame(bnd_1, columns=columns_wall_i)
df_bnd_2 = pd.DataFrame(bnd_2, columns=columns_wall_i)
df_bnd_3 = pd.DataFrame(bnd_3, columns=columns_wall_j)
df_bnd_4 = pd.DataFrame(bnd_4, columns=columns_wall_j)
df_bnd_5 = pd.DataFrame(bnd_5, columns=columns_wall_k)
df_bnd_6 = pd.DataFrame(bnd_6, columns=columns_wall_k)

# =============================================
# =============================================
# =============================================
count = 0
face_1 = []
face_2 = []
for line in lines:
    
    if line == '   1-1 BLOCKING DATA:\n':
        nbli = int(lines[count+2])
   
        sta_1 = count + 4
        end_1 = sta_1 + nbli
        face_1 = lines[sta_1: end_1]
        
        sta_2 = end_1 + 1
        end_2 = sta_2 + nbli
        face_2 = lines[sta_2: end_2]
        
    count += 1

face_1 = [line.strip() for line in face_1]
face_2 = [line.strip() for line in face_2]

data_1 = [line.split() for line in face_1]
data_2 = [line.split() for line in face_2]

columns = ['NUMBER', 'GRID', 'ISTA', 'JSTA', 'KSTA', 'IEND', 'JEND', 'KEND', 'ISVA1', 'ISVA2']

df_1 = pd.DataFrame(data_1, columns=columns)
df_2 = pd.DataFrame(data_2, columns=columns)

# =============================================
# ================ fort.11 ====================
# =============================================
with open('input files/fort.11', 'w') as f:
    
    # ================== group 1 ==================
    f.write(f"TITLE: {title}\n")
    f.write(f'IDIM,\n{3:>4},\n')
    
    # ================== group 2 ==================
    group2 = ['IZON', 'IZFACE', 'IBND', 'ID', 'ISNGL']
    for g2 in group2:
        f.write(f'{g2:>6},')
    f.write('\n')
    
    IZFACE = df_1.loc[-1, 'NUMBER']
    IBND = len(IBCZON)
    ID = IZON*6 - IZFACE*2 - IBND
    ISNGL = 0
    
    f.write(f'{IZON:>6},{IZFACE:>6},{IBND:>6},{ID:>6},{ISNGL:>6},\n')
    
    # ================== group 3 ==================
    group3 = ['IZT', 'JZT', 'KZT', 'LPROC', 'CBG1', 'CBG2', 'CBG3', 'CBV1', 'CBV2', 'CBV3']

    for g3 in group3:
        f.write(f'{g3:>6},')
    f.write('\n')

    processor = 0
    for blk3 in range(IZON):
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
    
    def pair(df, blk, isva):
        if isva == '1':
            return df.loc[blk, 'ISTA'], df.loc[blk, 'IEND']
        elif isva == '2':
            return df.loc[blk, 'JSTA'], df.loc[blk, 'JEND']
        elif isva == '3':
            return df.loc[blk, 'KSTA'], df.loc[blk, 'KEND']
    
    def direction(df, blk, isva_1, isva_2):
        if isva_1 == '1' and isva_2 == '2':
            if df.loc[blk, 'KSTA'] == '1':
                return '6'
            else:
                return '5'
        elif isva_1 == '2' and isva_2 == '3':
            if df.loc[blk, 'ISTA'] == '1':
                return '2'
            else:
                return '1'
        elif isva_1 == '1' and isva_2 == '3':
            if df.loc[blk, 'JSTA'] == '1':
                return '4'
            else:
                return '3'
        
    for i in range(nbli):
        line_1 = []
        line_2 = []
        
        ifcyc = df_1.loc[i, 'NUMBER']
        izb_1 = df_1.loc[i, 'GRID']
        izb_2 = df_2.loc[i, 'GRID']
        inonuf = '0'
        
        isva_11 = df_1.loc[i, 'ISVA1']
        isva_12 = df_1.loc[i, 'ISVA2']
        isva_21 = df_2.loc[i, 'ISVA1']
        isva_22 = df_2.loc[i, 'ISVA2']
        
        ijz_11, ijz_12 = pair(df_1, i, isva_11)
        jkz_11, jkz_12 = pair(df_1, i, isva_12)
        ijz_21, ijz_22 = pair(df_2, i, isva_21)
        jkz_21, jkz_22 = pair(df_2, i, isva_22)
        
        izf_1 = direction(df_1, i, isva_11, isva_12)
        izf_2 = direction(df_2, i, isva_21, isva_22)
        
        line_1.append(ifcyc)
        line_1.append(izb_1)
        line_1.append(izf_1)
        line_1.append(ijz_11)
        line_1.append(ijz_12)
        line_1.append(jkz_11)
        line_1.append(jkz_12)
        line_1.append(inonuf)
        
        line_2.append(izb_2)
        line_2.append(izf_2)
        line_2.append(ijz_21)
        line_2.append(ijz_22)
        line_2.append(jkz_21)
        line_2.append(jkz_22)
        
        for item in line_1:
            f.write(f'{item:>6},')
        f.write('\n')
        f.write(f'{' ':7}')
        for item in line_2:
            f.write(f'{item:>6},')
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