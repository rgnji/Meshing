import pandas as pd

processor_num = 3
# ===== group 1 =====
title = 'INJECTOR A'
IDIM = 3
# ===== group 4 =====
cyclic_column = ['grid1','segment1','face1','grid2','segment2','face2']
cyclic = [[1,1,5,1,1,6],
          [1,2,5,1,2,6],
          [1,3,5,1,3,6],
          [1,4,5,1,4,6],
          [2,1,5,2,1,6],
          [2,2,5,2,2,6],
          [3,1,5,3,1,6],
          [3,2,5,3,2,6]]
# ===== group 6 =====
IWTM = 1
HQDOX = 0
IWALL = 0
DENNX = 1.
VISWX = 1.
# ===== group 8 =====
IDATA = 1
IGEO = 9    # any number > 0 but 1 and 9, which are for special cases
ITT = 1E5
ITPNT = 100
ICOUP = 1   # for unstable flow, ICOUP > 1
NLIMT = 1
IAX = 1
ICYC = 3
# ===== group 9 =====
# reference time = XREF/UREF = 1.4E-5
DTT = 1E-5
IREC = 0
REC = 0.25
THETA = 1   # time marching scheme
BETAP = 1   # compressible for real fluid model
IEXX = 1 #
PRAT = 0 #
# ===== group 10 =====
IPC = 50625
JPC = 1
IPEX = 50705
JPEX = 1
IMN = 50705
JMN = 1
# ===== group 11 =====
VISC = 18.37e-6
IG = 2      # laminar or turbulent
ITURB = 1   # use low-Re model if calculating viscous boundary layer
AMC = 1     # compressible for real fluid model
GAMA = 1.455
CBE = 0
CBH = 0 #
EREXT = 1E-8
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
# ============== read cfl3d file ==============
# =============================================
with open('input files/injector A o grid.inp', 'r') as f:
    # skip
    for _ in range(15):
        f.readline()
    
    # read izon
    sw = 1
    while(sw):
        line = f.readline()
        a = line.split()
        if a[0] == 'NGRID':
            sw = 0
    
    izon = int(f.readline().split()[1])
    
    # skip
    for _ in range(izon + 1):
        f.readline()
    
    # read dimension
    col_dim = f.readline().strip().split()
    dim = []
    for _ in range(izon):
        dim.append(f.readline().strip().split())
    
    df_dim = pd.DataFrame(dim, columns=col_dim)
    
    # skip
    for _ in range(5*(izon+1)):
        f.readline()
    
    # read boundary
    bnd = []
    col_bnd = []
    bnd_head = ['I0:', 'IDIM:', 'J0:', 'JDIM:', 'K0:', 'KDIM:']
    
    line = f.readline().strip().split()
    col_bnd.append(line[1:])
    
    for i in range(6):
        if i == 5:
            end = 'MSEQ'
        else:
            end = bnd_head[i+1]
        
        sw = 1
        bnd_blk = []
        while(sw):
            
            line = f.readline().strip().split()
            
            if line[0] == end:
                sw = 0
                if i != 5:
                    col_bnd.append(line[1:])
            else:
                bnd_blk.append(line)
                if line[-1] != '0':
                    for _ in range(2):
                        f.readline()
        
        bnd.append(bnd_blk)
    
    df_bnd = []
    df_bnd.append(pd.DataFrame(bnd[1], columns=col_bnd[1])) #idim
    df_bnd.append(pd.DataFrame(bnd[0], columns=col_bnd[0])) #i0
    df_bnd.append(pd.DataFrame(bnd[3], columns=col_bnd[3])) #jdim
    df_bnd.append(pd.DataFrame(bnd[2], columns=col_bnd[2])) #j0
    df_bnd.append(pd.DataFrame(bnd[5], columns=col_bnd[5])) #kdim
    df_bnd.append(pd.DataFrame(bnd[4], columns=col_bnd[4])) #k0
    
    # skip
    sw = 1
    while(sw):
        if f.readline() == '   1-1 BLOCKING DATA:\n':
            sw = 0
    
    # read patch interface
    f.readline()
    izface = int(f.readline().strip().split()[0])
    
    col_patch = f.readline().strip().split()
    col_patch.remove(':')
    
    face_1 = []
    for _ in range(izface):
        line = f.readline().strip().split()
        face_1.append(line)
    
    f.readline()
    face_2 = []
    for _ in range(izface):
        line = f.readline().strip().split()
        face_2.append(line)
    
    df_1 = pd.DataFrame(face_1, columns=col_patch)
    df_2 = pd.DataFrame(face_2, columns=col_patch)

# =============================================
# ================ fort.11 ====================
# =============================================
filename = 'input files/fort.11'
with open(filename, 'w') as f:
    
    # ================== group 1 ==================
    f.write(f"TITLE: {title}\n")
    f.write(f'IDIM,\n{3:>4},\n')
    
    # ================== group 2 ==================
    group2 = ['IZON', 'IZFACE', 'IBND', 'ID', 'ISNGL']
    for g2 in group2:
        f.write(f'{g2:>6},')
    f.write('\n')
    
    """
    1001 -> general symmetry plane    -> symmetry
    1003 -> inflow/outflow            -> mass outlet
    2003 -> inflow w/ total condition -> air inlet
    2004 -> no-slip wall              -> wall
    2005 -> periodic in space         -> cyclic
    2009 -> inflow                    -> h2o inlet
    """
    
    IZON = izon
    IZFACE = int(izface) + len(cyclic)
    IBND = 0
    ID = 0
    ISNGL = 0
    
    for face2 in range(6):
        for i in range(len(df_bnd[face2].index)):
            if df_bnd[face2].loc[i, 'BCTYPE'] in ['1003', '2003', '2009']:
                IBND += 1
            elif df_bnd[face2].loc[i, 'BCTYPE'] == '1001':
                IBND += 1
                ISNGL += 1
            elif df_bnd[face2].loc[i, 'BCTYPE'] == '2004':
                ID += 1
    
    f.write(f'{IZON:>6},{IZFACE:>6},{IBND:>6},{ID:>6},{ISNGL:>6},\n')
    
    # ================== group 3 ==================
    group3 = ['IZT', 'JZT', 'KZT', 'LPROC', 'CBG1', 'CBG2', 'CBG3', 'CBV1', 'CBV2', 'CBV3']

    for g3 in group3:
        f.write(f'{g3:>6},')
    f.write('\n')

    processor = 0
    for blk3 in range(IZON):
        f.write(f'{df_dim.loc[blk3, "IDIM"]:>6},{df_dim.loc[blk3, "JDIM"]:>6},{df_dim.loc[blk3, "KDIM"]:>6},')

        f.write(f'{(processor % processor_num) + 1:>6},')
        processor += 1

        for g3 in range(6):
            f.write(f'{0.:>6},')

        f.write('\n')
    
    # ================== group 4 ==================
    group4_1 = ['IFCYC', 'IZB1', 'IZF1', 'IJZ11', 'IJZ12', 'JKZ11', 'JKZ12', 'INONUF']
    group4_2 = ['IZB2', 'IZF2', 'IJZ21', 'IJZ22', 'JKZ21', 'JKZ22']

    for g4 in group4_1:
        f.write(f'{g4:>6},')
    f.write('\n')
    f.write(f'{" ":>7} ')
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
        
    for i in range(izface):
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
        f.write(f'{" ":7}')
        for item in line_2:
            f.write(f'{item:>6},')
        f.write('\n')
    
    # cyclic boundary
    df_cyc = pd.DataFrame(cyclic, columns=cyclic_column)
    
    ifcyc = '-1'
    inonuf = '0'
    for i4 in range(len(cyclic)):
        line_1 = []
        line_2 = []
        
        b_1 = df_cyc.loc[i4, 'grid1']
        b_2 = df_cyc.loc[i4, 'grid2']
        s_1 = df_cyc.loc[i4, 'segment1']
        s_2 = df_cyc.loc[i4, 'segment2']
        f_1 = df_cyc.loc[i4, 'face1']
        f_2 = df_cyc.loc[i4, 'face2']
        
        for t1 in range(len(df_bnd[f_1-1].index)):
            df = df_bnd[f_1-1]
            if df.loc[t1, 'GRID'] == str(b_1) and df.loc[t1, 'SEGMENT'] == str(s_1):
                t1_col = df.columns.tolist()
                ijz_11 = df.loc[t1, t1_col[3]]
                ijz_12 = df.loc[t1, t1_col[4]]
                jkz_11 = df.loc[t1, t1_col[5]]
                jkz_12 = df.loc[t1, t1_col[6]]
        
        for t2 in range(len(df_bnd[f_2-1].index)):
            df = df_bnd[f_2-1]
            if df.loc[t2, 'GRID'] == str(b_2) and df.loc[t2, 'SEGMENT'] == str(s_2):
                t2_col = df.columns.tolist()
                ijz_21 = df.loc[t2, t2_col[3]]
                ijz_22 = df.loc[t2, t2_col[4]]
                jkz_21 = df.loc[t2, t2_col[5]]
                jkz_22 = df.loc[t2, t2_col[6]]
        
        izb_1 = str(b_1)
        izb_2 = str(b_2)
        izf_1 = str(f_1)
        izf_2 = str(f_2)
        
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
        f.write(f'{" ":7}')
        for item in line_2:
            f.write(f'{item:>6},')
        f.write('\n')
    
    # ================== group 5 ==================
    group5 = ['IBCZON','IDBC','ITYBC','IJBB','IJBS','IJBT','JKBS','JKBT']

    for g5 in group5:
        f.write(f'{g5:>6},')
    f.write('\n')
    
    """
    1001 -> general symmetry plane    -> symmetry
    1003 -> inflow/outflow            -> mass outlet
    2003 -> inflow w/ total condition -> air inlet
    2004 -> no-slip wall              -> wall
    2005 -> periodic in space         -> cyclic
    2009 -> inflow                    -> h2o inlet
    """
    
    def find_g5(block, face):
        if face == 1:
            ijbb = df_dim.loc[block-1, 'IDIM']
        elif face == 3:
            ijbb = df_dim.loc[block-1, 'JDIM']
        elif face == 5:
            ijbb = df_dim.loc[block-1, 'KDIM']
        else:
            ijbb = 1
        return ijbb
    
    for face5 in range(6):
        
        for i in range(len(df_bnd[face5].index)):
            
            if df_bnd[face5].loc[i, 'BCTYPE'] in ['1001', '1003', '2003', '2009']:
                
                ibczon = df_bnd[face5].loc[i, 'GRID']
                idbc = face5 + 1
                
                ijbb = find_g5(int(ibczon), idbc)
                
                i_col = df_bnd[face5].columns.tolist()
                ijbs = df_bnd[face5].loc[i, i_col[3]]
                ijbt = df_bnd[face5].loc[i, i_col[4]]
                jkbs = df_bnd[face5].loc[i, i_col[5]]
                jkbt = df_bnd[face5].loc[i, i_col[6]]
            
                if df_bnd[face5].loc[i, 'BCTYPE'] == '1003':
                    itybc = '2'
                    
                elif df_bnd[face5].loc[i, 'BCTYPE'] in ['2003', '2009']:
                    itybc = '-1'
                
                elif df_bnd[face5].loc[i, 'BCTYPE'] == '1001':
                    itybc = '3'
                    
                f.write(f'{ibczon:>6},{idbc:>6},{itybc:>6},')
                f.write(f'{ijbb:>6},{ijbs:>6},{ijbt:>6},{jkbs:>6},{jkbt:>6},\n')
    
    # ================== group 6 ==================
    group6 = ['IWBZON','L1','L2','M1','M2','N1','N2','IWTM','HQDOX','IWALL','DENNX','VISWX']
    for g6 in group6:
        f.write(f'{g6:>6},')
    f.write('\n')
    
    def find_g6(df, i, face):
        
        if face == 0: #idim
            ijbb = df.loc[i, 'IDIM']
            ijbt = df.loc[i, 'JDIM']
            jkbt = df.loc[i, 'KDIM']
        elif face == 1: #i0
            ijbb = '1'
            ijbt = df.loc[i, 'JDIM']
            jkbt = df.loc[i, 'KDIM']
        elif face == 2: #jdim
            ijbb = df.loc[i, 'JDIM']
            ijbt = df.loc[i, 'IDIM']
            jkbt = df.loc[i, 'KDIM']
        elif face == 3: #j0
            ijbb = '1'
            ijbt = df.loc[i, 'IDIM']
            jkbt = df.loc[i, 'KDIM']
        elif face == 4: #kdim
            ijbb = df.loc[i, 'KDIM']
            ijbt = df.loc[i, 'IDIM']
            jkbt = df.loc[i, 'JDIM']
        elif face == 5: #k0
            ijbb = '1'
            ijbt = df.loc[i, 'IDIM']
            jkbt = df.loc[i, 'JDIM']
        
        return ijbb, ijbt, jkbt
    
    for face6 in range(6):
        
        for i in range(len(df_bnd[face6].index)):
            
            if df_bnd[face6].loc[i, 'BCTYPE'] == '2004':
                
                iwbzon = df_bnd[face6].loc[i, 'GRID']
                iwbzon = int(iwbzon)
                
                if face6 in [0, 1]:
                    l1, m2, n2 = find_g6(df_dim, iwbzon-1, face6)
                    l2 = l1
                    m1 = 1
                    n1 = 1
                elif face6 in [2, 3]:
                    m1, l2, n2 = find_g6(df_dim, iwbzon-1, face6)
                    m2 = m1
                    l1 = 1
                    n1 = 1
                elif face6 in [4, 5]:
                    n1, l2, m2 = find_g6(df_dim, iwbzon-1, face6)
                    n2 = n1
                    l1 = 1
                    m1 = 1
                
                f.write(f'{iwbzon:>6},{l1:>6},{l2:>6},{m1:>6},{m2:>6},{n1:>6},{n2:>6},')
                f.write(f'{IWTM:>6},{HQDOX:>6},{IWALL:>6},{DENNX:>6},{VISWX:>6},\n')
    
    # ================== group 7 ==================
    group7 = ['ISNZON', 'ISNBC', 'ISNAX', 'ISNBS', 'ISNBT']
    for g7 in group7:
        f.write(f'{g7:>6},')
    f.write('\n')
    
    for face7 in range(6):
        for i in range(len(df_bnd[face7].index)):
            df = df_bnd[face7]
            if df.loc[i, 'BCTYPE'] == '1001':
                isnzon = df.loc[i, 'GRID']
                isnbc = face7 + 1
                # assume singular line is along i-axis (same as cyclic boundary)
                isnax = 1
                isnbs = df.loc[i, 'ISTA']
                isnbt = df.loc[i, 'IEND']
                
                f.write(f'{isnzon:>6},{isnbc:>6},{isnax:>6},{isnbs:>6},{isnbt:>6},')
                f.write('\n')
    
    # ================== group 8 ==================
    group8 = ['IDATA', 'IGEO', 'ITT', 'ITPNT', 'ICOUP', 'NLIMT', 'IAX', 'ICYC']
    for g8 in group8:
        f.write(f'{g8:>9},')
    f.write('\n')
    f.write(f'{IDATA:>9},{IGEO:>9},{ITT:>.3E},{ITPNT:>9},{ICOUP:>9},{NLIMT:>9},{IAX:>9},{ICYC:>9},\n')

    # ================== group 9 ==================
    group9 = ['DTT', 'IREC', 'REC', 'THETA', 'BETAP', 'IEXX', 'PRAT']
    for g9 in group9:
        f.write(f'{g9:>9},')
    f.write('\n')
    f.write(f'{DTT:>.3E},{IREC:>9},{REC:>9},{THETA:>9},{BETAP:>9},{IEXX:>9},{PRAT:>9},\n')

    # ================== group 10 ==================
    group10 = ['IPC', 'JPC', 'IPEX', 'JPEX', 'IMN', 'JMN']
    for g10 in group10:
        f.write(f'{g10:>9},')
    f.write('\n')
    f.write(f'{IPC:>9},{JPC:>9},{IPEX:>9},{JPEX:>9},{IMN:>9},{JMN:>9},\n')

    # ================== group 11 ==================
    group11 = ['VISC', 'IG', 'ITURB', 'AMC', 'GAMA', 'CBE', 'CBH', 'EREXT']
    for g11 in group11:
        f.write(f'{g11:>9},')
    f.write('\n')
    f.write(f'{VISC:>9},{IG:>9},{ITURB:>9},{AMC:>9},{GAMA:>9},{CBE:>9},{CBH:>9},{EREXT:>9},\n')

    # ================== group 12 ==================
    group12 = ['ISWU','ISWP','ISWK','ISKEW']
    for g12 in group12:
        f.write(f'{g12:>9},')
    f.write('\n')
    f.write(f'{ISWU:>9},{ISWP:>9},{ISWK:>9},{ISKEW:>9},\n')

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
    
    f.write('FLUID\n')
    f.write('c Species(a20) , ideal gas(=0), real fluid(=1)\n')
    f.write(f"{'H2O':<20}{'1':>5}\n")
    f.write(f"{'AIR':<20}{'0':>5}\n")
    f.write(f"{'DONE':<20}{'0':>5}\n")
    
    # ================== end ==================
    print(filename + ' established.')