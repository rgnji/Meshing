import numpy as np
from struct import unpack

#
#  read fort.12
#
with open('fort.12', 'rb') as f:
    buff = f.read(4*3)
    IZON = unpack('<3i', buff)[1]
    
    IZT, JZT, KZT = [], [], []
    for i in range(IZON):
        buff = f.read(4*5)
        arr = unpack('<5i', buff)
        IZT.append(arr[1])
        JZT.append(arr[2])
        KZT.append(arr[3])
    
    XT, YT, ZT = [], [], []
    for i in range(IZON):
        size = IZT[i] * JZT[i] * KZT[i]
        buff = f.read(4*(size+2))
        grid = np.array(unpack(f'<{size+2}f', buff)[1:-1])
        XT.append(grid)
        buff = f.read(4*(size+2))
        grid = np.array(unpack(f'<{size+2}f', buff)[1:-1])
        YT.append(grid)
        buff = f.read(4*(size+2))
        grid = np.array(unpack(f'<{size+2}f', buff)[1:-1])
        ZT.append(grid)

#
#  read fort.13
#
with open('fort.13', 'rb') as f:
    buff = f.read(4*7)
    INSO_1, INSO_4, INSO_5, INSO_7, NGAS = unpack('<7i', buff)[1:6]
    
    BLKDN, BLKU, BLKV, BLKW, BLKP = [], [], [], [], []
    BLKDK, BLKDE, BLKQ, BLKAM, BLKFM = [], [], [], [], []
    for i in range(IZON):
        size = IZT[i] * JZT[i] * KZT[i]
        
        buff = f.read(4*size)
        grid = np.array(unpack(f'<{size}f', buff))
        BLKDN.append(grid)
        
        buff = f.read(4*size)
        grid = np.array(unpack(f'<{size}f', buff))
        BLKU.append(grid)
        
        buff = f.read(4*size)
        grid = np.array(unpack(f'<{size}f', buff))
        BLKV.append(grid)
        
        buff = f.read(4*size)
        grid = np.array(unpack(f'<{size}f', buff))
        BLKW.append(grid)
        
        buff = f.read(4*size)
        grid = np.array(unpack(f'<{size}f', buff))
        BLKP.append(grid)
        
        buff = f.read(4*size)
        grid = np.array(unpack(f'<{size}f', buff))
        BLKDK.append(grid)
        
        buff = f.read(4*size)
        grid = np.array(unpack(f'<{size}f', buff))
        BLKDE.append(grid)
        
        buff = f.read(4*size)
        grid = np.array(unpack(f'<{size}f', buff))
        BLKQ.append(grid)
        
        buff = f.read(4*size)
        grid = np.array(unpack(f'<{size}f', buff))
        BLKAM.append(grid)
        
        buff = f.read(4*size)
        grid1 = np.array(unpack(f'<{size}f', buff))
        buff = f.read(4*size)
        grid2 = np.array(unpack(f'<{size}f', buff))
        BLKFM.append([grid1, grid2])

#
#  tecplot output
#
with open('tecout.dat', 'w') as f:
    VAR1 = 'VARIABLES=X,Y,Z,U,V,W,PRES,DEN,MACH,QUAL'
    f.write('TITLE = PROPERTIES\n')
    f.write(f"{VAR1:<45}{'H2O':<10}{'AIR':<10}\n")
    
    for i in range(IZON):
        f.write(f"ZONE I={IZT[i]:>5} J={JZT[i]:>5} K={KZT[i]:>5} F=BLOCK\n")
        
        size = IZT[i] * JZT[i] * KZT[i]
        count = 0
        for j in range(size):
            f.write(f"  {XT[i][j]:13.6E}")
            count += 1
            if count % 5 == 0:
                f.write('\n')
        if count % 5 != 0:
            f.write('\n')
        
        count = 0
        for j in range(size):
            f.write(f"  {YT[i][j]:13.6E}")
            count += 1
            if count % 5 == 0:
                f.write('\n')
        if count % 5 != 0:
            f.write('\n')
            
        count = 0
        for j in range(size):
            f.write(f"  {ZT[i][j]:13.6E}")
            count += 1
            if count % 5 == 0:
                f.write('\n')
        if count % 5 != 0:
            f.write('\n')
        
        count = 0
        for j in range(size):
            f.write(f"  {BLKU[i][j]:13.6E}")
            count += 1
            if count % 5 == 0:
                f.write('\n')
        if count % 5 != 0:
            f.write('\n')
        
        count = 0
        for j in range(size):
            f.write(f"  {BLKV[i][j]:13.6E}")
            count += 1
            if count % 5 == 0:
                f.write('\n')
        if count % 5 != 0:
            f.write('\n')
            
        count = 0
        for j in range(size):
            f.write(f"  {BLKW[i][j]:13.6E}")
            count += 1
            if count % 5 == 0:
                f.write('\n')
        if count % 5 != 0:
            f.write('\n')
            
        count = 0
        for j in range(size):
            f.write(f"  {BLKP[i][j]:13.6E}")
            count += 1
            if count % 5 == 0:
                f.write('\n')
        if count % 5 != 0:
            f.write('\n')
            
        count = 0
        for j in range(size):
            f.write(f"  {BLKDN[i][j]:13.6E}")
            count += 1
            if count % 5 == 0:
                f.write('\n')
        if count % 5 != 0:
            f.write('\n')
            
        count = 0
        for j in range(size):
            f.write(f"  {BLKAM[i][j]:13.6E}")
            count += 1
            if count % 5 == 0:
                f.write('\n')
        if count % 5 != 0:
            f.write('\n')
            
        count = 0
        for j in range(size):
            f.write(f"  {BLKQ[i][j]:13.6E}")
            count += 1
            if count % 5 == 0:
                f.write('\n')
        if count % 5 != 0:
            f.write('\n')
        
        for kk in range(NGAS):
            count = 0
            for j in range(size):
                f.write(f"  {BLKFM[i][kk][j]:13.6E}")
                count += 1
                if count % 5 == 0:
                    f.write('\n')
            if count % 5 != 0:
                f.write('\n')

print('tecout.dat established.')