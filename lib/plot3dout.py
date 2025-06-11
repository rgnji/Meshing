import numpy as np
import struct

#=================== unformatted plot3d format (with record marker) ===================
# X, Y, Z will be 4D matrices X(IZ, KZT, JZT, IZT)
def unformatted_fort12(X, Y, Z):
    filename = 'fort.12'
    decimals = 5
    
    IZON = len(X)
    IZT, JZT, KZT = [], [], []
    for i in range(IZON):
        IZT.append(X[i].shape[2])
        JZT.append(X[i].shape[1])
        KZT.append(X[i].shape[0])
    
    with open(filename, 'wb') as f:
        f.write(struct.pack('<i', 4))
        f.write(struct.pack('<i', IZON))
        f.write(struct.pack('<i', 4))
        
        for i in range(IZON):
            f.write(struct.pack('<i', 12))
            f.write(struct.pack('<i', IZT[i]))
            f.write(struct.pack('<i', JZT[i]))
            f.write(struct.pack('<i', KZT[i]))
            f.write(struct.pack('<i', 12))
        
        for i in range(IZON):
            size = IZT[i]*JZT[i]*KZT[i]
            
            XT = X[i].flatten()
            f.write(struct.pack('<i', 4*size))
            for j in range(size):    
                f.write(struct.pack('<f', round(XT[j], decimals)))    
            f.write(struct.pack('<i', 4*size))
            
            YT = Y[i].flatten()
            f.write(struct.pack('<i', 4*size))
            for j in range(size):    
                f.write(struct.pack('<f', round(YT[j], decimals)))    
            f.write(struct.pack('<i', 4*size))
            
            ZT = Z[i].flatten()
            f.write(struct.pack('<i', 4*size))
            for j in range(size):    
                f.write(struct.pack('<f', round(ZT[j], decimals)))    
            f.write(struct.pack('<i', 4*size))
    
    return filename + " established."

#=================== unformatted plot3d format (with record marker) ===================
def unformatted_fort13(INSO_1, INSO_4, INSO_5, INSO_7, NGAS, IZON, den, u, v, w, p, dk, de, am, q, fm):
    filename = 'fort.13'

    with open(filename, 'wb') as f:

        f.write(struct.pack('<i', 20))
        f.write(struct.pack('<i', INSO_1))
        f.write(struct.pack('<i', INSO_4))
        f.write(struct.pack('<i', INSO_5))
        f.write(struct.pack('<i', INSO_7))
        f.write(struct.pack('<i', NGAS))
        f.write(struct.pack('<i', 20))

        for i in range(IZON):
            len_den = np.array(den[i], dtype=np.float32).nbytes
            f.write(struct.pack('<i', len_den))
            den[i].astype(np.float32).tofile(f)
            f.write(struct.pack('<i', len_den))

            len_u = np.array(u[i], dtype=np.float32).nbytes
            f.write(struct.pack('<i', len_u))
            u[i].astype(np.float32).tofile(f)
            f.write(struct.pack('<i', len_u))

            len_v = np.array(v[i], dtype=np.float32).nbytes
            f.write(struct.pack('<i', len_v))
            v[i].astype(np.float32).tofile(f)
            f.write(struct.pack('<i', len_v))

            len_w = np.array(w[i], dtype=np.float32).nbytes
            f.write(struct.pack('<i', len_w))
            w[i].astype(np.float32).tofile(f)
            f.write(struct.pack('<i', len_w))

            len_p = np.array(p[i], dtype=np.float32).nbytes
            f.write(struct.pack('<i', len_p))
            p[i].astype(np.float32).tofile(f)
            f.write(struct.pack('<i', len_p))

            len_dk = np.array(dk[i], dtype=np.float32).nbytes
            f.write(struct.pack('<i', len_dk))
            dk[i].astype(np.float32).tofile(f)
            f.write(struct.pack('<i', len_dk))

            len_de = np.array(de[i], dtype=np.float32).nbytes
            f.write(struct.pack('<i', len_de))
            de[i].astype(np.float32).tofile(f)
            f.write(struct.pack('<i', len_de))

            len_am = np.array(am[i], dtype=np.float32).nbytes
            f.write(struct.pack('<i', len_am))
            am[i].astype(np.float32).tofile(f)
            f.write(struct.pack('<i', len_am))

            len_q = np.array(q[i], dtype=np.float32).nbytes
            f.write(struct.pack('<i', len_q))
            q[i].astype(np.float32).tofile(f)
            f.write(struct.pack('<i', len_q))

            # the same order what cec table in fort.11 is
            for kk in range(NGAS):
                len_fm = np.array(fm[i][kk], dtype=np.float32).nbytes
                f.write(struct.pack('<i', len_fm))
                fm[i][kk].astype(np.float32).tofile(f)
                f.write(struct.pack('<i', len_fm))
    
    return filename + ' established.'

#
#  ASCII fort.12
#
def ascii_fort12(X, Y, Z):
    filename = 'fort12.txt'
    decimals = 6
    IZON = len(X)
    IZT, JZT, KZT = [], [], []
    for i in range(IZON):
        IZT.append(X[i].shape[2])
        JZT.append(X[i].shape[1])
        KZT.append(X[i].shape[0])
    
    with open(filename, 'w') as f:
        f.write(f'{IZON:>5}\n')
        
        for i in range(IZON):
            f.write(f'{IZT[i]:>5}')
            f.write(f'{JZT[i]:>5}')
            f.write(f'{KZT[i]:>5}')
            f.write('\n')
        
        for i in range(IZON):
            size = IZT[i]*JZT[i]*KZT[i]
            XT = X[i].flatten() * 10
            YT = Y[i].flatten() * 10
            ZT = Z[i].flatten() * 10
            
            count = 0
            for j in range(size):
                f.write(f'{round(XT[j], decimals):>16.8E}')
                count += 1
                if count % 5 == 0:
                    f.write('\n')
            if count % 5 != 0:
                f.write('\n')
                
            count = 0
            for j in range(size):
                f.write(f'{round(YT[j], decimals):>16.8E}')
                count += 1
                if count % 5 == 0:
                    f.write('\n')
            if count % 5 != 0:
                f.write('\n')
            
            count = 0
            for j in range(size):
                f.write(f'{round(ZT[j], decimals):>16.8E}')
                count += 1
                if count % 5 == 0:
                    f.write('\n')
            if count % 5 != 0:
                f.write('\n')
            
    return filename + '(ASCII) established.'

#
#  ASCII fort.13
#
def ascii_fort13(INSO_1, INSO_4, INSO_5, INSO_7, NGAS, IZON, DN, U, V, W, P, DK, DE, AM, Q, FM):
    filename = 'fort13.txt'
    IZON = len(DN)
    IZT, JZT, KZT = [], [], []
    for i in range(IZON):
        IZT.append(DN[i].shape[2])
        JZT.append(DN[i].shape[1])
        KZT.append(DN[i].shape[0])
    
    with open(filename, 'w') as f:
        f.write(f'{INSO_1:>5}{INSO_4:>5}{INSO_5:>5}{INSO_7:>5}{NGAS:>5}\n')
        
        for i in range(IZON):
            size = IZT[i]*JZT[i]*KZT[i]
            
            DNT = DN[i].flatten() * 10
            count = 0
            for j in range(size):
                f.write(f'{DNT[j]:>16.8E}')
                count += 1
                if count % 5 == 0:
                    f.write('\n')
            if count % 5 != 0:
                f.write('\n')
                
            UT = U[i].flatten() * 10
            count = 0
            for j in range(size):
                f.write(f'{UT[j]:>16.8E}')
                count += 1
                if count % 5 == 0:
                    f.write('\n')
            if count % 5 != 0:
                f.write('\n')
                
            VT = V[i].flatten() * 10
            count = 0
            for j in range(size):
                f.write(f'{VT[j]:>16.8E}')
                count += 1
                if count % 5 == 0:
                    f.write('\n')
            if count % 5 != 0:
                f.write('\n')
            
            WT = W[i].flatten() * 10
            count = 0
            for j in range(size):
                f.write(f'{WT[j]:>16.8E}')
                count += 1
                if count % 5 == 0:
                    f.write('\n')
            if count % 5 != 0:
                f.write('\n')
                
            PT = P[i].flatten() * 10
            count = 0
            for j in range(size):
                f.write(f'{PT[j]:>16.8E}')
                count += 1
                if count % 5 == 0:
                    f.write('\n')
            if count % 5 != 0:
                f.write('\n')
            
            DKT = DK[i].flatten() * 10
            count = 0
            for j in range(size):
                f.write(f'{DKT[j]:>16.8E}')
                count += 1
                if count % 5 == 0:
                    f.write('\n')
            if count % 5 != 0:
                f.write('\n')
                
            DET = DE[i].flatten() * 10
            count = 0
            for j in range(size):
                f.write(f'{DET[j]:>16.8E}')
                count += 1
                if count % 5 == 0:
                    f.write('\n')
            if count % 5 != 0:
                f.write('\n')
            
            AMT = AM[i].flatten() * 10
            count = 0
            for j in range(size):
                f.write(f'{AMT[j]:>16.8E}')
                count += 1
                if count % 5 == 0:
                    f.write('\n')
            if count % 5 != 0:
                f.write('\n')
            
            QT = Q[i].flatten() * 10
            count = 0
            for j in range(size):
                f.write(f'{QT[j]:>16.8E}')
                count += 1
                if count % 5 == 0:
                    f.write('\n')
            if count % 5 != 0:
                f.write('\n')
                
            for kk in range(NGAS):
                FMT = FM[i][kk].flatten() * 10
                count = 0
                for j in range(size):
                    f.write(f'{FMT[j]:>16.8E}')
                    count += 1
                    if count % 5 == 0:
                        f.write('\n')
                if count % 5 != 0:
                    f.write('\n')
    
    return filename + '(ASCII) established.'

#=================== binary plot3d format (no record marker) ===================
# X, Y, Z will be 4D matrices X(IZ, KZT, JZT, IZT)
def binary_fort12(X, Y, Z):
    filename = "fort12.bin.xyz"
    decimals = 5
    
    IZON = len(X)
    IZT, JZT, KZT = [], [], []
    for i in range(IZON):
        IZT.append(X[i].shape[2])
        JZT.append(X[i].shape[1])
        KZT.append(X[i].shape[0])
    
    with open(filename, 'wb') as f:
        f.write(struct.pack('<i', IZON))
        
        for i in range(IZON):
            f.write(struct.pack('<i', IZT[i]))
            f.write(struct.pack('<i', JZT[i]))
            f.write(struct.pack('<i', KZT[i]))
        
        for i in range(IZON):
            size = IZT[i]*JZT[i]*KZT[i]
            
            XT = X[i].flatten()
            for j in range(size):    
                f.write(struct.pack('<f', round(XT[j], decimals)))    
            
            YT = Y[i].flatten()
            for j in range(size):    
                f.write(struct.pack('<f', round(YT[j], decimals)))    
            
            ZT = Z[i].flatten()
            for j in range(size):    
                f.write(struct.pack('<f', round(ZT[j], decimals)))    
    
    return filename + " established."

#=================== binary plot3d format (no record marker) ===================
def binary_fort13(INSO_1, INSO_4, INSO_5, INSO_7, NGAS, IZON, den, u, v, w, p, dk, de, am, q, fm):
    filename = 'fort13.bin.q'

    with open(filename, 'wb') as f:
        
        # nblocks
        f.write(struct.pack('<i', len(den)))
        
        for i in range(len(den)):
            ni = den[i].shape[2]
            nj = den[i].shape[1]
            nk = den[i].shape[0]
            f.write(struct.pack('<i', ni))
            f.write(struct.pack('<i', nj))
            f.write(struct.pack('<i', nk))
        
        for i in range(len(den)):
            den[i].tofile(f)
            u[i].tofile(f)
            v[i].tofile(f)
            w[i].tofile(f)
            p[i].tofile(f)
    
    return filename + ' established.'