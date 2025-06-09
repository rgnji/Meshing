import numpy as np
import struct

#=================== unformatted plot3d format (with record marker) ===================
# X, Y, Z will be 4D matrices X(IZ, KZT, JZT, IZT)
def unformatted_fort12(X, Y, Z):
    filename = "fort.12"

    with open(filename, "wb") as f:
        # write number of blocks
        blocks = np.array(len(X), dtype=np.int32) # len(X) = number of blocks

        f.write(struct.pack("<i", 4)) # start marker
        blocks.tofile(f)
        f.write(struct.pack("<i", 4)) # end marker

        # write maximum grid number in I, J and K direction for each block
        for i in range(len(X)):
            dim = []
            for j in range(3):
                dim.append(X[i].shape[2-j]) # IZT, JZT, KZT
            
            f.write(struct.pack('<i', 12))
            dim = np.array(dim, dtype=np.int32)
            dim.tofile(f)
            f.write(struct.pack('<i', 12))

        # write coordinates of all the points in a block
        # write all the "x", then all the "y", then "z"
        for block in range(len(X)):
            x_flat = np.array(X[block].flatten(), dtype=np.float32) # total element size
            len_x = x_flat.nbytes

            f.write(struct.pack("<i", len_x)) # write x 
            x = np.array(X[block], dtype=np.float32) # change dtype
            x.tofile(f) # row-major
            f.write(struct.pack("<i", len_x))

            f.write(struct.pack("<i", len_x)) # write y
            y = np.array(Y[block], dtype=np.float32)
            y.tofile(f)
            f.write(struct.pack("<i", len_x))

            f.write(struct.pack("<i", len_x)) # write z
            z = np.array(Z[block], dtype=np.float32)
            z.tofile(f)
            f.write(struct.pack("<i", len_x))
    
    return filename + " established."


#=================== binary plot3d format (no record marker) ===================
# X, Y, Z will be 4D matrices X(IZ, KZT, JZT, IZT)
def binary_fort12(X, Y, Z):
    filename = "fort12.bin.xyz"

    with open(filename, "wb") as f:
        # write number of blocks
        blocks = np.array(len(X), dtype=np.int32) # len(X) = number of blocks
        f.write(struct.pack("<i", blocks))

        # write maximum grid number in I, J and K direction for each block
        dim = []
        for i in range(len(X)):
            for j in range(3):
                dim.append(X[i].shape[2-j]) # IZT, JZT, KZT
        dim = np.array(dim, dtype=np.int32)
        dim.tofile(f)

        # write coordinates of all the points in a block
        # write all the "x", then all the "y", then "z"
        for block in range(len(X)):
            x = np.array(X[block], dtype=np.float64) # change dtype
            x.tofile(f) # row-major

            y = np.array(Y[block], dtype=np.float64)
            y.tofile(f)

            z = np.array(Z[block], dtype=np.float64)
            z.tofile(f)
    
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
            dk[i].tofile(f)
            de[i].tofile(f)
            am[i].tofile(f)
            q[i].tofile(f)
            for kk in range(len(fm[i])):
                fm[i][kk].tofile(f)
    
    return filename + ' established.'