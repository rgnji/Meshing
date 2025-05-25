import numpy as np
import struct

#=================== unformatted plot3d format (with record marker) ===================
# X, Y, Z will be 4D matrices X(IZ, KZT, JZT, IZT)
def unformatted_fort12(X, Y, Z):
    filename = "fort.12"

    with open(filename, "wb") as f:
        # write number of blocks
        blocks = np.array(len(X), dtype=np.int32) # len(X) = number of blocks
        len_blocks = blocks.nbytes

        f.write(struct.pack("<i", len_blocks)) # start marker
        f.write(struct.pack("<i", blocks))
        f.write(struct.pack("<i", len_blocks)) # end marker

        # write maximum grid number in I, J and K direction for each block
        dim = []
        for i in range(len(X)):
            for j in range(3):
                dim.append(X[i].shape[2-j]) # IZT, JZT, KZT
        dim = np.array(dim, dtype=np.int32)
        len_dim = dim.nbytes

        f.write(struct.pack("<i", len_dim))
        dim.tofile(f)
        f.write(struct.pack("<i", len_dim))

        # write coordinates of all the points in a block
        # write all the "x", then all the "y", then "z"
        for block in range(len(X)):
            x_flat = np.array(X[block].flatten(), dtype=np.float64) # total element size
            len_x = x_flat.nbytes

            f.write(struct.pack("<i", len_x)) # write x 
            x = np.array(X[block], dtype=np.float64) # change dtype
            x.tofile(f) # row-major
            f.write(struct.pack("<i", len_x))

            f.write(struct.pack("<i", len_x)) # write y
            y = np.array(Y[block], dtype=np.float64)
            y.tofile(f)
            f.write(struct.pack("<i", len_x))

            f.write(struct.pack("<i", len_x)) # write z
            z = np.array(Z[block], dtype=np.float64)
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
def unformatted_fort13(INSO_1, INSO_4, INSO_5, INSO_7, NGAS, IZON, den, u, v, w, p, dk, de, q, fm):
    filename = 'fort.13'

    with open(filename, 'wb') as f:
        len1 = np.array([INSO_1], dtype=np.int32).nbytes
        len2 = np.array([INSO_4], dtype=np.int32).nbytes
        len3 = np.array([INSO_5], dtype=np.int32).nbytes
        len4 = np.array([INSO_7], dtype=np.int32).nbytes
        len5 = np.array([NGAS], dtype=np.int32).nbytes

        f.write(struct.pack('<i', len1))
        f.write(struct.pack('<i', INSO_1))
        f.write(struct.pack('<i', len1))

        f.write(struct.pack('<i', len2))
        f.write(struct.pack('<i', INSO_4))
        f.write(struct.pack('<i', len2))

        f.write(struct.pack('<i', len3))
        f.write(struct.pack('<i', INSO_5))
        f.write(struct.pack('<i', len3))

        f.write(struct.pack('<i', len4))
        f.write(struct.pack('<i', INSO_7))
        f.write(struct.pack('<i', len4))

        f.write(struct.pack('<i', len5))
        f.write(struct.pack('<i', NGAS))
        f.write(struct.pack('<i', len5))

        for i in range(IZON):
            len_den = np.array(den[i], dtype=np.float64).nbytes
            f.write(struct.pack('<i', len_den))
            den[i].astype(np.float64).tofile(f)
            f.write(struct.pack('<i', len_den))

            len_u = np.array(u[i], dtype=np.float64).nbytes
            f.write(struct.pack('<i', len_u))
            u[i].astype(np.float64).tofile(f)
            f.write(struct.pack('<i', len_u))

            len_v = np.array(v[i], dtype=np.float64).nbytes
            f.write(struct.pack('<i', len_v))
            v[i].astype(np.float64).tofile(f)
            f.write(struct.pack('<i', len_v))

            len_w = np.array(w[i], dtype=np.float64).nbytes
            f.write(struct.pack('<i', len_w))
            w[i].astype(np.float64).tofile(f)
            f.write(struct.pack('<i', len_w))

            len_p = np.array(p[i], dtype=np.float64).nbytes
            f.write(struct.pack('<i', len_p))
            p[i].astype(np.float64).tofile(f)
            f.write(struct.pack('<i', len_p))

            len_dk = np.array(dk[i], dtype=np.float64).nbytes
            f.write(struct.pack('<i', len_dk))
            dk[i].astype(np.float64).tofile(f)
            f.write(struct.pack('<i', len_dk))

            len_de = np.array(de[i], dtype=np.float64).nbytes
            f.write(struct.pack('<i', len_de))
            de[i].astype(np.float64).tofile(f)
            f.write(struct.pack('<i', len_de))

            len_q = np.array(q[i], dtype=np.float64).nbytes
            f.write(struct.pack('<i', len_q))
            q[i].astype(np.float64).tofile(f)
            f.write(struct.pack('<i', len_q))

            # the same order what cec table in fort.11 is
            for kk in range(NGAS):
                len_fm = np.array(fm[i][kk], dtype=np.float64).nbytes
                f.write(struct.pack('<i', len_fm))
                fm[i][kk].astype(np.float64).tofile(f)
                f.write(struct.pack('<i', len_fm))
    
    return filename + ' established.'


#=================== binary plot3d format (no record marker) ===================
def binary_fort13(INSO_1, INSO_4, INSO_5, INSO_7, NGAS, IZON, den, u, v, w, p, dk, de, q, fm):
    filename = 'fort13.bin.xyz'

    with open(filename, 'wb') as f:
        f.write(struct.pack('<i', INSO_1))
        f.write(struct.pack('<i', INSO_4))
        f.write(struct.pack('<i', INSO_5))
        f.write(struct.pack('<i', INSO_7))
        f.write(struct.pack('<i', NGAS))

        for i in range(IZON):
            den[i].astype(np.float64).tofile(f)
            u[i].astype(np.float64).tofile(f)
            v[i].astype(np.float64).tofile(f)
            w[i].astype(np.float64).tofile(f)
            p[i].astype(np.float64).tofile(f)
            dk[i].astype(np.float64).tofile(f)
            de[i].astype(np.float64).tofile(f)
            q[i].astype(np.float64).tofile(f)
            # the same order what cec table in fort.11 is
            for kk in range(NGAS):
                fm[i][kk].astype(np.float64).tofile(f)
    
    return filename + ' established.'