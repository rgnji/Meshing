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
