import numpy as np
import struct

# X, Y, Z will be 4D matrices X(IZ, KZT, JZT, IZT)
def bin_fort12(X, Y, Z):
    with open("fort.12", "wb") as f:
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
    
    return "fort.12 established."
