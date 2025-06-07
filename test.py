import numpy as np
import struct

with open('fort.12', 'rb') as f:
    buff  = f.read(4*3)
    blocks = struct.unpack('<3i', buff)[1]
    print(blocks)

    dim = []
    for i in range(blocks):
        buff = f.read(4*5)
        dim.append(struct.unpack('<5i', buff)[1:4])
    print(len(dim))
    
    X, Y, Z = [], [], []
    for i in range(blocks):
        size = dim[i][0] * dim[i][1] * dim[i][2]
        print(size)
        for j in range(3):
            buff = f.read(4)
            buff = f.read(4*size)
            tmp = struct.unpack(f'<{size}f', buff)
            X.append(np.array(tmp, dtype=np.float32).reshape(dim[i][2], dim[i][1], dim[i][0]))
            buff = f.read(4)

            buff = f.read(4)
            buff = f.read(4*size)
            tmp = struct.unpack(f'<{size}f', buff)
            Y.append(np.array(tmp, dtype=np.float32).reshape(dim[i][2], dim[i][1], dim[i][0]))
            buff = f.read(4)

            buff = f.read(4)
            buff = f.read(4*size)
            tmp = struct.unpack(f'<{size}f', buff)
            Z.append(np.array(tmp, dtype=np.float32).reshape(dim[i][2], dim[i][1], dim[i][0]))
            buff = f.read(4)
    