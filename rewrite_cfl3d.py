from struct import unpack
import numpy as np
from lib.plot3dout import unformatted_fort12

with open('input files/injector A o grid.bin', 'rb') as f:
    buff = f.read(4*3)
    izon = unpack('<3i', buff)[1]
    
    dim = []
    f.read(4)
    for i in range(izon):
        buff = f.read(4*3)
        dim.append(unpack('<3i', buff))
    f.read(4)
    
    X, Y, Z = [], [], []
    for i in range(izon):
        izt = dim[i][0]
        jzt = dim[i][1]
        kzt = dim[i][2]
        size = izt * jzt * kzt
        
        f.read(4)
        buff = f.read(8 * size)
        x = np.array(unpack(f'{size}d', buff)).reshape((kzt, jzt, izt))
        x = np.float32(x)
        X.append(x)
        
        buff = f.read(8 * size)
        y = np.array(unpack(f'{size}d', buff)).reshape((kzt, jzt, izt))
        y = np.float32(y)
        Y.append(y)
        
        buff = f.read(8 * size)
        z = np.array(unpack(f'{size}d', buff)).reshape((kzt, jzt, izt))
        z = np.float32(z)
        Z.append(z)
        f.read(4)
    
txt = unformatted_fort12(X, Y, Z, 'input files/')
print(txt)
        