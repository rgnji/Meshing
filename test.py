"""
# gas inlet
for qua in range(1,5): 
    x, y, z = o_grid(qua, 14.11, 0.5*14.11, 20, 0, 6, 15, 28, 10)
    XT.append(x)
    YT.append(y)
    ZT.append(z)
    izon_count += 1
    ibnd_count += 1
    id_count += 1
    IZT.append(x.shape[2])
    JZT.append(x.shape[1])
    KZT.append(x.shape[0])

x, y, z = h_grid(0.5*14.11, 20, 0, 6, 15, 10)
XT.append(x)
YT.append(y)
ZT.append(z)
izon_count += 1
ibnd_count += 1
IZT.append(x.shape[2])
JZT.append(x.shape[1])
KZT.append(x.shape[0])

# liquid inlet
xx, yy, zz = liquid_cylinder(r_orifice, R_liquidout, R_liquidin, num_s, num_r, num_rl, num_arc, num_pipe, inlet_pipe_len)

XT += xx
YT += yy
ZT += zz

for i in range(len(xx)):
    izon_count += 1
    IZT.append(xx[i].shape[2])
    JZT.append(xx[i].shape[1])
    KZT.append(xx[i].shape[0])
ibnd_count += 10
id_count += (4 + 10) * 2

# liquid cylinder
for qua in range(1,5): 
    x, y, z = o_grid(qua, 15.61, 15.11, 90, -1, 5, 15, 2, 10)
    XT.append(x)
    YT.append(y)
    ZT.append(z)
    izon_count += 1
    id_count += 2
    IZT.append(x.shape[2])
    JZT.append(x.shape[1])
    KZT.append(x.shape[0])

# gas recess
for qua in range(1,5): 
    x, y, z = o_grid(qua, 14.11, 0.5*14.11, 20, -6, 28.22, 15, 28, 50)
    XT.append(x)
    YT.append(y)
    ZT.append(z)
    izon_count += 1
    IZT.append(x.shape[2])
    JZT.append(x.shape[1])
    KZT.append(x.shape[0])

x, y, z = h_grid(0.5*14.11, 20, -6, 28.22, 15, 50)
XT.append(x)
YT.append(y)
ZT.append(z)
izon_count += 1
IZT.append(x.shape[2])
JZT.append(x.shape[1])
KZT.append(x.shape[0])

# wall recess
for qua in range(1,5): 
    x, y, z = o_grid(qua, 15.11, 14.11, 90, -6, 28.22, 15, 4, 50)
    XT.append(x)
    YT.append(y)
    ZT.append(z)
    izon_count += 1
    id_count += 1
    IZT.append(x.shape[2])
    JZT.append(x.shape[1])
    KZT.append(x.shape[0])

# liquid recess
for qua in range(1,5): 
    x, y, z = o_grid(qua, 15.61, 15.11, 90, -6, 28.22, 15, 2, 50)
    XT.append(x)
    YT.append(y)
    ZT.append(z)
    izon_count += 1
    id_count += 1
    IZT.append(x.shape[2])
    JZT.append(x.shape[1])
    KZT.append(x.shape[0])

# gas out
for qua in range(1,5): 
    x, y, z = o_grid(qua, 14.11, 0.5*14.11, 20, -34.22, 50, 15, 28, 85)
    XT.append(x)
    YT.append(y)
    ZT.append(z)
    izon_count += 1
    ibnd_count += 1
    IZT.append(x.shape[2])
    JZT.append(x.shape[1])
    KZT.append(x.shape[0])

x, y, z = h_grid(0.5*14.11, 20, -34.22, 50, 15, 85)
XT.append(x)
YT.append(y)
ZT.append(z)
izon_count += 1
ibnd_count += 1
IZT.append(x.shape[2])
JZT.append(x.shape[1])
KZT.append(x.shape[0])

# wall out
for qua in range(1,5): 
    x, y, z = o_grid(qua, 15.11, 14.11, 90, -34.22, 50, 15, 4, 85)
    XT.append(x)
    YT.append(y)
    ZT.append(z)
    izon_count += 1
    ibnd_count += 1
    IZT.append(x.shape[2])
    JZT.append(x.shape[1])
    KZT.append(x.shape[0])

# liquid out
for qua in range(1,5): 
    x, y, z = o_grid(qua, 15.61, 15.11, 90, -34.22, 50, 15, 2, 85)
    XT.append(x)
    YT.append(y)
    ZT.append(z)
    izon_count += 1
    ibnd_count += 1
    IZT.append(x.shape[2])
    JZT.append(x.shape[1])
    KZT.append(x.shape[0])

# out
for qua in range(1,5): 
    x, y, z = o_grid(qua, 25, 15.61, 90, -34.22, 50, 15, 40, 85)
    XT.append(x)
    YT.append(y)
    ZT.append(z)
    izon_count += 1
    ibnd_count += 3
    IZT.append(x.shape[2])
    JZT.append(x.shape[1])
    KZT.append(x.shape[0])

"""
import numpy as np
ls = np.array([1,2,3])
a, b, c = ls * 3
print(a, b, c)