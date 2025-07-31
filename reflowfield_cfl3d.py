import numpy as np
import pandas as pd
from struct import unpack
from lib.plot3dout import unformatted_fort13

#
#  read cfl3d file
#
with open('input files/injector A quarter.bin', 'rb') as f:
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
        #x = np.float32(x)
        X.append(x)
        
        buff = f.read(8 * size)
        y = np.array(unpack(f'{size}d', buff)).reshape((kzt, jzt, izt))
        #y = np.float32(y)
        Y.append(y)
        
        buff = f.read(8 * size)
        z = np.array(unpack(f'{size}d', buff)).reshape((kzt, jzt, izt))
        #z = np.float32(z)
        Z.append(z)
        f.read(4)

with open('input files/injector A quarter.inp', 'r') as f:
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
    
    f.readline()
    
    # skip
    for _ in range(izon + 1):
        f.readline()
    
    # read dimension
    f.readline()
    for _ in range(izon):
        f.readline()
    
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

#
#  governing equation
#
INSO_1 = 1 # U
INSO_4 = 1 # TM
INSO_5 = 1 # DK
INSO_7 = 1 # FL
NGAS = 2

#
#  reference state
#
UREF = 69.593
XREF = 1E-3
PREF11 = 1E5
TREF = 300 
RMXBAR = 288.5939026
GAMAREF = 1.3985
DENREF = PREF11/(RMXBAR*TREF) # 1.155
PREF = DENREF * UREF**2 # 5594

#
#  geometry
#
inlet_orifice_rad = 1E-3    # m
inlet_gas_rad = 14.11E-3
swirl_inn_rad = 15.11E-3
swirl_out_rad = 15.61E-3
swirl_height = 6E-3

inlet_area = np.pi * (swirl_out_rad**2 - swirl_inn_rad**2)

#
#  liquid inlet condition
#
lq_den = 1000        # kg/m3

lq_flowrate = 29.8E-3 # kg/s
lq_vel_z = lq_flowrate / (inlet_area * lq_den) * -1 # m/s
lq_vel = lq_flowrate / ( (swirl_out_rad - swirl_inn_rad) * swirl_height * lq_den )
lq_vel_tot = np.sqrt(lq_vel**2 + lq_vel_z**2)

lq_pre = 101300  # Pa

lq_vsc = 0.89E-3 # kg/m s
lq_re = lq_vel_tot * lq_den * 0.5E-3 / lq_vsc
lq_i = 0.16 * lq_re**(-0.125)
lq_dk = 1.5 * (lq_vel_tot * lq_i)**2
lq_de = 0.09**0.75 * lq_dk**1.5 / lq_i

lq_fm = [1, 0]      # 1 -> H2O, 2 -> AIR

#
#  gas inlet condition
#
gs_den = 1

gs_flowrate = 23E-3
gs_vel_z = gs_flowrate / (np.pi * inlet_gas_rad**2 * gs_den) * -1

gs_pre = 101300

gs_vsc = 1.849E-5
gs_re = gs_vel_z * gs_den * inlet_gas_rad / gs_vsc
gs_i = 0.16 * gs_re**(-0.125)
gs_dk = 1.5 * (gs_vel_z * gs_i)**2
gs_de = 0.09**0.75 * gs_dk**1.5 / gs_i

gs_fm = [0, 1]

#
#  internal condition
#
int_den = gs_den
int_u = 0
int_v = 0
int_w = 0
int_pre = 101300
int_tm = 300
int_am = 0
int_q = 0

#
#  fort.13
#
den, u, v, w, p, tm = [],[],[],[],[],[]
dk, de, am, q, fm = [],[],[],[],[]

# internal
for i in range(izon):
    izt = dim[i][0]
    jzt = dim[i][1]
    kzt = dim[i][2]
    
    den.append( np.full((kzt, jzt, izt), int_den/DENREF) )
    u.append(   np.full((kzt, jzt, izt), int_u/UREF)     )
    v.append(   np.full((kzt, jzt, izt), int_v/UREF)     )
    w.append(   np.full((kzt, jzt, izt), int_w/UREF)     )
    p.append(   np.full((kzt, jzt, izt), int_pre/PREF)   )
    tm.append(  np.full((kzt, jzt, izt), int_tm/TREF)   )
    
    vel_int_up = 1E-4
    vel_int_dw = 1E-6
    re_int_up = int_den * vel_int_up * inlet_gas_rad / gs_vsc
    re_int_dw = int_den * vel_int_dw * inlet_gas_rad / gs_vsc
    i_int_up = 0.16 * re_int_up**(-0.125)
    i_int_dw = 0.16 * re_int_dw**(-0.125)
    dk_int_up = 1.5 * (vel_int_up * i_int_up)**2
    dk_int_dw = 1.5 * (vel_int_dw * i_int_dw)**2
    de_int_up = 0.09**0.75 * dk_int_up**1.5 / i_int_up
    de_int_dw = 0.09**0.75 * dk_int_dw**1.5 / i_int_dw
    dk.append(np.random.uniform(dk_int_dw/UREF**2,        dk_int_up/UREF**2,        (kzt, jzt, izt)))
    de.append(np.random.uniform(de_int_dw/(UREF**3*XREF), de_int_up/(UREF**3*XREF), (kzt, jzt, izt)))
    
    am.append( np.full((kzt, jzt, izt), int_am) )
    q.append(  np.full((kzt, jzt, izt), int_q)  )
    fm.append( [np.full((kzt, jzt, izt), gs_fm[0]),
                np.full((kzt, jzt, izt), gs_fm[1])] )

def find_face(face, block):
    if face == 1: # idim
        return block[:, :, -1]
    elif face == 2: # i0
        return block[:, :, 0]
    elif face == 3: # jdim
        return block[:, -1, :]
    elif face == 4: # j0
        return block[:, 0, :]
    elif face == 5: # kdim
        return block[-1, :, :]
    elif face == 6: # k0
        return block[0, :, :]

def replace_face(face, block, value):
    if face == 1: # idim
        block[:, :, -1] = value
    elif face == 2: # i0
        block[:, :, 0] = value
    elif face == 3: # jdim
        block[:, -1, :] = value
    elif face == 4: # j0
        block[:, 0, :] = value
    elif face == 5: # kdim
        block[-1, :, :] = value
    elif face == 6: # k0
        block[0, :, :] = value
    return

# inlet condition
for i in range(izon):
    for face in range(1,7):
        """
        1001 -> general symmetry plane    -> symmetry
        1003 -> inflow/outflow            -> mass outlet
        2003 -> inflow w/ total condition -> air inlet
        2004 -> no-slip wall              -> wall
        2005 -> periodic in space         -> cyclic
        2009 -> inflow                    -> h2o inlet
        """
        if df_bnd[face-1].loc[i, 'BCTYPE'] == '2009':
            
            # den
            den_i = find_face(face, den[i])
            replace_face(face, den[i], np.full_like(den_i, lq_den/DENREF))
            
            # u,v,w
            x_i = find_face(face, X[i])
            y_i = find_face(face, Y[i])
            theta = np.atan2(y_i, x_i) + np.pi/2
            u_i = np.cos(theta) * lq_vel/UREF
            v_i = np.sin(theta) * lq_vel/UREF
            replace_face(face, u[i], u_i)
            replace_face(face, v[i], v_i)
            
            w_i = find_face(face, w[i])
            replace_face(face, w[i], np.full_like(w_i, lq_vel_z/UREF))
            
            # dk,de
            dk_i = find_face(face, dk[i])
            de_i = find_face(face, de[i])
            replace_face(face, dk[i], np.full_like(dk_i, lq_dk/UREF**2)        )
            replace_face(face, de[i], np.full_like(de_i, lq_de/(UREF**3*XREF)) )
            
            # fm
            fm_i_1 = find_face(face, fm[i][0])
            fm_i_2 = find_face(face, fm[i][1])
            replace_face(face, fm[i][0], np.full_like(fm_i_1, lq_fm[0]))
            replace_face(face, fm[i][1], np.full_like(fm_i_2, lq_fm[1]))

        elif df_bnd[face-1].loc[i, 'BCTYPE'] == '2003':
            
            # u,v,w
            w_i = find_face(face, w[i])
            replace_face(face, w[i], np.full_like(w_i, gs_vel_z/UREF))
            
            # dk, de
            dk_i = find_face(face, dk[i])
            de_i = find_face(face, de[i])
            replace_face(face, dk[i], np.full_like(dk_i, gs_dk/UREF**2)        )
            replace_face(face, de[i], np.full_like(de_i, gs_de/(UREF**3*XREF)) )
"""
# no-slip condition
for i in range(izon):
    for face in range(1,7):
        if df_bnd[face-1].loc[i, 'BCTYPE'] == '2004':
            
            # u,v,w
            u_tmp = find_face(face, u[i])
            v_tmp = find_face(face, v[i])
            w_tmp = find_face(face, w[i])
            replace_face(face, u[i], np.full_like(u_tmp, 0))
            replace_face(face, v[i], np.full_like(v_tmp, 0))
            replace_face(face, w[i], np.full_like(w_tmp, 0))
            
            # dk,de
            dk_tmp = find_face(face, dk[i])
            de_tmp = find_face(face, de[i])
            replace_face(face, dk[i], np.full_like(dk_tmp, 0))
            replace_face(face, de[i], np.full_like(de_tmp, 0))
"""
#
#
#
txt = unformatted_fort13(INSO_1, INSO_4, INSO_5, INSO_7, NGAS, izon, den, u, v, w, p, tm, dk, de, am, q, fm, 'input files/')
print(txt)