import numpy as np
from struct import unpack, pack

from lib.plot3dout import unformatted_fort13, binary_fort13

# ============================================
# ============================================
# geometrical condition
with open("fort12.bin.xyz", "rb") as f12:
    # read fort.12 as double precision
    data = f12.read(4)
    blocks = unpack("<i", data)[0]  # number of blocks

    dim = []
    for i in range(blocks):
        data = f12.read(4*3)
        dim.append(unpack("<3i", data)) # (IZT, JZT, KZT)
    
    coor = np.frombuffer(f12.read(), dtype=np.float64)

    XT = []
    YT = []
    ZT = []
    start = 0
    end = 0
    for i in range(blocks):
        now = dim[i]
        size = now[0] * now[1] * now[2]

        start = end
        end += size # end index
        XT.append( coor[start:end].reshape(dim[i][::-1]) )
        start = end
        end += size
        YT.append( coor[start:end].reshape(dim[i][::-1]) )
        start = end
        end += size
        ZT.append( coor[start:end].reshape(dim[i][::-1]) )

# ============================================
# ============================================
# governing equation
INSO_1 = 1 # U
INSO_4 = 0 # TM
INSO_5 = 1 # DK
INSO_7 = 1 # FL
NGAS = 2

IZON = blocks

# constant
MW_liquid = 18 * 10**-3
MW_gas = (0.75 * 28 + 0.25 * 32) * 10**-3
RREF_liquid = 8.31446261815324 / MW_liquid
RREF_gas = 8.31446261815324 / MW_gas

inlet_liquid_r = 0.001
inlet_gas_r = 0.05
mass_rate_liquid = 0
mass_rate_gas = 0
density_liquid = 1000 
density_gas = 1 


# reference state
UREF = 1
TREF = 1
XREF = 1
PREF = 1
RHOREF_liquid = PREF / (RREF_liquid * TREF)
RHOREF_gas = PREF / (RREF_gas * TREF)

# internal initial condition
u = 0 
v = 0 
w = 0 
p = 101300 # Pa
q = 0 # quality, 1 for pure liquid and 0 for pure gas
density = 1


# boundary condition
# current bc: conserved mass flow rate
# mass flow rate => velocity => u, v, w
IBCZON = [1,2,3,4,5,
          42,43,48,49,54,55,60,61,66,67]
IDBC = [5,5,5,5,5,
        5,5,5,5,5,5,5,5,5,5]
# 1 for gas inlet, 2 for liquid inlet
ITYBC = [[1],[1],[1],[1],[1],
         [2],[2],[2],[2],[2],[2],[2],[2],[2],[2]]

fm_gas_inlet = [0, 1]
fm_liquid_inlet = [1, 0]

inlet_liquid_velocity = mass_rate_liquid / (RHOREF_liquid * inlet_liquid_r**2 * np.pi)
inlet_gas_velocity = mass_rate_gas / (RHOREF_gas * inlet_gas_r**2 * np.pi)

inlet_liquid_q = 1
inlet_gas_q = 0

inlet_liquid_density = 1000 # kg/m3
inlet_gas_density = 1

inlet_liquid_k = 3/2 * (0.05 * (mass_rate_liquid/density_liquid * inlet_liquid_r**2*np.pi))
inlet_gas_k = 3/2 *(0.05 * (mass_rate_gas/density_gas * inlet_gas_r**2*np.pi))

inlet_liquid_epsilon = 0.09**0.75 * inlet_liquid_k**1.5 / (2 * inlet_liquid_r)
inlet_gas_epsilon = 0.09**0.75 * inlet_gas_k**1.5 / (2 * inlet_gas_r)


# ============================================
# ============================================
def face_check(face, block, bc):
    if face == 1:
        block[:, :, -1] = bc
    elif face == 2:
        block[:, :, 0] = bc
    elif face == 3:
        block[:, -1, :] = bc
    elif face == 4:
        block[:, 0, :] = bc
    elif face == 5:
        block[-1, :, :] = bc
    elif face == 6:
        block[0, :, :] = bc
    return

def initial_uvw(block, face, velocity): # mass flow rate => velocity => u, v, w
    if face == 1: # I-max
        i, j, k = XT[block].shape[2]-1, 0, 0
        neighbor = np.array([XT[block][k][j][i-1], YT[block][k][j][i-1], ZT[block][k][j][i-1]])
    elif face == 2: # I-1
        i, j, k = 0, 0, 0
        neighbor = np.array([XT[block][k][j][i+1], YT[block][k][j][i+1], ZT[block][k][j][i+1]])
    elif face == 3: # J-max
        i, j, k = 0, XT[block].shape[1]-1, 0
        neighbor = np.array([XT[block][k][j-1][i], YT[block][k][j-1][i], ZT[block][k][j-1][i]])
    elif face == 4: # J-1
        i, j, k = 0, 0, 0
        neighbor = np.array([XT[block][k][j+1][i], YT[block][k][j+1][i], ZT[block][k][j+1][i]])
    elif face == 5: # K-max
        i, j, k = 0, 0, XT[block].shape[0]-1
        neighbor = np.array([XT[block][k-1][j][i], YT[block][k-1][j][i], ZT[block][k-1][j][i]])
    elif face == 6: # K-1
        i, j, k = 0, 0, 0
        neighbor = np.array([XT[block][k+1][j][i], YT[block][k+1][j][i], ZT[block][k+1][j][i]])

    boundary = np.array([XT[block][k][j][i], YT[block][k][j][i], ZT[block][k][j][i]])
    direction = (neighbor - boundary) / np.linalg.norm(neighbor - boundary)
    uu, vv, ww = velocity * direction
    return uu, vv, ww

# generate fort.13 
den_blocks, u_blocks, v_blocks, w_blocks, p_blocks = [], [], [], [], []
dk_blocks, de_blocks, am_blocks, q_blocks, fm_blocks = [], [], [], [], []

for i in range(IZON): # mass flow rate inlet => u, v, w, density, (quality), mass fraction need to be changed
    izt, jzt, kzt = dim[i]

    # internal flow field 
    blk_den = np.full((kzt, jzt, izt), density / RHOREF_gas)
    blk_u = np.full((kzt, jzt, izt), u / UREF)
    blk_v = np.full((kzt, jzt, izt), v / UREF)
    blk_w = np.full((kzt, jzt, izt), w / UREF)
    blk_p = np.full((kzt, jzt, izt), p / PREF)

    blk_dk = np.random.uniform(1e-6 / UREF**2, 
                                1e-4 / UREF**2, 
                                (kzt, jzt, izt))
    blk_de = np.random.uniform((0.09**0.75 * 1e-6**0.75 / inlet_gas_r) / (UREF**3 / XREF), 
                                (0.09**0.75 * 1e-4**0.75 / inlet_gas_r) / (UREF**3 / XREF), 
                                (kzt, jzt, izt))
    
    blk_am = np.zeros((kzt, jzt, izt))
    blk_q = np.full((kzt, jzt, izt), q)
    blk_fm = [np.full((kzt, jzt, izt), 0), 
                np.full((kzt, jzt, izt), 1)]

    # deal with boundary face
    if i+1 in IBCZON:
        blk_index = IBCZON.index(i+1)

        # if gas inlet
        if 1 in ITYBC[blk_index]:
            inlet_gas_u, inlet_gas_v, inlet_gas_w = initial_uvw(i, IDBC[blk_index], inlet_gas_velocity / UREF)
            face_check(IDBC[blk_index], blk_u, inlet_gas_u)
            face_check(IDBC[blk_index], blk_v, inlet_gas_v)
            face_check(IDBC[blk_index], blk_w, inlet_gas_w)

            face_check(IDBC[blk_index], blk_dk, inlet_gas_k / UREF**2)
            face_check(IDBC[blk_index], blk_de, inlet_gas_epsilon / (UREF**3 / XREF))

        # if liquid inlet
        if 2 in ITYBC[blk_index]:
            face_check(IDBC[blk_index], blk_den, inlet_liquid_density / RHOREF_liquid)

            inlet_liquid_u, inlet_liquid_v, inlet_liquid_w = initial_uvw(i, IDBC[blk_index], inlet_liquid_velocity / UREF)
            face_check(IDBC[blk_index], blk_u, inlet_liquid_u)
            face_check(IDBC[blk_index], blk_v, inlet_liquid_v)
            face_check(IDBC[blk_index], blk_w, inlet_liquid_w)

            face_check(IDBC[blk_index], blk_dk, inlet_liquid_k / UREF**2)
            face_check(IDBC[blk_index], blk_de, inlet_liquid_epsilon / (UREF**3 / XREF))

            for kk in range(NGAS):
                face_check(IDBC[blk_index], blk_fm[kk], fm_liquid_inlet[kk])

    den_blocks.append(blk_den)
    u_blocks.append(blk_u)
    v_blocks.append(blk_v)
    w_blocks.append(blk_w)
    p_blocks.append(blk_p)
    dk_blocks.append(blk_dk)
    de_blocks.append(blk_de)
    am_blocks.append(blk_am)
    q_blocks.append(blk_q)
    fm_blocks.append(blk_fm)

# establish fort.13
txt = binary_fort13(INSO_1, INSO_4, INSO_5, INSO_7, NGAS, IZON, den_blocks, u_blocks, v_blocks, w_blocks, p_blocks, dk_blocks, de_blocks, am_blocks, q_blocks, fm_blocks)
print(txt)

txt = unformatted_fort13(INSO_1, INSO_4, INSO_5, INSO_7, NGAS, IZON, den_blocks, u_blocks, v_blocks, w_blocks, p_blocks, dk_blocks, de_blocks, am_blocks, q_blocks, fm_blocks)
print(txt)