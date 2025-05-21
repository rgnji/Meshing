import numpy as np
from struct import unpack, pack

with open("fort12.bin.xyz", "rb") as f12:
    # read fort.12 as double precision
    data = f12.read(8)
    IZON = unpack("<i", data)[0]  # number of blocks

    dim = []
    for i in range(IZON):
        data = f12.read(8*3)
        dim.append(unpack("<i", data)) # (IZT, JZT, KZT)
    
    coor = np.frombuffer(f12.read(), dtype=np.float64)

    XT = []
    YT = []
    ZT = []
    start = 0
    end = 0
    for i in range(IZON):
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
NGAS = 1

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
IBCZON = [1,2,3,4,5,
          42,43,48,49,54,55,60,61,66,67]
IDBC = [6,6,6,6,6,
        2,2,2,2,2,2,2,2,6,6]
ITYBC = [1,1,1,1,1,
         2,2,2,2,2,2,2,2,2,2] # 1 for gas inlet, 2 for liquid inlet

inlet_liquid_velocity = mass_rate_liquid / (RHOREF_liquid * inlet_liquid_r**2 * np.pi)
inlet_gas_velocity = mass_rate_gas / (RHOREF_gas * inlet_gas_r**2 * np.pi)

inlet_liquid_q = 1
inlet_gas_q = 0

inlet_liquid_density = 1000 # kg/m3
inlet_gas_density = 1

k_liquid = 3/2 * (0.05 * (mass_rate_liquid/density_liquid * inlet_liquid_r**2*np.pi))
k_gas = 3/2 *(0.05 * (mass_rate_gas/density_gas * inlet_gas_r**2*np.pi))

epsilon_liquid = 0.09**0.75 * k_liquid**1.5 / (2 * inlet_liquid_r)
epsilon_gas = 0.09**0.75 * k_gas**1.5 / (2 * inlet_gas_r)


# ============================================
# ============================================
# generate initial flow field
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

with open('fort13.bin.xyz', 'rb') as f13:
    f13.write(pack('<i', INSO_1))
    f13.write(pack('<i', INSO_4))
    f13.write(pack('<i', INSO_5))
    f13.write(pack('<i', INSO_7))
    f13.write(pack('<i', NGAS))

    for i in range(IZON): # mass inlet => inlet density, quality need to be changed
        izt, jzt, kzt = dim[i]

        blk_den = np.full((kzt, jzt, izt), density / RHOREF_gas)
        blk_q = np.full((kzt, jzt, izt), q)

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

        if i+1 in IBCZON:
            blk_index = IBCZON.index(i+1)
            if ITYBC[blk_index] == 1:
                face_check(IDBC[blk_index], blk_u, )
                face_check(IDBC[blk_index], blk_den, inlet_gas_density)
                face_check(IDBC[blk_index], blk_q, inlet_gas_q)
            elif ITYBC[blk_index] == 2:
                face_check(IDBC[blk_index], blk_den, inlet_liquid_density)
                face_check(IDBC[blk_index], blk_q, inlet_liquid_q)

        blk_den.astype(np.float64).tofile(f13)
        blk_u.astype(np.float64).tofile(f13)
        blk_v.astype(np.float64).tofile(f13)
        blk_w.astype(np.float64).tofile(f13)
        blk_p.astype(np.float64).tofile(f13)

        blk_dk.astype(np.float64).tofile(f13)
        blk_de.astype(np.float64).tofile(f13)

        blk_q.astype(np.float64).tofile(f13)
