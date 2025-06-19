import numpy as np
from struct import unpack
from lib.plot3dout import unformatted_fort13, binary_fort13, ascii_fort13

# 
#  read fort.12
# 
with open('fort.12', 'rb') as f:
    buff = f.read(4*3)
    IZON = unpack('<3i', buff)[1]
    
    IZT, JZT, KZT = [], [], []
    for i in range(IZON):
        buff = f.read(4*5)
        arr = unpack('<5i', buff)
        IZT.append(arr[1])
        JZT.append(arr[2])
        KZT.append(arr[3])
    
    XT, YT, ZT =[], [], []
    for i in range(IZON):
        size = IZT[i]*JZT[i]*KZT[i]
        
        f.read(4)
        buff = f.read(4*size)
        xflat = np.array(unpack(f'<{size}f', buff))
        XT.append(xflat.reshape((KZT[i], JZT[i], IZT[i])))
        f.read(4)
        
        f.read(4)
        buff = f.read(4*size)
        yflat = np.array(unpack(f'<{size}f', buff))
        YT.append(yflat.reshape((KZT[i], JZT[i], IZT[i])))
        f.read(4)
        
        f.read(4)
        buff = f.read(4*size)
        zflat = np.array(unpack(f'<{size}f', buff))
        ZT.append(zflat.reshape((KZT[i], JZT[i], IZT[i])))
        f.read(4)
        
# 
#  inlet blocks
# 
IBCZON = [1,2,3,4,5,
          25,31,37,43,49,
          26,32,38,44,50]
IDBC = [5,5,5,5,5,
        5,5,5,5,5,5,5,5,5,5]

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
RINLQ = 1E-3    # m
RINGS = 14.11E-3
RSWINN = 15.11E-3
RSWOUT = 15.61E-3
ANGIN = np.arccos(15.11 / 15.61)
PIPDIRX = 0.01511 - 0.01549
PIPDIRY = 0.00392 - 0.00246
PIPLEN = np.sqrt(PIPDIRX**2 + PIPDIRY**2)

#
#  liquid inlet condition
#
DNLQ = 1000        # kg/m3
FLOLQ = 39.72E-3/4 # kg/s
VELINLQ = FLOLQ / (np.pi*RINLQ**2 * DNLQ)
UINLQ = -VELINLQ * np.cos(ANGIN)
VINLQ = VELINLQ * np.sin(ANGIN)
WINLQ = 0
PINLQ = 101300  # Pa
DKINLQ = 0.015*VELINLQ**2
DEINLQ = 0.09*DKINLQ**1.5/(2*RINLQ)
FMLQ = [1, 0]      # 1 -> H2O, 2 -> AIR

#
#  gas inlet condition
#
DNGS = 1
FLOGS = 297.9E-3/4
VELINGS = FLOGS / (np.pi*RINGS**2 * DNGS)
UINGS = 0
VINGS = 0
WINGS = -VELINGS
PINGS = 101300
DKINGS = 0.015*VELINGS**2
DEINGS = 0.09*DKINGS**1.5/(2*RINGS)
FMGS = [0, 1]

#
#  internal condition
#
DNIN = 1
UIN = 0
VIN = 0
WIN = 0
PIN = 101300
TM = 300
AMIN = 0
QIN = 0

#
#  liquid inlet pipe
#
UPIPLQ = UINLQ * PIPDIRX / PIPLEN
VPIPLQ = VINLQ * PIPDIRY / PIPLEN

#
#  liquid swirl
#
WSWLQ = -FLOLQ / (np.pi*RSWOUT**2 - np.pi*RSWINN**2) / DNGS

#
#  recess chamber
#
WRECESS = -(FLOGS + FLOLQ) / (np.pi*RSWOUT**2) / DNGS

#
#
#
BLKDN, BLKU, BLKV, BLKW, BLKP, BLKTM = [], [], [], [], [], []
BLKDK, BLKDE, BLKQ, BLKAM, BLKFM = [], [], [], [], []

#
#  internal
#
for i in range(IZON):
    BLKDN.append( np.full((KZT[i], JZT[i], IZT[i]), DNIN/DENREF) )
    BLKU.append( np.full((KZT[i], JZT[i], IZT[i]), UIN/UREF) )
    BLKV.append( np.full((KZT[i], JZT[i], IZT[i]), VIN/UREF) )
    BLKW.append( np.full((KZT[i], JZT[i], IZT[i]), WIN/UREF) )
    BLKP.append( np.full((KZT[i], JZT[i], IZT[i]), PIN/PREF) )
    BLKTM.append( np.full((KZT[i], JZT[i], IZT[i]), TM/TREF) )
    
    DKINUP = 1E-4
    DKINDW = 1E-6
    DEINUP = 0.09*DKINUP**1.5/(2*15.61)
    DEINDW = 0.09*DKINDW**1.5/(2*15.61)
    BLKDK.append( np.random.uniform(DKINDW/UREF**2, DKINUP/UREF**2, (KZT[i], JZT[i], IZT[i])) )
    BLKDE.append( np.random.uniform(DEINDW/(UREF**3*XREF), DEINUP/(UREF**3*XREF), (KZT[i], JZT[i], IZT[i])) )
    
    BLKAM.append( np.full((KZT[i], JZT[i], IZT[i]), AMIN) ) 
    BLKQ.append( np.full((KZT[i], JZT[i], IZT[i]), QIN) )
    BLKFM.append( [np.full((KZT[i], JZT[i], IZT[i]), FMGS[0]), 
                   np.full((KZT[i], JZT[i], IZT[i]), FMGS[1])] )
"""
#
#  liquid inlet pipe
#

for i in [23,27,31,35,39]:
    BLKU[i-1][:,:,:] = UPIPLQ / UREF
    BLKV[i-1][:,:,:] = VPIPLQ / UREF
    
for i in [24,28,32,36,40]:
    BLKU[i-1][:,:,:] = - UPIPLQ / UREF
    BLKV[i-1][:,:,:] = - VPIPLQ / UREF
    
#
#  liquid swirl
#
for i in [6,7,8,9,25,26,29,30,33,34,37,38]:
    BLKW[i-1][:,:,:] = WSWLQ / UREF

#
#  gas inlet pipe
#
for i in [1,2,3,4,5]:
    BLKW[i-1][:,:,:] = WINGS / UREF

#
#  recess chamber
#
for i in [10,11,12,13,14,15,16,17,18,19,20,21,22]:
    BLKW[i-1][:,:,:] = WRECESS / UREF
"""
#
#  inlet boundary condition
#
for i in IBCZON[:5]:
    BLKDN[i-1][-1, :, :] = DNGS/DENREF
    BLKU[i-1][-1, :, :] = UINGS/UREF
    BLKV[i-1][-1, :, :] = VINGS/UREF
    BLKW[i-1][-1, :, :] = WINGS/UREF
    BLKP[i-1][-1, :, :] = PINGS/PREF
    BLKDK[i-1][-1, :, :] = DKINGS/UREF**2
    BLKDE[i-1][-1, :, :] = DEINGS/(UREF**3*XREF)
    BLKFM[i-1][0][-1, :, :] = FMGS[0]
    BLKFM[i-1][1][-1, :, :] = FMGS[1]

for i in IBCZON[5:10]:
    BLKDN[i-1][-1, :, :] = DNLQ/DENREF
    BLKU[i-1][-1, :, :] = UINLQ/UREF
    BLKV[i-1][-1, :, :] = VINLQ/UREF
    BLKW[i-1][-1, :, :] = WINLQ/UREF
    BLKP[i-1][-1, :, :] = PINLQ/PREF
    BLKDK[i-1][-1, :, :] = DKINLQ/UREF**2
    BLKDE[i-1][-1, :, :] = DEINLQ/(UREF**3*XREF)
    BLKFM[i-1][0][-1, :, :] = FMLQ[0]
    BLKFM[i-1][1][-1, :, :] = FMLQ[1]

for i in IBCZON[10:]:
    BLKDN[i-1][-1, :, :] = DNLQ/DENREF
    BLKU[i-1][-1, :, :] = -UINLQ/UREF
    BLKV[i-1][-1, :, :] = -VINLQ/UREF
    BLKW[i-1][-1, :, :] = WINLQ/UREF
    BLKP[i-1][-1, :, :] = PINLQ/PREF
    BLKDK[i-1][-1, :, :] = DKINLQ/UREF**2
    BLKDE[i-1][-1, :, :] = DEINLQ/(UREF**3*XREF)
    BLKFM[i-1][0][-1, :, :] = FMLQ[0]
    BLKFM[i-1][1][-1, :, :] = FMLQ[1]

#
#  no-slip condition at inlet
#
for i in [1,2,3,4,
          25,31,37,43,
          26,32,38,44]:
    BLKU[i-1][-1, 0, :] = 0
    BLKV[i-1][-1, 0, :] = 0
    BLKW[i-1][-1, 0, :] = 0

#
#  fort.13
#
txt = unformatted_fort13(INSO_1, INSO_4, INSO_5, INSO_7, NGAS, IZON, BLKDN, BLKU, BLKV, BLKW, BLKP, BLKTM, BLKDK, BLKDE, BLKAM, BLKQ, BLKFM)
print(txt)
"""
txt = binary_fort13(INSO_1, INSO_4, INSO_5, INSO_7, NGAS, IZON, BLKDN, BLKU, BLKV, BLKW, BLKP, BLKTM, BLKDK, BLKDE, BLKAM, BLKQ, BLKFM)
print(txt)
txt = ascii_fort13(INSO_1, INSO_4, INSO_5, INSO_7, NGAS, IZON, BLKDN, BLKU, BLKV, BLKW, BLKP, BLKTM, BLKDK, BLKDE, BLKAM, BLKQ, BLKFM)
print(txt)
"""