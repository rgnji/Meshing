import numpy as np
from struct import unpack

from lib.plot3dout import unformatted_fort13, binary_fort13

# 
#  read fort.12
# 
with open("fort12.bin.xyz", "rb") as f12:
    # read fort.12 as double precision
    data = f12.read(4)
    IZON = unpack("<i", data)[0]  # number of blocks

    dim = []
    for i in range(IZON):
        data = f12.read(4*3)
        dim.append(unpack("<3i", data)) # (IZT, JZT, KZT)
    
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

# 
#  inlet blocks
# 
IBCZON = [1,2,3,4,5,
          42,43,48,49,54,55,60,61,66,67]
IDBC = [5,5,5,5,5,
        5,5,5,5,5,5,5,5,5,5]
# 1 for gas inlet, 2 for liquid inlet
ITYBC = [[1],[1],[1],[1],[1],
         [2],[2],[2],[2],[2],[2],[2],[2],[2],[2]]

#
#  governing equation
#
INSO_1 = 1 # U
INSO_4 = 0 # TM
INSO_5 = 1 # DK
INSO_7 = 1 # FL
NGAS = 2

#
#  constant
#
# molecular weight
MWLQ = 18 * 10**-3               
MWGS = (0.75 * 28 + 0.25 * 32) * 10**-3
# universal gas constant
RREFLQ = 8.31446261815324 / MWLQ 
RREFGS = 8.31446261815324 / MWGS

#
#  reference state
#
UREF = 1
TREF = 1
XREF = 1
PREF = 1
RHOREFLQ = PREF / (RREFLQ * TREF)
RHOREFGS = PREF / (RREFGS * TREF)

#
#  geometrt
#
RINLQ = 1E-3    # m
RINGS = 14.11E-3
ANGIN = np.arccos(15.11 / 15.61)

#
#  liquid inlet condition
#
DNLQ = 1000    # kg/m3
VELINLQ = 340  # m/s
FLOLQ = VELINLQ * (RINLQ**2*np.pi) * DNLQ # kg/s
UINLQ = VELINLQ * np.cos(ANGIN)
VINLQ = VELINLQ * np.sin(ANGIN)
WINLQ = 0
PINLQ = 101300 # Pa
DKINLQ = 0.015*VELINLQ**2
DEINLQ = 0.09*DKINLQ**1.5/(2*RINLQ)
FMLQ = [1, 0] # 1 -> H2O, 2 -> AIR

#
#  gas inlet condition
#
DNGS = 1
VELINGS = 340
FLOGS = VELINGS * (RINGS**2*np.pi) * DNGS
UINGS = 0
VINGS = 0
WINGS = VELINGS
PINGS = 101300
DKINGS = 0.015*VELINGS**2
DEINGS = 0.09*DKINGS**1.5/(2*RINGS)
FMGS = [0, 1]

#
#  internal initial condition
#
DNIN = 1
UIN = 0
VIN = 0
WIN = 0
PIN = 101300
QIN = 0

#
#
#
BLKDN, BLKU, BLKV, BLKW, BLKP = [], [], [], [], []
BLKDK, BLKDE, BLKQ, BLKAM, BLKFM = [], [], [], [], []

#
#
#
for i in range(IZON):
    IZT, JZT, KZT = dim[i]
    
    BLKDN = np.full((KZT, JZT, IZT), DNIN/RHOREFGS)
    BLKU = np.full((KZT, JZT, IZT), UIN/UREF)
    BLKV = np.full((KZT, JZT, IZT), VIN/VREF)
    BLKW = np.full((KZT, JZT, IZT), WIN/WREF)
    BLKP = np.full((KZT, JZT, IZT), PIN/PREF)
    
    BLKDK = np.full((KZT, JZT, IZT), )