import numpy as np

XT = np.array([[1,2],[3,4]])
YT = np.array([[2,4],[6,8]])
XFLT = XT.flatten()
YFLT = YT.flatten()
COOR = np.concatenate(( [XFLT], [YFLT] ))
VEC = np.transpose(COOR)

ANGCNT = 1
RTN = np.array( [[np.cos(ANGCNT), np.sin(ANGCNT)], [-np.sin(ANGCNT), np.cos(ANGCNT)]] )

print(COOR)
print(RTN @ COOR)