from lib.laplace import solve_grid_laplace_sor
import numpy as np

def pipe(PIPRIN, PIPLEN, CNT, LQROUT, D1, D3):
    """
    CNT = [X, Y, Z]
    """
    BNDX = []
    BNDZ = []
    
    for qua in range(1,5):
        ANGSTR = (qua - 1) * 0.5 * np.pi
        ANGEND = qua * 0.5 * np.pi
        ANG = np.linspace(ANGSTR, ANGEND, D1+1)
        
        BNDX.append(np.cos(ANG))
        BNDZ.append(np.sin(ANG))
        
    BNDX[0][0]  = 1
    BNDX[0][-1] = 0
    BNDX[1][0]  = 0
    BNDX[1][-1] = -1
    BNDX[2][0]  = -1
    BNDX[2][-1] = 0
    BNDX[3][0]  = 0
    BNDX[3][-1] = 1 
    
    BNDZ[0][0]  = 0
    BNDZ[0][-1] = 1 
    BNDZ[1][0]  = 1 
    BNDZ[1][-1] = 0
    BNDZ[2][0]  = 0
    BNDZ[2][-1] = -1 
    BNDZ[3][0]  = -1
    BNDZ[3][-1] = 0
    
    for i in range(4):
        BNDX[i] *= PIPRIN
        BNDZ[i] *= PIPRIN
    
    X2D = np.zeros((D1+1, D1+1), float)
    Z2D = np.zeros((D1+1, D1+1), float)
    X2D[:, 1:-1] = np.linspace(BNDX[0][1:-1], BNDX[2][1:-1][::-1], D1+1)
    Z2D[:, 1:-1] = np.linspace(BNDZ[0][1:-1], BNDZ[2][1:-1][::-1], D1+1)
    X2D[:, 0]  = BNDX[3][::-1]
    X2D[:, -1] = BNDX[1]
    Z2D[:, 0]  = BNDZ[3][::-1]
    Z2D[:, -1] = BNDZ[1]
    
    X2D, Z2D, _ = solve_grid_laplace_sor(X2D, Z2D)
    
    X2D += CNT[0]
    Z2D += CNT[2]
    Y2D = np.full_like(X2D, -PIPLEN)
    
    THETA = np.acos(X2D / LQROUT)
    Y2DCV = np.sin(THETA) * LQROUT * -1
    
    X = np.broadcast_to(X2D, (D3+1, *X2D.shape))
    Z = np.broadcast_to(Z2D, (D3+1, *Z2D.shape))
    Y = np.linspace(Y2DCV, Y2D, D3+1)
    
    return X, Y, Z

def wall(PIPRIN, PIPLEN, CNT, LQROUT, LQRINN, D1, D2, D3):
    XT = []
    YT = []
    ZT = []
    
    x, y, z = pipe(PIPRIN, PIPLEN, CNT, LQROUT, D1, D3)
    XT.append(x)
    YT.append(y)
    ZT.append(z)
    
    X2DOUT = XT[0][0, :, :]
    Y2DOUT = YT[0][0, :, :]
    Z2DOUT = ZT[0][0, :, :]
    
    THETA = np.acos(X2DOUT / LQROUT)
    X2DINN = np.cos(THETA) * LQRINN
    Y2DINN = np.sin(THETA) * LQRINN * -1
    
    X = np.linspace(X2DINN, X2DOUT, D2+1)
    Y = np.linspace(Y2DINN, Y2DOUT, D2+1)
    Z = np.broadcast_to(Z2DOUT, (D2+1, *Z2DOUT.shape))
    XT.append(X)
    YT.append(Y)
    ZT.append(Z)
    
    for i in range(2):
        xx = XT[i] * -1
        yy = YT[i] * -1
        zz = ZT[i]
        XT.append(xx)
        YT.append(yy)
        ZT.append(zz)
    
    return XT, YT, ZT

def liquid_inlet(PIPRIN, PIPLEN, CNT, LQROUT, LQRINN, D1, D2, D3, LQD1):
    XT = []
    YT = []
    ZT = []
    
    x, y, z = wall(PIPRIN, PIPLEN, CNT, LQROUT, LQRINN, D1, D2, D3)
    XT += x
    YT += y
    ZT += z

    # ========================    
    WALLX1 = XT[1]
    WALLY1 = YT[1]
    WALLZ1 = ZT[1]
    WALLX2 = XT[3]
    WALLY2 = YT[3]
    WALLZ2 = ZT[3]
        
    X1Q1 = WALLX1[:, 0, :]
    Y1Q1 = WALLY1[:, 0, :]
    Z1Q1 = WALLZ1[:, 0, :]
    X1Q4 = WALLX1[:, ::-1, 0]
    Y1Q4 = WALLY1[:, ::-1, 0]
    Z1Q4 = WALLZ1[:, ::-1, 0]
    
    X2Q2 = WALLX2[:, ::-1, -1]
    Y2Q2 = WALLY2[:, ::-1, -1]
    Z2Q2 = WALLZ2[:, ::-1, -1]
    X2Q3 = WALLX2[:, -1, :]
    Y2Q3 = WALLY2[:, -1, :]
    Z2Q3 = WALLZ2[:, -1, :]
    
    COMX1 = np.concatenate((X1Q4, X1Q1[:, 1:]), axis=1)
    COMY1 = np.concatenate((Y1Q4, Y1Q1[:, 1:]), axis=1)
    COMZ1 = np.concatenate((Z1Q4, Z1Q1[:, 1:]), axis=1)
    
    COMX2 = np.concatenate((X2Q3, X2Q2[:, 1:]), axis=1)
    COMY2 = np.concatenate((Y2Q3, Y2Q2[:, 1:]), axis=1)
    COMZ2 = np.concatenate((Z2Q3, Z2Q2[:, 1:]), axis=1)
    
    ANG1 = np.atan2(COMY1, COMX1)
    R1 = np.sqrt(COMX1**2 + COMY1**2)
    ANG2 = np.atan2(COMY2, COMX2)
    R2 = np.sqrt(COMX2**2 + COMY2**2)
    ANG = np.linspace(ANG1, ANG2, LQD1+1)
    R = np.linspace(R1, R2, LQD1+1)
    XX1 = R * np.cos(ANG)
    YY1 = R * np.sin(ANG)
    ZZ1 = np.linspace(COMZ1, COMZ2, LQD1+1)
    XT.append(XX1)
    YT.append(YY1)
    ZT.append(ZZ1)
    
    XX2 = R * np.cos(ANG + np.pi)
    YY2 = R * np.sin(ANG + np.pi)
    XT.append(XX2)
    YT.append(YY2)
    ZT.append(ZZ1)
    
    YCNT = np.sin(np.acos(CNT[0] / LQROUT)) * LQROUT * -1
    ANGCNT = np.atan2(YCNT, CNT[0]) * -1
    RTN = np.array( [[np.cos(ANGCNT), -np.sin(ANGCNT)], [np.sin(ANGCNT), np.cos(ANGCNT)]] )
    for zon in range(len(XT)):
        XFLT = XT[zon].flatten()
        YFLT = YT[zon].flatten()
        COOR = np.concatenate(( [XFLT], [YFLT] ))
        CORRT = RTN @ COOR
        XT[zon] = CORRT[0].reshape(XT[zon].shape)
        YT[zon] = CORRT[1].reshape(YT[zon].shape)
    
    return XT, YT, ZT