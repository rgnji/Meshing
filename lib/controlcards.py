def gridz_size(blocks):
    IZT = []
    JZT = []
    KZT = []
    cards = []

    for blk in blocks:
        IZT.append(blk.shape[2])
        JZT.append(blk.shape[1])
        KZT.append(blk.shape[0])
    
    for i in range(len(IZT)):
        temp = []

        temp.append(IZT[i])
        temp.append(JZT[i])
        temp.append(KZT[i])

        if i < 48:
            temp.append(1)
        else:
            temp.append(2)
        
        for i in range(6):
            temp.append(0.0)
        
        cards.append(temp)
    
    return cards

#=================================== group 4 ===================================
def patch(XT, block1, block2, izf1, izf2, dblk11, dblk12, dblk21, dblk22):
    """
    Parameters:
        XT: blocks
        block1: one-based
        block2: one-based
        izf1: 1~6 for block 1
        izf2: 1~6 for block 2
        dblk11: 1 for (1~XT[block1].shape[x]), 0 for reversed
        dblk12: 1 for (1~XT[block1].shape[x]), 0 for reversed
        dblk21: 1 for (1~XT[block2].shape[x]), 0 for reversed
        dblk22: 1 for (1~XT[block2].shape[x]), 0 for reversed
    """
    IFCYC = 1

    IZB1 = block1
    IZF1 = izf1
    IZB2 = block2
    IZF2 = izf2

    def running(izf):
        if izf == 1:
            first = 1
            second = 0
        elif izf == 2:
            first = 2
            second = 0
        else:
            first = 2
            second = 1
        return first, second
    
    def running2(direc, end):
        if direc:
            ijz11 = 1
            ijz12 = end
        else:
            ijz11 = end
            ijz12 = 1
        return ijz11, ijz12

    
    d11, d12 = running(izf1)
    IJZ11, IJZ12 = running2(dblk11, XT[block1-1].shape[d11])
    JKZ11, JKZ12 = running2(dblk12, XT[block1-1].shape[d12])
    d21, d22 = running(izf2)
    IJZ21, IJZ22 = running2(dblk21, XT[block2-1].shape[d21])
    JKZ21, JKZ22 = running2(dblk22, XT[block2-1].shape[d22])

    INONUF = 0

    return[[IFCYC, IZB1, IZF1, IJZ11, IJZ12, JKZ11, JKZ12, INONUF],
           [IZB2, IZF2, IJZ21, IJZ22, JKZ21, JKZ22]]

def o_patch(XT, startblock):
    """
    Parameters:
        startblocks: one-based, first quadrant
    Returns:
        output: 3-D list
    """
    output = []

    for i in range(3):
        output.append(patch(XT, startblock+i, startblock+1+i, 1, 2, 1, 1, 1, 1))
    output.append(patch(XT, startblock+3, startblock, 1, 2, 1, 1, 1, 1))

    return output

def h_patch(XT, startblock):
    """
    Parameters:
        startblocks: one-based, first quadrant
    Returns:
        output: 3-D list
    """
    output = []
    f = [4, 1, 3, 2]

    for i in range(4):
        if i+1 == 1 or 2:
            direc = 1
        else:
            direc = 0
        output.append(patch(XT, startblock+i, startblock+4, 3, f[i], 1, 1, direc, 1))
    
    return output

def h_patch_liquid(XT, startblock, endblock):
    output = []
    f = [4, 1, 3, 2]

    for i in range(4):
        if i+1 == 1 or 2:
            direc = 1
        else:
            direc = 0
        output.append(patch(XT, startblock+i*6, endblock, 1, f[i], 1, 1, direc, 1))
    
    return output

def o_inner_patch(XT, startblock_out, startblock_in):
    """
    Parameters:
        startblocks: one-based, first quadrant
    Returns:
        output: 3-D list
    """
    output = []
    
    for i in range(4):
        output.append(patch(XT, startblock_out+i, startblock_in+i, 3, 4, 1, 1, 1, 1))
    
    return output

