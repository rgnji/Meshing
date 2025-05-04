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
def patch(XT, startblock, endblock, direction):
    """
    Parameters:
        startblock: one-based
        endblock: one-based
        direction: 1~3
    """
    IFCYC = 1

    if direction == 1: # I
        f = 1
        first = 1
        second = 0
    elif direction == 2: # J
        f = 3
        first = 2
        second = 0
    elif direction == 3: # K
        f = 5
        first = 2
        second = 1

    IZB1 = startblock
    IZF1 = f
    IZB2 = endblock
    IZF2 = f + 1
    
    IJZ11 = 1
    IJZ12 = XT[startblock-1].shape[first] 
    JKZ11 = 1
    JKZ12 = XT[startblock-1].shape[second] 
    IJZ21 = 1
    IJZ22 = XT[startblock-1].shape[first] 
    JKZ21 = 1
    JKZ22 = XT[startblock-1].shape[second] 

    INONUF = 0

    return [[IFCYC, IZB1, IZF1, IJZ11, IJZ12, JKZ11, JKZ12, INONUF],
            [IZB2, IZF2, IJZ21, IJZ22, JKZ21, JKZ22]]

def o_patch(XT, startblock):
    """
    Parameters:
        startblocks: one-based
    Returns:
        output: 3-D list
    """
    output = []

    for i in range(3): # 0,1,2
        output.append(patch(XT, startblock+i, startblock+1+i, 1))
    
    output.append(patch(XT, startblock+3, startblock, 1))

    return output

def h_patch(XT, startblock):
    """
    Parameters:
        startblocks: one-based
    Returns:
        output: 3-D list
    """
    output = []

    for i in range(4):
        output.append(patch(XT, startblock+i, startblock+4, 2))
    
    return output

def o_inner_patch(XT, startblock_out, startblock_in):
    """
    Parameters:
        startblocks: one-based
    Returns:
        output: 3-D list
    """
    output = []
    
    for i in range(4):
        output.append(patch(XT, startblock_out+i, startblock_in+i, 2))
    
    return output