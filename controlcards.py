def group3(blocks):
    IZT = []
    JZT = []
    KZT = []

    for blk in blocks:
        IZT.append(blk.shape[2])
        JZT.append(blk.shape[1])
        KZT.append(blk.shape[0])