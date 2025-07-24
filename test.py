from struct import unpack

with open('input files/injector A.bin', 'rb') as f:
    buff=f.read(4*3)
    blocks=unpack('<3i', buff)[1]
    print(unpack('<3i', buff))
    
    buff=f.read(4)
    print(unpack('<i', buff))
    for i in range(blocks):
        buff=f.read(4*3)
        print(unpack('<3i', buff))
    
    buff=f.read(4)
    print(unpack('<i', buff))