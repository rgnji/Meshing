import numpy as np

a=np.array([2,3,4])
b=np.array([3,4,5])
c=[a,b]
print(len(c))
print(len(c[0]))
print(np.shape(c[0])[0])