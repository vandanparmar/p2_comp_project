import numpy as np


def fun(x,y):
	z = x
	a = y
	return 1

print(np.fromfunction(fun,(3,3)))

mass = np.array([[0.0,0.0],[0.0,0.0]])
print(mass.shape)
mass = np.reshape(mass,(4,-1))
mass
powed = np.power(mass,-1)
print(powed)
powed[powed>1e10]=1e10
print(powed)