import numpy as np
import matplotlib.pyplot as plt
from fractal import getP

def getMeasure(epsilon):
    d = 0.05
    k = 2.25
    nmax = (2*d/epsilon)**(1./k)
    print(epsilon, nmax)

    x, p, gradP = getP(d, k, nmax, fareyMethod = 'maxDen', getN = False)
    maxP = max(p) # indicative of total measure
    return maxP
    
epsilons = [1e-1, 5e-2, 2e-2, 1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4, 5e-5, 2e-5, 1e-5, 5e-6, 2e-6, 1e-6]
maxPs = []
for e in epsilons:
    maxPs.append(getMeasure(e))
    
plt.plot(epsilons,maxPs,'ro')
plt.xscale('log')
plt.show()
