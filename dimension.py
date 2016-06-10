from pylab import *
import numpy as np
from fractal import getP #(d, k, nmax), returns x, p, gradp
from scipy.interpolate import interp1d

x, p, gradp = getP(0.1, 2, 200)

def getLength(x,p,r):
    newX = np.arange(0,1,0.0001)
    newP = interp1d(x,p)
    xs, ys = [x[0]], [p[0]]
    theta = np.arange(0,np.pi/2.,0.001)
    length = 0
    while xs[-1] + r < 1:
        xNow = xs[-1] + r*np.cos(theta)
        yNow = ys[-1] + r*np.sin(theta)
        pNow = newP(xNow)
        match = np.argwhere(np.isclose(yNow,pNow, atol=0.0005))
        match = round(np.mean(match))
        xs.append(xNow[match])
        ys.append(yNow[match])
        length += r
    
    remaining = sqrt((xs[-1]-1)**2+(ys[-1]-p[-1])**2)
    length += remaining
        
    return length
        
lengths = []
rs = [1,.1, .01, .001, .0001, 1e-5]
for r in rs:
    lengths.append(getLength(x,p,r))
    
print(lengths)

for i in range(len(rs)-1):
    D = -np.log(lengths[i+1]/rs[i+1])/np.log(rs[i+1])
    print(D)