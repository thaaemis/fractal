import matplotlib, math, itertools
import matplotlib.pyplot as plt
import numpy as np



def CFtoVal(CFarray):
    CF = [x for x in CFarray]
    val = 1/float(CF.pop())
    for i in range(len(CF)):
        element = CF.pop()
        val = 1/(float(element)+val)
    return val

def makeTrees(nmax = 7, length = 5):

    iteration = itertools.product(np.arange(1,nmax+1), repeat=length)
    valBoundArray = []
    for i in iteration:
        elements = [float(g) for g in i]
        val = CFtoVal(elements)
        maxBound = max(elements)
        valBoundArray.append((val,maxBound))
    
    cmap = matplotlib.cm.get_cmap('spectral')
    for j in valBoundArray:
        plt.plot([j[0]], [0], '.', color=cmap(j[1]/7.))
        
    plt.show()
    
makeTrees()
