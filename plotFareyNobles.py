import matplotlib.pyplot as plt
import math
from fractal import getP
from selfSimilar import rational, makeFractionList, getFareyPath

phi = (math.sqrt(5)-1)/2.

def CFtoVal(CFarray):
    CF = [x for x in CFarray]
    val = 1/float(CF.pop())
    for i in range(len(CF)):
        element = CF.pop()
        val = 1/(element+val)
    return val

def valtoCF(val, stepsMax = 100):
    CF = []
    while True:
        CF.append(math.floor(1/val))
        val = 1/val - CF[-1]
        if (abs(val) < 1e-10) or (len(CF) > stepsMax):
            break
    return CF

def plotLine(value,axisoffset=0):
    try:
        len(value)
        value = CFtoVal(value)
    except TypeError:
        value = float(value)
        
    maxN = 100
    d = 0.2
    k = 2.1
    convergents = False
    farey, center, pathOut = getFareyPath(value, maxN, d, k, convergents)
    
    cmap = plt.get_cmap('spectral')
    for j in range(len(farey)-1):
        maxBound = max(valtoCF(farey[j+1].val))
        color = cmap(maxBound/10)
        plt.plot([i.val-axisoffset for i in farey[j:j+2]],
                [i.den for i in farey[j:j+2]],color=color)
        (diophantineMin, diophantineMax) = farey[j+1].diophantine
        plt.plot([diophantineMin-axisoffset,diophantineMax-axisoffset],
                [farey[j+1].den]*2,
                color = color, linewidth=10)
        plt.yscale('log')
        
goldenMean = valtoCF((math.sqrt(5)-1)/2)
# print([index for index,value in enumerate(goldenMean) if value > 1])


        
def plotAssortment(offset = 0):
    plotLine([1,400,1,1,1,1,2,1,1,1],axisoffset=offset)
    plotLine([2,3,1,9,16,3,1,1,1,1,4]*3,axisoffset=offset)
    for i in range(20):
        a = [1]*20
        a[i] = 2
        plotLine(a,axisoffset=offset)
        b = [1]*20
        b[i] = 3
        plotLine(b,axisoffset=offset)
        c = [1,2] * 10
        c[i] = 3
        plotLine(c,axisoffset=offset)
    plotLine([1,15,1,1,1,1,2,1,1,1],axisoffset=offset)
    plotLine([1]*100,axisoffset=offset)



plt.subplot(2,1,1)
plotAssortment(0)
plt.axis([0,1,1,1e3])
    
plt.subplot(2,1,2)
plotAssortment(phi)
plt.xscale('symLog',linthreshx=1e-7)
plt.axis([-0.5,0.5,1,1e3])
plt.show()
