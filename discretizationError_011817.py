import numpy as np
from scipy.interpolate import interp1d
from pylab import *
from fractal import getP
import pandas as pd

pExactFilename = 'highRes_pApprox.pkl'

def getExactP(fareyLevel, save = False, filename = pExactFilename):
    d, k = 0.1, 2.1
    x, p, gradp = getP(d, k, nmax = fareyLevel, fareyMethod = 'treeSteps', getN = False)
    pData = pd.DataFrame(columns=('x', 'p', 'gradp'))
    pData['x'] = x
    pData['p'] = p
    pData['gradp'] = gradp
    
    if save:
        
        pData.to_pickle(filename)
    
    return x, p, gradp
            
def retrieveP(filename = pExactFilename, plotAll = False):
    pData = pd.read_pickle(filename)
    
    if plotAll:
        plot(pData['x'],pData['gradp'])
        show()
    
    x = pData['x']
    p = pData['p']
    gradp = pData['gradp'] 
    
    return x, p, gradp
    
def retrieveDicts(low, high):
    x, p, gradp = {}, {}, {}
    interpFunctions = {}
    for i in range(low,high+1):
        x[i], p[i], gradp[i] = getExactP(i)
    return x, p, gradp, interpFunctions    

def analyzeP():
    def getHist(x, p, gradp):

        # analyze largest and smallest gaps between coordinates in exact grid
        xFilt = x[1::2]
        diffX = np.abs(np.diff(xFilt))
        print(max(diffX), min(diffX), np.std(diffX))
        hist, binEdges = np.histogram(np.log(diffX), bins=200)
        figure(1)
        plot(binEdges[1:], hist)
        yscale('log')
        
        return max(diffX), min(diffX), np.mean(diffX), np.std(diffX)
    
    low, high = 1,9
    x, p, gradp, interpFunctions = retrieveDicts(low, high)
    for i in range(high, low-1, -1):
        upper, lower, mean, std = getHist(x[i], p[i], gradp[i])
        
        figure(2)
        f = interp1d(x[i], gradp[i])
        interpFunctions[i] = f
        
        # plot(i, upper, 'g')
        # plot(i, lower, 'r')
        # errorbar(i, mean, yerr=std, color='b')
        plot(x[i], interpFunctions[i](x[i]),'o', markersize=12-i, label=i)
    legend(loc='best')    
    show()
    
def saveManyP():
    for i in range(13,19):
        print('getting ',i)
        getExactP(i, save=True)
        
def measureConvergence():
    x, p, gradp = retrieveP()
    # need to shift x-points for interpolation
    xNew = np.asarray(x)
    xMax = np.max(x)
    for i, xPt in enumerate(x):
        if i % 2 == 1:
            xNew[i] = xPt*1.000000000000001*(1.+np.sqrt(5.))/2./xMax
        else:
            xNew[i] = xPt*0.999999999999999*(1.+np.sqrt(5.))/2./xMax
    gradpExact = interp1d(xNew, gradp)
    numPts = len(x)
    #plot(x, [gradpExact(i) for i in x],'.')
    print(p[numPts-1])
    
    # look at dumb grid way
    nVals = np.arange(10,5000,10)
    integrals = []
    for n in nVals:
        gridPts = np.arange(0.0,max(x),max(x)/float(n))
        # get Riemann integral
        riemannInt = 0
        for i in range(0,n):
            point = float(i)/float(n)
            riemannInt += gradpExact(point)/float(n)
            color = 'r' if n < 20 else 'b'
            # plot(point, gradpExact(point),'o',markersize=np.log(n), color=color)
        integrals.append(riemannInt)
    plot(nVals,integrals,'ro')
    plot([0,1000],[p[numPts-1],p[numPts-1]],'g--')
    
    # get fractal grid approximations
    low, high = 1, 10
    gridSize, fractPs = [],[]
    xDict, pDict, gradpDict, interpFunctions = retrieveDicts(low, high)
    for i in range(1,11):
        gridSize.append(len(pDict[i]))
        fractPs.append(pDict[i][-1])
    plot(gridSize,fractPs,'b*')
    xscale('log')
    show()
measureConvergence()
# saveManyP()
