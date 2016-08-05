import numpy as np
import matplotlib.pyplot as plt
from fractions import Fraction
import itertools, sys
sys.path.insert(0, '/home/brian/GitHub/colormap/')
import colormaps as cmaps

def getConvergents(CFin,addNoble = False):
    convergents, CF = [], [x for x in CFin]
    if addNoble:
        for j in range(20):
            CF.append(1)
    for i in range(1,len(CF)+1):
        out = Fraction(0)
        for number in reversed(CF[:i]):
            result = number + out
            out = Fraction(1,result)
        convergents.append(out)
    return convergents
    
def getDcrit(convergents,alpha, k): #list of convergents and mean value (CF.append([1,1,1,1,1])
    dCrits = []
    for i in range(len(convergents)):
        separation = abs(convergents[i] - alpha)
        dCritHere = separation * convergents[i].denominator**k
        dCrits.append(dCritHere)
    dCrit = min(dCrits)
    return dCrit, dCrits

def setDiffsPlot(CF,d0,ySym = True):
    CFtext = [str(j)+',' for j in CF]
    bText = [CFtext[0],r'...,$a_i$+c,...',CFtext[-1]]
    CFtext = ''.join(CFtext)
    bText= ''.join(bText)
    CFtext = '[' + CFtext[:-1] + ']'
    bText = '[' + bText[:-1] + ']'
    print(CFtext)

    plt.ylabel(r'$d^{crit}_b - d^{crit}_a$',fontsize=20)
    plt.xlabel(r'Element changed',fontsize=20)
    xmin, xmax, ymin, ymax = plt.axis()
    if ySym:
        plt.yscale('symlog',linthreshy=1e-15)
        yLoc = [y*ymax for y in [.1, .01, .001]]
    else:
        yLoc = [y*(ymax-ymin)+ymin for y in [0.95, 0.85, 0.75]]
    plt.plot([0,xmax],[0,0],'k--',label='_')
    plt.text((xmax-xmin)*0.15,yLoc[0],r'$a = [a_i] =$'+CFtext,fontsize=15)
    plt.text((xmax-xmin)*0.15,yLoc[1],r'$b_i = $'+bText,fontsize=15)
    plt.text((xmax-xmin)*0.15,yLoc[2],r'$d_a^{crit} = $'+str(float(d0)),fontsize=15)
    plt.legend(loc='best')
    plt.show()

def plotDiff(base, increase, k=2):
    
    alpha = getConvergents(base,True)[-1]
    increase = int(increase)
    d0, d0s = getDcrit(getConvergents(base),alpha,k)
    ds2 = []
    
    for i in range(len(base)):
        newCF = [x for x in base]
        newCF[i] += increase
        convergents = getConvergents(newCF)
        alpha = getConvergents(newCF, True)[-1]
        d, array = getDcrit(convergents, alpha, k)
        ds2.append(d)

    labelStr = r'$c = $' + str(increase)
    plt.plot([(x-d0) for x in ds2],label=labelStr)
    return d0

def makeDiffPlots():
    alpha = [1,2]*12
    d0 = plotDiff(alpha, 1)
    plotDiff(alpha, 2)
    plotDiff(alpha, 5)
    plotDiff(alpha, 9)
    plotDiff(alpha, 20)
    setDiffsPlot(alpha, d0)
    
def dCritIndexPlot(nmax = 5, length = 6, k = 2.1):
    iteration = itertools.product(np.arange(1,nmax+1), repeat=length)
    dCrits, avgElements, firstInds = [], [], []
    plt.figure()
    cMap = cmaps.magma
    for i in iteration:
        i = list(i)
        maxEl = max(i)
        plt.subplot(int(nmax), 1, int(maxEl))
        convergents = getConvergents(i,addNoble = False)
        alpha = getConvergents(i, True)[-1]
        dCrit, ds = getDcrit(convergents, alpha, k)
        avgElement = np.mean(i)
        firstInd = i.index(maxEl)
        dCrits.append(dCrit)
        avgElements.append(avgElement)
        firstInds.append(firstInd)
        plt.plot(firstInd,dCrit,'r.',color=cMap(avgElement/maxEl),markersize=10)
        plt.axis([-0.2, length-0.8, 0, 0.5])
    plt.subplot(nmax,1,nmax)
    # plt.colorbar(orientation='horizontal')
    plt.show()
    
dCritIndexPlot()
    
    
