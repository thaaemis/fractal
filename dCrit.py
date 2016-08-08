import numpy as np
import matplotlib.pyplot as plt
from fractions import Fraction
import itertools, sys
sys.path.insert(0, '/home/brian/GitHub/colormap/')
import colormaps as cmaps

def getConvergents(CFin):
    convergents, CF = [], [x for x in CFin]
    # Need irrational which these convergents are approaching.
    # an easy way is to pick the noble number nearby, ie.,
    # by appending a lot of 1s to the CF.
    for j in range(20):
        CF.append(1)
    for i in range(0,len(CF)+1):
        out = Fraction(0)
        for number in reversed(CF[:i]):
            result = number + out
            out = Fraction(1,result)
        convergents.append(out)
    return convergents[0:-20], convergents[-1]
    
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
    plt.xscale('symlog',linthreshx=1e-14)
    # plt.yscale('log')
    plt.show()

def plotDiff(base, increase, k=2.0, symlog = True):
    
    increase = int(increase)
    convergents0, alpha0 = getConvergents(base)
    d0, d0s = getDcrit(convergents0,alpha0,k)
    if symlog:
        offset = alpha0
    else:
        offset = 0
    plt.plot(alpha0-offset,d0,'k*',markersize=20)
    
    ds, alphas = [], []
    for i in range(len(base)-1):
        newCF = [x for x in base]
        newCF[i] += increase
        convergents, alpha = getConvergents(newCF)
        d, array = getDcrit(convergents, alpha, k)
        ds.append(d)
        alphas.append(alpha)
        plt.text(alpha-offset, d, str(i))

    labelStr = r'$c = $' + str(increase)

    plt.plot([x-offset for x in alphas],ds,'o',label=labelStr)
    return d0

def makeDiffPlots():
    alpha = [1,2,3,4,5,4,3,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    d0 = plotDiff(alpha, 1)
    for i in range(2,15):
       plotDiff(alpha,i)
    setDiffsPlot(alpha, d0, False)
    
def dCritIndexPlot(nmax = 5, length = 6, k = 2):
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
    
def compareCFs():
    convergents, alpha = getConvergents([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
    dNoble, array = getDcrit(convergents, alpha, 2)
    convergents, alpha = getConvergents([1,1,1,4,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
    dSudden, array = getDcrit(convergents, alpha, 2)
    convergents, alpha = getConvergents([1,2,3,4,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
    dGradUp, array = getDcrit(convergents, alpha, 2)
    convergents, alpha = getConvergents([1,2,2,4,3,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
    dGradRoughBoth, array = getDcrit(convergents, alpha, 2)
    convergents, alpha = getConvergents([1,2,3,4,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
    dGradRougherBoth, array = getDcrit(convergents, alpha, 2)
    convergents, alpha = getConvergents([1,2,3,4,3,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
    dGradBoth, array = getDcrit(convergents, alpha, 2)
    convergents, alpha = getConvergents([1,1,1,4,3,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
    dGradDown, array = getDcrit(convergents, alpha, 2)
    convergents, alpha = getConvergents([1,2,3,4,5,4,3,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
    dBigger, array = getDcrit(convergents, alpha, 2)

    print('Noble: ', float(dNoble))
    print('Sudden 4: ', float(dSudden))
    print('Gradual Up: ', float(dGradUp))
    print('Gradual Down: ', float(dGradDown))
    print('Gradual Both: ', float(dGradBoth))
    print('Gradual Both, rough sooner: ', float(dGradRoughBoth))
    print('Gradual Both, rough later: ', float(dGradRougherBoth))
    print('Higher bound: ', float(dBigger))

convergents, alpha = getConvergents([1,1,1,1,1,1,1,50,1,1,1,1,1,1,1])
d50Soon, array = getDcrit(convergents, alpha, 2)
convergents, alpha = getConvergents([1,1,1,1,1,1,1,1,50,1,1,1,1,1,1])
d50Late, array = getDcrit(convergents, alpha, 2)

print('Sooner: ', float(d50Soon))
print('Later : ', float(d50Late))

makeDiffPlots()
    
    
