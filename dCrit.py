import numpy as np
import matplotlib.pyplot as plt
from fractions import Fraction
import itertools, sys, matplotlib
sys.path.insert(0, '/home/brian/GitHub/colormap/')
import colormaps as cmaps

#Colorbar info
vMin = 3.5
vMax = 6.5

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
    CFtext = [str(j)+',' for j in CF[:5]]
    bText = [CFtext[0],r'...,$a_i$+c,...',CFtext[-1]]
    CFtext = ''.join(CFtext)
    bText= ''.join(bText)
    CFtext = '[' + CFtext[:-1] + ',...]'
    bText = '[' + bText[:-1] + ']'
    print(CFtext)
    # goldenText = plt.text(0.64, 0.38, r'Golden mean, $\varphi$',fontsize=20)

    plt.ylabel(r'$d_{crit}$',fontsize=30)
    plt.xlabel(r'Value',fontsize=30)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.axis([0,0.99,0,0.5])
    xmin, xmax, ymin, ymax = plt.axis()
    if ySym:
        plt.yscale('symlog',linthreshy=1e-15)
        yLoc = [y*ymax for y in [.1, .01, .001]]
    else:
        yLoc = [y*(ymax-ymin)+ymin for y in [0.95, 0.8, 0.65]]
    plt.plot([0,xmax],[0,0],'k--',label='_')
    # plt.text((xmax-xmin)*0.15,yLoc[0],r'$a = [a_i] =$'+CFtext,fontsize=25)
    # plt.text((xmax-xmin)*0.15,yLoc[1],r'$b_i = $'+bText,fontsize=25)
    # plt.text((xmax-xmin)*0.15,yLoc[2],r'$d_a^{crit} = $'+str(float(d0)),fontsize=25)
    plt.text(0.35, 0.44, r'$\gamma(\alpha) = \sum_{i=1}^{n} (\alpha_{i+1} - \alpha_i)^2$',
        fontsize=20)
    # plt.legend(loc='best',numpoints=1)
    # plt.xscale('symlog',linthreshx=1e-14)
    # plt.yscale('log')
    fig = plt.gcf()
    ax1 = fig.add_axes([0.87, 0.10, 0.03, 0.8])

    cmap = cmaps.magma
    norm = matplotlib.colors.Normalize(vmin=vMin, vmax=vMax) 
    
    cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='vertical')
    cb1.set_label(r'Element abruptness, $\gamma(\alpha)$',fontsize=20)
    plt.show()

def measure(CF):
    convergents, tmp = getConvergents(CF)
    total = 0
    for i in range(len(convergents)-2):
        brjuno = np.log(convergents[i+2].denominator)/float(convergents[i].denominator)
        total += brjuno
        # print(convergents[i].denominator,brjuno)
    print(CF, ', Total: ', total)
    return(total)

def plotDiff(base, increase, k=2, symlog = False):
    cmap = cmaps.magma
    colorVal = measure(base)
    increase = int(increase)
    convergents0, alpha0 = getConvergents(base)
    d0, d0s = getDcrit(convergents0,alpha0,k)
    if symlog:
        offset = alpha0
    else:
        offset = 0
    plt.plot(alpha0-offset,d0,'o',markersize=15,color=cmap((colorVal-vMin)/(vMax-vMin)))
    
    ds, alphas = [], []
    for i in range(1,1): #len(base)-1):
        newCF = [x for x in base]
        newCF[i] += increase
        convergents, alpha = getConvergents(newCF)
        d, array = getDcrit(convergents, alpha, k)
        ds.append(d)
        alphas.append(alpha)
        if (i < 6) or (i%4 == 0):
            plt.text(alpha-offset+0.006, d-0.005, str(i),fontsize=15)

    labelStr = r'$M = $' + str(increase+1)

    # plt.plot([x-offset for x in alphas],ds,'go',label=labelStr,markersize=15)
    # plt.plot(alpha0-offset,d0,'k*',markersize=25)
    return d0

def makeDiffPlots():
#     alpha = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

#    for i in range(2,3):
#       plotDiff(alpha,i)
    nmax, length = 5,5
    iteration = itertools.product(np.arange(1,nmax+1), repeat=length)
    for x in iteration:
        # maxEl = max(x)
        # plt.subplot(5,1,maxEl)
        d0 = plotDiff(x, 0)
    setDiffsPlot(x, d0, False)
    
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
    CFs = [ [4,4,4,4,4,4,4],
            [1,2,3,4,3,2,1],
            [1,2,3,4,1,1,1],
            [1,1,1,4,3,2,1],
            [1,1,1,1,1,1,4],
            [4,1,1,1,1,1,1]]

    for CF in CFs:
        convergents, alpha = getConvergents(CF)
        d, array = getDcrit(convergents, alpha, 2)
        print(CF,float(d))



makeDiffPlots()
    
    
