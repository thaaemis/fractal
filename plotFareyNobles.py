# To do in this code:
# put all PLOTTING functions in Tkinter module description.Regularize everything.s
# Plot convergent lines VS. direct Farey lines 
# How to finally fix color normalization all the way through (L/R asymmetry)

import math, itertools, matplotlib
matplotlib.use("QT4Agg")
import matplotlib.pyplot as plt
import numpy as np
from fractal import getP, rationalList
from selfSimilar import rational, makeFractionList, getFareyPath
# from Tkinter import *
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# from matplotlib.figure import Figure
import sys
sys.path.insert(0, '/home/brian/GitHub/colormap/')
import colormaps as cmaps
import pickle

phi = (math.sqrt(5)-1)/2.
# with open('fareySequences.pkl','r') as f:
#    fareySequences = pickle.load(f)

def CFtoVal(CFarray):
    CF = [x for x in CFarray]
    val = 1/float(CF.pop())
    for i in range(len(CF)):
        element = CF.pop()
        val = 1/(float(element)+val)
    return val

def valtoCF(val, stepsMax = 1E5):
    CF = []
    while True:
        CF.append(math.floor(1/val))
        val = 1/val - CF[-1]
        if (abs(val) < 1e-10) or (len(CF) > stepsMax):
            break
    return CF

def Brjuno(fareyList, j):
    fareyList = fareyList[1:]
    val, n = 0, 0
    while j+n < len(fareyList):
        val += math.log(fareyList[n+j].den)/fareyList[n].den
        n += 1
    return val

    

def plotLine(CF,axisoffset=0., d = 0.15,k = 2., 
    cmap = matplotlib.cm.ScalarMappable(cmap=cmaps.viridis), axisLog = False):

    def getNextRational(element, fareyList, d, k):
        n = element * fareyList[-1].num + fareyList[-2].num
        m = element * fareyList[-1].den + fareyList[-2].den
        return rational(n, m, d, k)
        
        
    CF = [float(c) for c in CF]
    farey = [rational(1,0,d,k), rational(0,1,d,k)] # initiate
    if (CF[-1] != 1.0):
        CF[-1] = float(CF[-1])-1
        CF.append(1.)

    neighbors = []
    for i in range(len(CF)):
        a = getNextRational(int(CF[i]),farey,d,k)
        farey.append(a)

    farey = farey[1:]
#    for i in range(len(farey)-1):
#        level = 0
#        while True:
#            fareySequence = fareySequences[level]
#            fareyNum = [x[0] for x in fareySequence]
#            fareyDen = [x[1] for x in fareySequence]
#            ind = []
#            for f in [farey[i],farey[i+1]]:
#                indNum = np.where(np.asarray(fareyNum) == f.num)[0]
#                indDen = np.where(np.asarray(fareyDen) == f.den)[0]
#                if (len(list(np.intersect1d(indNum, indDen))) == 0):
#                    level += 1
#                    continue
#                else:
#                    ind.append(list(np.intersect1d(indNum, indDen))[0])
#            if len(ind) == 2:
#                neighboring = True if abs(ind[0]-ind[1]) == 1 else False
#                neighbors.append(neighboring)
#                break
    plt.plot([d-axisoffset,-d-axisoffset],
            [farey[0].den]*2, color = cmap.to_rgba(1), linewidth=4)
    for j in range(len(farey)-1):
        changed = False
        if CF[0] > 1:
            CF[0] -= 1
            changed = True
        maxBound = max(CF[:j]) if j > 0 else 1
        if changed:
            CF[0] += 1
        color = cmap.to_rgba(maxBound) # if cmap.norm == None else cmap(maxBound/10)
        # plot line only if partners are Farey neighbors
        # if neighbors[j]:
        plt.plot([i.val-axisoffset for i in farey[j:j+2]],
                [i.den for i in farey[j:j+2]],color=color,linewidth=1)
        (diophantineMin, diophantineMax) = farey[j+1].diophantine
        plt.plot([diophantineMin-axisoffset,diophantineMax-axisoffset],
                [farey[j+1].den]*2,
                color = color, linewidth=4)
    
    # print(CF, farey, neighbors)
    # print(Brjuno(farey,1), Brjuno(farey,2), Brjuno(farey,3))
    
    return farey
        
goldenMean = valtoCF((math.sqrt(5)-1)/2)

def makeAlltoN(d=0.2,kHere=2., nmax=3, length=2, axisOffset = 0, 
        colorMap=matplotlib.cm.ScalarMappable(norm=matplotlib.colors.LogNorm(vmin = 1, 
        vmax = 10), cmap=cmaps.viridis)):

    # plt.imshow([[-5,-2],[-1,1]],norm = norm, extent = [-1, -0.5, 1000e3, 10001e3])
    # ax1 = plt.colorbar(label='Max Element in CF')
    # colorMap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmaps.viridis)
    
    iteration = itertools.product(np.arange(1,nmax+1), repeat=length)
    els, fareyVals = [], []
    for i in reversed(list(iteration)):
        elements = [float(g) for g in i]
        els.append(elements)
        fareyVals.append(plotLine(elements,d=d,k=kHere,cmap = colorMap, 
            axisoffset = axisOffset))
        elements.insert(0,1.0)
        els.append(elements)
        fareyVals.append(plotLine(elements,d=d,k=kHere,cmap = colorMap, 
            axisoffset = axisOffset))
        elements = elements[1:]
        elements[0] = float(nmax+1)
        els.append(elements)
        fareyVals.append(plotLine(elements,d=d,k=kHere,cmap = colorMap, 
            axisoffset = axisOffset))

def plotAssortment(offset = 0, d = 0.15, k = 2, path = None, axisLog = False):
    if axisLog == True:
        offset = (math.sqrt(5)-1)/2.
        
    if path == None:
        plt.clf()
        norm = matplotlib.colors.LogNorm(vmin = 2, vmax = 20)
        plt.imshow([[-5,-2],[-1,1]],norm = norm, extent = [-1, -0.5, 1000e3, 10001e3])
           
        ax1 = plt.colorbar(label='Max Element in CF')
        
        levels = 18
        for i in range(levels):
            a = [1]*levels
            a[i] = 2
            plotLine(a,axisoffset=offset, d = d, k = k, axisLog = axisLog)
            b = [1]*levels
            b[i] = 3
            plotLine(b,axisoffset=offset, d = d, k = k, axisLog = axisLog)
            c = [1,2,1] * (levels/3)
            c[i] = 3
            plotLine(c,axisoffset=offset, d = d, k = k, axisLog = axisLog)
            f = [1] * levels
            f[i] = 4
            plotLine(f,axisoffset=offset, d = d, k = k, axisLog = axisLog)
            e = [2] * levels
            e[i] = 4
            plotLine(e,axisoffset=offset, d = d, k = k, axisLog = axisLog)
        plotLine([1,15,1,1,1,1,2,1,1,1],axisoffset=offset, d = d, k = k, axisLog = axisLog)
        plotLine([1]*100,axisoffset=offset, d = d, k = k, axisLog = axisLog)
    else:
        plotLine(path, axisoffset = offset, d = d, k = k, axisLog = axisLog)

def restartPlot(restart = True,symlog = False):
    
    fig1 = plt.figure(1,figsize=(10,10))
    if restart:
        plt.clf()
    norm = matplotlib.colors.LogNorm(vmin = 1, vmax = 2)
    cMap = cmaps.viridis
    cMap.set_under('r')
    colorMap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cMap)
    
    # plt.imshow([[-5,-2],[-1,1]],norm = norm, extent = [-1, -0.5, 1e-7, 2e-7],cmap=cmaps.viridis)
    # ax1 = plt.colorbar(label='Max Element in CF',extend='both')
    
    #legend 
    plt.plot([-1,-1],[1e-3,1e-3],'k-',linewidth=1,label="Farey steps")
    plt.plot([-2,-1],[1e-3,1e-3],'k-',linewidth=4,label="Diophantine widths")
    plt.legend(loc=1)
    
    axisOffset = 0
    plt.ylabel('Denominator',fontsize=20)
    plt.yscale('log',basey=10)
    plt.xlabel('Value',fontsize=20)
    plt.axis([0,1,3e-1,1e3])
    plt.gca().invert_yaxis()
    
    return fig1, colorMap, axisOffset


restartPlot(True, symlog=True)
CF = [1]*30
# makeAlltoN(axisOffset = phi, nmax = 3, length = 6, d = 1-phi, kHere = 1.999)
plotLine(CF,axisoffset=phi, d = 1-phi,k = 2.)
CF = [1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
plotLine(CF,axisoffset=CFtoVal(CF), d = 0.3828,k = 2.)
plt.xlabel(r'Value - $\phi$',fontsize=20)
plt.xscale('symlog',linthreshx=1e-7)
plt.axis([-1,1,3e-1,1e3])

plt.show()








def module():

    fig1, colorMap, axisOffset = restartPlot(True)
        
    top = Tk()
    allLabels = Label(top,text="Insert parameters d and k")
    allLabels.grid(row=0, column=0,
                rowspan=2)

    top.title('Farey Tree Generator')
    top.minsize(600,600)
    
    dInput = Entry(top,width=10)
    dInput.grid(row=0,column=1)
    kInput = Entry(top,width=10)
    kInput.grid(row=1,column=1)

    dInput.insert(0, 0.15)

    kInput.insert(0,2.0)
    
    def makeLog():
        plt.yscale('log')
        plt.gca().invert_yaxis()
        plt.xscale('symlog',linthreshx = 1e-7)
        axisOffset = (math.sqrt(5)-1)/2.

    def makeLin():
        plt.yscale('log')
        plt.gca().invert_yaxis()
        plt.xscale('linear')
        axisOffset = 0


    axisLogSubmit = Button(top, text="Xscale: symlog", command=makeLog)
    axisLogSubmit.grid(row=7,column=1)
    axisLinSubmit = Button(top, text="Xscale: linear", command=makeLog)
    axisLinSubmit.grid(row=7,column=2)

        
    makeButton = Button(top, text="Generate many values",
        command = lambda: plotAssortment(d=float(dInput.get()), 
                                         k=float(kInput.get()),
                                         axisLog=axisLogBool.get()))
    # makeButton.grid(row=0,column=2, rowspan=2)

    dividerLabel1 = Label(top,text="  ")
    dividerLabel1.grid(row=2,column=0)

    dividerLabel = Label(top,text="Insert continued fraction: ")
    dividerLabel.grid(row=3,column=0)
    CFinput = Entry(top, width=30)
    CFinput.grid(row=3,column=1)
    submitPathButton = Button(top,text="Add Value", command=
        lambda: plotLine(CFinput.get().split(','),
            d = float(dInput.get()), k=float(kInput.get()),
            cmap = colorMap,axisoffset=axisOffset))
    submitPathButton.grid(row=3,column=2)
    
    dividerLabel2 = Label(top,text="  ")
    dividerLabel2.grid(row=4,column=0)

    
    treeLabel = Label(top, text="Add tree for nMax, length: ")
    treeLabel.grid(row=5,column=0)
    nmaxInput = Entry(top,width = 10)
    nmaxInput.grid(row=5,column=1)
    lengthInput = Entry(top, width = 10)
    lengthInput.grid(row=5,column=2)
    submitTreeButton = Button(top, text="Add Tree",
        command=lambda: makeAlltoN(d = float(dInput.get()),
            kHere = float(kInput.get()), nmax = int(nmaxInput.get()),
            length = int(lengthInput.get()),axisOffset = axisOffset,
            colorMap = colorMap))
    submitTreeButton.grid(row=6, column=0,columnspan=3)

    
    
    
    

    # canvas = FigureCanvasTkAgg(fig1, top)
    # toolbar = NavigationToolbar2TkAgg(canvas, top)
    
    def includeFig():
        # canvas.get_tk_widget().delete("all")
        canvas.show()
        canvas.get_tk_widget().grid(row=8,column=0,columnspan=3)


        # toolbar.update()
        # canvas._tkcanvas.grid(row=4,column=0,columnspan=3)

    showButton = Button(top, text="Show",
        command = includeFig)
    showButton.grid(row=9,column=0,columnspan=3)
        
    Button(top, text="Quit", command = top.quit).grid(row=9,column=1,columnspan=3)
    top.mainloop()

# module()

#plt.figure(1)
#d, k = 0.2, 1.9
#plt.subplot(2,1,1)
#plotAssortment(0, d = d, k = k)
#plt.axis([0,1,1,1e3])
#plt.xlabel('Value',fontsize=20)
#plt.ylabel('Denominator',fontsize=20)
#    
#plt.subplot(2,1,2)
#plotAssortment(phi, d = d, k = k)
#plt.xscale('symLog',linthreshx=1e-7)
#plt.xlabel(r'Value - $\varphi$',fontsize=20)
#plt.ylabel('Denominator',fontsize=20)
#    

#plt.axis([-0.5,0.5,1,1e3])


#plt.figure(2)
#d, k = 0.2, 2.1
#plt.subplot(2,1,1)
#plotAssortment(0, d = d, k = k)
#plt.axis([0,1,1,1e3])
#plt.xlabel('Value',fontsize=20)
#plt.ylabel('Denominator',fontsize=20)
#    
#plt.subplot(2,1,2)
#plotAssortment(phi, d = d, k = k)
#plt.xscale('symLog',linthreshx=1e-7)
#plt.axis([-0.5,0.5,1,1e3])
#plt.xlabel(r'Value - $\varphi$',fontsize=20)
#plt.ylabel('Denominator',fontsize=20)
#plt.show()
