import matplotlib, math, itertools
import matplotlib.pyplot as plt
import numpy as np
from fractal import getP
from selfSimilar import rational, makeFractionList, getFareyPath
from Tkinter import *
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

phi = (math.sqrt(5)-1)/2.

def CFtoVal(CFarray):
    CF = [x for x in CFarray]
    val = 1/float(CF.pop())
    for i in range(len(CF)):
        element = CF.pop()
        val = 1/(float(element)+val)
    return val

def valtoCF(val, stepsMax = 100):
    CF = []
    while True:
        CF.append(math.floor(1/val))
        val = 1/val - CF[-1]
        if (abs(val) < 1e-10) or (len(CF) > stepsMax):
            break
    return CF

def plotLine(value,axisoffset=0, d = 0.15,k = 2, cmap = plt.get_cmap('spectral')):
    try:
        len(value)
        value = CFtoVal(value)
    except TypeError:
        value = float(value)
        
    maxN = 100
    convergents = False
    farey, center, pathOut = getFareyPath(value, maxN, d, k, convergents)
    
    for j in range(len(farey)-1):
        maxBound = max(valtoCF(farey[j+1].val))
        color = cmap(maxBound/20)
        plt.plot([i.val-axisoffset for i in farey[j:j+2]],
                [i.den for i in farey[j:j+2]],color=color)
        (diophantineMin, diophantineMax) = farey[j+1].diophantine
        plt.plot([diophantineMin-axisoffset,diophantineMax-axisoffset],
                [farey[j+1].den]*2,
                color = color, linewidth=4)
        plt.yscale('log',basey=2)
        plt.axis([0,1,0.9,1e3])
        
goldenMean = valtoCF((math.sqrt(5)-1)/2)
# print([index for index,value in enumerate(goldenMean) if value > 1])

def makeAlltoN(d=0.2,k=2., nmax=10, length=10):
    lists = itertools.combinations_with_replacement(range(1,nmax+1),length)
    for i in lists:
        elements = [float(k) for k in i]
        plotLine(CFtoVal(elements),d=d,k=k)
        
def plotAssortment(offset = 0, d = 0.15, k = 2, path = None):
    if path == None:
        plt.clf()
        norm = matplotlib.colors.LogNorm(vmin = 1, vmax = 20)
        plt.imshow([[-5,-2],[-1,1]],norm = norm, extent = [-1, -0.5, 1000e3, 10001e3])
           
        ax1 = plt.colorbar(label='Max Element in CF')
        
        levels = 18
        for i in range(levels):
            a = [1]*levels
            a[i] = 2
            plotLine(a,axisoffset=offset, d = d, k = k)
            b = [1]*levels
            b[i] = 3
            plotLine(b,axisoffset=offset, d = d, k = k)
            c = [1,2,1] * (levels/3)
            c[i] = 3
            plotLine(c,axisoffset=offset, d = d, k = k)
            f = [1] * levels
            f[i] = 4
            plotLine(f,axisoffset=offset, d = d, k = k)
            e = [2] * levels
            e[i] = 4
            plotLine(e,axisoffset=offset, d = d, k = k)
        plotLine([1,15,1,1,1,1,2,1,1,1],axisoffset=offset, d = d, k = k)
        plotLine([1]*100,axisoffset=offset, d = d, k = k)
    else:
        plotLine(path, d = d, k = k)
def module():

    
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

        
    makeButton = Button(top, text="Generate many values",
        command = lambda: plotAssortment(d=float(dInput.get()), 
                                         k=float(kInput.get())))
    makeButton.grid(row=0,column=2, rowspan=2)

    
    dividerLabel = Label(top,text="Insert continued fraction: ")
    dividerLabel.grid(row=2,column=0)
    CFinput = Entry(top, width=30)
    CFinput.grid(row=2,column=1)
    submitPathButton = Button(top,text="Add Value", command=
        lambda: plotLine(CFinput.get().split(','),
            d = float(dInput.get()), k=float(kInput.get())))
    submitPathButton.grid(row=2,column=2)
    
    treeLabel = Label(top, text="Add tree for nMax, length: ")
    treeLabel.grid(row=3,column=0)
    nmaxInput = Entry(top,width = 10)
    nmaxInput.grid(row=3,column=1)
    lengthInput = Entry(top, width = 10)
    lengthInput.grid(row=3,column=2)
    submitTreeButton = Button(top, text="Add Tree",
        command=lambda: makeAlltoN(d = float(dInput.get()),
            k = float(kInput.get()), nmax = int(nmaxInput.get()),
            length = int(lengthInput.get())))
    submitTreeButton.grid(row=4, column=0,columnspan=3)

    fig1 = plt.figure(1)
    canvas = FigureCanvasTkAgg(fig1, top)
    # toolbar = NavigationToolbar2TkAgg(canvas, top)
    
    def includeFig():
        canvas.get_tk_widget().delete("all")
        canvas.show()
        canvas.get_tk_widget().grid(row=6,column=0,columnspan=3)


        # toolbar.update()
        # canvas._tkcanvas.grid(row=4,column=0,columnspan=3)

    showButton = Button(top, text="Show",
        command = includeFig)
    showButton.grid(row=7,column=0,columnspan=3)
        
    Button(top, text="Quit", command = top.quit).grid(row=8,column=0,columnspan=3)
    top.mainloop()

module()

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
