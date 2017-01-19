from pylab import *
import numpy as np
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.integrate import quad
from scipy.interpolate import interp1d

from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.colors as colors

from fractal import getP

def gridPlots3D():
    kList = np.arange(1.3, 2.5, 0.1)
    iList = range(1,9)

    # Mesh in 3d
    kMesh, iMesh = np.meshgrid(kList, iList)
    zMesh = kMesh*0

    # more info on grid spacing
    avgSpacingMesh = kMesh*0
    minSpacingMesh = kMesh*0

    #triangle list
    kTri, iTri, zTri = [], [], []
    minSpacingTri, avgSpacingTri, colTri = [], [], []
    for kEl, k in enumerate(kList): 
        for iEl, i in enumerate(iList):
            x, p, gradP = getP(0.1, k, i, fareyMethod = 'treeSteps')
            
            print('k=', str(k), 'Farey ', i, len(x))
            print(kEl, iEl)
            zMesh[iEl, kEl] = len(x)
            kTri.append(k)
            iTri.append(i)
            zTri.append(len(x)/8./2.**i)
            
            distTmp, dist = np.diff(x), []
            for j in distTmp:
                if j != 0.:
                    dist.append(j)
            minSpacingMesh[iEl, kEl] = np.max(dist)*-1
            minSpacingTri.append(np.max(dist)*-1)
            avgSpacingMesh[iEl, kEl] = np.mean(dist)*-1
            avgSpacingTri.append(np.mean(dist)*-1)
            if np.mean(dist) == -1:
                colTri.append([1,0,0,0])
            else:
                colTri.append([0,0,1,0])
            
    # figure 1
    fig = figure()
    ax = fig.add_subplot(111, projection = '3d') 
    # surface = ax.plot_surface(kMesh, iMesh, zMesh)
    surface = ax.plot_trisurf(kTri, iTri, zTri)
    ax.set_xlabel('k',fontsize=20)
    ax.set_ylabel('Farey Level', fontsize=20)
    ax.set_zlabel('Number gridpoints',fontsize=20)

    #figure 2
    fig2 = figure(2)
    ax2 = fig2.add_subplot(111, projection = '3d') 
    surface2 = ax2.plot_trisurf(kMesh, iMesh, np.log10(avgSpacingTri))
    ax2.set_xlabel('k',fontsize=20)
    ax2.set_ylabel('Farey Level', fontsize=20)
    ax2.set_zlabel('Mean grid spacing',fontsize=20)
    show()

def showGrid():
    maxI = 11
    for i in range(1,maxI):
        x, p, gradP = getP(0.1, 2.1, maxI - i, fareyMethod = 'treeSteps')
        plot(x,np.asarray(x)*0.+i,'k|',markersize=15,linewidth=4)
        print(i, len(x))
    axis([-0.01,1.01,0.5,maxI-0.5])
    show()
    
def showNewGrid():
    maxI = 11
    xOld = []
    for i in range(1,maxI):
        x, p, gradP = getP(0.1, 2.1, i, fareyMethod = 'treeSteps')
        for xi in x:
            col = 'k' if xi in xOld else 'r'
            plot(xi,maxI-i,'|',color=col,markersize=15,linewidth=4,mew=3)
        xOld = x
        print(i, len(x))
    axis([-0.01,1.01,0.5,maxI-0.5])
    show()
    
showNewGrid()

