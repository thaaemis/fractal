from pylab import *
import time
import matplotlib as mpl
from fractal import cylinderB
import numpy as np
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import Axes3D
sys.path.insert(0, '/home/brian/GitHub/colormap/')
import colormaps as cmaps

r, Bz, Btheta, Jtheta, Jz, x, p, gradp = cylinderB(.15, 2, 8, 0.00001, 
        fareyMethod='treeSteps')

gridSpace = 0.02
xMesh, yMesh = np.meshgrid(np.arange(0,1,gridSpace),  \
                           np.arange(0,1,gridSpace))

JzFxn = interp1d(r, Jz)
JthetaFxn = interp1d(r, Jtheta)
BzFxn = interp1d(r, Bz)
BthetaFxn = interp1d(r, Btheta)
pFxn = interp1d(x,p)
gradPFxn = interp1d(x, gradp)

JzMesh = [JzFxn(x) for x in yMesh[:,:]]
JthetaMesh = [JthetaFxn(x) for x in yMesh[:,]]
BzMesh = [BzFxn(x) for x in yMesh[:,:]]
BthetaMesh = [BthetaFxn(x) for x in yMesh[:,]]

maxJz, minJz = np.max(Jz), np.min(Jz)
cmap = cmaps.viridis
norm = mpl.colors.Normalize(vmin = minJz, vmax = maxJz)
thetaMesh, rMesh = np.meshgrid(np.arange(0,3*math.pi/2.,gridSpace),
                               np.arange(0,1,gridSpace))
pMesh = [pFxn(x) for x in rMesh]
gradPMesh = [gradPFxn(x) for x in rMesh]
side1xMesh = rMesh*0+1.
side1yMesh = rMesh * np.cos(thetaMesh)
side1zMesh = -1 * rMesh * np.sin(thetaMesh)

color2Mesh = rMesh * 0 + gradPFxn(0)
side2rMesh = 0 * rMesh + 1
side2thetaMesh = thetaMesh
side2xMesh = rMesh
side2yMesh = side2rMesh * np.cos(side2thetaMesh)
side2zMesh = side2rMesh * -1 * np.sin(side2thetaMesh)

fig = figure()
ax = fig.add_subplot(111, projection = '3d')
color = (JzMesh - minJz)/(maxJz - minJz)

# close cylinder
print(size(side2xMesh), size(side2yMesh), size(side2zMesh), size(color2Mesh))
time.sleep(1)
ax.scatter(side2xMesh, side2yMesh, side2zMesh, 
    c = color2Mesh, cmap = cmap, linewidth = 0)# Jz on z = 0 plane
ax.scatter(xMesh, yMesh, yMesh*0., 
    c = JzMesh, cmap = cmap, linewidth = 0)
# Jtheta on y = 0
ax.scatter(xMesh, yMesh*0, yMesh, 
    c = JthetaMesh, cmap = cmap, linewidth = 0)
# pressure on x = 1
print(size(side1xMesh), size(side1yMesh), size(side1zMesh), size(pMesh))
time.sleep(1)
ax.scatter(side1xMesh, side1yMesh, side1zMesh, 
    c = gradPMesh, cmap = cmap, linewidth = 0)

ax.set_ylabel('y')
ax.set_xlabel('x')
ax.set_zlabel('z')
ax.view_init(elev = 25, azim = 40)
show()


