from pylab import *
import numpy as np
from fractal import getP

nmax = 100
drange = np.arange(0,.51,0.005)
krange = np.arange(1, 4.02,0.01)
lenD, lenK = len(drange), len(krange)

storeP = np.zeros((lenD,lenK,nmax))
clf()
cmap = get_cmap('RdBu')
cmap.set_under('black')
cmap.set_over('green')

for n in range(1,nmax):
    for d in range(lenD):
        for k in range(lenK):
            storeP[d,k,n-1] = max(getP(drange[d],krange[k],n)[1])
            
    pcolor(krange,drange,storeP[:,:,n-1],cmap=cmap,vmin=0.001,vmax=.999)
    title('Max pressure at Farey level '+str(n))
    colorbar(extend='both')
    ylim(0,.5)
    xlim(1,4)
    savefig('C:/Users/bkraus/Dropbox/Hudson/Images_dkMaps_020516/dk'+str(n)+'.jpg')
    clf()