# Generate fractal diagrams with proper r

import numpy as np
from pylab import *
from fractal import cylinderB, getP

def plotPressures():
    # iota, Bz, Bth, Jth, Jz, x, p, gradp = cylinderB(0.1,2.1,5,.00001,
    #     R=3.,fareyMethod = 'treeSteps')

    x1, p1, gradp1 = getP(0.1, 2.1, 10, fareyMethod = 'treeSteps', getN = False)
    x2, p2, gradp2 = getP(0.35, 2.1, 10, fareyMethod = 'treeSteps', getN = False)
    x3, p3, gradp3 = getP(0.38, 2.1, 10, fareyMethod = 'treeSteps', getN = False)

    fig = figure()
    plot(x1, [p/max(p1) for p in p1], 'r', label=r'$d = 0.25$')
    plot(x2, [p/max(p2) for p in p2], 'b', label=r'$d = 0.35$')
    plot(x3, [p/max(p3) for p in p3], 'g', label=r'$d = 0.38$')
    xlabel(r'$r/a$', fontsize=30)
    ylabel(r'$p / \max(p)$', fontsize=30)

    xticks(fontsize=20)
    yticks(fontsize=20)
    legend(loc='best',fontsize=20)
    text(0.6, 0.2, r'$k = 2.1,$ Farey level 10', fontsize=30)
    axis([0,1.1,-0.05,1.05])
    show()

    # subplot(212)
    # plot(r, Jz, r, Jth)
    # show()
    
def plotB():
    iota, Bz, Bth, Jth, Jz, x, p, gradp = cylinderB(0.1,2.1,5,.00001,
        R=3.,fareyMethod = 'treeSteps')
        
    plot(iota, Jz, iota, Jth)
    show()
    
plotB()
