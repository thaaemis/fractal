from pylab import *
import numpy as np
import pandas as pd
from scipy.ndimage.filters import gaussian_filter1d
from scipy.interpolate import interp1d
from savitzky_golay import savitzky_golay
from fractal import getP, cylinderB
import pickle

def getProfiles(nmax = 4, smoothWindow = 51):
    d, k, fareyMethod = 0.1, 2.2, 'treeSteps'
    pickleFile = 'standardPressure%d.pkl' % nmax
    try:
        with open(pickleFile, 'r') as f:
            x, p, gradp = pickle.load(f)
    except IOError:
        x, p, gradp = getP(d, k, nmax, fareyMethod)
        with open(pickleFile, 'w') as f:
            pickle.dump((x, p, gradp), f)
    # plot(x, p, 'b-')
    # show()
    p, gradp = list(p), list(gradp)
    p.reverse()
    gradp.reverse()
    p, gradp = np.asarray(p), np.asarray(gradp)
    pInterp = interp1d(x, p)
    xDiffs = np.diff(x)
    positiveDiffs = [element for element in xDiffs if element != 0]
    minxDiff = np.min(np.abs(positiveDiffs))
    print(minxDiff)
    denseRange = np.arange(min(x), max(x), minxDiff)

    evenGridP = pInterp(denseRange)
    newx = denseRange
    newp = savitzky_golay(evenGridP, smoothWindow, 2)
    newgradp = np.gradient(newp)/np.gradient(newx)
    # plot(x, p, 'g-', denseRange, newp, 'r-')
    # show()

    # input('')
    oldr, oldBz, oldBtheta, oldJtheta, oldJz, x1, p1, gradp1 = cylinderB(d, k, nmax, 1e-5, 1., 
                                                       fareyMethod, list(x), list(p), list(gradp))

    # subplot(311)
    # plot(x, p, 'b-', newx, newp, 'r-', markersize=1)
    newx, newp, newgradp = list(newx), list(newp), list(newgradp)
    newx.reverse()
    newp.reverse()
    newgradp.reverse()
    newr, newBz, newBtheta, newJtheta, newJz, tmp, newp, newgradp = cylinderB(d, k, nmax, 1e-5, 1., 
                                                      fareyMethod, newx, newp, newgradp)

    oldProfiles = (oldr, oldBz, oldBtheta, oldJz, oldJtheta, x, p, gradp)
    newProfiles = (newr, newBz, newBtheta, newJz, newJtheta, newx, newp, newgradp)
    profilePickle = 'profiles-n%d-order2.pkl' % nmax
    with open(profilePickle, 'w') as f:
        pickle.dump((oldProfiles, newProfiles), f)

    return oldProfiles, newProfiles

def comparePlots(nmax = 4):
    
    try:
        profilePickle = 'profiles-n%d-order2.pkl' % nmax
        with open(profilePickle, 'r') as f:
            oldProfiles, newProfiles = pickle.load(f)
    except IOError:
        oldProfiles, newProfiles = getProfiles(nmax, smoothWindow = 51)

    for profiles, color, label in [(oldProfiles, 'b', 'Ideal'), (newProfiles, 'r', 'Smoothed')]:
        r, Bz, Btheta, Jz, Jtheta, x, p, gradp = profiles
        
        subplot(211)
        plot(x, max(p)-p, color + '-', markersize=1, label = label)
        subplot(212)
        plot(r, Jtheta, color + '.', markersize=.5)
    subplot(211)    
    text(0.1, 0.4, r'$p(r)$')
        # plot(newx, max(newp) - newp, 'r-', markersize=1, label = 'Smoothed')
    legend(loc='best')
    gca().get_xaxis().set_visible(False)
    subplot(212)
    text(0.1, 1.0, r'$J_\theta(r)$')
    xlabel('Radius')

    show()

comparePlots(8)
