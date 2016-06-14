from pylab import *
import numpy as np

def main(LR=True, plotOn=False):
    fareyNum = [1., 1.] # initial conditions
    fareyDen = [1., 2.] # for left side of tree
    diophantine = [] # edge of relevant excluded region
    fraction = [] # size of remaining region / size between rationals
    d, k = 0.2, 2

    for n in range(2, 45):
        if LR:
            fareyNum.append(np.float64(fareyNum[n-2]+fareyNum[n-1]))
            fareyDen.append(np.float64(fareyDen[n-2]+fareyDen[n-1]))
        else:
            fareyNum.append(np.float64(fareyNum[n-1]+1.))
            fareyDen.append(np.float64(fareyDen[n-1]+3.))
            
    rationals = np.array(fareyNum)/np.array(fareyDen)
    center = np.median(rationals[-5:-1])
    # center = 2/(1+math.sqrt(5))
    print(center)

    # Get Diophantine windows
    for i in range(len(fareyNum)):
        window = np.float64(d/fareyDen[i]**k)
        edge = rationals[i] + window if rationals[i] < center else \
            rationals[i] - window
        diophantine.append(edge)
        if i >= 1:
            wholeSize = np.float64(abs(rationals[i] - rationals[i-1]))
            leftSize = abs(diophantine[i]-diophantine[i-1])
            fraction.append(np.float64(leftSize/wholeSize))
    fraction = np.array(fraction)
    if plotOn:
        figure(1)
        subplot(2,1,1)
        if LR:
            semilogy(abs(fraction-np.median(fraction)))
        else:
            semilogy(abs(fraction-np.median(fraction)),'r')
    diophantine = np.array(diophantine)
    # Plot excluded windows visually
    if plotOn:
        subplot(2,1,2)
        if LR:
            plot(abs(rationals-center),range(len(rationals)),'ko')
            plot(abs(diophantine-center),range(len(rationals)),'bx')
            
        else:
            plot(abs(rationals-center),range(len(rationals)),'ro')
            plot(abs(diophantine-center),range(len(rationals)),'rx')
        xscale('log')
        
main(True,True)
main(False,True)
# main(7,8)

show()