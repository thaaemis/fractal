from pylab import *
import numpy as np

def main(LR=True, a = 1, plotOn=False):
    fareyNum = [1., 1.] # initial conditions
    fareyDen = [1., 2.] # for left side of tree
    diophantine = [] # edge of relevant excluded region
    fraction = [] # size of remaining region / size between rationals
    d, k = 0.2, 2

    for n in range(2, 100):
        if LR:
            fareyNum.append(np.float64(fareyNum[n-2]+fareyNum[n-1]))
            fareyDen.append(np.float64(fareyDen[n-2]+fareyDen[n-1]))
        else:
            fareyNum.append(np.float64(fareyNum[n-1]+1.))
            fareyDen.append(np.float64(fareyDen[n-1]+np.float64(a)))
            
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
    # if plotOn:
        # figure(1)
        # subplot(2,1,1)
        # if LR:
            # plot(fraction)
        # else:
            # plot(fraction,'r')
        # text(3,1.02,'Remaining fraction between two rational excluded regions')
        # xticks([])
        # ylabel(r"$f_L$")
        
        # subplot(2,1,2)
        # if LR:
            # semilogy(abs(fraction-np.mean(fraction)))
        # else:
            # semilogy(abs(fraction-np.mean(fraction)),'r')
        # legend(['Alternating noble steps',r"Steps approaching 1/3"],loc='best')
        # xlabel('Steps down Farey tree')
        # ylabel(r"Error from median $f_L$")
        # subplots_adjust(hspace=0)
        
    diophantine = np.array(diophantine)
    # Plot excluded windows visually
    if plotOn:
        figure(2)
        if LR:
            plot(abs(rationals-center),range(len(rationals)),'ko')
            plot(abs(diophantine-center),range(len(rationals)),'bx')
            
        else:
            plot(abs(rationals-center),range(len(rationals)),'ro')
            plot(abs(diophantine-center),range(len(rationals)),'rx')
        xscale('log')
        legend(['Rational center','Diophantine window edge'],loc='best')
        xlabel('Distance from asymptotic limit')
        ylabel('Steps down Farey tree')
        
main(True,1,True)
main(False,1,True)
main(False,2,True)
main(False,3,True)

# main(7,8)

show()