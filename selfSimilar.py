from pylab import *
import numpy as np

def main(path, plotOn=False, d = 0.2, k = 2, col = 'k'):

    class rational:
        def __init__(self,num,den):
            self.num = int(num)
            self.den = int(den)
            self.val = float(num)/float(den)
            self.width = d/self.den**k
            self.diophantine = (self.val - self.width, self.val + self.width)
        
        def __lt__(self,other):
            return self.val < other.val
            
        def __sub__(self,other):
            return self.val - other.val
            
        def __repr__(self):
            return str(self.num)+'/'+str(self.den)

    def getFraction(farey1, farey2):
        separation = abs(farey1 - farey2)
        widthLeft = min(abs(farey1.diophantine[0] - farey2.diophantine[1]), 
                        abs(farey1.diophantine[1] - farey2.diophantine[0]))
        try:
            fraction = widthLeft / separation
        except ZeroDivisionError:
            fraction = 0
        return fraction

    farey = [rational(0,1), rational(1, 1), rational(1,2)]
    fraction = [] # size of remaining region / size between rationals
    maxN = 45
    if type(path) == str:
        path = path*int((maxN+1)/len(path))
        target = None
    else:
        target = path
        path = ''
    
    def nextLevel(farey1, farey2):
        num, den = farey1.num + farey2.num, farey1.den + farey2.den
        return rational(num,den)

    # Generate Farey tree in appropriate path/number
    for n in range(2, maxN):
        if type(target) == (float or int):
            if target > farey[-1].val:
                path = path + 'R' 
            else:
                path = path + 'L'
        if path[n-2] == 'L':
            relevantFarey = [x for x in farey if x < farey[-1]]
        elif path[n-2] == 'R':
            relevantFarey = [x for x in farey if x > farey[-1]]
        winner, minDist = relevantFarey[0], abs(farey[-1]-relevantFarey[0])
        for i in range(1,len(relevantFarey)):
            if abs(farey[-1]-relevantFarey[i]) < minDist:
                winner = relevantFarey[i]
                minDist = abs(farey[-1]-relevantFarey[i])
        
        farey.append(nextLevel(farey[-1],winner))

    if type(target) == (float or int):
        center = target
    else:
        center = np.mean([x.val for x in farey[-5:-1]])
    
    print(center)
    
    # Get irrational fraction between two rationals
    fraction = []
    for i in range(1,len(farey)):
        fraction.append(getFraction(farey[i-1],farey[i]))
    fraction = np.asarray(fraction)
    
    if plotOn:
        figure(1)
        subplot(2,1,1)
        plot(fraction,color=col)
        text(3,1.02,'Remaining fraction between two rational excluded regions')
        xticks([])
        ylabel(r"$f_L$")
        
        subplot(2,1,2)
        semilogy(abs(fraction-np.median(fraction)),color=col)
        # legend(['Alternating noble steps',r"Steps approaching 1/3"],loc='best')
        xlabel('Steps down Farey tree')
        ylabel(r"Error from asymptotic $f_L$")
        subplots_adjust(hspace=0)
        
    # Plot excluded windows visually
    if plotOn:
        figure(2)
        for i in range(len(farey)):
            plot(i,abs(farey[i].val-center),color=col,marker='o')
            plot(i,min(abs(farey[i].diophantine[0]-center),abs(farey[i].diophantine[1]-center)),
                    color=col,marker='x')
        yscale('log')
        legend(['Rational','Diophantine window edge'],loc='best')
        ylabel('Distance from asymptotic limit')
        xlabel('Steps down Farey tree')

dmin, dmax = 0.15, 0.35
cmap = cm.get_cmap('gnuplot2')
junk = contourf([[0,0],[0,0]],np.arange(dmin,dmax+0.01,0.001),cmap=cmap)
clf()

# for d in np.arange(dmin, dmax, 0.01): 
#    main(True,1,True, d, cmap((d-dmin)/(dmax-dmin)))
main('RLLL',True)
main('RLL',True)
main('RL',True)
# main(1/math.pi,True)
fig = figure(1)
fig.subplots_adjust(right=0.8)
cbarAx = fig.add_axes([0.85,0.15,0.05,0.7])
colorbar(junk,cax=cbarAx)
show()