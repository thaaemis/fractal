from pylab import *
import numpy as np

class rational:
    def __init__(self,num,den,d,k):
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

def getFareyPath(path,maxN,d,k):    
    farey = [rational(0,1,d,k), rational(1, 1,d,k), rational(1,2,d,k)]
    if type(path) == str:
        path = path*int((maxN)/len(path)+1)
        target = None
    else:
        target = path
        path = ''
    def nextLevel(farey1, farey2):
        num, den = farey1.num + farey2.num, farey1.den + farey2.den
        return rational(num,den,d,k)

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
        
    return(farey, center, path)


def makeFractionList(path, plotOn=False, d = 0.1, k = 2, col = 'k'):

    def getFraction(fareyList):
        farey1, farey2 = fareyList.pop(0),fareyList.pop()
        separation = abs(farey1 - farey2)
        widthLeft = min(abs(farey1.diophantine[0] - farey2.diophantine[1]), 
                        abs(farey1.diophantine[1] - farey2.diophantine[0]))
        while len(fareyList) > 0:
            next = fareyList.pop()
            widthLeft = widthLeft - abs(next.diophantine[1]-next.diophantine[0])
        if separation > 1e-15:
            fraction = widthLeft / separation
        else:
            fraction = 0
        return fraction
    
    pathOriginal = path
    farey, center, path = getFareyPath(path, 200, d, k)
    
    # Get irrational fraction between two rationals
    rawFraction, fraction, fracOrder = [], [], []
    skip = 1 if type(pathOriginal) == (int or float) else len(pathOriginal)-1
    for i in range(1,len(farey)):
        if skip > 1 and i > skip:
            rawFraction.append(getFraction(farey[i-skip:i+1]))
        else:
            rawFraction.append(getFraction(farey[i-1:i+1]))
    for i in range(len(rawFraction)):
        if rawFraction[i] > 0:
            fracOrder.append(i)
            fraction.append(rawFraction[i])
    fraction = np.asarray(fraction)
    
    if plotOn:
        figure(1)
        subplot(2,1,1)
        plot(fracOrder,fraction,color=col,marker='o')
        text(3,1.02,'Remaining fraction between two rational excluded regions')
        xticks([])
        ylabel(r"$f_L$")
        
        subplot(2,1,2)
        semilogy(fracOrder,abs(fraction-np.median(fraction)),color=col,marker='o') 
        # legend(['Alternating noble steps',r"Steps approaching 1/3"],loc='best')
        xlabel('Steps down Farey tree')
        ylabel(r"Error from asymptotic $f_L$")
        subplots_adjust(hspace=0)
        
    # Plot excluded windows visually
    if plotOn:
        figure(2)
        for i in range(len(farey)):
            plot(i,abs(farey[i].val-center),color=col,marker='o')
            # plot(i,min(abs(farey[i].diophantine[0]-center),abs(farey[i].diophantine[1]-center)),
            #       color=col,marker='x')
        yscale('log')
        # legend(['Rational','Diophantine window edge'],loc='best')
        ylabel('Distance from asymptotic limit',fontsize=20)
        xlabel('Steps down Farey tree',fontsize=20)
        
    
    if pathOriginal == 'RRLL':
        return([np.median(fraction[0:][::2]), np.median(fraction[1:][::2])])
    else:
        return(np.median(fraction))

def plotFractions():
    makeFractionList('RL',plotOn = True, col='c')
    makeFractionList('RLL',plotOn = True,col='k')
    makeFractionList('RLLL',plotOn = True,col='r')
    makeFractionList('RRLLL',plotOn = True,col='m')
    show()
    if False:
        dmin, dmax = 0.0, 0.5
        cmap = cm.get_cmap('gnuplot2')
        junk = contourf([[0,0],[0,0]],np.arange(dmin,dmax+0.01,0.001),cmap=cmap)
        clf()
        figure(2)
        plot([0],[0],'co',[0],[0],'bo',[0],[0],'go',[0],[0],'ko',[0],[0],'mo',[0],[0],'ro')
        ks = np.arange(1.7,5,0.01)
        makeFractionList(1/2.,True, col='b')
        makeFractionList(1/math.pi,True, col='g')

        
        plot([2**((float(x)-1)/2) for x in range(0,50)],'b')
        yscale('log')
        legend(["1",r"$1/2$", r"$1/\pi$", r"$1/e$",'RRLL (silver mean)', 'RL (golden mean)'],loc='best',numpoints=1)
        show()

        fracs1 = []
        fracse = []
        fracsG = []
        fracsSeven = []
        fracsSodd = []
        makeFractionList('RRLL',True, col='m',d=0.2)
        show()

        ds = np.arange(0,0.5,0.01)
        for d in ds:
            asymptote = []
            for k in ks: 
               makeFractionList(True,1,True, d, cmap((d-dmin)/(dmax-dmin)))
               asymptote.append(makeFractionList(1/1.,False,d,k))
            fracs1.append(makeFractionList(1.,False,col='c',d=d))
            fracse.append(makeFractionList(1/math.e,False, col='k',d=d))
            pair = makeFractionList('RRLL',False, col='m',d=d)
            fracsSeven.append(pair[0])
            fracsSodd.append(pair[1])
            fracsG.append(makeFractionList('RL',False, col='r',d=d))

            makeFractionList(1/2.,True, col='c')
            makeFractionList(19512/33201.,True, col='c')
            makeFractionList(1/math.e,True, col='m')
            makeFractionList('RRL',True, col='g')

            makeFractionList(1/math.pi,True)
            fig = figure(1)
            fig.subplots_adjust(right=0.8)
            cbarAx = fig.add_axes([0.85,0.15,0.05,0.7])
            colorbar(junk,cax=cbarAx)

            plot(ks,asymptote,color=cmap(2*d))
            
        plot(ds, fracs1, 'co', ds, fracse, 'ko', ds, fracsG, 'ro', ds, fracsSeven, 'mo', ds, fracsSodd, 'mx',[0,0.5],[0.5,0.5],'k--')
        legend(['1',r"$1/e$",'RL (golden mean)','RRLL (odd)', 'RRLL (even)'],numpoints=1,loc='best')
        ylabel(r"Asymptotic $f_L$ remaining",fontsize=20)
        xlabel("d",fontsize=20)
        show()

        show()

def plotDmin(d = 0.2, k = 2., leg = False):
    alpha = {(math.sqrt(5)-1)/2.:r'$1-1/\varphi$',
              1/math.pi:r'$1/\pi$', 
              1/math.sqrt(2):r'$1/\sqrt{2}$', 
              1/math.sqrt(10):r'$1/\sqrt{10}$',
              1/55.:r'$1/55$'}
    dminLookup = dict()
    colors = ['k','g','r','b','m','c']
    for a in alpha.keys():
        farey, center, path = getFareyPath(a,50,d,k)
        dmin = [1]
        dImmediate = []
        for f in farey:
            dNow = abs(a*f.den**k - f.num*f.den**(k-1))
            dmin.append(min(dNow,dmin[-1]))
            dImmediate.append(dNow)
        dmin.pop(0)
        dminLookup[a] = dmin
        col = colors.pop(0)
        plot(dmin,'-',c=col,label=alpha[a])
        plot(dImmediate,'.',c=col,label="_",markersize=6)
    if leg:
        legend()
    yscale('log')
    [xmin, xmax, ymin, ymax] = axis()
    text(5,math.exp(log(ymax)-(log(ymax)-log(ymin))*0.15),r'$k = $'+str(k),fontsize=25)
 
figure()
subplot(2,2,1) 
plotDmin(k=1.96)
ylabel(r'$d_{min}$')
subplot(2,2,2)
plotDmin(k=2,leg=True)
subplot(2,2,3)
plotDmin(k=2.1)
ylabel(r'$d_{min}$')
xlabel(r'Farey level $\mathscr{F}_j$')
subplot(2,2,4)
plotDmin(k=3)
xlabel(r'Farey level $\mathscr{F}_j$')
currentAxes = axis()
legend()
ax = gca()
ax.legend_.remove()
plot([0,0],'-',[0,0],'.',c='k',markersize=10,linewidth=5)
legend([r'$d_{min}$',r'$d,$ this step'])
axis(currentAxes)
matplotlib.rcParams.update({'font.size': 20})
show()
