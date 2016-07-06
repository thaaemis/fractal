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
        if other.__class__.__name__ == 'rational':
            return self.val - other.val
        else:
            return self.val - other
        
    def __repr__(self):
        return str(self.num)+'/'+str(self.den)

def getFareyPath(path,maxN,d,k, convergents = False):    
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
        if target != None:
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
        nextUp = nextLevel(farey[-1],winner)
        if abs(nextUp - farey[-1]) < 1e-14 or (target != None and nextUp.val == target): 
            break        # So there aren't problems with computer accuracy
        farey.append(nextUp)
        
    if type(target) == (float or int):
        center = target
    else:
        center = np.mean([x.val for x in farey[-5:-1]])
    if convergents == False:
        return(farey, center, path)
    else:
        convergents = [farey[0]] if center < 0.5 else [farey[1]]
        indices = [0]
        for i in range(1,len(farey)-1):
            # last mediant on one side of target is a convergent
            if ( (farey[i+1] - center) > 0 ) != ( (farey[i] - center) > 0 ):
                convergents.append(farey[i])
                # indices.append(i)
        return(convergents,center,path) # ,indices)

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

def getDmax(d = 0.2, k = 2., leg = False, plotOn = False):
    class number:
        def __init__(self, val, string, type):
            self.val = val
            self.string = string
            self.type = type
            
    alpha = [number((math.sqrt(5)-1)/2., r'$1-1/\varphi$','noble-like'),
             number(1/math.pi, r'$1/\pi$', 'irrational'),
             number(1/(math.e+math.pi), r'$1/(\pi+e)$', 'irrational'),
             number(7331/10321., r'$7331/10321$', 'rational') ,
             number(15831/19321., r'$15831/19321$', 'rational') ,
             number(1/math.sqrt(5), r'$1/\sqrt{5}$', 'bound4'),
             number(1/math.sqrt(6), r'$1/\sqrt{6}$', 'bound4'),
             number(1/math.sqrt(8), r'$1/\sqrt{8}$', 'bound4'),
             number(1/math.sqrt(10), r'$1/\sqrt{10}$', 'bound6'),
             number(1/math.sqrt(11), r'$1/\sqrt{11}$', 'bound6'),
             number(1/math.sqrt(15), r'$1/\sqrt{15}$', 'bound6'),
             number(1/math.sqrt(17), r'$1/\sqrt{17}$', 'irrational'),
             number(1/math.sqrt(51), r'$1/\sqrt{51}$', 'irrational'),
             number(1/math.sqrt(101), r'$1/\sqrt{101}$', 'irrational')]
    markerLookup = {'noble-like':'D', 'irrational':'*', 'rational':'^', 'bound4':"8", 'bound6':"."}
    # sizeLookup = {'noble-like':7, 'irrational':8, 'rational':8}
    dMaxLookup = dict()
    colors = ['k','g','r','b','m','c','y']
    for a in alpha:
        farey, center, path = getFareyPath(a.val,600,d,k,convergents=True)
        dmax = [1]
        dImmediate = []
        for f in farey:
            dNow = abs(a.val*f.den**k - f.num*f.den**(k-1))
            dmax.append(min(dNow,dmax[-1]))
            dImmediate.append(dNow)
        dmax.pop(0)
        dMaxLookup[a.val] = dImmediate
        try:
            col = colors.pop(0)
        except IndexError:
            colors = ['k','g','r','b','m','c','y']
            col = '0.75'

        if plotOn:
            # plot(dmax,'-',c=col,label="_")
            plot(dImmediate,'.',c=col,label=a.string,
                markersize=6,marker=markerLookup[a.type])
    if plotOn and leg:
        legend(loc='best',numpoints=1)
    if plotOn:
        yscale('log')
        [xmin, xmax, ymin, ymax] = axis()
        # axis([xmin, 12, ymin, ymax])
        text(8,math.exp(log(ymax)-(log(ymax)-log(ymin))*0.85),r'$k = $'+str(k),fontsize=20)
    return(alpha, dMaxLookup)

def plotDmax():
     
    figure()
    subplot(2,2,1) 
    getDmax(k=1.9,plotOn=True)
    ylabel(r'$d_{max}$',fontsize=30)
    subplot(2,2,2)
    getDmax(k=2,leg=True,plotOn=True)
    [xmin, xmax, ymin, ymax] = axis()
    text(1,math.exp(log(ymax)-(log(ymax)-log(ymin))*0.85),r'$k = 2$',fontsize=20)
    subplot(2,2,3)
    getDmax(k=2.1,plotOn=True)
    ylabel(r'$d_{max}$',fontsize=30)
    xlabel(r'Convergent #',fontsize=30)
    subplot(2,2,4)
    getDmax(k=3,plotOn=True)
    xlabel(r'Convergent #',fontsize=30)
    currentAxes = axis()
    legend()
    ax = gca()
    ax.legend_.remove()
    plot([0,0],'-',[0,0],'.',c='k',markersize=10,linewidth=5)
    # legend([r'$d_{max}$'])# ,r'$d,$ this step'])
    axis(currentAxes)
    matplotlib.rcParams.update({'font.size': 16})
    show()

def testConvergents(): # needs uncommenting indices in getFareyPath above
    d, k, a = 0.2, 1.9, 1/math.sqrt(2)
    convergents, center, path, indices = getFareyPath(a,300,d,k,convergents=True)
    farey, center, path = getFareyPath(a,100,d,k,convergents=False)
    dFarey, dConvergents = [], []
    for f in farey:
        dNow = abs(a*f.den**k - f.num*f.den**(k-1))
        dFarey.append(dNow)
    for f in convergents:
        dNow = abs(a*f.den**k - f.num*f.den**(k-1))
        dConvergents.append(dNow)
    plot(dFarey,'bx',indices,dConvergents,'ko')
    yscale('log')
    show()
    
def nobleSlopes():

    def findAvgSlope(points,step): # gets exponential slope between points of dmax
        slopes = [(log(points[i])-log(points[i-1]))/step for i in range(1,len(points))]
        return np.mean(slopes), np.std(slopes)

    markerLookup = {'noble-like':'D', 'irrational':'*', 'rational':'^', 'bound4':"8", 'bound6':"."}
    figure()
    cmap = cm.spectral
    allSlopes = dict()
    kstep = 0.05
    for k in np.arange(1,3,kstep):
        nums, dMaxes = getDmax(d=0.2,k=k)
        slopeLookup = dict()
        maximum = 0
        for num in nums:
            maximum = max(abs(num.val - 2./(1+sqrt(5))), maximum)
        for num in nums:
            dmax = dMaxes[num.val]
            dmax = dmax[4:]
            all = findAvgSlope(dmax,1)
            even = findAvgSlope(dmax[::2],2)
            odd = findAvgSlope(dmax[1::2],2)
            if all[1] > abs(0.1*all[0]):
                slopeLookup[num.val] = np.mean([even[0],odd[0]])
            else:
                slopeLookup[num.val] = all[0]
            label = '_' if k != 2.5 else num.string
            plot(k,slopeLookup[num.val],'o',c=cmap(abs(num.val- 2./(1+sqrt(5)))/maximum),
                    marker=markerLookup[num.type],label=label,markersize=10)
        allSlopes[k] = slopeLookup
    for num in nums:
        slopevsK = []
        ks = sort(allSlopes.keys())
        for k in ks:
            x = allSlopes[k][num.val]
            if x != NAN:
                slopevsK.append(x)
        slope = mean([(slopevsK[i] - slopevsK[i-1])/kstep for i in range(1,len(slopevsK))])
        print(num.string, num.val, '%.2f' %slope)
            
    legend(loc='best',numpoints=1)
    
    show()
        
        
    
    
# plotDmax()

nobleSlopes()

def examineDenominators(target = 1./sqrt(7)):

    farey, center, path = getFareyPath(target,500,0.2,2,convergents=True)    
    print([x for x in farey])
    subplot(2,1,1)
    plot([abs(x.val-target) for x in farey],'bo')
    yscale('log')
    subplot(2,1,2)
    denominatorRatio = []
    for i in range(1,len(farey)):
        denominatorRatio.append(float(farey[i].den)/float(farey[i-1].den))
    plot(denominatorRatio,'ro')
    show()    
        
    
    
    
    
    
    
    
    
    