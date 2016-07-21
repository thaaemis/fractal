from pylab import *
import numpy as np
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.integrate import quad
    
def integFractal(x,fractal):
    fxn = []
    for i in range(1,len(x)+1):
        fxn.append(trapz(fractal[0:i],x[0:i]))
    return fxn
    
def getP(d, k, nmax, fareyMethod = 'maxDen', getN = False): # returns x, p, gradp

    def integStepFxn(x,fxn):
        integral = [0]
        for i in range(1,len(x),2):
            integral.append(integral[i-1]+(x[i]-x[i-1])*fxn[i])
            integral.append(integral[i-1]+(x[i]-x[i-1])*fxn[i])
        integral.pop()
        return integral    

    def stepFxn(xL,xR,epsilon):
        xL = xL[1:] # trim irrelevant start point
        xR = xR[0:-1] #trim irrelevant end point
        x, fxn = [0], [0]
        for i in range(0,len(xR)):
            # Stay zero
            fxn.append(0)
            x.append(xR[i]) #-epsilon)
            # Jump to 1
            fxn.append(1)
            x.append(xR[i]) #+epsilon)
            #Stay 1
            fxn.append(1)
            x.append(xL[i]) #-epsilon)
            # Jump to 0
            fxn.append(0)
            x.append(xL[i]) #+epsilon)
        x.append(1)
        fxn.append(0)
        return x,fxn

    def farey( n, asc=True ): # give all rationals for nth farey tree
        """Python function to print the nth Farey sequence, either ascending or descending."""
        row = set()
        if asc: 
            a, b, c, d = 0, 1,  1 , n     # (*)
        else:
            a, b, c, d = 1, 1, n-1, n     # (*)
        row.add((a,b))
        while (asc and c <= n) or (not asc and a > 0):
            k = int((n + b)/d)
            a, b, c, d = c, d, k*c - a, k*d - b
            row.add((a,b))
        return row

    class window:
        def __init__(self,n,m,d,k):
            self.r = float(n)/float(m) # rational center
            self.num = n # numerator
            self.den = m # denominator
            self.low  = (n * m**(k-1) - d)/m**k # excluded upper bound
            self.high = (n * m**(k-1) + d)/m**k # excluded low bound
            
        def __lt__(self,other):
            return self.r < other.r
            
        def getLowDen(self,list):
            lowDen = min(x.den for x in list)
            lowWindow = (y for y in list if y.den == lowDen)
            return min(lowWindow)
            
        def overlap(self,other):
            if self.r < other.r:
                return self.high >= other.low
            else:
                return self.low <= other.high
                    
        def merge(self,overlaps):
            overlaps.append(self)
            outWindow = self.getLowDen(overlaps)
            outWindow.low = min(x.low for x in overlaps)
            outWindow.high = max(x.high for x in overlaps)
            return outWindow

    # Define initial regions: clean = no overlaps.
    cleanRegions = [window(0,1,d,k),window(1,1,d,k)]
    
    # Make list of all rationals to consider.
    fareyLevels = []
    if fareyMethod == 'maxDen':
        if getN != False:
            n = nmax
            while len(farey(n)) < getN:
                n = n + 1
            nmax = n
        for n in range(2, nmax+1):
            fareyLevels.append(farey(n)-farey(n-1))

    elif fareyMethod == 'treeSteps':
        fareyAll = [(0,1),(1,1)]
        fareyLevels.append(set(fareyAll))
        for n in range(nmax+1):
            newLevel, count = [], 1
            fareyAllCopy = [x for x in fareyAll]
            # newLevel.append(fareyAll[i])
            for i in range(len(fareyAll)-1):
                newLevel.append((fareyAllCopy[i][0]+fareyAllCopy[i+1][0],fareyAllCopy[i][1]+fareyAllCopy[i+1][1]))
                fareyAll.insert(i+count,newLevel[-1])
                count += 1
            fareyLevels.append(set(newLevel))
    nTot = sum([len(x) for x in fareyLevels]) # How many rationals are considered
                
    # Loop over levels of Farey tree from 2 until nmax:
    for newregions in fareyLevels:
        # Check what to do with new window around this rational.
        for x in newregions:
            # Set up new window with zero overlaps so far
            newWindow = window(x[0],x[1],d,k)
            # Remake old region list to iterate over
            regionsMod = cleanRegions
            cleanRegions, overlapRegions = [], []
            for y in regionsMod:
                if not newWindow.overlap(y):
                    cleanRegions.append(y)
                else:
                    overlapRegions.append(y)
            if len(overlapRegions) == 0: # no overlaps
                cleanRegions.append(newWindow)
            else:
                cleanRegions.append(newWindow.merge(overlapRegions))
        cleanRegions.sort()
        
        # for a in cleanRegions:
            # plot(a.r,n,'ro',a.low,n,'k<',
                # a.high,n,'k>',ms=15./sqrt(a.den))    
    # ylim(1,nmax+1)
    # show()
    
    xL = [x.low  for x in cleanRegions]
    xR = [x.high for x in cleanRegions]
    x, gradp = stepFxn(xL,xR,.0000001)
    p = integStepFxn(x,gradp)
    return x, p, gradp # , nTot

def deriv(y,x):
    return np.gradient(y)/np.gradient(x)
    
def cylinderB(d,k,nmax,RKstep,R=2.,fareyMethod = 'maxDen'): # returns iota, Bz, Bth, Jth, Jz, x, p, gradp
    
    x, p, gradp = getP(d, k, nmax, fareyMethod = fareyMethod)
    p.reverse()
    gradp.reverse()
    
    def iota(r):
        return 1-7*r**2/8.
    def iotaPrime(r):
        return -7*r/4.
    
    def BzPrime(r, Bz, gradp):
        BzP = -(R**2 + iota(r)**2 * r**2)**-1 * (gradp * R**2 / Bz + \
            Bz * (r**2 * iota(r) * iotaPrime(r) + 2 * r * iota(r)**2) )
        return BzP
    
    # set initial conditions
    B0 = 3.
    Bz, r, i = [B0], [0], 0
    BzNo, BzYes = [B0], [B0]
    pFull, pPrimeFull = [p[0]],[0]
    
    
    
    while r[-1] <= 1:
        # do RK4 method for solving for Bz(r)
        h = min([x[i+1]-x[i],RKstep])
        rNow, BzNow = r[-1], Bz[-1]
        BzNoNow, BzYesNow = BzNo[-1], BzYes[-1]
        k1 = BzPrime(rNow,BzNow,gradp[i])
        k2 = BzPrime(rNow+h/2,BzNow+h/2*k1,gradp[i])
        k3 = BzPrime(rNow+h/2,BzNow+h/2*k2,gradp[i])
        k4 = BzPrime(rNow+h,BzNow+h*k3,gradp[i])
        
        # find Bz(r) in instances of all gradp == 0, gradp == 1
        # gradp == 0
        k1no = BzPrime(rNow,BzNoNow,0)
        k2no = BzPrime(rNow+h/2,BzNoNow+h/2*k1,0)
        k3no = BzPrime(rNow+h/2,BzNoNow+h/2*k2,0)
        k4no = BzPrime(rNow+h,BzNoNow+h*k3,0)
        # gradp == 1
        k1yes = BzPrime(rNow,BzYesNow,1)
        k2yes = BzPrime(rNow+h/2,BzYesNow+h/2*k1,1)
        k3yes = BzPrime(rNow+h/2,BzYesNow+h/2*k2,1)
        k4yes = BzPrime(rNow+h,BzYesNow+h*k3,1)
        
        r.append(rNow+h)
        Bz.append(BzNow+h/6*(k1+2*k2+2*k3+k4))
        BzNo.append(BzNoNow+h/6*(k1no+2*k2no+2*k3no+k4no))
        BzYes.append(BzYesNow+h/6*(k1yes+2*k2yes+2*k3yes+k4yes))
        pFull.append(p[i])
        pPrimeFull.append(gradp[i])
        
        # Update counter if we pass an index for x, p, gradp
        if r[-1] >= x[i+1]:
            try:
                i = i+1
                x[i+1]
            except IndexError:
                break
    
    r = np.array(r)
    Bz = np.array(Bz)
    Btheta = r*np.array(Bz)*iota(r)/R
    Jtheta = -deriv(Bz,r)
    with np.errstate(divide='ignore'):
        Jz = 1/r*deriv(r*Btheta,r)
    
    return(r.tolist(), Bz.tolist(), Btheta.tolist(), 
           Jtheta.tolist(), Jz.tolist(), x, p, gradp, pFull, pPrimeFull,
           BzNo, BzYes)

def makePlots():           
    # Make many plots of B, J, p for different parameters.
    params = [('maxDen',200),('treeSteps',11),('treeSteps',12)]
    params.sort(reverse=False)
    fareyLevel = 15
    fig = figure(1) # subplot(gs[0,0:2])
    ax = fig.add_subplot(111)
    # axins = inset_axes(ax, 3,3, loc=3)

    figure(2)
    gs = gridspec.GridSpec(2, 2, height_ratios = [1,1]) # gridspec.GridSpec(3, 2, height_ratios = [1.5,1,1])
    for param in params:
        r, Bz, Btheta, Jtheta, Jz, x, p, gradp = cylinderB(.15,2,param[1],0.00001,fareyMethod=param[0])
        magB = sqrt(np.array(Bz)**2+np.array(Btheta)**2).tolist()
        figure(1)
        ax.plot(x,np.array(p))
        # axins.plot(x,np.array(p))
        figure(2)
        subplot(gs[0,0]) # subplot(gs[1,0])
        plot(r,Bz)
        subplot(gs[1,0])
        plot(r,Btheta)
        subplot(gs[1,1])
        plot(r,Jtheta)
        subplot(gs[0,1])
        plot(r,Jz)
        
    figure(1) # subplot(gs[0,0:2])
    sca(ax)
    text(0.8, 0.5, r"$p(r)$", fontsize=30)
    xlabel(r'r/a', fontsize=20)
    legend(['R = '+str(g) for g in params])
    [xmin, xmax, ymin, ymax] = axis()
    plot([2/(1+math.sqrt(5)),2/(1+math.sqrt(5))],[ymin,ymax],'k--')
    # sca(axins)
    # plot([2/(1+math.sqrt(5)),2/(1+math.sqrt(5))],[ymin,ymax],'k--')
    # axis([0.618,0.61808,0.,0.51])
    # xticks([])
    # yticks([])
    # text(0.617965,.33445,r"55/89",color="blue")
    # text(.618046, 0.3341,r"89/144",color="green")
    # annotate(r"144/233",xy=(.618026,0.333781),xytext=(.61798,.3334),color="red",
    #     arrowprops=dict(facecolor='red'))
    figure(2)
    subplot(gs[0,0])
    text(0.1, 0.85, r"$B_z(r)$", fontsize=26)
    subplot(gs[1,0])
    text(0.9, max(Btheta)*0.8, r"$B_\theta(r)$", fontsize=26)
    xlabel(r'r/a',fontsize=20)
    subplot(gs[1,1])
    text(0.95, max(Jtheta)*0.8, r"$J_\theta(r)$", fontsize=26)
    xlabel(r'r/a',fontsize=20)
    subplot(gs[0,1])
    text(0.1, -0.1, r"$J_z(r)$", fontsize=26)
    legend(['d = '+str(g) for g in params])
    show()
    
def pressureAsymptotes():
    #Plot saturation of p as Farey tree saturates
    n = []
    pressuresBigWindow = []
    pressuresQuickDecay = []
    pressuresSmallWindow = []

    for i in range(10,101,10):
        n.append(i)
        pressuresBigWindow.append(max(getP(.25,2,i)[1]))
        pressuresSmallWindow.append(max(getP(.1,2,i)[1]))
        pressuresQuickDecay.append(max(getP(.25,3,i)[1]))

    plot(n,pressuresBigWindow,'r',label='d = 0.25, k = 2')
    plot(n,pressuresSmallWindow,'g',label='d = 0.1, k = 2')
    plot(n,pressuresQuickDecay,'b',label='d = 0.25, k = 3')
    legend()
    xlabel('Farey tree level')
    ylabel('Peak pressure')
    show()

# Show convergence to some max(p) with Farey ordering, maxDen ordering
def rationalOrdering():
    d = 0.15
    k = 2
    maxLevel = 15
    nFarey, fareyMax, nDen, denMax = [], [], [], []
    for nmax in range(1,maxLevel):
        x, p, gradp, nTotFarey = getP(d, k, nmax, 'treeSteps')
        fareyMax.append(max(p))
        nFarey.append(nTotFarey)
        x, p, gradp, nTotDen = getP(d, k, nmax, 'maxDen', nTotFarey)
        denMax.append(max(p))
        nDen.append(nTotDen)
        print(nTotFarey, nTotDen)
    plot(nFarey,fareyMax,'bx', nDen, denMax, 'rx')
    xscale('log')
    show()
    
def streamsBz():

    R = 2

    def iota(r):
        return 1-7*r**2/8.
        
    def iotaPrime(r):
        return -7*r/4.

    def BzPrime(r, Bz, gradp):
        BzP = -(R**2 + iota(r)**2 * r**2)**-1 * (gradp * R**2 / Bz + \
            Bz * (r**2 * iota(r) * iotaPrime(r) + 2 * r * iota(r)**2) )
        return BzP
    
    def Sprime(r, S):
        Sp = -(r**2 * iota(r) * iotaPrime(r) + 2*r*iota(r)**2)/(R**2 + r**2 * iota(r)**2) * S
        return Sp
        
    def Fprime(r, S, F, gradp):
        Fp = -R**2*gradp/((R**2 + iota(r)**2 * r**2)*S**2)/F
        return Fp
    
    size = 15
    r = np.linspace(0,1,size)
    Bz = np.linspace(0,3,size)
    S = Bz
    F = Bz
    Rad, BZgrid = np.meshgrid(r,Bz)
    
    BzPNo = BzPrime(Rad,BZgrid,0)
    BzPYes = BzPrime(Rad,BZgrid,1)
    dx = np.ones((size,size))
    SP = Sprime(Rad, BZgrid)
    FP = Fprime(Rad, BZgrid, BZgrid, 1)
    
    subplot(3,1,1)
    quiver(Rad,BZgrid,dx,BzPNo, color='k', headwidth=2.3, label=r'$\nabla p = 0$')
    quiver(Rad,BZgrid,dx,BzPYes,color='r', headwidth=2.3, label=r'$\nabla p = 1$')
    plot([0],[0],'k',[0],[0],'r')
    leg = legend([r'$\nabla p = 0$',r'$\nabla p = 1$'],numpoints=1,fontsize=20)
    tmp = leg.get_texts()[1]
    setp(tmp, color='r')
    xlabel(r'$r$',fontsize=20)
    ylabel(r'$B_z$',fontsize=20)
    
    subplot(3,1,2)
    streamplot(Rad,BZgrid,dx,SP, color='g') # headwidth=2.3, label=r'$\nabla p = 0$')
    subplot(3,1,3)
    streamplot(Rad,BZgrid,dx,FP,color='b') #headwidth=2.3, label=r'$\nabla p = 1$')

    
    show()
    
    
# makePlots()

# Show "honing in" on golden mean as d increases towards critical value
# 
# for d in [0.3, 0.38, 0.381, 0.3819, 0.38196, 0.381966, (3-math.sqrt(5))/2]:
    # r, Bz, Btheta, Jtheta, Jz, x, p, gradp = cylinderB(d,2,10,0.00001,fareyMethod='treeSteps')
    # subplot(1,2,1)
    # plot(x, [-1*i for i in gradp])
    # xlabel('x')
    # axis([0,1,-1.1,0.1])
    # subplot(1,2,2)
    # plot(x,np.asarray(p)/max(p))
    # xlabel('x')
# # axis([0,1,-0.01,0.2])
# legend([0.3, 0.38, 0.381, 0.3819, 0.38196, 0.381966],loc='best')
# show()

# Shows how Bz is bound between gradP = 1 case and gradP = 0

r, Bz, Btheta, Jtheta, Jz, x, p, gradp, pFull, pPrimeFull, BzNo, BzYes = cylinderB(0.12,2, 10, 0.00001, 
                fareyMethod='treeSteps')

rOn, BzOn, rOff, BzOff = [], [], [], []
rF, F = [], [1]
for i in range(len(r)):
    if pPrimeFull[i] == 0:
        rOff.append(r[i])
        BzOff.append(Bz[i])
        rF.append(r[i])
        F.append(F[-1])
    else:
        rOn.append(r[i])
        BzOn.append(Bz[i])

def iota(r):
    return 1-7*r**2/8.
        
def iotaPrime(r):
    return -7*r/4.

def smoothFunction(x):
    def integrand(r):
        R = 3.
        return ((r**2 * iota(r) * iotaPrime(r) + 2*r*iota(r)**2)/(R**2 + r**2 * iota(r)**2))
    R = 3.
    math.exp(-1*float(quad(lambda r:(r**2*iota(r)*iotaPrime(r)+2*r*iota(r)**2)/(R**2+r**2*iota(r)**2), 0., x)[0]))

R = 3
print(quad(lambda k:(k**2*iota(k)*iotaPrime(k)+2*k*iota(k)**2)/(R**2+k**2*iota(k)**2), 0., x))
plot(r,smoothFunction(r))
show()
    
F = F[1:]
plot(rOn,np.asarray(BzOn)+0.02,'b.',markersize=1,label=r"$\nabla p = 1$")
plot(rOff,BzOff,'r.',markersize=1,label=r"$\nabla p = 0$")
plot(rF, F, 'g.',markersize=1,label='Fractal part')
legend()
show()

                
# subplot(2,1,1)
# plot(r,Bz,'b',r,Btheta,'r')
# plot(r,BzNo,'k--',r,BzYes,'g--')
# subplot(2,1,2)
# scaledBz, scaledBtheta, scaledJz, scaledJtheta = [], [], [], []
# for i in range(len(Bz)):
    # if pFull[i] == 0:
        # scaledBz.append(0)
        # scaledBtheta.append(0)
    # else:
        # scaledBz.append(Bz[i]/(pFull[i]))
        # scaledBtheta.append(Btheta[i]/pFull[i])
# plot(r,scaledBz,'b',r,scaledBtheta,'r')
# plot(r,BzNo,'k--',r,BzYes,'g--')
# show()