from pylab import *
import numpy as np

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
    
def stepFxn(xL,xR,epsilon):
    xL = xL[1:] # trim irrelevant start point
    xR = xR[0:-1] #trim irrelevant end point
    x, fxn = [0], [0]
    for i in range(0,len(xR)):
        # Stay zero
        fxn.append(0)
        x.append(xR[i]-epsilon)
        # Jump to 1
        fxn.append(1)
        x.append(xR[i]+epsilon)
        #Stay 1
        fxn.append(1)
        x.append(xL[i]-epsilon)
        # Jump to 0
        fxn.append(0)
        x.append(xL[i]+epsilon)
    x.append(1)
    fxn.append(0)
    return x,fxn

def integStepFxn(x,fxn):
    integral = [0]
    for i in range(1,len(x),2):
        integral.append(integral[i-1]+(x[i]-x[i-1])*fxn[i])
        integral.append(integral[i-1]+(x[i]-x[i-1])*fxn[i])
    integral.pop()
    return integral    
    
def integFractal(x,fractal):
    fxn = []
    for i in range(1,len(x)+1):
        fxn.append(trapz(fractal[0:i],x[0:i]))
    return fxn
    
def deriv(y,x):
    return np.gradient(y)/np.gradient(x)
    
def getP(d,k,nmax):
    # Define initial regions: clean = no overlaps.
    cleanRegions = [window(0,1,d,k),window(1,1,d,k)]
    
    # Loop over levels of Farey tree from 2 until nmax:
    for n in range(2,nmax+1):
        # Compute new Farey rationals
        newregions = list(farey(n)-farey(n-1))
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
    # return x, p, gradp
    return max(p)

def cylinderB(d,k,nmax,RKstep):
    x, p, gradp = getP(d, k, nmax)
    
    p.reverse()
    gradp.reverse()
    
    # for iota = r = Bz/Bth, make function for Bz' = ...
    def f(iota, Bz, gradp):
        return (-gradp/Bz - 2*iota*Bz)/(1+iota**2)
    
    # set initial conditions
    Bz, iota, i = [1], [0], 0
    
    while iota[-1] <= 1:
        # do RK4 method for solving for Bz(iota)
        h = min([x[i+1]-x[i],RKstep])
        iotaNow, BzNow = iota[-1], Bz[-1]
        k1 = f(iotaNow,BzNow,gradp[i])
        k2 = f(iotaNow+h/2,BzNow+h/2*k1,gradp[i])
        k3 = f(iotaNow+h/2,BzNow+h/2*k2,gradp[i])
        k4 = f(iotaNow+h,BzNow+h*k3,gradp[i])
        iota.append(iotaNow+h)
        Bz.append(BzNow+h/6*(k1+2*k2+2*k3+k4))
        # Update counter if we pass an index for x, p, gradp
        if iota[-1] >= x[i+1]:
            try:
                i = i+1
                x[i+1]
            except IndexError:
                break
    
    iota = np.array(iota)
    Bz = np.array(Bz)
    Btheta = np.array(iota)*np.array(Bz)
    Jtheta = -deriv(Bz,iota)
    Jz = 1/iota*deriv(iota*Btheta,iota)
    
    return(iota.tolist(), Bz.tolist(), Btheta.tolist(), 
           Jtheta.tolist(), Jz.tolist(), x, p, gradp)
    

# params = [0.1,0.2,0.25]
# params.sort(reverse=False)
# fareyLevel = 200
# for param in params:
    # iota, Bz, Btheta, Jtheta, Jz, x, p, gradp = cylinderB(param,2,fareyLevel,0.00001)
    # magB = sqrt(np.array(Bz)**2+np.array(Btheta)**2).tolist()
    # figure(1)
    # subplot(211)
    # plot(iota,Bz)
    # subplot(212)
    # plot(iota,Btheta)
    # # subplot(313)
    # # plot(iota,magB)
    # figure(2)
    # subplot(211)
    # plot(iota,Jtheta)
    # subplot(212)
    # plot(iota,Jz)
    # figure(3)
    # plot(x,p)
# figure(1)
# subplot(211)
# legend(['d = '+str(g) for g in params])
# title('Bz(r)')
# subplot(212)
# title('Btheta(r)')
# figure(2)
# subplot(211)
# legend(['d = '+str(g) for g in params])
# title('Jtheta(r)')
# subplot(212)
# title('Jz(r)')
# figure(3)
# title('p(r)')
# legend(['d = '+str(g) for g in params])
# show()


#Plot saturation of p as Farey tree saturates
n = []
pressuresBigWindow = []
pressuresQuickDecay = []
pressuresSmallWindow = []

for i in range(10,201,10):
    n.append(i)
    pressuresBigWindow.append(getP(.25,2,i))
    pressuresSmallWindow.append(getP(.1,2,i))
    pressuresQuickDecay.append(getP(.25,3,i))

plot(n,pressuresBigWindow,'r',label='d = 0.25, k = 2')
plot(n,pressuresSmallWindow,'g',label='d = 0.1, k = 2')
plot(n,pressuresQuickDecay,'b',label='d = 0.25, k = 3')
legend()
xlabel('Farey tree level')
ylabel('Peak pressure')
show()