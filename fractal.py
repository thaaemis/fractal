from pylab import *

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
        
    def printFrac(self):
        print self.num,"/",self.den
        

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
    
def plotRegions(d,k,nmax,plotYN=True): # plot regions for each n
    windowSets = []
    
    def getWindows(n,d,k,plotYN=True):
        rationals = farey(n)
        windows = []
        for x in rationals:
            a = window(x[0],x[1],d,k)
            windows.append(a)
            if plotYN:
                plot(a.r,n,'ro',a.low,n,'k<',
                a.high,n,'k>',ms=15./sqrt(x[1]))
        return windows
    
    for n in range(1,nmax+1):
        a = getWindows(n,d,k,plotYN)
        windowSets.append(a)
    if plotYN:
        xlim(xmin=-0.1,xmax=1.1)
        ylim(ymin=0,ymax=nmax+2)

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
    
def main(d,k,nmax):
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
        
        # plotRegions(d,k,nmax,True)
        # for a in cleanRegions:
            # plot(a.r,nmax+1,'ro',a.low,nmax+1,'k<',
                # a.high,nmax+1,'k>',ms=15./sqrt(a.den))
        # show()
    
    xL = [x.low  for x in cleanRegions]
    xR = [x.high for x in cleanRegions]
    x, gradp = stepFxn(xL,xR,.0000001)
    plot(x,gradp)
    show()
    p = integStepFxn(x,gradp)
    plot(x,p)
    show()
    plot(x,integFractal(x,p))
    show()
    return p[-1]

n = []
pressuresBigWindow = []
pressuresQuickDecay = []
pressuresSmallWindow = []

for i in range(150,151,25):
    n.append(i)
    pressuresBigWindow.append(main(.25,2,i))
    pressuresSmallWindow.append(main(.1,2,i))
    pressuresQuickDecay.append(main(.25,3,i))

# plot(n,pressuresBigWindow,'r',label='BigWindow')
# plot(n,pressuresSmallWindow,'g',label='SmallWindow')
# plot(n,pressuresQuickDecay,'b',label='QuickDecay')
# legend()
# show()