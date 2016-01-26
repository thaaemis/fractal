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
        
    def checkOverlap(self,other):
        if self.r < other.r:
            return self.high >= other.low
        else:
            return self.low <= other.high
        
    def totalOverlap(self,other):
        return (self.high < other.high) and (other.low < self.low)
        
    def merge(self,other):
        low = min(self.low,other.low)
        high = max(self.high,other.high)
        other.low = low
        other.high = high
        return other
        
    def printFrac(self):
        print self.num,"/",self.den
        

def farey( n, asc=True ):
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

def plotRegions(d,k,nmax,plotYN=True):
    windowSets = []
    for n in range(1,nmax+1):
        a = getWindows(n,d,k,plotYN)
        windowSets.append(a)
    if plotYN:
        xlim(xmin=-0.1,xmax=1.1)
        ylim(ymin=0,ymax=nmax+2)
        # show()

def main(d,k,nmax):
    # Define initial regions.
    regions = [window(0,1,d,k),window(1,1,d,k)]
    
    # Loop over levels of Farey tree from 2 until nmax:
    for n in range(2,nmax+1):
        # Compute new Farey rationals
        newregions = list(farey(n)-farey(n-1))
        regionsCopy = regions
        # Check what to do with new window around this rational.
        for x in newregions:
            addBool = True
            newWindow = window(x[0],x[1],d,k)
            for y in regionsCopy:
                if not newWindow.checkOverlap(y):
                    continue
                elif newWindow.totalOverlap(y):
                    addBool = False
                    break
                else:
                    regions.remove(y)
                    newWindow = newWindow.merge(y)
                    break
            if addBool:
                regions.append(newWindow)
        regions.sort()
        regionsCopy = regions
        regionsCopy2 = regionsCopy
        regions = []
        for x in regionsCopy:
            new = x
            for y in regionsCopy2:
                if not y.checkOverlap(new):
                    continue
                else:
                    if y.den < new.den:
                        new = x.merge(new)
                    else:
                        new = new.merge(x)
            regions.append(new)
        regions = set(regions)
        regions = list(regions)
                
        print(len(regions))
        plotRegions(d,k,nmax,True)
        for a in regions:
            plot(a.r,nmax+1,'ro',a.low,nmax+1,'k<',
                a.high,nmax+1,'k>',ms=15./sqrt(a.den))
        show()


#main(0.1,2,5)
#main(1,2,5)
main(.35,1.5,8)
                
