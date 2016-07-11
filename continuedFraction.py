from pylab import *
import numpy as np

class CF:
    
    def __init__(self,list):
        self.elements = list
        els = list
        answer = 0
        for i in reversed(range(len(els))):
            answer = 1/(float(els[i])+answer)
        self.val = answer
        self.maxEl = max(list)
        self.weightedMaxEl = max(list)- \
            float(list.index(max(list)))/float(len(list))
        

CFs = []
maxEl = []
maxint = 2
ints = range(1,maxint+1)
for i in ints:
    for j in ints:
        for k in ints:
            for l in ints:
                for m in ints:
                    for n in ints:
                        for o in ints:
                            for p in ints:
                                elements = [i, j, k, l, m, n, o,p]
                                CFs.append(CF(elements))

# plot(np.asarray(CFs)-(math.sqrt(5)-1)/2.,maxEl,'r.',markersize=1)
# plot(CFs[0]-(math.sqrt(5)-1)/2.,maxEl[0],'bD',markersize=3)
# axis([-1,1,0.8,maxint+0.2])
# xscale('symlog',linthreshx = 1e-3)
# show()

histogram = []
step = 0.001
for x in np.arange(0,1,step):
    count = sum([2./2.**i.maxEl for i in CFs if (i.val >= x and i.val < x+step)])
    histogram.append(count)
    
plot(np.arange(0,1,step),histogram)
show()