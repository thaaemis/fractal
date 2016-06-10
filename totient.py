import fractions
import numpy as np
from pylab import *

def totient(x):
    n = 0
    for k in range(0,x+1):
        if fractions.gcd(k,x) == 1:
            n += 1
    return n
    
x = np.arange(1,5001)
out = [totient(y) for y in x]
plot(x,out,'r.')
title('Totient Function')
xlabel('denominator m')
ylabel('# of new excluded regions')
show()