from pylab import *
import numpy as np
from fractal import getP #(d, k, nmax), returns x, p, gradp
import pickle

# for a given d, want to iterate to find:
#  - Critical k at which pressure goes to 0
#  - last irrational window which survives

# Above is done, but I need to do some iteration with increasing n
# to find if the "last irrational" will survive up to n -> inf

def critK(d,farey,wantPlot=False):
    kGuess = -math.log(1/(2*d)-1)/math.log(2) # if 0, 1/2 overlap
    step = 1.
    while True:
        if max(getP(d,kGuess,farey)[1]) == 0:
            kGuess = kGuess + step
        else:
            if step < 1e-5:
                break
            else:
                kGuess = kGuess - step
                step = step/10.
    if wantPlot:
        xK, pK, gradpK = getP(d,kGuess+step,farey)
        plot(xK,gradpK,label='A little pressure',linewidth=2)
        xM, pM, gradpM = getP(d,kGuess-step,farey)
        plot(xM,gradpM,'r',label='No Pressure')
        ylim((-0.1,1.1))
        legend()
        show()
    x, p, gradp = getP(d, kGuess+step, farey)
    return kGuess, (x[3]+x[2])/2 # last irrational < 0.5

def getDX(farey):
    ds = np.arange(0.38,0.44,0.001)
    ks = []
    xirrs = []
    for d in ds:
        k, xirr = critK(d,farey)
        ks.append(k)
        xirrs.append(xirr)

    # for i in range(len(ds)):
        # print("%.2f, %.5f" % (ds[i],ks[i]))

    fig = figure(figsize=(6,10),dpi=300)    
    
    subplot(211)
    title('d vs. k_Critical, Farey = '+str(farey))
    plot(ks, ds)
    xlabel('k Critical')
    ylim((0.38,0.44))
    xlim((2,3))

    ylabel('d')
    subplot(212)
    plot(ds,xirrs)
    xlim((0.38,0.44))
    #ylim((0.38,0.44))
    xlabel('d')
    ylabel('Last surviving irrational')
    savefig('C:/Users/bkraus/Dropbox/Hudson/criticalK_images_021516/critKvsD_n'+str(farey)+'.png')
    with open('C:/Users/bkraus/Dropbox/Hudson/criticalK_images_021516/critKvsD_n'+str(farey)+'.pkl','w') as f:
        pickle.dump([ds,ks,xirrs], f)
    clf()
 
def replotFromPkl():
    kArray = []
    for farey in range(2,200,1):
        with open('C:/Users/bkraus/Dropbox/Hudson/criticalK_images_021216/critKvsD_n'+str(farey)+'.pkl','r') as f:
            [ds,ks,xirrs] = pickle.load(f)
        kArray.append(ks)
        fig = figure(figsize=(6,10),dpi=200)    

        subplot(211)
        title('d vs. k_Critical, Farey = '+str(farey))
        plot(ks, ds)
        xlabel('k Critical')
        ylabel('d')
        ylim((0,0.5))
        xlim((0,6))
        subplot(212)
        plot(ks,xirrs)
        xlabel('k')
        ylabel('Last surviving irrational')
        ylim((0,0.5))
        xlim((0,6))
        savefig('C:/Users/bkraus/Dropbox/Hudson/criticalK_images_021216/plotsKvsXirr/KvsxIrr_n'+str(farey)+'.png')
        close()
    
    for i in range(0,197,1):
        plot(kArray[i], ds, 'k', alpha = (i+100)*1.0/400.)

    plot(kArray[197], ds, 'r',linewidth=2)
    plot(-np.log(1/(2*np.array(ds))-1)/math.log(2),ds,'g',linewidth=2)

    title('d vs. k_Critical, Farey from 2 to 200')
    xlabel('k Critical')
    ylabel('d')
    ylim((0,0.5))
    xlim((0,6))
    show()
         
replotFromPkl()
         
# for farey in range(2,100,1):
   # getDX(farey)

