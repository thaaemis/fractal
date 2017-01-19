from pylab import *
import numpy as np
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.integrate import quad
from scipy.interpolate import interp1d
from fractal import cylinderB

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

def makePlots():           
    # Make many plots of B, J, p for different parameters.
    params = [(0.1, 2.2)] #, (0.05, 2.2), (0.1, 2.2)] #[('maxDen',200),('treeSteps',11),('treeSteps',12)]
    # params.sort(reverse=False)
    fareyLevel = 15
    fig = figure(1) # subplot(gs[0,0:2])
    axP = gca()
    figure(4)
    axPprime = gca()
    axPprime.axis([0,1,-0.15,1.15])

    
    def plotPprime(x,gradp,ax):
        sca(ax)
        for i in np.arange(0,len(gradp),2):
            plot(x[i:i+2],gradp[i:i+2],'g-',linewidth=2)
        print('ok')

    for param in params:
        r, Bz, Btheta, Jtheta, Jz, x, p, gradp = cylinderB(param[0],param[1],9,0.00001,R=4.,fareyMethod='treeSteps')
        magB = sqrt(np.array(Bz)**2+np.array(Btheta)**2).tolist()
        figure(1)
        p = np.asarray(p)/max(p)

        axP.plot(x,np.array(p),'g',linewidth=4)
        figure(4)
        plotPprime(x,gradp,axPprime)
        #axPprime.plot(x,np.array(gradp),'g.',markersize=4)
        figure(2)
        plot(r,Bz,'.',markersize=4)
#         plot(r,Btheta,'.',markersize=4)
        figure(3)
        axJ = gca()
        axJ.plot(r,Jtheta,'r.',markersize=4)
        axJ.plot(r,Jz,'b.',markersize=4)
        
    figure(1) # subplot(gs[0,0:2])
    text(0.1, 0.35, r"$p(r)$", fontsize=50, color='g')
    
    figure(4)
    text(0.1, 0.35, r"$\nabla p(r)$", fontsize=50, color='g')

#     xlabel(r'r/a', fontsize=20)
    # legend(['R = '+str(g) for g in params])
#     [xmin, xmax, ymin, ymax] = axis()
#     plot([2/(1+math.sqrt(5)),2/(1+math.sqrt(5))],[ymin,ymax],'k--')

    figure(2)
    text(0.1, 1.1, r"$B_z(r)$", fontsize=50,color='b')
    # legend([r'$(d, k)$ = '+str(g) for g in params],loc='best')
#     text(0.9, max(Btheta)*0.8, r"$B_\theta(r)$", fontsize=50)
    xlabel(r'r/a',fontsize=40)
    xticks(fontsize=40)
    yticks(np.arange(0, 1.3, 0.4), fontsize=40)
    axis([0,1.2, -0.1, 1.25])
    axB = gca()
    
    figure(3)
    text(0.35, max(Jtheta)/1.3, r"$J_\theta(r)$", fontsize=50,color='r')
    xlabel(r'r/a',fontsize=40)
    text(0.05, 0.35, r"$J_z(r)$", fontsize=50,color='b')
    xticks(fontsize=40)
    yticks(np.arange(-0.4,1.4,0.4),fontsize=40)
    # legend([r'$(d, k)$ = '+str(g) for g in params])

    if True:
        axins = zoomed_inset_axes(axJ, 5, loc=5,) # zoom = 6
        axins.plot(r,Jz,'b.',r,Jtheta,'r.',markersize=4)
        
        axinsP = inset_axes(axP, 2,2 , loc='best',
            bbox_to_anchor=(0.2, 0.2),bbox_transform=gcf().transFigure)
                                                                       
        #axinsP = zoomed_inset_axes(axP, 15, loc=3) # zoom = 6
        axinsP.plot(x,np.array(p),linewidth=4)

        axinsPprime = zoomed_inset_axes(axPprime, 40, loc=0) # zoom = 6
        # axinsPprime.plot(x,np.array(gradp),'g.',markersize=4)
        plotPprime(x,gradp,axinsPprime)
                
        axins2 = zoomed_inset_axes(axB, 20, loc=5) # zoom = 6
        axins2.plot(r,Bz,'b.',markersize=4)

        figure(3)
        x1, x2, y1, y2 = 0.35, .45, .32, 0.465
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        xticks(visible=False)
        yticks(visible=False)

        figure(1)
        sca(axinsP)
        x1, x2, y1, y2 = 1.02134, 1.0407, 0.09, 0.1107
        axinsP.set_xlim(x1, x2)
        axinsP.set_ylim(y1, y2)
        xticks(visible=False)
        yticks(visible=False)

        figure(4)
        x1, x2, y1, y2 = 0.83, 0.84, 0.995, 1.005
        axinsPprime.set_xlim(x1, x2)
        axinsPprime.set_ylim(y1, y2)
        xticks(visible=False)
        yticks(visible=False)


        figure(2)
        x1, x2, y1, y2 = 0.698, .721, .944, 0.962
        axins2.set_xlim(x1, x2)
        axins2.set_ylim(y1, y2)
        xticks(visible=False)
        yticks(visible=False)

    # draw a bbox of the region of the inset axes in the parent axes and
    # connecting lines between the bbox and the inset axes area
        mark_inset(axJ, axins, loc1=2, loc2=3, fc="none", ec="0.5", lw=1)
        mark_inset(axB, axins2, loc1=1, loc2=2, fc="none", ec="0.5", lw=1)
        mark_inset(axP, axinsP, loc1=1, loc2=4, fc="none", ec="0.5", lw=1)
        mark_inset(axPprime, axinsPprime, loc1=2, loc2=3, fc="none", ec="0.5", lw=1)
        draw()


    show()
    
makePlots()

