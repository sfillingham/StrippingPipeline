import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table


def work():

    data = Table.read('orbitalproperties_ELVIS_2Gyrbound.dat', format = 'ascii')
    ID = np.array(data['ID'])
    fperi = np.array(data['first_peri'])
    mperi = np.array(data['min_peri'])
    fperi_tau = np.array(data['first_peri_tau'])
    mperi_tau = np.array(data['min_peri_tau'])
    r50_tau = np.array(data['r50_tau'])
    
    r50cut = r50_tau > 0.0

    pericut = (fperi > 0.5) & (fperi < 1.0)
    pericut2 = (fperi > 0.5)#&(r50_tau < 1.0)
    print len(fperi_tau[pericut])
    print len(fperi_tau[pericut2])
    print len(fperi_tau)
    #print np.mean(fperi_tau[pericut])
    #print np.median(fperi_tau[pericut])

####################################
#Plot work
####################################

    axwidth = 3
    axlength = 10
    fontsize=28

    #############################################################################################
    ## First figure, this is 'direct comparison' plot for the Spekkens 14 paper
    #############################################################################################
    
    plt.rc('axes',linewidth=axwidth)
    plt.figure(figsize=(13,7))

    plt.subplot2grid((2,4),(0,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([-0.1,1.1,-2.1,9.1])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm First\ Pericenter\ Distance\ ({\it R}_{vir})$', fontsize = 28)
    ax.set_ylabel(r'$\rm \tau_{{\it R}50} - \tau_{peri}\ (Gyr)$', fontsize = 28)

    #plt.text(8.3,-2.5, r'$\rm No\ Stripping$', fontsize = 22, color = 'b')
    #plt.text(8.3,-3.0, r'$\rm RPS\ Only$', fontsize = 22, color = 'r')
    #plt.text(8.3,-3.5, r'$\rm RPS\ and\ KH$', fontsize = 22, color = 'g')

    diff = fperi_tau[r50cut] - r50_tau[r50cut]
    diffcut = diff < 0.0
    #print diff[diffcut]

    #p1a = plt.errorbar(np.log10(ABmstar6), np.log10(rps_blue6), xerr = None, yerr=None, fmt='bo', markersize = 10)
    plt.plot(fperi[r50cut], r50_tau[r50cut] - fperi_tau[r50cut], color = 'DarkBlue', marker = 'o', linestyle = 'None', markersize = 16)

    #print np.where(diff == np.min(diff))
    #print r50_tau[r50cut2]
    #print fperi_tau[r50cut2]
        
    xtickloc = [0.0,0.5,1.0]
    xtickstr = ['$'+str(kk)+'$' for kk in xtickloc]
    ytickloc = [-2,0,2,4,6,8]
    ytickstr = ['$'+str(kk)+'$' for kk in ytickloc]
    
    ax.set_xticks(xtickloc)
    ax.set_xticklabels(xtickstr, position = (0,-0.01))
    ax.set_yticks(ytickloc)
    ax.set_yticklabels(ytickstr)
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(18)
        line.set_markeredgewidth(3)
    for tick in ax.xaxis.get_minor_ticks():
        tick.label1.set_fontsize(fontsize/2)
    for tick in ax.yaxis.get_minor_ticks():
        tick.label1.set_fontsize(fontsize/2)
    
    ax.tick_params(which='major',width=axwidth,length=axlength+5)
    ax.tick_params(which='minor',width=axwidth,length=axlength)
####################################################################################################################
    plt.subplot2grid((2,4),(0,2), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([-0.1,5.1,-2.1,9.1])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm \tau_{peri}\ (Gyr)$', fontsize = 28)
    #ax.set_ylabel(r'$\rm \tau \ difference\ (Gyr)$', fontsize = 28)

    #plt.text(8.3,-2.5, r'$\rm No\ Stripping$', fontsize = 22, color = 'b')
    #plt.text(8.3,-3.0, r'$\rm RPS\ Only$', fontsize = 22, color = 'r')
    #plt.text(8.3,-3.5, r'$\rm RPS\ and\ KH$', fontsize = 22, color = 'g')

    #diff = fperi_tau[r50cut] - r50_tau[r50cut]
    #diffcut = diff < 0.0
    #print diff[diffcut]

    #p1a = plt.errorbar(np.log10(ABmstar6), np.log10(rps_blue6), xerr = None, yerr=None, fmt='bo', markersize = 10)
    plt.plot(fperi_tau[r50cut], r50_tau[r50cut] - fperi_tau[r50cut], color = 'DarkBlue', marker = 'o', linestyle = 'None', markersize = 16)

    #print np.where(diff == np.min(diff))
    #print r50_tau[r50cut2]
    #print fperi_tau[r50cut2]
        
    xtickloc = [0,2,4]
    xtickstr = ['$'+str(kk)+'$' for kk in xtickloc]
    ytickloc = [-2,0,2,4,6,8]
    ytickstr = ['$'+str(kk)+'$' for kk in ytickloc]
    
    ax.set_xticks(xtickloc)
    ax.set_xticklabels(xtickstr, position = (0,-0.01))
    ax.set_yticks(ytickloc)
    ax.set_yticklabels(ytickstr)
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(18)
        line.set_markeredgewidth(3)
    for tick in ax.xaxis.get_minor_ticks():
        tick.label1.set_fontsize(fontsize/2)
    for tick in ax.yaxis.get_minor_ticks():
        tick.label1.set_fontsize(fontsize/2)
    
    ax.tick_params(which='major',width=axwidth,length=axlength+5)
    ax.tick_params(which='minor',width=axwidth,length=axlength)

    #plt.savefig('tau_peri_R50_diff.pdf')

    #############################################################################################
    ## Second figure, cumulative histogram work
    #############################################################################################

    plt.rc('axes',linewidth=axwidth)
    plt.figure(figsize=(8,13))

    plt.subplot2grid((4,2),(0,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.13, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([-0.1,1.1,-0.1,250])
    #plt.axis([0,10,0,1.1])
    ax = plt.gca()
    #ax.set_xlabel(r'$\rm Distance\ ({\it R}_{vir})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Number\ of\ Galaxies$', fontsize = 28)

    plt.text(0.3,200,r'$\rm Pericenter\ Distance$', fontsize = 22, color = 'k')
    plt.text(0.35,180,r'$\rm mean\ =\ $'+np.str(np.mean(fperi)), fontsize = 22, color = 'b')
    plt.text(0.35,160,r'$\rm median\ =\ $'+np.str(np.median(fperi)), fontsize = 22, color = 'b')

    n, bins, patches = plt.hist(fperi, 20, color = 'b', linewidth = 4.0, histtype = 'step')
    #plt.plot(fperi, fperi_tau, color = 'DarkBlue', marker = 'o', linestyle = 'None', markersize = 16)
    #print mperi
          
    xtickloc = [0.0,0.5,1.0]
    xtickstr = ['$'+str(kk)+'$' for kk in xtickloc]
    ytickloc = [0,50,100,150,200]
    ytickstr = ['$'+str(kk)+'$' for kk in ytickloc]
    
    ax.set_xticks(xtickloc)
    ax.set_xticklabels(xtickstr, position = (0,-0.01))
    ax.set_yticks(ytickloc)
    ax.set_yticklabels(ytickstr)
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(18)
        line.set_markeredgewidth(3)
    for tick in ax.xaxis.get_minor_ticks():
        tick.label1.set_fontsize(fontsize/2)
    for tick in ax.yaxis.get_minor_ticks():
        tick.label1.set_fontsize(fontsize/2)
    
    ax.tick_params(which='major',width=axwidth,length=axlength+5)
    ax.tick_params(which='minor',width=axwidth,length=axlength)
#############################################################################################
    plt.subplot2grid((4,2),(2,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.13, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([-0.1,1.1,-0.1,5.1])
    #plt.axis([0,10,0,1.1])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Distance\ ({\it R}_{vir})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Infall\ Timescale\ (Gyr)$', fontsize = 28)

    #plt.text(0.4,200,r'$\rm Pericenter\ Distance$', fontsize = 22, color = 'k')
    #plt.text(0.5,180,r'$\rm mean\ =\ $'+np.str(np.mean(fperi)), fontsize = 22, color = 'b')
    #plt.text(0.5,160,r'$\rm median\ =\ $'+np.str(np.median(fperi)), fontsize = 22, color = 'b')

    #n, bins, patches = plt.hist(fperi, 20, color = 'b', linewidth = 4.0, histtype = 'step')
    plt.plot(fperi, fperi_tau, color = 'DarkBlue', marker = 'o', linestyle = 'None', markersize = 16)
    #print mperi
          
    xtickloc = [0.0,0.5,1.0]
    xtickstr = ['$'+str(kk)+'$' for kk in xtickloc]
    ytickloc = [0,1,2,3,4]
    ytickstr = ['$'+str(kk)+'$' for kk in ytickloc]
    
    ax.set_xticks(xtickloc)
    ax.set_xticklabels(xtickstr, position = (0,-0.01))
    ax.set_yticks(ytickloc)
    ax.set_yticklabels(ytickstr)
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(18)
        line.set_markeredgewidth(3)
    for tick in ax.xaxis.get_minor_ticks():
        tick.label1.set_fontsize(fontsize/2)
    for tick in ax.yaxis.get_minor_ticks():
        tick.label1.set_fontsize(fontsize/2)
    
    ax.tick_params(which='major',width=axwidth,length=axlength+5)
    ax.tick_params(which='minor',width=axwidth,length=axlength)

    #plt.savefig('Rdistributions_peri.pdf')

    #############################################################################################
    ## Third figure, histogram work
    #############################################################################################

    plt.rc('axes',linewidth=axwidth)
    plt.figure(figsize=(10,10))

    plt.subplot2grid((2,2),(0,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([-0.1,10.1,0,650])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Timescale\ (Gyr)$', fontsize = 28)
    ax.set_ylabel(r'$\rm Number\ of\ Galaxies$', fontsize = 28)

    plt.text(4,600,r'$\rm Pericenter\ Timescale$', fontsize = 22, color = 'k')
    plt.text(4.25,550,r'$\rm mean\ =\ $'+np.str(np.mean(fperi_tau)), fontsize = 22, color = 'b')
    plt.text(4.25,500,r'$\rm median\ =\ $'+np.str(np.median(fperi_tau)), fontsize = 22, color = 'b')
    plt.text(4,450,r'$\rm 0.5 {\it R}_{vir}\ Timescale$', fontsize = 22, color = 'k')
    plt.text(4.25,400,r'$\rm mean\ =\ $'+np.str(np.mean(r50_tau)), fontsize = 22, color = 'g')
    plt.text(4.25,350,r'$\rm median\ =\ $'+np.str(np.median(r50_tau)), fontsize = 22, color = 'g')

    n, bins, patches = plt.hist(fperi_tau, 15, color = 'b', linewidth = 4.0, histtype = 'step')
    n, bins, patches = plt.hist(r50_tau, 20, color = 'g', linewidth = 4.0, histtype = 'step')
          
    xtickloc = [0,2,4,6,8,10]
    xtickstr = ['$'+str(kk)+'$' for kk in xtickloc]
    ytickloc = [0,200,400,600]
    ytickstr = ['$'+str(kk)+'$' for kk in ytickloc]
    
    ax.set_xticks(xtickloc)
    ax.set_xticklabels(xtickstr, position = (0,-0.01))
    ax.set_yticks(ytickloc)
    ax.set_yticklabels(ytickstr)
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(18)
        line.set_markeredgewidth(3)
    for tick in ax.xaxis.get_minor_ticks():
        tick.label1.set_fontsize(fontsize/2)
    for tick in ax.yaxis.get_minor_ticks():
        tick.label1.set_fontsize(fontsize/2)
    
    ax.tick_params(which='major',width=axwidth,length=axlength+5)
    ax.tick_params(which='minor',width=axwidth,length=axlength)

    #plt.savefig('taudistributions_peri_R50.pdf')

###################################################
# End of Plotting, time to display figure
###################################################

    plt.show()

########################################################################################
