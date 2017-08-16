#plotting routine for RPS paper based on get_h1data_fitparams script
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table


def work(hot,vel,tq,check):

    if check == 'pearson':
        ABdata = Table.read('RPS_strippedfrac_hot'+np.str(hot)+'_v'+vel+'_pearsoniterate100.dat', format = 'ascii')
        ABdata2 = Table.read('/Users/Sean/research/elvis/RamPressure/h1profiles/modeldata/RPS_strippedfrac_hot'+np.str(hot)+'_v'+vel+'_Model.dat', format = 'ascii')
    else:
        ABdata = Table.read('stripdata/RPS_strippedfrac_hot'+np.str(hot)+'_v'+vel+'_ABiterate10.dat', format = 'ascii')
        ABdata2 = Table.read('modeldata/RPS_strippedfrac_hot'+np.str(hot)+'_v'+vel+'_Model.dat', format = 'ascii')
        
    ABmstar = np.array(ABdata['Mstar'])
    KH = np.array(ABdata['KH_1Gyr'])
    print len(ABmstar)
    ABfrac = np.array(ABdata['NFW'])
    
    ABmstar2 = np.array(ABdata2['Mstar'])
    KH2 = np.array(ABdata2['KH_1Gyr'])
    print len(ABmstar2)
    ABfrac2 = np.array(ABdata2['NFW'])
    #upperhot = hot-5
    #ABupperdata = Table.read('stripdata/RPS_strippedfrac_hot'+np.str(upperhot)+'_v'+vel+'_AB.dat', format = 'ascii')
    #ABupper = np.array(ABupperdata['NFW'])
    #lowerhot = hot+5
    #ABlowerdata = Table.read('stripdata/RPS_strippedfrac_hot'+np.str(lowerhot)+'_v'+vel+'_AB.dat', format = 'ascii')
    #ABlower = np.array(ABlowerdata['NFW'])
    
    #ABerror = np.empty([2,len(ABmstar)])
    #ABerror[1] = ABupper - ABfrac
    #ABerror[0] = ABfrac - ABlower

    if check == 'kh_hist':
        frac1 = tq*KH
        frac2 = tq*KH2

    elif check == 'iterate':
        frac1 = ABfrac + tq*KH
        frac2 = ABfrac2 + tq*KH2

    elif check == 'model':
        frac1 = ABfrac
        frac2 = ABfrac2

    elif check == 'pearson':
        frac1 = ABfrac
        frac2 = ABfrac2

    elif check == 'scatter':
        frac1 = ABfrac
        frac2 = ABfrac2

    elif check == 'rps_hist':
        frac1 = ABfrac
        frac2 = ABfrac2

    else:
        frac1 = ABfrac + tq*KH
        frac2 = ABfrac2 + tq*KH2
        
    frac1[np.where(frac1 > 1.0)] = 1.00
    frac2[np.where(frac2 > 1.0)] = 1.00

    #bin data and determine the mean and std in each bin
    binsize = 0.3
    iterations = 13
    #binlist = np.array([6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0])
    binlist = np.linspace(6.0, 12.0, 50)
    plotlist = np.array([6.25,6.75,7.25,7.75,8.5,9.5,10.25,10.75,11.25,11.75])
    
    mean1 = np.empty(len(binlist))
    upper1 = np.empty(len(binlist))
    lower1 = np.empty(len(binlist))
    mean2 = np.empty(len(binlist))
    upper2 = np.empty(len(binlist))
    lower2 = np.empty(len(binlist))

    for i in range(len(binlist)):
        cut1 = (np.log10(ABmstar) > (binlist[i]-binsize)) & (np.log10(ABmstar) < (binlist[i]+binsize))
        cut2 = (np.log10(ABmstar2) > (binlist[i]-binsize)) & (np.log10(ABmstar2) < (binlist[i]+binsize))

        mean1[i] = np.mean(frac1[cut1])
        upper1[i] = np.mean(frac1[cut1]) + np.std(frac1[cut1])
        lower1[i] = np.mean(frac1[cut1]) - np.std(frac1[cut1])

        mean2[i] = np.mean(frac2[cut2])
        upper2[i] = np.mean(frac2[cut2]) + np.std(frac2[cut2])
        lower2[i] = np.mean(frac2[cut2]) - np.std(frac2[cut2])

    mean1[np.where(mean1 > 1.00)] = 1.00
    upper1[np.where(upper1 > 1.00)] = 1.00
    lower1[np.where(lower1 > 1.00)] = 1.00

    mean2[np.where(mean2 > 1.00)] = 1.00
    upper2[np.where(upper2 > 1.00)] = 1.00
    lower2[np.where(lower2 > 1.00)] = 1.00
    

############
#Median and mean stripped fraction for object below 10^8 Msun
############

    cut = ABmstar < 10**8.0
    lowfrac = ABfrac[cut]
    print np.mean(lowfrac)
    print np.median(lowfrac)
    print np.std(lowfrac)


####################################
#Plot work
####################################

    axwidth = 3
    axlength = 10
    fontsize=28
    
    plt.rc('axes',linewidth=axwidth)
    
    plt.figure(figsize=(11,8))

############
#FIGURE 
############
    plt.subplot2grid((1,1),(0,0), colspan = 1, rowspan = 1)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.87, wspace=0.6, hspace=0.05)
    plt.axis([5.8,11.2,-0.05,1.05])
    ax = plt.gca()
    ax2 = ax.twiny()
    ax2.set_xlabel(r'$\rm Halo\ Mass\ (M_{\odot})$', fontsize = 28, labelpad = 20)
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    #plt.text(9.3,0.95,r'$\rm RPS:$', fontsize = 24)
    plt.text(9.4,0.9,r'$\rm {\it n}_{halo}\ = \ 10^{-'+np.str(0.1*hot)+'} cm^{-3}$', fontsize = 22)
    plt.text(9.4,0.8,r'$\rm {\it V}_{sat}\ = \ '+vel+'\ km\ s^{-1}$', fontsize = 22)

    if check == 'scatter':
        plt.plot(binlist, mean1, color = 'Black', marker = 'None', linestyle = '-', linewidth = 5.0, alpha = 0.3)
        p1a = plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=None, fmt='yo', markersize = 10, alpha = 0.0)
        #plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=ABerror, fmt='yo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        plt.plot(np.log10(ABmstar), ABfrac, color = 'Goldenrod', marker ='o', markersize = 16, linestyle = 'None')
        #plt.fill_between(binlist, upper1, lower1, facecolor = 'k', alpha = 0.2)
        #plt.plot((10.5-binsize/2.,10.5+binsize/2.),(0.6,0.6), color = 'k', linestyle = '-', linewidth = 2.0)

    elif check == 'model':
        plt.plot(binlist, mean1, color = 'Black', marker = 'None', linestyle = '-', linewidth = 5.0, alpha = 0.3)
        #p1a = plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=None, fmt='yo', markersize = 10, alpha = 0.0)
        #plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=ABerror, fmt='yo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        plt.plot(np.log10(ABmstar2), ABfrac2, color = 'Goldenrod', marker ='o', markersize = 16, linestyle = 'None')
        #plt.fill_between(binlist, upper1, lower1, facecolor = 'k', alpha = 0.2)
        #plt.plot((10.5-binsize/2.,10.5+binsize/2.),(0.6,0.6), color = 'k', linestyle = '-', linewidth = 2.0)

    elif check == 'pearson':
        plt.plot(binlist, mean1, color = 'Black', marker = 'None', linestyle = '-', linewidth = 5.0, alpha = 0.3)
        plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=None, fmt='yo', markersize = 10, alpha = 1.0)
        #plt.plot(np.log10(ABmstar), ABfrac, marker = 'o', color = 'y', markersize = 10, alpha = 1.0, linestyle = 'None')
        #plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=ABerror, fmt='yo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        #plt.plot(np.log10(ABmstar2), ABfrac2, color = 'Goldenrod', marker ='o', markersize = 16, linestyle = 'None')
        #plt.fill_between(binlist, upper1, lower1, facecolor = 'k', alpha = 0.2)
        #plt.plot((10.5-binsize/2.,10.5+binsize/2.),(0.6,0.6), color = 'k', linestyle = '-', linewidth = 2.0)

    elif check == 'kh_hist':
        #plotting routine for KH only histogram
        plt.fill_between(binlist, upper1, lower1, facecolor = 'y', alpha = 0.4)
        plt.plot(binlist, mean1, color = 'Goldenrod', marker = 'None', linestyle = '-', linewidth = 4.0)
        plt.plot((10.5-binsize/2.,10.5+binsize/2.),(0.4,0.4), color = 'k', linestyle = '-', linewidth = 2.0)

    elif check == 'iterate':
        #plotting routine for KH only histogram
        plt.fill_between(binlist, upper1, lower1, facecolor = 'y', alpha = 0.4)
        plt.plot(binlist, mean1, color = 'Goldenrod', marker = 'None', linestyle = '-', linewidth = 4.0)
        plt.plot((10.5-binsize/2.,10.5+binsize/2.),(0.4,0.4), color = 'k', linestyle = '-', linewidth = 2.0)

    else:
        #plotting routine for alt
        plt.fill_between(binlist, upper1, lower1, facecolor = 'y', alpha = 0.4)
        plt.plot(binlist, mean1, color = 'Goldenrod', marker = 'None', linestyle = '-', linewidth = 4.0)
        plt.plot((10.5-binsize/2.,10.5+binsize/2.),(0.4,0.4), color = 'k', linestyle = '-', linewidth = 2.0)
    
    xtickloc = [6,7,8,9,10,11]
    xtickstr = [r'$\rm 10^{6}$',r'$\rm 10^{7}$',r'$\rm 10^{8}$',r'$\rm 10^{9}$',r'$\rm 10^{10}$',r'$\rm 10^{11}$']
    ytickloc = [0.0,0.2,0.4,0.6,0.8,1.0]
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


    ax2.set_xlim(5.8,11.2)
    ax2.set_ylim(-0.05,1.05)
    x2tickloc = [6.487,8.442,10.429]
    x2tickstr = [r'$\rm 10^{10}$',r'$\rm 10^{11}$',r'$\rm 10^{12}$']
    ax2.set_xticks(x2tickloc)
    ax2.set_xticklabels(x2tickstr)#, position = (0,1.02))
    for tick in ax2.xaxis.get_major_ticks():
        tick.label2.set_fontsize(fontsize)
    for tick in ax2.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for line in ax2.get_xticklines() + ax2.get_yticklines():
        line.set_markersize(18)
        line.set_markeredgewidth(3)
    for tick in ax2.xaxis.get_minor_ticks():
        tick.label1.set_fontsize(fontsize/2)
    for tick in ax2.yaxis.get_minor_ticks():
        tick.label1.set_fontsize(fontsize/2)

    ax2.tick_params(which='major',width=axwidth,length=axlength+5)
    ax2.tick_params(which='minor',width=axwidth,length=axlength)

    if check == 'scatter':
        plt.savefig('RPS_motivation_hot'+np.str(hot)+'_v'+vel+'_AB_meanline.pdf')
        plt.savefig('RPS_motivation_hot'+np.str(hot)+'_v'+vel+'_AB_meanline.png')

    elif check == 'model':
        plt.savefig('RPS_motivation_hot'+np.str(hot)+'_v'+vel+'_model_meanline.pdf')

    elif check == 'pearson':
        plt.savefig('RPS_motivation_hot'+np.str(hot)+'_v'+vel+'_pearson_meanline.pdf')

    #elif check == 'kh_hist':
        #plt.savefig('RPS_motivation_hot'+np.str(hot)+'_v'+vel+'_AB_KHhist.pdf')

    elif check == 'iterate':
        plt.savefig('RPS_motivation_hot'+np.str(hot)+'_v'+vel+'_AB_iterate10.pdf')

    else:
        plt.savefig('RPS_motivation_hot'+np.str(hot)+'_v'+vel+'_AB_hist.pdf')

    plt.show()

    
####################################################################################################################################
## Motivation plot with two panels to show both the RPS and KH
####################################################################################################################################

def double(hot,vel,tq,check):
    
    ABdata = Table.read('stripdata/RPS_strippedfrac_hot'+np.str(hot)+'_v'+vel+'_AB.dat', format = 'ascii')
    ABmstar = np.array(ABdata['Mstar'])
    KH = np.array(ABdata['KH_1Gyr'])
    print len(ABmstar)
    ABfrac = np.array(ABdata['NFW'])
    upperhot = hot-5
    ABupperdata = Table.read('stripdata/RPS_strippedfrac_hot'+np.str(upperhot)+'_v'+vel+'_AB.dat', format = 'ascii')
    ABupper = np.array(ABupperdata['NFW'])
    lowerhot = hot+5
    ABlowerdata = Table.read('stripdata/RPS_strippedfrac_hot'+np.str(lowerhot)+'_v'+vel+'_AB.dat', format = 'ascii')
    ABlower = np.array(ABlowerdata['NFW'])
    
    ABerror = np.empty([2,len(ABmstar)])
    ABerror[1] = ABupper - ABfrac
    ABerror[0] = ABfrac - ABlower

    fracRPS = ABfrac
    fracKH = tq*KH

    #if check == 'kh_hist':
        #fracRPS = ABfrac
        #fracKH = tq*KH

    #else:
        #fracKH = ABfrac + tq*KH
        
    fracKH[np.where(fracKH > 1.0)] = 1.00

    #bin data and determine the mean and std in each bin
    binsize = 0.5
    iterations = 13
    #binlist = np.array([6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0])
    binlist = np.linspace(6.0, 12.0, 30)
    plotlist = np.array([6.25,6.75,7.25,7.75,8.5,9.5,10.25,10.75,11.25,11.75])
    
    meanRPS = np.empty(len(binlist))
    upperRPS = np.empty(len(binlist))
    lowerRPS = np.empty(len(binlist))
    meanKH = np.empty(len(binlist))
    upperKH = np.empty(len(binlist))
    lowerKH = np.empty(len(binlist))

    for i in range(len(binlist)):
        cut1 = (np.log10(ABmstar) > (binlist[i]-binsize)) & (np.log10(ABmstar) < (binlist[i]+binsize))

        meanRPS[i] = np.mean(fracRPS[cut1])
        upperRPS[i] = np.mean(fracRPS[cut1]) + np.std(fracRPS[cut1])
        lowerRPS[i] = np.mean(fracRPS[cut1]) - np.std(fracRPS[cut1])

        meanKH[i] = np.mean(fracKH[cut1])
        upperKH[i] = np.mean(fracKH[cut1]) + np.std(fracKH[cut1])
        lowerKH[i] = np.mean(fracKH[cut1]) - np.std(fracKH[cut1])

    meanRPS[np.where(meanRPS > 1.00)] = 1.00
    upperRPS[np.where(upperRPS > 1.00)] = 1.00
    lowerRPS[np.where(lowerRPS > 1.00)] = 1.00

    meanKH[np.where(meanKH > 1.00)] = 1.00
    upperKH[np.where(upperKH > 1.00)] = 1.00
    lowerKH[np.where(lowerKH > 1.00)] = 1.00
    


############
#Median and mean stripped fraction for object below 10^8 Msun
############

    cut = ABmstar < 10**8.0
    lowfrac = ABfrac[cut]
    print np.mean(lowfrac)
    print np.median(lowfrac)
    print np.std(lowfrac)


####################################
#Plot work
####################################

    axwidth = 3
    axlength = 10
    fontsize=28
    
    plt.rc('axes',linewidth=axwidth)
    plt.figure(figsize=(26,8))

####################
#FIGURE 1 - RPS
####################
    plt.subplot2grid((2,4),(0,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.6, hspace=0.1)
    plt.axis([5.8,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(9.3,0.9,r'$\rm {\it n}_{halo}\ = \ 10^{-'+np.str(0.1*hot)+'} cm^{-3}$', fontsize = 22)
    plt.text(9.3,0.8,r'$\rm {\it V}_{sat}\ = \ '+vel+'\ km\ s^{-1}$', fontsize = 22)

    ###Plotting routine for RPS scatter###
    ###
    if check == 'scatter':
        p1a = plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=None, fmt='yo', markersize = 10)
        #plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=ABerror, fmt='yo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        plt.plot(np.log10(ABmstar), ABfrac, color = 'Goldenrod', marker ='o', markersize = 16, linestyle = 'None')

    ###Plotting routine for RPS histogram###
    ###
    else:
        plt.fill_between(binlist, upperRPS, lowerRPS, facecolor = 'c', alpha = 0.4)
        plt.plot(binlist, meanRPS, color = 'DarkCyan', marker = 'None', linestyle = '-', linewidth = 4.0)
    
    xtickloc = [6,7,8,9,10,11]
    xtickstr = [r'$\rm 10^{6}$',r'$\rm 10^{7}$',r'$\rm 10^{8}$',r'$\rm 10^{9}$',r'$\rm 10^{10}$',r'$\rm 10^{11}$']
    ytickloc = [0.0,0.2,0.4,0.6,0.8,1.0]
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

####################
#FIGURE 2 - KH
####################
    plt.subplot2grid((2,4),(0,2), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.6, hspace=0.1)
    plt.axis([5.8,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(9.3,0.9,r'$\rm {\it n}_{halo}\ = \ 10^{-'+np.str(0.1*hot)+'} cm^{-3}$', fontsize = 22)
    plt.text(9.3,0.8,r'$\rm {\it V}_{sat}\ = \ '+vel+'\ km\ s^{-1}$', fontsize = 22)

    ###plotting routine for KH only scatter###
    ###
    if check == 'scatter':
        p1a = plt.errorbar(np.log10(ABmstar), fracKH, xerr = None, yerr=None, fmt='yo', markersize = 10)
        #plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=ABerror, fmt='yo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        plt.plot(np.log10(ABmstar), fracKH, color = 'Goldenrod', marker ='o', markersize = 16, linestyle = 'None')

    ###plotting routine for KH only histogram###
    ###
    else:
        plt.fill_between(binlist, upperKH, lowerKH, facecolor = 'c', alpha = 0.4)
        plt.plot(binlist, meanKH, color = 'DarkCyan', marker = 'None', linestyle = '-', linewidth = 4.0)
    
    xtickloc = [6,7,8,9,10,11]
    xtickstr = [r'$\rm 10^{6}$',r'$\rm 10^{7}$',r'$\rm 10^{8}$',r'$\rm 10^{9}$',r'$\rm 10^{10}$',r'$\rm 10^{11}$']
    ytickloc = [0.0,0.2,0.4,0.6,0.8,1.0]
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


####################################
# Plotting done: save, show and exit
####################################

    if check == 'normal':
        
        plt.savefig('RPS_motivation_hot'+np.str(hot)+'_v'+vel+'_AB_double.pdf')

    else:
        plt.savefig('RPS_motivation_hot'+np.str(hot)+'_v'+vel+'_AB_hist_double.pdf')

    plt.show()


###############################
## Hybrid
######################

def hybrid(hot1,vel1,hot2,vel2,tq,check):
    
    ABdata1 = Table.read('stripdata/RPS_strippedfrac_hot'+np.str(hot1)+'_v'+vel1+'_AB.dat', format = 'ascii')
    ABmstar1 = np.array(ABdata1['Mstar'])
    KH1 = np.array(ABdata1['KH_1Gyr'])
    ABdata2 = Table.read('stripdata/RPS_strippedfrac_hot'+np.str(hot2)+'_v'+vel2+'_AB.dat', format = 'ascii')
    ABmstar2 = np.array(ABdata2['Mstar'])
    KH2 = np.array(ABdata2['KH_1Gyr'])

    ABfrac1 = np.array(ABdata1['NFW'])
    ABfrac2 = np.array(ABdata2['NFW'])
    
    upperhot = hot1-5
    ABupperdata = Table.read('stripdata/RPS_strippedfrac_hot'+np.str(upperhot)+'_v'+vel1+'_AB.dat', format = 'ascii')
    ABupper = np.array(ABupperdata['NFW'])
    lowerhot = hot1+5
    ABlowerdata = Table.read('stripdata/RPS_strippedfrac_hot'+np.str(lowerhot)+'_v'+vel1+'_AB.dat', format = 'ascii')
    ABlower = np.array(ABlowerdata['NFW'])
    
    ABerror = np.empty([2,len(ABmstar1)])
    ABerror[1] = ABupper - ABfrac1
    ABerror[0] = ABfrac1 - ABlower

    if check == 'kh_hist':
        frac1 = tq*KH

    #elif check == 'scatter':
        #frac1 = ABfrac

    elif check == 'rps_hist':
        frac1 = ABfrac1
        frac2 = ABfrac2

    else:
        frac1 = ABfrac1 + tq*KH2
        
    frac1[np.where(frac1 > 1.0)] = 1.00

    #bin data and determine the mean and std in each bin
    binsize = 0.5
    iterations = 13
    #binlist = np.array([6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0])
    binlist = np.linspace(6.0, 12.0, 30)
    plotlist = np.array([6.25,6.75,7.25,7.75,8.5,9.5,10.25,10.75,11.25,11.75])
    
    mean1 = np.empty(len(binlist))
    upper1 = np.empty(len(binlist))
    lower1 = np.empty(len(binlist))

    for i in range(len(binlist)):
        cut1 = (np.log10(ABmstar1) > (binlist[i]-binsize)) & (np.log10(ABmstar1) < (binlist[i]+binsize))

        mean1[i] = np.mean(frac1[cut1])
        upper1[i] = np.mean(frac1[cut1]) + np.std(frac1[cut1])
        lower1[i] = np.mean(frac1[cut1]) - np.std(frac1[cut1])

    mean1[np.where(mean1 > 1.00)] = 1.00
    upper1[np.where(upper1 > 1.00)] = 1.00
    lower1[np.where(lower1 > 1.00)] = 1.00
    


############
#Median and mean stripped fraction for object below 10^8 Msun
############

    cut = ABmstar1 < 10**8.0
    lowfrac = ABfrac1[cut]
    print np.mean(lowfrac)
    print np.median(lowfrac)
    print np.std(lowfrac)


####################################
#Plot work
####################################

    axwidth = 3
    axlength = 10
    fontsize=28
    
    plt.rc('axes',linewidth=axwidth)
    
    plt.figure(figsize=(11,8))

############
#FIGURE 
############
    plt.subplot2grid((1,1),(0,0), colspan = 1, rowspan = 1)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.6, hspace=0.1)
    plt.axis([5.8,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(9.3,0.95,r'$\rm RPS:$', fontsize = 24)
    plt.text(9.4,0.9,r'$\rm {\it n}_{halo}\ = \ 10^{-'+np.str(0.1*hot1)+'} cm^{-3}$', fontsize = 22)
    plt.text(9.4,0.8,r'$\rm {\it V}_{sat}\ = \ '+vel1+'\ km\ s^{-1}$', fontsize = 22)

    plt.text(9.3,0.7,r'$\rm KH:$', fontsize = 24)
    plt.text(9.4,0.65,r'$\rm {\it n}_{halo}\ = \ 10^{-'+np.str(0.1*hot2)+'} cm^{-3}$', fontsize = 22)
    plt.text(9.4,0.55,r'$\rm {\it V}_{sat}\ = \ '+vel2+'\ km\ s^{-1}$', fontsize = 22)

    if check == 'scatter':
        p1a = plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=None, fmt='yo', markersize = 10)
        #plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=ABerror, fmt='yo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        plt.plot(np.log10(ABmstar), ABfrac, color = 'Goldenrod', marker ='o', markersize = 16, linestyle = 'None')

    elif check == 'kh_hist':
        #plotting routine for KH only histogram
        plt.fill_between(binlist, upper1, lower1, facecolor = 'y', alpha = 0.4)
        plt.plot(binlist, mean1, color = 'Goldenrod', marker = 'None', linestyle = '-', linewidth = 4.0)

    else:
        #plotting routine for alt
        plt.fill_between(binlist, upper1, lower1, facecolor = 'y', alpha = 0.4)
        plt.plot(binlist, mean1, color = 'Goldenrod', marker = 'None', linestyle = '-', linewidth = 4.0)
    
    xtickloc = [6,7,8,9,10,11]
    xtickstr = [r'$\rm 10^{6}$',r'$\rm 10^{7}$',r'$\rm 10^{8}$',r'$\rm 10^{9}$',r'$\rm 10^{10}$',r'$\rm 10^{11}$']
    ytickloc = [0.0,0.2,0.4,0.6,0.8,1.0]
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

    #if check == 'scatter':
        #plt.savefig('RPS_motivation_hot'+np.str(hot1)+'_v'+vel1+'_AB_hybrid.pdf')

    #elif check == 'kh_hist':
        #plt.savefig('RPS_motivation_hot'+np.str(hot1)+'_v'+vel1+'_AB_KHhist_hybrid.pdf')

    #else:
        #plt.savefig('RPS_motivation_hot'+np.str(hot1)+'_v'+vel1+'_AB_hist_hybrid.pdf')

    plt.show()
