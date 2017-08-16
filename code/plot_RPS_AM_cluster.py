#plotting routine for RPS paper based on get_h1data_fitparams script
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

########################################################################################

def work(tq,check):
    
    ABdata1 = Table.read('stripdata/RPS_strippedfrac_hot30_v400_AB.dat', format = 'ascii')
    ABmstar1 = np.array(ABdata1['Mstar'])
    ABfrac1 = np.array(ABdata1['NFW'])
    KH1 = np.array(ABdata1['KH_1Gyr'])

    ABdata2 = Table.read('stripdata/RPS_strippedfrac_hot30_v800_AB.dat', format = 'ascii')
    ABmstar2 = np.array(ABdata2['Mstar'])
    ABfrac2 = np.array(ABdata2['NFW'])
    KH2 = np.array(ABdata2['KH_1Gyr'])

    ABdata3 = Table.read('stripdata/RPS_strippedfrac_hot35_v400_AB.dat', format = 'ascii')
    ABmstar3 = np.array(ABdata3['Mstar'])
    ABfrac3 = np.array(ABdata3['NFW'])
    KH3 = np.array(ABdata3['KH_1Gyr'])

    ABdata4 = Table.read('stripdata/RPS_strippedfrac_hot35_v800_AB.dat', format = 'ascii')
    ABmstar4 = np.array(ABdata4['Mstar'])
    ABfrac4 = np.array(ABdata4['NFW'])
    KH4 = np.array(ABdata4['KH_1Gyr'])

    if check == 'kh_hist':
        frac1 = tq*KH1
        frac2 = tq*KH2
        frac3 = tq*KH3
        frac4 = tq*KH4

    else:
        frac1 = ABfrac1 + tq*KH1
        frac2 = ABfrac2 + tq*KH2
        frac3 = ABfrac3 + tq*KH3
        frac4 = ABfrac4 + tq*KH4

    frac1[np.where(frac1 > 1.0)] = 1.00
    frac2[np.where(frac2 > 1.0)] = 1.00
    frac3[np.where(frac3 > 1.0)] = 1.00
    frac4[np.where(frac4 > 1.0)] = 1.00

    #################################
    #bin the data by Mstar for alternative figures
    ###############
    binsize = 0.5
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
    mean3 = np.empty(len(binlist))
    upper3 = np.empty(len(binlist))
    lower3 = np.empty(len(binlist))
    mean4 = np.empty(len(binlist))
    upper4 = np.empty(len(binlist))
    lower4 = np.empty(len(binlist))

    for i in range(len(binlist)):
        cut1 = (np.log10(ABmstar1) > (binlist[i]-binsize)) & (np.log10(ABmstar1) < (binlist[i]+binsize))
        cut2 = (np.log10(ABmstar2) > (binlist[i]-binsize)) & (np.log10(ABmstar2) < (binlist[i]+binsize))
        cut3 = (np.log10(ABmstar3) > (binlist[i]-binsize)) & (np.log10(ABmstar3) < (binlist[i]+binsize))
        cut4 = (np.log10(ABmstar4) > (binlist[i]-binsize)) & (np.log10(ABmstar4) < (binlist[i]+binsize))

        mean1[i] = np.mean(frac1[cut1])
        upper1[i] = np.mean(frac1[cut1]) + np.std(frac1[cut1])
        lower1[i] = np.mean(frac1[cut1]) - np.std(frac1[cut1])
        mean2[i] = np.mean(frac2[cut2])
        upper2[i] = np.mean(frac2[cut2]) + np.std(frac2[cut2])
        lower2[i] = np.mean(frac2[cut2]) - np.std(frac2[cut2])
        mean3[i] = np.mean(frac3[cut3])
        upper3[i] = np.mean(frac3[cut3]) + np.std(frac3[cut3])
        lower3[i] = np.mean(frac3[cut3]) - np.std(frac3[cut3])
        mean4[i] = np.mean(frac4[cut4])
        upper4[i] = np.mean(frac4[cut4]) + np.std(frac4[cut4])
        lower4[i] = np.mean(frac4[cut4]) - np.std(frac4[cut4])

    mean1[np.where(mean1 > 1.00)] = 1.00
    mean2[np.where(mean2 > 1.00)] = 1.00
    mean3[np.where(mean3 > 1.00)] = 1.00
    mean4[np.where(mean4 > 1.00)] = 1.00
    upper1[np.where(upper1 > 1.00)] = 1.00
    lower1[np.where(lower1 > 1.00)] = 1.00
    upper2[np.where(upper2 > 1.00)] = 1.00
    lower2[np.where(lower2 > 1.00)] = 1.00
    upper3[np.where(upper3 > 1.00)] = 1.00
    lower3[np.where(lower3 > 1.00)] = 1.00
    upper4[np.where(upper4 > 1.00)] = 1.00
    lower4[np.where(lower4 > 1.00)] = 1.00


####################################
#Plot work
####################################

    axwidth = 3
    axlength = 10
    fontsize=28
    
    plt.rc('axes',linewidth=axwidth)
    plt.figure(figsize=(16,16))

    ####################
    #FIGURE 2-1
    ####################
    plt.subplot2grid((4,4),(0,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([5.8,11.5,-0.05,1.05])
    ax = plt.gca()
    #ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(9.1,0.9,r'$\rm {\it n}_{halo}\ = \ 10^{-3.0} cm^{-3}$', fontsize = 22)
    plt.text(9.1,0.8,r'$\rm {\it V}_{sat}\ = \ 400\ km\ s^{-1}$', fontsize = 22)

    if check == 'scatter':
        p1b = plt.errorbar(np.log10(ABmstar1), ABfrac1, xerr = None, yerr=None, fmt='co', markersize = 10)
        #plt.errorbar(np.log10(ABmstar1), ABfrac1, xerr = None, yerr=LTerror9, fmt='co', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        plt.plot(np.log10(ABmstar1), ABfrac1, color = 'DarkCyan', marker = 'o', linestyle = 'None', markersize = 16)

    elif check == 'kh_hist':
        #plotting routine for KH only histogram
        plt.fill_between(binlist, upper1, lower1, facecolor = 'c', alpha = 0.4)
        plt.plot(binlist, mean1, color = 'DarkCyan', marker = 'None', linestyle = '-', linewidth = 3.0)

    else:
        #plotting routine for alt
        plt.fill_between(binlist, upper1, lower1, facecolor = 'c', alpha = 0.4)
        plt.plot(binlist, mean1, color = 'DarkCyan', marker = 'None', linestyle = '-', linewidth = 3.0)
    
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
    #for tick in ax.xaxis.get_minor_ticks():
        #tick.label1.set_fontsize(fontsize/2)
    for tick in ax.yaxis.get_minor_ticks():
        tick.label1.set_fontsize(fontsize/2)
    
    ax.tick_params(which='major',width=axwidth,length=axlength+5)
    ax.tick_params(which='minor',width=axwidth,length=axlength)
    
    ####################
    #FIGURE 2-2
    ####################
    plt.subplot2grid((4,4),(0,2), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([5.8,11.5,-0.05,1.05])
    ax = plt.gca()
    #ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(6.3,0.12,r'$\rm {\it n}_{halo}\ = \ 10^{-3.0} cm^{-3}$', fontsize = 22)
    #plt.text(6.3,0.12,r'$\rm \rho_{halo}\ = \ 10^{-3.0} cm^{-3}$', fontsize = 22)
    plt.text(6.3,0.02,r'$\rm {\it V}_{sat}\ = \ 800\ km\ s^{-1}$', fontsize = 22)

    if check == 'scatter':
        p1b = plt.errorbar(np.log10(ABmstar2), ABfrac2, xerr = None, yerr=None, fmt='mo', markersize = 10)
        #plt.errorbar(np.log10(ABmstar2), ABfrac2, xerr = None, yerr=LTerror9, fmt='mo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        plt.plot(np.log10(ABmstar2), ABfrac2, color = 'DarkMagenta', marker = 'o', linestyle = 'None', markersize = 16)

    elif check == 'kh_hist':
        #plotting routine for KH only histogram
        plt.fill_between(binlist, upper2, lower2, facecolor = 'm', alpha = 0.4)
        plt.plot(binlist, mean2, color = 'DarkMagenta', marker = 'None', linestyle = '-', linewidth = 3.0)

    else:
        #plotting routine for alt
        plt.fill_between(binlist, upper2, lower2, facecolor = 'm', alpha = 0.4)
        plt.plot(binlist, mean2, color = 'DarkMagenta', marker = 'None', linestyle = '-', linewidth = 3.0)

    xtickloc = [6,7,8,9,10,11,11.5]
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
    #FIGURE 2-3
    ####################
    plt.subplot2grid((4,4),(2,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([5.8,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(9.1,0.9,r'$\rm {\it n}_{halo}\ = \ 10^{-3.5} cm^{-3}$', fontsize = 22)
    #plt.text(9.1,0.9,r'$\rm \rho_{halo}\ = \ 10^{-3.5} cm^{-3}$', fontsize = 22)
    plt.text(9.1,0.8,r'$\rm {\it V}_{sat}\ = \ 400\ km\ s^{-1}$', fontsize = 22)

    if check == 'scatter':
        p1b = plt.errorbar(np.log10(ABmstar3), ABfrac3, xerr = None, yerr=None, fmt='co', markersize = 10)
        #plt.errorbar(np.log10(ABmstar3), ABfrac3, xerr = None, yerr=LTerror9, fmt='co', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        plt.plot(np.log10(ABmstar3), ABfrac3, color = 'DarkCyan', marker = 'o', linestyle = 'None', markersize = 16)

    elif check == 'kh_hist':
        #plotting routine for KH only histogram
        plt.fill_between(binlist, upper3, lower3, facecolor = 'c', alpha = 0.4)
        plt.plot(binlist, mean3, color = 'DarkCyan', marker = 'None', linestyle = '-', linewidth = 3.0)

    else:
        #plotting routine for alt
        plt.fill_between(binlist, upper3, lower3, facecolor = 'c', alpha = 0.4)
        plt.plot(binlist, mean3, color = 'DarkCyan', marker = 'None', linestyle = '-', linewidth = 3.0)
    
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
    #FIGURE 2-4
    ####################
    plt.subplot2grid((4,4),(2,2), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([5.8,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(6.3,0.12,r'$\rm {\it n}_{halo}\ = \ 10^{-3.5} cm^{-3}$', fontsize = 22)
    #plt.text(6.3,0.12,r'$\rm \rho_{halo}\ = \ 10^{-3.5} cm^{-3}$', fontsize = 22)
    plt.text(6.3,0.02,r'$\rm {\it V}_{sat}\ = \ 800\ km\ s^{-1}$', fontsize = 22)

    if check == 'scatter':
        p1b = plt.errorbar(np.log10(ABmstar4), ABfrac4, xerr = None, yerr=None, fmt='mo', markersize = 10)
        #plt.errorbar(np.log10(ABmstar4), ABfrac4, xerr = None, yerr=LTerror9, fmt='mo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        plt.plot(np.log10(ABmstar4), ABfrac4, color = 'DarkMagenta', marker = 'o', linestyle = 'None', markersize = 16)

    elif check == 'kh_hist':
        #plotting routine for KH only histogram
        plt.fill_between(binlist, upper4, lower4, facecolor = 'm', alpha = 0.4)
        plt.plot(binlist, mean4, color = 'DarkMagenta', marker = 'None', linestyle = '-', linewidth = 3.0)

    else:
        #plotting routine for alt
        plt.fill_between(binlist, upper4, lower4, facecolor = 'm', alpha = 0.4)
        plt.plot(binlist, mean4, color = 'DarkMagenta', marker = 'None', linestyle = '-', linewidth = 3.0)
    
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


###################################################
# End of Plotting, time to save and display figure
###################################################

    if check == 'scatter':
        plt.savefig('RPS_AM_cluster.pdf')

    elif check == 'kh_hist':
        plt.savefig('RPS_AM_cluster_KHonly_hist'+np.str(tq)+'Gyr.pdf')

    else:
        plt.savefig('RPS_AM_cluster_KHonly_hist'+np.str(tq)+'Gyr.pdf')
        

    plt.show()


########################################################################################