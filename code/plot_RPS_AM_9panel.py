#plotting routine for RPS paper based on get_h1data_fitparams script
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

########################################################################################

def figure_9panel():
    
    ABdata1 = Table.read('stripdata/RPS_strippedfrac_hot40_v200_AB.dat', format = 'ascii')
    ABmstar1 = np.array(ABdata1['Mstar'])
    ABfrac1 = np.array(ABdata1['NFW'])

    ABdata2 = Table.read('stripdata/RPS_strippedfrac_hot40_v250_AB.dat', format = 'ascii')
    ABmstar2 = np.array(ABdata2['Mstar'])
    ABfrac2 = np.array(ABdata2['NFW'])

    ABdata3 = Table.read('stripdata/RPS_strippedfrac_hot40_v300_AB.dat', format = 'ascii')
    ABmstar3 = np.array(ABdata3['Mstar'])
    ABfrac3 = np.array(ABdata3['NFW'])

    ABdata4 = Table.read('stripdata/RPS_strippedfrac_hot35_v200_AB.dat', format = 'ascii')
    ABmstar4 = np.array(ABdata4['Mstar'])
    ABfrac4 = np.array(ABdata4['NFW'])

    ABdata5 = Table.read('stripdata/RPS_strippedfrac_hot35_v250_AB.dat', format = 'ascii')
    ABmstar5 = np.array(ABdata5['Mstar'])
    ABfrac5 = np.array(ABdata5['NFW'])

    ABdata6 = Table.read('stripdata/RPS_strippedfrac_hot35_v300_AB.dat', format = 'ascii')
    ABmstar6 = np.array(ABdata6['Mstar'])
    ABfrac6 = np.array(ABdata6['NFW'])

    ABdata7 = Table.read('stripdata/RPS_strippedfrac_hot45_v200_AB.dat', format = 'ascii')
    ABmstar7 = np.array(ABdata7['Mstar'])
    ABfrac7 = np.array(ABdata7['NFW'])

    ABdata8 = Table.read('stripdata/RPS_strippedfrac_hot45_v250_AB.dat', format = 'ascii')
    ABmstar8 = np.array(ABdata8['Mstar'])
    ABfrac8 = np.array(ABdata8['NFW'])

    ABdata9 = Table.read('stripdata/RPS_strippedfrac_hot45_v300_AB.dat', format = 'ascii')
    ABmstar9 = np.array(ABdata9['Mstar'])
    ABfrac9 = np.array(ABdata9['NFW'])

####################################
#Plot work
####################################

    axwidth = 3
    axlength = 10
    fontsize=28
    
    plt.rc('axes',linewidth=axwidth)
    plt.figure(figsize=(26,24))

    ####################
    #FIGURE 3-1
    ####################
    plt.subplot2grid((6,6),(0,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    #ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(8.3,0.9,r'$\rm \rho_{hot}\ \sim \ 10^{-4.0} particles/cm^{3}$', fontsize = 22)
    plt.text(8.3,0.8,r'$\rm V_{host}\ \sim \ 200\ km/s$', fontsize = 22)

    p1b = plt.errorbar(np.log10(ABmstar1), ABfrac1, xerr = None, yerr=None, fmt='co', markersize = 10)
    #plt.errorbar(np.log10(ABmstar1), ABfrac1, xerr = None, yerr=LTerror9, fmt='co', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(ABmstar1), ABfrac1, color = 'DarkCyan', marker = 'o', linestyle = 'None', markersize = 16)
    
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
    #FIGURE 3-2
    ####################
    plt.subplot2grid((6,6),(0,2), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    #ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(8.3,0.9,r'$\rm \rho_{hot}\ \sim \ 10^{-4.0} particles/cm^{3}$', fontsize = 22)
    plt.text(8.3,0.8,r'$\rm V_{host}\ \sim \ 250\ km/s$', fontsize = 22)
    
    p1b = plt.errorbar(np.log10(ABmstar2), ABfrac2, xerr = None, yerr=None, fmt='mo', markersize = 10)
    #plt.errorbar(np.log10(ABmstar2), ABfrac2, xerr = None, yerr=LTerror9, fmt='mo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(ABmstar2), ABfrac2, color = 'DarkMagenta', marker = 'o', linestyle = 'None', markersize = 16)

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
    #FIGURE 3-3
    ####################
    plt.subplot2grid((6,6),(0,4), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    #ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(8.3,0.9,r'$\rm \rho_{hot}\ \sim \ 10^{-4.0} particles/cm^{3}$', fontsize = 22)
    plt.text(8.3,0.8,r'$\rm V_{host}\ \sim \ 300\ km/s$', fontsize = 22)

    p1b = plt.errorbar(np.log10(ABmstar3), ABfrac3, xerr = None, yerr=None, fmt='yo', markersize = 10)
    #plt.errorbar(np.log10(ABmstar3), ABfrac3, xerr = None, yerr=LTerror9, fmt='yo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(ABmstar3), ABfrac3, color = 'Goldenrod', marker = 'o', linestyle = 'None', markersize = 16)
    
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
    #FIGURE 3-4
    ####################
    plt.subplot2grid((6,6),(2,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    #ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(8.3,0.9,r'$\rm \rho_{hot}\ \sim \ 10^{-3.5} particles/cm^{3}$', fontsize = 22)
    plt.text(8.3,0.8,r'$\rm V_{host}\ \sim \ 200\ km/s$', fontsize = 22)

    p1b = plt.errorbar(np.log10(ABmstar4), ABfrac4, xerr = None, yerr=None, fmt='co', markersize = 10)
    #plt.errorbar(np.log10(ABmstar4), ABfrac4, xerr = None, yerr=LTerror9, fmt='co', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(ABmstar4), ABfrac4, color = 'DarkCyan', marker = 'o', linestyle = 'None', markersize = 16)
    
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
    #FIGURE 3-5
    ####################
    plt.subplot2grid((6,6),(2,2), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    #ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(8.3,0.9,r'$\rm \rho_{hot}\ \sim \ 10^{-3.5} particles/cm^{3}$', fontsize = 22)
    plt.text(8.3,0.8,r'$\rm V_{host}\ \sim \ 250\ km/s$', fontsize = 22)

    p1b = plt.errorbar(np.log10(ABmstar5), ABfrac5, xerr = None, yerr=None, fmt='mo', markersize = 10)
    #plt.errorbar(np.log10(ABmstar5), ABfrac5, xerr = None, yerr=LTerror9, fmt='mo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(ABmstar5), ABfrac5, color = 'DarkMagenta', marker = 'o', linestyle = 'None', markersize = 16)
    
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
    #FIGURE 3-6
    ####################
    plt.subplot2grid((6,6),(2,4), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    #ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(8.3,0.9,r'$\rm \rho_{hot}\ \sim \ 10^{-3.5} particles/cm^{3}$', fontsize = 22)
    plt.text(8.3,0.8,r'$\rm V_{host}\ \sim \ 300\ km/s$', fontsize = 22)

    p1b = plt.errorbar(np.log10(ABmstar6), ABfrac6, xerr = None, yerr=None, fmt='yo', markersize = 10)
    #plt.errorbar(np.log10(ABmstar6), ABfrac6, xerr = None, yerr=LTerror9, fmt='yo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(ABmstar6), ABfrac6, color = 'Goldenrod', marker = 'o', linestyle = 'None', markersize = 16)
    
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
    #FIGURE 3-7
    ####################
    plt.subplot2grid((6,6),(4,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(8.3,0.9,r'$\rm \rho_{hot}\ \sim \ 10^{-4.5} particles/cm^{3}$', fontsize = 22)
    plt.text(8.3,0.8,r'$\rm V_{host}\ \sim \ 200\ km/s$', fontsize = 22)

    p1b = plt.errorbar(np.log10(ABmstar7), ABfrac7, xerr = None, yerr=None, fmt='co', markersize = 10)
    #plt.errorbar(np.log10(ABmstar7), ABfrac7, xerr = None, yerr=LTerror9, fmt='co', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(ABmstar7), ABfrac7, color = 'DarkCyan', marker = 'o', linestyle = 'None', markersize = 16)
    
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
    #FIGURE 3-8
    ####################
    plt.subplot2grid((6,6),(4,2), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(8.3,0.9,r'$\rm \rho_{hot}\ \sim \ 10^{-4.5} particles/cm^{3}$', fontsize = 22)
    plt.text(8.3,0.8,r'$\rm V_{host}\ \sim \ 250\ km/s$', fontsize = 22)

    p1b = plt.errorbar(np.log10(ABmstar8), ABfrac8, xerr = None, yerr=None, fmt='mo', markersize = 10)
    #plt.errorbar(np.log10(ABmstar8), ABfrac8, xerr = None, yerr=LTerror9, fmt='mo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(ABmstar8), ABfrac8, color = 'DarkMagenta', marker = 'o', linestyle = 'None', markersize = 16)
    
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
    #FIGURE 3-9
    ####################
    plt.subplot2grid((6,6),(4,4), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(8.3,0.9,r'$\rm \rho_{hot}\ \sim \ 10^{-4.5} particles/cm^{3}$', fontsize = 22)
    plt.text(8.3,0.8,r'$\rm V_{host}\ \sim \ 300\ km/s$', fontsize = 22)

    p1b = plt.errorbar(np.log10(ABmstar9), ABfrac9, xerr = None, yerr=None, fmt='yo', markersize = 10)
    #plt.errorbar(np.log10(ABmstar9), ABfrac9, xerr = None, yerr=LTerror9, fmt='yo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(ABmstar9), ABfrac9, color = 'Goldenrod', marker = 'o', linestyle = 'None', markersize = 16)
    
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
    
    plt.savefig('RPS_AM_9panel.pdf')

    plt.show()


########################################################################################
