#plotting routine for RPS paper based on get_h1data_fitparams script
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

########################################################################################

def figure_9panel():
    
    ABdata = Table.read('stripdata/RPS_strippedfrac_hot40_v250_AB.dat', format = 'ascii')
    ABmstar = np.array(ABdata['Mstar'])
    ABfrac = np.array(ABdata['NFW'])

    Bdata = Table.read('stripdata/RPS_strippedfrac_hot40_v250_Apace.dat', format = 'ascii')
    Bmstar = np.array(Bdata['Mstar'])
    Bfrac = np.array(Bdata['Burkert'])
    
    LTdata1 = Table.read('stripdata/RPS_strippedfrac_hot40_v250_LT.dat', format = 'ascii')
    LTmstar1 = np.array(LTdata1['Mstar'])
    LTfrac1 = np.array(LTdata1['NFW'])
    LTupper1 = np.array(LTdata1['upper_error'])
    LTlower1 = np.array(LTdata1['lower_error'])

    LTerror1 = np.empty([2,len(LTmstar1)])
    LTerror1[1] = LTupper1 - LTfrac1
    LTerror1[0] = LTfrac1 - LTlower1

    LTdata4 = Table.read('stripdata/RPS_strippedfrac_hot45_v200_LT.dat', format = 'ascii')
    LTmstar4 = np.array(LTdata4['Mstar'])
    LTfrac4 = np.array(LTdata4['NFW'])
    LTupper4 = np.array(LTdata4['upper_error'])
    LTlower4 = np.array(LTdata4['lower_error'])
    
    LTerror4 = np.empty([2,len(LTmstar4)])
    LTerror4[1] = LTupper4 - LTfrac4
    LTerror4[0] = LTfrac4 - LTlower4

    LTdata5 = Table.read('stripdata/RPS_strippedfrac_hot45_v250_LT.dat', format = 'ascii')
    LTmstar5 = np.array(LTdata5['Mstar'])
    LTfrac5 = np.array(LTdata5['NFW'])
    LTupper5 = np.array(LTdata5['upper_error'])
    LTlower5 = np.array(LTdata5['lower_error'])
    
    LTerror5 = np.empty([2,len(LTmstar5)])
    LTerror5[1] = LTupper5 - LTfrac5
    LTerror5[0] = LTfrac5 - LTlower5

    LTdata6 = Table.read('stripdata/RPS_strippedfrac_hot45_v150_LT.dat', format = 'ascii')
    LTmstar6 = np.array(LTdata6['Mstar'])
    LTfrac6 = np.array(LTdata6['NFW'])
    LTupper6 = np.array(LTdata6['upper_error'])
    LTlower6 = np.array(LTdata6['lower_error'])
    
    LTerror6 = np.empty([2,len(LTmstar6)])
    LTerror6[1] = LTupper6 - LTfrac6
    LTerror6[0] = LTfrac6 - LTlower6

    LTdata7 = Table.read('stripdata/RPS_strippedfrac_hot35_v200_LT.dat', format = 'ascii')
    LTmstar7 = np.array(LTdata7['Mstar'])
    LTfrac7 = np.array(LTdata7['NFW'])
    LTupper7 = np.array(LTdata7['upper_error'])
    LTlower7 = np.array(LTdata7['lower_error'])
    
    LTerror7 = np.empty([2,len(LTmstar7)])
    LTerror7[1] = LTupper7 - LTfrac7
    LTerror7[0] = LTfrac7 - LTlower7

    LTdata8 = Table.read('stripdata/RPS_strippedfrac_hot35_v250_LT.dat', format = 'ascii')
    LTmstar8 = np.array(LTdata8['Mstar'])
    LTfrac8 = np.array(LTdata8['NFW'])
    LTupper8 = np.array(LTdata8['upper_error'])
    LTlower8 = np.array(LTdata8['lower_error'])
    
    LTerror8 = np.empty([2,len(LTmstar8)])
    LTerror8[1] = LTupper8 - LTfrac8
    LTerror8[0] = LTfrac8 - LTlower8

    LTdata9 = Table.read('stripdata/RPS_strippedfrac_hot35_v150_LT.dat', format = 'ascii')
    LTmstar9 = np.array(LTdata9['Mstar'])
    LTfrac9 = np.array(LTdata9['NFW'])
    LTupper9 = np.array(LTdata9['upper_error'])
    LTlower9 = np.array(LTdata9['lower_error'])
    
    LTerror9 = np.empty([2,len(LTmstar9)])
    LTerror9[1] = LTupper9 - LTfrac9
    LTerror9[0] = LTfrac9 - LTlower9
    

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

    p1b = plt.errorbar(np.log10(LTmstar1), LTfrac1, xerr = None, yerr=None, fmt='c^', markersize = 10)
    plt.errorbar(np.log10(LTmstar1), LTfrac1, xerr = None, yerr=LTerror1, fmt='c^', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(LTmstar1), LTfrac1, 'c^', markersize = 12)
    
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
    
    p1a = plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=None, fmt='co', markersize = 10)
    #plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=ABerror, fmt='co', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(ABmstar), ABfrac, 'co', markersize = 12)

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

    p1a = plt.errorbar(np.log10(Bmstar), Bfrac, xerr = None, yerr=None, fmt='cd', markersize = 10)
    #plt.errorbar(np.log10(Bmstar), Bfrac, xerr = None, yerr=Berror, fmt='cd', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(Bmstar), Bfrac, 'cd', markersize = 12)
    
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

    p1b = plt.errorbar(np.log10(LTmstar4), LTfrac4, xerr = None, yerr=None, fmt='c^', markersize = 10)
    plt.errorbar(np.log10(LTmstar4), LTfrac4, xerr = None, yerr=LTerror4, fmt='c^', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(LTmstar4), LTfrac4, 'c^', markersize = 12)
    
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

    p1b = plt.errorbar(np.log10(LTmstar5), LTfrac5, xerr = None, yerr=None, fmt='c^', markersize = 10)
    plt.errorbar(np.log10(LTmstar5), LTfrac5, xerr = None, yerr=LTerror5, fmt='c^', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(LTmstar5), LTfrac5, 'c^', markersize = 12)
    
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

    p1b = plt.errorbar(np.log10(LTmstar6), LTfrac6, xerr = None, yerr=None, fmt='c^', markersize = 10)
    plt.errorbar(np.log10(LTmstar6), LTfrac6, xerr = None, yerr=LTerror6, fmt='c^', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(LTmstar6), LTfrac6, 'c^', markersize = 12)
    
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

    p1b = plt.errorbar(np.log10(LTmstar7), LTfrac7, xerr = None, yerr=None, fmt='c^', markersize = 10)
    plt.errorbar(np.log10(LTmstar7), LTfrac7, xerr = None, yerr=LTerror7, fmt='c^', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(LTmstar7), LTfrac7, 'c^', markersize = 12)
    
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

    p1b = plt.errorbar(np.log10(LTmstar8), LTfrac8, xerr = None, yerr=None, fmt='c^', markersize = 10)
    plt.errorbar(np.log10(LTmstar8), LTfrac8, xerr = None, yerr=LTerror8, fmt='c^', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(LTmstar8), LTfrac8, 'c^', markersize = 12)
    
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

    p1b = plt.errorbar(np.log10(LTmstar9), LTfrac9, xerr = None, yerr=None, fmt='c^', markersize = 10)
    plt.errorbar(np.log10(LTmstar9), LTfrac9, xerr = None, yerr=LTerror9, fmt='c^', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(LTmstar9), LTfrac9, 'c^', markersize = 12)
    
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
    
    plt.savefig('RPS_figure3_MWboxes.pdf')

    plt.show()

########################################################################################

def figure4():
        
    LTdata1 = Table.read('stripdata/RPS_strippedfrac_hot30_v400_AB.dat', format = 'ascii')
    LTmstar1 = np.array(LTdata1['Mstar'])
    LTfrac1 = np.array(LTdata1['NFW'])
    #LTupper1 = np.array(LTdata1['upper_error'])
    #LTlower1 = np.array(LTdata1['lower_error'])
    
    #LTerror1 = np.empty([2,len(LTmstar1)])
    #LTerror1[1] = LTupper1 - LTfrac1
    #LTerror1[0] = LTfrac1 - LTlower1

    LTdata2 = Table.read('stripdata/RPS_strippedfrac_hot30_v500_AB.dat', format = 'ascii')
    LTmstar2 = np.array(LTdata2['Mstar'])
    LTfrac2 = np.array(LTdata2['NFW'])
    #LTupper2 = np.array(LTdata2['upper_error'])
    #LTlower2 = np.array(LTdata2['lower_error'])
    
    #LTerror2 = np.empty([2,len(LTmstar2)])
    #LTerror2[1] = LTupper2 - LTfrac2
    #LTerror2[0] = LTfrac2 - LTlower2

    LTdata3 = Table.read('stripdata/RPS_strippedfrac_hot30_v600_AB.dat', format = 'ascii')
    LTmstar3 = np.array(LTdata3['Mstar'])
    LTfrac3 = np.array(LTdata3['NFW'])
    #LTupper3 = np.array(LTdata3['upper_error'])
    #LTlower3 = np.array(LTdata3['lower_error'])
    
    #LTerror3 = np.empty([2,len(LTmstar3)])
    #LTerror3[1] = LTupper3 - LTfrac3
    #LTerror3[0] = LTfrac3 - LTlower3

    LTdata4 = Table.read('stripdata/RPS_strippedfrac_hot35_v400_AB.dat', format = 'ascii')
    LTmstar4 = np.array(LTdata4['Mstar'])
    LTfrac4 = np.array(LTdata4['NFW'])
    #LTupper4 = np.array(LTdata4['upper_error'])
    #LTlower4 = np.array(LTdata4['lower_error'])
    
    #LTerror4 = np.empty([2,len(LTmstar4)])
    #LTerror4[1] = LTupper4 - LTfrac4
    #LTerror4[0] = LTfrac4 - LTlower4

    LTdata5 = Table.read('stripdata/RPS_strippedfrac_hot35_v500_AB.dat', format = 'ascii')
    LTmstar5 = np.array(LTdata5['Mstar'])
    LTfrac5 = np.array(LTdata5['NFW'])
    #LTupper5 = np.array(LTdata5['upper_error'])
    #LTlower5 = np.array(LTdata5['lower_error'])
    
    #LTerror5 = np.empty([2,len(LTmstar5)])
    #LTerror5[1] = LTupper5 - LTfrac5
    #LTerror5[0] = LTfrac5 - LTlower5

    LTdata6 = Table.read('stripdata/RPS_strippedfrac_hot35_v600_AB.dat', format = 'ascii')
    LTmstar6 = np.array(LTdata6['Mstar'])
    LTfrac6 = np.array(LTdata6['NFW'])
    #LTupper6 = np.array(LTdata6['upper_error'])
    #LTlower6 = np.array(LTdata6['lower_error'])
    
    #LTerror6 = np.empty([2,len(LTmstar6)])
    #LTerror6[1] = LTupper6 - LTfrac6
    #LTerror6[0] = LTfrac6 - LTlower6


####################################
#Plot work
####################################

    axwidth = 3
    axlength = 10
    fontsize=28
    
    plt.rc('axes',linewidth=axwidth)
    plt.figure(figsize=(26,16))

    ####################
    #FIGURE 3-1
    ####################
    plt.subplot2grid((4,6),(0,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    #ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(6.3,0.15,r'$\rm \rho_{hot}\ \sim \ 10^{-3.0} particles/cm^{3}$', fontsize = 22)
    plt.text(6.3,0.05,r'$\rm V_{host}\ \sim \ 400\ km/s$', fontsize = 22)

    p1b = plt.errorbar(np.log10(LTmstar1), LTfrac1, xerr = None, yerr=None, fmt='b^', markersize = 10)
    #plt.errorbar(np.log10(LTmstar1), LTfrac1, xerr = None, yerr=LTerror1, fmt='b^', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(LTmstar1), LTfrac1, 'b^', markersize = 12)
    
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
    #FIGURE 3-2
    ####################
    plt.subplot2grid((4,6),(0,2), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    #ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(6.3,0.15,r'$\rm \rho_{hot}\ \sim \ 10^{-3.0} particles/cm^{3}$', fontsize = 22)
    plt.text(6.3,0.05,r'$\rm V_{host}\ \sim \ 500\ km/s$', fontsize = 22)
    
    p1a = plt.errorbar(np.log10(LTmstar2), LTfrac2, xerr = None, yerr=None, fmt='b^', markersize = 10)
    #plt.errorbar(np.log10(LTmstar2), LTfrac2, xerr = None, yerr=LTerror2, fmt='b^', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(LTmstar2), LTfrac2, 'b^', markersize = 12)

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
    plt.subplot2grid((4,6),(0,4), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    #ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(6.3,0.15,r'$\rm \rho_{hot}\ \sim \ 10^{-3.0} particles/cm^{3}$', fontsize = 22)
    plt.text(6.3,0.05,r'$\rm V_{host}\ \sim \ 600\ km/s$', fontsize = 22)

    p1a = plt.errorbar(np.log10(LTmstar3), LTfrac3, xerr = None, yerr=None, fmt='b^', markersize = 10)
    #plt.errorbar(np.log10(LTmstar3), LTfrac3, xerr = None, yerr=LTerror3, fmt='b^', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(LTmstar3), LTfrac3, 'b^', markersize = 12)
    
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
    plt.subplot2grid((4,6),(2,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(8.3,0.9,r'$\rm \rho_{hot}\ \sim \ 10^{-3.5} particles/cm^{3}$', fontsize = 22)
    plt.text(8.3,0.8,r'$\rm V_{host}\ \sim \ 400\ km/s$', fontsize = 22)

    p1b = plt.errorbar(np.log10(LTmstar4), LTfrac4, xerr = None, yerr=None, fmt='b^', markersize = 10)
    #plt.errorbar(np.log10(LTmstar4), LTfrac4, xerr = None, yerr=LTerror4, fmt='b^', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(LTmstar4), LTfrac4, 'b^', markersize = 12)
    
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
    plt.subplot2grid((4,6),(2,2), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(8.3,0.9,r'$\rm \rho_{hot}\ \sim \ 10^{-3.5} particles/cm^{3}$', fontsize = 22)
    plt.text(8.3,0.8,r'$\rm V_{host}\ \sim \ 500\ km/s$', fontsize = 22)

    p1b = plt.errorbar(np.log10(LTmstar5), LTfrac5, xerr = None, yerr=None, fmt='b^', markersize = 10)
    #plt.errorbar(np.log10(LTmstar5), LTfrac5, xerr = None, yerr=LTerror5, fmt='b^', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(LTmstar5), LTfrac5, 'b^', markersize = 12)
    
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
    plt.subplot2grid((4,6),(2,4), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(8.3,0.9,r'$\rm \rho_{hot}\ \sim \ 10^{-3.5} particles/cm^{3}$', fontsize = 22)
    plt.text(8.3,0.8,r'$\rm V_{host}\ \sim \ 600\ km/s$', fontsize = 22)

    p1b = plt.errorbar(np.log10(LTmstar6), LTfrac6, xerr = None, yerr=None, fmt='b^', markersize = 10)
    #plt.errorbar(np.log10(LTmstar6), LTfrac6, xerr = None, yerr=LTerror6, fmt='b^', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(LTmstar6), LTfrac6, 'b^', markersize = 12)
    
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
    
    plt.savefig('RPS_clusters.pdf')

    plt.show()
