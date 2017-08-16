#plotting routine for RPS paper based on get_h1data_fitparams script
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

########################################################################################

def figure1():
    
    ABdata = Table.read('stripdata/RPS_strippedfrac_hot40_v200_AB.dat', format = 'ascii')
    ABmstar = np.array(ABdata['Mstar'])
    ABfrac = np.array(ABdata['NFW'])
    ABupperdata = Table.read('stripdata/RPS_strippedfrac_hot35_v200_AB.dat', format = 'ascii')
    ABupper = np.array(ABupperdata['NFW'])
    ABlowerdata = Table.read('stripdata/RPS_strippedfrac_hot45_v200_AB.dat', format = 'ascii')
    ABlower = np.array(ABlowerdata['NFW'])
    
    ABerror = np.empty([2,len(ABmstar)])
    ABerror[1] = ABupper - ABfrac
    ABerror[0] = ABfrac - ABlower
    
    LTdata = Table.read('stripdata/RPS_strippedfrac_hot40_v200_LT.dat', format = 'ascii')
    LTmstar = np.array(LTdata['Mstar'])
    LTfrac = np.array(LTdata['NFW'])
    LTupper = np.array(LTdata['upper_error'])
    LTlower = np.array(LTdata['lower_error'])
    
    LTerror = np.empty([2,len(LTmstar)])
    LTerror[1] = LTupper - LTfrac
    LTerror[0] = LTfrac - LTlower
    
    Bdata = Table.read('stripdata/RPS_strippedfrac_hot40_v200_Apace.dat', format = 'ascii')
    Bmstar = np.array(Bdata['Mstar'])
    Bfrac = np.array(Bdata['Burkert'])
    Bupperdata = Table.read('stripdata/RPS_strippedfrac_hot35_v200_Apace.dat', format = 'ascii')
    Bupper = np.array(Bupperdata['Burkert'])
    Blowerdata = Table.read('stripdata/RPS_strippedfrac_hot45_v200_Apace.dat', format = 'ascii')
    Blower = np.array(Blowerdata['Burkert'])
    
    Berror = np.empty([2,len(Bmstar)])
    Berror[1] = Bupper - Bfrac
    Berror[0] = Bfrac - Blower

####################################
#Plot work
####################################

    axwidth = 3
    axlength = 10
    fontsize=28
    
    plt.rc('axes',linewidth=axwidth)
    
    plt.figure(figsize=(26,8))
    
    ####################
    #FIGURE 1a
    ####################
    plt.subplot2grid((2,6),(0,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.6, hspace=0.1)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)
    
    p1a = plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=None, fmt='co', markersize = 10)
    plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=ABerror, fmt='co', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(ABmstar), ABfrac, 'co', markersize = 10)

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
#FIGURE 1b
####################
    plt.subplot2grid((2,6),(0,2), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.6, hspace=0.1)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    p1b = plt.errorbar(np.log10(LTmstar), LTfrac, xerr = None, yerr=None, fmt='c^', markersize = 10)
    plt.errorbar(np.log10(LTmstar), LTfrac, xerr = None, yerr=LTerror, fmt='c^', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(LTmstar), LTfrac, 'c^', markersize = 10)
    
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
#FIGURE 1c
####################
    plt.subplot2grid((2,6),(0,4), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.6, hspace=0.1)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    p1a = plt.errorbar(np.log10(Bmstar), Bfrac, xerr = None, yerr=None, fmt='cd', markersize = 10)
    plt.errorbar(np.log10(Bmstar), Bfrac, xerr = None, yerr=Berror, fmt='cd', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(Bmstar), Bfrac, 'cd', markersize = 10)
    
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

    plt.savefig('RPS_figure1.pdf')

    plt.show()

########################################################################################
    
def figure2():
    
    ABdata = Table.read('stripdata/RPS_strippedfrac_hot30_v300_AB.dat', format = 'ascii')
    ABmstar = np.array(ABdata['Mstar'])
    ABfrac = np.array(ABdata['NFW'])
    ABupperdata = Table.read('stripdata/RPS_strippedfrac_hot25_v300_AB.dat', format = 'ascii')
    ABupper = np.array(ABupperdata['NFW'])
    ABlowerdata = Table.read('stripdata/RPS_strippedfrac_hot35_v300_AB.dat', format = 'ascii')
    ABlower = np.array(ABlowerdata['NFW'])
    
    ABerror = np.empty([2,len(ABmstar)])
    ABerror[1] = ABupper - ABfrac
    ABerror[0] = ABfrac - ABlower
    
    LTdata = Table.read('stripdata/RPS_strippedfrac_hot30_v300_LT.dat', format = 'ascii')
    LTmstar = np.array(LTdata['Mstar'])
    LTfrac = np.array(LTdata['NFW'])
    LTupper = np.array(LTdata['upper_error'])
    LTlower = np.array(LTdata['lower_error'])
    
    LTerror = np.empty([2,len(LTmstar)])
    LTerror[1] = LTupper - LTfrac
    LTerror[0] = LTfrac - LTlower
    
    Bdata = Table.read('stripdata/RPS_strippedfrac_hot30_v300_Apace.dat', format = 'ascii')
    Bmstar = np.array(Bdata['Mstar'])
    Bfrac = np.array(Bdata['Burkert'])
    Bupperdata = Table.read('stripdata/RPS_strippedfrac_hot25_v300_Apace.dat', format = 'ascii')
    Bupper = np.array(Bupperdata['Burkert'])
    Blowerdata = Table.read('stripdata/RPS_strippedfrac_hot35_v300_Apace.dat', format = 'ascii')
    Blower = np.array(Blowerdata['Burkert'])
    
    Berror = np.empty([2,len(Bmstar)])
    Berror[1] = Bupper - Bfrac
    Berror[0] = Bfrac - Blower
    
    ####################################
    #Plot work
    ####################################
    
    axwidth = 3
    axlength = 10
    fontsize=28
    
    plt.rc('axes',linewidth=axwidth)
    
    plt.figure(figsize=(26,8))
    
    ####################
    #FIGURE 2a
    ####################
    plt.subplot2grid((2,6),(0,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.6, hspace=0.1)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)
    
    p1a = plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=None, fmt='mo', markersize = 10)
    plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=ABerror, fmt='mo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(ABmstar), ABfrac, 'mo', markersize = 10)
    
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
    #FIGURE 2b
    ####################
    plt.subplot2grid((2,6),(0,2), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.6, hspace=0.1)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)
    
    p1b = plt.errorbar(np.log10(LTmstar), LTfrac, xerr = None, yerr=None, fmt='m^', markersize = 10)
    plt.errorbar(np.log10(LTmstar), LTfrac, xerr = None, yerr=LTerror, fmt='m^', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(LTmstar), LTfrac, 'm^', markersize = 10)
    
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
    #FIGURE 2c
    ####################
    plt.subplot2grid((2,6),(0,4), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.6, hspace=0.1)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)
    
    p1a = plt.errorbar(np.log10(Bmstar), Bfrac, xerr = None, yerr=None, fmt='md', markersize = 10)
    plt.errorbar(np.log10(Bmstar), Bfrac, xerr = None, yerr=Berror, fmt='md', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(Bmstar), Bfrac, 'md', markersize = 10)
    
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

    plt.savefig('RPS_figure2.pdf')
    
    plt.show()


########################################################################################

def figure3():
    
    ABdata = Table.read('stripdata/RPS_strippedfrac_hot40_v200_AB.dat', format = 'ascii')
    ABmstar = np.array(ABdata['Mstar'])
    ABfrac = np.array(ABdata['NFW'])
    ABupperdata = Table.read('stripdata/RPS_strippedfrac_hot40_v250_AB.dat', format = 'ascii')
    ABupper = np.array(ABupperdata['NFW'])
    ABlowerdata = Table.read('stripdata/RPS_strippedfrac_hot40_v150_AB.dat', format = 'ascii')
    ABlower = np.array(ABlowerdata['NFW'])
    
    ABerror = np.empty([2,len(ABmstar)])
    ABerror[1] = ABupper - ABfrac
    ABerror[0] = ABfrac - ABlower
    
    LTdata = Table.read('stripdata/RPS_strippedfrac_hot40_v200_LT.dat', format = 'ascii')
    LTmstar = np.array(LTdata['Mstar'])
    LTfrac = np.array(LTdata['NFW'])
    LTupper = np.array(LTdata['upper_error'])
    LTlower = np.array(LTdata['lower_error'])
    
    LTerror = np.empty([2,len(LTmstar)])
    LTerror[1] = LTupper - LTfrac
    LTerror[0] = LTfrac - LTlower
    
    Bdata = Table.read('stripdata/RPS_strippedfrac_hot40_v200_Apace.dat', format = 'ascii')
    Bmstar = np.array(Bdata['Mstar'])
    Bfrac = np.array(Bdata['Burkert'])
    Bupperdata = Table.read('stripdata/RPS_strippedfrac_hot40_v250_Apace.dat', format = 'ascii')
    Bupper = np.array(Bupperdata['Burkert'])
    Blowerdata = Table.read('stripdata/RPS_strippedfrac_hot40_v150_Apace.dat', format = 'ascii')
    Blower = np.array(Blowerdata['Burkert'])
    
    Berror = np.empty([2,len(Bmstar)])
    Berror[1] = Bupper - Bfrac
    Berror[0] = Bfrac - Blower

####################################
#Plot work
####################################

    axwidth = 3
    axlength = 10
    fontsize=28
    
    plt.rc('axes',linewidth=axwidth)
    
    plt.figure(figsize=(26,8))
    
    ####################
    #FIGURE 3a
    ####################
    plt.subplot2grid((2,6),(0,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.6, hspace=0.1)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)
    
    p1a = plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=None, fmt='bo', markersize = 10)
    plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=ABerror, fmt='bo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(ABmstar), ABfrac, 'bo', markersize = 10)

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
#FIGURE 3b
####################
    plt.subplot2grid((2,6),(0,2), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.6, hspace=0.1)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    p1b = plt.errorbar(np.log10(LTmstar), LTfrac, xerr = None, yerr=None, fmt='b^', markersize = 10)
    plt.errorbar(np.log10(LTmstar), LTfrac, xerr = None, yerr=LTerror, fmt='b^', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(LTmstar), LTfrac, 'b^', markersize = 10)
    
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
#FIGURE 3c
####################
    plt.subplot2grid((2,6),(0,4), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.6, hspace=0.1)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    p1a = plt.errorbar(np.log10(Bmstar), Bfrac, xerr = None, yerr=None, fmt='bd', markersize = 10)
    plt.errorbar(np.log10(Bmstar), Bfrac, xerr = None, yerr=Berror, fmt='bd', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(Bmstar), Bfrac, 'bd', markersize = 10)
    
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

    plt.savefig('RPS_figure3.pdf')

    plt.show()


########################################################################################

def figure4():
    
    ABdata = Table.read('stripdata/RPS_strippedfrac_hot35_v200_AB.dat', format = 'ascii')
    ABmstar = np.array(ABdata['Mstar'])
    ABfrac = np.array(ABdata['NFW'])
    ABupperdata = Table.read('stripdata/RPS_strippedfrac_hot30_v200_AB.dat', format = 'ascii')
    ABupper = np.array(ABupperdata['NFW'])
    ABlowerdata = Table.read('stripdata/RPS_strippedfrac_hot40_v200_AB.dat', format = 'ascii')
    ABlower = np.array(ABlowerdata['NFW'])
    
    ABerror = np.empty([2,len(ABmstar)])
    ABerror[1] = ABupper - ABfrac
    ABerror[0] = ABfrac - ABlower
    
    LTdata = Table.read('stripdata/RPS_strippedfrac_hot35_v200_LT.dat', format = 'ascii')
    LTmstar = np.array(LTdata['Mstar'])
    LTfrac = np.array(LTdata['NFW'])
    LTupper = np.array(LTdata['upper_error'])
    LTlower = np.array(LTdata['lower_error'])
    
    LTerror = np.empty([2,len(LTmstar)])
    LTerror[1] = LTupper - LTfrac
    LTerror[0] = LTfrac - LTlower
    
    Bdata = Table.read('stripdata/RPS_strippedfrac_hot35_v200_Apace.dat', format = 'ascii')
    Bmstar = np.array(Bdata['Mstar'])
    Bfrac = np.array(Bdata['Burkert'])
    Bupperdata = Table.read('stripdata/RPS_strippedfrac_hot30_v200_Apace.dat', format = 'ascii')
    Bupper = np.array(Bupperdata['Burkert'])
    Blowerdata = Table.read('stripdata/RPS_strippedfrac_hot40_v200_Apace.dat', format = 'ascii')
    Blower = np.array(Blowerdata['Burkert'])
    
    Berror = np.empty([2,len(Bmstar)])
    Berror[1] = Bupper - Bfrac
    Berror[0] = Bfrac - Blower

####################################
#Plot work
####################################

    axwidth = 3
    axlength = 10
    fontsize=28
    
    plt.rc('axes',linewidth=axwidth)
    
    plt.figure(figsize=(26,8))
    
    ####################
    #FIGURE 4a
    ####################
    plt.subplot2grid((2,6),(0,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.6, hspace=0.1)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)
    
    p1a = plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=None, fmt='ro', markersize = 10)
    plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=ABerror, fmt='ro', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(ABmstar), ABfrac, 'ro', markersize = 10)

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
#FIGURE 4b
####################
    plt.subplot2grid((2,6),(0,2), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.6, hspace=0.1)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    p1b = plt.errorbar(np.log10(LTmstar), LTfrac, xerr = None, yerr=None, fmt='r^', markersize = 10)
    plt.errorbar(np.log10(LTmstar), LTfrac, xerr = None, yerr=LTerror, fmt='r^', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(LTmstar), LTfrac, 'r^', markersize = 10)
    
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
#FIGURE 4c
####################
    plt.subplot2grid((2,6),(0,4), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.6, hspace=0.1)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    p1a = plt.errorbar(np.log10(Bmstar), Bfrac, xerr = None, yerr=None, fmt='rd', markersize = 10)
    plt.errorbar(np.log10(Bmstar), Bfrac, xerr = None, yerr=Berror, fmt='rd', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(Bmstar), Bfrac, 'rd', markersize = 10)
    
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

    plt.savefig('RPS_figure4.pdf')

    plt.show()
########################################################################################

def motivate(hot,vel,tq,check):
    
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

    if check == 'kh_hist':
        frac1 = tq*KH1

    else:
        frac1 = ABfrac + tq*KH
        
    frac1[np.where(frac1 > 1.0)] = 1.00

    binlist = np.array([6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 10.5, 11.0, 11.5, 12.0])
    plotlist = np.array([6.25,6.75,7.25,7.75,8.5,9.5,10.25,10.75,11.25,11.75])
    
    mean1 = np.empty(len(binlist)-1)
    upper1 = np.empty(len(binlist)-1)
    lower1 = np.empty(len(binlist)-1)

    for i in range(len(binlist)-1):
        cut1 = (np.log10(ABmstar) > binlist[i]) & (np.log10(ABmstar) < binlist[i+1])

        mean1[i] = np.mean(frac1[cut1])
        upper1[i] = np.mean(frac1[cut1]) + np.std(frac1[cut1])
        lower1[i] = np.mean(frac1[cut1]) - np.std(frac1[cut1])

    mean1[np.where(mean1 > 1.00)] = 1.00
    upper1[np.where(upper1 > 1.00)] = 1.00
    lower1[np.where(lower1 > 1.00)] = 1.00
    


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

####################
#FIGURE 
####################
    plt.subplot2grid((1,1),(0,0), colspan = 1, rowspan = 1)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.6, hspace=0.1)
    plt.axis([5.8,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(9.3,0.9,r'$\rm {\it n}_{halo}\ = \ 10^{-3.5} cm^{-3}$', fontsize = 22)
    plt.text(9.3,0.8,r'$\rm {\it V}_{sat}\ = \ 300\ km\ s^{-1}$', fontsize = 22)

    if check == 'normal':
        p1a = plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=None, fmt='yo', markersize = 10)
        #plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=ABerror, fmt='yo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        plt.plot(np.log10(ABmstar), ABfrac, color = 'Goldenrod', marker ='o', markersize = 16, linestyle = 'None')

    elif check == 'kh_hist':
        #plotting routine for KH only histogram
        plt.fill_between(plotlist, upper1, lower1, facecolor = 'y', alpha = 0.4)
        plt.plot(plotlist, mean1, color = 'Goldenrod', marker = 'None', linestyle = '-', linewidth = 4.0)

    else:
        #plotting routine for alt
        plt.fill_between(plotlist, upper1, lower1, facecolor = 'y', alpha = 0.4)
        plt.plot(plotlist, mean1, color = 'Goldenrod', marker = 'None', linestyle = '-', linewidth = 4.0)
    
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



    #plt.savefig('RPS_motivation_hot'+np.str(hot)+'_v'+vel+'_AB.pdf')
    #plt.savefig('RPS_motivation_hot'+np.str(hot)+'_v'+vel+'_AB_n.pdf')

    plt.show()

