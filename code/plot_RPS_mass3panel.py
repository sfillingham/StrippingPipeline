#plotting routine for RPS paper based on get_h1data_fitparams script
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

########################################################################################

def figure_3panel():
    
    ABdata = Table.read('stripdata/RPS_strippedfrac_hot35_v300_AB.dat', format = 'ascii')
    ABmstar = np.array(ABdata['Mstar'])
    ABfrac = np.array(ABdata['NFW'])
    ABfracKH = np.array(ABdata['KH_1Gyr'])

    ABfrac = ABfrac + ABfracKH

    Bdata = Table.read('stripdata/RPS_strippedfrac_hot35_v300_Apace.dat', format = 'ascii')
    Bmstar = np.array(Bdata['Mstar'])
    Bfrac = np.array(Bdata['Burkert'])
    BfracKH = np.array(Bdata['KH_1Gyr'])

    Bfrac = Bfrac + BfracKH
    
    LTdata = Table.read('stripdata/RPS_strippedfrac_hot35_v300_LT.dat', format = 'ascii')
    LTmstar = np.array(LTdata['Mstar'])
    LTfrac = np.array(LTdata['NFW'])
    #LTupper = np.array(LTdata['upper_error'])
    #LTlower = np.array(LTdata['lower_error'])

    #LTerror = np.empty([2,len(LTmstar)])
    #LTerror[1] = LTupper - LTfrac
    #LTerror[0] = LTfrac - LTlower


    ABname = np.array(ABdata['NAME'])
    APacename = np.array(Bdata['NAME'])
    LTname = np.array(LTdata['NAME'])
    ABmass = np.array(ABdata['Mstar'])
    APacemass = np.array(Bdata['Mstar'])
    LTmass = np.array(LTdata['Mstar'])
    ABstrip = np.array(ABdata['NFW'])
    APacestrip = np.array(Bdata['Burkert'])
    LTstrip = np.array(LTdata['NFW'])

    print len(LTname)
    print len(ABname)
    print len(APacename)

    diff1strip = np.array([])
    diff1mass = np.array([])
    diff2strip = np.array([])
    diff2mass = np.array([])
    diff3strip = np.array([])
    diff3mass = np.array([])

#####################
# Diff #1: AB - Apace
#####################

    for i in range(len(APacename)):

        cut1 = APacename[i] == ABname
        diff = ABstrip[cut1] - APacestrip[i]
        mass = APacemass[i]

        diff1strip = np.append(diff1strip,diff)
        diff1mass = np.append(diff1mass,mass)

##################
# Diff #2: AB - LT
##################

    for i in range(len(LTname)):

        cut2 = LTname[i] == ABname
        diff = ABstrip[cut2] - LTstrip[i]
        mass = LTmass[i]

        diff2strip = np.append(diff2strip,diff)
        diff2mass = np.append(diff2mass,mass)


####################################
#Plot work
####################################

    axwidth = 3
    axlength = 10
    fontsize=28
    
    plt.rc('axes',linewidth=axwidth)
    plt.figure(figsize=(26,8))

    ####################
    #FIGURE 3-1
    ####################
    plt.subplot2grid((1,3),(0,0), colspan = 1, rowspan = 1)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.2, hspace=0.1)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Stripped\ Fraction\ ({\it f}_{stripped})$', fontsize = 28)

    plt.text(8.2,0.92,r'$\rm NFW:\ Abundance\ Matching$', fontsize = 20)
    #plt.text(8.3,0.85,r'$\rm \rho_{halo}\ \sim \ 10^{-4.0} particles/cm^{3}$', fontsize = 20)
    #plt.text(8.3,0.75,r'$V_{\rm sat}\ \sim \ 250\ {\rm km/s}$', fontsize = 20)

    p1a = plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=None, fmt='co', markersize = 10)
    #plt.errorbar(np.log10(ABmstar), ABfrac, xerr = None, yerr=ABerror, fmt='co', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(ABmstar), ABfrac, color = 'Goldenrod', marker = 'o', markersize = 16, linestyle = 'None')
    
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
    plt.subplot2grid((1,3),(0,1), colspan = 1, rowspan = 1)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.2, hspace=0.1)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(8.5,0.92,r'$\rm NFW:\ Oh\ et\ al.\ (2015)$', fontsize = 20)
    #plt.text(8.3,0.85,r'$\rm \rho_{halo}\ \sim \ 10^{-3.5} particles/cm^{3}$', fontsize = 20)
    #plt.text(8.3,0.75,r'$V_{\rm sat}\ \sim \ 300\ {\rm km/s}$', fontsize = 20)
    
    p1b = plt.errorbar(np.log10(LTmstar), LTfrac, xerr = None, yerr=None, fmt='c^', markersize = 10)
    #plt.errorbar(np.log10(LTmstar1), LTfrac1, xerr = None, yerr=LTerror1, fmt='c^', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(LTmstar), LTfrac, color = 'Goldenrod', marker = '^', markersize = 16, linestyle = 'None')

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
    plt.subplot2grid((1,3),(0,2), colspan = 1, rowspan = 1)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.2, hspace=0.1)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(8.0,0.92,r'$\rm Burkert:\ Pace\ et\ al.\ (in\ prep)$', fontsize = 20)
    #plt.text(8.3,0.85,r'$\rm \rho_{halo}\ \sim \ 10^{-3.5} particles/cm^{3}$', fontsize = 20)
    #plt.text(8.3,0.75,r'$V_{\rm sat}\ \sim \ 300\ {\rm km/s}$', fontsize = 20)

    p1a = plt.errorbar(np.log10(Bmstar), Bfrac, xerr = None, yerr=None, fmt='cd', markersize = 10)
    #plt.errorbar(np.log10(Bmstar), Bfrac, xerr = None, yerr=Berror, fmt='cd', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(Bmstar), Bfrac, color = 'Goldenrod', marker = 'd', markersize = 16, linestyle = 'None')
    
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
    
    #plt.savefig('RPS_mass3panel_check.pdf')

    plt.show()
