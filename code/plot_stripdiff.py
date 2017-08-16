#This routine will take as input a single file from each of the 3 dark matter profile options
#and output the difference in stripped fraction for each

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

def work(hot, vel):

    ABdata = Table.read('stripdata/RPS_strippedfrac_hot'+hot+'_v'+vel+'_AB.dat', format = 'ascii')
    APacedata = Table.read('stripdata/RPS_strippedfrac_hot'+hot+'_v'+vel+'_Apace.dat', format = 'ascii')
    LTdata = Table.read('stripdata/RPS_strippedfrac_hot'+hot+'_v'+vel+'_LT.dat', format = 'ascii')

    ABname = np.array(ABdata['NAME'])
    APacename = np.array(APacedata['NAME'])
    LTname = np.array(LTdata['NAME'])
    ABmass = np.array(ABdata['Mstar'])
    APacemass = np.array(APacedata['Mstar'])
    LTmass = np.array(LTdata['Mstar'])
    ABstrip = np.array(ABdata['NFW'])
    APacestrip = np.array(APacedata['Burkert'])
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

#####################
# Diff #3: LT - Apace
#####################

    for i in range(len(APacename)):

        cut3 = APacename[i] == LTname
        diff = LTstrip[cut3] - APacestrip[i]
        mass = APacemass[i]

        if len(diff) != 0:
            diff3strip = np.append(diff3strip,diff)
            diff3mass = np.append(diff3mass,mass)

        else:
            continue

    print diff3strip
    
            

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
    plt.axis([6,11.5,-0.8,1.0])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Stripped\ Fraction\ Difference$', fontsize = 28)

    plt.text(8.5,0.8,r'$\rm Ab.Match - APace$', fontsize = 22)
    #plt.text(8.3,0.85,r'$\rm \rho_{hot}\ \sim \ 10^{-4.0} particles/cm^{3}$', fontsize = 22)
    #plt.text(8.3,0.75,r'$\rm V_{host}\ \sim \ 250\ km/s$', fontsize = 22)

    p1a = plt.errorbar(np.log10(diff1mass), diff1strip, xerr = None, yerr=None, fmt='yo', markersize = 10)
    #plt.errorbar(np.log10(diff1mass), diff1strip, xerr = None, yerr=ABerror, fmt='yo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(diff1mass), diff1strip, color = 'Goldenrod', marker = 'o', markersize = 16, linestyle = 'None')
    
    xtickloc = [6,7,8,9,10,11]
    xtickstr = [r'$\rm 10^{6}$',r'$\rm 10^{7}$',r'$\rm 10^{8}$',r'$\rm 10^{9}$',r'$\rm 10^{10}$',r'$\rm 10^{11}$']
    ytickloc = [-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0]
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
    plt.subplot2grid((1,3),(0,1), colspan = 1, rowspan = 1)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.2, hspace=0.1)
    plt.axis([6,11.5,-0.8,1.0])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction\ Difference$', fontsize = 28)

    plt.text(8.0,0.8,r'$\rm Ab.Match - LittleTHINGS$', fontsize = 22)
    #plt.text(8.3,0.85,r'$\rm \rho_{hot}\ \sim \ 10^{-4.0} particles/cm^{3}$', fontsize = 22)
    #plt.text(8.3,0.75,r'$\rm V_{host}\ \sim \ 250\ km/s$', fontsize = 22)

    p2a = plt.errorbar(np.log10(diff2mass), diff2strip, xerr = None, yerr=None, fmt='yo', markersize = 10)
    #plt.errorbar(np.log10(diff2mass), diff2strip, xerr = None, yerr=ABerror, fmt='yo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(diff2mass), diff2strip, color = 'Goldenrod', marker = 'o', markersize = 16, linestyle = 'None')
    
    xtickloc = [6,7,8,9,10,11]
    xtickstr = [r'$\rm 10^{6}$',r'$\rm 10^{7}$',r'$\rm 10^{8}$',r'$\rm 10^{9}$',r'$\rm 10^{10}$',r'$\rm 10^{11}$']
    ytickloc = [-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0]
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
    plt.axis([6,11.5,-0.8,1.0])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction\ Difference$', fontsize = 28)

    plt.text(8.5,0.8,r'$\rm LittleTHINGS - APace$', fontsize = 22)
    #plt.text(8.3,0.85,r'$\rm \rho_{hot}\ \sim \ 10^{-4.0} particles/cm^{3}$', fontsize = 22)
    #plt.text(8.3,0.75,r'$\rm V_{host}\ \sim \ 250\ km/s$', fontsize = 22)

    p3a = plt.errorbar(np.log10(diff3mass), diff3strip, xerr = None, yerr=None, fmt='yo', markersize = 10)
    #plt.errorbar(np.log10(diff3mass), diff3strip, xerr = None, yerr=ABerror, fmt='yo', markersize = 5, elinewidth = 2.5, capsize = 10.0)
    plt.plot(np.log10(diff3mass), diff3strip, color = 'Goldenrod', marker = 'o', markersize = 16, linestyle = 'None')
    
    xtickloc = [6,7,8,9,10,11]
    xtickstr = [r'$\rm 10^{6}$',r'$\rm 10^{7}$',r'$\rm 10^{8}$',r'$\rm 10^{9}$',r'$\rm 10^{10}$',r'$\rm 10^{11}$']
    ytickloc = [-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0]
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
    
    plt.savefig('stripdiff_fid.pdf')

    plt.show()
    
