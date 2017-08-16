#plotting routine for RPS paper based on get_h1data_fitparams script
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

########################################################################################

def figure(check):
    
    ABdata = Table.read('stripdata/RPS_strippedfrac_hot35_v300_AB.dat', format = 'ascii')
    ABmstar = np.array(ABdata['Mstar'])
    ABfrac = np.array(ABdata['NFW'])
    ABKH = np.array(ABdata['KH_1Gyr'])

    Bdata = Table.read('stripdata/RPS_strippedfrac_hot35_v300_Apace.dat', format = 'ascii')
    Bmstar = np.array(Bdata['Mstar'])
    Bfrac = np.array(Bdata['Burkert'])
    BKH = np.array(Bdata['KH_1Gyr'])    
    
    LTdata1 = Table.read('stripdata/RPS_strippedfrac_hot35_v300_LT.dat', format = 'ascii')
    LTmstar1 = np.array(LTdata1['Mstar'])
    LTfrac1 = np.array(LTdata1['NFW'])
    LTKH = np.array(LTdata1['KH_1Gyr'])
    #LTupper1 = np.array(LTdata1['upper_error'])
    #LTlower1 = np.array(LTdata1['lower_error'])

    #LTerror1 = np.empty([2,len(LTmstar1)])
    #LTerror1[1] = LTupper1 - LTfrac1
    #LTerror1[0] = LTfrac1 - LTlower1

    #ABfrac = ABfrac + ABKH
    #Bfrac = Bfrac + BKH
    #LTfrac1 = LTfrac1 + LTKH

    ABname = np.array(ABdata['NAME'])
    APacename = np.array(Bdata['NAME'])
    LTname = np.array(LTdata1['NAME'])
    ABmass = np.array(ABdata['Mstar'])
    APacemass = np.array(Bdata['Mstar'])
    LTmass = np.array(LTdata1['Mstar'])
    ABstrip = np.array(ABdata['NFW'])
    APacestrip = np.array(Bdata['Burkert'])
    LTstrip = np.array(LTdata1['NFW'])
    ABstrip_KH = np.array(ABdata['KH_1Gyr'])
    APacestrip_KH = np.array(Bdata['KH_1Gyr'])
    LTstrip_KH = np.array(LTdata1['KH_1Gyr'])

    #ABstrip = ABstrip + ABstrip_KH
    #LTstrip = LTstrip + LTstrip_KH
    #APacestrip = APacestrip + APacestrip_KH

    print len(LTname)
    print len(ABname)
    print len(APacename)

    diff1strip = np.array([])
    diff1mass = np.array([])
    diff2strip = np.array([])
    diff2mass = np.array([])

##################
# Diff #1: AB - LT
##################

    for i in range(len(LTname)):

        cut1 = LTname[i] == ABname
        diff = (ABstrip[cut1] - LTstrip[i])
        mass = LTmass[i]

        diff1strip = np.append(diff1strip,diff)
        diff1mass = np.append(diff1mass,mass)


#####################
# Diff #2: AB - Apace
#####################

    for i in range(len(APacename)):

        cut2 = APacename[i] == ABname
        diff = (ABstrip[cut2] - APacestrip[i])
        mass = APacemass[i]

        diff2strip = np.append(diff2strip,diff)
        diff2mass = np.append(diff2mass,mass)


######################
## Extra calculations
######################

    cut1 = diff1mass < 1e9
    cut2 = diff2mass < 1e9
    d1s = diff1strip[cut1]
    d2s = diff2strip[cut2]
    mean1 = np.mean(diff1strip[cut1])
    mean2 = np.mean(diff2strip[cut2])
    sig1 = np.std(diff1strip[cut1])
    sig2 = np.std(diff2strip[cut2])
    rms1 = np.sqrt(np.sum(d1s**2)/len(d1s))
    rms2 = np.sqrt(np.sum(d2s**2)/len(d2s))
    print 'mean1 = '+np.str(mean1)
    print 'mean2 = '+np.str(mean2)
    print 'rms1 = '+np.str(rms1)
    print 'rms2 = '+np.str(rms2)
    print 'sig1 = '+np.str(sig1)
    print 'sig2 = '+np.str(sig2)

####################################
#bin data and determine the mean and std in each bin
####################################
    binsize = 0.5
    iterations = 13
    #binlist = np.array([6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0])
    binlist = np.linspace(6.0, 12.0, 20)
    plotlist = np.array([6.25,6.75,7.25,7.75,8.5,9.5,10.25,10.75,11.25,11.75])

    LTmean1 = np.empty(len(binlist))
    LTupper1 = np.empty(len(binlist))
    LTlower1 = np.empty(len(binlist))
    Bmean1 = np.empty(len(binlist))
    Bupper1 = np.empty(len(binlist))
    Blower1 = np.empty(len(binlist))
    mean1 = np.empty(len(binlist))
    upper1 = np.empty(len(binlist))
    lower1 = np.empty(len(binlist))
    mean2 = np.empty(len(binlist))
    upper2 = np.empty(len(binlist))
    lower2 = np.empty(len(binlist))

    for i in range(len(binlist)):
        cutLT = (np.log10(LTmass) > (binlist[i]-binsize)) & (np.log10(LTmass) < (binlist[i]+binsize))
        cutB = (np.log10(APacemass) > (binlist[i]-binsize)) & (np.log10(APacemass) < (binlist[i]+binsize))
        cut1 = (np.log10(diff1mass) > (binlist[i]-binsize)) & (np.log10(diff1mass) < (binlist[i]+binsize))
        cut2 = (np.log10(diff2mass) > (binlist[i]-binsize)) & (np.log10(diff2mass) < (binlist[i]+binsize))

        LTmean1[i] = np.mean(LTstrip[cutLT])
        LTupper1[i] = np.mean(LTstrip[cutLT]) + np.std(LTstrip[cutLT])
        LTlower1[i] = np.mean(LTstrip[cutLT]) - np.std(LTstrip[cutLT])
        
        Bmean1[i] = np.mean(APacestrip[cutB])
        Bupper1[i] = np.mean(APacestrip[cutB]) + np.std(APacestrip[cutB])
        Blower1[i] = np.mean(APacestrip[cutB]) - np.std(APacestrip[cutB])
        
        mean1[i] = np.mean(diff1strip[cut1])
        upper1[i] = np.mean(diff1strip[cut1]) + np.std(diff1strip[cut1])
        lower1[i] = np.mean(diff1strip[cut1]) - np.std(diff1strip[cut1])
        
        mean2[i] = np.mean(diff2strip[cut2])
        upper2[i] = np.mean(diff2strip[cut2]) + np.std(diff2strip[cut2])
        lower2[i] = np.mean(diff2strip[cut2]) - np.std(diff2strip[cut2])

    LTmean1[np.where(LTmean1 > 1.00)] = 1.00
    LTupper1[np.where(LTupper1 > 1.00)] = 1.00
    LTlower1[np.where(LTlower1 > 1.00)] = 1.00
    
    Bmean1[np.where(Bmean1 > 1.00)] = 1.00
    Bupper1[np.where(Bupper1 > 1.00)] = 1.00
    Blower1[np.where(Blower1 > 1.00)] = 1.00
    
    mean1[np.where(mean1 > 1.00)] = 1.00
    upper1[np.where(upper1 > 1.00)] = 1.00
    lower1[np.where(lower1 > 1.00)] = 1.00
    
    mean2[np.where(mean2 > 1.00)] = 1.00
    upper2[np.where(upper2 > 1.00)] = 1.00
    lower2[np.where(lower2 > 1.00)] = 1.00



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
    #plt.subplot2grid((1,3),(0,0), colspan = 2, rowspan = 1)
    plt.subplot2grid((1,7),(0,0), colspan = 2, rowspan = 1)
    plt.subplots_adjust(left=0.07, bottom=0.13, right=0.93, top=0.97, wspace=0.2, hspace=0.1)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(9.5,0.92,r'$\rm NFW\ fits$', fontsize = 28)
    #plt.text(8.5,0.92,r'$\rm NFW:\ Oh\ et\ al.\ (2015)$', fontsize = 24)
    #plt.text(8.3,0.85,r'$\rm \rho_{halo}\ \sim \ 10^{-3.5} particles/cm^{3}$', fontsize = 20)
    #plt.text(8.3,0.75,r'$V_{\rm sat}\ \sim \ 300\ {\rm km/s}$', fontsize = 20)

    if check == 'scatter':
        p1b = plt.errorbar(np.log10(LTmstar1), LTfrac1, xerr = None, yerr=None, fmt='mh', markersize = 10)
        #plt.errorbar(np.log10(LTmstar1), LTfrac1, xerr = None, yerr=LTerror1, fmt='mh', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        plt.plot(np.log10(LTmstar1), LTfrac1, color = 'DarkMagenta', marker = 'H', markersize = 16, linestyle = 'None')

    else:
        plt.fill_between(binlist, LTupper1, LTlower1, facecolor = 'm', alpha = 0.4)
        plt.plot(binlist, LTmean1, color = 'DarkMagenta', marker = 'None', linestyle = '--', linewidth = 4.0)

    xtickloc = [7,8,9,10,11,11.5]
    xtickstr = [r'$\rm 10^{7}$',r'$\rm 10^{8}$',r'$\rm 10^{9}$',r'$\rm 10^{10}$',r'$\rm 10^{11}$']
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
    plt.subplot2grid((1,7),(0,2), colspan = 2, rowspan = 1)
    plt.subplots_adjust(left=0.07, bottom=0.13, right=0.93, top=0.97, wspace=0.0, hspace=0.1)
    plt.axis([6,11.5,-0.05,1.05])
    ax = plt.gca()
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Stripped\ Fraction$', fontsize = 28)

    plt.text(9.2,0.92,r'$\rm Burkert\ fits$', fontsize = 28)
    #plt.text(8.0,0.92,r'$\rm Burkert:\ Pace\ et\ al.\ (in\ prep)$', fontsize = 24)
    #plt.text(8.3,0.85,r'$\rm \rho_{halo}\ \sim \ 10^{-3.5} particles/cm^{3}$', fontsize = 20)
    #plt.text(8.3,0.75,r'$V_{\rm sat}\ \sim \ 300\ {\rm km/s}$', fontsize = 20)

    if check == 'scatter':

        p1a = plt.errorbar(np.log10(Bmstar), Bfrac, xerr = None, yerr=None, fmt='cd', markersize = 10)
        #plt.errorbar(np.log10(Bmstar), Bfrac, xerr = None, yerr=Berror, fmt='cd', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        plt.plot(np.log10(Bmstar), Bfrac, color = 'DarkCyan', marker = 'D', markersize = 12, linestyle = 'None')

    else:
        plt.fill_between(binlist, Bupper1, Blower1, facecolor = 'c', alpha = 0.4)
        plt.plot(binlist, Bmean1, color = 'DarkCyan', marker = 'None', linestyle = '-', linewidth = 4.0)
    
    xtickloc = [7,8,9,10,11]
    xtickstr = [r'$\rm 10^{7}$',r'$\rm 10^{8}$',r'$\rm 10^{9}$',r'$\rm 10^{10}$',r'$\rm 10^{11}$']
    ytickloc = [0.0,0.2,0.4,0.6,0.8,1.0]
    #ytickstr = ['$'+str(kk)+'$' for kk in ytickloc]
    
    ax.set_xticks(xtickloc)
    ax.set_xticklabels(xtickstr, position = (0,-0.01))
    ax.set_yticks(ytickloc)
    #ax.set_yticklabels(ytickstr)
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    #for tick in ax.yaxis.get_major_ticks():
        #tick.label1.set_fontsize(fontsize)
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(18)
        line.set_markeredgewidth(3)
    for tick in ax.xaxis.get_minor_ticks():
        tick.label1.set_fontsize(fontsize/2)
    #for tick in ax.yaxis.get_minor_ticks():
        #tick.label1.set_fontsize(fontsize/2)
    
    ax.tick_params(which='major',width=axwidth,length=axlength+5)
    ax.tick_params(which='minor',width=axwidth,length=axlength)


    ####################
    #FIGURE 3-3
    ####################

    plt.subplot2grid((1,7),(0,5), colspan = 2, rowspan = 1)
    plt.subplots_adjust(left=0.07, bottom=0.13, right=0.93, top=0.97, wspace=0.0, hspace=0.1)
    plt.axis([6,11.5,-0.95,0.35])
    ax = plt.gca()
    #ax.yaxis.tick_right()
    #ax.yaxis.set_label_position("right")
    #plt.axis([6,11.5,-0.95,0.45])
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm \Delta {\it f}_{stripped}$', fontsize = 28)

    #plt.text(8.5,0.8,r'$\rm Ab.Match - Oh et al.$', fontsize = 22)
    #plt.text(8.3,0.85,r'$\rm \rho_{hot}\ \sim \ 10^{-4.0} particles/cm^{3}$', fontsize = 22)
    #plt.text(8.3,0.75,r'$\rm V_{host}\ \sim \ 250\ km/s$', fontsize = 22)

    #plt.text(8.0,0.8,r'$\rm Ab.Match - APace$', fontsize = 22)
    #plt.text(8.3,0.85,r'$\rm \rho_{hot}\ \sim \ 10^{-4.0} particles/cm^{3}$', fontsize = 22)
    #plt.text(8.3,0.75,r'$\rm V_{host}\ \sim \ 250\ km/s$', fontsize = 22)


    plt.axhline(y = 0.0, xmin = 0.0, xmax = 1.0, linestyle = '-.', color = 'k', linewidth = 3.0)

    if check == 'scatter':

        p1a = plt.errorbar(np.log10(diff1mass), diff1strip, xerr = None, yerr=None, fmt='mh', markersize = 10)
        #plt.errorbar(np.log10(diff1mass), diff1strip, xerr = None, yerr=ABerror, fmt='mh', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        plt.plot(np.log10(diff1mass), diff1strip, color = 'DarkMagenta', marker = 'H', markersize = 16, linestyle = 'None')

        p2a = plt.errorbar(np.log10(diff2mass), diff2strip, xerr = None, yerr=None, fmt='cd', markersize = 10)
        #plt.errorbar(np.log10(diff2mass), diff2strip, xerr = None, yerr=ABerror, fmt='cd', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        plt.plot(np.log10(diff2mass), diff2strip, color = 'DarkCyan', marker = 'D', markersize = 12, linestyle = 'None')

    else:

        plt.fill_between(binlist, upper1, lower1, facecolor = 'm', alpha = 0.4)
        plt.plot(binlist, mean1, color = 'DarkMagenta', marker = 'None', linestyle = '--', linewidth = 4.0)

        plt.fill_between(binlist, upper2, lower2, facecolor = 'c', alpha = 0.4)
        plt.plot(binlist, mean2, color = 'DarkCyan', marker = 'None', linestyle = '-', linewidth = 4.0)

    xtickloc = [7,8,9,10,11]
    xtickstr = [r'$\rm 10^{7}$',r'$\rm 10^{8}$',r'$\rm 10^{9}$',r'$\rm 10^{10}$',r'$\rm 10^{11}$']
    ytickloc = [-0.8,-0.6,-0.4,-0.2,0.0,0.2]
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
        plt.savefig('RPS_3panel_massplusdiff.pdf')
        plt.savefig('RPS_3panel_massplusdiff.png')

    else:
        plt.savefig('RPS_3panel_massplusdiff_hist.pdf')
        plt.savefig('RPS_3panel_massplusdiff_hist.png')

    plt.show()
