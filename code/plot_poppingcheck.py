import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table


def work(check):

    sdata = Table.read('SHIELD_list.dat', format = 'ascii')
    #smstar = np.array(sdata['Mstar'])
    sHI = np.array(sdata['HI'])
    svmag = np.array(sdata['Mv'])
    smstar = 10**((4.83-(svmag))/2.5)

    LTdata = Table.read('AB_galaxylist.dat', format = 'ascii')
    #smstar = np.array(sdata['Mstar'])
    LTHI = np.array(LTdata['HI'])
    LTvmag = np.array(LTdata['Mv'])
    LTmstar = 10**((4.83-(LTvmag))/2.5)

    Tdata = Table.read('THINGS_galaxylist.dat', format = 'ascii')
    Tmstar = np.array(Tdata['Mstar'])
    THI = np.array(Tdata['HI'])
    #Tvmag = np.array(Tdata['Mv'])
    #Tmstar = 10**((4.83-(Tvmag))/2.5)

    ABdata = Table.read('stripdata/RPS_strippedfrac_hot35_v300_AB.dat', format = 'ascii')
    ABmstar = np.array(ABdata['Mstar'])
    ABfrac = np.array(ABdata['NFW'])
    rps_blue = np.array(ABdata['Bgasfrac'])
    rps_red = np.array(ABdata['Rgasfrac'])
    KH = np.array(ABdata['KH_1Gyr'])
    kh_red = np.array(ABdata['KH_red'])

    ABdata6 = Table.read('modeldata/RPS_strippedfrac_hot35_v300_Model.dat', format = 'ascii')
    ABmstar6 = np.array(ABdata6['Mstar'])
    ABfrac6 = np.array(ABdata6['NFW'])
    rps_blue6 = np.array(ABdata6['Bgasfrac'])
    rps_red6 = np.array(ABdata6['Rgasfrac'])
    KH6 = np.array(ABdata6['KH_1Gyr'])
    kh_red6 = np.array(ABdata6['KH_red'])

    pdata = Table.read('/Users/Sean/research/elvis/RamPressure/poppingdata/gas_candels.data/mstar.vs.mHI.dat', format = 'ascii')
    z = np.array(pdata['z'])
    mstar = np.array(pdata['Mstar'])
    mHI = np.array(pdata['HI50'])
    HIupper = np.array(pdata['HI84'])
    HIlower = np.array(pdata['HI16'])
    zcut = z == 0.75

    pdata1 = Table.read('/Users/Sean/research/elvis/RamPressure/poppingdata/gas_sham.data/mstar.vs.gas.dat', format = 'ascii')
    z1 = np.array(pdata1['z'])
    mstar1 = np.array(pdata1['Mstar'])
    mHI1 = np.array(pdata1['HI'])
    sigHI = np.array(pdata1['sigHI'])
    HIupper1 = mHI1+sigHI
    HIlower1 = mHI1-sigHI
    zcut0 = z1 == 0.0
    zcut1 = z1 == 0.5
    zcut2 = z1 == 1.0

    masterHI = np.append(LTHI, THI)
    masterHI = np.append(masterHI, sHI)

    #################################
    #bin the data by Mstar for alternative figures
    ###############
    binsize = 0.6
    iterations = 13
    #binlist = np.array([6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0])
    binlist = np.linspace(6.0, 11.0, 40)
    plotlist = np.array([6.25,6.75,7.25,7.75,8.5,9.5,10.25,10.75,11.25,11.75])
    
    median1 = np.empty(len(binlist))
    upper1 = np.empty(len(binlist))
    lower1 = np.empty(len(binlist))
    median2 = np.empty(len(binlist))
    upper2 = np.empty(len(binlist))
    lower2 = np.empty(len(binlist))

    gasmass = np.log10(rps_blue*ABmstar)

    for i in range(len(binlist)):
        cut1 = (np.log10(ABmstar) > (binlist[i]-binsize)) & (np.log10(ABmstar) < (binlist[i]+binsize))

        #mean1[i] = np.mean(rps_blue[cut1]*ABmstar[cut1])
        #upper1[i] = np.mean(rps_blue[cut1]*ABmstar[cut1]) + np.std(rps_blue[cut1]*ABmstar[cut1])
        #lower1[i] = np.mean(rps_blue[cut1]*ABmstar[cut1]) - np.std(rps_blue[cut1]*ABmstar[cut1])

        median2[i] = np.median(masterHI[cut1])
        upper2[i] = np.mean(masterHI[cut1]) + np.std(masterHI[cut1])
        lower2[i] = np.mean(masterHI[cut1]) - np.std(masterHI[cut1])

        median1[i] = np.median(gasmass[cut1])
        upper1[i] = np.mean(gasmass[cut1]) + np.std(gasmass[cut1])
        lower1[i] = np.mean(gasmass[cut1]) - np.std(gasmass[cut1])

    #mean1[np.where(mean1 > 1.00)] = 1.00
    #upper1[np.where(upper1 > 1.00)] = 1.00
    #lower1[np.where(lower1 > 1.00)] = 1.00

####################################
#Plot work
####################################

    axwidth = 3
    axlength = 10
    fontsize=28
    
    plt.rc('axes',linewidth=axwidth)
    plt.figure(figsize=(13,13))

    plt.subplot2grid((2,2),(0,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([4.8,11.2,4.8,11.2])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm log_{10}(M_{HI}/M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Gas\ Mass,\ (M_{\odot}))$', fontsize = 28)

    #p1a = plt.errorbar(np.log10(ABmstar6), np.log10(rps_blue6*ABmstar6), xerr = None, yerr=None, fmt='ro', markersize = 10)
    plt.plot(np.log10(ABmstar6), np.log10(rps_blue6*ABmstar6), color = 'Red', marker = 'None', linestyle = '-', markersize = 16, linewidth = 4.0)
    #plt.plot(np.log10(ABmstar1), ABfrac1, color = 'Goldenrod', marker = 'o', linestyle = 'None', markersize = 16)

    #plt.plot(np.log10(ABmstar), np.log10(rps_blue*ABmstar), color = 'DarkMagenta', marker = 'D', linestyle = 'None', markersize = 16)
    plt.plot(binlist, median1, color = 'DarkMagenta', linestyle = '-', linewidth = 4.0)
    plt.fill_between(binlist, lower1, upper1, color = 'DarkMagenta', edgecolor = 'DarkMagenta', alpha = 0.5)

    #plt.plot(np.log10(ABmstar), masterHI, color = 'Magenta', marker = 'D', linestyle = 'None', markersize = 16)
    plt.plot(binlist, median2, color = 'Green', linestyle = '-', linewidth = 4.0)
    plt.fill_between(binlist, lower2, upper2, color = 'Green', edgecolor = 'DarkMagenta', alpha = 0.5)

    #plt.plot(mstar[zcut], mHI[zcut], color = 'Black', marker = 'None', linestyle = '-', linewidth = 3.0)
    #plt.fill_between(mstar[zcut], HIlower[zcut], HIupper[zcut], color = 'Black', edgecolor = 'Black', alpha = 0.5)

    plt.plot(mstar1[zcut0], mHI1[zcut0], color = 'Blue', marker = 'None', linestyle = '-', linewidth = 3.0)
    plt.fill_between(mstar1[zcut0], HIlower1[zcut0], HIupper1[zcut0], color = 'Blue', edgecolor = 'Blue', alpha = 0.5)

    #plt.plot(mstar1[zcut1], mHI1[zcut1], color = 'Green', marker = 'None', linestyle = '-', linewidth = 3.0)
    #plt.fill_between(mstar1[zcut1], HIlower1[zcut1], HIupper1[zcut1], color = 'Green', edgecolor = 'Green', alpha = 0.5)

    #plt.plot(mstar1[zcut2], mHI1[zcut2], color = 'Yellow', marker = 'None', linestyle = '-', linewidth = 3.0)
    #plt.fill_between(mstar1[zcut1], HIlower1[zcut1], HIupper1[zcut1], color = 'Yellow', edgecolor = 'Yellow', alpha = 0.5)

    #LT data
    #plt.plot(np.log10(LTmstar), LTHI, color = 'DarkOrange', marker = 'o', linestyle = 'None', linewidth = 3.0, markersize = 16)
    #plt.fill_between(mstar1[zcut1], HIlower1[zcut1], HIupper1[zcut1], color = 'DarkOrange', edgecolor = 'Yellow', alpha = 0.5)

    #SHIELD data
    #plt.plot(np.log10(smstar), sHI, color = 'DarkCyan', marker = 'o', linestyle = 'None', linewidth = 3.0, markersize = 16)
    #plt.fill_between(mstar1[zcut1], HIlower1[zcut1], HIupper1[zcut1], color = 'DarkCyan', edgecolor = 'Yellow', alpha = 0.5)

    #THINGS data
    #plt.plot(Tmstar, THI, color = 'DarkGreen', marker = 'o', linestyle = 'None', linewidth = 3.0, markersize = 16)
    #plt.fill_between(mstar1[zcut1], HIlower1[zcut1], HIupper1[zcut1], color = 'DarkGreen', edgecolor = 'Yellow', alpha = 0.5)

    xtickloc = [6,7,8,9,10,11]
    xtickstr = [r'$\rm 10^{6}$',r'$\rm 10^{7}$',r'$\rm 10^{8}$',r'$\rm 10^{9}$',r'$\rm 10^{10}$',r'$\rm 10^{11}$']#,r'$\rm 10^{9}$',r'$\rm 10^{10}$',r'$\rm 10^{11}$']
    ytickloc = [5,6,7,8,9,10,11]
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

    plt.savefig('model_poppingcheck.pdf')

    masterHI = np.append(LTHI, THI)
    masterHI = np.append(masterHI, sHI)

    print len(masterHI)

    diff = gasmass - masterHI

    plt.figure(2)
    plt.plot(np.log10(ABmstar), diff, 'bo')
    plt.show()
