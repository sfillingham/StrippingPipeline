import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from astropy.table import Table


def work(check):

    ABdata1 = Table.read('stripdata/RPS_strippedfrac_hot40_v200_AB.dat', format = 'ascii')
    ABmstar1 = np.array(ABdata1['Mstar'])
    ABfrac1 = np.array(ABdata1['NFW'])
    rps_blue1 = np.array(ABdata1['Bgasfrac'])
    rps_red1 = np.array(ABdata1['Rgasfrac'])
    kh_red1 = np.array(ABdata1['KH_red'])

    ABdata2 = Table.read('stripdata/RPS_strippedfrac_hot40_v250_AB.dat', format = 'ascii')
    ABmstar2 = np.array(ABdata2['Mstar'])
    ABfrac2 = np.array(ABdata2['NFW'])
    rps_blue2 = np.array(ABdata2['Bgasfrac'])
    rps_red2 = np.array(ABdata2['Rgasfrac'])
    kh_red2 = np.array(ABdata2['KH_red'])

    ABdata3 = Table.read('stripdata/RPS_strippedfrac_hot40_v300_AB.dat', format = 'ascii')
    ABmstar3 = np.array(ABdata3['Mstar'])
    ABfrac3 = np.array(ABdata3['NFW'])
    rps_blue3 = np.array(ABdata3['Bgasfrac'])
    rps_red3 = np.array(ABdata3['Rgasfrac'])
    kh_red3 = np.array(ABdata3['KH_red'])

    ABdata4 = Table.read('stripdata/RPS_strippedfrac_hot35_v200_AB.dat', format = 'ascii')
    ABmstar4 = np.array(ABdata4['Mstar'])
    ABfrac4 = np.array(ABdata4['NFW'])
    rps_blue4 = np.array(ABdata4['Bgasfrac'])
    rps_red4 = np.array(ABdata4['Rgasfrac'])
    kh_red4 = np.array(ABdata4['KH_red'])

    ABdata5 = Table.read('stripdata/RPS_strippedfrac_hot35_v250_AB.dat', format = 'ascii')
    ABmstar5 = np.array(ABdata5['Mstar'])
    ABfrac5 = np.array(ABdata5['NFW'])
    rps_blue5 = np.array(ABdata5['Bgasfrac'])
    rps_red5 = np.array(ABdata5['Rgasfrac'])
    kh_red5 = np.array(ABdata5['KH_red'])

    ABdata6 = Table.read('stripdata/RPS_strippedfrac_hot30_v300_AB.dat', format = 'ascii')
    ABmstar6 = np.array(ABdata6['Mstar'])
    ABfrac6 = np.array(ABdata6['NFW'])
    rps_blue6 = np.array(ABdata6['Bgasfrac'])
    rps_red6 = np.array(ABdata6['Rgasfrac'])
    KH6 = np.array(ABdata6['KH_1Gyr'])
    kh_red6 = np.array(ABdata6['KH_red'])

    #######################################
    # TCBFlash data
    #######################################

    data = Table.read('dwarf_data.dat', format = 'ascii')
    fielddata = Table.read('fielddwarf_data.dat', format = 'ascii')

    galname = np.array(data['Galname'])
    stellarmass = np.array(data['M_star'])
    HImass = np.array(data['M_HI'])

    fieldname = np.array(fielddata['Galname'])
    fieldstellarmass = np.array(fielddata['M_star'])
    fieldHImass = np.array(fielddata['M_HI'])

    if check == 'gasfrac':
        limitcut = (HImass < 0) & (HImass > -1)
        limit_h1frac = -1*HImass[limitcut]/stellarmass[limitcut]
        #limit_dist = dist[limitcut]
        limit_mstar = stellarmass[limitcut]

        #MWdist = dist[1:5]
        MWmstar = stellarmass[1:5]
        MW_h1frac = HImass[1:5]/stellarmass[1:5]
        #Adist = dist[5:12]
        Amstar = stellarmass[5:12]
        A_h1frac = HImass[5:12]/stellarmass[5:12]
        MWlabel = galname[:5]
        Alabel = galname[5:]
        Aname = galname[5:]
        starcut = (galname == 'IC10')
        starcut1 = (galname == 'NGC6822')
        starcut2 = (galname == 'PegDIG')

        #Newdist = dist[13:16]
        New_mstar = stellarmass[13:16]
        New_h1frac = HImass[13:16]/stellarmass[13:16]

        #MWLowdist = dist[16]
        MWLowmstar = stellarmass[16]
        MWLow_h1frac = HImass[16]/stellarmass[16]
        #ALowdist = dist[17:19]
        ALowmstar = stellarmass[17:19]
        ALow_h1frac = HImass[17:19]/stellarmass[17:19]
        #StarLowdist = dist[19]
        StarLowmstar = stellarmass[19]
        StarLow_h1frac = HImass[19]/stellarmass[19]
        #Phoenixdist = dist[20]
        Phoenixmstar = stellarmass[20]
        Phoenix_h1frac = HImass[20]/stellarmass[20]

    else:
        limitcut = (HImass < 0) & (HImass > -1)
        limit_h1frac = -1*HImass[limitcut]
        #limit_dist = dist[limitcut]
        limit_mstar = stellarmass[limitcut]

        #MWdist = dist[1:5]
        MWmstar = stellarmass[1:5]
        MW_h1frac = HImass[1:5]
        #Adist = dist[5:12]
        Amstar = stellarmass[5:12]
        A_h1frac = HImass[5:12]
        MWlabel = galname[:5]
        Alabel = galname[5:]
        Aname = galname[5:]
        starcut = (galname == 'IC10')
        starcut1 = (galname == 'NGC6822')
        starcut2 = (galname == 'PegDIG')

        #Newdist = dist[13:16]
        New_mstar = stellarmass[13:16]
        New_h1frac = HImass[13:16]

        #MWLowdist = dist[16]
        MWLowmstar = stellarmass[16]
        MWLow_h1frac = HImass[16]
        #ALowdist = dist[17:19]
        ALowmstar = stellarmass[17:19]
        ALow_h1frac = HImass[17:19]
        #StarLowdist = dist[19]
        StarLowmstar = stellarmass[19]
        StarLow_h1frac = HImass[19]
        #Phoenixdist = dist[20]
        Phoenixmstar = stellarmass[20]
        Phoenix_h1frac = HImass[20]

    
    #########################################
    # Calculations for Section 4 in RPS paper
    #########################################
    cut6 = ABmstar6 <= 100000000.
    ABfrac = ABfrac6[cut6]
    KHfrac = KH6[cut6]
    rps_blue6 = rps_blue6[cut6]
    rps_red6 = rps_red6[cut6]
    kh_red6 = kh_red6[cut6]
    ABmstar6 = ABmstar6[cut6]

    totalrem = ABfrac + KHfrac
    totalrem[np.where(totalrem > 1.0)[0]] = 1.0
    
    fracHI = (1-(totalrem))*rps_blue6
    print 'RPS+KH'
    print np.mean(fracHI)
    print np.std(fracHI)
    kh_stripped = fracHI[np.where(fracHI <= 0.136)[0]]
    print len(kh_stripped)
    print len(fracHI)
    fracHI[np.where(fracHI < 0.0001)[0]] = 0.001
    print np.min(fracHI)

    print 'RPS only'
    print np.mean(rps_red6)
    #print np.std(rps_red6)
    rps_stripped = rps_red6[np.where(rps_red6 <= 0.136)[0]]
    print len(rps_stripped)
    print len(rps_red6)

####################################
#Plot work
####################################

    axwidth = 4
    axlength = 10
    fontsize=28

    #############################################################################################
    ## First figure, this is 'direct comparison' plot for the Spekkens 14 paper
    #############################################################################################
    
    plt.rc('axes',linewidth=axwidth)
    plt.figure(figsize=(13,13))

    if check == 'gasfrac':

        plt.subplot2grid((2,2),(0,0), colspan = 2, rowspan = 2)
        plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
        plt.axis([5.8,8.2,-5,2.1])
        ax = plt.gca()
        ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
        ax.set_ylabel(r'$\rm Gas\ Fraction,\ M_{ISM}/M_{*}$', fontsize = 28)
        #ax.set_ylabel(r'$\rm Gas\ Mass,\ (M_{\odot}))$', fontsize = 28)

    

        plt.plot(np.log10(MWmstar*10**(6)), np.log10(MW_h1frac), linestyle = 'None', marker='D', color = 'white', markersize=16.0, mec='DarkMagenta', mew=3.0)
        plt.plot(np.log10(Amstar*10**(6)), np.log10(A_h1frac), linestyle = 'None', marker='D', color = 'white', markersize=16.0, mec='DarkMagenta', mew=3.0)
        #plt.plot(Newdist, np.log10(New_h1frac), linestyle = 'None', marker='D', color = 'black', markersize=8.0, mec='black', mew=2.0)
        plt.plot(np.log10(stellarmass[starcut]*10**(6)), np.log10(HImass[starcut]/stellarmass[starcut]), linestyle = 'None', marker='o', color = 'white', markersize=16.0, label = 'IC10', mec='DarkCyan', mew=3.0)
        plt.plot(np.log10(stellarmass[starcut1]*10**(6)), np.log10(HImass[starcut1]/stellarmass[starcut1]), linestyle = 'None', marker='o', color = 'white', markersize=16.0, label = 'NGC6822', mec='DarkCyan', mew=3.0)
        plt.plot(np.log10(stellarmass[starcut2]*10**(6)), np.log10(HImass[starcut2]/stellarmass[starcut2]), linestyle = 'None', marker='o', color = 'white', markersize=16.0, label = 'PegDIG', mec='DarkCyan', mew=3.0)
        plt.plot(np.log10(fieldstellarmass*10**(6)), np.log10(fieldHImass/fieldstellarmass), linestyle = 'None', marker='o', color = 'white', markersize=16.0, label = 'Field', mec='DarkCyan', mew=3.0)
        #plt.plot(MWLowdist, np.log10(MWLow_h1frac), linestyle = 'None', marker='D', fillstyle='none', color = 'red', markersize=8.0, mew=2.0)
        #plt.plot(ALowdist, np.log10(ALow_h1frac), linestyle = 'None', marker='D', fillstyle='none', color = 'red', markersize=8.0, mew=2.0)
        #plt.plot(StarLowdist, np.log10(StarLow_h1frac), linestyle = 'None', marker='D', fillstyle='none', color = 'blue', markersize=8.0, mew=2.0)
        #plt.plot(Phoenixdist, np.log10(Phoenix_h1frac), linestyle = 'None', marker='o', fillstyle='none', color = 'blue', markersize=8.0, mew=2.0)

        plt.plot(np.log10(limit_mstar[0]*10**(6)), np.log10(limit_h1frac[0]), linestyle = 'None', marker='D', color = 'white', markersize=16.0, mec='DarkMagenta', mew=3.0)
        plt.plot(np.log10(limit_mstar[1:]*10**(6)), np.log10(limit_h1frac[1:]), linestyle = 'None', marker='D', color = 'white', markersize=16.0, mec='DarkMagenta', mew=3.0)

        plt.arrow(np.log10(21*10**(6)), np.log10(0.007/21), 0.0, -0.5, head_width=0.05, head_length=0.1, fc='DarkMagenta', ec='DarkMagenta')
        plt.arrow(np.log10(3.9*10**(6)), np.log10(0.35/3.9), 0.0, -0.5, head_width=0.05, head_length=0.1, fc='DarkMagenta', ec='DarkMagenta')
        plt.arrow(np.log10(1.1*10**(6)), np.log10(0.1/1.1), 0.0, -0.5, head_width=0.05, head_length=0.1, fc='DarkMagenta', ec='DarkMagenta')
        plt.arrow(np.log10(62*10**(6)), np.log10(0.005/62), 0.0, -0.5, head_width=0.05, head_length=0.1, fc='DarkMagenta', ec='DarkMagenta')
        plt.arrow(np.log10(7.6*10**(6)), np.log10(0.27/7.6), 0.0, -0.5, head_width=0.05, head_length=0.1, fc='DarkMagenta', ec='DarkMagenta')
        plt.arrow(np.log10(9.5*10**(6)), np.log10(0.36/9.5), 0.0, -0.5, head_width=0.05, head_length=0.1, fc='DarkMagenta', ec='DarkMagenta')
        plt.arrow(np.log10(2.8*10**(6)), np.log10(0.015/2.8), 0.0, -0.5, head_width=0.05, head_length=0.1, fc='DarkMagenta', ec='DarkMagenta')

        plt.errorbar(np.log10(ABmstar6), np.log10(rps_blue6), xerr = None, yerr=None, fmt='co', markersize = 10)
        #plt.errorbar(np.log10(ABmstar1), ABfrac1, xerr = None, yerr=LTerror9, fmt='co', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        p1 = plt.plot(np.log10(ABmstar6), np.log10(rps_blue6), color = 'DarkCyan', marker = 'o', linestyle = 'None', markersize = 16)
        #plt.plot(np.log10(ABmstar1), ABfrac1, color = 'Goldenrod', marker = 'o', linestyle = 'None', markersize = 16)

        #p1b = plt.errorbar(np.log10(ABmstar6), np.log10(rps_red6), xerr = None, yerr=None, fmt='ro', markersize = 10)
        #plt.errorbar(np.log10(ABmstar1), ABfrac1, xerr = None, yerr=LTerror9, fmt='co', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        #plt.plot(np.log10(ABmstar6), np.log10(rps_red6), color = 'DarkRed', marker = 'o', linestyle = 'None', markersize = 16)
        #plt.plot(np.log10(ABmstar1), ABfrac1, color = 'Goldenrod', marker = 'o', linestyle = 'None', markersize = 16)

        plt.errorbar(np.log10(ABmstar6), np.log10(fracHI), xerr = None, yerr=None, fmt='mD', markersize = 10)
        #plt.errorbar(np.log10(ABmstar1), ABfrac1, xerr = None, yerr=LTerror9, fmt='co', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        plt.plot(np.log10(ABmstar6), np.log10(fracHI), color = 'DarkMagenta', marker = 'D', linestyle = 'None', markersize = 16)
        #plt.plot(np.log10(ABmstar1), ABfrac1, color = 'Goldenrod', marker = 'o', linestyle = 'None', markersize = 16)

        xtickloc = [6,7,8]
        xtickstr = [r'$\rm 10^{6}$',r'$\rm 10^{7}$',r'$\rm 10^{8}$']#,r'$\rm 10^{9}$',r'$\rm 10^{10}$',r'$\rm 10^{11}$']
        ytickloc = [-4,-3,-2,-1,0,1,2]
        ytickstr = ['$'+str(kk)+'$' for kk in ytickloc]

    else:
        plt.subplot2grid((2,2),(0,0), colspan = 2, rowspan = 2)
        plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
        plt.axis([5.8,8.2,3.5,9.5])
        ax = plt.gca()
        ax.set_xlabel(r'$\rm Stellar\ Mass\ (M_{\odot})$', fontsize = 28)
        #ax.set_xlabel(r'$\rm log_{10}(M_{*}/M_{\odot})$', fontsize = 28)
        #ax.set_ylabel(r'$\rm log_{10}(M_{HI}/M_{\odot})$', fontsize = 28)
        ax.set_ylabel(r'$\rm HI \ Mass\ (M_{\odot})$', fontsize = 28)

        p1a = plt.plot(np.log10(stellarmass[starcut]*10**(6)), np.log10(HImass[starcut]*10**(6)), linestyle='None', marker='o', color='white', markersize=16.0, mec='DarkCyan', mew=3.0, label=r"$\rm Local\ Field$")
        p2a = plt.plot(np.log10(MWmstar*10**(6)), np.log10(MW_h1frac*10**(6)), linestyle='None', marker='D', color='white', markersize=16.0, mec='DarkMagenta', mew=3.0, label=r"$\rm Local\ Group\ Satellites$")
        plt.plot(np.log10(Amstar*10**(6)), np.log10(A_h1frac*10**(6)), linestyle='None', marker='D', color='white', markersize=16.0, mec='DarkMagenta', mew=3.0)
        #plt.plot(Newdist, np.log10(New_h1frac), linestyle = 'None', marker='D', color = 'black', markersize=8.0, mec='black', mew=2.0)
        
        plt.plot(np.log10(stellarmass[starcut1]*10**(6)), np.log10(HImass[starcut1]*10**(6)), linestyle='None', marker='o', color='white', markersize=16.0, mec='DarkCyan', mew=3.0)
        plt.plot(np.log10(stellarmass[starcut2]*10**(6)), np.log10(HImass[starcut2]*10**(6)), linestyle='None', marker='o', color='white', markersize=16.0, mec='DarkCyan', mew=3.0)
        plt.plot(np.log10(fieldstellarmass*10**(6)), np.log10(fieldHImass*10**(6)), linestyle='None', marker='o', color='white', markersize=16.0, mec='DarkCyan', mew=3.0)
        #plt.plot(MWLowdist, np.log10(MWLow_h1frac), linestyle = 'None', marker='D', fillstyle='none', color = 'red', markersize=8.0, mew=2.0)
        #plt.plot(ALowdist, np.log10(ALow_h1frac), linestyle = 'None', marker='D', fillstyle='none', color = 'red', markersize=8.0, mew=2.0)
        #plt.plot(StarLowdist, np.log10(StarLow_h1frac), linestyle = 'None', marker='D', fillstyle='none', color = 'blue', markersize=8.0, mew=2.0)
        #plt.plot(Phoenixdist, np.log10(Phoenix_h1frac), linestyle = 'None', marker='o', fillstyle='none', color = 'blue', markersize=8.0, mew=2.0)

        plt.plot(np.log10(limit_mstar[0]*10**(6)), np.log10(limit_h1frac[0]*10**(6)), linestyle = 'None', marker='D', color = 'white', markersize=16.0, mec='DarkMagenta', mew=3.0)
        plt.plot(np.log10(limit_mstar[1:]*10**(6)), np.log10(limit_h1frac[1:]*10**(6)), linestyle = 'None', marker='D', color = 'white', markersize=16.0, mec='DarkMagenta', mew=3.0)

        plt.arrow(np.log10(21*10**(6)), np.log10(0.007*10**(6)), 0.0, -0.5, head_width=0.05, head_length=0.1, fc='DarkMagenta', ec='DarkMagenta')
        plt.arrow(np.log10(3.9*10**(6)), np.log10(0.35*10**(6)), 0.0, -0.5, head_width=0.05, head_length=0.1, fc='DarkMagenta', ec='DarkMagenta')
        plt.arrow(np.log10(1.1*10**(6)), np.log10(0.1*10**(6)), 0.0, -0.5, head_width=0.05, head_length=0.1, fc='DarkMagenta', ec='DarkMagenta')
        plt.arrow(np.log10(62*10**(6)), np.log10(0.005*10**(6)), 0.0, -0.5, head_width=0.05, head_length=0.1, fc='DarkMagenta', ec='DarkMagenta')
        plt.arrow(np.log10(7.6*10**(6)), np.log10(0.27*10**(6)), 0.0, -0.5, head_width=0.05, head_length=0.1, fc='DarkMagenta', ec='DarkMagenta')
        plt.arrow(np.log10(9.5*10**(6)), np.log10(0.36*10**(6)), 0.0, -0.5, head_width=0.05, head_length=0.1, fc='DarkMagenta', ec='DarkMagenta')
        plt.arrow(np.log10(2.8*10**(6)), np.log10(0.015*10**(6)), 0.0, -0.5, head_width=0.05, head_length=0.1, fc='DarkMagenta', ec='DarkMagenta')
        #plt.arrow(np.log10(ABmstar6), np.log10(fracHI*ABmstar6), 0.0, -0.5, head_width=0.05, head_length=0.1, fc='DarkMagenta', ec='DarkMagenta')

        plt.errorbar(np.log10(ABmstar6), np.log10(rps_blue6*ABmstar6), xerr = None, yerr=None, fmt='co', markersize = 10)
        #plt.errorbar(np.log10(ABmstar1), ABfrac1, xerr = None, yerr=LTerror9, fmt='co', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        p1 = plt.plot(np.log10(ABmstar6), np.log10(rps_blue6*ABmstar6), color = 'DarkCyan', marker = 'o', linestyle = 'None', markersize = 16, label=r"$\rm Pre-Stripping\ Sample$")
        #plt.plot(np.log10(ABmstar1), ABfrac1, color = 'Goldenrod', marker = 'o', linestyle = 'None', markersize = 16)

        mstarlist = np.linspace(5,9,50)
        plt.plot(mstarlist, mstarlist+np.log10(0.136), color = 'k', linestyle = '--', linewidth = 4)

        #p1b = plt.errorbar(np.log10(ABmstar6), np.log10(rps_red6), xerr = None, yerr=None, fmt='ro', markersize = 10)
        #plt.errorbar(np.log10(ABmstar1), ABfrac1, xerr = None, yerr=LTerror9, fmt='co', markersize = 5, elinewidth = 2.5, capsize = 10.0)
        #plt.plot(np.log10(ABmstar6), np.log10(rps_red6), color = 'DarkRed', marker = 'o', linestyle = 'None', markersize = 16)
        #plt.plot(np.log10(ABmstar1), ABfrac1, color = 'Goldenrod', marker = 'o', linestyle = 'None', markersize = 16)

        

        for i in range(len(ABmstar6)):
            if fracHI[i] < 0.01:
                #plotHI = fracHI[i]*ABmstar6[i]
                plotHI = 1.e3
                plt.arrow(np.log10(ABmstar6[i]), np.log10(plotHI), 0.0, -0.5, head_width=0.05, head_length=0.1, fc='DarkMagenta', ec='DarkMagenta')
                plt.errorbar(np.log10(ABmstar6[i]), np.log10(plotHI), xerr = None, yerr=None, fmt='mD', markersize = 10)
                p2 = plt.plot(np.log10(ABmstar6[i]), np.log10(plotHI), color = 'DarkMagenta', marker = 'D', linestyle = 'None', markersize = 16)
            else:
                plt.errorbar(np.log10(ABmstar6[i]), np.log10(fracHI[i]*ABmstar6[i]), xerr = None, yerr=None, fmt='mD', markersize = 10)
                p2 = plt.plot(np.log10(ABmstar6[i]), np.log10(fracHI[i]*ABmstar6[i]), color = 'DarkMagenta', marker = 'D', linestyle = 'None', markersize = 16)

        p2 = plt.plot(np.log10(ABmstar6[len(ABmstar6)-1]), np.log10(fracHI[i]*ABmstar6[i]), color = 'DarkMagenta', marker = 'D', linestyle = 'None', markersize = 16, label=r"$\rm Post-Stripping\ Sample$")
            
        
        xtickloc = [6,7,8]
        xtickstr = [r'$\rm 10^{6}$',r'$\rm 10^{7}$',r'$\rm 10^{8}$']
        #xtickstr = ['$'+str(kk)+'$' for kk in xtickloc]
        x2tickloc = [6,7,8]
        x2tickstr = [r'$\rm 10^{6}$',r'$\rm 10^{7}$',r'$\rm 10^{8}$']
        ytickloc = [2,3,4,5,6,7,8,9]
        #ytickstr = ['$'+str(kk)+'$' for kk in ytickloc]
        ytickstr = [r'$\rm 10^{2}$',r'$\rm 10^{3}$',r'$\rm 10^{4}$',r'$\rm 10^{5}$',r'$\rm 10^{6}$',r'$\rm 10^{7}$',r'$\rm 10^{8}$',r'$\rm 10^{9}$']
    
        # Make legend##################
        plt.legend(frameon = False, numpoints = 1, prop={'size':24}, loc = (0.05,0.8))
    
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


    if check == 'gasfrac':
        plt.savefig('RPS_AM_spekkenscomparison_halo30.pdf')

    else:
        plt.savefig('RPS_AM_spekkens_gasmass_halo30.pdf')

    #############################################################################################
    ## Second figure, cumulative histogram work
    #############################################################################################

    rps_blue6[np.where(rps_blue6 == 0.)[0]] = 0.01
    rps_red6[np.where(rps_red6 == 0.)[0]] = 0.01
    fracHI[np.where(fracHI <= 0.01)[0]] = 0.01

    plt.rc('axes',linewidth=axwidth)
    plt.figure(figsize=(13,13))

    plt.subplot2grid((2,2),(0,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([-2.1,1.1,0.0,1.05])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Gas\ Fraction$', fontsize = 28)
    ax.set_ylabel(r'$\rm Cumulative\ Gas\ Fraction$', fontsize = 28)

    plt.text(0.3,1.0, r'$\rm No\ Stripping$', fontsize = 22, color = 'b')
    plt.text(0.3,0.95, r'$\rm RPS\ Only$', fontsize = 22, color = 'r')
    plt.text(0.3,0.9, r'$\rm RPS\ and\ KH$', fontsize = 22, color = 'g')

    n, bins, patches = plt.hist(np.log10(rps_blue6), color = 'b', linewidth = 4.0, cumulative = -1, normed = 1, histtype = 'step')
    n, bins, patches = plt.hist(np.log10(rps_red6), color = 'r', linewidth = 4.0, cumulative = -1, normed = 1, histtype = 'step')
    n, bins, patches = plt.hist(np.log10(fracHI), color = 'g', linewidth = 4.0, cumulative = -1, normed = 1, histtype = 'step')
          
    xtickloc = [-2,-1,0,1]
    xtickstr = [r'$\rm 10^{-2}$',r'$\rm 10^{-1}$',r'$\rm 10^{0}$',r'$\rm 10^{1}$']
    ytickloc = [0.0,0.2,0.4,0.6,0.8,0.9,1.0]
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

    #plt.savefig('RPS_AM_gasfrac_cumulativehist.pdf')

    #############################################################################################
    ## Third figure, histogram work
    #############################################################################################

    plt.rc('axes',linewidth=axwidth)
    plt.figure(figsize=(13,13))

    plt.subplot2grid((2,2),(0,0), colspan = 2, rowspan = 2)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.4, hspace=0.4)
    plt.axis([-2.1,1.1,0,16])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Gas\ Fraction$', fontsize = 28)
    ax.set_ylabel(r'$\rm Number\ of\ Galaxies$', fontsize = 28)

    plt.text(0.3,15, r'$\rm No\ Stripping$', fontsize = 22, color = 'b')
    plt.text(0.3,14, r'$\rm RPS\ Only$', fontsize = 22, color = 'r')
    plt.text(0.3,13, r'$\rm RPS\ and\ KH$', fontsize = 22, color = 'g')

    n, bins, patches = plt.hist(np.log10(rps_blue6), color = 'b', linewidth = 4.0, histtype = 'step')
    n, bins, patches = plt.hist(np.log10(rps_red6), color = 'r', linewidth = 4.0, histtype = 'step')
    n, bins, patches = plt.hist(np.log10(fracHI), color = 'g', linewidth = 4.0, histtype = 'step')
          
    xtickloc = [-2,-1,0,1]
    xtickstr = [r'$\rm 10^{-2}$',r'$\rm 10^{-1}$',r'$\rm 10^{0}$',r'$\rm 10^{1}$']
    ytickloc = [0,5,10,15]
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

    #plt.savefig('RPS_AM_gasfrac_hist.pdf')

###################################################
# End of Plotting, time to display figure
###################################################

    plt.show()

########################################################################################
