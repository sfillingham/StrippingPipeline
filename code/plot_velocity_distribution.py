import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

def plot3panel(check='None'):
################################################################
##This function was originally used to plot the velocity distribution in 3 panels
##It has now been modified to adress point #2 in the referee report
##Checking the satellite mass dependence on Vsat and Rperi
################################################################

    data1 = Table.read('velocity_distribution_virial.dat', format = 'ascii')
    data2 = Table.read('velocity_distribution_halfvirial.dat', format = 'ascii')
    data3 = Table.read('velocity_distribution_pericenter.dat', format = 'ascii')

    masterdata = Table.read('/Users/Sean/research/elvis/ELVIS_Data_Clean/elvis_alldwarfs_clean_m59.dat', format = 'ascii')
    mhalo = np.array(masterdata['Mpeak(Msun)'])

    vel1a = np.array(data1['vel_host1'])
    vel1b = np.array(data1['vel_host2'])
    vel2a = np.array(data2['vel_host1'])
    vel2b = np.array(data2['vel_host2'])
    vel3a = np.array(data3['vel_host1'])
    vel3b = np.array(data3['vel_host2'])

    clean1a = vel1a[vel1a > 0.0]
    clean1b = vel1b[vel1b > 0.0]
    clean2a = vel2a[vel2a > 0.0]
    clean2b = vel2b[vel2b > 0.0]
    clean3a = vel3a[vel3a > 0.0]
    clean3b = vel3b[vel3b > 0.0]
    halo1a = mhalo[vel1a > 0.0]
    halo1b = mhalo[vel1b > 0.0]
    halo2a = mhalo[vel2a > 0.0]
    halo2b = mhalo[vel2b > 0.0]
    halo3a = mhalo[vel3a > 0.0]
    halo3b = mhalo[vel3b > 0.0]

    mvel1 = np.append(clean1a, clean1b)
    mvel2 = np.append(clean2a, clean2b)
    mvel3 = np.append(clean3a, clean3b)
    mhalo1 = np.append(halo1a, halo1b)
    mhalo2 = np.append(halo2a, halo2b)
    mhalo3 = np.append(halo3a, halo3b)

    print 'mean1 = '+np.str(np.mean(mvel1))
    print 'mean2 = '+np.str(np.mean(mvel2))
    print 'mean3 = '+np.str(np.mean(mvel3))
    print 'median1 = '+np.str(np.median(mvel1))
    print 'median2 = '+np.str(np.median(mvel2))
    print 'median3 = '+np.str(np.median(mvel3))
    print 'stddev1 = '+np.str(np.std(mvel1))
    print 'stddev2 = '+np.str(np.std(mvel2))
    print 'stddev3 = '+np.str(np.std(mvel3))

    ##############################################
    #bin the data by Mstar for alternative figures
    ##############################################
    binsize = 0.6
    iterations = 13
    #binlist = np.array([6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0])
    binlist = np.linspace(8.0, 11.5, 50)
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

    for i in range(len(binlist)):
        cut1 = (np.log10(mhalo1) > (binlist[i]-binsize)) & (np.log10(mhalo1) < (binlist[i]+binsize))
        cut2 = (np.log10(mhalo2) > (binlist[i]-binsize)) & (np.log10(mhalo2) < (binlist[i]+binsize))
        cut3 = (np.log10(mhalo3) > (binlist[i]-binsize)) & (np.log10(mhalo3) < (binlist[i]+binsize))

        mean1[i] = np.mean(mvel1[cut1])
        upper1[i] = np.mean(mvel1[cut1]) + np.std(mvel1[cut1])
        lower1[i] = np.mean(mvel1[cut1]) - np.std(mvel1[cut1])
        mean2[i] = np.mean(mvel2[cut2])
        upper2[i] = np.mean(mvel2[cut2]) + np.std(mvel2[cut2])
        lower2[i] = np.mean(mvel2[cut2]) - np.std(mvel2[cut2])
        mean3[i] = np.mean(mvel3[cut3])
        upper3[i] = np.mean(mvel3[cut3]) + np.std(mvel3[cut3])
        lower3[i] = np.mean(mvel3[cut3]) - np.std(mvel3[cut3])


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
    plt.subplot2grid((1,3),(0,0), colspan = 1, rowspan = 1)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.2, hspace=0.1)
    #plt.axis([0,600,0,170])
    plt.axis([8.8,12,0,600])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Halo\ Mass\ (M_{\odot})$', fontsize = 28)
    ax.set_ylabel(r'$\rm Velocity\ Relative\ to\ Host\ (km/s)$', fontsize = 28)

    #plt.text(250,155,r'$V_{\rm infall}$', fontsize = 20)
    #plt.text(250,148,r'$\rm mean\ =\ 184\ km/s$', fontsize = 20)
    #plt.text(250,140,r'$\rm \sigma\ =\ 65\ km/s$', fontsize = 20)

    #n, bins, patches = plt.hist(mvel1, 50, histtype='step', color = 'DarkCyan', linewidth = 4.0)
    if check == 'scatter':
        plt.plot(np.log10(mhalo1), mvel1, color = 'c', marker = 'o', markersize = 14, linestyle = 'None')
    else:
        plt.fill_between(binlist, upper1, lower1, facecolor = 'c', alpha = 0.4)
        plt.plot(binlist, mean1, color = 'DarkCyan', marker = 'None', linestyle = '-', linewidth = 3.0)

    #xtickloc = [0,100,200,300,400,500,600]
    xtickloc = [9,10,11,12]
    xtickstr = ['$'+str(kk)+'$' for kk in xtickloc]
    #ytickloc = [0, 50, 100, 150]
    ytickloc = [0,100,200,300,400,500,600]
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
    plt.subplot2grid((1,3),(0,1), colspan = 1, rowspan = 1)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.2, hspace=0.1)
    #plt.axis([0,600,0,170])
    plt.axis([8.8,12,0,600])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Halo\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Velocity\ Relative\ to\ Host\ (km/s)$', fontsize = 28)

    #plt.text(250,155,r'$V_{\rm 0.5R_{vir}}$', fontsize = 20)
    #plt.text(250,148,r'$\rm mean\ =\ 240\ km/s$', fontsize = 20)
    #plt.text(250,140,r'$\rm \sigma\ =\ 65\ km/s$', fontsize = 20)

    #n, bins, patches = plt.hist(mvel2, 50, histtype='step', color = 'DarkMagenta', linewidth = 4.0)
    if check == 'scatter':
        plt.plot(np.log10(mhalo2), mvel2, color = 'm', marker = 'o', markersize = 14, linestyle = 'None')
    else:
        plt.fill_between(binlist, upper2, lower2, facecolor = 'm', alpha = 0.4)
        plt.plot(binlist, mean2, color = 'DarkMagenta', marker = 'None', linestyle = '-', linewidth = 3.0)

    #xtickloc = [0,100,200,300,400,500,600]
    xtickloc = [9,10,11,12]
    xtickstr = ['$'+str(kk)+'$' for kk in xtickloc]
    #ytickloc = [0, 50, 100, 150]
    ytickloc = [0,100,200,300,400,500,600]
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
    plt.subplot2grid((1,3),(0,2), colspan = 1, rowspan = 1)
    plt.subplots_adjust(left=0.1, bottom=0.13, right=0.97, top=0.97, wspace=0.2, hspace=0.1)
    #plt.axis([0,600,0,170])
    plt.axis([8.8,12,0,600])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm Halo\ Mass\ (M_{\odot})$', fontsize = 28)
    #ax.set_ylabel(r'$\rm Velocity\ Relative\ to\ Host\ (km/s)$', fontsize = 28)

    #plt.text(250,155,r'$V_{\rm pericenter}$', fontsize = 20)
    #plt.text(250,148,r'$\rm mean\ =\ 318\ km/s$', fontsize = 20)
    #plt.text(250,140,r'$\rm \sigma\ =\ 74\ km/s$', fontsize = 20)

    #n, bins, patches = plt.hist(mvel3, 50, histtype='step', color = 'Goldenrod', linewidth = 4.0)
    if check == 'scatter':
        plt.plot(np.log10(mhalo3), mvel3, color = 'Goldenrod', marker = 'o', markersize = 14, linestyle = 'None')
    else:
        plt.fill_between(binlist, upper3, lower3, facecolor = 'y', alpha = 0.4)
        plt.plot(binlist, mean3, color = 'Goldenrod', marker = 'None', linestyle = '-', linewidth = 3.0)

    #xtickloc = [0,100,200,300,400,500,600]
    xtickloc = [9,10,11,12]
    xtickstr = ['$'+str(kk)+'$' for kk in xtickloc]
    #ytickloc = [0,50,100,150]
    ytickloc = [0,100,200,300,400,500,600]
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

    if check == 'scatter':
        plt.show()

    else:
        plt.savefig('refreport_infallvelocity.pdf')
        plt.show()

def plot1panel():

    data1 = Table.read('velocity_distribution_virial.dat', format = 'ascii')
    data2 = Table.read('velocity_distribution_halfvirial.dat', format = 'ascii')
    data3 = Table.read('velocity_distribution_pericenter.dat', format = 'ascii')

    vel1a = np.array(data1['vel_host1'])
    vel1b = np.array(data1['vel_host2'])
    vel2a = np.array(data2['vel_host1'])
    vel2b = np.array(data2['vel_host2'])
    vel3a = np.array(data3['vel_host1'])
    vel3b = np.array(data3['vel_host2'])

    clean1a = vel1a[vel1a > 0.0]
    clean1b = vel1b[vel1b > 0.0]
    clean2a = vel2a[vel2a > 0.0]
    clean2b = vel2b[vel2b > 0.0]
    clean3a = vel3a[vel3a > 0.0]
    clean3b = vel3b[vel3b > 0.0]

    mvel1 = np.append(clean1a, clean1b)
    mvel2 = np.append(clean2a, clean2b)
    mvel3 = np.append(clean3a, clean3b)

    cut300 = mvel3[mvel3 > 300.0]
    cut250 = mvel3[mvel3 > 250.0]
    cut200 = mvel3[mvel3 > 200.0]

    print 'mean1 = '+np.str(np.mean(mvel1))
    print 'mean2 = '+np.str(np.mean(mvel2))
    print 'mean3 = '+np.str(np.mean(mvel3))
    print 'median1 = '+np.str(np.median(mvel1))
    print 'median2 = '+np.str(np.median(mvel2))
    print 'median3 = '+np.str(np.median(mvel3))
    print 'stddev1 = '+np.str(np.std(mvel1))
    print 'stddev2 = '+np.str(np.std(mvel2))
    print 'stddev3 = '+np.str(np.std(mvel3))

    print np.float(len(cut300)) / len(mvel3)
    print np.float(len(cut250)) / len(mvel3)
    print np.float(len(cut200)) / len(mvel3)

####################################
#Plot work
####################################

    axwidth = 3
    axlength = 10
    fontsize=30
    plt.rc('axes',linewidth=axwidth)
    plt.figure(figsize=(11,11))

####################
#FIGURE 1a
####################
    plt.subplot2grid((1,1),(0,0), colspan = 1, rowspan = 1)
    plt.subplots_adjust(left=0.12, bottom=0.13, right=0.97, top=0.97, wspace=0.2, hspace=0.1)
    #plt.axis([0,600,0,170])
    plt.axis([0,600,0,1.2])
    ax = plt.gca()
    ax.set_ylabel(r'$\rm Number\ of\ Subhalos$', fontsize = 30)
    ax.set_xlabel(r'$\rm Velocity\ Relative\ to\ Host\ (km/s)$', fontsize = 30)

    #plt.text(50,195,r'$V_{\rm infall}$', fontsize = 20, color = 'DarkCyan')
    #plt.text(50,188,r'$\rm mean\ =\ 184\ km/s$', fontsize = 20)
    #plt.text(50,180,r'$\rm \sigma\ =\ 65\ km/s$', fontsize = 20)

    #plt.text(150,195,r'$V_{\rm 0.5R_{vir}}$', fontsize = 20, color = 'DarkMagenta')
    #plt.text(150,188,r'$\rm mean\ =\ 240\ km/s$', fontsize = 20)
    #plt.text(150,180,r'$\rm \sigma\ =\ 65\ km/s$', fontsize = 20)

    #plt.text(350,195,r'$V_{\rm pericenter}$', fontsize = 20, color = 'Goldenrod')
    #plt.text(350,188,r'$\rm mean\ =\ 318\ km/s$', fontsize = 20)
    #plt.text(350,180,r'$\rm \sigma\ =\ 74\ km/s$', fontsize = 20)

    n, bins, patches = plt.hist(mvel1, 50, histtype='step', normed = 'True', cumulative = 'True', color = 'DarkCyan', linewidth = 4.0, label = r'$\rm velocity\ at\ infall$')
    n, bins, patches = plt.hist(mvel2, 50, histtype='step', normed = 'True', cumulative = 'True', color = 'DarkMagenta', linewidth = 4.0, label = r'$\rm velocity\ at\ 0.5{\it R}_{vir}$')
    n, bins, patches = plt.hist(mvel3, 50, histtype='step', normed = 'True', cumulative = 'True', color = 'Goldenrod', linewidth = 4.0, label = r'$\rm velocity\ at\ pericenter$')

    #n, bins, patches = plt.hist(mvel1, 50, histtype='step', color = 'DarkCyan', linewidth = 4.0, label = r'$\rm velocity\ at\ infall$')
    #n, bins, patches = plt.hist(mvel2, 50, histtype='step', color = 'DarkMagenta', linewidth = 4.0, hatch = '/', label = r'$\rm velocity\ at\ 0.5{\it R}_{vir}$')
    #n, bins, patches = plt.hist(mvel3, 50, histtype='step', color = 'Goldenrod', linewidth = 4.0, hatch = '|',label = r'$\rm velocity\ at\ pericenter$')

    xtickloc = [0,100,200,300,400,500,600]
    xtickstr = ['$'+str(kk)+'$' for kk in xtickloc]
    #ytickloc = [0, 50, 100, 150]
    ytickloc = [0, 0.5, 1.0]
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

    plt.legend(frameon = False, numpoints = 1, prop={'size':24}, loc = 2)

    #plt.savefig('infallvelocity_distribution_1panel.pdf')
    #plt.savefig('infallvelocity_distribution_1panel.png')

    plt.savefig('infallvelocity_cdistribution_1panel.pdf')
    plt.savefig('infallvelocity_cdistribution_1panel.png')

    plt.show()
