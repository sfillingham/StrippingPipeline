import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as cb
import matplotlib.cm as cmx
from astropy.io import ascii
from astropy.table import Table

def profilefig(check, colorcode='cool'):

    LTdata = Table.read('LT_galaxylist.txt', format = 'ascii')
    Tdata = Table.read('THINGS_galaxylist.dat', format = 'ascii')
    Sdata = Table.read('SHIELD_list.dat', format = 'ascii')

    modeldata = Table.read('model_input.dat', format = 'ascii')
    modelname = np.array(modeldata['NAME'])
    modelmstar = np.array(modeldata['Mstar'])

    LTname = np.array(LTdata['NAME'])
    Tname = np.array(Tdata['NAME'])
    Sname = np.array(Sdata['NAME'])
    num = len(LTname)
    num2 = len(Tname) + num
    
    Mv = np.array(LTdata['Mv'])
    LTmstar = 10**((4.83-(Mv))/2.5)
    Tmstar = np.array(Tdata['Mstar'])
    Smstar = np.array(Sdata['Mstar'])
    clean = Tmstar != -1

    Tname_clean = Tname[clean]
    num2 = len(Tname_clean) + num
    Tmstar_clean = Tmstar[clean]
    
    master_Mstar1 = np.append(np.log10(LTmstar), Tmstar_clean)
    master_Mstar = np.append(master_Mstar1, Smstar)
    print len(master_Mstar)
    master_name1 = np.append(LTname, Tname_clean)
    master_name2 = np.append(master_name1, Sname)
    print master_name2
    print len(master_name2)

    if check == 'model':
        cool = plt.get_cmap(colorcode)
        cNorm = colors.Normalize(np.min(modelmstar), np.max(modelmstar))
        sm = plt.cm.ScalarMappable(norm = cNorm, cmap = cool)
        
    else:
        cool = plt.get_cmap(colorcode)
        cNorm = colors.Normalize(np.min(master_Mstar), np.max(master_Mstar))
        sm = plt.cm.ScalarMappable(norm = cNorm, cmap = cool)

####################################
#Plot work
####################################

    axwidth = 3
    axlength = 10
    fontsize=28
    
    plt.rc('axes',linewidth=axwidth)
    
    fig = plt.figure(figsize=(8,9))
    plt.subplot2grid((1,1),(0,0), colspan = 1, rowspan = 1)
    plt.subplots_adjust(left=0.15, bottom=0.13, right=0.97, top=0.83, wspace=0.6, hspace=0.1)
    #plt.axis([0.1,1.4,15,22])
    plt.axis([0,24,15,22])
    ax = plt.gca()
    ax.set_xlabel(r'$\rm {\it R}\ (kpc)$', fontsize = 28)
    ax.set_ylabel(r'$\rm log_{10}\ {\it N}_{HI}\ (cm^{-2})$', fontsize = 28)

    if check == 'model':
        
        for i in range(len(modelname)):

            dwarfdata = Table.read('data/'+modelname[i]+'_clean.dat', format = 'ascii')
            gas = np.array(dwarfdata['N_HI'])
            r = np.array(dwarfdata['R(kpc)'])

            #gas_clean = gas[gas < 50] + 0.13354
            gas_clean = gas[gas < 50] 
            r_clean = r[gas < 50]

            colorval = sm.to_rgba(modelmstar[i])

            #plt.plot(np.log10(r_clean), gas_clean, color = colorval, linestyle = '-', linewidth = 2.0)
            plt.plot(r_clean, gas_clean, color = colorval, linestyle = '-', linewidth = 2.0)
        
    else:

        for i in range(len(LTname)):

            dwarfdata = Table.read('data/'+LTname[i]+'_clean.dat', format = 'ascii')
            gas = np.array(dwarfdata['N_HI'])
            r = np.array(dwarfdata['R(kpc)'])

            gas_clean = gas[gas < 50] + 0.13354
            r_clean = r[gas < 50]

            colorval = sm.to_rgba(master_Mstar[i])

            #plt.plot(np.log10(r_clean), gas_clean, color = colorval, linestyle = '-', linewidth = 2.0)
            plt.plot(r_clean, gas_clean, color = colorval, linestyle = '-', linewidth = 2.0)

        for i in range(len(Tname_clean)):

            dwarfdata = Table.read('data/'+Tname_clean[i]+'_clean.dat', format = 'ascii')
            gas = np.array(dwarfdata['N_HI'])
            r = np.array(dwarfdata['R(kpc)'])

            gas_clean = gas[gas < 50] 
            r_clean = r[gas < 50]

            colorval = sm.to_rgba(master_Mstar[num+i])

            #plt.plot(np.log10(r_clean), gas_clean, color = colorval, linestyle = '-', linewidth = 2.0)
            plt.plot(r_clean, gas_clean, color = colorval, linestyle = '-', linewidth = 2.0)

        for i in range(len(Sname)):

            dwarfdata = Table.read('data/'+Sname[i]+'_clean.dat', format = 'ascii')
            gas = np.array(dwarfdata['N_HI'])
            r = np.array(dwarfdata['R(kpc)'])

            gas_clean = gas[gas < 50] + 0.13354
            r_clean = r[gas < 50]

            colorval = sm.to_rgba(master_Mstar[num2+i])

            #plt.plot(np.log10(r_clean), gas_clean, color = colorval, linestyle = '-', linewidth = 2.0)
            plt.plot(r_clean, gas_clean, color = colorval, linestyle = '-', linewidth = 2.0)

        for i in range(len(Sname)):

            dwarfdata = Table.read('data/'+Sname[i]+'_clean_kpc.dat', format = 'ascii')
            gas = np.array(dwarfdata['N_HI'])
            r = np.array(dwarfdata['R(kpc)'])

            gas_clean = gas[gas < 50] + 0.13354
            r_clean = r[gas < 50]

            colorval = sm.to_rgba(master_Mstar[num2+i])

            #plt.plot(r_clean, gas_clean, color = 'Red', linestyle = '-', linewidth = 2.0)

    
    xtickloc = [0,3,6,9,12,15,18,21,24]
    #xtickloc = [0.2,0.5,1,5,10,20]
    xtickstr = ['$'+str(kk)+'$' for kk in xtickloc]
    ytickloc = [16,17,18,19,20,21,22]
    ytickstr = ['$'+str(kk)+'$' for kk in ytickloc]
    
    ax.set_xticks(xtickloc)
    #ax.set_xticks(np.log10(xtickloc))
    ax.set_xticklabels(xtickstr, position = (0,-0.01))
    ax.set_yticks(ytickloc)
    ax.set_yticklabels(ytickstr)
    #ax.set_xscale('log')
    
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
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(2)
    
    ax.tick_params(which='major',width=axwidth,length=axlength+5)
    ax.tick_params(which='minor',width=axwidth,length=axlength)
    
    sm._A = []
    cbaxes = fig.add_axes([0.15, 0.86, 0.82, 0.04]) 
    #cb = plt.colorbar(ax1, cax = cbaxes)  
    cbar = plt.colorbar(sm, cax = cbaxes, orientation = 'horizontal')
    cbartickloc = [6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5]
    cbartickstr = ['$'+str(kk)+'$' for kk in cbartickloc]
    cbar.set_ticks(cbartickloc)
    cbar.set_ticklabels(cbartickstr)
    cbar.ax.tick_params(labelsize=24, pad=-51, width = 2, size = 5, top = 'on')
    #cbaxes.set_xlabel(r'$\rm log_{10}\ {\it M}_{\star}\ (M_{\odot})$', fontsize=24, labelpad=-58)
    #cbar.set_label(r'$\rm log\ {\it M}_{\star}\ (M_{\odot})$', fontsize = 24)

    #if check == 'model':
        #plt.savefig('HI_profiles_model.pdf')
        #plt.savefig('HI_profiles_model.png')
    #else:
        #plt.savefig('HI_profiles.pdf')
        #plt.savefig('HI_profiles.png')
        #plt.savefig('HI_profiles_loglog.pdf')
        #plt.savefig('HI_profiles_loglog.png')
        
    #plt.show()
    
