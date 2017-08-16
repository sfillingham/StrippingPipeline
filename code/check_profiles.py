import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
import massprofile as mass

def profiles():
    
    dist = np.arange(0,10,0.1)
    ab_match = Table.read('GK14AM.txt', format = 'ascii')
    abMpeak = np.array(ab_match['Mpeak(Msun)'])
    abMstar = np.array(ab_match['Mstar(Msun)'])
    
    #A_Pace data
    pdata = Table.read('best_fit_parameters_mlfixed.dat', format = 'ascii')
    pname = np.array(pdata['NAME'])
    r0 = 10**np.array(pdata['r0'])
    rho0 = 10**np.array(pdata['rho0'])
    rs = 10**np.array(pdata['rs'])
    rhos = 10**np.array(pdata['rhos'])
    
    #LT data
    LTdata = Table.read('LT_galaxylist.txt', format = 'ascii')
    LTname = np.array(LTdata['NAME'])
    vmag = np.array(LTdata['Mv'])
    massinput = 10**np.array(LTdata['virmass'])
    LTconcentration = np.array(LTdata['c'])
    
    #AB data
    ABdata = Table.read('AB_galaxylist.txt', format = 'ascii')
    ABname = np.array(ABdata['NAME'])
    ABvmag = np.array(ABdata['Mv'])
    
    ABmstar = 10**((4.83-(ABvmag))/2.5)
    ABmassinput = np.interp(ABmstar, abMstar, abMpeak)
    ABconcentration = 9.6*(ABmassinput/1.429e12)**(-0.075)

    for i in range(len(pdata)):

        LTmask = LTname == pname[i]
        clean_LTname = LTname[LTmask]
        clean_LTmass = massinput[LTmask]
        clean_LTc = LTconcentration[LTmask]
        
        ABmask = ABname == pname[i]
        clean_ABname = ABname[ABmask]
        clean_ABmass = ABmassinput[ABmask]
        clean_ABc = ABconcentration[ABmask]

        burkertmass = mass.burkert(r0[i], rho0[i], dist)
        
        pNFWmass = mass.NFW_pace(rs[i], rhos[i], dist)
        
        NFWmassLT = mass.NFW(clean_LTmass, clean_LTc, dist)
        
        NFWmassAB = mass.NFW(clean_ABmass, clean_ABc, dist)

        plt.figure(i)
        plt.plot(dist, np.log10(burkertmass), 'b')
        plt.plot(dist, np.log10(NFWmassLT), 'r')
        plt.plot(dist, np.log10(pNFWmass), 'g')
        plt.plot(dist, np.log10(NFWmassAB), 'DarkOrange')
        plt.text(6.0,6.4,r'$\rm Burkert$', color = 'b', fontsize = 20)
        plt.text(6.0,6.1,r'$\rm Pace\ NFW$', color = 'g', fontsize = 20)
        plt.text(6.0,5.8,r'$\rm LT\ NFW$', color = 'r', fontsize = 20)
        plt.text(6.0,5.5,r'$\rm AB\ NFW$', color = 'DarkOrange', fontsize = 20)

        plt.ylabel(r'$\rm log\left(\frac{M_{*}}{M_{\odot}}\right)$', fontsize = 20)
        plt.xlabel(r'$\rm radial\ distance\ (kpc)$', fontsize = 20)

        plt.savefig('profilecheck/'+pname[i]+'_profilecheck.pdf')
        plt.show()



def concentration():
    
    ab_match = Table.read('GK14AM.txt', format = 'ascii')
    abMpeak = np.array(ab_match['Mpeak(Msun)'])
    abMstar = np.array(ab_match['Mstar(Msun)'])

    #LT data
    LTdata = Table.read('LT_galaxylist.txt', format = 'ascii')
    LTname = np.array(LTdata['NAME'])
    LTvmag = np.array(LTdata['Mv'])
    LTmassinput = 10**np.array(LTdata['virmass'])
    LTconcentration = np.array(LTdata['c'])

    #AB data
    ABdata = Table.read('AB_galaxylist.txt', format = 'ascii')
    ABname = np.array(ABdata['NAME'])
    ABvmag = np.array(ABdata['Mv'])
            
    ABmstar = 10**((4.83-(ABvmag))/2.5)
    ABmassinput = np.interp(ABmstar, abMstar, abMpeak)
    ABconcentration = 9.6*(ABmassinput/1.429e12)**(-0.075)
    
    c_final = np.empty(len(LTname))

    for i in range(len(LTname)):

        mask = ABname == LTname[i]
        clean_name = ABname[mask]
        clean_c = ABconcentration[mask]

        c_final[i] = clean_c

    print c_final
    print LTconcentration
    
    plt.figure(1)
    plt.axis([5,15,5,15])
    plt.plot(LTconcentration, c_final, 'bo')
    plt.plot(np.arange(20), np.arange(20), 'k--')
    plt.ylabel(r'$\rm Klypin\ Concentration$', fontsize = 20)
    plt.xlabel(r'$\rm Little\ THINGS\ Concentration$', fontsize = 20)
        
    plt.savefig('profilecheck/concentrationcheck.pdf')
    plt.show()





