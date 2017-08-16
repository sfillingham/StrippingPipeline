import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import code.test_stripping_individual as tsi
import astropy.units as u
import astropy.constants as consts
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline as IUSpline


## This routine will determine the amount of gas removed from a sample of galaxy HI profiles iteratively.
## The total time, and number of timesteps will need to be specified

def HIloss_trapz(dwarf, Rlimit):

    testdata = Table.read('THINGS_galaxylist.dat', format = 'ascii')
    testname = np.array(testdata['NAME'])

    data = Table.read('data/'+dwarf+'_clean.dat', format = 'ascii')
    radius = np.array(data['R(kpc)'])
    gas = np.array(data['N_HI'])

    check = dwarf == testname
    test = testname[check]

    #This accounts for the 1.36 multiplicative factor based on He abundance
    if len(test) == 0:
        gas = gas + 0.13354
    else:
        gas = gas

    clean_radius = radius[gas < 50]
    clean_gas = gas[gas < 50]

    radiusfinal = np.linspace(min(clean_radius), max(clean_radius), len(clean_radius)*10)
    gasinterp = IUSpline(clean_radius, clean_gas, k=3)
    gasfinal = gasinterp(radiusfinal)
    

    func = 2*np.pi*(radiusfinal*((3.08E21)**2))*(10**gasfinal)
    mass = np.trapz(func,radiusfinal)

    striplimit = Rlimit #np.max([Rlimit,np.min(radiusfinal[gasfinal < gaslimit])])
    stripradius = radiusfinal[radiusfinal < striplimit]
    stripfunc = func[radiusfinal < striplimit]
    nonstripmass = np.trapz(stripfunc, stripradius)

    frac = (mass - nonstripmass) / mass

    gasmass = nonstripmass*(consts.m_p / consts.M_sun)
    
    return frac, gasmass, mass*(consts.m_p / consts.M_sun)


def KH(totaltime, timestep, hotdensity, velinfall):

    ab_match = Table.read('GK14AM.txt', format = 'ascii')
    abMpeak = np.array(ab_match['Mpeak(Msun)'])
    abMstar = np.array(ab_match['Mstar(Msun)'])

    check = 'AB'

    print 'GK14 read in'

    data = Table.read('AB_galaxylist.txt', format = 'ascii')
    name = np.array(data['NAME'])
    vmag = np.array(data['Mv'])

    print 'AB data read in'

    tdata = Table.read('THINGS_galaxylist.dat', format = 'ascii')
    tname = np.array(tdata['NAME'])
    tmstar = np.array(tdata['Mstar'])

    print 'THINGS read in'

    sdata = Table.read('SHIELD_list.dat', format = 'ascii')
    sname = np.array(sdata['NAME'])
    svmag = np.array(sdata['Mv'])
    smstar = 10**((4.83-(svmag))/2.5)
    #smstar = np.array(10**sdata['Mstar'])

    print 'SHIELD read in'
            
    tcut = tmstar > 0.0
    tmstar_clean = 10**tmstar[tcut]
    tname_clean = tname[tcut]
                        
    mstar = 10**((4.83-(vmag))/2.5)
    massinput = np.interp(mstar, abMstar, abMpeak)
    concentration = 9.6*(massinput/1.429e12)**(-0.075)

    tmassinput = np.interp(tmstar_clean, abMstar, abMpeak)
    tconcentration = 9.6*(tmassinput/1.429e12)**(-0.075)

    smassinput = np.interp(smstar, abMstar, abMpeak)
    sconcentration = 9.6*(smassinput/1.429e12)**(-0.075)

    tmaster_name = np.append(name, tname_clean)
    tmaster_mass = np.append(massinput, tmassinput)
    tmaster_con = np.append(concentration, tconcentration)

    master_name = np.append(tmaster_name, sname)
    master_mass = np.append(tmaster_mass, smassinput)
    master_con = np.append(tmaster_con, sconcentration)

    print 'final data collected'

    r0 = 1.0
    rho0 = 1.0
    rs = 1.0
    rhos = 1.0

    fracNFW = np.empty(len(master_name))
    gasNFW = np.empty(len(master_name))
    gmassNFW = np.empty(len(master_name))
    upperfracNFW = np.empty(len(master_name))
    lowerfracNFW = np.empty(len(master_name))
    stripRNFW = np.empty(len(master_name))
    visc1 = np.empty(len(master_name))
    DMmassO = np.empty(len(master_name))
        
    Mstar1 = np.append(mstar, tmstar_clean)
    Mstar = np.append(Mstar1, smstar)

    time, dt = np.linspace(0, totaltime, num=timestep, retstep = True)

    for i in range(1):#len(master_name)):

        #print master_name[i]
        stripNFW, DMmass, gasmass, radius, DarkMass = tsi.strip(master_name[i], master_mass[i], master_con[i], r0, rho0, rs, rhos, hotdensity, velinfall, check)

        totalmass = gasmass/(consts.M_sun/consts.m_p)*(3.086*10**21)**2

        DMmassO[i] = DMmass
        #print np.log10(DMmassO[i])
        #print np.log10(DarkMass[-1])
        
        stripRNFW[i] = stripNFW
        print stripNFW
            
        fracNFW[i],gasNFW[i],gmassNFW[i] = HIloss_trapz(master_name[i], stripNFW)

        #print fracNFW

        #Test whether KH instability will strip
        strip_pressure = (10**(-1*hotdensity))*(velinfall)**2 #units: cm**(-3) (km/s)**2
        print strip_pressure

        Ddepth = stripNFW #This might need to be changed to something less obvious
        h1_avg = gasNFW[i]/(np.pi*(stripNFW**2)*Ddepth)
        KH_restore = (((h1_avg)*consts.G*DMmass[0])/(2*np.pi*stripNFW))*((1e-4*consts.M_sun)*(u.s**2/u.m**3)*(consts.M_sun/consts.m_p)/(3.086*10**21)**4)
        print 'KH1 ='
        print KH_restore

        tgasmass = totalmass
        gasradius = radius
        gasrem = gasNFW[i]
        dmass = DarkMass

        for j in range(len(time)):

            print 'time = '
            print time[j]
            
            if (KH_restore < strip_pressure):
                new_visc1 = (2e10)*((stripNFW/20.)**2)*((10**(-1*hotdensity))/(10**(-3)))*(velinfall/1000.0)*dt
            else:
                new_visc1 = 0.0

            visc1[i] = visc1[i] + new_visc1

            print new_visc1

            # Determine the new strip radius based on the amount of ISM removed via KH
            mass_check = np.trapz(2*np.pi*gasradius*tgasmass,gasradius) - new_visc1
            m_check = np.cumsum(tgasmass)
            print len(gasradius[np.where(m_check < mass_check)[0]])
            stripNFW = gasradius[np.where(m_check < mass_check)[0]][-1]

            ## Reset the input values and total gas remaining based on the new stripping radius
            gasradius = gasradius[np.where(m_check < mass_check)[0]]
            tgasmass = tgasmass[np.where(m_check < mass_check)[0]]
            dmass = dmass[np.where(m_check < mass_check)[0]]
            gasrem = np.trapz(2*np.pi*gasradius*tgasmass, gasradius)
            print gasrem
            Ddepth = stripNFW 
            h1_avg = gasrem/(np.pi*(stripNFW**2)*Ddepth)
            KH_restore = (((h1_avg)*consts.G*dmass[-1])/(2*np.pi*stripNFW))*((1e-4*consts.M_sun)*(u.s**2/u.m**3)*(consts.M_sun/consts.m_p)/(3.086*10**21)**4)
            print KH_restore
            
            #print visc1[0]
            print gasradius[-1]
            
        #print 'visc1 = '
        #print visc1[0]


    ratio_red = gasNFW / Mstar
    ratio_blue = gmassNFW / Mstar

    KH1 = np.empty(len(Mstar))
    for i in range(len(Mstar)):
        if gasNFW[i] == 0:
            KH1[i] = 0.0
        else:
            KH1[i] = visc1[i] / gasNFW[i]
               
            
    mastertable = Table([master_name, Mstar, fracNFW, ratio_red, ratio_blue, KH1], names=('NAME', 'Mstar', 'NFW', 'Rgasfrac', 'Bgasfrac', 'KH_1Gyr'), meta={'strippedfrac': 'RPS'})
            
    masteroutputfile = 'stripdata/RPS_strippedfrac_hot'+np.str(np.int(10*hotdensity))+'_v'+np.str(np.int(velinfall))+'_'+check+'.dat'

    print 'gasmass'
    return gasmass, radius

    #mastertable.write(masteroutputfile, format = 'ascii')

    #return Mstar, fracNFW, gasNFW, gmassNFW, DMmassO, visc1
