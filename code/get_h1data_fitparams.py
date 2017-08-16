import numpy as np
#import matplotlib.pyplot as plt
from astropy.table import Table
import code.test_stripping_individual as tsi
import code.get_HImassloss_trapz as get
import code.massprofile as massprofile
import astropy.units as u
import astropy.constants as consts
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline as IUSpline



def frac_vs_mstar(hotdensity, velinfall, inputfile, check):
    
    ab_match = Table.read('GK14AM.txt', format = 'ascii')
    abMpeak = np.array(ab_match['Mpeak(Msun)'])
    abMstar = np.array(ab_match['Mstar(Msun)'])

    if check == 'Apace':
        
        data = Table.read(inputfile, format = 'ascii')
        name = np.array(data['NAME'])
        r0 = 10**np.array(data['r0'])
        rho0 = 10**np.array(data['rho0'])
        rs = 10**np.array(data['rs'])
        rhos = 10**np.array(data['rhos'])

        stardata = Table.read('AB_galaxylist.txt', format = 'ascii')
        tstardata = Table.read('THINGS_galaxylist.dat', format = 'ascii')
        starname = np.array(stardata['NAME'])
        tstarname = np.array(tstardata['NAME'])
        master_name = np.append(starname, tstarname)
        
        vmag = np.array(stardata['Mv'])
        mstar = 10**((4.83-(vmag))/2.5)
        tmstar = 10**np.array(tstardata['Mstar'])
        master_mstar = np.append(mstar, tmstar)
        
        massinput = 1.0
        concentration = 1.0
        
        ########################################################
        # define the empty arrays which will be filled with data
        ########################################################
        
        fracNFW = np.empty(len(name))
        gasNFW = np.empty(len(name))
        gmassNFW = np.empty(len(name))
        fracBurk = np.empty(len(name))
        gasBurk = np.empty(len(name))
        gmassBurk = np.empty(len(name))
        stripRNFW = np.empty(len(name))
        stripRBurk = np.empty(len(name))
        Mstar = np.empty(len(name))

        visc1_NFW = np.empty(len(name))
        visc1_Burk = np.empty(len(name))
        
        ############################################################################
        #loop through each list of dwarf names and determine the relevant properties
        ############################################################################
        
        for i in range(len(name)):
            
            cut = master_name == name[i]
            Mstar[i] = master_mstar[cut]
        
            stripNFW, stripBurk, NFWmass, Burkmass = tsi.strip(name[i], massinput, concentration, r0[i], rho0[i], rs[i], rhos[i], hotdensity, velinfall, check)
        
            stripRNFW[i] = stripNFW
            stripRBurk[i] = stripBurk

            fracNFW[i], gasNFW[i], gmassNFW[i] = get.HIloss_trapz(name[i], stripNFW)
            fracBurk[i], gasBurk[i], gmassBurk[i] = get.HIloss_trapz(name[i], stripBurk)

            #Test whether KH instability will strip
            strip_pressure = (10**(-1*hotdensity))*(velinfall)**2 #units: cm**(-3) (km/s)**2

            DdepthNFW = stripNFW #This might need to be changed to something less obvious
            h1_avgNFW = gasNFW[i]/(np.pi*(stripNFW**2)*DdepthNFW)
            KH_restoreNFW = (((h1_avgNFW)*consts.G*NFWmass)/(2*np.pi*stripNFW))*((1e-4*consts.M_sun)*(u.s**2/u.m**3)*(consts.M_sun/consts.m_p)/(3.086*10**21)**4)
            
            DdepthBurk = stripBurk #This might need to be changed to something less obvious
            h1_avgBurk = gasBurk[i]/(np.pi*(stripBurk**2)*DdepthBurk)
            KH_restoreBurk = (((h1_avgBurk)*consts.G*Burkmass)/(2*np.pi*stripBurk))*((1e-4*consts.M_sun)*(u.s**2/u.m**3)*(consts.M_sun/consts.m_p)/(3.086*10**21)**4)

            if (KH_restoreNFW < strip_pressure):
                visc1_NFW[i] = (2e10)*((stripNFW/20.)**2)*((10**(-1*hotdensity))/(10**(-3)))*(velinfall/1000.0)
            else:
                visc1_NFW[i] = 0.0

            if (KH_restoreBurk < strip_pressure):
                visc1_Burk[i] = (2e10)*((stripBurk/20.)**2)*((10**(-1*hotdensity))/(10**(-3)))*(velinfall/1000.0)
            else:
                visc1_Burk[i] = 0.0

        ratio_redBurk = gasBurk / Mstar
        ratio_blueBurk = gmassBurk / Mstar

        KH1_Burk = np.empty(len(Mstar))
        #KH_blue = np.empty(len(Mstar))
        KH_red_Burk = np.empty(len(Mstar))
        #for i in range(len(Mstar)):
            #if gasNFW[i]==0:
                #KH1[i] = 0.0
            #else:
                #KH1[i] = visc1[i] / gasNFW[i]

        KH1_Burk = visc1_Burk / gmassBurk
        KH_red_Burk = visc1_Burk / Mstar


        mastertable = Table([name, Mstar, fracNFW, fracBurk, gasNFW, gasBurk, ratio_redBurk, ratio_blueBurk, KH1_Burk, KH_red_Burk], names=('NAME', 'Mstar', 'NFW', 'Burkert', 'gasNFW', 'gasBurk', 'Rgasfrac', 'Bgasfrac', 'KH_1Gyr', 'KH_red'), meta={'strippedfrac': 'RPS'})
    
        masteroutputfile = 'stripdata/RPS_strippedfrac_hot'+np.str(np.int(10*hotdensity))+'_v'+np.str(np.int(velinfall))+'_'+check+'.dat'
        mastertable.write(masteroutputfile, format = 'ascii')
###################################################################################################################################
               
    
    else:
        
        if check == 'LT':
        
            data = Table.read(inputfile, format = 'ascii')
            name = np.array(data['NAME'])
            vmag = np.array(data['Mv'])
            mstar = 10**((4.83-(vmag))/2.5)
            LTmassinput = 10**np.array(data['virmass'])
            con = np.array(data['c'])

            tdata = Table.read('THINGS_galaxylist.dat', format = 'ascii')
            tname = np.array(tdata['NAME'])
            tmstar = np.array(tdata['Mstar'])
            tmassinput = np.array(tdata['virmass'])
            tcon = np.array(tdata['c'])

            sdata = Table.read('SHIELD_list.dat', format = 'ascii')
            sname = np.array(sdata['NAME'])
            svmag = np.array(sdata['Mv'])
            smstar = 10**((4.83-(svmag))/2.5)
            #smstar = np.array(10**sdata['Mstar'])

            clean = (tmstar > 0.0)&(tmassinput > 0.0)&(tcon > 0.0)
            tcon_clean = tcon[clean]
            tmstar_clean = 10**tmstar[clean]
            tname_clean = tname[clean]
            tmassinput_clean = 10**tmassinput[clean]

            smassinput = np.interp(smstar, abMstar, abMpeak)
            sconcentration = 9.6*(smassinput/1.429e12)**(-0.075)
            
            concentration = np.append(con,tcon_clean)
            massinput = np.append(LTmassinput, tmassinput_clean)
            master_name = np.append(name, tname_clean)
            master_mstar = np.append(mstar, tmstar_clean)
            
            upper = concentration + 2.0
            lower = concentration - 2.0
        
        elif check == 'AB':
        
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

        elif check == 'pearson':
        
            data = Table.read('pearson_galaxylist.dat', format = 'ascii')
            name = np.array(data['NAME'])
            mstar = np.array(data['Mstar'])

            print 'Pearson data read in'
            mstar = 10**mstar
            massinput = np.interp(mstar, abMstar, abMpeak)
            concentration = 9.6*(massinput/1.429e12)**(-0.075)

            master_name = name
            master_mass = massinput
            master_con = concentration

            print 'Pearson data collected'

        else:
        
            data = Table.read('model_input.dat', format = 'ascii')
            name = np.array(data['NAME'])
            mstar = np.array(data['Mstar'])

            print 'Model data read in'
            mstar = 10**mstar
            massinput = np.interp(mstar, abMstar, abMpeak)
            concentration = 9.6*(massinput/1.429e12)**(-0.075)

            master_name = name
            master_mass = massinput
            master_con = concentration

            print 'final data collected'

        r0 = 1.0
        rho0 = 1.0
        rs = 1.0
        rhos = 1.0
    

        ########################################################
        # define the empty arrays which will be filled with data
        ########################################################
        fracNFW = np.empty(len(master_name))
        gasNFW = np.empty(len(master_name))
        gmassNFW = np.empty(len(master_name))
        upperfracNFW = np.empty(len(master_name))
        lowerfracNFW = np.empty(len(master_name))
        stripRNFW = np.empty(len(master_name))
        visc1 = np.empty(len(master_name))
        DMmassO = np.empty(len(master_name))

        if check == 'Model':
            Mstar = mstar
        elif check == 'pearson':
            Mstar = mstar
        else:
            Mstar1 = np.append(mstar, tmstar_clean)
            Mstar = np.append(Mstar1, smstar)

    
        ############################################################################
        #loop through each list of dwarf names and determine the relevant properties
        ############################################################################
        
        if check == 'LT':
        
            for i in range(len(master_name)):
            
                stripNFW, DMmass = tsi.strip(master_name[i], massinput[i], concentration[i], r0, rho0, rs, rhos, hotdensity, velinfall, check)
                #upperNFW = tsi.strip(master_name[i], massinput[i], upper[i], r0, rho0, rs, rhos, hotdensity, velinfall, check)
                #lowerNFW = tsi.strip(master_name[i], massinput[i], lower[i], r0, rho0, rs, rhos, hotdensity, velinfall, check)
            
                stripRNFW[i] = stripNFW
                
                fracNFW[i], gasNFW[i], gmassNFW[i] = get.HIloss_trapz(master_name[i], stripNFW)
                #upperfracNFW[i] = get.HIloss_trapz(master_name[i], upperNFW)
                #lowerfracNFW[i] = get.HIloss_trapz(master_name[i], lowerNFW)

                #Test whether KH instability will strip
                strip_pressure = (10**(-1*hotdensity))*(velinfall)**2 #units: cm**(-3) (km/s)**2
                print 'SP = '+np.str(strip_pressure)
                
                Ddepth = stripNFW #This might need to be changed to something less obvious
                h1_avg = gasNFW[i]/(np.pi*(stripNFW**2)*Ddepth)
                KH_restore = (((h1_avg)*consts.G*DMmass)/(2*np.pi*stripNFW))*((1e-4*consts.M_sun)*(u.s**2/u.m**3)*(consts.M_sun/consts.m_p)/(3.086*10**21)**4)
                print 'KH = '+np.str(KH_restore)
                
                if (KH_restore < strip_pressure):
                    visc1[i] = (2e10)*((stripNFW/20.)**2)*((10**(-1*hotdensity))/(10**(-3)))*(velinfall/1000.0)
                else:
                    visc1[i] = 0.0

            ratio_red = gasNFW / master_mstar
            ratio_blue = gmassNFW / master_mstar

            KH1 = np.empty(len(master_mstar))
            #KH_blue = np.empty(len(Mstar))
            KH_red = np.empty(len(master_mstar))
            #for i in range(len(Mstar)):
                #if gasNFW[i]==0:
                    #KH1[i] = 0.0
                #else:
                    #KH1[i] = visc1[i] / gasNFW[i]

            KH1 = visc1 / gmassNFW
            KH_red = visc1 / master_mstar

    
            mastertable = Table([master_name, master_mstar, fracNFW, ratio_red, ratio_blue, KH1, KH_red], names=('NAME', 'Mstar', 'NFW', 'Rgasfrac', 'Bgasfrac', 'KH_1Gyr', 'KH_red'), meta={'strippedfrac': 'RPS'})
        
            masteroutputfile = 'stripdata/RPS_strippedfrac_hot'+np.str(np.int(10*hotdensity))+'_v'+np.str(np.int(velinfall))+'_'+check+'.dat'
            mastertable.write(masteroutputfile, format = 'ascii')

        elif check == 'AB':

            for i in range(len(master_name)):

                #############################################################################################################################
                ## This determines the radial location in each dwarf where the restoring force can no longer overcome the ram pressure
                ## The outputs are the radial distance and the dark matter mass enclosed inside that radius.
                #############################################################################################################################
                print master_name[i]
                stripNFW, DMmass = tsi.strip(master_name[i], master_mass[i], master_con[i], r0, rho0, rs, rhos, hotdensity, velinfall, check)

                DMmassO[i] = DMmass
                print np.log10(DMmassO[i])
        
                stripRNFW[i] = stripNFW
                print stripNFW

                
                #############################################################################################################################
                ## This determines the amount of HI that is removed based on the radial distance determined above.
                ## The outputs are the gas fraction, gas mass removed, and the total gas mass that the dwarf began with.
                #############################################################################################################################
                fracNFW[i],gasNFW[i],gmassNFW[i] = get.HIloss_trapz(master_name[i], stripNFW)

                #############################################################################################################################
                ## This determines the amount of HI removed via viscous stripping, one single round
                #############################################################################################################################
                if fracNFW[i] == 1.0:
                    visc1[i] = 0.0
                    
                else:
                ##################################################################################################################
                ###Test whether KH instability will strip
                ## The strip pressure is constant and set by whatever MW properties we adopt
                ###############################################################################################
                    strip_pressure = (10**(-1*hotdensity))*(velinfall)**2 #units: cm**(-3) (km/s)**2
                    #print strip_pressure
                
                    Ddepth = stripNFW #This might need to be changed to something less obvious
                    h1_avg = gasNFW[i]/(np.pi*(stripNFW**2)*Ddepth)
                    KH_restore = (((h1_avg)*consts.G*DMmass)/(2*np.pi*stripNFW))*((1e-4*consts.M_sun)*(u.s**2/u.m**3)*(consts.M_sun/consts.m_p)/(3.086*10**21)**4)
                    #print KH_restore

                ###############################################################################################
                ## Not the best way to test this
                ###############################################################################################
                #F = gasNFW[i] / DMmassO[i]
                #M_kh = (1.6*10**(12))*((F/0.1)**(-7./2))*(((10**(-1*hotdensity))/(10**(-4.)))**(7./2))*((velinfall/1000.)**7)

                ###############################################################################################
                ## Check the inequality once, then viscous strip for 1 Gyr if True
                ###############################################################################################
                #if (KH_restore < strip_pressure):
                    #visc1[i] = (2e10)*((stripNFW/20.)**2)*((10**(-1*hotdensity))/(10**(-3)))*(velinfall/1000.0)
                #else:
                    #visc1[i] = 0.0

                #############################################################################################################################
                ## Below is the iterative viscous stripping routine
                ## Begin by selecting the timesteps, i.e. the total time and interval spacing
                ## Loop through each timestep and check the stripping inequality, determining the stripped fraction at each step
                ## Current problem is accurately determining the HI mass left after each step...go time
                #############################################################################################################################
                    totaltime = 1000.0
                    timeinterval = np.linspace(1, totaltime, 1000)
                    masslost = 0.0 #initialize the masslost to zero
                    visc1[i] = 0.0 #initialize the viscous stripping to zero
                    HImassremain = gasNFW[i]
                    print 'log HImassremain = '+np.str(np.log10(HImassremain))
                    print 'log DM mass = '+np.str(np.log10(DMmass))
                    print 'Rstrip = '+np.str(Ddepth)
                    print 'HI avg = '+np.str(h1_avg)

                    testdata = Table.read('THINGS_galaxylist.dat', format = 'ascii')
                    testname = np.array(testdata['NAME'])
                #############################################################################################################################
                ## These are arrays of values, specifically the mass lost for each additional timestep. This is cumulative...
                #############################################################################################################################
                
                    for j in range(len(timeinterval)):

                        print 'j = '+np.str(j)
                    
                        dt = 0.001
                    
                        if (KH_restore < strip_pressure):
                            masslost = ((2e10)*((stripNFW/20.)**2)*((10**(-1*hotdensity))/(10**(-3)))*(velinfall/1000.0))*dt

                            if ((masslost > HImassremain) & (j>0)):
                                masslost = HImassremain
                            else: 
                                masslost = masslost
                        else:
                            masslost = 0.0
                        
                        print 'log masslost = '+np.str(np.log10(masslost))
                        
                    #############################################################################################################################
                    ## Remove the mass lost from the total HI mass, and determine the new strip radius (stripNFW), HI remaining (gasNFW),
                    ##   and total dark matter mass inside stripNFW (DMmass)
                    #############################################################################################################################
                        massremain = HImassremain - masslost
                    
                        if massremain <= 0:
                            visc1[i] = visc1[i] + HImassremain
                            print 'log mass removed = ' + np.str(np.log10(visc1[i]))
                            print 'mass remaining is zero, halo no longer exists'
                            #break
                        else:
                            visc1[i] = visc1[i] + masslost
                            print 'log mass remain = ' + np.str(np.log10(massremain))
                            print 'log mass removed = ' + np.str(np.log10(visc1[i]))
                            print 'halo still exist...carry on!'
                            #continue

                            data = Table.read('data/'+master_name[i]+'_clean.dat', format = 'ascii')
                            radius = np.array(data['R(kpc)'])
                            gas = np.array(data['N_HI'])

                            H2check = master_name[i] == testname
                            test = testname[H2check]

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
                            Rlimit = np.max(radiusfinal)

                            all_mass = np.log10(integrate.cumtrapz(func,radiusfinal)*(consts.m_p / consts.M_sun))
                            #print massremain
                    
                            cut = np.where(all_mass < np.log10(massremain))[0]
                            #print all_mass[cut]

                            #### This ensures the new radius is always a real value and not an empty array
                            if (len(radiusfinal[cut]) == 0):
                                newR = np.min(radiusfinal)
                                stripNFW = newR
                                print 'new Rism = '+ np.str(newR)
                                #print np.max(radiusfinal[cut])
                                HImassremain = 10**np.min(all_mass)
                                print 'new log HImassremain = '+ np.str(np.log10(HImassremain))
                            else:
                                newR = np.max(radiusfinal[cut])
                                stripNFW = newR
                                print 'new Rism = '+ np.str(newR)
                                #print np.max(radiusfinal[cut])
                                HImassremain = np.max(10**all_mass[cut])
                                print 'new log HImassremain = '+ np.str(np.log10(HImassremain))

                    #############################################################################################################################
                    ## New Dark Matter mass definition
                    #############################################################################################################################
                            NFWmass = massprofile.NFW(master_mass[i], master_con[i], radiusfinal)
                            DMcut = np.where(radiusfinal == newR)[0]
                            newDMmass = NFWmass[DMcut]

                            print 'new log DMmass = ' + np.str(np.log10(newDMmass))

                    #############################################################################################################################
                    ## Define new variables after the HI mass has been removed.
                    #############################################################################################################################
                            Ddepth = newR 
                            h1_avg = HImassremain/(np.pi*(newR**2)*Ddepth)
                            print 'new HI avg = '+np.str(h1_avg)
                            KH_restore = (((h1_avg)*consts.G*newDMmass)/(2*np.pi*newR))*((1e-4*consts.M_sun)*(u.s**2/u.m**3)*(consts.M_sun/consts.m_p)/(3.086*10**21)**4)
                    

            ratio_red = gasNFW / Mstar
            ratio_blue = gmassNFW / Mstar

            KH1 = np.empty(len(Mstar))
            #KH_blue = np.empty(len(Mstar))
            KH_red = np.empty(len(Mstar))
            #for i in range(len(Mstar)):
                #if gasNFW[i]==0:
                    #KH1[i] = 0.0
                #else:
                    #KH1[i] = visc1[i] / gasNFW[i]

            KH1 = visc1 / gmassNFW
            KH_red = visc1 / Mstar
               
            
            mastertable = Table([master_name, Mstar, fracNFW, ratio_red, ratio_blue, KH1, KH_red], names=('NAME', 'Mstar', 'NFW', 'Rgasfrac', 'Bgasfrac', 'KH_1Gyr', 'KH_red'), meta={'strippedfrac': 'RPS'})
            
            masteroutputfile = 'stripdata/RPS_strippedfrac_hot'+np.str(np.int(10*hotdensity))+'_v'+np.str(np.int(velinfall))+'_'+check+'iterate1000.dat'
            mastertable.write(masteroutputfile, format = 'ascii')

        elif check == 'pearson':

            for i in range(len(master_name)):

                #############################################################################################################################
                ## This determines the radial location in each dwarf where the restoring force can no longer overcome the ram pressure
                ## The outputs are the radial distance and the dark matter mass enclosed inside that radius.
                #############################################################################################################################
                print master_name[i]
                stripNFW, DMmass = tsi.strip(master_name[i], master_mass[i], master_con[i], r0, rho0, rs, rhos, hotdensity, velinfall, check)

                DMmassO[i] = DMmass
                print np.log10(DMmassO[i])
        
                stripRNFW[i] = stripNFW
                print stripNFW

                
                #############################################################################################################################
                ## This determines the amount of HI that is removed based on the radial distance determined above.
                ## The outputs are the gas fraction, gas mass removed, and the total gas mass that the dwarf began with.
                #############################################################################################################################
                fracNFW[i],gasNFW[i],gmassNFW[i] = get.HIloss_trapz(master_name[i], stripNFW, check)

                #############################################################################################################################
                ## This determines the amount of HI removed via viscous stripping, one single round
                #############################################################################################################################
                if fracNFW[i] == 1.0:
                    visc1[i] = 0.0
                    
                else:
                ##################################################################################################################
                ###Test whether KH instability will strip
                ## The strip pressure is constant and set by whatever MW properties we adopt
                ###############################################################################################
                    strip_pressure = (10**(-1*hotdensity))*(velinfall)**2 #units: cm**(-3) (km/s)**2
                    #print strip_pressure
                
                    Ddepth = stripNFW #This might need to be changed to something less obvious
                    h1_avg = gasNFW[i]/(np.pi*(stripNFW**2)*Ddepth)
                    KH_restore = (((h1_avg)*consts.G*DMmass)/(2*np.pi*stripNFW))*((1e-4*consts.M_sun)*(u.s**2/u.m**3)*(consts.M_sun/consts.m_p)/(3.086*10**21)**4)
                    #print KH_restore

                ###############################################################################################
                ## Not the best way to test this
                ###############################################################################################
                #F = gasNFW[i] / DMmassO[i]
                #M_kh = (1.6*10**(12))*((F/0.1)**(-7./2))*(((10**(-1*hotdensity))/(10**(-4.)))**(7./2))*((velinfall/1000.)**7)

                ###############################################################################################
                ## Check the inequality once, then viscous strip for 1 Gyr if True
                ###############################################################################################
                #if (KH_restore < strip_pressure):
                    #visc1[i] = (2e10)*((stripNFW/20.)**2)*((10**(-1*hotdensity))/(10**(-3)))*(velinfall/1000.0)
                #else:
                    #visc1[i] = 0.0

                #############################################################################################################################
                ## Below is the iterative viscous stripping routine
                ## Begin by selecting the timesteps, i.e. the total time and interval spacing
                ## Loop through each timestep and check the stripping inequality, determining the stripped fraction at each step
                ## Current problem is accurately determining the HI mass left after each step...go time
                #############################################################################################################################
                    totaltime = 1000.0
                    timeinterval = np.linspace(1, totaltime, 100)
                    masslost = 0.0 #initialize the masslost to zero
                    visc1[i] = 0.0 #initialize the viscous stripping to zero
                    HImassremain = gasNFW[i]
                    print 'log HImassremain = '+np.str(np.log10(HImassremain))
                    print 'log DM mass = '+np.str(np.log10(DMmass))
                    print 'Rstrip = '+np.str(Ddepth)
                    print 'HI avg = '+np.str(h1_avg)

                    testdata = Table.read('THINGS_galaxylist.dat', format = 'ascii')
                    testname = np.array(testdata['NAME'])
                #############################################################################################################################
                ## These are arrays of values, specifically the mass lost for each additional timestep. This is cumulative...
                #############################################################################################################################
                
                    for j in range(len(timeinterval)):

                        print 'j = '+np.str(j)
                    
                        dt = 0.001
                    
                        if (KH_restore < strip_pressure):
                            masslost = ((2e10)*((stripNFW/20.)**2)*((10**(-1*hotdensity))/(10**(-3)))*(velinfall/1000.0))*dt

                            if ((masslost > HImassremain) & (j>0)):
                                masslost = HImassremain
                            else: 
                                masslost = masslost
                        else:
                            masslost = 0.0
                        
                        print 'log masslost = '+np.str(np.log10(masslost))
                        
                    #############################################################################################################################
                    ## Remove the mass lost from the total HI mass, and determine the new strip radius (stripNFW), HI remaining (gasNFW),
                    ##   and total dark matter mass inside stripNFW (DMmass)
                    #############################################################################################################################
                        massremain = HImassremain - masslost
                    
                        if massremain <= 0:
                            visc1[i] = visc1[i] + HImassremain
                            print 'log mass removed = ' + np.str(np.log10(visc1[i]))
                            print 'mass remaining is zero, halo no longer exists'
                            #break
                        else:
                            visc1[i] = visc1[i] + masslost
                            print 'log mass remain = ' + np.str(np.log10(massremain))
                            print 'log mass removed = ' + np.str(np.log10(visc1[i]))
                            print 'halo still exist...carry on!'
                            #continue

                            data = Table.read(master_name[i]+'.txt', format = 'ascii')
                            radius = np.array(data['col1'])
                            gas = np.array(np.log10(data['col2']))

                            H2check = master_name[i] == testname
                            test = testname[H2check]

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

                            func = 0.5*np.pi*(radiusfinal*((3.08E21)**2))*(10**gasfinal)
                            mass = np.trapz(func,radiusfinal)
                            Rlimit = np.max(radiusfinal)

                            all_mass = np.log10(integrate.cumtrapz(func,radiusfinal)*(consts.m_p / consts.M_sun))
                            #print massremain
                    
                            cut = np.where(all_mass < np.log10(massremain))[0]
                            #print all_mass[cut]

                            #### This ensures the new radius is always a real value and not an empty array
                            if (len(radiusfinal[cut]) == 0):
                                newR = np.min(radiusfinal)
                                stripNFW = newR
                                print 'new Rism = '+ np.str(newR)
                                #print np.max(radiusfinal[cut])
                                HImassremain = 10**np.min(all_mass)
                                print 'new log HImassremain = '+ np.str(np.log10(HImassremain))
                            else:
                                newR = np.max(radiusfinal[cut])
                                stripNFW = newR
                                print 'new Rism = '+ np.str(newR)
                                #print np.max(radiusfinal[cut])
                                HImassremain = np.max(10**all_mass[cut])
                                print 'new log HImassremain = '+ np.str(np.log10(HImassremain))

                    #############################################################################################################################
                    ## New Dark Matter mass definition
                    #############################################################################################################################
                            NFWmass = massprofile.NFW(master_mass[i], master_con[i], radiusfinal)
                            DMcut = np.where(radiusfinal == newR)[0]
                            newDMmass = NFWmass[DMcut]

                            print 'new log DMmass = ' + np.str(np.log10(newDMmass))

                    #############################################################################################################################
                    ## Define new variables after the HI mass has been removed.
                    #############################################################################################################################
                            Ddepth = newR 
                            h1_avg = HImassremain/(np.pi*(newR**2)*Ddepth)
                            print 'new HI avg = '+np.str(h1_avg)
                            KH_restore = (((h1_avg)*consts.G*newDMmass)/(2*np.pi*newR))*((1e-4*consts.M_sun)*(u.s**2/u.m**3)*(consts.M_sun/consts.m_p)/(3.086*10**21)**4)
                    

            ratio_red = gasNFW / Mstar
            ratio_blue = gmassNFW / Mstar

            KH1 = np.empty(len(Mstar))
            #KH_blue = np.empty(len(Mstar))
            KH_red = np.empty(len(Mstar))
            #for i in range(len(Mstar)):
                #if gasNFW[i]==0:
                    #KH1[i] = 0.0
                #else:
                    #KH1[i] = visc1[i] / gasNFW[i]

            KH1 = visc1 / gmassNFW
            KH_red = visc1 / Mstar
               
            
            mastertable = Table([master_name, Mstar, fracNFW, ratio_red, ratio_blue, KH1, KH_red], names=('NAME', 'Mstar', 'NFW', 'Rgasfrac', 'Bgasfrac', 'KH_1Gyr', 'KH_red'), meta={'strippedfrac': 'RPS'})
            
            masteroutputfile = 'RPS_strippedfrac_hot'+np.str(np.int(10*hotdensity))+'_v'+np.str(np.int(velinfall))+'_'+check+'iterate100.dat'
            mastertable.write(masteroutputfile, format = 'ascii')

        else:

            for i in range(len(master_name)):

                print master_name[i]
                stripNFW, DMmass = tsi.strip(master_name[i], master_mass[i], master_con[i], r0, rho0, rs, rhos, hotdensity, velinfall, check)

                DMmassO[i] = DMmass
                print np.log10(DMmassO[i])
        
                stripRNFW[i] = stripNFW
                print stripNFW
            
                fracNFW[i],gasNFW[i],gmassNFW[i] = get.HIloss_trapz(master_name[i], stripNFW)

                #Test whether KH instability will strip
                strip_pressure = (10**(-1*hotdensity))*(velinfall)**2 #units: cm**(-3) (km/s)**2
                print strip_pressure
                
                Ddepth = stripNFW #This might need to be changed to something less obvious
                h1_avg = gasNFW[i]/(np.pi*(stripNFW**2)*Ddepth)
                KH_restore = (((h1_avg)*consts.G*DMmass)/(2*np.pi*stripNFW))*((1e-4*consts.M_sun)*(u.s**2/u.m**3)*(consts.M_sun/consts.m_p)/(3.086*10**21)**4)
                print KH_restore
                
                #F = gasNFW[i] / DMmassO[i]
                #M_kh = (1.6*10**(12))*((F/0.1)**(-7./2))*(((10**(-1*hotdensity))/(10**(-4.)))**(7./2))*((velinfall/1000.)**7)

                if (KH_restore < strip_pressure):
                    visc1[i] = (2e10)*((stripNFW/20.)**2)*((10**(-1*hotdensity))/(10**(-3)))*(velinfall/1000.0)
                else:
                    visc1[i] = 0.0

            ratio_red = gasNFW / Mstar
            ratio_blue = gmassNFW / Mstar

            KH1 = np.empty(len(Mstar))
            #KH_blue = np.empty(len(Mstar))
            KH_red = np.empty(len(Mstar))
            #for i in range(len(Mstar)):
                #if gasNFW[i]==0:
                    #KH1[i] = 0.0
                #else:
                    #KH1[i] = visc1[i] / gasNFW[i]

            KH1 = visc1 / gmassNFW
            KH_red = visc1 / Mstar
               
            
            mastertable = Table([master_name, Mstar, fracNFW, ratio_red, ratio_blue, KH1, KH_red], names=('NAME', 'Mstar', 'NFW', 'Rgasfrac', 'Bgasfrac', 'KH_1Gyr', 'KH_red'), meta={'strippedfrac': 'RPS'})
            
            masteroutputfile = 'modeldata/RPS_strippedfrac_hot'+np.str(np.int(10*hotdensity))+'_v'+np.str(np.int(velinfall))+'_'+check+'.dat'
            mastertable.write(masteroutputfile, format = 'ascii')

        
        print 'balls'

