import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
import code.test_stripping_individual as tsi
import code.get_HImassloss_trapz as get



def frac_vs_mstar(hotdensity, velinfall):
    
    ab_match = Table.read('GK14AM.txt', format = 'ascii')
    abMpeak = np.array(ab_match['Mpeak(Msun)'])
    abMstar = np.array(ab_match['Mstar(Msun)'])
    
#############################################
# Define inputs for each plot
#############################################

    #use for the ABmatching analysis and Isothermal sphere
    highmass = np.array(['ddo47', 'ddo50', 'ddo52', 'ddo165', 'ddo168', 'ngc2366', 'ngc3738', 'ngc4214'])
    highmassinput = np.array([-15.5, -16.6, -15.4, -15.6, -15.7, -16.8, -17.1, -17.6])
    #use for the LITTLE THINGS based NFW profile, data from Oh et al. 2015
    LThighmass = np.array(['ddo47', 'ddo50', 'ddo52', 'ddo168', 'ngc2366', 'ngc3738'])#no data for ddo165 and ngc4214
    LThighmass_star = np.array([-15.5, -16.6, -15.4, -15.7, -16.8, -17.1])
    LThighmassinput = np.array([9.930, 9.463, 9.664, 9.520, 9.841, 9.838])#no data for ddo165 and ngc4214
    LThighconcentration = np.array([9.2, 10.3, 9.4, 9.2, 9.4, 8.1]) #no data for ddo165 and ngc4214
    
    #use for the ABmatching analysis and Isothermal sphere
    lowmass7 = np.array(['ddo43', 'ddo46', 'ddo53', 'ddo63', 'ddo70', 'ddo75', 'ddo87', 'ddo101', 'ddo126', 'ddo133', 'ddo154', 'ddo167', 'ddo187', 'ddo216', 'f564-v3', 'ngc4163', 'ugc8508', 'wlm'])
    lowmass7input = np.array([-15.1, -14.7, -13.8, -14.8, -14.1, -13.9, -15.0, -15.0, -14.9, -14.8, -14.2, -13.0, -12.7, -13.7, -14.0, -14.4, -13.6, -14.4])
    #use for the LITTLE THINGS based NFW profile, data from Oh et al. 2015
    LTlowmass7 = np.array(['ddo43', 'ddo46', 'ddo53', 'ddo70','ddo87', 'ddo101', 'ddo126', 'ddo133', 'ddo154', 'ddo216', 'f564-v3', 'ugc8508', 'wlm'])
    LTlowmass7_star = np.array([-15.1, -14.7, -13.8, -14.1, -15.0, -15.0, -14.9, -14.8, -14.2, -13.7, -14.0, -13.6, -14.4])
    LTlowmass7input = np.array([9.159, 9.609, 8.567, 8.768, 9.734, 9.279, 9.182, 9.261, 9.614, 7.807, 9.003, 8.913, 9.002])
    LTlowconcentration7 = np.array([10.2, 9.0, 10.4, 10.3, 9.5, 9.4, 10.1, 9.8, 9.8, 11.7, 10.4, 10.0, 10.2])
    
    #use for the ABmatching analysis and Isothermal sphere
    lowmass6 = np.array(['CVnIdwA', 'ddo69', 'ddo155', 'ddo210', 'M81dwA', 'SagDIG'])
    lowmass6input = np.array([-12.4, -11.7, -12.5, -10.9, -11.7, -12.5])
    #use for the LITTLE THINGS based NFW profile, data from Oh et al. 2015
    LTlowmass6 = np.array(['CVnIdwA', 'ddo210'])
    LTlowmass6_star = np.array([-12.4, -10.9])
    LTlowmass6input = np.array([8.529, 6.964])
    LTlowconcentration6 = np.array([11.0, 12.2]) #missing ddo69, ddo155, M81dwA, and SagDIG: assume all have the same concentration as CVnIdwA
    
    lowmass = np.append(lowmass7,lowmass6)
    lowmassinput = np.append(lowmass7input, lowmass6input)
    LTlowmass = np.append(LTlowmass7, LTlowmass6)
    LTlowmass_star = np.append(LTlowmass7_star, LTlowmass6_star)
    LTlowmassinput = np.append(LTlowmass7input, LTlowmass6input)
    LTlowconcentration = np.append(LTlowconcentration7, LTlowconcentration6)

########################################################
# define the empty arrays which will be filled with data for figures
########################################################
    highfrac = np.empty(len(highmass))
    highfracISO = np.empty(len(highmass))
    highfrac2 = np.empty(len(highmass))
    highfracB = np.empty(len(highmass))
    highstripR = np.empty(len(highmass))
    highstripRISO = np.empty(len(highmass))
    highstripR2 = np.empty(len(highmass))
    highstripRBurk = np.empty(len(highmass))
    highMstar = np.empty(len(highmass))
    
    LThighfrac = np.empty(len(LThighmass))
    LThighfracISO = np.empty(len(LThighmass))
    LThighfrac2 = np.empty(len(LThighmass))
    LThighstripR = np.empty(len(LThighmass))
    LThighstripRISO = np.empty(len(LThighmass))
    LThighstripR2 = np.empty(len(LThighmass))
    LThighMstar = np.empty(len(LThighmass))
    
    lowfrac = np.empty(len(lowmass))
    lowfracISO = np.empty(len(lowmass))
    lowfrac2 = np.empty(len(lowmass))
    lowfracB = np.empty(len(lowmass))
    lowstripR = np.empty(len(lowmass))
    lowstripRISO = np.empty(len(lowmass))
    lowstripR2 = np.empty(len(lowmass))
    lowstripRBurk = np.empty(len(lowmass))
    lowMstar = np.empty(len(lowmass))
    
    LTlowfrac = np.empty(len(LTlowmass))
    LTlowfracISO = np.empty(len(LTlowmass))
    LTlowfrac2 = np.empty(len(LTlowmass))
    LTlowstripR = np.empty(len(LTlowmass))
    LTlowstripRISO = np.empty(len(LTlowmass))
    LTlowstripR2 = np.empty(len(LTlowmass))
    LTlowMstar = np.empty(len(LTlowmass))
    
##############################################################
#loop through each list of dwarf names and determine the relevant properties
##############################################################

    for i in range(len(highmass)):
        
        mstar = 10**((4.83-(highmassinput[i]))/2.5)
        massinput = np.interp(mstar, abMstar, abMpeak)
        highconcentration = 9.6*(massinput/1.429e12)**(-0.075)
        hot = hotdensity
        vel = velinfall
        #print massinput
        
        stripradius2, stripradiusISO = tsi.strip(highmass[i], massinput, highconcentration, hot, vel)
        
        #highstripR[i] = stripradius
        highstripRISO[i] = stripradiusISO
        highstripR2[i] = stripradius2
        #highstripRBurk = stripradiusB
        
        highMstar[i] = mstar
        
        #highfrac[i] = get.HIloss_trapz(highmass[i], stripradius)
        highfracISO[i] = get.HIloss_trapz(highmass[i], stripradiusISO)
        highfrac2[i] = get.HIloss_trapz(highmass[i], stripradius2)
        #highfracB[i] = get.HIloss_trapz(highmass[i], stripradiusB)
        
    for i in range(len(LThighmass)):
        
        mstar = 10**((4.83-(LThighmass_star[i]))/2.5)
        LThighMstar[i] = mstar
        
        LTstripradius2, LTstripradiusISO = tsi.strip(LThighmass[i], 10**LThighmassinput[i], LThighconcentration[i], hotdensity, velinfall)
        
        #LThighstripR[i] = LTstripradius
        LThighstripRISO[i] = LTstripradiusISO
        LThighstripR2[i] = LTstripradius2

    
        #LThighfrac[i] = get.HIloss_trapz(LThighmass[i], LTstripradius)
        LThighfracISO[i] = get.HIloss_trapz(LThighmass[i], LTstripradiusISO)
        LThighfrac2[i] = get.HIloss_trapz(LThighmass[i], LTstripradius2)

    for i in range(len(lowmass)):
        
        mstar = 10**((4.83-(lowmassinput[i]))/2.5)
        massinput = np.interp(mstar, abMstar, abMpeak)
        lowconcentration = 9.6*(massinput/1.429e12)**(-0.075)
        hot = hotdensity
        vel = velinfall
        #print massinput
        
        stripradius2, stripradiusISO = tsi.strip(lowmass[i], massinput, lowconcentration, hot, vel)

        
        #lowstripR[i] = stripradius
        lowstripRISO[i] = stripradiusISO
        lowstripR2[i] = stripradius2
        #lowstripRBurk[i] = stripradiusB
        
        lowMstar[i] = mstar
        
        #lowfrac[i] = get.HIloss_trapz(lowmass[i], stripradius)
        lowfracISO[i] = get.HIloss_trapz(lowmass[i], stripradiusISO)
        lowfrac2[i] = get.HIloss_trapz(lowmass[i], stripradius2)
        #lowfracB[i] = get.HIloss_trapz(lowmass[i], stripradiusB)
        
    for i in range(len(LTlowmass)):
    
        mstar = 10**((4.83-(LTlowmass_star[i]))/2.5)
        LTlowMstar[i] = mstar
        
        LTstripradius2, LTstripradiusISO = tsi.strip(LTlowmass[i], 10**LTlowmassinput[i], LTlowconcentration[i], hotdensity, velinfall)
        
        #LTlowstripR[i] = LTstripradius
        LTlowstripRISO[i] = LTstripradiusISO
        LTlowstripR2[i] = LTstripradius2
        
        #LTlowfrac[i] = get.HIloss_trapz(LTlowmass[i], LTstripradius)
        LTlowfracISO[i] = get.HIloss_trapz(LTlowmass[i], LTstripradiusISO)
        LTlowfrac2[i] = get.HIloss_trapz(LTlowmass[i], LTstripradius2)


##################################################
## Begin plotting routine
##################################################

#begin plotting routine for ABmatching and Isothermal sphere
#plt.figure(1)
#plt.axis([6.0, 9.0, 0.0, 1.1])
#plt.xlabel(r'$\rm Stellar\ Mass$')
#plt.ylabel(r'$f_{\rm stripped}$')
    
#ph1 = plt.plot(np.log10(highMstar), highfrac, 'go', markersize = 10.0)
    #phISO = plt.plot(np.log10(highMstar), highfracISO, 'gs', markersize = 10.0)
#ph2 = plt.plot(np.log10(highMstar), highfrac2, 'g*', markersize = 10.0)
#pl1 = plt.plot(np.log10(lowMstar), lowfrac, 'bo', markersize = 10.0)
    #plISO = plt.plot(np.log10(lowMstar), lowfracISO, 'bs', markersize = 10.0)
#pl2 = plt.plot(np.log10(lowMstar), lowfrac2, 'b*', markersize = 10.0)
    #plt.plot(np.log10(lowMstar7), lowfrac7, 'co', linewidth = 4.0, markersize = 10.0)
    #plt.plot(np.log10(lowMstar7), lowfrac72, 'c*', linewidth = 4.0, markersize = 10.0)
    #plt.plot(np.log10(lowMstar6), lowfrac6, 'mo', linewidth = 4.0, markersize = 10.0)
    #plt.plot(np.log10(lowMstar6), lowfrac62, 'm*', linewidth = 4.0, markersize = 10.0)

    #plt.legend([p1, p2, p3, p4], [r"$\rm Isothermal,\ \sigma = 50\ km/s$", r"$\rm NFW,\ c=12$",r"$\rm Isothermal,\ \sigma = 30\ km/s$", r"$\rm NFW,\ c=15$"], loc = 4, frameon=False, numpoints=1, prop={'size':20}, handletextpad = 0.1)
    
    #plt.savefig('frac_vs_mstar_IsoNFW_hot35.pdf')

#begin plotting routine for ABmatching NFW profile and LITTLE THINGS NFW profile
    plt.figure(2)
    plt.axis([6.0, 10.0, 0.0, 1.1])
    plt.xlabel(r'$\rm Stellar\ Mass$')
    plt.ylabel(r'$f_{\rm stripped}$')

    #pLTh2 = plt.plot(np.log10(LThighMstar), LThighfrac2, 'g^', markersize = 10.0)
    ph2 = plt.plot(np.log10(highMstar), highfrac2, 'g*', markersize = 10.0)
    #pLTl2 = plt.plot(np.log10(LTlowMstar), LTlowfrac2, 'b^', markersize = 10.0)
    pl2 = plt.plot(np.log10(lowMstar), lowfrac2, 'g*', markersize = 10.0)
    #plt.plot(np.log10(lowMstar7), lowfrac7, 'co', linewidth = 4.0, markersize = 10.0)
    #plt.plot(np.log10(lowMstar7), lowfrac72, 'c*', linewidth = 4.0, markersize = 10.0)
    #plt.plot(np.log10(lowMstar6), lowfrac6, 'mo', linewidth = 4.0, markersize = 10.0)
    #plt.plot(np.log10(lowMstar6), lowfrac62, 'm*', linewidth = 4.0, markersize = 10.0)

    #plt.legend([p1, p2, p3, p4], [r"$\rm Isothermal,\ \sigma = 50\ km/s$", r"$\rm NFW,\ c=12$",r"$\rm Isothermal,\ \sigma = 30\ km/s$", r"$\rm NFW,\ c=15$"], loc = 4, frameon=False, numpoints=1, prop={'size':20}, handletextpad = 0.1)

    plt.savefig('frac_vs_mstar_NFWab_hot'+np.str(np.int(hotdensity*10))+'.pdf')

    plt.figure(3)
    plt.axis([6.0, 10.0, 0.0, 1.1])
    plt.xlabel(r'$\rm Stellar\ Mass$')
    plt.ylabel(r'$f_{\rm stripped}$')
    
    pLTh2 = plt.plot(np.log10(LThighMstar), LThighfrac2, 'g^', markersize = 10.0)
    #ph2 = plt.plot(np.log10(highMstar), highfrac2, 'g*', markersize = 10.0)
    pLTl2 = plt.plot(np.log10(LTlowMstar), LTlowfrac2, 'g^', markersize = 10.0)
    #pl2 = plt.plot(np.log10(lowMstar), lowfrac2, 'b*', markersize = 10.0)
    #plt.plot(np.log10(lowMstar7), lowfrac7, 'co', linewidth = 4.0, markersize = 10.0)
    #plt.plot(np.log10(lowMstar7), lowfrac72, 'c*', linewidth = 4.0, markersize = 10.0)
    #plt.plot(np.log10(lowMstar6), lowfrac6, 'mo', linewidth = 4.0, markersize = 10.0)
    #plt.plot(np.log10(lowMstar6), lowfrac62, 'm*', linewidth = 4.0, markersize = 10.0)
    
    #plt.legend([p1, p2, p3, p4], [r"$\rm Isothermal,\ \sigma = 50\ km/s$", r"$\rm NFW,\ c=12$",r"$\rm Isothermal,\ \sigma = 30\ km/s$", r"$\rm NFW,\ c=15$"], loc = 4, frameon=False, numpoints=1, prop={'size':20}, handletextpad = 0.1)
    
    plt.savefig('frac_vs_mstar_LT_hot'+np.str(np.int(hotdensity*10))+'.pdf')

    
    plt.show()


    #return highstripR, highstripR2, highfrac, highfrac2, highMstar, lowstripR, lowstripR2, lowfrac, lowfrac2, lowMstar












