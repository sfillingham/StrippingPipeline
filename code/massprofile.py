import numpy as np
import astropy.constants as consts
import astropy.units as u

def Rvir(input):

    delvir = 337.0
    omegaM = 0.3
    h = 0.7
    massconversion = consts.M_sun*(1e3)/u.kg
    rho_crit = (1.879e-29)*h**2 #units: g cm^-3
    rho_universe = omegaM*rho_crit #average mass density of the universe
    #virialradius = ((3.0*input*massconversion)/(4.0*np.pi*delvir*rho_universe))**(1.0/3)/consts.kpc*u.m*1e-2
    
    #this virial radius uses James's simpler normalization.
    virialradius = 300.0*((input/(1.0e12))**(1/3.0))

    return virialradius #units: kpc

######################################
def Rvir2(input):

    delvir = 337.0
    omegaM = 0.3
    h = 0.7
    massconversion = consts.M_sun*(1e3)/u.kg
    rho_crit = (1.879e-29)*h**2 #units: g cm^-3
    rho_universe = omegaM*rho_crit #average mass density of the universe
    #virialradius = ((3.0*input*massconversion)/(4.0*np.pi*delvir*rho_universe))**(1.0/3)/consts.kpc*u.m*1e-2
    
    virialradius = 300.0*(input/200.0)

    return virialradius #units: kpc

######################################
def M200(input):

    delvir = 337.0
    omegaM = 0.3
    h = 0.7
    massconversion = consts.M_sun*(1e3)/u.kg
    rho_crit = (1.879e-29)*h**2 #units: g cm^-3
    rho_universe = omegaM*rho_crit #average mass density of the universe
    #virialradius = ((3.0*input*massconversion)/(4.0*np.pi*delvir*rho_universe))**(1.0/3)/consts.kpc*u.m*1e-2
    
    virialmass = (1e12)*((input/200.0)**2)*(Rvir2(input)/300.0)

    return virialmass #units: kpc

######################################


def NFW(virmass, concentration, x):
    
    delvir = 337.0
    omegaM = 0.3
    h = 0.7
    massconversion = consts.M_sun*(1e3)/u.kg #takes Msun to grams
    rho_crit = (1.879e-29)*h**2 #units: g cm^-3
    rho_universe = omegaM*rho_crit*1e3*(consts.kpc/u.m)**3/(consts.M_sun/u.kg) #average mass density of the universe, convert to Msun kpc^-3
    rho_scale = delvir*rho_universe*(concentration**3)/(3*(np.log(1+concentration)-concentration/(1+concentration))) #units: Msun kpc^-3
    
    rs = Rvir(virmass)/concentration #units: kpc
    #print rs

    #mass = 4*np.pi*rho_scale*(rs**3)*(np.log(1+(x/rs)) - (x/rs)/(1+(x/rs)))
    
    #This mass function uses James's simple approximation
    mass = virmass*((np.log(1+(x/rs)) - (x/rs)/(1+(x/rs)))/(np.log(1+concentration) - concentration/(1+concentration)))

    return mass #units of Msun

######################################

def NFW_pace(rs, rhos, x):
    
    mass = 4*np.pi*rhos*(rs**3.)*((np.log(1+(x/rs))) - (x/rs)/(1+(x/rs)))

    return mass

########################################

def burkert(r0, rho0, x):
    
    mass = 4*np.pi*rho0*(r0**3.)*((np.log(1+(x/r0))/2.) + (np.log(1+(x/r0)**2.)/4.) - np.arctan(x/r0)/2.)
    
    return mass

######################################

def Isothermal(virmass,x):
    
    sigmasquare = (consts.G*virmass)/(2*Rvir(virmass))

    #mass = (2*(sigmasquare)*x/consts.G)*(1e6*consts.kpc)*((u.m/u.s)**2)/(consts.M_sun)
    
    mass = (virmass/Rvir(virmass))*x

    return mass

########################################


