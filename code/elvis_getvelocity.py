#This script take the interpolated output files from the ELVIS suite for a given mass range and determines the distribution of satellite
#velocities relative to the host in which they reside today.
#
#Note: This can also be used to determine the distribution of pericenter distances. This feature should be commented out
#      when the radial distance cut is anything other than pericenter.
#
import numpy as np
from astropy.table import Table
from astropy.cosmology import LambdaCDM
import astropy.units as u
cosmo = LambdaCDM(H0=71,Om0=0.266,Ode0=0.734)

def work(vel):

    host1data = np.loadtxt('/Users/Sean/research/elvis/FieldPaper/bigbox/mass59/elvis_alldwarfs_v'+vel+'_distance_time_interp_host1.dat')
    host2data = np.loadtxt('/Users/Sean/research/elvis/FieldPaper/bigbox/mass59/elvis_alldwarfs_v'+vel+'_distance_time_interp_host2.dat')
    vel1data = np.loadtxt('/Users/Sean/research/elvis/FieldPaper/bigbox/mass59/elvis_alldwarfs_v'+vel+'_velocity_time_interp_host1.dat')
    vel2data = np.loadtxt('/Users/Sean/research/elvis/FieldPaper/bigbox/mass59/elvis_alldwarfs_v'+vel+'_velocity_time_interp_host2.dat')
    scale = np.loadtxt('/Users/Sean/research/elvis/FieldPaper/bigbox/mass59/elvis_alldwarfs_v'+vel+'_scalefactor_interp.dat')
    vir1data = np.loadtxt('/Users/Sean/research/elvis/FieldPaper/bigbox/mass59/elvis_alldwarfs_v'+vel+'_virial1_interp.dat')
    vir2data = np.loadtxt('/Users/Sean/research/elvis/FieldPaper/bigbox/mass59/elvis_alldwarfs_v'+vel+'_virial2_interp.dat')
    subhalolist = Table.read('/Users/Sean/research/elvis/FieldPaper/bigbox/inputdata/elvis_alldwarfs_Etracks_distanceinput_v'+vel+'_m59.dat', format = 'ascii')
    totalsublist = Table.read('/Users/Sean/research/elvis/ELVIS_Data_Clean/elvis_alldwarfs_clean_m59.dat', format = 'ascii')

    print len(totalsublist)
    
    subhaloID = np.array(subhalolist['ID'])
    hostlist = np.array(subhalolist['ELVISname'])
    subhalotest = np.array(totalsublist['ID'])

    zinfall = np.array([])
    hostnum = np.array([])
    distance = np.array([])
    velocity1 = np.array([])
    velocity2 = np.array([])
    minvel1 = np.array([])
    minvel2 = np.array([])
    virial1 = np.array([])
    virial2 = np.array([])


    print len(host1data)

    for i in range(len(host1data)):

        #read in host virial radius data 
        hostname = hostlist[i]
        print len(host1data)-i
        Rvir1 = vir1data[i]
        Rvir2 = vir2data[i]

        #select distance data for subhalo, multiplying by 1000 converts units from Mpc to kpc
        dist1_data = host1data[i]*1000 
        dist2_data = host2data[i]*1000
        vel1_data = vel1data[i]
        vel2_data = vel2data[i]
        distance1 = dist1_data[0]
        distance2 = dist2_data[0]
        a_data = scale[i]

        #eliminate data from before subhalo formed (ie a = 0) or only data from late times, (a > 0.25)
        acut = a_data > 0.25
        a = a_data[acut]
        dist1 = dist1_data[acut]
        dist2 = dist2_data[acut]
        vel1 = vel1_data[acut]
        vel2 = vel2_data[acut]
        vir1 = Rvir1[acut]
        vir2 = Rvir2[acut]

        #determine when the subhalo crosses the virial radius, this is the infall redshift
        redcut1 = (dist1 > 0) & (dist1 < vir1) & (vir1 > 0) & (a < 0.862)
        redcut2 = (dist2 > 0) & (dist2 < vir2) & (vir2 > 0) & (a < 0.862)

        z = (1-a)/a
        t = cosmo.lookback_time(z)/u.Gyr
        
        zinfall1_cut = z[redcut1]
        zinfall2_cut = z[redcut2]
        tinfall1_cut = t[redcut1]
        tinfall2_cut = t[redcut2]
        vinfall1_cut = vel1[redcut1]
        vinfall2_cut = vel2[redcut2]
        dist1_cut = dist1[redcut1]
        dist2_cut = dist2[redcut2]
        vir1_cut = vir1[redcut1]
        vir2_cut = vir2[redcut2]

        frac1 = dist1_cut / vir1_cut
        frac2 = dist2_cut / vir2_cut

        check1 = np.diff(np.sign(np.diff(dist1_cut)))
        check2 = np.diff(np.sign(np.diff(dist2_cut)))

        peri1 = frac1[check1 > 0.0]
        peri2 = frac2[check2 > 0.0]
        t1 = tinfall1_cut[check1 > 0.0]
        t2 = tinfall2_cut[check2 > 0.0]
        v1 = vinfall1_cut[check1 > 0.0]
        v2 = vinfall2_cut[check2 > 0.0]

        if len(peri1) != 0:
            p1 = peri1[0]
            m1 = np.min(peri1)
            vp1 = v1[0]
            vm1 = v1[np.where(peri1 == np.min(peri1))][0]
            tp1 = t1[0]
            tm1 = t1[np.where(peri1==np.min(peri1))][0]

            velocity1 = np.append(velocity1, vp1)
            minvel1 = np.append(minvel1, vm1)

        elif len(peri2) != 0:
            p2 = peri2[0]
            m2 = np.min(peri2)
            vp2 = v2[0]
            vm2 = v2[np.where(peri2 == np.min(peri2))][0]
            tp2 = t2[0]
            tm2 = t2[np.where(peri2==np.min(peri2))][0]

            velocity2 = np.append(velocity2, vp2)
            minvel2 = np.append(minvel2, vm2)

        else:
            continue



    return velocity1, minvel1, velocity2, minvel2
    #return pericenter1, time1, dt1, pericenter2, time2, dt2
