#This script take the interpolated output files from the ELVIS suite for a given mass range and determines the distribution of satellite
#pericenters relative to the host in which they reside today.
#
#
import numpy as np
from astropy.table import Table
from astropy.cosmology import LambdaCDM
import astropy.units as u
cosmo = LambdaCDM(H0=71,Om0=0.266,Ode0=0.734)

def work(zcut):

    host1data = np.loadtxt('/Users/Sean/research/elvis/FieldPaper/bigbox/mass59/elvis_alldwarfs_v0.3_distance_time_interp_host1.dat')
    host2data = np.loadtxt('/Users/Sean/research/elvis/FieldPaper/bigbox/mass59/elvis_alldwarfs_v0.3_distance_time_interp_host2.dat')
    scale = np.loadtxt('/Users/Sean/research/elvis/FieldPaper/bigbox/mass59/elvis_alldwarfs_v0.3_scalefactor_interp.dat')
    vir1data = np.loadtxt('/Users/Sean/research/elvis/FieldPaper/bigbox/mass59/elvis_alldwarfs_v0.3_virial1_interp.dat')
    vir2data = np.loadtxt('/Users/Sean/research/elvis/FieldPaper/bigbox/mass59/elvis_alldwarfs_v0.3_virial2_interp.dat')
    subhalolist = Table.read('/Users/Sean/research/elvis/FieldPaper/bigbox/inputdata/elvis_alldwarfs_Etracks_distanceinput_v0.3_m59.dat', format = 'ascii')
    totalsublist = Table.read('/Users/Sean/research/elvis/ELVIS_Data_Clean/elvis_alldwarfs_clean_m59.dat', format = 'ascii')
    
    subhaloID = np.array(subhalolist['ID'])
    hostlist = np.array(subhalolist['ELVISname'])
    subhalotest = np.array(totalsublist['ID'])

    satID1 = np.array([])
    satID2 = np.array([])
    time1 = np.array([])
    time2 = np.array([])
    dt1 = np.array([])
    dt2 = np.array([])
    tperi = np.array([])
    tdiff = np.array([])
    distance = np.array([])
    pericenter1 = np.array([])
    pericenter2 = np.array([])
    min1 = np.array([])
    min2 = np.array([])
    mtime1 = np.array([])
    mtime2 = np.array([])
    mdt1 = np.array([])
    mdt2 = np.array([])
    t50_1 = np.array([])
    t50_2 = np.array([])

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
        distance1 = dist1_data[0]
        distance2 = dist2_data[0]
        a_data = scale[i]

        #eliminate data from before subhalo formed (ie a = 0) or only data from late times, (a > 0.25)
        acut = a_data > 0.25
        a = a_data[acut]
        dist1 = dist1_data[acut]
        dist2 = dist2_data[acut]
        vir1 = Rvir1[acut]
        vir2 = Rvir2[acut]

        #determine when the subhalo crosses the virial radius
        #this is the redshift which corresponds to when the subhalo crosses the chosen distance
        acut = 1./(1+zcut)
        redcut1 = (dist1 > 0) & (dist1 < vir1) & (vir1 > 0) & (a < acut)
        redcut2 = (dist2 > 0) & (dist2 < vir2) & (vir2 > 0) & (a < acut)

        z = (1-a)/a
        t = cosmo.lookback_time(z)/u.Gyr
        
        zinfall1_cut = z[redcut1]
        zinfall2_cut = z[redcut2]
        tinfall1_cut = t[redcut1]
        tinfall2_cut = t[redcut2]
        dist1_cut = dist1[redcut1]
        dist2_cut = dist2[redcut2]
        vir1_cut = vir1[redcut1]
        vir2_cut = vir2[redcut2]

        frac1 = dist1_cut / vir1_cut
        frac2 = dist2_cut / vir2_cut

        #determine when the subhalo crosses some fraction of the virial radius
        rfrac = 0.5
        fraccut1 = (frac1 <= rfrac)
        fraccut2 = (frac2 <= rfrac)
        timeR50_1 = tinfall1_cut[fraccut1]
        timeR50_2 = tinfall2_cut[fraccut2]

        check1 = np.diff(np.sign(np.diff(dist1_cut)))
        check2 = np.diff(np.sign(np.diff(dist2_cut)))
        #print check1

        peri1 = frac1[check1 > 0.0]
        peri2 = frac2[check2 > 0.0]
        t1 = tinfall1_cut[check1 > 0.0]
        t2 = tinfall2_cut[check2 > 0.0]

        if len(peri1) != 0:
            p1 = peri1[0]
            #print p1
            m1 = np.min(peri1)
            #print m1
            tp1 = t1[0]
            #print tp1
            tm1 = t1[np.where(peri1==np.min(peri1))][0]
            #print tm1
            diff1 = np.max(tinfall1_cut) - tp1
            #print diff1
            mdiff1 = np.max(tinfall1_cut) - tm1
            #print mdiff1
            
            satID1 = np.append(satID1, subhaloID[i])
            pericenter1 = np.append(pericenter1, p1)
            time1 = np.append(time1, tp1)
            dt1 = np.append(dt1, diff1)
            min1 = np.append(min1, m1)
            mtime1 = np.append(mtime1, tm1)
            mdt1 = np.append(mdt1, mdiff1)
            if len(timeR50_1) > 0:
                time501 = np.max(tinfall1_cut) - np.max(timeR50_1)
                t50_1 = np.append(t50_1, time501)
            else:
                t50_1 = np.append(t50_1, -1)


        elif len(peri2) != 0:
            p2 = peri2[0]
            m2 = np.min(peri2)
            tp2 = t2[0]
            tm2 = t2[np.where(peri2==np.min(peri2))][0]
            diff2 = np.max(tinfall2_cut) - tp2
            mdiff2 = np.max(tinfall2_cut) - tm2

            satID2 = np.append(satID2, subhaloID[i])
            pericenter2 = np.append(pericenter2, p2)
            time2 = np.append(time2, tp2)
            dt2 = np.append(dt2, diff2)
            min2 = np.append(min2, m2)
            mtime2 = np.append(mtime2, tm2)
            mdt2 = np.append(mdt2, mdiff2)
            if len(timeR50_2) > 0:
                time502 = np.max(tinfall2_cut) - np.max(timeR50_2)
                t50_2 = np.append(t50_2, time502)
            else:
                t50_2 = np.append(t50_2, -1)

        else:
            continue

    satID = np.append(satID1,satID2)
    p = np.append(pericenter1,pericenter2)
    m = np.append(min1,min2)
    tp = np.append(time1,time2)
    tm = np.append(mtime1,mtime2)
    diff = np.append(dt1,dt2)
    mdiff = np.append(mdt1,mdt2)
    t50 = np.append(t50_1,t50_2)

    mastertable = Table([satID, p, m, tp, tm, diff, mdiff, t50], names=('ID', 'first_peri', 'min_peri', 'first_peri_lookback', 'min_peri_lookback', 'first_peri_tau', 'min_peri_tau', 'r50_tau'), meta={'orbital properties': 'ELVIS'})
        
    masteroutputfile = 'orbitalproperties_ELVIS_2Gyrbound.dat'
    mastertable.write(masteroutputfile, format = 'ascii')

    #return min1, mtime1, mdt1, min2, mtime2, mdt2
    #return t50_1, t50_2
    #return pericenter1, time1, dt1, pericenter2, time2, dt2




