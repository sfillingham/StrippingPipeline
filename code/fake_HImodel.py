import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

def work():

    data = Table.read('AB_galaxylist.dat', format = 'ascii')
    name = np.array(data['NAME'])
    Mv = np.array(data['Mv'])
    mstar = 10**((4.83-(Mv))/2.5)
    slope = np.empty(len(name))
    norm = np.empty(len(name))
    maxdist = np.empty(len(name))

    pdata1 = Table.read('/Users/Sean/research/elvis/RamPressure/poppingdata/gas_sham.data/mstar.vs.gas.dat', format = 'ascii')
    z1 = np.array(pdata1['z'])
    cut = z1 == 0.0
    mstar1 = np.array(pdata1['Mstar'])
    mHI1 = np.array(pdata1['HI'])
    mstar1 = mstar1[cut]
    mHI1 = mHI1[cut]

    


    for i in range(len(name)):
        print i

        dwarfdata = Table.read('data/'+np.str(name[i])+'_clean.dat', format = 'ascii')
        n = np.array(dwarfdata['N_HI'])
        r = np.array(dwarfdata['R(kpc)'])

        N = n[n<50]
        R = r[n<50]

        fit = np.polyfit(R,N,1)

        slope[i] = fit[0]
        norm[i] = fit[1]
        maxdist[i] = np.max(R)
        print maxdist[i]

    tdata = Table.read('THINGS_galaxylist.dat', format = 'ascii')
    tname = np.array(tdata['NAME'])
    tmstar = np.array(tdata['Mstar'])
    tslope = np.empty(len(tname))
    tnorm = np.empty(len(tname))
    tmaxdist = np.empty(len(tname))

    for j in range(len(tname)):
        print j

        dwarfdata = Table.read('data/'+np.str(tname[j])+'_clean.dat', format = 'ascii')
        n = np.array(dwarfdata['N_HI'])
        r = np.array(dwarfdata['R(kpc)'])

        N = n[n<50]
        R = r[n<50]

        fit = np.polyfit(R,N,1)

        tslope[j] = fit[0]
        tnorm[j] = fit[1]
        tmaxdist[j] = np.max(R)
        print tmaxdist[j]

    sdata = Table.read('SHIELD_list.dat', format = 'ascii')
    sname = np.array(sdata['NAME'])
    smstar = np.array(sdata['Mstar'])
    sslope = np.empty(len(sname))
    snorm = np.empty(len(sname))
    smaxdist = np.empty(len(sname))

    for k in range(len(sname)):
        print k

        dwarfdata = Table.read('data/'+np.str(sname[k])+'_clean.dat', format = 'ascii')
        n = np.array(dwarfdata['N_HI'])
        r = np.array(dwarfdata['R(kpc)'])

        N = n[n<50]
        R = r[n<50]

        fit = np.polyfit(R,N,1)

        sslope[k] = fit[0]
        snorm[k] = fit[1]
        smaxdist[k] = np.max(R)
        print smaxdist[k]

    slope = np.append(slope,tslope)
    slope = np.append(slope,sslope)
    
    norm = np.append(norm,tnorm)
    norm = np.append(norm,snorm)
    print np.mean(norm)
    print np.median(norm)

    mstar = np.append(mstar,10**tmstar)
    mstar = np.append(mstar,10**smstar)

    maxdist = np.append(maxdist, tmaxdist)
    maxdist = np.append(maxdist, smaxdist)

    #######################################################################
    slopefit = np.polyfit(np.log10(mstar),slope,2)
    #normfit = np.polyfit(np.log10(mstar),norm,1)
    Rdistfit = np.polyfit(np.log10(mstar), maxdist, 1)

    mstarfit = np.polyfit(mstar1, mHI1,1)

    x = np.linspace(5,11,50)
    Rmax = np.empty(len(x))
    y = np.empty(len(x))
    z = np.empty(len(x))
    for i in range(len(x)):
        if x[i] <= 10.0:
            y[i] = slopefit[0]*x[i]**2 + slopefit[1]*x[i] + slopefit[2]
        else:
            y[i] = 0.0

    #z = -7*normfit[0]*x + normfit[1]+0.35
    Rmax = Rdistfit[0]*x + Rdistfit[1]
    HImass = mstarfit[0]*x + mstarfit[1]
    #get HImass from Popping MstarvsMHI relation
    z = HImass - np.log10(2*np.pi) - (y*Rmax**2)/2 - Rmax*np.log10(Rmax) + Rmax

    plt.figure(4)
    plt.plot(mstar1, mHI1, 'ko')
    plt.plot(x, HImass, 'b-')

    plt.figure(5)
    plt.plot(np.log10(mstar), maxdist, 'ko')
    plt.plot(x, Rmax, 'b-')

            
    #make model galaxy names
    numofgal = 50
    fakename = np.chararray(numofgal, itemsize = 5)
    for i in range(numofgal):
        fakename[i] = 'gal'+np.str(i)
        star = x[i]

        #rmax = 3.333*star - 11.667
        rmax = 23
        r = np.linspace(0,rmax,100)
        n = y[np.where(x == star)[0]]*r + z[np.where(x == star)[0]]

        if np.min(n) < 16:
            r = r[n>16]
            n = n[n>16]
        else:
            r=r
            n=n

        modeltable = Table([n,r], names = ('N_HI', 'R(kpc)'), meta={'HI Surface Density': 'Model_HI'})
        outputfile = 'data/'+fakename[i]+'_clean.dat'
        modeltable.write(outputfile, format = 'ascii')

        plt.figure(3)
        plt.plot(r,n,'b-')

    #plt.savefig('model_HI.pdf')
        

    print fakename
    
    inputtable = Table([fakename, x], names=('NAME', 'Mstar'), meta={'HI Surface Density': 'Model_HI'})
    inputfile = 'model_input.dat'
    inputtable.write(inputfile, format = 'ascii')


    plt.figure(1)
    plt.plot(np.log10(mstar),slope,'ko')
    plt.plot(x,y,'g-', linewidth = 2.0)

    #plt.savefig('model_slope.pdf')

    norm_mean = np.empty(len(norm))
    norm_mean[:] = np.mean(norm)
    
    plt.figure(2)
    plt.plot(np.log10(mstar),norm,'ko')
    plt.plot(x,z, 'b-', linewidth = 2.0)

    #plt.savefig('model_centraldensity.pdf')

    plt.show()

    
