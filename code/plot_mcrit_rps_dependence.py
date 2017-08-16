import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

def mcrit():
    
    nhalo_list = np.array([4.5,4.0,3.5,3.25,3.0])#,2.75,2.5])
    nhalo_list1 = np.array([4.0,3.5,3.25])
    nhalo_list2 = np.array([3.5,3.25,3.0])
    color_list = np.array(['blue','green','red','magenta','cyan'])
    vsat_list = np.array([250,300,350,400,450,500,550,600,650,700,800,900,1000])
    vsat_list1 = np.array([200,250,300,350,400,450,500])#,550,600])
    vsat_list2 = np.array([400,450,500,550,600,650,700,800,900,1000])

    rps = np.empty(len(nhalo_list)*len(vsat_list))
    mcrit = np.empty(len(nhalo_list)*len(vsat_list))

    for i in range(len(nhalo_list)):

        nhalo = nhalo_list[i]

        for j in range(len(vsat_list)):

            vsat = vsat_list[j]

            k = i*len(vsat_list)+j
            rps[k] = 10**(-1*nhalo)*vsat**2
            print rps[k]

            data = Table.read('modeldata/RPS_strippedfrac_hot'+np.str(np.int(10*nhalo))+'_v'+np.str(np.int(vsat))+'_Model.dat', format = 'ascii')
            NFWstrip = np.array(data['NFW'])
            mstar = np.array(data['Mstar'])

            cut = NFWstrip > 0.05
            if len(mstar[cut] > 0):
                mcrit[k] = np.max(mstar[cut])
            else:
                mcrit[k] = 0
            print mcrit[k]

    #### Local Group-like objects ##################
    MWrps = np.empty(len(nhalo_list1)*len(vsat_list1))
    MWmcrit = np.empty(len(nhalo_list1)*len(vsat_list1))

    for i in range(len(nhalo_list1)):

        nhalo = nhalo_list1[i]

        for j in range(len(vsat_list1)):

            vsat = vsat_list1[j]

            k = i*len(vsat_list1)+j
            MWrps[k] = 10**(-1*nhalo)*vsat**2
            print MWrps[k]

            data = Table.read('modeldata/RPS_strippedfrac_hot'+np.str(np.int(10*nhalo))+'_v'+np.str(np.int(vsat))+'_Model.dat', format = 'ascii')
            NFWstrip = np.array(data['NFW'])
            mstar = np.array(data['Mstar'])

            cut = NFWstrip > 0.05
            if len(mstar[cut] > 0):
                MWmcrit[k] = np.max(mstar[cut])
            else:
                MWcrit[k] = 0
            print MWmcrit[k]

    #### Cluster-like objects ##################
    Crps = np.empty(len(nhalo_list2)*len(vsat_list2))
    Cmcrit = np.empty(len(nhalo_list2)*len(vsat_list2))

    for i in range(len(nhalo_list2)):

        nhalo = nhalo_list2[i]

        for j in range(len(vsat_list2)):

            vsat = vsat_list2[j]

            k = i*len(vsat_list2)+j
            Crps[k] = 10**(-1*nhalo)*vsat**2
            print Crps[k]

            data = Table.read('modeldata/RPS_strippedfrac_hot'+np.str(np.int(10*nhalo))+'_v'+np.str(np.int(vsat))+'_Model.dat', format = 'ascii')
            NFWstrip = np.array(data['NFW'])
            mstar = np.array(data['Mstar'])

            cut = NFWstrip > 0.05
            if len(mstar[cut] > 0):
                Cmcrit[k] = np.max(mstar[cut])
            else:
                Cmcrit[k] = 0
            print Cmcrit[k]

    finalRPS = np.append(MWrps,Crps)
    finalMcrit = np.append(MWmcrit,Cmcrit)

    plt.figure(1)
    plt.axis([0,3.5,8,11.5])
    plt.plot(np.log10(rps),np.log10(mcrit), marker = 'o', markersize = 10+i, color = 'cyan', linestyle = 'None')
    #plt.xlabel(r'$\rm log({\it n}_{halo}{\it V}_{sat}^{\ 2})\ (cm^{-3}\ (km\ s^{-1})^{2})$')
    plt.xlabel(r'$\rm log_{10}(P_{ram})\ [cm^{-3}\ (km\ s^{-1})^{2}]$')
    plt.ylabel(r'$\rm log_{10}(M_{crit})\ [M_{\odot}]$')

    plt.figure(2)
    plt.axis([0,3.5,8,11.5])
    plt.plot(np.log10(MWrps),np.log10(MWmcrit), marker = 'o', markersize = 10+i, color = 'magenta', linestyle = 'None')
    #plt.xlabel(r'$\rm log({\it n}_{halo}{\it V}_{sat}^{\ 2})\ (cm^{-3}\ (km\ s^{-1})^{2})$')
    plt.xlabel(r'$\rm log_{10}(P_{ram})\ [cm^{-3}\ (km\ s^{-1})^{2}]$')
    plt.ylabel(r'$\rm log_{10}(M_{crit})\ [M_{\odot}]$')

    plt.figure(3)
    plt.axis([0,3.5,8,11.5])
    plt.plot(np.log10(Crps),np.log10(Cmcrit), marker = 'o', markersize = 10+i, color = 'Goldenrod', linestyle = 'None')
    #plt.xlabel(r'$\rm log({\it n}_{halo}{\it V}_{sat}^{\ 2})\ (cm^{-3}\ (km\ s^{-1})^{2})$')
    plt.xlabel(r'$\rm log_{10}(P_{ram})\ [cm^{-3}\ (km\ s^{-1})^{2}]$')
    plt.ylabel(r'$\rm log_{10}(M_{crit})\ [M_{\odot}]$')

    plt.figure(4)
    plt.axis([0,3.5,8,11.5])
    plt.plot(np.log10(finalRPS),np.log10(finalMcrit), marker = 'o', markersize = 10+i, color = 'blue', linestyle = 'None')
    #plt.xlabel(r'$\rm log({\it n}_{halo}{\it V}_{sat}^{\ 2})\ (cm^{-3}\ (km\ s^{-1})^{2})$')
    plt.xlabel(r'$\rm log_{10}(P_{ram})\ [cm^{-3}\ (km\ s^{-1})^{2}]$')
    plt.ylabel(r'$\rm log_{10}(M_{crit})\ [M_{\odot}]$')

    return rps, mcrit
            

    
