#This code will check how close to pericenter

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.cosmology as cosmo

def work():

    host1data = np.loadtxt('/Users/Sean/research/elvis/FieldPaper/bigbox/mass59/elvis_alldwarfs_v0.3_distance_time_interp_host1.dat')
    host2data = np.loadtxt('/Users/Sean/research/elvis/FieldPaper/bigbox/mass59/elvis_alldwarfs_v0.3_distance_time_interp_host2.dat')
    scale = np.loadtxt('/Users/Sean/research/elvis/FieldPaper/bigbox/mass59/elvis_alldwarfs_v0.3_scalefactor_interp.dat')


    for i in range(100):

        dist = host1data[i]*1000
        test = dist[749] < 350.

        if dist[749] < 350:
            plt.figure(1)
            plt.plot(dist)
            plt.show()

        else:
            plt.figure(2)
            plt.plot(dist)
            plt.show()
            

    
