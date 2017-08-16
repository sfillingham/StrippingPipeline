        if ((len(zinfall1_cut) !=0) & (len(zinfall2_cut) !=0)):
            zinfall1 = np.max(zinfall1_cut)
            zinfall2 = np.max(zinfall2_cut)
            vinfall1 = vinfall1_cut[np.where(zinfall1_cut == zinfall1)]
            vinfall2 = vinfall2_cut[np.where(zinfall2_cut == zinfall2)]
            redin1 = zinfall1_cut[np.where(dist1_cut == np.min(dist1_cut))]
            redin2 = zinfall2_cut[np.where(dist2_cut == np.min(dist2_cut))]

        elif ((len(zinfall1_cut) ==0) & (len(zinfall2_cut) !=0)):
            zinfall1 = -1
            zinfall2 = np.max(zinfall2_cut)
            vinfall1 = -1
            vinfall2 = vinfall2_cut[np.where(zinfall2_cut == zinfall2)]
            redin1 = -1
            redin2 = zinfall2_cut[np.where(dist2_cut == np.min(dist2_cut))]

        elif ((len(zinfall1_cut) !=0) & (len(zinfall2_cut) ==0)):
            zinfall1 = np.max(zinfall1_cut)
            zinfall2 = -1
            vinfall1 = vinfall1_cut[np.where(zinfall1_cut == zinfall1)]
            vinfall2 = -1
            redin1 = zinfall1_cut[np.where(dist1_cut == np.min(dist1_cut))]
            redin2 = -1

        else:
            zinfall1 = -1
            zinfall2 = -1
            vinfall1 = -1
            vinfall2 = -1
            redin1 = -1
            redin2 = -1

        velocity1 = np.append(velocity1, vinfall1)
        velocity2 = np.append(velocity2, vinfall2)
        #red1 = np.append(red1, redin1)
        #red2 = np.append(red2, redin2)

        output = Table([velocity1, velocity2], names=('vel_host1', 'vel_host2'), meta={'Velocity': 'Distribution'})
        outputfile = 'velocity_distribution_virial.dat'
        output.write(outputfile, format = 'ascii')
