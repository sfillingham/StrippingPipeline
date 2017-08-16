        if ((len(zinfall1_cut) !=0) & (len(zinfall2_cut) !=0)):
            zinfall1 = np.amax(zinfall1_cut)
            zinfall2 = np.amax(zinfall2_cut)
            peri1 = np.min(dist1_cut)
            peri2 = np.min(dist2_cut)
            redin1 = zinfall1_cut[np.where(dist1_cut == np.min(dist1_cut))]
            redin2 = zinfall2_cut[np.where(dist2_cut == np.min(dist2_cut))]

        elif ((len(zinfall1_cut) ==0) & (len(zinfall2_cut) !=0)):
            zinfall1 = -1
            zinfall2 = np.amax(zinfall2_cut)
            peri1 = -1
            peri2 = np.min(dist2_cut)
            redin1 = -1
            redin2 = zinfall2_cut[np.where(dist2_cut == np.min(dist2_cut))]

        elif ((len(zinfall1_cut) !=0) & (len(zinfall2_cut) ==0)):
            zinfall1 = np.amax(zinfall1_cut)
            zinfall2 = -1
            peri1 = np.min(dist1_cut)
            peri2 = -1
            redin1 = zinfall1_cut[np.where(dist1_cut == np.min(dist1_cut))]
            redin2 = -1

        else:
            zinfall1 = -1
            zinfall2 = -1
            peri1 = -1
            peri2 = -1
            redin1 = -1
            redin2 = -1

        pericenter1 = np.append(pericenter1, peri1)
        pericenter2 = np.append(pericenter2, peri2)

        output = Table([pericenter1, pericenter2], names=('peri_host1', 'peri_host2'), meta={'Pericenter': 'Distribution'})
        outputfile = 'pericenter_distribution.dat'
        output.write(outputfile, format = 'ascii')

    return pericenter1, pericenter2
