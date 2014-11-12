__author__ = 'cheng109'

import matplotlib.pylab as plt
import numpy as np
from matplotlib.patches import Ellipse
import cor
def plotEllipticityMap(starListDict, spacing, freq=1):
    n = 30  # enlarge buy 'n'
    ells = []

    cor.update_grid(starListDict, spacing)
    mergedDict = cor.mergeAllStarList(starListDict)
    #for key, starList in starListDict.items():

    for (RA_grid, DEC_grid), star in mergedDict["both"].items():
        star.update(Moffat=True,Ellipticity=False)
        if np.remainder(RA_grid, freq)==0 and np.remainder(DEC_grid, freq)==0:
            ells.append(Ellipse([star.sex_RA, star.sex_DEC],
                      width=star.moffat_smallest_fwhm/18000*n,
                      height = star.moffat_largest_fwhm/18000*n,
                      angle = star.moffat_phi    #in degree
                      ))
    fig = plt.figure()

    ax = fig.add_subplot(111, aspect='equal')
    for e in ells:
        ax.add_artist(e)
        #e.set_clip_box(ax.bbox)
        #e.set_alpha(plt.rand())
        #e.set_facecolor(plt.rand(3))
    centerRA = 180
    centerDEC = 0
    side = 0.125
    ax.set_xlim(centerRA-side, centerRA+side*15)
    ax.set_ylim(centerDEC-side, centerDEC+side)
    plt.show()








