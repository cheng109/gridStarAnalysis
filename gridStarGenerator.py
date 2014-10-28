"""
  @package 
  @file gridCatalog.py
  @brief generate stars on square grid
 
  @brief Created by:
  @author Jun Cheng (Purdue)
 
  @warning This code is not fully validated
  and not ready for full release.  Please
  treat results with caution.

Usage: python gridCatalog.py
 
Notes: 1. You could set the 'magnitude', 'density'

"""
import sys,os
import numpy as np
import matplotlib.pyplot as plt
import math
#from random import gauss
header ="""Unrefracted_Azimuth 0
Unrefracted_Altitude 89
Slalib_date 1994/7/19/0.298822999997
Opsim_rotskypos 0
Opsim_rottelpos 0
Opsim_moondec -90
Opsim_moonra 180
Opsim_expmjd 49552.3
Opsim_moonalt -90
Opsim_sunalt -90
Opsim_filter 2
Opsim_dist2moon 180.0
Opsim_moonphase 10.0
Opsim_obshistid 99999999
Opsim_rawseeing 0.67
SIM_SEED     1000
SIM_MINSOURCE 1
SIM_TELCONFIG 0
SIM_CAMCONFIG 1
SIM_VISTIME 15.0
SIM_NSNAP 1
"""

def getGrid(minRA, maxRA, minDEC, maxDEC, deltaRA, deltaDEC):
    RA=[]
    DEC=[]
    rangeRA=int(math.floor((maxRA-minRA)/deltaRA))
    rangeDEC=int(math.floor((maxDEC-minDEC)/deltaDEC))
    num = rangeRA*rangeDEC     # the total number of stars
    for i in range(rangeRA):
        RA.append(minRA+i*deltaRA)
    for i in range(rangeDEC):
        DEC.append(minDEC+i*deltaDEC)
    position=[]
    for ra in RA:
        for dec in DEC:
            position.append((ra, dec))
    return position

def writeCatalog(catalogName,header, mag, position):
    f = open(catalogName,"w")
    f.write(header)
    counter = 0
    for ra, dec in position:
        f.write('object ' + str(counter) +' '+ str(position[0]) +' '+ str(position[1]) + ' ' + str(mag[counter]))
        f.write(' ../sky/sed_flat.txt 0 0 0 0 0 0 star none none\n')
        counter+=1
    f.close()

def genMagnitudeList(num,mu,sigma):
    mag =[]
    for i in range(num):
        mag.append(gauss(mu,sigma))
    return mag

def plotCheck(position):
    for (x, y) in position:
        plt.plot(x,y, '*')
    plt.xlabel("RA (degree)")
    plt.ylabel("DEC (degree)")
    plt.savefig("foo.png")
    #plt.show()

if __name__=='__main__':
    pointing_RA = 180
    pointing_DEC = 0
    header = "Unrefracted_RA_deg " +str(pointing_RA)+ "\nUnrefracted_Dec_deg " + str(pointing_DEC) +"\n" + header

    position = getGrid(180-2.5,180+2.5, 0-2.5, 0+2.5, 0.003,0.003)
    plotCheck(position)
    #magnitude = genMagnitudeList(num, 18,3.0)
    #magnitude = np.random.uniform(18,20,num)
    magnitude = [18]*len(position)
    writeCatalog("gridCatalog",header, magnitude, position)














