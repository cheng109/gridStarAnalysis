__author__ = 'cheng109'

import star_list
import sys
import numpy as np
import matplotlib.pyplot as plt
import astrometric_error as astro_error

def update_grid(starListDict, spacing):
    firstStarList = starListDict.values()[0]
    temp_RA = firstStarList[0].sex_RA
    temp_DEC = firstStarList[0].sex_DEC

    for key, starList in starListDict.items():
        for i in range(len(starList)):
            starList[i].RA_grid=int(np.round((starList[i].sex_RA-temp_RA)/spacing))
            starList[i].DEC_grid=int(np.round((starList[i].sex_DEC-temp_DEC)/spacing))
            #plt.plot(starList[i].RA_grid, starList[i].DEC_grid, '*b')

def mergeAllStarList(starListDict):
    mergedDict_RA = {}
    mergedDict_DEC = {}
    newDict = {}
    for key, starList in starListDict.items():
        for star in starList:
            if star.RA_grid not in mergedDict_RA:
                mergedDict_RA[star.RA_grid]=[]
            if star.DEC_grid not in mergedDict_DEC:
                mergedDict_DEC[star.DEC_grid]=[]
            mergedDict_RA[star.RA_grid].append(star)
            mergedDict_DEC[star.DEC_grid].append(star)
            newDict[(star.RA_grid, star.DEC_grid)] = star


    mergedDict = {}
    mergedDict['ra']=mergedDict_RA
    mergedDict['dec']=mergedDict_DEC
    mergedDict['both']=newDict
    return mergedDict


def createDistanceList(indexDict, direction, freq=4):

    cov = []
    unit = []
    if direction=='ra':
        for RA_grid, starList in indexDict.items():
            if np.remainder(RA_grid, freq)==0:
                for i in range(len(starList)-1):
                    for j in range(i+1, len(starList)):
                        A = starList[i]
                        B = starList[j]
                        cov.append(A.sex_e1*B.sex_e1 + A.sex_e2*B.sex_e2)
                        unit.append(int(abs(B.DEC_grid - A.DEC_grid)))
    if direction=='dec':
        for DEC_grid, starList in indexDict.items():
            if np.remainder(DEC_grid, freq)==0:
                for i in range(len(starList)-1):
                    for j in range(i+1, len(starList)):
                        A = starList[i]
                        B = starList[j]
                        cov.append(A.sex_e1*B.sex_e1 + A.sex_e2*B.sex_e2)
                        unit.append(int(abs(B.RA_grid - A.RA_grid)))


    x=list(set(unit))
    mean = [0]*len(x)
    d = [[] for i in range(len(x))]
    for i in range(len(unit)):
        d[unit[i]-1].append(cov[i])

    for j in range(len(d)):
        mean[j] = np.mean(d[j])

    return x, mean

def plotCorrelation(x, cov):

    plt.plot(x,cov, '*-')
    plt.xlabel(r"$\theta$ (arcmin)")
    plt.ylabel("Correlation")
    plt.show()

def writeToFile(x, cov, outFileName):
    f=open(outFileName, 'w')
    for i in range(len(x)):
        f.write(str(x[i])+"\t"+str(cov[i])+"\n")
    f.close()


def do_writeOutAllStars(starListDict, outFileName):
    f=open(outFileName, 'w')
    for key, value in starListDict.items():
        for star in value:
            f.write(str(star.sex_RA) + "\t"
            +str(star.sex_DEC) + "\t"
            #+str(star.psfSize) + "\t"
            #+str(star.e1) + "\t"
            #+str(star.e2) +"\t"
            +str(star.sex_ellipticity) + "\t"
            #+str(star.moffat_smallest_fwhm) + "\t"
            #+str(star.sex_FWHM) + "\t"
            #+str(star.moffat_largest_fwhm) +"\t"
            + "\n")
    f.close()


def do_correlatin(starListDict, spacing, freq):
    update_grid(starListDict, spacing)
    print "updating grids done ! "
    mergedDict = mergeAllStarList(starListDict,)
    x,mean = createDistanceList(mergedDict, direction="dec",freq=freq)
    x[:]=[a*spacing*60.0 for a in x]   # convert the unit to be arcmin.
    plotCorrelation(x, mean)
    writeToFile(x, mean, "outputFile.txt")


if __name__=="__main__":
    freq = 4
    spacing = 0.003    # unit degree
    winX=10
    winY=10
    xlim = [20,4000-20]
    ylim = [20,4000-20]
    starListDict = star_list.generateStarListDictionary(sys.argv[1:], (winX, winY), (xlim, ylim), sex=False, Moffat=False, dir="data/defect_OFF_atm_ON/") #with ending '/'
    #writeOutAllStars(starListDict, "star_info.txt")

    do_correlatin(starListDict, spacing, freq=freq)

    print "completely done ! "


