__author__ = 'cheng109'

import star_list
import sys
import numpy as np
import matplotlib.pyplot as plt
import astrometric_error as astro_error


def genFileNameList(listFile):
    fileList = []
    infile = open(listFile,'r')
    for line in infile.readlines():
        if line[0]!='#':
            temp=line.split()
            for fits in temp:
                fileList.append(fits)

    infile.close()
    return fileList


def update_grid(starListDict, spacing):
    firstStarList = starListDict.values()[0]
    temp_RA = firstStarList[0].sex_RA
    temp_DEC = firstStarList[0].sex_DEC

    for key, starList in starListDict.items():
        for i in range(len(starList)):
            starList[i].RA_grid=int(np.round((starList[i].sex_RA-temp_RA)/spacing))
            starList[i].DEC_grid=int(np.round((starList[i].sex_DEC-temp_DEC)/spacing))
            #plt.plot(starList[i].RA_grid, starList[i].DEC_grid, '*b')
    #plt.show()

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
        for RA_grid, starList in indexDict['ra'].items():
            if np.remainder(RA_grid, freq)==0:
                for i in range(len(starList)-1):
                    for j in range(i+1, len(starList)):
                        A = starList[i]
                        B = starList[j]
                        cov.append(A.sex_e1*B.sex_e1 + A.sex_e2*B.sex_e2)
                        unit.append(int(abs(B.DEC_grid - A.DEC_grid)))
    count =0
    if direction=='dec':
        for DEC_grid, starList in indexDict['dec'].items():
            if np.remainder(DEC_grid, freq)==0:
                for i in range(len(starList)-1):
                    for j in range(i+1, len(starList)):
                        A = starList[i]
                        B = starList[j]
                        cov.append(A.sex_e1*B.sex_e1 + A.sex_e2*B.sex_e2)
                        #cov.append(A.mine_e1*B.mine_e1 + A.mine_e2*B.mine_e2)
                        unit.append(int(abs(B.RA_grid - A.RA_grid)))
                        if unit[-1]==191:
                            count +=1
                            print "haha"

    print "count = ", count
    x=sorted(list(set(unit)))

    min_d = min(unit)
    max_d = max(unit)
    d = [[] for i in range(len(x))]
    mean = [0 for i in range(len(x))]
    for i, u in enumerate(unit):
 #       d[unit[i]-min_d].append(cov[i])
       d[x.index(u)].append(cov[i])
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

def do_correlatin(starListDict, outputFile, spacing, freq):
    update_grid(starListDict, spacing)
    print "updating grids done ! "
    mergedDict = mergeAllStarList(starListDict)
    x,cov = createDistanceList(mergedDict, direction="dec",freq=freq)
    x[:]=[a*spacing*60.0/4 for a in x]   # convert the unit to be arcmin.
    #plotCorrelation(x, cov)
    writeToFile(x, cov, outputFile)


if __name__=="__main__":
    freq = 4
    spacing = 0.003    # unit degree
    winX=100
    winY=100
    xlim = [100,40000-100]
    ylim = [100,40000-100]
    starListDict = star_list.generateStarListDictionary(sys.argv[1:], (winX, winY), (xlim, ylim), sex=False, Moffat=False, dir="data/defect_OFF_atm_ON/") #with ending '/'
    #writeOutAllStars(starListDict, "star_info.txt")

    do_correlatin(starListDict, "correlation.txt", spacing, freq=freq)

    print "completely done ! "


