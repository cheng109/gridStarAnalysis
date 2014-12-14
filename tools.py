__author__ = 'cheng109'

import star_list
import sys
import numpy as np
import matplotlib.pyplot as plt
import astrometric_error as astro_error


def poly_correct_ellipticity(mergedDict,freq):
    e1 = []
    e2 = []
    X =[]
    Y =[]
    newDict = {}
    for (RA_grid, DEC_grid), star in mergedDict["both"].items():
        if np.remainder(RA_grid, freq)==0 and np.remainder(DEC_grid, freq)==0:
            newDict[(RA_grid,DEC_grid)]=star
            e1.append(star.sex_e1)
            e2.append(star.sex_e2)
            X.append(star.sex_RA)
            Y.append(star.sex_DEC)
    num=len(X)
    X=np.array(X)
    Y=np.array(Y)
    v = np.array([np.ones(num), X, Y, X*X, X*Y, Y*Y, X*X*X, X*X*Y,X*Y*Y, Y*Y*Y, X*X*X*X, X*X*X*Y, X*X*Y*Y,X*Y*Y*Y,Y*Y*Y*Y])
    coX, resid_1, rank_1, singval_1 = np.linalg.lstsq(v.T, e1)
    coY, resid_2, rank_2, singval_2 = np.linalg.lstsq(v.T, e2)
    return coX, coY


def get_correct_sex(star, coX, coY):
    X = star.sex_RA
    Y = star.sex_DEC
    star.sex_correct_e1 = star.sex_e1 - np.dot(coX, [1, X, Y, X*X, X*Y, Y*Y, X*X*X, X*X*Y,X*Y*Y, Y*Y*Y, X*X*X*X, X*X*X*Y, X*X*Y*Y,X*Y*Y*Y,Y*Y*Y*Y])
    star.sex_correct_e2 = star.sex_e2 - np.dot(coY, [1, X, Y, X*X, X*Y, Y*Y, X*X*X, X*X*Y,X*Y*Y, Y*Y*Y, X*X*X*X, X*X*X*Y, X*X*Y*Y,X*Y*Y*Y,Y*Y*Y*Y])
    return star


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


def createDistanceList(mergedDict, direction, freq=4, polyCorrect=True):
    cov = []
    unit = []
    if polyCorrect==True:
        coX, coY= poly_correct_ellipticity(mergedDict, freq)

    if direction=='ra':
        for RA_grid, starList in mergedDict['ra'].items():
            if np.remainder(RA_grid, freq)==0:
                for i in range(len(starList)-1):
                    for j in range(i+1, len(starList)):
                        A = starList[i]
                        B = starList[j]
                        if polyCorrect==True:
                            new_A = get_correct_sex(A, coX, coY)
                            new_B = get_correct_sex(B, coX, coY)
                            A_e1 = new_A.sex_correct_e1
                            A_e2 = new_A.sex_correct_e2
                            B_e1 = new_B.sex_correct_e1
                            B_e2 = new_B.sex_correct_e2
                        else:
                            A_e1 = A.sex_correct_e1
                            A_e2 = A.sex_correct_e2
                            B_e1 = B.sex_correct_e1
                            B_e2 = B.sex_correct_e2
                        cov.append(A_e1*B_e1 + A_e2*B_e2)
                        unit.append(int(abs(B.DEC_grid - A.DEC_grid)))
    count =0
    if direction=='dec':
        for DEC_grid, starList in mergedDict['dec'].items():
            if np.remainder(DEC_grid, freq)==0:
                for i in range(len(starList)-1):
                    for j in range(i+1, len(starList)):
                        A = starList[i]
                        B = starList[j]
                        if polyCorrect==True:
                            new_A = get_correct_sex(A, coX, coY)
                            new_B = get_correct_sex(B, coX, coY)
                            A_e1 = new_A.sex_correct_e1
                            A_e2 = new_A.sex_correct_e2
                            B_e1 = new_B.sex_correct_e1
                            B_e2 = new_B.sex_correct_e2
                        else:
                            A_e1 = A.sex_correct_e1
                            A_e2 = A.sex_correct_e2
                            B_e1 = B.sex_correct_e1
                            B_e2 = B.sex_correct_e2
                        cov.append(A_e1*B_e1 + A_e2*B_e2)
                        unit.append(int(abs(B.RA_grid - A.RA_grid)))


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
    x,cov = createDistanceList(mergedDict, direction="correct",freq=freq, polyCorrect=True)
    x[:]=[a*spacing*60.0/4 for a in x]   # convert the unit to be arcmin.
    plotCorrelation(x, cov)
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


