"""
  @package 
  @file moffatFitting.py
  @brief compute the fwhm from the 'segmentation.fits' and 'image.fits'
 
  @brief Created by:
  @author Jun Cheng (Purdue)
  
  @warning This code is not fully validated
  and not ready for full release.  Please
  treat results with caution.

Usage: python sex_validation.py <image.fits> <test.cat> <result.txt>
 
Notes:
1. This code is mainly used to validate that central value from sextractor is CORRECT.
2. It takes <image.fits> and <test.cat> from sextractor as input file, and output the results.
3. This code can save all the moffat fitting results of stars throughout the whole image.

"""
import star_list
import numpy as np
import sys
import matplotlib.pyplot as plt

def angle_distance(A, B):
    A[:] = [x/180*np.pi for x in A]
    B[:] = [x/180*np.pi for x in B]
    return np.arccos(np.cos(A[1])*np.cos(B[1])*np.cos(A[0]-B[0]) + np.sin(A[1])*np.sin(B[1]))*180/np.pi*3600000    # in arcsec

def update_grid(starList, spacing):
    temp_CenterX = 10000
    temp_CenterY = 10000
    temp_RA = 0
    temp_DEC = 0

    for i in range(len(starList)):
        if (starList[i].sex_CenterX+starList[i].sex_CenterY)<(temp_CenterX+temp_CenterY):
             temp_CenterX = starList[i].sex_CenterX
             temp_CenterY = starList[i].sex_CenterY
             temp_RA = starList[i].sex_RA
             temp_DEC = starList[i].sex_DEC

    for i in range(len(starList)):
        starList[i].RA_grid=np.round((starList[i].sex_RA-temp_RA)/spacing)
        starList[i].DEC_grid=np.round((starList[i].sex_DEC-temp_DEC)/spacing)
        #plt.plot(starList[i].RA_grid, starList[i].DEC_grid, '*b')
    #plt.show()


def createDistanceList(starList):
    global unit_dist
    distance = [];
    unit = [];
    for i in range(len(starList)-1):
        for j in range(i+1, len(starList)):
            A = starList[i]
            B = starList[j]
            if int(A.RA_grid)==int(B.RA_grid):
                distance.append(angle_distance([A.sex_RA, A.sex_DEC], [B.sex_RA, B.sex_DEC]))
                unit.append(int(abs(B.DEC_grid - A.DEC_grid)))

    x=list(set(unit))
    std = [0]*len(x)
    d = [[] for i in range (len(x))]
    for i in range(len(unit)):
        d[unit[i]-1].append(distance[i])

    for j in range(len(d)):
        std[j] = np.std(d[j])
    return x, std

def plotAstroError(x,std,labelName):

    plt.plot(x,std, '*--', label='lines')
    plt.xlabel("Seperation (arcsec)")
    plt.ylabel("STD (MAS)")
    #patch = mpatches.Patch(color='red', label=labelName)
    #plt.legend(handles=[patch])
    #plt.show()

def writeToFile(x, std, fileName):
    f = open(fileName[:-5]+"_output.txt",'w')
    for i in range(len(x)):
        f.write(str(x[i]) + "\t" + str(std[i])+"\n")
    f.close()

if __name__=="__main__":
    spacing = 0.0077
    winX=10
    winY=10
    xlim = [20,4000-20]
    ylim = [20,4000-20]
    starListDict = star_list.generateStarListDictionary(sys.argv[1:], (winX, winY), (xlim, ylim), sex=True, Moffat=False)
    for key in starListDict:
        update_grid(starListDict[key], spacing)
        x, std = createDistanceList(starListDict[key])
        x=[spacing*i*3600 for i in x]
        writeToFile(x, std, key)

        plotAstroError(x,std, key)
    plt.show()
