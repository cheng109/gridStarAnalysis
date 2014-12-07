__author__ = 'cheng109'


import sys
import star_list
import tools
import numpy as np
import matplotlib.pylab as plt
import time, threading
import pickle

def do_compareSEX(starListDict):
    i=0
    index =[]
    error = []
    mine = []
    sex = []
    for key, starList in starListDict.items():
        for star in starList:
            i=i+1
            index.append(i)
            mine.append(1-star.mine_ellipticity)
            sex.append(star.sex_ellipticity)
            error.append(mine[-1]-sex[-1])
    print "sex", sex
    print "mine",mine
    print "error",error
    plt.plot(index, mine, 'g*-', index, sex,'y*-', index, error,'b*-')

   # plt.hist(error)
    plt.show()



def writeStarInfo(starListDict, outFile):
    f = open(outFile, 'w')
    for key, starList in starListDict.items():
        for star in starList:
            info= star.star_info()
            for i in info:
                f.write(i +"\t")
            f.write("\n")
    f.close()




if __name__=='__main__':
    starListDict ={}
    threadList = []
    spacing = 0.003*4    # unit degree
    freq = 1
    winX=80
    winY=80
    xlim = [100,16000-100]
    ylim = [100,16288-100]

    #xlim = [10,4000-10]
    #ylim = [10,4000-10]

    DIR="/Volumes/HD3/data/CCD_ATM_OPTICS/ccd_OFF_atm_ON_perturbation/"
    #DIR="/Users/cheng109/work/source/gridStarAnalysis/gridstaranalysis/data/defect_ON_atm_ON/"
    start_time=time.time()
    fileList = tools.genFileNameList("OFF_ON_perturbation_fileList.txt")

    starListDict =star_list.generateStarListDictionary(fileList, (winX, winY), (xlim, ylim),
                                                       sex=True,
                                                       Moffat=False,
                                                       Ellipticity=False,
                                                       dir=DIR) #with ending '/'


    tools.do_correlatin(starListDict,DIR+"correlation.txt", spacing, freq=freq)
    writeStarInfo(starListDict, DIR+"star-info-out.txt")
    print "Initial Time usage: " + str(time.time()-start_time)
