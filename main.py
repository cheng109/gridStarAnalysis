__author__ = 'cheng109'


import sys
import star_list
import cor
import e_map
import matplotlib.pylab as plt
import plotLib

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




if __name__=='__main__':
    spacing = 0.003    # unit degree
    freq = 4
    winX=12
    winY=12
    xlim = [20,4000-20]
    ylim = [20,4000-20]
    starListDict = star_list.generateStarListDictionary(sys.argv[1:], (winX, winY), (xlim, ylim),
                                                        sex=False,
                                                        Moffat=False,
                                                        Ellipticity=False,
                                                        dir="data/defect_OFF_atm_ON/") #with ending '/'
    #cor.do_writeOutAllStars(starListDict, "star_info_OFF_ON.txt")

    #cor.do_correlatin(starListDict, spacing, freq=freq)
    #do_compareSEX(starListDict)

    plotLib.plotEllipticityMap(starListDict,spacing,freq=6)
    print "completely done ! "