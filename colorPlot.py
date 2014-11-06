__author__ = 'cheng109'


import matplotlib.pylab as plt
import numpy as np
def colorPlot():
    RA=[]
    DEC=[]
    FWHM=[]
    e1=[]
    e2=[]
    e=[]
    FWHM_max=0
    FWHM_min=1
    f=open("star_info.txt")
    i=0
    for line in f.readlines():
        temp= line.split()
        RA.append(float(temp[0]))
        DEC.append(float(temp[1]))
        e1.append(float(temp[2]))
        e2.append(float(temp[3]))
        e.append(np.sqrt(e1[-1]**2+e2[-1]**2))
        FWHM.append(float(temp[6]))
        if e[-1]>FWHM_max:
            FWHM_max=e[-1]
        if e[-1]<FWHM_min:
            FWHM_min=e[-1]
        i=i+1
       # if i>3000:
       #     break


    for i in range(len(RA)):
        plt.scatter(RA[i],DEC[i], color=str(float(e[i]-FWHM_min)/FWHM_max))
     #   if i>3000:
     #       break

    plt.show()


colorPlot()
print "haha"