"""
  @package 
  @file star_list.py
  @brief compute the fwhm from the 'segmentation.fits' and 'image.fits'
 
  @brief Created by:
  @author Jun Cheng (Purdue)
  
  @warning This code is not fully validated
  and not ready for full release.  Please
  treat results with caution.

Usage: python sex_validation.py <image.fits> <test.cat> <result.txt>
 
Notes:

1. This is a procedure to create a dictionary {'fitsFileName01': StarList01, 'fitsFileName02': StarList02,......}
2. If you don't have a .cat file from SExtractor, then you can create one by setting 'sex=True'
3. If you want to used Moffat fitted value instead of the one from SExtrator, specify 'Moffat=True'

"""
import sys, subprocess
from astropy.io import fits
from astropy.wcs import WCS
from aspylib import astro
import numpy as np
import matplotlib.pyplot as plt
import time
import threading


class Star:
    def __init__(self, ID, sex_CenterX, sex_CenterY, sex_RA, sex_DEC, data, wcs, origin, windows, e):
        self.winX =windows[0]
        self.winY =windows[1]

        self.ID = ID
        self.sex_CenterX = sex_CenterX
        self.sex_CenterY = sex_CenterY
        self.sex_RA = sex_RA
        self.sex_DEC = sex_DEC
        self.data = data
        self.wcs = wcs
        self.origin = origin

        self.RA_grid = 0
        self.DEC_grid = 0

        self.psfSize=0

        self.sex_X2=e[0]
        self.sex_Y2=e[1]
        self.sex_XY=e[2]

        self.sex_e1 = (self.sex_X2-self.sex_Y2)/(self.sex_X2+self.sex_Y2)
        self.sex_e2 = 2*self.sex_XY/(self.sex_X2+self.sex_Y2)

        self.sex_correct_e1=0
        self.sex_correct_e2=0
        #self.sex_e1 = (self.sex_X2-self.sex_Y2+self.sex_CenterX**2-self.sex_CenterY**2)/(self.sex_X2+self.sex_Y2+self.sex_CenterX**2+self.sex_CenterY**2)
        #self.sex_e2 = 2*(self.sex_XY+self.sex_CenterY*self.sex_CenterX)/(self.sex_X2+self.sex_Y2+self.sex_CenterX**2+self.sex_CenterY**2)

        self.sex_ellipticity= np.sqrt(self.sex_e1**2+self.sex_e2**2)
        self.sex_FWHM = e[4]

        self.mine_X2=0
        self.mine_Y2=0
        self.mine_XY=0
        self.mine_e1 =0
        self.mine_e2 =0
        self.mine_ellipticity= 0


        self.moffat_A = 0
        self.moffat_B = 0
        self.moffat_C = 0
        self.moffat_max = 0
        self.moffat_background = 0
        self.moffat_amplitude = 0
        self.moffat_CenterX = 0
        self.moffat_CenterY = 0
        self.moffat_smallest_fwhm =0
        self.moffat_largest_fwhm=0
        self.moffat_phi=0
        self.moffat_ellipticity =0
        self.moffat_red_chi_square = 0
        self.moffat_beta = 0

    def star_info(self):
        info = [self.ID, self.RA_grid, self.DEC_grid,
                self.sex_X2, self.sex_Y2, self.sex_XY,
                self.sex_e1, self.sex_e2,self.sex_ellipticity,self.sex_FWHM,
                self.sex_correct_e1, self.sex_correct_e2,
                self.moffat_CenterX, self.moffat_CenterY, self.moffat_smallest_fwhm, self.moffat_largest_fwhm, self.moffat_phi, self.moffat_ellipticity,
                ]
        info_str =  [str(i) for i in info]
        return info_str

    def update_mine_ellpticity(self):
        sumAll = np.sum(self.data)
        for i in range(self.data.shape[0]):
            for j in range(self.data.shape[1]):
                self.mine_X2 += self.data[i][j]*i*i
                self.mine_Y2 += self.data[i][j]*j*j
                self.mine_XY += self.data[i][j]*i*j
        self.mine_X2 =self.mine_X2/sumAll
        self.mine_Y2 =self.mine_Y2/sumAll
        self.mine_XY =self.mine_XY/sumAll
        self.mine_e1 = (self.mine_X2-self.mine_Y2)/(self.mine_X2+self.mine_Y2)
        self.mine_e2 = (2*self.mine_XY)/(self.mine_X2+self.mine_Y2)
        self.mine_ellipticity = np.sqrt(self.mine_e1**2 + self.mine_e2**2)

    def update_wcs(self):
        self.moffat_RA, self.moffat_DEC = self.wcs.all_pix2world(self.moffat_global_CenterX,self.moffat_global_CenterY,0)

    def update_moffat(self):
        moffat_results = astro.fit_moffat_elliptical([0,0],self.data)

        self.moffat_max = moffat_results[0]
        self.moffat_background = moffat_results[1]
        self.moffat_amplitude = moffat_results[2]
        self.moffat_CenterX = moffat_results[3]
        self.moffat_CenterY = moffat_results[4]
        self.moffat_smallest_fwhm = moffat_results[5]
        self.moffat_largest_fwhm  = moffat_results[6]
        self.moffat_phi  = moffat_results[7]
        self.moffat_beta = moffat_results[8]

        phi = moffat_results[7]
        beta = moffat_results[8]

        denom = 2*np.sqrt(np.power(2,1/self.moffat_beta)-1)
        alpha1=self.moffat_smallest_fwhm/denom
        alpha2=self.moffat_largest_fwhm/denom

        self.moffat_A = (np.cos(phi)/alpha1)**2 + (np.sin(phi)/alpha2)**2
        self.moffat_B = (np.sin(phi)/alpha1)**2 + (np.cos(phi)/alpha2)**2
        self.moffat_C = 2*np.sin(phi)*np.cos(phi)*(1/alpha1**2-1/alpha2**2)
        self.moffat_global_CenterX = self.origin[1] + self.moffat_CenterY
        self.moffat_global_CenterY = self.origin[0] + self.moffat_CenterX


    def update(self, Moffat, Ellipticity):
        #self.update_mine_ellpticity()
        if Ellipticity==True:
            self.update_mine_ellpticity()
        if Moffat==True:
            self.update_moffat()

    def plotFitting(self):

        data = self.data
        xmargin = np.sum(data,axis=1)
        ymargin = np.sum(data,axis=0)

        xx = np.arange(0,len(xmargin))-self.moffat_CenterX-0.5
        yy = np.arange(0,len(ymargin))-self.moffat_CenterY-0.5

        f, (ax3,ax4)=plt.subplots(1,2,sharex='col', sharey='row')
        xx3=np.arange(-1*self.winX, self.winX,.1)
        yy3=np.arange(-1*self.winY, self.winY,.1)
        X, Y = np.meshgrid(xx3, yy3)
        moffat_Z = self.moffat_background + self.moffat_amplitude/np.power(1+self.moffat_A*X**2
                                                                    +self.moffat_B*Y**2
                                                                    +self.moffat_C*X*Y,
                                                                    self.moffat_beta)
        moffat_xmargin = np.sum(moffat_Z, axis=1)/10
        moffat_ymargin = np.sum(moffat_Z, axis=0)/10

        ax3.set_xlabel("X")
        ax4.set_xlabel("Y")
        ax3.bar(xx,xmargin)
        ax3.plot(xx3,moffat_xmargin,'-r', linewidth=3.0)
        ax3.set_xlim([-10,10])
        ax4.bar(xx,ymargin)
        ax4.plot(xx3,moffat_ymargin,'-r', linewidth=3.0)
        ax4.set_xlim([-10,10])
        imgName = "fitting_chip_"+"_ID_"+str(self.ID) +".eps"
        plt.show()
        #plt.savefig(imgName, format='eps', dpi=1000)
        #plt.close()



def generateStarList(fitsFileName, catFileName, (winX, winY), (limX, limY), Moffat, Ellipticity):
    starList=[]
    catFile=open(catFileName,'r')
    imgFile = fits.open(fitsFileName)
    wcs = WCS(fitsFileName)
    image_map = imgFile[0].data


    for line in catFile.readlines():
        if line[0]!='#':
            temp=line.split()
            ID = float(temp[0])
            sex_CenterX = float(temp[1])
            sex_CenterY = float(temp[2])
            sex_RA  = float(temp[3])
            sex_DEC = float(temp[4])
            sex_X2 = float(temp[5])
            sex_Y2 = float(temp[6])
            sex_XY = float(temp[7])
            sex_elipticity = float(temp[8])
            sex_FWHM = float(temp[10])
            e = (sex_X2, sex_Y2, sex_XY, sex_elipticity, sex_FWHM)
            if sex_CenterY>limY[0] and sex_CenterY< limY[1] and sex_CenterX>limX[0] and sex_CenterX < limX[1]:

                data = image_map[np.floor(sex_CenterY-winY):np.floor(sex_CenterY+winY)+1,np.floor(sex_CenterX-winX):np.floor(sex_CenterX+winX)+1]
                origin = (np.floor(sex_CenterY)-winY, np.floor(sex_CenterX)-winX)


                star = Star(ID, sex_CenterX, sex_CenterY, sex_RA, sex_DEC, data, wcs, origin, (winX, winY), e)
                star.update(Moffat, Ellipticity)
                starList.append(star)

    catFile.close()
    return starList




def generateStarListDictionary(fitsNameList, (winX, winY), (xlim, ylim), sex, Moffat, Ellipticity, dir=None ):
    keys = fitsNameList
    starDict = dict.fromkeys(keys)
    threadList =[]
    for i in range(len(fitsNameList)):
        fitsName = fitsNameList[i]
        if dir!=None:
            fullFitsName = str(dir)+fitsNameList[i]
            catName = str(dir)+fitsName[:-5] +"_sex.cat"
        else:
            fullFitsName = fitsNameList[i]
            catName = fitsName[:-5] +"_sex.cat"
        if sex ==True:
            subprocess.call("sex " + fullFitsName + " -CATALOG_NAME " + catName, shell=True)
        starList= generateStarList(fullFitsName, catName, (winX, winY), (xlim, ylim), Moffat, Ellipticity)
        starDict[fitsName] = starList


    return starDict






