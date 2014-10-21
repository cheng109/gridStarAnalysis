"""
  @package 
  @file moffatFitting.py
  @brief compute the fwhm from the 'segmentation.fits' and 'image.fits'
 
  @brief Created by:
  @author Jun Cheng (Purdue)
  
  @warning This code is not fully validated
  and not ready for full release.  Please
  treat results with caution.

Usage: python moffatFitting.py <image.fits> 
Notes: 1.

"""
from aspylib import astro
import sys,os,commands
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
#import mayavi.mlab as mlab

class Star:
    def __init__(self, ID, globalPos, centerX, centerY):
        self.ID = ID
        self.globalPos = globalPos
        self.empty = 1
        self.boundary = 0
        self.global_R = []

        self.fwhm = 0
        self.ellipticity = 0
        self.centerX = centerX
        self.centerY = centerY
        self.sigma = 0

        self.sigmaX=0
        self.sigmaY=0
       
        self.gauss_fwhmX=0
        self.gauss_fwhmY=0
        self.moffat_fwhmX=0
        self.moffat_fwhmY=0

        self.winSizeX=100
        self.winSizeY=100

        self.gauss_A = 0
        self.gauss_B = 0
        self.gauss_C = 0
        self.gauss_max = 0
        self.gauss_background = 0
        self.gauss_amplitude = 0
        self.gauss_centerX = 0
        self.gauss_centerY = 0
        self.gauss_smallest_fwhm =0
        self.gauss_largest_fwhm=0
        self.gauss_phi=0
        self.gauss_ellipticity =0
        self.gauss_red_chi_square = 0

        self.moffat_A = 0
        self.moffat_B = 0
        self.moffat_C = 0
        self.moffat_max = 0
        self.moffat_background = 0
        self.moffat_amplitude = 0
        self.moffat_centerX = 0
        self.moffat_centerY = 0
        self.moffat_smallest_fwhm =0
        self.moffat_largest_fwhm=0
        self.moffat_phi=0
        self.moffat_ellipticity =0
        self.moffat_red_chi_square = 0
        self.moffat_beta = 0

    def update(self, image):
        
        self.global_R = np.sqrt((self.centerX-20000)**2+(self.centerY-600)**2)
        data = image[int(self.centerY-self.winSizeY):int(self.centerY+self.winSizeY),int(self.centerX-self.winSizeX):int(self.centerX+self.winSizeX)]

        gauss_results=astro.fit_gauss_elliptical([0,0],data)

        self.gauss_max = gauss_results[0]
        self.gauss_back = gauss_results[1]
        self.gauss_amplitude = gauss_results[2]
        self.gauss_centerX = gauss_results[3]
        self.gauss_centerY = gauss_results[4]
        self.gauss_smallest_fwhm =fwhm1= gauss_results[5]
        self.gauss_largest_fwhm =fwhm2 = gauss_results[6]
        self.gauss_phi = phi =gauss_results[7]
       
       
        self.gauss_fwhmX=np.sqrt((fwhm2*np.cos(phi))**2+(fwhm1*np.sin(phi))**2)      
        self.gauss_fwhmY=np.sqrt((fwhm2*np.sin(phi))**2+(fwhm1*np.cos(phi))**2)
        self.gauss_ellipticity = 1-fwhm1/fwhm2

        # gauss_chi_square
        ssum =0
        S0 = gauss_results[1]
        S1 = gauss_results[2]
        x0 = gauss_results[3]
        y0 = gauss_results[4]
        smallest_fwhm = gauss_results[5]
        largest_fwhm = gauss_results[6]
        phi = gauss_results[7]
        
        sigma1=self.gauss_smallest_fwhm/2.3548
        sigma2=self.gauss_largest_fwhm/2.3548
        self.gauss_A = (np.cos(phi)/sigma1)**2 + (np.sin(phi)/sigma2)**2
        self.gauss_B = (np.sin(phi)/sigma1)**2 + (np.cos(phi)/sigma2)**2
        self.gauss_C = 2*np.sin(phi)*np.cos(phi)*(1/sigma1**2-1/sigma2**2)
        self.n=0
        for x in range(data.shape[0]):
            for y in range(data.shape[1]):
                if data[x][y]>0:
                    self.n+=1
                    delta = S0+S1*np.exp(-0.5*(self.gauss_A*(x-x0)**2+self.gauss_B*(y-y0)**2+self.gauss_C*(x-x0)*(y-y0)))-data[x][y]
                    ssum = ssum + float(delta**2)/(data[x][y])
        self.gauss_red_chi_square = ssum/(self.n-1-7)


        moffat_results=astro.fit_moffat_elliptical([0,0],data)

        self.moffat_max = moffat_results[0]
        self.moffat_background = moffat_results[1]
        self.moffat_amplitude = moffat_results[2]
        self.moffat_centerX = moffat_results[3]
        self.moffat_centerY = moffat_results[4]
        self.moffat_smallest_fwhm =fwhm1= moffat_results[5]
        self.moffat_largest_fwhm =fwhm2 = moffat_results[6]
        self.moffat_phi  = moffat_results[7]
        self.moffat_beta = moffat_results[8]
        
        self.moffat_fwhmX=np.sqrt((fwhm2*np.cos(phi))**2+(fwhm1*np.sin(phi))**2)
        self.moffat_fwhmY=np.sqrt((fwhm2*np.sin(phi))**2+(fwhm1*np.cos(phi))**2)
        self.moffat_ellipticity = 1-fwhm1/fwhm2
        ssum =0
        S0 = moffat_results[1]
        S1 = moffat_results[2]
        x0 = moffat_results[3]
        y0 = moffat_results[4]
        smallest_fwhm = moffat_results[5]
        largest_fwhm = moffat_results[6]
        phi = moffat_results[7]
        beta = moffat_results[8]

        denom = 2*np.sqrt(np.power(2,1/self.moffat_beta)-1)
        alpha1=self.moffat_smallest_fwhm/denom
        alpha2=self.moffat_largest_fwhm/denom
        self.moffat_A = (np.cos(phi)/alpha1)**2 + (np.sin(phi)/alpha2)**2
        self.moffat_B = (np.sin(phi)/alpha1)**2 + (np.cos(phi)/alpha2)**2
        self.moffat_C = 2*np.sin(phi)*np.cos(phi)*(1/alpha1**2-1/alpha2**2)
        self.n=0
        for x in range(data.shape[0]):
            for y in range(data.shape[1]):
                if data[x][y]>0:
                    self.n+=1
                    delta = S0+S1/np.power(1+self.moffat_A*(x-x0)**2+self.moffat_B*(y-y0)**2+self.moffat_C*(x-x0)*(y-y0),beta)-data[x][y]
                    ssum = ssum + float(delta**2)/(data[x][y])
        self.moffat_red_chi_square = ssum/(self.n-1-8)


    def fitPlot(self, image):
         
        data = image[int(self.centerY-self.winSizeY):int(self.centerY+self.winSizeY),int(self.centerX-self.winSizeX):int(self.centerX+self.winSizeX)]
        self.sigmaX = self.gauss_fwhmX/2.3548
        self.sigmaY = self.gauss_fwhmY/2.3548
        
        xmargin = np.sum(data,axis=1)
        ymargin = np.sum(data,axis=0)
        
        xx = np.arange(0,len(xmargin))-self.gauss_centerX
        yy = np.arange(0,len(ymargin))-self.gauss_centerY

        f, ((ax1,ax2),(ax3,ax4))=plt.subplots(2,2,sharex='col', sharey='row')
        
        ax1.plot(xx,xmargin,'*')#,label="FWHM_X=%.3f"%self.gauss_fwhmX)
        ax2.plot(yy,ymargin,'*')#,label="FWHM_Y=%.3f"%self.gauss_fwhmY)
        xx2=np.arange(-1*self.winSizeX, self.winSizeX,0.1) #-self.gauss_centerX
        yy2=np.arange(-1*self.winSizeY, self.winSizeY,0.1) #-self.gauss_centerY
        alongX = np.sum(data)*np.sqrt(1/(2*np.pi*self.sigmaX**2))*np.exp(-xx2**2/(2*self.sigmaX**2)) 
        alongY = np.sum(data)*np.sqrt(1/(2*np.pi*self.sigmaY**2))*np.exp(-yy2**2/(2*self.sigmaY**2)) 
        
        xx3=np.arange(-1*self.winSizeX, self.winSizeX,1) #-self.gauss_centerX
        yy3=np.arange(-1*self.winSizeY, self.winSizeY,1) #-self.gauss_centerY
        
        X, Y = np.meshgrid(xx3, yy3)
        gauss_Z = self.gauss_background + self.gauss_amplitude*np.exp(-0.5*(self.gauss_A*X**2
                                                                            +self.gauss_B*Y**2
                                                                            +self.gauss_C*X*Y))
        moffat_Z = self.moffat_background + self.moffat_amplitude/np.power(1+self.moffat_A*X**2
                                                                    +self.moffat_B*Y**2
                                                                    +self.moffat_C*X*Y,
                                                                    self.moffat_beta)
        gauss_xmargin = np.sum(gauss_Z, axis=1)
        gauss_ymargin = np.sum(gauss_Z, axis=0)
        moffat_xmargin = np.sum(moffat_Z, axis=1)
        moffat_ymargin = np.sum(moffat_Z, axis=0)
        ax1.plot(xx3,gauss_xmargin,'-')
        ax2.plot(xx3,gauss_xmargin,'-')
        
        ax1.set_ylabel("counts")
        ax1.set_xlabel("X")
        ax2.set_xlabel("Y")
        ax3.set_xlabel("X")
        ax4.set_xlabel("Y")

        ax1.set_title("Gassian: A=%.3f"%self.gauss_largest_fwhm + ", B=%.3f"%self.gauss_smallest_fwhm)
        ax3.set_title("Moffat: A=%.3f"%self.moffat_largest_fwhm + ", B=%.3f"%self.moffat_smallest_fwhm)        

        ax3.plot(xx,xmargin,'*')#,label="FWHM_X=%.3f"%self.moffat_fwhmX)
        ax3.plot(xx3,moffat_xmargin,'-')
        ax4.plot(yy,ymargin,'*')#, label="FWHM_Y=%.3f"%self.moffat_fwhmY)
        ax4.plot(yy3,moffat_ymargin,'-')
        
        imgName = "fitting_chip_"+str(self.globalPos)+"_ID_"+str(self.ID) +".eps"
        plt.savefig(imgName, format='eps', dpi=1000)
        plt.close()

    def surfPlot(self,image):
        data = image[int(self.centerY-self.winSizeY):int(self.centerY+self.winSizeY),int(self.centerX-self.winSizeX):int(self.centerX+self.winSizeX)]
        mlab.clf()
        x,y=np.mgrid[0:data.shape[0],0:data.shape[1]]
        print data.shape
        print x.shape,y.shape
        print x,y,data
        mlab.surf(data*0.1,warp_scale='auto')
        mlab.surf(data*0.12,warp_scale='auto')
        mlab.show()
        mlab.close()

def writeTofile(starList, fileName):
    file = open(fileName,'w')
    for i in range(len(starList)):
        #if starList[i].empty==0 and starList[i].boundary==0:
        file.write(str(starList[i].ID)+"\t"
                   +str(starList[i].globalPos)+"\t"
                   # Gaussian results
                   +str(starList[i].centerX-20000+40000*starList[i].globalPos)+"\t"
                   +str(starList[i].gauss_smallest_fwhm)+"\t"
                   +str(starList[i].gauss_largest_fwhm)+"\t"
                   +str(starList[i].gauss_smallest_fwhm*starList[i].gauss_largest_fwhm*np.pi/4)+"\t"
                   +str(starList[i].gauss_ellipticity)+"\t"
                   +str(starList[i].gauss_phi)+"\t"
                   +str(starList[i].gauss_red_chi_square)+"\t"
                   # Moffat results
                   +str(starList[i].moffat_smallest_fwhm)+"\t"
                   +str(starList[i].moffat_largest_fwhm)+"\t"
                   +str(starList[i].moffat_smallest_fwhm*starList[i].moffat_largest_fwhm*np.pi/4)+"\t"
                   +str(starList[i].moffat_ellipticity)+"\t"
                   +str(starList[i].moffat_phi)+"\t"
                   +str(starList[i].moffat_red_chi_square)+"\t"
                   +"\n")
    file.close()
    
def readCat(filename):
    f=open(filename,'r')
    centerX=[]
    centerY=[]
    ID =[]
    lowerLim = 200#2030
    upperLim = 900#2080

    leftLim = 200
    rightLim = 40000-200
    for line in f.readlines():
        if line[0]!='#':
            temp=line.split()
            if float(temp[5])>lowerLim and float(temp[5])<upperLim and float(temp[4])>leftLim and float(temp[4])<rightLim: 
                centerX.append(temp[4])
                centerY.append(temp[5])
                ID.append(int(temp[0]))
    f.close()
    return ID, centerX, centerY


def generateStarList(imageFile, catFile, globalPos):
    f1 =fits.open(imageFile)
    image =f1[0].data    
    ID, centerX, centerY = readCat(catFile)
    numStar=len(centerX)
    starList=[]
    for i in range(numStar):
        starList.append(Star(ID[i], globalPos, float(centerX[i]),float(centerY[i])))
        starList[i].update(image)
        starList[i].fitPlot(image)
    return starList, image
      
starList, image = generateStarList(sys.argv[1],sys.argv[2], int(sys.argv[3]))
writeTofile(starList,"x_y_elipticity_fwhm.txt")
#starList[int(sys.argv[4])-1].surfPlot(image)

