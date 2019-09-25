#!/usr/bin/env python
'''
GOAL:
  - first program to analyze WISEsize data
    
  
'''
import pylab as pl
import numpy as np
import os
import mystuff as my
from astropy.io import fits
from astropy.table import Table
from astropy.table import Column
#import astropy.cosmology.funcs.angular_diameter_distance as DA
import astrofuncs
from LCScommon import *
minsize_kpc = 6.


def drawbox(data,style):#feed in center x,y,dx,dy,rotation E of N
    #xcoords of unrotated box, going around CCW
    xl=array([data[0]-0.5*data[2],data[0]+0.5*data[2],data[0]+0.5*data[2],data[0]-0.5*data[2],data[0]-0.5*data[2]],'d')
    yl=array([data[1]-0.5*data[3],data[1]-0.5*data[3],data[1]+0.5*data[3],data[1]+0.5*data[3],data[1]-0.5*data[3] ],'d')

    xl=array([-0.5*data[2],+0.5*data[2],+0.5*data[2],-0.5*data[2],-0.5*data[2]],'d')
    yl=array([-0.5*data[3],-0.5*data[3],+0.5*data[3],+0.5*data[3],-0.5*data[3] ],'d')

    ang=data[4]*pi/180.*-1.#convert rotation to radians
    #rotate coordinates
    xp=cos(ang)*xl-sin(ang)*yl
    yp=sin(ang)*xl+cos(ang)*yl

    #put back on absolute scale
    xp=data[0]+xp
    yp=data[1]+yp
    #draw rotated box
    plot(xp,yp,style)


class wisegal():
    def __init__(self):
        hdulist=fits.open(homedir+'research/WISEsize/LCS_WISE_Zoo_spirals.fits')
        self.s=hdulist[1].data
        hdulist.close()
        #self.DA = self.s.SERSIC_TH50*astrofuncs.DA(self.s.ZDIST,H0/100.)
        self.cullsample()
        self.sdssspecflag=(self.s.ISDSS > -1)
        self.calcagn()
        # print fraction of AGN
        print 'AGN  = ',sum(self.agnflag & self.sampleflag & self.sdssspecflag),' out of ',sum(self.sampleflag & self.sdssspecflag)
        self.agnfrac=1.*sum(self.agnflag & self.sampleflag & self.sdssspecflag)/sum(self.sampleflag & self.sdssspecflag)
        print 'AGN fraction = %5.2f'%(self.agnfrac)
        print 'Remaining sample = ',sum(self.sampleflag & self.sdssspecflag)*(1.-self.agnfrac)
    def cullsample(self):
        self.sampleflag =   (self.s.W3MAG < 11.) & (self.s.SERSIC_TH50 > 5.5) & (self.s.ZDIST > 0.015) & (self.s.ZDIST < 0.043)
    def calcstellarmass(self):
        self.logstellarmassTaylor=1.15+0.70*(self.s.ABSMAG[:,3]-self.s.ABSMAG[:,5]) -0.4*(self.s.ABSMAG[:,5]+ 5.*log10(h))

    def calcagn(self):
        self.AGNKAUFF= ((log10(self.s.O3FLUX/self.s.HBFLUX) > (.61/(log10(self.s.N2FLUX/self.s.HAFLUX)-.05)+1.3)) | (log10(self.s.N2FLUX/self.s.HAFLUX) > 0.))
#y=(.61/(x-.47)+1.19)
        self.AGNKEWLEY= ((log10(self.s.O3FLUX/self.s.HBFLUX) > (.61/(log10(self.s.N2FLUX/self.s.HAFLUX)-.47)+1.19)) | (log10(self.s.N2FLUX/self.s.HAFLUX) > 0.3))
        self.agnflag = self.AGNKAUFF
        self.agnflag = self.AGNKEWLEY
    def plotagn(self):
         figure()
         clf()
         keepflag=self.sampleflag & self.sdssspecflag
         #keepflag=self.sdssspecflag
         x=np.log10(self.s.N2FLUX/self.s.HAFLUX)
         y=np.log10(self.s.O3FLUX/self.s.HBFLUX)

         plot(x[keepflag],y[keepflag],'k.')
         pl.plot(x[self.AGNKAUFF & keepflag],y[self.AGNKAUFF & keepflag],'k*',mec='None',markersize=4, label='_nolabel_')
         pl.plot(x[self.agnflag & keepflag],y[self.agnflag & keepflag],'k*',mec='None',markersize=14, label='_nolabel_')
         #pl.plot(x[self.agnflag],y[self.agnflag],'ro',markersize=4, label='_nolabel_')
         #draw AGN diagnostic lines
         x=arange(-3,.4,.01)
         y=(.61/(x-.47)+1.19)
         #Kewley
         plot(x,y,'c',label='Kewley & Dopita 2002')
         x=arange(-3,.0,.01)
         y =(.61/(x-.05)+1.3)#Kauffman 2003?
         plot(x,y,'g',label='Kauffmann et al. 2003')
         y = ((-30.787+(1.1358*x)+((.27297)*(x)**2))*tanh(5.7409*x))-31.093 #Stasinska 2006	    
         plot(x,y,'r',label='Stasinska et al. 2006')

         pl.axis([-1.5,.49,-1.,1.5])
         pl.xlabel(r'$\log_{10}(NII/H\alpha)$',fontsize=20)
         pl.ylabel(r'$\log_{10}(OIII/H\beta)$',fontsize=20)
         pl.legend(loc='upper left',prop={'size':12})
         pl.savefig(homedir+'research/LocalClusters/SamplePlots/AGNclassification.png')
    def plotredshifthist(self):
        figure()
        hist(self.s.ZDIST[self.sampleflag],bins=30)
        xlabel('$Redshift$',fontsize=20)
        ylabel('$Number \ of \ Galaxies$',fontsize=20)
    def plotcolormag(self):
        figure()
        hexbin(self.s.ABSMAG[:,4][self.sampleflag],self.s.ABSMAG[:,1][self.sampleflag]-self.s.ABSMAG[:,4][self.sampleflag],cmap='gray_r')#,'k.')#,color='0.5',alpha='0.5')
        axis([-23,-17,0,7])
    def plotpositions(self):
        figure(figsize=(10,4))
        subplots_adjust(bottom=.15,top=.95,left=.1,right=.95)
        scatter(self.s.RA[self.sampleflag],self.s.DEC[self.sampleflag],s=10,alpha=.25,marker='.')#,label='This proposal')#,color='0.5',alpha='0.5')
        legend(['This proposal'],loc='upper left')#,'LCS'])
        for i in range(len(clusternames)):
            s=' '+clusternames[i]
            text(clusterRA[clusternames[i]],clusterDec[clusternames[i]],s,fontsize=12,color='r')
            #cra.append(clusterRA[clusternames[i]])
            #cdec.append(clusterDec[clusternames[i]])
            drawbox(cluster24Box[clusternames[i]],'r-')
        xlim(100,275)
        ylim(-5,70)
        xlabel('$RA \ (deg) $',fontsize=14)
        ylabel('$DEC \ (deg) $',fontsize=14)
        savefig('/Users/rfinn/proposals/NASA2017/WISEsize-LCS.pdf')

        

w=wisegal()
