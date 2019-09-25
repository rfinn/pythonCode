#!/usr/bin/env python
'''
writing to help run sextractor and make some common plots

working on reduction of mosaic Halpha imaging.
Becky has done basic reductions.  We still need to
convolve images to a common seeing and stack.
Then do continuum subtraction.

GOAL:
- measure FWHM of an image using sextractor


PROCEDURE:

USAGE:

REQUIRED MODULES:


written by Rose A. Finn
01-March-2014

'''

import os
import pyfits
import atpy
from pylab import *


# flux and class_star cuts to pick out non-saturated stars
minflux=1.e5
maxflux=1.e6
minclass_star=0.97
maxclass_star=1.0

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

class image:
    def __init__(self,image,config_file='default.sex.mosaic'):
        self.image=image
        self.prefix,t=self.image.split('.fits')
        self.config_file=config_file
        try:
            self.readcat_atpy()
        except IOError:
            self.runsextractor()
            self.readcat_atpy()
        self.get_fwhm()
        print self.prefix,': FWHM = ',self.fwhm,'+/-',self.fwhm_std
        self.plot_classstar_fluxauto()
        self.plot_fwhm()
    def runsextractor(self):
        s='sex '+self.image+' -c '+self.config_file+' -CATALOG_NAME '+self.prefix+'_test_cat.fits -CATALOG_TYPE FITS_LDAC'
        os.system(s)

    def readcat(self):
        #self.cat=atpy.Table(self.prefix+'_test_cat.fits',hdu=2)
        print self.prefix+'_test_cat.fits'
        f=pyfits.open(self.prefix+'_test_cat.fits')
        self.cat=f[2].data
        f.close()
        self.starflag=(self.cat.FLUX_AUTO > minflux) & (self.cat.FLUX_AUTO < maxflux) & (self.cat.CLASS_STAR > minclass_star)
        
    def readcat_atpy(self):
        self.cat=atpy.Table(self.prefix+'_test_cat.fits',hdu=2)
        #print self.prefix+'_test_cat.fits'
        #f=pyfits.open(self.prefix+'_test_cat.fits')
        #self.cat=f[2].data
        #f.close()
        self.starflag=(self.cat.FLUX_AUTO > minflux) & (self.cat.FLUX_AUTO < maxflux) & (self.cat.CLASS_STAR > minclass_star)
    def plot_classstar_fluxauto(self):
        figure()
        plot(self.cat.FLUX_AUTO,self.cat.CLASS_STAR,'k.')
        plot(self.cat.FLUX_AUTO[self.starflag],self.cat.CLASS_STAR[self.starflag],'b.')
        gca().set_xscale('log')
        xlabel('FLUX_AUTO')
        ylabel('CLASS_STAR')
        title(self.prefix)
        boxparams=array([(minflux+maxflux)/2.,(minclass_star+maxclass_star)/2.,(maxflux-minflux),(maxclass_star-minclass_star),0],'f')
        drawbox(boxparams,'r-')
        savefig(self.prefix+'_classstar_fluxauto.png')
    def get_fwhm(self):
        m=mean(self.cat.FWHM_IMAGE[self.starflag])
        s=std(self.cat.FWHM_IMAGE[self.starflag])
        keepflag=(abs(m-self.cat.FWHM_IMAGE) < s) & self.starflag #& (self.cat.FLAGS < 1)
        m=self.cat.FWHM_IMAGE[keepflag]
        self.fwhm = mean(self.cat.FWHM_IMAGE[keepflag])
        self.fwhm_std = std(self.cat.FWHM_IMAGE[keepflag])
        #self.fwhm =(m)
        #self.fwhm_std = s
        
    def plot_classstar_magauto(self):
        figure()
        plot(self.cat.MAG_AUTO,self.cat.CLASS_STAR,'k.')
        #gca().set_xscale('log')
        xlabel('MAG_AUTO')
        ylabel('CLASS_STAR')
        title(self.prefix)
        savefig(self.prefix+'_classstar_magauto.png')
    def plot_fwhm(self):
        figure()
        hbins=arange(self.fwhm-2,self.fwhm+2.,.05)
        hbins=arange(0,30,.1)
        hist(self.cat.FWHM_IMAGE,bins=hbins)
        hist(self.cat.FWHM_IMAGE[self.starflag],bins=hbins)
        axvline(self.fwhm,color='r')
        axvline(self.fwhm-self.fwhm_std,color='r',ls='--')
        axvline(self.fwhm+self.fwhm_std,color='r',ls='--')
        xlabel('FWHM_IMAGE')
        s='fwhm = %5.2f +/- %5.2f'%(self.fwhm,self.fwhm_std)
        ax=gca()
        text(.1,.9,s,transform=ax.transAxes,horizontalalignment='left',fontsize=12)
        title(self.prefix)
        axis([2,10,0,200])
        savefig(self.prefix+'_fwhm_hist.png')
