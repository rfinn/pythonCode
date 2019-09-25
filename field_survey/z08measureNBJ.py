#!/usr/bin/env python
'''
writing code to run sextractor on aligned NB and J images

'''
from z08common import *
from pylab import *
import os
#import asciidata
import atpy
from pyraf import iraf


RAlist=['13:17','17:16','13:50']

NB='rcl1317-F1187N-image-2013-04-01T05h.fits'
J='rcl1317-JX-image-2013-02-02T10h.fits'
NBJratio=0.073
newfirmpixelscale=0.4
def getheader(im,keyword):
    iraf.imgets(image=im,param=keyword)
    t=float(iraf.imgets.value)
    return t
def getgain(im):
    iraf.imgets(image=im,param='GAIN')
    t=float(iraf.imgets.value)
    return t

def getgain(im):
    iraf.imgets(image=im,param='GAIN')
    t=float(iraf.imgets.value)
    return t


def runsextractor():
    zp=getheader(NB,'MAGZERO')
    #gain=getheader(NB,'GAIN')
    #ncombine=getheader(NB,'NCOMBINE')
    #effective_gain=gain*ncombine
#s='sex '+NB+' -c /Users/rfinn/research/LocalClusters/sextractor/default.sex.newfirm -FILTER_NAME /Users/rfinn/research/LocalClusters/sextractor/default.conv -GAIN '+str(effective_gain)+' -MAG_ZEROPOINT '+str(zp)
    s='sex '+NB+' -c /Users/rfinn/research/LocalClusters/sextractor/default.sex.newfirm -FILTER_NAME /Users/rfinn/research/LocalClusters/sextractor/default.conv -MAG_ZEROPOINT '+str(zp)+' -CATALOG_NAME testNB.fits -CATALOG_TYPE FITS_1.0 '
    print s
    os.system(s)
    zp=getheader(J,'MAGZERO')
    s='sex '+NB+','+J+' -c /Users/rfinn/research/LocalClusters/sextractor/default.sex.newfirm -FILTER_NAME /Users/rfinn/research/LocalClusters/sextractor/default.conv -MAG_ZEROPOINT '+str(zp)+' -CATALOG_NAME testJ.fits -CATALOG_TYPE FITS_1.0 '
    print s
    os.system(s)

#runsextractor()

class cluster:
    def __init__(self,clustername):
        self.radeg=clusterRAdeg[clustername]
        self.decdeg=clusterDecdeg[clustername]
        self.prefix=clustername

        self.nbimage='rcl1317-F1187N-image-2013-04-01T05h.fits'
        self.jimage='rcl1317-JX-image-2013-02-02T10h.fits'

        self.nbzp=getheader(self.nbimage,'MAGZERO')
        self.jzp=getheader(self.jimage,'MAGZERO')
        #self.readSextractor('testJ.cat')
        #self.readSextractorNB('testNB.cat')
        infile='testNB.fits'
        self.nb=atpy.Table(infile)
        infilej='testJ.fits'
        self.j=atpy.Table(infilej)
        self.ximage=self.nb.X_IMAGE
        self.yimage=self.nb.Y_IMAGE
        self.ra=self.nb.ALPHA_J2000
        self.dec=self.nb.DELTA_J2000
        self.nbfwhm=self.nb.FWHM_IMAGE
        self.jfwhm=self.j.FWHM_IMAGE
        self.nbmag=self.nb.MAG_APER[:,0]
        self.jmag=self.j.MAG_APER[:,0]
        self.nbmag1=self.nb.MAG_APER[:,1]
        self.jmag1=self.j.MAG_APER[:,1]
        self.nbflux1=10.**((self.nbmag1-self.nbzp)/-2.5)
        self.jflux1=10.**((self.jmag1-self.jzp)/-2.5)
        #self.nbmag2=self.nb.MAG_APER2'],'f')
        #self.jmag2=array(self.j['MAG_APER2'],'f')

        self.goodflag=(self.ximage >120) & (self.yimage < 3750) & (self.yimage > 70) & ( ((self.ximage-295)**2 + (self.yimage-56.)**2) > 438.**2) 

        self.supersource= (self.nbmag < 18) & ((self.jmag1-self.nbmag1) > 2)
        self.nbsource=self.goodflag & (self.nbmag1 < 22) & ((self.jmag1-self.nbmag1) > .9) & ~self.supersource
        self.contsubflux=self.nbflux1-NBJratio*self.jflux1
        self.contsubfluxerr=sqrt(self.nb.MAGERR_APER[:,1]**2-(NBJratio)**2*(self.j.MAGERR_APER[:,1])**2)
        self.contsubsnr=abs(self.contsubflux/self.contsubfluxerr)
    def writenbregfile(self):
        f=homedir+'research/z08clusters/RegionsFiles/'+self.prefix+'_NBexcess.reg'
        outfile=open(f,'w')
        outfile.write('global color=green width=2 \n')
        outfile.write('fk5\n')

        ra=self.ra[self.nbsource]
        dec=self.dec[self.nbsource]
        snr=self.contsubsnr[self.nbsource]
        for i in range(len(ra)):
            outfile.write('circle(%12.8f, %12.8f, 3\") # color=red \n'%(ra[i],dec[i]))
            if snr[i] > 2.5:
                outfile.write('circle(%12.8f, %12.8f, 15\") # color=blue \n'%(ra[i],dec[i]))
        ra=self.ra[self.supersource]
        dec=self.dec[self.supersource]
        snr=self.contsubsnr[self.supersource]
        for i in range(len(ra)):
            outfile.write('circle(%12.8f, %12.8f, 10\") # color=green \n'%(ra[i],dec[i]))
            if snr[i] > 2.5:
                outfile.write('circle(%12.8f, %12.8f, 15\") # color=blue \n'%(ra[i],dec[i]))

        outfile.close()
    def comparefwhm(self):
        mybins=arange(1,5,.05)
        figure()
        hist(self.nb.FWHM_IMAGE,bins=mybins,histtype='step',color='b',label='NB')
        hist(self.j.FWHM_IMAGE,bins=mybins,histtype='step',color='r',label='J')
                         
    def plotNBJflux(self):
        figure(figsize=(8,6))
        m1=self.nbflux1
        m2=self.jflux1
        subplots_adjust(left=.1,right=.95)
        y=self.nbflux1-NBJratio*self.jflux1
        ery=sqrt(self.nb.MAGERR_APER[:,1]**2-(NBJratio)**2*(self.j.MAGERR_APER[:,1])**2)
        flag=self.goodflag & (self.j.MAG_APER[:,1] > 17.)
        for i in range(2):
            subplot(2,1,i+1)
            if i == 0:
                x=m1
                xl='NB Flux'
            elif i == 1:
                x=m2
                xl='J Flux'
            sp=scatter(x[flag],y[flag],c=self.contsubsnr[flag],vmin=1,vmax=15)
            colorbar(sp)
            #errorbar(x[self.goodflag],y[self.goodflag],yerr=ery[self.goodflag],fmt=None)
            xlabel(xl)
            ylabel('J-NB Flux')
            axhline(y=0,color='k',ls='--')
            gca().set_xscale('log')
            axis([.3,3000,-1,4])

    def plotclassmag(self):
        figure(figsize=(8,6))

        for i in range(2):
            subplot(2,1,i+1)
            if i == 0:
                x=self.nb.MAG_APER[:,1]
                y=self.nb.CLASS_STAR
                xl='NB MAG_APER'
            elif i == 1:
                x=self.j.MAG_APER[:,1]
                y=self.j.CLASS_STAR
                xl='J MAG_APER'
            plot(x[self.goodflag],y[self.goodflag],'k.')
            #errorbar(x[self.goodflag],y[self.goodflag],yerr=ery[self.goodflag],fmt=None)
            xlabel(xl)
            ylabel('CLASS_STAR')
            #axhline(y=0,color='k',ls='--')
            #gca().set_xscale('log')
            axis([10, 30,-0.1,1.1])

    def plotfwhmmag(self):
        figure(figsize=(10,4))
        subplots_adjust(left=.1,right=.95,bottom=.15,top=0.95,wspace=.3)

        for i in range(2):
            subplot(1,2,i+1)
            if i == 0:
                x=self.nb.MAG_APER[:,1]
                y=self.nb.FWHM_IMAGE
                xl='NB MAG_APER'
                axl=3.
                dq=3.65
            elif i == 1:
                x=self.j.MAG_APER[:,1]
                y=self.j.FWHM_IMAGE
                xl='J MAG_APER'
                axl=4.
                dq=4.07
            plot(x[self.goodflag],y[self.goodflag],'k.')
            #errorbar(x[self.goodflag],y[self.goodflag],yerr=ery[self.goodflag],fmt=None)
            xlabel(xl)
            ylabel('FWHM_IMAGE')
            axhline(y=axl,color='k',ls='--')
            axhline(y=dq,color='r',ls='--')
            #gca().set_xscale('log')
            axis([10, 25,2,6])
            yticks(arange(2,7))
            y1, y2=gca().get_ylim()
            ax2=twinx()
            ax2.set_ylim(newfirmpixelscale*y1,newfirmpixelscale*y2)


    def plotcolormag(self):
        m1=self.nb.MAG_APER
        m2=self.j.MAG_APER

        y=m2-m1
        figure()
        for i in range(2):
            subplot(2,1,i+1)
            if i == 0:
                x=m1
                xl='NB MAG_APER'
            elif i == 1:
                x=m2
                xl='J MAG_APER'
            plot(x[self.goodflag],y[self.goodflag],'k.')
            xlabel(xl)
            ylabel('J-NB MAG_APER')
            axhline(y=0.35,color='k',ls='--')
            axis([12,30,-2,9])



    def plotcolormagap1(self):
        m1=self.nb.MAG_APER1
        m2=self.j.MAG_APER1
        y=m2-m1
        figure()
        for i in range(2):
            subplot(2,1,i+1)
            if i == 0:
                x=m1
                xl='NB MAG_APER1'
            elif i == 1:
                x=m2
                xl='J MAG_APER1'
            plot(x[self.goodflag],y[self.goodflag],'k.')
            xlabel(xl)
            ylabel('J-NB MAG_APER1')
            axhline(y=0.35,color='k',ls='--')
            axis([12,30,-2,9])

    def plotcolormagap2(self):
        m1=self.nb.MAG_APER2
        m2=self.j.MAG_APER2
        y=m2-m1
        figure()
        for i in range(2):
            subplot(2,1,i+1)
            if i == 0:
                x=m1
                xl='NB MAG_APER2'
            elif i == 1:
                x=m2
                xl='J MAG_APER2'
            plot(x[self.goodflag],y[self.goodflag],'k.')
            xlabel(xl)
            ylabel('J-NB MAG_APER2')
            axhline(y=0.35,color='k',ls='--')
            axis([12,30,-2,9])
cl1317=cluster('RDCSJ1317')
