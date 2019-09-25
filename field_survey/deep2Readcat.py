#!/usr/bin/env python
import pyfits
from pylab import *
import os
import mystuff as my

class deep2:
    def __init__(self):
        infile='/home/rfinn/research/DEEP2/deep2.dr4.groupcat.fits'
        tb=pyfits.open(infile)
        tbdata=tb[1].data
        self.galObjno=tbdata.field('OBJNO')
        self.galRA=tbdata.field('RA')
        self.galDec=tbdata.field('DEC')
        self.galZ=tbdata.field('Z')
        self.galClustnum=tbdata.field('CLUSTNUM')

        tbdata2=tb[2].data
        self.clustnum=tbdata2.field('CLUSTNUM')
        self.field=tbdata2.field('FIELD')
        self.ra=tbdata2.field('RA')
        self.dec=tbdata2.field('DEC')
        self.z=tbdata2.field('Z')
        self.vdisp=tbdata2.field('VDISP')
        self.nmembers=tbdata2.field('N_MEMBERS')
        self.highpurity=tbdata2.field('HIGHPURITY')


        tb.close()
        self.newfirmSample()
        
    def newfirmSample(self):#select those clusters that lie w/in newfirm z window
        centerwave=1.187#center wavelength in microns
        maxwave=1.18825#max wavelength with Transmission > 90%
        minwave = 1.18125#min wavelength w/T>90%
        hawave=.6563 #rest wavelength of Halpha in microns
        self.zmin=minwave/hawave-1.#min redshift where Halpha lies in filter
        self.zmax=maxwave/hawave-1.#min redshift where Halpha lies in filter
        self.filterflag = (self.z > self.zmin) & (self.z < self.zmax)
        self.springflag= (self.ra > 150.) & (self.ra < 260.)
        self.winterflag=(self.ra < 50 ) | (self.ra > 300.)
        self.nmembflag = (self.nmembers > 1)
        self.sigmaflag=(self.vdisp > 250.)
        #flag = self.filterflag & self.springflag & self.nmembflag & self.sigmaflag
        flag = self.filterflag  & self.nmembflag & self.sigmaflag

        self.nfclustnum=self.clustnum[flag]#only clusters that lie in NF filter
        self.nffield=self.field[flag]#only clusters that lie in NF filtera
        self.nfra=self.ra[flag]#only clusters that lie in NF filtera
        self.nfdec=self.dec[flag]#only clusters that lie in NF filtera
        self.nfz=self.z[flag]#only clusters that lie in NF filtera
        self.nfvdisp=self.vdisp[flag]#only clusters that lie in NF filtera
        self.nfnmembers=self.nmembers[flag]#only clusters that lie in NF filtera
        self.nfhighpurity=self.highpurity[flag]#only clusters that lie in NF filtera

        self.nfclustnum=self.clustnum[flag]#only clusters that lie in NF filter

        #individual galaxies
        gflag=(self.galZ < self.zmax) & (self.galZ > self.zmin)
        self.nfgalra=self.galRA[gflag]#only clusters that lie in NF filtera
        self.nfgaldec=self.galDec[gflag]#only clusters that lie in NF filtera
        self.nfgalz=self.galZ[gflag]#only clusters that lie in NF filtera

    def plotNFzdist(self):
        figure()
        hist(self.nfz,20,histtype='step')
        axvline(self.zmin,color='k',ls='--')
        axvline(self.zmax,color='k',ls='--')
        
        title('Redshift Distribution')
        xlabel('Redshift')
        ylabel('Number of Groups')
    def plotNFvdispdist(self):
        figure()
        hist(self.nfvdisp,20,histtype='step')
        title('Velocity Dispersion Distribution')
        xlabel('Velocity Dispersion (km/s)')
        ylabel('Number of Groups')
    def plotNFpositions(self):
        figure()
        subplot(1,2,1)
        plot(self.nfra/15.,self.nfdec,'k.')
        axis([14.2,14.45,51.5,54.])
        subplot(1,2,2)
        plot(self.nfra/15.,self.nfdec,'k.')
        axis([16.78,16.9,34.6,35.5])
        ax=gca()
        text(0,1.05,'z=0.8 Deep2 Group Positions (spring targets)',transform=ax.transAxes,horizontalalignment='center')

        text(-1.4,.5,'Dec (deg)',transform=ax.transAxes,verticalalignment='center',rotation=90)
        text(-.1,-.1,'RA (hr)',transform=ax.transAxes,horizontalalignment='center')

    def plotNFpositionsdeg(self):
        figure(figsize=(12,10))
        dx=.5
        dy=.5
        rot=0
        xy1=array([215.3,53.3,dx,dy,rot],'f')
        
        xy2=array([213.8,52.2,dx,dy,rot],'f')
        xy3=array([253.1,34.97,dx,dy,rot],'f')
        
        subplot(2,2,1)
        scatter(self.nfra,self.nfdec,s=self.nfvdisp/10,color='b',alpha=.7,marker='o')
        plot(self.nfgalra,self.nfgaldec,'k.')
        my.drawbox(xy1,'r:')
        my.drawbox(xy2,'r:')
        axis('equal')
        axis([213.4,216.1,51.6,54.])
        subplot(2,2,2)
        scatter(self.nfra,self.nfdec,s=self.nfvdisp/10,color='b',alpha=0.7,marker='o')
        plot(self.nfgalra,self.nfgaldec,'k.')
        #plot(self.nfra,self.nfdec,'r.')
        my.drawbox(xy3,'r:')
        axis('equal')
        axis([251.4,253.5,34.6,35.5])

        ax=gca()
        text(0,1.05,'z=0.8 Deep2 Group Positions (spring targets)',transform=ax.transAxes,horizontalalignment='center')
        text(-.1,-.1,'RA (hr)',transform=ax.transAxes,horizontalalignment='center')
        text(-1.4,.5,'Dec (deg)',transform=ax.transAxes,verticalalignment='center',rotation=90)

        subplot(2,2,3)
        scatter(self.nfra,self.nfdec,s=self.nfvdisp/10,color='b',alpha=0.7,marker='o')
        plot(self.nfgalra,self.nfgaldec,'k.')
        axis('equal')
        axis([36.4,38.1,0,1.])
        subplot(2,2,4)
        scatter(self.nfra,self.nfdec,s=self.nfvdisp/10,color='b',alpha=0.7,marker='o')
        plot(self.nfgalra,self.nfgaldec,'k.')
        axis('equal')
        axis([351.2,353.7,-.1,.5])

    def plotNFnmembersvdisp(self):
        figure()
        plot(self.nfvdisp,self.nfnmembers,'k.')
        
    def sampleplots(self):
        self.plotNFzdist()
        self.plotNFvdispdist()
        self.plotNFpositionsdeg()


d2=deep2()

        
