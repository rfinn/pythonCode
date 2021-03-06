#!/usr/bin/env python

'''
GOAL:
  compare the results from sextractor photometry with the results from wise

'''

import pyfits
from LCScommon import *
from pylab import *
import os
import mystuff as my
from LCSReadmasterBaseNSA import *
import matplotlib.cm as cm
import chary_elbaz_24um as chary
Lsol=3.826e33#normalize by solar luminosity
bellconv=9.8e-11#converts Lir (in L_sun) to SFR/yr
bellconv=4.5e-44#Kenn 98 conversion fro erg/s to SFR/yr

#usefwhm24=1
fieldColor='0.7'
Remin=1.
mypath=os.getcwd()
if mypath.find('Users') > -1:
    print "Running on Rose's mac pro"
    homedir='/Users/rfinn/'
elif mypath.find('home') > -1:
    print "Running on coma"
    homedir='/home/rfinn/'

figuredir=homedir+'Dropbox/Research/LocalClusters/SamplePlots/'
figuredir=homedir+'research/LocalClusters/SamplePlots/'


class cluster(baseClusterNSA):
    def __init__(self,clustername):
        baseClusterNSA.__init__(self,clustername)
        mypath=os.getcwd()
        if mypath.find('Users') > -1:
            print "Running on Rose's mac pro"
            infile='/Users/rfinn/research/LocalClusters/MasterTables/'+clustername+'mastertable.WithProfileFits.fits'
        elif mypath.find('home') > -1:
            print "Running on coma"
            infile=homedir+'research/LocalClusters/MasterTables/'+clustername+'mastertable.WithProfileFits.fits'

        self.mipssnrflag = self.mipssnr > 6.
        try:
            self.readsnr24NSA()
        except:
            print self.prefix,": couldn't read SNR24 file"
        try:
            self.readGalfitSersicResults()
        except:
            print self.prefix,": couln't read galfit sersic results"
        try:
            self.readGalfitResults()
        except:
            print self.prefix,": couldn't read galfit results"
        self.fwhm24=self.sex24.FWHM_DEG*3600.
        self.fwhm24=self.sex24.FLUX_RADIUS1*mipspixelscale
        self.member=self.membflag
        self.member=self.dvflag
        self.sample24flag= self.sex24.MATCHFLAG24 & ~self.agnflag & (self.n.SERSIC_TH50 > Remin) & (log10(self.stellarmass) > 9.5) & (log10(self.stellarmass) < 12)
        self.blueclustersample=self.member & self.blueflag & self.sample24flag
        self.bluefieldsample=~self.member & self.blueflag & self.sample24flag
        self.greenclustersample=self.member & self.greenflag & self.sample24flag
        self.greenfieldsample=~self.member & self.greenflag & self.sample24flag
        self.redclustersample=self.member & self.redflag & self.sample24flag
        self.redfieldsample=~self.member & self.redflag & self.sample24flag
        self.varlookup={'stellarmass':log10(self.stellarmass),'Re':self.n.SERSIC_TH50,'R24':self.sex24.FLUX_RADIUS1*mipspixelscale,'NUV':self.n.ABSMAG[:,1],'r':self.n.ABSMAG[:,4],'m24':self.sex24.MAG_BEST,'redshift':self.n.ZDIST,'NUVr':(self.n.ABSMAG[:,1]-self.n.ABSMAG[:,4]),'NUV24':(self.n.ABSMAG[:,1]-self.sex24.MAG_BEST),'24mass':(self.sex24.MAG_BEST-log10(self.stellarmass)),'ratioR':self.sex24.FLUX_RADIUS1*mipspixelscale/self.n.SERSIC_TH50}

        self.fw4=8.363*10.**(-1.*self.wise.W4MAG_3/2.5)*1.e6 # uJy
        self.fw4err=8.363*10.**(-1.*(self.wise.W4MAG_3 - self.wise.W4SIGM_3)/2.5)*1.e6 -self.fw4# uJy


    def plotmipswise(self,plotsingle=1):
        #zeropoint fluxes from http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#conv2ab
        #self.fw1=309.54*10.**(-1.*wtable.w1mpro/2.5)
        #self.fw2=171.787*10.**(-1.*wtable.w2mpro/2.5)
        #self.fw3=31.674*10.**(-1.*wtable.w3mpro/2.5)
        #self.fw4=8.363*10.**(-1.*wtable.w4mpro/2.5)
        if plotsingle:
            figure()
        flag=self.n.On24ImageFlag_1 & self.sex24.MATCHFLAG24
        x=self.sex24.FLUX_BEST[flag]*mipsflux2umJyconv/1e3 # convert to mJy
        y=self.fw4[flag]/1e3
        plot(x,y,'ko',markersize=4)
        xl=arange(-1,2.5,.1)
        xl=10.**xl
        plot(xl,xl,'k--')
        axhline(y=5.1,color='k',ls=':')
        gca().set_yscale('log')
        gca().set_xscale('log')
        if plotsingle:
            xlabel('MIPS 24um flux (mJy)')
            ylabel('WISE 22um flux (mJy)')
        ax=gca()
        axis([15/1e3,5.e2,.2,5.e2])

        text(.1,.85,'$'+self.clustername+'$',fontsize=14,transform=ax.transAxes,horizontalalignment='left')

    def plotmipssnrflux(self,plotsingle=1):
        if plotsingle:
            figure()
        flag=self.n.On24ImageFlag_1 & self.sex24.MATCHFLAG24
        x=self.sex24.FLUX_BEST[flag]*mipsflux2umJyconv/1e3 # convert to mJy
        xerr=self.sex24.FLUXERR_BEST[flag]*mipsflux2umJyconv/1e3 # convert to mJy
        y=abs(x/xerr)
        plot(x,y,'ro',markersize=4)
        y=abs(self.fw4/self.fw4err)
        plot(self.fw4/1e3,y,'ko',mfc='None',markersize=4)
        xl=arange(-1,2.,.1)
        xl=10.**xl
        #plot(xl,xl,'k--')
        axvline(x=5.1,color='k',ls=':')
        axvline(x=2,color='r',ls='--')
        axhline(y=5.,color='k',ls='--')
        axhline(y=3.,color='k',ls='--')
        axhline(y=10.,color='b',ls='--')
        gca().set_yscale('log')
        gca().set_xscale('log')
        if plotsingle:
            xlabel('MIPS 24um flux (mJy)')
            ylabel('SNR')
        ax=gca()
        axis([15/1e3,5.e2,.5,90])

        text(.1,.85,'$'+self.clustername+'$',fontsize=14,transform=ax.transAxes,horizontalalignment='left')

        
    def plotF24hist(self):
        mybins=arange(0,5000,100)
	y=hist(self.mipsflux[self.apexflag],bins=mybins,histtype='step')
	ngal=y[0]
	x=y[1]
	xbin=zeros(len(ngal))
	for i in range(len(xbin)):
            xbin[i]=0.5*(x[i]+x[i+1])
	#clf()                                                                                       
        self.xbin=xbin
        self.ngal=ngal
	plot(xbin,ngal,'ro')
	errorbar(xbin,ngal,sqrt(ngal),fmt=None)
        axis([0,2000,0,20])
        title(self.clustername)
	ax=gca()
	#ax.set_yscale('log')
	#ax.set_xscale('log')
        #xlabel('24um Flux')

    def plotL24hist(self):
        #figure(1)
        
	y=hist(log10(self.L24[self.apexflag]),bins=10,histtype='stepfilled')
	ax=gca()
	#ax.set_yscale('log')
        axis([6.5,11.,1,22])
        xmin,xmax=xlim()
        xticks(arange(round(xmin),xmax+1,1,'i'),fontsize=10)
        ymin,ymax=ylim()
        title(self.clustername)
        #yticks(arange(round(ymin),ymax+1,2,'i'),fontsize=10)


def plotF24histall():
    figure(figsize=[10,10])
    clf()
    subplots_adjust(wspace=.25,hspace=.35)
    i=0
    for cl in mylocalclusters:
        i += 1
        subplot(3,3,i)
        cl.plotF24hist()
    ax=gca()
    text(-.75,-.35,'$log_{10}(F_{24} \  (\mu Jy))$',fontsize=22,horizontalalignment='center',transform=ax.transAxes)
    subplot(3,3,4)
    text(-2.8,1.9,'$N_{gal}$',fontsize=22,verticalalignment='center',rotation=90,transform=ax.transAxes)
    savefig(figuredir+'PlotF24histAll.png')


def plotwisemips24all():
    figure(figsize=[10,8])
    clf()
    subplots_adjust(wspace=.02,hspace=.02)
    i=0
    for cl in mylocalclusters:
        i +=1
        subplot(3,3,i)
        print cl.prefix
        cl.plotmipswise(plotsingle=0)
        multiplotaxes(i)

    multiplotlabels('$MIPS \ F_{24} \ (mJy) $','$ WISE \ F_{22} \ (mJy) $')
    savefig(figuredir+'PlotWiseMipsall.eps')
def plotmipssnrfluxall():
    figure(figsize=[10,8])
    clf()
    subplots_adjust(wspace=.02,hspace=.02)
    i=0
    for cl in mylocalclusters:
        i +=1
        subplot(3,3,i)
        print cl.prefix
        cl.plotmipssnrflux(plotsingle=0)
        multiplotaxes(i)

    multiplotlabels('$MIPS \ F_{24} \ (mJy) $','$ SNR_{24} $')
    savefig(figuredir+'PlotWiseMipsall.eps')


mkw11=cluster('MKW11')
coma=cluster('Coma')

mkw8=cluster('MKW8')
awm4=cluster('AWM4')
a2052=cluster('A2052')
a2063=cluster('A2063')
ngc=cluster('NGC6107')
herc=cluster('Hercules')
a1367=cluster('A1367')
clustersbymass=[mkw8,mkw11,ngc,awm4,a2052,a2063,herc,a1367,coma]
clustersbydistance=[a1367,coma,mkw11,mkw8,ngc,awm4,a2063,a2052]
clustersbylx=[mkw11,ngc,mkw8,awm4,herc,a1367,a2063,a2052,coma]
mylocalclusters=clustersbymass
