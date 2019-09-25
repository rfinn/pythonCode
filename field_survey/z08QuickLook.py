#!/usr/bin/env python
'''
writing code to do a quick analysis of the depth of the 09/2009 data

RXJ1821 looks pretty good

RXJ1716 is compromised by two bright stars on the field

'''
from z08common import *
from pylab import *
import atpy, os
from z08baseCluster import *

class cluster(baseCluster):
    def __init__(self,clustername):
        baseCluster.__init__(self,clustername)
        self.starflag=self.classstar > 0.8
        self.prefix=clustername
    def plotsnrVSmag(self):
        figure()
        plot(self.magbest,self.snr,'k.',label='MAG_BEST')
        plot(self.magauto,self.snrauto,'b.',label='MAG_AUTO')
        axis([18.5,22.,0.,12.])
        axhline(y=10,ls='--',color='k')
        axvline(x=21.6,ls=':',color='k')
        legend(loc='upper right',numpoints=1)
        title(self.prefix+': SNR vs Mag, J-band')
        xlabel('J-band Mag')
        ylabel('S/N Ratio (J-band)')
        savefig(self.prefix+'snrvmagJ.png')
    def plotsnrVSmagNB(self):
        figure()
        plot(self.NBmagbest,self.NBsnr,'k.',label='MAG_BEST')
        plot(self.NBmagauto,self.NBsnrauto,'b.',label='MAG_AUTO')

        axhline(y=10,ls='--',color='k')
        axvline(x=21.6,ls=':',color='k')
        legend(loc='upper right',numpoints=1)
        ylabel('S/N Ratio (Narrow-band)')
        xlabel('NB Magnitude')
        title(self.prefix+': SNR vs Mag, Narrow-band')
        axis([18.5,22.,0.,3.])
        savefig(self.prefix+'snrvmagNB.png')
    def runSextractor2image(self):
        # run sextractor in 2-image mode
        os.system('sex rRXJ1821J.fits -CATALOG_NAME testJ.cat')
        os.system('sex rRXJ1821J.fits, rRXJ1821NB.fits -MAG_ZEROPOINT 18.405 -CATALOG_NAME testNB.cat')
    def readcats(self):
        self.readSextractor('testJ.cat')
        self.readSextractorNB('testNB.cat')
    def plotNBvsJ(self):
        figure()
        keepflag=self.NBmagauto < 30


        errorbar(self.magauto[keepflag],self.NBmagauto[keepflag],yerr=self.NBmagautoerr[keepflag],xerr=self.magautoerr[keepflag],fmt=None,color='k')
        plot(self.magauto[keepflag],self.NBmagauto[keepflag],'k.',label='All')
        plot(self.magauto[keepflag & self.starflag],self.NBmagauto[keepflag & self.starflag],'r.', label='Stars')
        legend(numpoints=1, loc='upper left')
        xlabel('J-band Mag')
        ylabel('Narrow-band Mag')
        xl=arange(10,23)
        plot(xl,xl,'k--')
        axis([10,22,8,30])
        savefig(self.prefix+'NBvJmag.png')
        figure()
        color=self.magauto - self.NBmagauto
        plot(self.magauto[keepflag],color[keepflag],'k.',label='All')
        plot(self.magauto[keepflag & self.starflag],color[keepflag & self.starflag],'r.', label='Stars')
        legend(numpoints=1, loc='upper left')
        errorbarflag=keepflag & (self.NBsnrauto > 3)
        errorbar(self.magauto[errorbarflag],color[errorbarflag],fmt=None,ecolor='k')
        title(self.prefix+': J-NB vs J')
        xlabel('J Mag')
        ylabel('J - NB')
        axhline(y=0,color='k',ls='--')
        savefig(self.prefix+'ColorvJmag.png')
    def plotclassstarJ(self):
        figure()
        plot(self.magauto,self.classstar,'k.')
        xlabel('J-band Mag Auto')
        ylabel('SE Class Star (1=Star)')
        title('SE Classifier Index vs J')
        axhline(y=0.8,color='r',ls='--')
        savefig(self.prefix+'StarClassvJmag.png')
rxj18=cluster('RXJ1821')
#rxj18.runSextractor2image()
rxj18.plotsnrVSmag()
rxj18.plotsnrVSmagNB()
rxj18.plotNBvsJ()
rxj18.plotclassstarJ()
#rxj18.readSextractor('test.cat')
#rxj18.plotsnrVSmag()
