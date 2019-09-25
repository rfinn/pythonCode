#!/usr/bin/env python
import numarray as N
from pylab import *

import matplotlib
matplotlib.rc('xtick',labelsize=12)     # fontsize of the tick labels
matplotlib.rc('ytick',labelsize=12)     # fontsize of the tick labels

class cluster:	    	
    def __init__(self,prefix):
	self.prefix=prefix

	file='/Users/rfinn/clusters/spitzer/MasterTables/'+str(prefix)+'mosaic_extract_final.tbl'

        input=open(file,'r')
        #get number of galaxies

        ngal=0
	for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue
            ngal=ngal+1
        input.close()


        self.id24 = N.zeros(ngal,'f')
        self.imagex24 = N.zeros(ngal,'f')
        self.imagey24  = N.zeros(ngal,'f')
        self.ra24 = N.zeros(ngal,'f')
        self.dec24 = N.zeros(ngal,'f')
        self.f24 = N.zeros(ngal,'d')#flux
        self.errf24 = N.zeros(ngal,'d')
        self.fap1 = N.zeros(ngal,'d')#flux in aperture 1 (1,1.5,2,2.6,3,3.5,4,4.5,5.,5.5) pixels
        self.fap2 = N.zeros(ngal,'d')#flux
        self.fap3 = N.zeros(ngal,'d')#flux
        self.fap4 = N.zeros(ngal,'d')#flux in ap 4 - this is one w/ap cor of 1.67 (Calzetti et al 2007)
        self.fap5 = N.zeros(ngal,'d')#flux
        self.fap6 = N.zeros(ngal,'d')#flux
        self.fap7 = N.zeros(ngal,'d')#flux
        self.fap8 = N.zeros(ngal,'d')#flux
        self.fap9 = N.zeros(ngal,'d')#flux
        self.fap10 = N.zeros(ngal,'d')#flux
        self.errfap1 = N.zeros(ngal,'d')#flux in aperture 1 (1,1.5,2,2.6,3,3.5,4,4.5,5.,5.5) pixels
        self.errfap2 = N.zeros(ngal,'d')#flux
        self.errfap3 = N.zeros(ngal,'d')#flux
        self.errfap4 = N.zeros(ngal,'d')#flux in ap 4 - this is one w/ap cor of 1.67 (Calzetti et al 2007)
        self.errfap5 = N.zeros(ngal,'d')#flux
        self.errfap6 = N.zeros(ngal,'d')#flux
        self.errfap7 = N.zeros(ngal,'d')#flux
        self.errfap8 = N.zeros(ngal,'d')#flux
        self.errfap9 = N.zeros(ngal,'d')#flux
        self.errfap10 = N.zeros(ngal,'d')#flux
        self.snr24 = N.zeros(ngal,'d')#SNR calculated by mopex
        self.ndeblend = N.zeros(ngal,'f')#SNR calculated by mopex


        input=open(file,'r')
        i=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue

            t=line.split()
            #(self.id24[i],self.imagex24[i],self.imagey24[i],self.ra24[i],self.dec24[i],self.f24[i],self.errf24[i])=(float(t[0]),float(t[1]),float(t[2]),float(t[3]),float(t[4]),float(t[5]),float(t[6]))
	    #print i,'length of t = ',len(t)
	    #print t
	    (self.id24[i],self.ndeblend[i],self.ra24[i],self.dec24[i],self.imagex24[i],self.imagey24[i],self.f24[i],self.errf24[i],self.snr24[i],self.fap1[i],self.fap2[i],self.fap3[i],self.fap4[i],self.fap5[i],self.fap6[i],self.fap7[i],self.fap8[i],self.fap9[i],self.fap10[i],self.errfap1[i],self.errfap2[i],self.errfap3[i],self.errfap4[i],self.errfap5[i],self.errfap6[i],self.errfap7[i],self.errfap8[i],self.errfap9[i],self.errfap10[i])=(float(t[0]),float(t[2]),float(t[3]),float(t[5]),float(t[8]),float(t[10]),float(t[13]),float(t[14]),float(t[18]),float(t[25]),float(t[26]),float(t[27]),float(t[28]),float(t[29]),float(t[30]),float(t[31]),float(t[32]),float(t[33]),float(t[34]),float(t[45]),float(t[46]),float(t[47]),float(t[48]),float(t[49]),float(t[50]),float(t[51]),float(t[52]),float(t[53]),float(t[54]))
            i=i+1
        input.close()#44 -> 43


	print prefix,' number of galaxies = ',i,len(self.f24),len(self.fap4)



def makeplot():
    subplot(4,4,1)
    c=c1
    makeplotsub(c)
    subplot(4,4,2)
    c=c2
    makeplotsub(c)
    subplot(4,4,3)
    c=c3
    makeplotsub(c)
    subplot(4,4,4)
    c=c4
    makeplotsub(c)
    subplot(4,4,5)
    c=c5
    makeplotsub(c)
    subplot(4,4,6)
    c=c6
    makeplotsub(c)
    subplot(4,4,7)
    #c=c7
    #makeplotsub(c)
    #subplot(4,4,8)
    c=c8
    makeplotsub(c)
    subplot(4,4,8)
    c=c9
    makeplotsub(c)
    #savefig('fluxcompalla.eps')

    #cla()
    #clf()
    subplot(4,4,9)
    c=c10
    makeplotsub(c)
    subplot(4,4,10)
    c=c11
    makeplotsub(c)
    subplot(4,4,11)
    c=c12
    makeplotsub(c)
    subplot(4,4,12)
    c=c13
    makeplotsub(c)
    subplot(4,4,13)
    c=c14
    makeplotsub(c)
    subplot(4,4,14)
    c=c15
    makeplotsub(c)
    subplot(4,4,15)
    c=c16
    makeplotsub(c)
    subplot(4,4,16)
    c=c17
    makeplotsub(c)
    ax=gca()
    s='MOPEX Flux'
    text(-1.5,-.5,s,fontsize=32,horizontalalignment='center',transform=ax.transAxes)
    s='Ap flux * 1.67'
    text(-4.2,2.,s,fontsize=32,verticalalignment='center',rotation='vertical',transform=ax.transAxes)
    savefig('fluxcompall.eps')



def makeplotsub(c):
    plot(c.f24,c.fap4*1.67,'bo',markersize=4)
    #plot(c.f24,c.fap4,'ro',markersize=4)
    #plot(c.f24,c.fap3,'go',markersize=2)
    #plot(c.f24,c.fap2,'yo',markersize=2)
    #plot(c.f24,c.fap1,'co',markersize=2)
    x=N.arange(1.,max(c.f24),100.)
    y=x
    plot(x,y,'k-')
    y2=y+N.log10(2.)
    plot(x,y2,'k--')
    y2=y+N.log10(.5)
    plot(x,y2,'k--')
    y=3.*x
    #plot(x,y,'k-')
    y=4.*x
    #plot(x,y,'k-')
    y=5.*x
    #plot(x,y,'k-')
    ax=gca()
    text(.1,.8,c.prefix,horizontalalignment='left',verticalalignment='center',transform=ax.transAxes)
    #xlabel('Mopex F(24)',fontsize=fsize)
    #ylabel('Ap Flux ',fontsize=fsize)
    axis([20.,4000.,20.,4000.])
    ax.set_xscale('log')
    ax.set_yscale('log')


def oldstuff():
    subplot(442)
    plot(fap1,fap2,'co')
    xlabel('fap1',fontsize=fsize)
    ylabel('fap2',fontsize=fsize)
    subplot(443)
    plot(fap1,fap4,'co')
    xlabel('fap1',fontsize=fsize)
    ylabel('fap4',fontsize=fsize)
    subplot(444)
    plot(fap2,fap4,'co')
    xlabel('fap2',fontsize=fsize)
    ylabel('fap4',fontsize=fsize)
    subplot(445)
    r=fap3/fap2
    plot(fap2,r,'co')
    xlabel('fap2',fontsize=fsize)
    ylabel('3/2',fontsize=fsize)
    subplot(446)
    r=fap4/fap2
    plot(fap2,r,'co')
    xlabel('fap2',fontsize=fsize)
    ylabel('4/2',fontsize=fsize)
    subplot(449)
    r=fap4/fap3
    plot(fap3,r,'co')
    xlabel('fap2',fontsize=fsize)
    ylabel('4/3',fontsize=fsize)
    
    r=array([1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5],'f')
    
    s=str(prefix)+'fluxcomp.eps'
    savefig(s)


c1=cluster('cl1018')
c2=cluster('cl1037')
c3=cluster('cl1040')
c4=cluster('cl105411')
c5=cluster('cl105412')
c6=cluster('cl1059')
#c7=cluster('cl1103')
c8=cluster('cl1138')
c9=cluster('cl1202')
c10=cluster('cl1216')
c11=cluster('cl1227')
c12=cluster('cl1232')
c13=cluster('cl1301')
c14=cluster('cl1353')
c15=cluster('cl1354')
c16=cluster('cl1411')
c17=cluster('cl1420')

makeplot()
