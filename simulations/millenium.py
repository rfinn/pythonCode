#!/usr/bin/env python
#"""
#useage:
#millenium.py switch1 switch2
#switch1
#   0=calculate quantities
#   1=update plots only
#switch2
#   0=use mini millenium catalog
#   1=use full millenium catalog
#"""
import sys
import sets
import copy
import os
import Numeric as N
import RandomArray as Rand
#import scipy
#import pylab #use matplotlib instead of scipy pylab.std
import mystuff as my
import ppgplot
import time
import matplotlib
from matplotlib import rc
rc('font',family='serif', style='normal', variant='normal',weight='bold', stretch='normal', size='large')
import pylab

starttime=time.clock()
print "start time = ",starttime

schdef=1.2
slwdef=4.
omega0=0.3
omegaL=0.7
h=.7
H0=h*100.#km/s/Mpc

G=4.31e-9#gravitational constant in (km/s)^2 Mpc Msolar^-1
cl=3.e5#speed of light in km/s
Mclmin=70.#min cluster mass in units of 10^12 Msun

mabscut=-20.38
DH=cl/H0
#nr=1. #number of virial radii away from cluster 
#nv=3. #number of sigma away from cluster
#read in catalog of halos, keeping those w/M>=Mclmin
def pylabsubplot1(x,y05,y1,y2,fx,fy05,fy1,fy2):
    nbin=5
    nbinscale=3
    xl=N.arange(100.,5000.,100.)
    #print 'x = ',x
    #print 'y = ',y05
    

    pylab.cla()
    pylab.clf()
    xminsig=2.6
    xmaxsig=3.1
    symheight=3.
    sigmamin=400.
    sigmamax=1200.
    nbin=7

    ymin=0.
    ymax=2.5
    xmin=xminsig
    xmax=xmaxsig

    pylab.subplots_adjust(left=0.125, right=.9,bottom=.1,top=0.9,wspace=.3,hspace=.1)
    p1=pylab.subplot(311)
    (xbin,ybin,ybinerr)=my.binitbins(sigmamin,sigmamax,nbin,x,y05)
    (fxbin,fybin,fybinerr)=my.binitbins(sigmamin,sigmamax,nbin,fx,fy05)
    ybin=ybin
    ybinerr=ybinerr
    fybin=fybin
    fybinerr=fybinerr
    try:
	scale=ybin[nbinscale]/fybin[nbinscale]
    except ZeroDivisionError:
	scale=1.
    fybin=fybin*scale
    fybinnerr=fybinerr*scale
    clabel="Control*%3.1f"%(scale)
    pylab.plot(xbin,ybin,'ko',label="Cluster",markersize=6.)
    pylab.plot(fxbin,fybin,'wo',label=clabel,markersize=6.)
    pylab.hold(True)
    t=pylab.legend(loc='upper left')
    pylab.errorbar(fxbin,fybin,fybinerr,fmt=None,ecolor='k')
    pylab.errorbar(xbin,ybin,ybinerr,fmt='ko')
    c1=ybin[nbinscale]/(xbin[nbinscale]/1000.)**3
    yl=c1*(xl/1000.)**3 
    pylab.plot(xl,yl,'k-')

    ymin=.5*min(fybin)
    ymax=1.5*max(fybin)
    pylab.axis([400.,1300.,ymin,ymax])
    ax=pylab.gca()
    ax.set_yscale("log")
    ax.set_xscale("log")
    pylab.setp(ax,xticklabels=[])

    pylab.text(.43,.85,r'$0.5\times R_{200}$',horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=20)
    pylab.subplot(312)
    (xbin,ybin,ybinerr)=my.binitbins(sigmamin,sigmamax,nbin,x,y1)
    (fxbin,fybin,fybinerr)=my.binitbins(sigmamin,sigmamax,nbin,fx,fy1)
    #print 'xbin=',xbin
    #print 'ybin=',ybin 
    ybin=ybin
    ybinerr=ybinerr
    fybin=fybin
    fybinerr=fybinerr
    try:
	scale=ybin[nbinscale]/fybin[nbinscale]
    except ZeroDivisionError:
	scale=1.


    fybin=fybin*scale
    fybinnerr=fybinerr*scale
    clabel="Control*%3.1f"%(scale)
    pylab.plot(xbin,ybin,'ko',label="Cluster",markersize=6.)
    pylab.plot(fxbin,fybin,'wo',label=clabel,markersize=6.)
    pylab.hold(True)
    pylab.legend(loc='upper left')
    pylab.errorbar(fxbin,fybin,fybinerr,fmt=None,ecolor='k')
    pylab.errorbar(xbin,ybin,ybinerr,fmt='ko')
    c1=ybin[nbinscale]/(xbin[nbinscale]/1000.)**3
    yl=c1*(xl/1000.)**3 
    pylab.plot(xl,yl,'k-')

    ymin=.5*min(fybin)
    ymax=1.5*max(fybin)
    pylab.axis([400.,1300.,ymin,ymax])

    ax=pylab.gca()
    ax.set_yscale("log")
    ax.set_xscale("log")
    pylab.setp(ax,xticklabels=[])
    pylab.text(.43,.85,r'$1.0\times R_{200}$',horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=20)

    pylab.subplot(313)

    (xbin,ybin,ybinerr)=my.binitbins(sigmamin,sigmamax,nbin,x,y2)
    (fxbin,fybin,fybinerr)=my.binitbins(sigmamin,sigmamax,nbin,fx,fy2)
    ybin=ybin
    ybinerr=ybinerr
    fybin=fybin
    fybinerr=fybinerr
    try:
	scale=ybin[nbinscale]/fybin[nbinscale]
    except ZeroDivisionError:
	scale=1.
    fybin=fybin*scale
    fybinnerr=fybinerr*scale
    clabel="Control*%3.1f"%(scale)
    pylab.plot(xbin,ybin,'ko',label="Cluster",markersize=6.)
    pylab.plot(fxbin,fybin,'wo',label=clabel,markersize=6.)
    pylab.hold(True)
    pylab.legend(loc='upper left')
    pylab.errorbar(fxbin,fybin,fybinerr,fmt=None,ecolor='k')
    pylab.errorbar(xbin,ybin,ybinerr,fmt='ko')
    c1=ybin[nbinscale]/(xbin[nbinscale]/1000.)**3
    yl=c1*(xl/1000.)**3 

    pylab.plot(xl,yl,'k-')
    ymin=.5*min(fybin)
    ymax=1.5*max(fybin)
    pylab.axis([400.,1300.,ymin,ymax])

    ticks=N.array([400.,600.,1000.],'d')
    ax=pylab.gca()
    ax.set_yscale("log")
    ax.set_xscale("log",)
    pylab.text(.43,.85,r'$2.0\times R_{200}$',horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=20)
    pylab.xlabel(r'$\sigma \ \rm{(km/s)}$',fontsize=24)
    pylab.setp(ax,xticklabels=['1000.'])

def pylabsubplot2(x,y05,y1,y2,fx,fy05,fy1,fy2):
    nbin=5
    nbinscale=3
    xl=N.arange(100.,5000.,100.)
    #print 'x = ',x
    #print 'y = ',y05
    

    pylab.cla()
    pylab.clf()
    xminsig=2.6
    xmaxsig=3.1
    symheight=3.
    sigmamin=400.
    sigmamax=1200.
    nbin=7

    ymin=0.
    ymax=2.5
    xmin=xminsig
    xmax=xmaxsig

    pylab.subplots_adjust(left=0.125, right=.9,bottom=.1,top=0.9,wspace=.3,hspace=.1)
    p1=pylab.subplot(311)
    (xbin,ybin,ybinerr)=my.binitbins(sigmamin,sigmamax,nbin,x,y05)
    (fxbin,fybin,fybinerr)=my.binitbins(sigmamin,sigmamax,nbin,fx,fy05)
    ybin=ybin
    ybinerr=ybinerr
    fybin=fybin
    fybinerr=fybinerr
    scale=ybin[nbinscale]#/fybin[nbinscale]
    fybin=fybin*scale
    fybinnerr=fybinerr*scale
    clabel="Control*%3.1f"%(scale)
    pylab.plot(xbin,ybin,'ko',label="Cluster",markersize=6.)
    pylab.plot(fxbin,fybin,'wo',label=clabel,markersize=6.)
    pylab.hold(True)
    t=pylab.legend(loc='upper left')
    #pylab.errorbar(fxbin,fybin,fybinerr,fmt=None,ecolor='k')
    #pylab.errorbar(xbin,ybin,ybinerr,fmt='ko')
    #c1=ybin[nbinscale]/(xbin[nbinscale]/1000.)**3
    #yl=c1*(xl/1000.)**3 
    #pylab.plot(xl,yl,'k-')

    #ymin=.5*min(fybin)
    #ymax=1.5*max(fybin)
    ymin=0.
    ymax=50.
    pylab.axis([400.,1300.,ymin,ymax])
    #ax=pylab.gca()
    #ax.set_yscale("log")
    #ax.set_xscale("log")
    #pylab.setp(ax,xticklabels=[])

    #pylab.text(.43,.85,r'$0.5\times R_{200}$',horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=20)
    pylab.subplot(312)
    (xbin,ybin,ybinerr)=my.binitbins(sigmamin,sigmamax,nbin,x,y1)
    (fxbin,fybin,fybinerr)=my.binitbins(sigmamin,sigmamax,nbin,fx,fy1)
    ybin=ybin
    ybinerr=ybinerr
    fybin=fybin
    fybinerr=fybinerr
    scale=ybin[nbinscale]#/fybin[nbinscale]
    fybin=fybin*scale
    fybinnerr=fybinerr*scale
    clabel="Control*%3.1f"%(scale)
    pylab.plot(xbin,ybin,'ko',label="Cluster",markersize=6.)
    pylab.plot(fxbin,fybin,'wo',label=clabel,markersize=6.)
    pylab.hold(True)
    pylab.legend(loc='upper left')
    #pylab.errorbar(fxbin,fybin,fybinerr,fmt=None,ecolor='k')
    #pylab.errorbar(xbin,ybin,ybinerr,fmt='ko')
    #c1=ybin[nbinscale]/(xbin[nbinscale]/1000.)**3
    #yl=c1*(xl/1000.)**3 
    #pylab.plot(xl,yl,'k-')

    ymin=.5*min(fybin)
    ymax=1.5*max(fybin)
    ymin=0.
    ymax=50.

    pylab.axis([400.,1300.,ymin,ymax])

    #ax=pylab.gca()
    #ax.set_yscale("log")
    #ax.set_xscale("log")
    #pylab.setp(ax,xticklabels=[])
    #pylab.text(.43,.85,r'$1.0\times R_{200}$',horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=20)

    pylab.subplot(313)

    (xbin,ybin,ybinerr)=my.binitbins(sigmamin,sigmamax,nbin,x,y2)
    (fxbin,fybin,fybinerr)=my.binitbins(sigmamin,sigmamax,nbin,fx,fy2)
    ybin=ybin
    ybinerr=ybinerr
    fybin=fybin
    fybinerr=fybinerr
    scale=ybin[nbinscale]#/fybin[nbinscale]
    fybin=fybin*scale
    fybinnerr=fybinerr*scale
    clabel="Control*%3.1f"%(scale)
    pylab.plot(xbin,ybin,'ko',label="Cluster",markersize=6.)
    pylab.plot(fxbin,fybin,'wo',label=clabel,markersize=6.)
    pylab.hold(True)
    pylab.legend(loc='upper left')
    #pylab.errorbar(fxbin,fybin,fybinerr,fmt=None,ecolor='k')
    #pylab.errorbar(xbin,ybin,ybinerr,fmt='ko')
    #c1=ybin[nbinscale]/(xbin[nbinscale]/1000.)**3
    #yl=c1*(xl/1000.)**3 

    #pylab.plot(xl,yl,'k-')
    ymin=.5*min(fybin)
    ymax=1.5*max(fybin)
    ymin=0.
    ymax=50.

    pylab.axis([400.,1300.,ymin,ymax])

    ticks=N.array([400.,600.,1000.],'d')
    #ax=pylab.gca()
    #ax.set_yscale("log")
    #ax.set_xscale("log",)
    #pylab.text(.43,.85,r'$2.0\times R_{200}$',horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=20)
    #pylab.xlabel(r'$\sigma \ \rm{(km/s)}$',fontsize=24)
    #pylab.setp(ax,xticklabels=['1000.'])

def plotstellarsigmapylab():
    print "stellar mass"
    x=(c.sigma)
    y05=c.ostellmass05/(10.)
    y1=c.ostellmass1/(10.)
    y2=c.ostellmass2/(10.)
    fx=(c.sigma)
    fy05=c.stellmass05/(10.)
    fy1=c.stellmass1/(10.)
    fy2=c.stellmass2/(10.)

    pylabsubplot1(x,y05,y1,y2,fx,fy05,fy1,fy2)
    pylab.subplot(312)
    pylab.ylabel(r'$\rm{\Sigma \ M_* \  (10^{11} M_\odot)}$',fontsize=24)
    pylab.savefig('stellmassmembsigma3.eps')

def plotsfrsigmapylab():
    print "SFR"
    x=(c.sigma)
    y05=c.osfr05#sumitmemb(3.,0.5,0)
    y1=c.osfr1#sumitmemb(3.,1.,0)
    y2=c.osfr2#sumitmemb(3.,2.,0)
    fx=(c.sigma)
    fy05=c.sfr05#sumitmemb(3.,0.5,0)
    fy1=c.sfr1#sumitmemb(3.,1.,0)
    fy2=c.sfr2#sumitmemb(3.,2.,0)
    #for i in range(len(c.sigma)):
	#print i, c.osfr1[i],c.sfr1[i],c.ongal1[i],c.ngal1[i]
    pylabsubplot1(x,y05,y1,y2,fx,fy05,fy1,fy2)

    pylab.subplot(312)
    pylab.ylabel(r'$\rm{\Sigma SFR \  (M_\odot/yr)}$',fontsize=24)
    pylab.savefig("sfrmembsigma3.eps")

def plotngalsigmapylab():
    print "Ngal"
    #pylab.show(True)
    x=(c.sigma)
    y05=c.ongal05#sumitmemb(3.,0.5,0)
    y1=c.ongal1#sumitmemb(3.,1.,0)
    y2=c.ongal2#sumitmemb(3.,2.,0)
    fx=(c.sigma)
    fy05=c.ngal05#sumitmemb(3.,0.5,0)
    fy1=c.ngal1#sumitmemb(3.,1.,0)
    fy2=c.ngal2#sumitmemb(3.,2.,0)
    pylabsubplot1(x,y05,y1,y2,fx,fy05,fy1,fy2)

    pylab.subplot(312)
    pylab.ylabel(r'$\rm{N_{gal}}$',fontsize=24,fontweight='bold')
    pylab.savefig("ngalmembsigma3.eps")


def plotdVdz():
    nv=3.
    nr=1.
    ppgplot.pgbeg("dVdz.ps/vcps",1,1)  #color port.
    ppgplot.pgpap(8.,1.25)
    ppgplot.pgpage
    ppgplot.pgsch(1.2) #font size
    ppgplot.pgslw(3)   #line width

    # 1st panel with symbols w/ stddev errorbars


    x1=.15
    x2=.45
    x3=.6
    x4=.95
    y1=.15
    y2=.425
    y3=.575
    y4=.85
    xlabel=14.1-14.
    ylabel=1.15
    schdef=1.2
    slwdef=4
    ppgplot.pgsch(schdef)
    xmin=0.
    xmax=1.1
    ymin=0.
    ymax=1.2

    ppgplot.pgsvp(x1,x4,y1,y4)  #sets viewport
    ppgplot.pgslw(slwdef)   #line width
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('bcnst',.2,2,'bcvnst',.2,2)  #tickmarks and labeling
    ppgplot.pgmtxt('b',2.5,0.5,0.5,"z")    #xlabel
    ppgplot.pgmtxt('l',2.6,0.5,0.5,"(1/DH)\u3\d c dV\dc\u/dv/d\gW")

    z=N.arange(0.,5.,.1)
    beta=((1+z)**2 - 1)/((1+z)**2 + 1)
    dV=N.zeros(len(z),'d')
    for i in range(len(z)):
        #dz=dv/(1+z[i])*(1- ((1+z[i])**2 -1)/((1+z[i])**2+1))**(-2)
        #z1=z[i]-0.5*dz
        #z2=z[i]+0.5*dz
        #dV[i]=my.dL(z2,h) - my.dL(z1,h)
        dA=my.DA(z[i],h)*206264./1000.
        dV[i]=DH*(1+z[i])*(dA)**2/(my.E(z[i]))/(1-beta[i])**2/DH**3
        #dV[i]=DH*(1+z[i])**2*(dA)**2/(my.E(z[i]))/DH**3#for comparison w/Hogg
        if z[i] < 1:
            print i,z[i],dV[i],dV[i]**(1./3.)
            
    ppgplot.pgline(z,dV) 

    ppgplot.pgend()





def plotngalsigmaradcuts():
    nr=1.
    nv=3.
    bbJmax=-18.
    ppgplot.pgbeg("ngalmhalo-radcut.ps/vcps",1,1)  #color port.
    ppgplot.pgpap(8.,1.25)
    ppgplot.pgpage
    ppgplot.pgsch(1.2) #font size
    ppgplot.pgslw(3)   #line width

    # 1st panel with symbols w/ stddev errorbars

    str1="R\dp\u < "
    str2=" R\dv\u"
    x1=.1
    x2=.45
    x3=.6
    x4=.95
    y1=.15
    y2=.425
    y3=.575
    y4=.85
    xlabel=14.25-14.
    ylabel=1.14
    ppgplot.pgsvp(x1,x2,y3,y4)  #sets viewport
    g.cutonlbj(bbJmax)
    #print "within plotradcuts, after cutonlbj, len(g.x1) = ",len(g.x1)
    nr=1.
    c.measurengalcontam(nv,nr,g)
    #print "nr = ",nr, " ave contam = ",N.average(c.contam)
    sub1plotngalmcl(c.mass,c.membincut,c.obsmembincut)
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    #label="R\dp\u < "+str(nr)+"R\dv\u"
    label=str1+str(nr)+str2
    ppgplot.pgtext(xlabel,ylabel,label) 


    nr=.5
    ppgplot.pgsvp(x1,x2,y1,y2)  #sets viewport
    #ppgplot.pgpanl(1,1)
    c.measurengalcontam(nv,nr,g)
    #print "nr = ",nr, " ave contam = ",N.average(c.contam)
    sub1plotngalmcl(c.mass,c.membincut,c.obsmembincut)
    label=str1+str(nr)+str2
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    ppgplot.pgtext(xlabel,ylabel,label) 


    ppgplot.pgend()

def sub1plotngalmcl(x,y1,y2):
    schdef=1.2
    slwdef=4
    ppgplot.pgsch(schdef)
    xmin=-.5
    xmax=1.
    ymin=0. 
    ymax=60.
    #nbin=5
    ppgplot.pgslw(slwdef)   #line width
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('bcnlst',1.0,0,'bcvnst',10.,2)  #tickmarks and labeling
    ppgplot.pgmtxt('b',2.5,0.5,0.5,"log\d10\u(10\u14\d M\dhalo\u/h\u-1\d M\d\(2281)\u)")    #xlabel
    ppgplot.pgmtxt('l',2.1,0.5,0.5,"N\dgal\u")

    (xbin,ybin,ybinerr)=my.biniterr(x,y1,nbin)
    print y1
    print 'ybin for y1 = ',ybin

    xbin=N.log10(xbin)+12.-14.#add back 10^12 Msun, then divide by 10^14
    ppgplot.pgsch(1.5)
    ppgplot.pgpt(xbin,ybin,7)
    ppgplot.pgslw(3)
    ppgplot.pgerrb(6,xbin,ybin,ybinerr,2.)
    #my.errory(xbin,ybin,ybinerr)
    (xbin,ybin,ybinerr)=my.biniterr(x,y2,nbin)
    print y2
    print 'ybin = ',ybin

    xbin=N.log10(xbin)+12.-14.#add back 10^12 Msun, then divide by 10^14
    ppgplot.pgsch(1.75)
    ppgplot.pgpt(xbin,ybin,17)
    ppgplot.pgsch(schdef)
    ppgplot.pgerrb(6,xbin,ybin,ybinerr,2.)
    ppgplot.pgslw(slwdef)


def linelabel(xlabel,ylabel,dx,ystep,dxl,dyl,label):#draw key
    schdef=ppgplot.pgqch()
    ppgplot.pgsch(1.1)
    ppgplot.pgtext(xlabel,ylabel,label)
    ppgplot.pgsch(1.1)
    ylabel=ylabel-2.*ystep
    xs=N.arange(xlabel,(xlabel+dx),.01)
    ys=ylabel*N.ones(len(xs),'d')
    ppgplot.pgsls(1)
    ppgplot.pgline(xs,ys)
    ppgplot.pgtext((xlabel+dxl+dx),(ylabel-dyl),"Memb")
    ylabel=ylabel-ystep
    ys=ylabel*N.ones(len(xs),'d')
    ppgplot.pgsls(4)
    ppgplot.pgline(xs,ys)
    ppgplot.pgtext((xlabel+dxl+dx),(ylabel-dyl),"Contam")
    ppgplot.pgsch(schdef)
    ppgplot.pgsls(1)

class Cluster:
    def __init__(self):
        self.id = []
        self.r200  = []
        self.sigma = []
        self.z  = []
        self.x1  = []
        self.x2  = []
        self.x3  = []

    def readfiles(self):
	ncopy=2
	for j in range(ncopy):
	    infile=open('/Users/rfinn/SDSS/fieldDR4/myclusters.cat','r')
	    for line in infile:
		if line.find('#') > -1:
		    continue
		t=line.split()
		self.id.append(float(t[0]))#C4 id name
		self.r200.append(float(t[3]))#R200 in Mpc
		self.sigma.append(float(t[4]))#sigma in km/s
		self.z.append(float(t[1]))
		c1=Rand.randint(0,len(g.x1all))#center on random galaxy in gal catalog
		#c2=Rand.randint(0,len(g.x1))#center on random galaxy in gal catalog
		#c3=Rand.randint(0,len(g.x1))#center on random galaxy in gal catalog
		self.x1.append(g.x1all[c1])
		self.x2.append(g.x2all[c1])
		self.x3.append(g.x3all[c1])

		#c1=Rand.random()#center on random position in simulation
		#c2=Rand.random()#center on random position in simulation
		#c3=Rand.random()#center on random 
		#self.x1.append(c1*simL)
		#self.x2.append(c2*simL)
		#self.x3.append(c3*simL)
	    infile.close()


	self.r200=N.array(self.r200,'d')
	self.sigma=N.array(self.sigma,'d')
	self.z=N.array(self.z,'d')
	self.x1=N.array(self.x1,'d')
	self.x2=N.array(self.x2,'d')
	self.x3=N.array(self.x3,'d')

	

    def getvobs(self):
        vh=((1.+self.z)**2 - 1)/((1+self.z)**2 +1)*3.e5
        #self.vobs=(vh + self.v3)/(1+vh*self.v3/(3.e5)**2)#Hubble flow
        #self.z=N.sqrt((1+self.vobs/3.e5)/(1-self.vobs/3.e5)) - 1.#observed redshift
        #print "vobs, z",self.vobs[0:5],self.z[0:5]
        #self.z=self.vobs/3.e5#observed redshift


    def getmemberids(self):#get galaxy members
        self.membids=N.zeros(len(self.x1),'d')
        self.membflag=N.ones(len(self.x1),'f')#flag to indicate if no members are found
        self.membincut=N.zeros(len(self.x1),'d')
        self.membids=list(self.membids)
        self.membdr=N.zeros(len(self.x1),'d')
        self.membdr=list(self.membdr)
	nomemb=0.
        for i in range(len(self.x1)):
            delta=nr#tolerance for matching
            temp=[]
	    mdv=[]
	    mdr=[]
	    distance=N.sqrt((self.x1[i]-g.x1)**2+(self.x2[i]-g.x2)**2+(self.x3[i]-g.x3)**2)/self.r200[i]
            (dsort,dindex)=my.sortwindex(distance)
            (temp,matchflag)=my.findmatch(0.,dsort,delta)
            #templist1=temp
            #for j in range(len(temp)):
            #    templist1[j]=dindex[temp[j]]

	    #membid=copy.copy(temp)
	    membid=[]
	    dr=[]
	    #print 'cdtheta = ',cdtheta
	    #print "temp=",temp
            try:
		for j in range(len(temp)):
		    membid.append(dindex[temp[j]])
		    #print j,temp[j],dindex[temp[j]],'temp=',temp
		    dr.append(distance[membid[j]])
		self.membids[i]=membid
		self.membincut[i]=len(membid)
		self.membdr[i]=dr
            except TypeError:
                #if (self.membids[i]+999) < 1.:
		nomemb=nomemb+1.
                #print "no members in cluster ",i,self.membids[i]
                self.membincut[i]=0
		self.membflag[i]=0.
                continue
	print nomemb," clusters with no members",nomemb/len(self.x1)
	    #print "cluster ",i," self.membincut[i] = ",self.membincut[i],self.membids[i],'temp=',temp
    def getobsmemberids(self):#get galaxies within nr*dr and nv*dv
        #print "within getobsmemberids(), nv,nr,len(x1) = ",nv,nr,len(g.x1)
        #print "within measurecontam(), nv,nr,len(x1) = ",nv,nr,len(g.x1),len(self.mass)
        #print self.mass
        self.obsmembids=N.zeros(len(self.x1),'d')
        self.obsmembflag=N.ones(len(self.x1),'d')
        self.obsmembids=list(self.obsmembids)
        self.obsdr=N.zeros(len(self.x1),'d')
        self.obsdr=list(self.obsdr)
        self.obsdv=N.zeros(len(self.x1),'d')
        self.obsdv=list(self.obsdv)
        self.obsmembincut=N.zeros(len(self.x1),'d')
	nomemb=0.
        vbox=my.vofz(zbox)#redshift of box center given epoch of observation
        #i.e. recession velocity corresponding to z=0.2, 0.4, 0.6, etc of GIF observations
        for i in range(len(self.x1)):
            #self.obsmembids[i]=[]
            vcl=H0*self.x3[i]
            deltar=nr*self.r200[i]#tolerance for matching
            deltav=nv*self.sigma[i]#tolerance for matching
            #print i,"cluster id, mass, r200, sigma = ",self.id[i],self.mass[i],self.r200[i],self.sigma[i],Mclmin
            temp=[]
            #put cluster in middle of box
            x1=g.x1-self.x1[i]+0.5*simL
            #flip coordinates so galaxy x1 ranges from 0 to simL 
            for j in range(len(x1)):
                if x1[j] > simL:
                    x1[j]=x1[j]-simL
                if x1[j] < 0:
                    x1[j]=x1[j]+simL
            x2=g.x2-self.x2[i]+0.5*simL
            for j in range(len(x2)):
                if x2[j] > simL:
                    x2[j]=x2[j]-simL
                if x2[j] < 0:
                    x2[j]=x2[j]+simL
            x3=g.x3-self.x3[i]+0.5*simL
            for j in range(len(x3)):
                if x3[j] > simL:
                    x3[j]=x3[j]-simL
                if x3[j] < 0:
                    x3[j]=x3[j]+simL
            gvobs=x3*H0 + g.v3
            gvobs=my.relsumv(vbox,gvobs)
            gz=my.zofv(gvobs)
            (x1sort,x1index)=my.sortwindex(x1)
            (x2sort,x2index)=my.sortwindex(x2)
            (temp,matchflag1)=my.findmatch(0.5*simL,x1sort,deltar)
	    if (matchflag1 > 0):
		templist1=[]
		for j in range(len(temp)):
		    templist1.append(x1index[temp[j]])
	    temp=[]            
            (temp,matchflag2)=my.findmatch(0.5*simL,x2sort,deltar)
	    if (matchflag2 > 0):
		templist2=[]
		for j in range(len(temp)):
		    templist2.append(x2index[temp[j]])
	    if (matchflag1 + matchflag2) < 2.:
		nomemb=nomemb+1.
		self.obsmembflag[i]=0.
		print "no observed members in cluster ",i
		continue
	    a=sets.Set(templist1)
	    b=sets.Set(templist2)
	    members=a&b
	    members=list(members)
		
            #templist1=N.argsort(templist1)
            #templist2=N.argsort(templist2)
            subset=[]
	    odv=[]
	    odr=[]

            ccvobs=my.relsumv(vbox,0.5*simL*H0)
            ccz=my.zofv(ccvobs)
            cdtheta=self.r200[i]/my.DA(ccz,h)
            ccx1=0.5*simL
            ccx2=0.5*simL
            for k in members:
                dr=(N.sqrt((x1[k]-ccx1)**2+(x2[k]-ccx2)**2))/self.r200[i]
                if dr < nr:
                    dv = (gvobs[k]-ccvobs)/(1-gvobs[k]*ccvobs/(9.e10))/self.sigma[i]
                    if abs(dv) <= nv:#deltav:
                        gdtheta=dr*self.r200[i]/my.DA(gz[k],h)
                        if gdtheta <= cdtheta:#make sure galaxy is in cone field-of-view
                            subset.append(k)
			    odv.append(dv)
			    odr.append(dr)
            self.obsmembids[i]=subset
	    self.obsdr[i]=odr
	    self.obsdv[i]=odv
            self.obsmembincut[i]=len(subset)
	    if self.obsmembincut[i] < 1.:
		nomemb=nomemb+1.
		self.obsmembflag[i]=0.
	print nomemb," clusters with no observed members",nomemb/len(self.x1)

    def sumitmemb(self,sig,r):
	ntot=N.zeros(len(self.x1),'d')
	sfr=N.zeros(len(self.x1),'d')
	stellmass=N.zeros(len(self.x1),'d')
	yave=N.zeros(len(self.x1),'d')
	for i in range(len(sfr)):
	    if self.membflag[i] < 1.:
		continue
	    membids=list(self.membids[i])
	    dr=list(self.membdr[i])
	    #print "HEY",i,r,dr
	    for j in range(len(membids)):
		k=int(membids[j])
		if (dr[j] < r):
		    ntot[i]=ntot[i]+1.
		    sfr[i]=sfr[i]+g.sfr[k]
		    stellmass[i]=stellmass[i]+g.StellarMass[k]
		    #print "check",i,j,k,membids[j],sfr[i],g.sfr[k],stellmass[i],g.StellarMass[k],self.membids[i]
	return ntot,sfr,stellmass

    def sumitobsmemb(self,sig,r):
	ntot=N.zeros(len(self.x1),'d')
	sfr=N.zeros(len(self.x1),'d')
	stellmass=N.zeros(len(self.x1),'d')
	yave=N.zeros(len(self.x1),'d')
	for i in range(len(sfr)):
	    if self.membflag[i] < 1.:
		continue
	    membids=list(self.obsmembids[i])
	    dr=list(self.obsdr[i])
	    for j in range(len(membids)):
		k=int(membids[j])
		if (dr[j] < r):
		    ntot[i]=ntot[i]+1.
		    sfr[i]=sfr[i]+g.sfr[k]
		    stellmass[i]=stellmass[i]+g.StellarMass[k]
	return ntot,sfr,stellmass


    def calctotalsfr(self):
	self.osfr05=N.zeros(len(self.x1),'d')#observed quantities
	self.osfr1=N.zeros(len(self.x1),'d')
	self.osfr2=N.zeros(len(self.x1),'d')
	self.ongal05=N.zeros(len(self.x1),'d')
	self.ongal1=N.zeros(len(self.x1),'d')
	self.ongal2=N.zeros(len(self.x1),'d')
	self.ostellmass05=N.zeros(len(self.x1),'d')#observed quantities
	self.ostellmass1=N.zeros(len(self.x1),'d')
	self.ostellmass2=N.zeros(len(self.x1),'d')
	
	(self.ongal05,self.osfr05,self.ostellmass05)=self.sumitobsmemb(nv,0.5)
	(self.ongal1,self.osfr1,self.ostellmass1)=self.sumitobsmemb(nv,1.)
	(self.ongal2,self.osfr2,self.ostellmass2)=self.sumitobsmemb(nv,2.)
	self.sfr05=N.zeros(len(self.x1),'d')#actual quantities
	self.sfr1=N.zeros(len(self.x1),'d')
	self.sfr2=N.zeros(len(self.x1),'d')
	self.ngal05=N.zeros(len(self.x1),'d')
	self.ngal1=N.zeros(len(self.x1),'d')
	self.ngal2=N.zeros(len(self.x1),'d')
	self.stellmass05=N.zeros(len(self.x1),'d')#actual quantities
	self.stellmass1=N.zeros(len(self.x1),'d')
	self.stellmass2=N.zeros(len(self.x1),'d')
	
	(self.ngal05,self.sfr05,self.stellmass05)=self.sumitmemb(nv,0.5)
	(self.ngal1,self.sfr1,self.stellmass1)=self.sumitmemb(nv,1.)
	(self.ngal2,self.sfr2,self.stellmass2)=self.sumitmemb(nv,2.)
	output=open('millen.dat','w')
	output.write("#sigma ngal05 ngal1 ngal2 sfr05 sfr1 sfr2 stellmass05 stellmass1 stellmass2 ongal05 ongal1 ongal2 osfr05 osfr1 osfr2 ostellmass05 ostellmass1 ostellmass2 \n")
	for i in range(len(self.ngal05)):
	    s="%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f \n"%(self.sigma[i],self.ngal05[i],self.ngal1[i],self.ngal2[i],self.sfr05[i],self.sfr1[i],self.sfr2[i],self.stellmass05[i],self.stellmass1[i],self.stellmass2[i],self.ongal05[i],self.ongal1[i],self.ongal2[i],self.osfr05[i],self.osfr1[i],self.osfr2[i],self.ostellmass05[i],self.ostellmass1[i],self.ostellmass2[i])
	    output.write(s)
	output.close()

    def readdatafile(self):
	input=open('millen.dat','r')
	i=0
	ngal=0
	for line in input:
            if my.beginsWith('#',line):
                continue
	    ngal = ngal+1
	input.close()
	self.sigma=N.zeros(ngal,'d')
	self.ngal05=N.zeros(ngal,'d')
	self.ngal1=N.zeros(ngal,'d')
	self.ngal2=N.zeros(ngal,'d')
	self.sfr05=N.zeros(ngal,'d')
	self.sfr1=N.zeros(ngal,'d')
	self.sfr2=N.zeros(ngal,'d')
	self.stellmass05=N.zeros(ngal,'d')
	self.stellmass1=N.zeros(ngal,'d')
	self.stellmass2=N.zeros(ngal,'d')
	self.ongal05=N.zeros(ngal,'d')
	self.ongal1=N.zeros(ngal,'d')
	self.ongal2=N.zeros(ngal,'d')
	self.osfr05=N.zeros(ngal,'d')
	self.osfr1=N.zeros(ngal,'d')
	self.osfr2=N.zeros(ngal,'d')
	self.ostellmass05=N.zeros(ngal,'d')
	self.ostellmass1=N.zeros(ngal,'d')
	self.ostellmass2=N.zeros(ngal,'d')


	input=open('millen.dat','r')
	i=0
	for line in input:
            if my.beginsWith('#',line):
                continue
	    t=line.split()
	    for j in range(len(t)):
		t[j]=float(t[j])
	    (self.sigma[i],self.ngal05[i],self.ngal1[i],self.ngal2[i],self.sfr05[i],self.sfr1[i],self.sfr2[i],self.stellmass05[i],self.stellmass1[i],self.stellmass2[i],self.ongal05[i],self.ongal1[i],self.ongal2[i],self.osfr05[i],self.osfr1[i],self.osfr2[i],self.ostellmass05[i],self.ostellmass1[i],self.ostellmass2[i])=t
	    #s="%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f \n"%(self.sigma[i],self.ngal05[i],self.ngal1[i],self.ngal2[i],self.sfr05[i],self.sfr1[i],self.sfr2[i],self.stellmass05[i],self.stellmass1[i],self.stellmass2[i],self.ongal05[i],self.ongal1[i],self.ongal2[i],self.osfr05[i],self.osfr1[i],self.osfr2[i],self.ostellmass05[i],self.ostellmass1[i],self.ostellmass2[i])
	    #print i,s
	    i=i+1
	input.close()


class Galaxy:
    def __init__(self):
        self.x1 = []#in Mpc/h
        self.x2  = []
        self.x3  = []
        self.x1all = []#in Mpc/h
        self.x2all  = []
        self.x3all  = []
        self.v1 = []#in km/s
        self.v2  = []
        self.v3  = []
	self.mu = []#-5logh
	self.mg = []
	self.mr = []
	self.mi = []
	self.mz = []
	self.bulgemu=[]#u mag of bulge, in 10^10 Msun/h
	self.bulgemg=[]
	self.bulgemr=[]
	self.bulgemi=[]
	self.bulgemz=[]
	self.StellarMass=[]#in 10^10 Msun/h
	self.BulgeMass=[]
	self.ColdGas=[]
	self.HotGas=[]
	self.EjectedMass=[]
	self.BlackHoleMass=[]
	self.sfr=[]#in Msun/yr


    def readfiles(self):
	mfile=int(sys.argv[2])
	if (mfile < 1.):

	    infile=open('/Users/rfinn/clusters/millenium/croton_etal.ugriz.mini.ascii/croton_etal.ugriz.mini.ascii','r')
	if (mfile > 0.):
	    infile=open('/Users/rfinn/clusters/millenium/croton_etal.ugriz.ascii/croton_etal.ugriz.ascii','r')

        i=0
        for line in infile:
            if line.find('#') > -1:
                continue
            fields=line.split()
            if float(fields[8]) < 99.:
		self.x1all.append(float(fields[0])/h)
		self.x2all.append(float(fields[1])/h)
		self.x3all.append(float(fields[2])/h)
            if (float(fields[8])+5.*N.log10(h)) < mabscut:
                self.x1.append(float(fields[0])/h)
                self.x2.append(float(fields[1])/h)
                self.x3.append(float(fields[2])/h)
                self.v1.append(float(fields[3]))
                self.v2.append(float(fields[4]))
                self.v3.append(float(fields[5]))
                self.mu.append(float(fields[6])+5.*N.log10(h))
                self.mg.append(float(fields[7])+5.*N.log10(h))
                self.mr.append(float(fields[8])+5.*N.log10(h))
                self.mi.append(float(fields[9])+5.*N.log10(h))
                self.mz.append(float(fields[10])+5.*N.log10(h))
                self.StellarMass.append(float(fields[16])/h)
                self.BulgeMass.append(float(fields[17])/h)
                self.sfr.append(float(fields[22]))

	self.x1=N.array(self.x1,'d')
	self.x2=N.array(self.x2,'d')
	self.x3=N.array(self.x3,'d')
	self.v1=N.array(self.v1,'d')
	self.v2=N.array(self.v2,'d')
	self.v3=N.array(self.v3,'d')
	self.mu=N.array(self.mu,'d')
	self.mg=N.array(self.mg,'d')
	self.mr=N.array(self.mr,'d')
	self.mi=N.array(self.mi,'d')
	self.mz=N.array(self.mz,'d')
	self.StellarMass=N.array(self.StellarMass,'d')
	self.BulgeMass=N.array(self.BulgeMass,'d')
	self.sfr=N.array(self.sfr,'d')
	print "all galaxies, no mag cut, = ",len(self.x1all)
        
    def getvobs(self):
        self.vobs=self.x3*H0+self.v3#Hubble flow
        #should evolve H0
        Hz=H0
        v=self.x3*Hz
        self.vobs=(v+self.v3)/(1+v*self.v3/(3.e5)**2)#Hubble flow
        self.ovobs=copy.copy(self.vobs)
        self.z=N.sqrt((1+self.vobs/3.e5)/(1-self.vobs/3.e5))-1.#observed redshift
zbox=0.07#center of cluster distribution
nv=2.
nr=4.
mode=int(sys.argv[1])
if (mode < 1.):
    g=Galaxy()
    print "Reading galaxy file"
    g.readfiles()
    print "got ",len(g.x1)," galaxies"
    simL=max(g.x1all)
    g.getvobs()

    c=Cluster()
    print "Reading cluster file"
    c.readfiles()
    print "got ",len(c.sigma)," clusters"
    c.getvobs()
    print "Getting true members"
    #pylab.plot(g.x1all,g.x2all,'k.')
    #r=pylab.median(c.r200)
    #pylab.scatter_classic(g.x1,g.x2,s=r,c='r')

    #print "Median r200 = ",r,min(c.r200),max(c.r200)
    #pylab.savefig('galpos.eps')
    c.getmemberids()
    print "elapsed time = ",time.clock()-starttime
    print "Getting observed members"
    c.getobsmemberids()
    print "elapsed time = ",time.clock()-starttime
    print "Calculating ngal,sfr,stellmass"
    c.calctotalsfr()
if (mode >0.):
    c=Cluster()
    c.readdatafile()
print "Making plots"
plotngalsigmapylab()
plotsfrsigmapylab()
plotstellarsigmapylab()
#pylab.plot(c.sigma,c.ongal2,'bo')
#pylab.xlabel('velocity dispersion')
#pylab.savefig("test.eps")
    
endtime=time.clock()
print "end time = ",endtime
print "elapsed time = ",endtime-starttime

