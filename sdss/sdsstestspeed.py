#!/usr/bin/env python
import sys, glob
import Numeric as N
import scipy, time
from math import *
import ppgplot
starttime=time.clock()
print "start time = ",starttime
omega0=0.3
omegaL=0.7
zminp=.04#min z for plots
zmaxp=.11#max z for plots
H0=.7
c=3.e5
superramin=12.*15.
superramax=14.5*15.
superdecmin=-10.
superdecmax=10.
superzmin=0.075
superzmax=0.085

#parameters for high-z clusters
sumsfrhiz=N.array([71.7,56,205,28.6,40.1,47,7.7,10.7,(253.*(70./100.)**2)],'f')
sfrcor=N.array([1.,1,1,1,1,1,2.8,2.8,1.],'f')#aperture correction
volcor=N.array([1.,(1./.94),(1./.86),1,1,1,(1./.99),(1./.74),1.],'f')#aperture correction
sumsfrhiz=sumsfrhiz*sfrcor*volcor
#sfrmhiz=N.log10(N.array([70,20,26,58,5.8,4,.3,.8],'f'))
zhiz=N.array([.704,.748,.794,.845,0.833,.228,.32,.183,0.39],'f')
sigmahiz=N.array([418,504,1018,415,1000,1023,1390,1274,560],'f')
kodamamasshiz=N.array([.8,1.4,11.4,(2.3*.7),(11.3*.7),(13.6*.7),(.7*7.3),(8.1*.7),(5.7*.7)],'f')#lensing masses
masshiz=12*(sigmahiz/1000.)**3.*1./N.sqrt(.7+.3*(1+zhiz)**2)
mratio=masshiz/kodamamasshiz
sfrmhiz=(sumsfrhiz/kodamamasshiz)*(.64)#correct for AGN contamination
#completeness cuts
version=1
if version == 1:
    mabsorig=-20.36#store original
    mabscut=-20.36 #corresponds roughly to r=17.5 at z=0.17
    sigmamin=0. #vel cut that corresponds to M200=2x10^14
    fHamin=35.
    lHamin=2.14#min lHa detected corresponding to fHa=35X10^-17erg/s/cm^2 at z=0.15 
    zmin=0.05 #min redshift of sample
    zmax=0.09 #min redshift of sample
    ngalmin=1. #min number of members
    zbin=N.arange(0.02,0.16,.04,'f')
if version == 2:
    mabsorig=-20.4#store original
    mabscut=-20.4 #corresponds roughly to r=17.5 at z=0.17
    sigmamin=0. #vel cut that corresponds to M200=2x10^14
    fHamin=20.
    lHamin=fHamin*4.*3.1415*(dL(zmax))**2#min lHa detected corresponding to fHa=35X10^-17erg/s/cm^2 at z=0.15 
    zmin=0.05 #min redshift of sample
    zmax=0.09 #min redshift of sample
    ngalmin=1. #min number of members
    zbin=N.arange(0.02,0.1,.04,'f')
sigmamin=600.
sigmamax=1200.
Nmin=35 #minimun # so that richness corresponds to M200=2x10^14
ewmin=4.#2A min equivalent width
nsig=3.
nr=1.
ngalmin=30.
z08min=0.07#min z cut for .08 superstructure
z08max=0.09
#g=Galaxy()
#c=Cluster()
def E(z):
    E=sqrt(omega0*(1+z)**3 + omegaL)
    return E

def dL(zmax): #luminosity distance
    zstep=10000.
    dz=zmax/(zstep-1)
    t=0.
    s=0.
    for i in range(zstep):
        zi=dz*(i-1)
        t=t+dz/E(zi)/(1+zi)
        s=s+dz/E(zi)
    dL=c/H0*(1+zi)*s
    return dL
def psplotinit(output):
    file=output+"/vcps"
    ppgplot.pgbeg(file,1,1)
    ppgplot.pgpap(8.,1.)
    ppgplot.pgsch(1.7) #font size
    ppgplot.pgslw(6)  #line width

#def makeps(name,file):
#    ppgplot.pgbeg(file".ps/vcps",1,1)
#    ppgplot.pgpap(8.,1.25)
#    ppgplot.pgsch(1.7) #font size
#    ppgplot.pgslw(4)  #line width
#    #plotrichnessz()
#    name
    

def binit(x,y,n):#bin arrays x, y into n bins, returning xbin,ybin
    nx=len(x)
    #x=N.array(x,'f')
    #y=N.array(y,'f')
    y=N.take(y,N.argsort(x))
    x=N.take(x,N.argsort(x))
    xbin=N.zeros(n,'f')
    ybin=N.zeros(n,'f')
    #ybinerr=N.zeros(n,'f')
    for i in range(n):
        nmin=i*int(float(nx)/float(n))
        nmax=(i+1)*int(float(nx)/float(n))
        xbin[i]=scipy.stats.stats.median(x[nmin:nmax])
        ybin[i]=scipy.stats.stats.median(y[nmin:nmax])
        #xbin[i]=N.average(x[nmin:nmax])
        #ybin[i]=N.average(y[nmin:nmax])
        #ybinerr[i]=scipy.stats.std(y[nmin:nmax])
    return xbin, ybin#, ybinerr

def drawbinned(x,y,nbin):
    xbin,ybin=binit(x,y,nbin)
    #ppgplot.pgsci(2)
    ppgplot.pgline(xbin,ybin)

def drawhist(x,y): #draw histogram of binned data, x=left side of bins, y=number per bin
    x1=[]
    y1=[]
    n=len(x)
    dx=x[1]-x[0]
    for i in range(n):
        x1.append(x[i])
        x1.append(x[i])
    x1.append(x[(n-1)]+dx)
    x1.append(x[(n-1)]+dx)
    y1.append(0.)
    for i in range(n):
        y1.append(y[i])
        y1.append(y[i])
    y1.append(0.)
    x1=N.array(x1)
    y1=N.array(y1)
    ppgplot.pgline(x1,y1)

def cumulative(input):
    x=N.sort(input)
    n=len(input)
    y=N.arange(0,1,(1./n))
    return x,y
def median(x):
    x=N.sort(x)
    n=float(len(x))
    nmed=int(n/2.)
    return x[nmed]

def sumit(sig,r,agn):
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    if agn < 1:
        for i in range(len(y)):
            y[i]=N.sum(N.compress((c.member[i] > 0) & (g.ew > ewmin)  & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.rvir[i]) & (g.Mabs < mabscut) & (g.agn < 1) & (g.lHa > lHamin), g.sfr))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
    if agn > 0:#include agn in summed sfr
        for i in range(len(y)):
            y[i]=N.sum(N.compress((c.member[i] > 0) & (g.ewr > ewmin)  & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.rvir[i]) & (g.Mabs < mabscut) & (g.lHa > lHamin), g.sfr))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr

    return y

def sumitfixedr(sig,r,agn):
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    if agn < 1:
        for i in range(len(y)):
            y[i]=N.sum(N.compress((c.member[i] > 0) & (g.ew > ewmin)  & (g.dv < sig*c.sigma[i]) & (g.dr < r) & (g.Mabs < mabscut) & (g.agn < 1) & (g.lHa > lHamin), g.sfr))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
    if agn > 0:#include agn in summed sfr
        for i in range(len(y)):
            y[i]=N.sum(N.compress((c.member[i] > 0) & (g.ewr > ewmin)  & (g.dv < sig*c.sigma[i]) & (g.dr < r) & (g.Mabs < mabscut) & (g.lHa > lHamin), g.sfr))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr

    return y

def sumstellmass(sig,r):
    y=N.zeros(len(c.z),'f')
    for i in range(len(y)):
        y[i]=N.sum(N.compress((c.member[i] > 0) & (g.stellarmass > 0) & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.rvir[i]) & (g.Mabs < mabscut), g.stellarmass))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
    return y

def richness(sig,r,agn,mabs):
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    if agn < 1:
        for i in range(len(y)):
            y[i]=len(N.compress((c.member[i] > 0) & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.rvir[i]) & (g.Mabs < mabs) & (g.agn < 1), g.sfr))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
    if agn > 0:#include agn in summed sfr
        for i in range(len(y)):
            y[i]=len(N.compress((c.member[i] > 0) & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.rvir[i]) & (g.Mabs < mabs), g.sfr))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr

    return y

def plotmclz():
    x=c.z
    y=N.log10(c.mass)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,.2,3,0,20)
    ppgplot.pglab("z","M\dcl \u (10\u14 \d M\d\(2281)\u)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)

def plotsigmaz():
    x=c.z
    y=c.sigma
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,0,3000,0,0)
    ppgplot.pglab("z","\gs (km/s)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)
def plotrichnesssigmaz():
    #mabs=-21.36
    mabs=mabscut
    x=c.z
    y=richness(6,1,0,mabs)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,0.,100,0)
    ppgplot.pglab("z","Richness(\gs\dmax\u)","")
    #ppgplot.pgpt(x,y,3)
    drawbinned(x,y,5)
    ppgplot.pgsci(2)
    y=richness(3,1,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(3)
    y=richness(2,1,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)
    y=richness(1,1,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)
    ppgplot.pgsls(3)#dashed line for 3 rvir cuts
    y=richness(6,3,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(2)
    y=richness(3,3,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(3)
    y=richness(2,3,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)
    y=richness(1,3,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsls(1)#solid line
    ppgplot.pgsci(1)
def plotrichnessr():
    mabs=-21.36
    x=c.z
    sig=3.
    y=richness(sig,3,0,mabs)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,0.,100,0)
    ppgplot.pglab("z","Richness(r\dmax\u)","")
    #ppgplot.pgpt(x,y,3)
    drawbinned(x,y,5)
    ppgplot.pgsci(2)
    y=richness(sig,2,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(3)
    y=richness(sig,1,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)
    y=richness(sig,.5,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)
def plotmagcutz():
    mabs=-21.5
    x=c.z
    y=richness(3,1,0,mabs)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,0.,80,0)
    ppgplot.pglab("z","Richness(M\dR\u)","")
    #ppgplot.pgpt(x,y,3)
    drawbinned(x,y,5)
    ppgplot.pgsci(2)
    mabs=mabs+0.5
    y=richness(3,1,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(3)
    mabs=mabs+0.5
    y=richness(3,1,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)
    mabs=mabs+0.5
    y=richness(3,1,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(5)
    mabs=mabs+0.5
    y=richness(3,1,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(6)
    mabs=mabs+0.5
    y=richness(3,1,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(7)
    mabs=mabs+0.5
    y=richness(3,1,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(8)
    mabs=mabs+0.5
    y=richness(3,1,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)
def plotcummag():
    mmin=-21.36
    magdist=N.compress((g.clusterz < .1) & (g.Mabs < mmin),g.Mabs)
    (x,y)=cumulative(magdist)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmax,xmin,0,1.02,0)
    ppgplot.pglab("M\dR\u","Cumulative Fraction","")
    ppgplot.pgline(x,y)
    ppgplot.pgsci(2)
    magdist=N.compress((g.clusterz > .1) & (g.Mabs < mmin),g.Mabs)
    (x,y)=cumulative(magdist)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(1)
def plotkauffmannz():#fraction with kauffmann mass versus z
    x=c.z
    y=(c.kauffmann)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,0.85,1.,0)
    ppgplot.pglab("z","Fraction w/ M\d* \u","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(1)

def plotcompapcors():
    x=N.compress((g.ewr > ewmin) & (g.clusterz > zmin),g.myapcor)
    y=N.compress((g.ewr > ewmin) & (g.clusterz > zmin),g.apcor)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,0)
    ppgplot.pglab("L\dB\u/L\dmine\u","Ap cor","")
    ppgplot.pgpt(x,y,1)
    ppgplot.pgsci(1)

def plotsfrstellarmass():
    x=(c.sumstellarmass)/(1.e12)
    y=c.sumsfr2
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,7,0,100,0)
    ppgplot.pglab("\gSM\d* \u (10\u12 \d M\d\(2281)\u)","\gSSFR (M\d\(2281)\u yr\u-1\d)","")
    ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(2)
    #y=N.log10(c.sumsfrbalogh)
    #x=N.log10(c.sumstellarmassbalogh)-12
    #ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(4)
    #y=N.log10(c.sumsfrhiz)
    #x=N.log10(c.sumstellarmasshiz)-12
    #ppgplot.pgpt(x,y,20)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)

def plotstellarmassz():
    #y=N.log10(N.compress(c.sumstellarmass > 0,c.sumstellarmass))-12
    y=c.sumstellarmass/(1.e12)
    x=(c.z)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,0,4.5,0)
    ppgplot.pglab("z","\gS M\d*\u (10\u12 \d M\d\(2281)\u)","")
    ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)

def plotgalsfrz():#plot sfr for indiv galaxies versus z
    y=(g.sfr +.1)
    x=(g.z)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,-1,100,0)
    ppgplot.pglab("z","SFR (M\d\(2281)\u yr\u-1\d)","")
    ppgplot.pgpt(x,y,1)
    ppgplot.pgsci(2)
    drawbinned(x,y,10)
    ppgplot.pgsci(1)

def plotgalfHaz():#plot sfr for indiv galaxies versus z
    y=N.compress(g.ewr > ewmin,g.fHa)
    x=N.compress(g.ewr > ewmin,g.z)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,-1,50,0)
    ppgplot.pglab("z","Flux Halpha","")
    ppgplot.pgpt(x,y,1)
    ppgplot.pgsci(2)
    drawbinned(x,y,10)
    ppgplot.pgsci(1)

def plotsfrz():
    y=(c.sumsfr2)
    x=(c.z)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,0,200,0)
    ppgplot.pglab("z","\gSSFR (M\d\(2281)\u yr\u-1\d)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    #y=N.log10(c.sumsfrbalogh)
    #ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(4)
    #y=N.log10(c.sumsfrhiz)
    #ppgplot.pgpt(x,y,20)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)

def plotsfrzallsub(y,n,ptype):
    x,y=binit(c.z,y,n)    
    #y=y/y[0]
    #y=N.log10(y)
    ppgplot.pgline(x,y)
    ppgplot.pgpt(x,y,ptype)    
def plotsfrzallr(): #plot all cuts on SFR versus z
    #x=(c.z)
    nbin=5
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,0,120,0)
    ppgplot.pglab("z","\gSSFR(r\dmax\u) (M\d\(2281)\u yr\u-1\d)","")
    y=sumit(3.,0.5,0)
    plotsfrzallsub(y,nbin,3)
    ppgplot.pgsci(2)
    y=(sumit(3.,1.,0))
    plotsfrzallsub(y,nbin,20)
    ppgplot.pgsci(4)
    y=(sumit(3.,2.,0))
    plotsfrzallsub(y,nbin,11)
    ppgplot.pgsci(3)
    y=(sumit(3.,3.,0))
    plotsfrzallsub(y,nbin,5)
    ppgplot.pgsci(1)
def plotsfrzalls(): #plot all cuts on SFR versus z
    #x=(c.z)
    nbin=5
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,0,80,0)
    ppgplot.pglab("z","\gSSFR(\gs\dmax\u) (M\d\(2281)\u yr\u-1\d)","")
    y=sumit(1,1,0)
    plotsfrzallsub(y,nbin,3)
    ppgplot.pgsci(2)
    y=(sumit(2.,1.,0))
    plotsfrzallsub(y,nbin,20)
    ppgplot.pgsci(4)
    y=(sumit(3.,1.,0))
    plotsfrzallsub(y,nbin,11)
    ppgplot.pgsci(3)
    y=(sumit(6.,1.,0))
    plotsfrzallsub(y,nbin,5)
    ppgplot.pgsci(1)

def plotsfrmclzall(): #plot all cuts on SFR versus z
    #x=(c.z)
    nbin=5
    y=sumit(0.5,nr,0)
    y=y/c.mass
    x,y=binit(c.z,y,nbin)
    #y=N.log10(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0.,.15,0,20,0)
    ppgplot.pglab("z","\gSSFR/M\dcl \u (M\d\(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(2)
    y=(sumit(1.,nr,0))
    y=y/c.mass
    x,y=binit(c.z,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,20)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(4)
    y=(sumit(2.,nr,0))
    y=y/c.mass
    x,y=binit(c.z,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,11)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(3)
    y=(sumit(3.,nr,0))
    y=y/c.mass
    x,y=binit(c.z,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,5)
    ppgplot.pgline(x,y)

    ppgplot.pgsci(1)

def plotsfrz08():
    y=(N.compress((c.z > z08min) & (c.z < z08max),c.sumsfr2))
    x=N.compress((c.z > z08min) & (c.z < z08max),c.z)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,0,2.6,0)
    ppgplot.pglab("z","\gSSFR (M\d\(2281)\u yr\u-1\d)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    #y=N.log10(c.sumsfrbalogh)
    #ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(4)
    #y=N.log10(c.sumsfrhiz)
    #ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(1)



def plotsfrmcl():
    y=c.sumsfr2
    x=N.log10(c.mass)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,1.6,-5,100,10)
    ppgplot.pglab("M\dcl \u (10\u14 \d M\d\(2281)\u)","\gSSFR (M\d\(2281)\u yr\u-1\d)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    #y=N.log10(c.sumsfrbalogh)
    #ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(4)
    #y=N.log10(c.sumsfrhiz)
    #ppgplot.pgpt(x,y,20)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)

def plotavesfrmcl():
    y=average(g.sfr)
    x=N.log10(c.mass)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,1.6,0,2.6,0,10)
    ppgplot.pglab("M\dcl \u (10\u14 \d M\d\(2281)\u)","\gSSFR (M\d\(2281)\u yr\u-1\d)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    #y=N.log10(c.sumsfrbalogh)
    #ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(4)
    #y=N.log10(c.sumsfrhiz)
    #ppgplot.pgpt(x,y,20)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)

def plotrichnessmcl():
    y=(c.richness)
    x=N.log10(c.mass)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,2,0,150,0,10)
    ppgplot.pglab("M\dcl \u (10\u14 \d M\d\(2281)\u)","Richness","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    #y=N.log10(c.sumsfrbalogh)
    #ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(4)
    #y=N.log10(c.sumsfrhiz)
    #ppgplot.pgpt(x,y,20)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)

def plotrichnessz():
    y=(c.richness)
    x=(c.z)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,0,150,0)
    ppgplot.pglab("z","Richness","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)
    ppgplot.pgsci(1)

def plotrichnesssigma():
    y=(c.richness)
    x=(c.sigma)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0)
    ppgplot.pglab("\gs (km/s)","Richness","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)
    ppgplot.pgsci(1)
def plottotallumprichnessz():
    y=(c.totlum/c.richness)
    x=(c.z)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,ymin,ymax,0)
    ppgplot.pglab("z","Total Lum/Richness","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    ppgplot.pgsci(4)
    ppgplot.pgsci(1)

def plottotallumz():
    y=(c.totlum)
    x=(c.z)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,ymin,ymax,0)
    ppgplot.pglab("z","R-band lum","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    ppgplot.pgsci(4)
    ppgplot.pgsci(1)

def plottotallumrichness():
    y=(c.totlum)
    x=(c.richness)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0.,175,ymin,ymax,0)
    ppgplot.pglab("Richness","R-band lum","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    ppgplot.pgsci(4)
    ppgplot.pgsci(1)

def plothiststellarmass():
    x=N.compress((g.stellarmass > 1.E6),g.stellarmass)
    x=N.log10(x)
    xmax=max(x)
    xmin=min(x)
    nx=len(x)
    #print "stellar mass max min = ",xmax,xmin
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(xmin,12,0,2100,0)
    ppgplot.pglab("M\d*\u/M\d\(2281)\u","N\dgal \u","")
    ppgplot.pghist(nx,x,xmin,xmax,20,1)
    ppgplot.pgsci(2)
    ppgplot.pgsci(4)
    ppgplot.pgsci(1)

def plothistnom():
    ppgplot.pgbox("",0.0,0,"L",0.0,0)

    ppgplot.pglab("M\dR\u","N\dgal \u w/o M\d*\u \u","")
    x=N.compress((g.Mabs < 0) & (g.stellarmass < -998.),g.Mabs)
    xmax=max(x)
    xmin=min(x)
    nx=len(x)
    ppgplot.pgenv(xmin,xmax,0,250,0)
    ppgplot.pghist(nx,x,xmin,xmax,10,1)
    ppgplot.pgsci(1)

def plothistmabs():
    z1=0.075
    z2=0.1
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    x=N.compress((g.dv < 3.*g.clustersigma) & (g.dr < g.clusterrvir) & (g.Mabs < 0) & (g.clusterz < z1),g.Mabs)
    ncluster=1.*len(N.compress(c.z < z1,c.z))
    #print g.clusterz
    xmax=max(x)
    xmin=min(x)
    ppgplot.pgenv(xmin,xmax,0,40,0)
    ppgplot.pglab("M\dR\u","N\dgal \u per cluster","")

    bins=N.arange(xmin,xmax,0.5)#make array of bins in step of 0.5mag
    nbin=len(bins)
    y=scipy.stats.histogram2(x,bins)
    y=y/ncluster
    ppgplot.pgpt(bins,y,3)
    drawhist(bins,y)
    #ppgplot.pgbin(nbin,bins,y)

    ppgplot.pgsci(2)
    x=N.compress((g.dv < 3.*g.clustersigma) & (g.dr < g.clusterrvir) & (g.Mabs < 0)  & (g.clusterz > z1) & (g.clusterz < z2),g.Mabs)
    ncluster=len(N.compress((c.z > z1) & (c.z < z2),c.z))
    y=scipy.stats.histogram2(x,bins)
    y=y/ncluster
    ppgplot.pgpt(bins,y,4)
    drawhist(bins,y)
    #ppgplot.pgbin(nbin,bins,y)
    ppgplot.pgsci(4)
    x=N.compress((g.dv < 3.*g.clustersigma) & (g.dr < g.clusterrvir) & (g.Mabs < 0)  & (g.clusterz > z2),g.Mabs)
    ncluster=len(N.compress(c.z > z2,c.z))
    y=scipy.stats.histogram2(x,bins)
    y=y/ncluster
    #ppgplot.pgbin(nbin,bins,y)
    ppgplot.pgpt(bins,y,5)
    drawhist(bins,y)
    ppgplot.pgsci(1)
    x=N.array([-21.36,-21.36])
    y=N.array([0,300])
    ppgplot.pgsls(4)#dashed line
    ppgplot.pgline(x,y)
    ppgplot.pgsls(1)#solid line
    

def plotstellarmassmabs():
    x=N.compress((g.Mabs < 0) & (g.stellarmass > 1.),g.Mabs)
    y=N.log10(N.compress((g.Mabs < 0) & (g.stellarmass > 1.),g.stellarmass))
    xmax=max(x)
    xmin=min(x)
    ymin=min(y)
    ymax=max(y)
    nx=len(x)
    print "log stellar mass",ymin,ymax
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,20)
    ppgplot.pglab("M\dR\u","M\d*\u/M\d\(2281)\u","")
    ppgplot.pgpt(x,y,1)
    ppgplot.pgsci(2)
    ppgplot.pgsci(4)
    ppgplot.pgsci(1)

def plotazmabs():
    x=N.compress((g.Mabs < 0) & (g.stellarmass > 1.),g.Mabs)
    y=(N.compress((g.Mabs < 0) & (g.stellarmass > 1.),g.az))
    xmax=max(x)
    xmin=min(x)
    ymin=min(y)
    ymax=max(y)
    nx=len(x)
    print "log stellar mass",ymin,ymax
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,20)
    ppgplot.pglab("M\dR\u","A\dz\u","")
    ppgplot.pgpt(x,y,1)
    ppgplot.pgsci(2)
    ppgplot.pgsci(4)
    ppgplot.pgsci(1)


def plothistabsmag():
    x=N.compress((g.Mabs < 0),g.Mabs)
    bins=N.arange(min(x),max(x),0.5,'f')
    #print bins
    ntot=scipy.stats.histogram2(x,bins)

    x=N.compress((g.stellarmass < -900.) & (g.Mabs < 0),g.Mabs)
    npart=scipy.stats.histogram2(x,bins)

    nfrac=N.zeros(len(npart),'f')
    nfracerr=N.zeros(len(npart),'f')
    for i in range(len(nfrac)):
        nfrac[i]=float(npart[i])/float(ntot[i])
        nfracerr[i]=sqrt(float(npart[i])+float(npart[i]**2)/float(ntot[i]))/float(ntot[i])
    #print npart
    #for i in range(len(bins)):
    #    print bins[i],npart[i],ntot[i],nfrac[i],nfracerr[i]
    xmax=max(x)
    xmin=min(x)
    nx=len(bins)
    bins=N.array(bins,'f')
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(xmin,xmax,0,0.5,0)
    ppgplot.pglab("M\dR\u","Fraction w/o M\d*\u","")
    #print nx, bins, nfrac
    ppgplot.pgpt(bins,nfrac,3)
    y1=nfrac+nfracerr
    y2=nfrac-nfracerr
    n=len(bins)
    ppgplot.pgerrb(6,bins,nfrac,nfracerr,1)
    #ppgplot.pgerrb(4,n,bins,nfrac,nfracerr,1)

def plotcumulativesfrlowz():
    sfr=N.compress((g.clusterz < .075) & (g.Mabs < 0),g.sfr)
    totsfr=N.sum(sfr)
    mabs=N.compress((g.clusterz < .075) & (g.Mabs < 0),g.Mabs)
    bins=N.arange(min(mabs),max(mabs),0.1,'f')
    nbins=len(bins)
    sfrpart=N.zeros(len(bins),'f')
    for i in range(nbins):
        sfrpart[i]=N.sum(N.compress(mabs < bins[i],sfr))/totsfr
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(-24.,-14.,-0.03,1.03,0)
    ppgplot.pglab("M\dR\u","SFR Fraction","")
    bins=N.array(bins,'f')
    sfrpart=N.array(sfrpart,'f')
    ppgplot.pgsci(2)
    ppgplot.pgline(bins,sfrpart)
    #print bins
    #print sfrpart
    #ppgplot.pgpt(bins,sfrpart,3)
    #print nx, bins, nfrac
    ppgplot.pgsci(1)
    x=N.array([-21.36,-21.36])
    y=N.array([0,300])
    ppgplot.pgsls(4)#dashed line
    ppgplot.pgline(x,y)
    ppgplot.pgsls(1)#solid line

def plotcumulativesfrewlowz():
    sfr=N.compress((g.clusterz < .075) & (g.Mabs < 0) & (g.ew > -100),g.sfr)
    totsfr=N.sum(sfr)
    mabs=N.compress((g.clusterz < .075) & (g.Mabs < 0) & (g.ew > -100),g.ew)
    bins=N.arange(min(mabs),max(mabs),5,'f')
    nbins=len(bins)
    sfrpart=N.zeros(len(bins),'f')
    for i in range(nbins):
        sfrpart[i]=N.sum(N.compress(mabs < bins[i],sfr))/totsfr
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(min(mabs),max(mabs),-0.03,1.03,0)
    ppgplot.pglab("EW (\(2078))","Cumulative SFR","")
    bins=N.array(bins,'f')
    sfrpart=N.array(sfrpart,'f')
    ppgplot.pgsci(2)
    ppgplot.pgline(bins,sfrpart)
    #print bins
    #print sfrpart
    #ppgplot.pgpt(bins,sfrpart,3)
    #print nx, bins, nfrac
    ppgplot.pgsci(1)
    x=N.array([4,4])
    y=N.array([-50,500])
    ppgplot.pgsls(4)#dashed line
    ppgplot.pgline(x,y)
    ppgplot.pgsls(1)#solid line

def plot3sigmas():
    ratio=N.array(len(c.sumsfra),'f')
    ratio=c.sumsfra/c.sumsfrb #total sfr at 1.5 sigma relative to 3 sigma
    y=N.log10(N.clip(ratio,.001,1E7))
    med=median(ratio)
    ave=N.average(ratio)
    print "Median 2/3\gs = ",med, ave
    x=N.log10(c.mass)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,2,-2,1.5,0,30)
    ppgplot.pglab("M\dcl \u (10\u14 \d M\d\(2281)\u)","\gSSFR / \gSSFR (within 3\gs)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    ratio=c.sumsfrc/c.sumsfrb
    y=N.log10(N.clip(ratio,1.E-2,1.E7))
    med=median(ratio)
    ave=N.average(ratio)
    print "Median 6/3\gs = ",med, ave
    ppgplot.pgpt(x,y,4)
    ppgplot.pgsci(1)

def plot3r():
    ratio=N.array(len(c.sumsfra),'f')
    ratio=c.sumsfrd/c.sumsfre #total sfr at 1.5 sigma relative to 3 sigma
    med=median(ratio)
    ave=N.average(ratio)
    print "Median 0.5/1 rvir = ",med, ave
    y=N.clip(ratio,1.E-2,1.E7)
    y=N.log10(y)
    x=N.log10(c.mass)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,2,-2,1.5,0,30)
    ppgplot.pglab("M\dcl \u (10\u14 \d M\d\(2281)\u)","\gSSFR / \gSSFR (within R\dv\u)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    ratio2=c.sumsfrf/c.sumsfre
    y=N.clip(ratio2,1.E-2,1E7)
    y=N.log10(y)
    ppgplot.pgpt(x,y,4)
    ppgplot.pgsci(1)
    med=median(ratio2)
    ave=N.average(ratio2)
    rmax=max(ratio2)
    print "Median, mean, max 2/1 rvir = ",med,ave,rmax

def plotdvvsr():
    y=g.dv/g.clustersigma
    x=(g.dr)/g.clusterrvir
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,3.1,-6.5,6.5,0,0)
    ppgplot.pglab("R / R\dv \u","\gD v/\gs","")
    ppgplot.pgpt(x,y,1)
    y=N.compress(g.sfr > 3, y)
    x=N.compress(g.sfr > 3, x)
    ppgplot.pgsci(2)
    x=N.array([0,3],'f')
    y=N.array([3,0],'f')
    ppgplot.pgline(x,y)
    y=N.array([-3,0],'f')
    ppgplot.pgline(x,y)
    #ppgplot.pgpt(x,y,1)
    ppgplot.pgsci(1)
def plotstellarmassmcl():
#    (a,b)=scipy.stats.spearmanr(c.sumstellarmass,c.mass)
#    print "rank correlation b/w tot Stellar mass and Mcl"
#    print a, b

    y=N.log10(N.clip(c.sumstellarmass,.00001,(max(c.sumstellarmass)+5.)))-12.
    x=N.log10(c.mass)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,1.6,-1.3,1.5,0,30)
    ppgplot.pglab("M\dcl \u (10\u14 \d M\d\(2281)\u)","\gSM\d* \u(10\u12 \d M\d\(2281)\u)","")
    ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(2)
    #y=N.log10(c.sumstellarmassbalogh)-12.
    #ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(4)
    #y=N.log10(c.sumstellarmasshiz)-12.
    #ppgplot.pgpt(x,y,20)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)
    

def plotsfrpstellarmassmcl():
    #y=N.log10(c.sumsfr2/c.sumstellarmass)+10.
    y=(c.sumsfr2/c.sumstellarmass)*1.e10
    x=N.log10(c.mass)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,2,0,1,0,10)
    ppgplot.pglab("M\dcl \u (10\u14 \d M\d\(2281)\u)","\gSSFR/\gSM\d* \u(M\d\(2281)\u yr\u-1\d/10\u10 \d M\d\(2281)\u)","")
    ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)


def plotsfrpmclsigma():
    x=c.sigma
    y=(c.sfrmass)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(300,1500,0,1000,0)#match hiz axes
    ppgplot.pglab("\gs (km s\u-1\d)","\gSSFR/M\dcl \u (M\d\(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    #y=N.log10(c.sfrmassbalogh)
    #ppgplot.pgpt(x,y,3)
    #ppgplot.pgsci(4)
    #y=N.log10(c.sfrmasshiz)
    #ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(1)

def plotsfrtotlumsigma():
    x=c.sigma
    y=(c.sfrtotlum)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(300,1500,0,50,0)#match hiz axes
    ppgplot.pglab("\gs (km s\u-1\d)","\gSSFR/L\dR tot\u(M\d\(2281)\u yr\u-1\d/)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)

def plotsfr1Mpctotlumsigma():
    x=c.sigma
    y=(c.sfr1Mpctotlum)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(300,1500,0,50,0)#match hiz axes
    ppgplot.pglab("\gs (km s\u-1\d)","\gSSFR/L\dR tot\u (R < 1Mpc) (M\d\(2281)\u yr\u-1\d/)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)

def plotsfr1Mpctotlumz():
    x=c.z
    y=(c.sfr1Mpctotlum)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,0,50,0)#match hiz axes
    ppgplot.pglab("\gs (km s\u-1\d)","\gSSFR/L\dR tot\u (R < 1Mpc) (M\d\(2281)\u yr\u-1\d/)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)

def plotsfrpmclmcl():
    #y=N.log10(c.sfrmass)
    y=(c.sfrmass)
    x=N.log10(c.mass)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,2,-.5,30.,0,10)#match hiz axes
    ppgplot.pglab("M\dcl \u (10\u14 \d M\d\(2281) \u)","\gSSFR/M\dcl \u (M\d\(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","")
    ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(2)
    #y=N.log10(c.sfrmassbalogh)
    #ppgplot.pgpt(x,y,3)
    #ppgplot.pgsci(4)
    #y=N.log10(c.sfrmasshiz)
    #ppgplot.pgpt(x,y,3)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)


def plotsfrpmclmcl08():
    #y=N.log10(c.sfrmass)
    y=(N.compress((c.z > z08min) & (c.z < z08max),c.sfrmass))
    x=N.log10(N.compress((c.z > z08min) & (c.z < z08max),c.mass))
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,2,-1,ymax,0,10)#match hiz axes
    ppgplot.pglab("M\dcl \u (10\u14 \d M\d\(2281) \u)","\gSSFR/M\dcl \u (M\d\(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    #y=N.log10(c.sfrmassbalogh)
    #ppgplot.pgpt(x,y,3)
    #ppgplot.pgsci(4)
    #y=N.log10(c.sfrmasshiz)
    #ppgplot.pgpt(x,y,3)
    drawbinned(x,y,3)
    x=N.log10(kodamamasshiz)
    y=sfrmhiz
    ppgplot.pgsci(4)
    ppgplot.pgpt(x,y,13)
    ppgplot.pgsci(1)

def plotsffracmcl08():
    #y=N.log10(c.sfrmass)
    y=(N.compress((c.z > z08min) & (c.z < z08max),c.sffrac))
    x=N.log10(N.compress((c.z > z08min) & (c.z < z08max),c.mass))
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,2,-.01,1,0,10)#match hiz axes
    ppgplot.pglab("M\dcl \u (10\u14 \d M\d\(2281) \u)","Fraction of Star-Forming Galaxies","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    #y=N.log10(c.sfrmassbalogh)
    #ppgplot.pgpt(x,y,3)
    #ppgplot.pgsci(4)
    #y=N.log10(c.sfrmasshiz)
    #ppgplot.pgpt(x,y,3)
    drawbinned(x,y,3)
    x=N.log10(kodamamasshiz)
    y=sfrmhiz
    ppgplot.pgsci(4)
    ppgplot.pgpt(x,y,13)
    ppgplot.pgsci(1)

def plotsfrpstellarmassz():
    y=c.sumsfr2/c.sumstellarmass*1.e10
    x=(c.z)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,0,1,0)
    ppgplot.pglab("z","\gSSFR/\gSM\d* \u(M\d\(2281)\u yr\u-1\d/10\u10 \d M\d\(2281)\u)","")
    ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(2)
    #y=N.log10(c.sumsfrbalogh/c.sumstellarmassbalogh)+10.
    #ppgplot.pgpt(x,y,20)
    #ppgplot.pgsci(4)
    #y=N.log10(c.sumsfrhiz/c.sumstellarmasshiz)+10.
    #ppgplot.pgpt(x,y,20)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)

def plotsfrpmclz():
    #y=(c.sfrmass)
    y=sumit(nsig,nr,0)
    y=y/c.mass
    x=c.z
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,-1,ymax,0,0)
    ppgplot.pglab("z","\gS SFR/M\d cl \u (M\d \(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","")
    ppgplot.pgpt(x,y,20)
    for i in range(len(x)):
        if y[i] > 30:
            print "sfrmcl > 30, z = ",c.z[i],c.ra[i],c.dec[i],c.super[i],superramin,superramax,superdecmin,superdecmax,superzmin,superzmax
    ppgplot.pgsci(2)
    #y=N.log10(c.sfrmassbalogh)
    #ppgplot.pgpt(x,y,3)
    #ppgplot.pgsci(4)
    #y=N.log10(c.sfrmasshiz)
    #ppgplot.pgpt(x,y,3)
    drawbinned(x,y,5)
    x=N.compress(c.super > 0,x)
    y=N.compress(c.super > 0,y)
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(1)

def plotsffracsfrpmcl():
    y=c.sffrac
    #x=N.log10(c.sfrmass)
    x=sumit(nsig,nr,0)
    x=x/c.mass
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(-2.,2.,-1,xmax,0)
    ppgplot.pglab("\gS SFR/M\d cl \u (M\d \(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","SF FRAC","")
    ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(2)
    #y=N.log10(c.sfrmassbalogh)
    #ppgplot.pgpt(x,y,3)
    #ppgplot.pgsci(4)
    #y=N.log10(c.sfrmasshiz)
    #ppgplot.pgpt(x,y,3)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)
def plotavesfrz(): #butcher oemler plot
    y=(c.avesfr)
    x=c.z
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,-0.2,15,0)
    ppgplot.pglab("z","Median SFR of SF Galaxies","")
    ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)

def plotavesfrmcl(): #butcher oemler plot
    y=(c.avesfr)
    x=N.log10(c.mass)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,2,-0.2,15,0)
    ppgplot.pglab("M\dcl\u","Median SFR of SF Galaxies","")
    ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)

def plotsffracz(): #butcher oemler plot
    y=(c.sffrac)
    x=c.z
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,-0.02,1,0)
    ppgplot.pglab("z","Fraction of SF Galaxies","")
    ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)
def plotsffracsigma(): #butcher oemler plot
    y=(c.sffrac)
    x=c.sigma
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(300,1500,-0.02,1,0)
    ppgplot.pglab("\gs (km s\u-1\d)","Fraction of SF Galaxies","")
    ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)

def plotsffracmcl(): #butcher oemler plot
    y=(c.sffrac)
    x=N.log10(c.mass)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,2,-0.02,1,0,10)
    ppgplot.pglab("M\dcl\u (10\u14\d M\d\(2281)\u)","Fraction of SF Galaxies","")
    ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(2)#red
    drawbinned(x,y,5)
    ppgplot.pgsci(4)#blue
    y=N.compress(c.sub < 1,c.sffrac)
    x=N.log10(N.compress(c.sub < 1,c.mass))
    drawbinned(x,y,5)
    ppgplot.pgsci(1)

def plotsffracsfrmcl(): #butcher oemler plot
    #y=N.compress((c.z < 0.09) & (c.z > .07),c.sffrac)
    #x=N.log10(N.compress((c.z < 0.09) & (c.z > .07),c.sfrmass))
    #print "SF FRAC VS SFR/Mcl for 0.07 < z < 0.09"
    #dospear(x,y)
    y=(c.sffrac)
    x=(c.sfrmass)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(-1,xmax,-0.02,1,0)
    ppgplot.pglab("\gS SFR/M\d cl \u (M\d \(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","Fraction of SF Galaxies","")
    ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)#blue
    y=N.compress(c.sub < 1,c.sffrac)
    x=(N.compress(c.sub < 1,c.sfrmass))
    drawbinned(x,y,5)

    ppgplot.pgsci(1)

def plotsffracrichness(): #butcher oemler plot
    #y=N.compress((c.z < 0.09) & (c.z > .07),c.sffrac)
    #x=N.log10(N.compress((c.z < 0.09) & (c.z > .07),c.sfrmass))
    #print "SF FRAC VS SFR/Mcl for 0.07 < z < 0.09"
    #dospear(x,y)
    y=(c.sffrac)
    x=(c.richness)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,150.,-0.02,1,0)
    ppgplot.pglab("Richness","Fraction of SF Galaxies","")
    ppgplot.pgpt(x,y,3)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)

def plotr200rv():
    x=c.rvir
    y=c.r200
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(.5,4,.5,4,0)
    ppgplot.pglab("R\dv\u (Mpc)","R\d200\u (Mpc)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(1)
    x=N.array([0,10],'f')
    y=N.array([0,10],'f')
    ppgplot.pgline(x,y)

def plothizz():
    #y=N.log10(c.sfrmass)
    y=(c.sfrmass)
    x=c.z
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,1,-1,50,0)
    ppgplot.pglab("z","\gS SFR/M\d cl \u (M\d \(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","")
    #ppgplot.pgpt(x,y,3)
    #ppgplot.pgsci(2)
    #y=N.log10(N.clip(c.sfrmassbalogh,.001,1000000))
    y=(c.sfrmassbalogh)
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(4)
    #y=N.log10(N.clip(c.sfrmasshiz,.001,10000000))
    #y=N.log10(N.clip(sumit(3.,.5,0),.001,10000000)/c.mass)
    y=(sumit(3.,.5,0))/c.mass
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(1)

    ppgplot.pgsci(3)
    #y=N.log10(sfrmhiz)
    y=(sfrmhiz)
    ppgplot.pgpt(zhiz,y,-3)
    ppgplot.pgsci(1)

def plothizsigma():
    #y=N.log10(c.sfrmass)
    y=(c.sfrmass)
    x=c.sigma
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(300,1500,-1,50,0)
    ppgplot.pglab("\gs (km s\u-1 \d)","\gS SFR/M\d cl \u (M\d \(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","")
    #ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    #y=N.log10(N.clip(c.sfrmassbalogh,.001,1000000))
    y=c.sfrmassbalogh
    #y=N.log10(c.sfrmassbalogh)
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(4)
    #y=N.log10(N.clip(c.sfrmasshiz,.001,1000000))
    y=c.sfrmasshiz
    #y=N.log10(c.sfrmasshiz)
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(1)

    ppgplot.pgsci(3)
    y=(sfrmhiz)
    ppgplot.pgpt(sigmahiz,y,-3)
    ppgplot.pgsci(1)

def plothizmcl():
    y=sumit(6.,0.5,0)/c.mass
    ymin=min(y)
    print "min C4 within 0.5R = ",ymin
    #y=(N.clip(y,.01,1.E7))
    x=(c.mass)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,30,-1,50,0)
    ppgplot.pglab("M\dcl\u (10\u14\d M\d\(2281)\u yr\u-1\d)","\gS SFR/M\d cl \u (M\d \(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","")
    ppgplot.pgpt(x,y,1)
    drawbinned(x,y,5)
    ppgplot.pgsci(1)
    ppgplot.pgsci(3)
    x=N.log10(N.compress(zhiz < 0.6,kodamamasshiz))
    y=N.compress(zhiz < 0.6,sfrmhiz)
    #ppgplot.pgpt(x,y,-3)
    x=N.log10(N.compress(zhiz > 0.6,kodamamasshiz))
    y=N.compress(zhiz > 0.6,sfrmhiz)
    ppgplot.pgsci(3)
    #ppgplot.pgpt(x,y,-4)
    y=(sfrmhiz)
    ppgplot.pgpt((kodamamasshiz),y,-3)
    ppgplot.pgsci(1)

def plotsfrpstellarmassfracrvir():
    x=N.compress(g.dv < 2*g.clustersigma,g.fracrvir)
    y=N.compress(g.dv < 2*g.clustersigma,g.sfrperstellarmass)
    a=zip(x,y)
    a=N.array(a,'f')
    nx=len(x)
    ny=len(y)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,1)
    ppgplot.pgenv(0,3.2,ymin,ymax,0)
    ppgplot.pglab("r/r\dvirial \u","SFR/M\d* \u (M\d \(2281)\u yr\u-1\d\u)","")
    ppgplot.pgpt(x,y,1)

def plotsfrmabs():
    #x=N.compress((g.Mabs > -999.) & (g.Mabs < 0.) & (g.sfr > 0.) & (g.agn < 1),g.Mabs)
    #y=N.log10(N.compress((g.Mabs > -999.) & (g.Mabs < 0.) & (g.sfr > 0.) & (g.agn < 1),g.sfr))
    x=N.compress((g.Mabs > -999.) & (g.Mabs < 0.) &  (g.agn < 1),g.Mabs)
    y=(N.compress((g.Mabs > -999.) & (g.Mabs < 0.)  & (g.agn < 1),g.sfr))
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,1)
    ppgplot.pgenv(-24.,-14.,-.5,10,0,0)
    ppgplot.pglab("M\dR\u","SFR (M\d \(2281)\u yr\u-1\d)","")
    ppgplot.pgpt(x,y,1)

def plotewmabs():
    x=N.compress((g.Mabs > -999.) & (g.Mabs < 0.)& (g.sfr > -998.) & (g.ew > -100.) & (g.agn < 1),g.Mabs)
    y=N.compress((g.Mabs > -999.) & (g.Mabs < 0.) & (g.sfr > -998.) & (g.ew > -100.) & (g.agn < 1),g.ew)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,1)
    ppgplot.pgenv(-24.,-14.,-100,420.,0)
    ppgplot.pglab("M\dR\u","EW(H\ga)","")
    ppgplot.pgpt(x,y,1)

def plotewsfr():
    x=N.compress((g.Mabs > -999.) & (g.Mabs < 0.) & (g.sfr > -998.) & (g.ew > -998.) & (g.agn < 1),g.sfr)
    y=(N.compress((g.Mabs > -999.) & (g.Mabs < 0.) & (g.sfr > -998.) & (g.ew > -998.)& (g.agn < 1),g.ew))
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,1)
#    ppgplot.pgenv(xmin,xmax,ymin,ymax,0)
    ppgplot.pgenv(xmin,20,-10,200,0)
    ppgplot.pglab("SFR","EW(H\ga)","")
    ppgplot.pgpt(x,y,1)


def dospear(x,y):
    (a,b)=scipy.stats.spearmanr(x,y)
    print "rank correl = %6.3f %6.5f" % (a,b)
    #(a,b)=scipy.stats.kendalltau(x,y)
    #print "kendall tau = %6.3f %6.5f" % (a,b)

def plotewvssfrperstellarmass():
    ppgplot.pgbox("",0.0,0,"",0.0,1)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0)
    ppgplot.pglab("SFR per Stellar Mass","EW","")
    ppgplot.pgpt(x,y,1)
    #plot sumsfr (calculated while reading catalogs) versus sumsfr2 (calculated in cluster class)
    x=c.sumsfr
    y=c.sumsfr2
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    linex=N.array([xmin,xmax])
    liney=N.array([xmin,xmax])

def plotsumsfrvsumsfr2():
    ppgplot.pgbox("",0.0,0,"",0.0,1)
    ppgplot.pgenv(xmin,xmax,xmin,xmax,0)
    ppgplot.pglab("\gS SFR (M\d \(2281)\u yr\u-1\d\u)","\gS SFR2 (M\d \(2281)\u yr\u-1\d\u)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgline(linex,liney)

def plotewvsfracrvir():
    x=(N.compress((g.dv < 2*g.clustersigma),g.fracrvir))
    y=(N.compress((g.dv < 2*g.clustersigma),g.ew))
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,1)
    ppgplot.pgenv(-0.1,3.2,-20.,100.,0)
    ppgplot.pglab("r/r\dvirial \u","LHalpha/L\dR \u","")
    ppgplot.pgpt(x,y,1)
    #plot L-Ha/LR versus sfr/stellar mass
    x=N.log10(N.compress((g.dv < 6*g.clustersigma) & (g.sfrperstellarmass > 0),g.sfrperstellarmass))
    y=(N.compress((g.dv < 6*g.clustersigma) & (g.sfrperstellarmass > 0),g.ew))
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    print xmin,xmax,ymin,ymax

def dostats():
    print "total number of clusters in final sample = ",len(c.z)
    print "fraction per cluster with Kauffmann stellar mass"
    med=scipy.median(c.kauffmann)
    ave=N.average(c.kauffmann)
    std=scipy.stats.std(c.kauffmann)
    print "ave = %6.3f %6.3f, med = %6.3f" % (ave,std,med)

    #agn contamination

    s1=sumit(3.,2.,0) #total sfr w/out agn
    s2=sumit(3.,2.,1) #total sfr w/agn
    agnfrac=(s2-s1)/(s1)
    contamination=N.average(agnfrac)
    std=scipy.stats.std(agnfrac)
    med=scipy.median(agnfrac)
    n1=N.sum(g.agn)
    n2=len(g.agn)
    frac=float(n1)/float(n2)
    print "AGN contamination by number w/in 3sig, 2R = ",frac
    print "AGN contamination in terms of total SFR w/in 3sig,2R= ",contamination,std,med

    s1=sumit(3.,0.5,0) #total sfr w/out agn
    s2=sumit(3.,0.5,1) #total sfr w/agn
    agnfrac=(s2-s1)/(s1)
    contamination=N.average(agnfrac)
    std=scipy.stats.std(agnfrac)
    med=scipy.median(agnfrac)
    print "AGN contamination in terms of total SFR w/in 3sig,0.5R= ",contamination,std,med
    y=N.take(c.sumsfr2,N.argsort(c.sumsfr2))
    x=N.take(c.id,N.argsort(c.sumsfr2))
    z=N.take(c.z,N.argsort(c.sumsfr2))
#    print "Clusters w/lowest tot SFR"
#    for i in range(9):
#        print x[i],z[i],y[i]
    temp=N.array(len(c.sumsfrhiz),'f')
    temp=c.sumsfrbalogh/c.sumsfrhiz
    med=median(temp)
    print "0.5Rvir, 6/3 sigma, med = ",med
    med=N.average(temp)
    std=scipy.stats.std(temp)
    print "0.5Rvir, 6/3 sigma, ave ",med, "+/-",std

    print "Mcl VS Z"
    x=c.z
    y=c.mass
    dospear(x,y)
    

    print "M* VS Mcl"
    x=c.mass
    y=c.sumstellarmass
    dospear(x,y)
    

    print "TOT SFR VS Mcl"
    x=c.mass
    #y=c.sumsfr2
    y=sumit(nsig,nr,0)
    dospear(x,y)
    

    #########################
    print "SF FRAC VS Z"
    x=c.z
    y=c.sffrac
    dospear(x,y)
    

    #########################
    print "SF FRAC VS Mcl"
    x=c.mass
    y=c.sffrac
    dospear(x,y)

    #########################
    print "MEDIAN SFR VS Z"
    x=c.z
    y=c.avesfr
    dospear(x,y)


    print "TOT SFR/Mcl VS Z"
    x=c.z
    y=sumit(nsig,nr,0)
    #y=c.sfrmass
    dospear(x,y)
    

    print "TOT SFR/Mcl VS Mcl"
    x=c.mass
    #y=c.sfrmass
    y=sumit(nsig,nr,0)
    dospear(x,y)

    print "TOT SFR/Mcl VS Z, 0.7 < z < 0.9"
    x=N.compress((c.z > z08min) & (c.z < z08max),c.z)
    y=N.compress((c.z > z08min) & (c.z < z08max),c.sfrmass)
    print "ncluster = ",len(x)
    dospear(x,y)
    

    print "TOT SFR/Mcl VS Mcl, 0.7 < z < 0.9"
    x=N.compress((c.z > z08min) & (c.z < z08max),c.mass)
    y=N.compress((c.z > z08min) & (c.z < z08max),c.sfrmass)
    dospear(x,y)

    print "SF FRAC VS Mcl, 0.7 < z < 0.9"
    x=N.compress((c.z > z08min) & (c.z < z08max),c.mass)
    y=N.compress((c.z > z08min) & (c.z < z08max),c.sffrac)
    dospear(x,y)


    print "TOT SFR/M* VS Z"
    x=c.z
    y=c.sumsfr2/c.sumstellarmass
    dospear(x,y)
    

    print "TOT SFR/M* VS Mcl"
    x=c.mass
    y=c.sumsfr2/c.sumstellarmass
    dospear(x,y)
    
    print "RICHNESS VS Mcl"
    x=c.mass
    y=c.richness
    dospear(x,y)

    print "RICHNESS VS Z"
    x=c.z
    y=c.richness
    dospear(x,y)

    print "SF FRAC VS RICHNESS"
    x=c.richness
    y=c.sffrac
    dospear(x,y)

    print "SF FRAC VS SFR/Mcl"
    x=c.sffrac
    y=c.sfrmass
#    for i in range(len(x)):
#        print i, x[i], y[i]
    dospear(x,y)


    x=len(N.compress((g.clusterz < 0.06) & (g.dv < 0) & (abs(g.dv/g.clustersigma) > 2),g.dv))
    y=len(N.compress((g.clusterz < 0.06) & (abs(g.dv/g.clustersigma) > 2),g.dv))
    frac=float(x)/float(y)
    print "FRACTION with dv < 0 for z < .06 clusters = ",frac

class Cluster:
    def __init__(self):
        self.id = []
        self.ra = []
        self.dec  = []
        self.z = []
        self.rvir = []
        self.sigma = []
        self.sumsfr = []
        self.mass = []
        self.sumsfr2 = []
    def creadfiles(self,clusters):
        print "number of clusters = ",len(clusters)
        for i in range(len(clusters)):
            prefix=clusters[i].split('.')
            j=0
            for line in open(clusters[i]):
                if line.find('#') > -1:
                    continue
                if j == 0:
                    fields=line.split()
                    sum=0
                    if float(fields[5]) < 0:#skip clusters with masses below sigmamin
                        break
                    self.id.append(int(fields[0]))
                    self.ra.append(float(fields[1]))
                    self.dec.append(float(fields[2]))
                    self.z.append(float(fields[3]))
                    self.rvir.append(float(fields[4]))#Rvirial in Mpc
                    self.sigma.append(float(fields[5]))#velocity dispersion in km/s
                    j=1            
                    continue
                if line.find('*') > -1:#skip any galaxies w/ *** in one or more fields
                    continue
                fields=line.split()
                if float(fields[8]) > (-100.): #get rid of large negative values
                    sum += float(fields[8])
            self.sumsfr.append(float(sum))

    def convarray(self):
        self.id = N.array(self.id,'i')
        self.ra = N.array(self.ra,'f')
        self.dec  = N.array(self.dec,'f')
        self.z = N.array(self.z,'f')
        self.rvir = N.array(self.rvir,'f')
        self.sigma = N.array(self.sigma,'f')
        self.sumsfr = N.array(self.sumsfr,'f')
        self.mass = N.zeros(len(self.sigma),'f')
        self.sfrmass = N.zeros(len(self.sigma),'f')

    def limitz(self): #limit cluster sample to z > zmin
        print len(self.id),len(self.z)
        #for i in range(len(self.id)):
        #    print i,self.id[i],self.z[i]
        self.id = N.compress((self.z > zmin) & (self.z < zmax),self.id)
        self.ra = N.compress((self.z > zmin) & (self.z < zmax),self.ra)
        self.dec = N.compress((self.z > zmin) & (self.z < zmax),self.dec)
        self.rvir = N.compress((self.z > zmin) & (self.z < zmax),self.rvir)
        self.sigma = N.compress((self.z > zmin) & (self.z < zmax),self.sigma)
        self.sumsfr = N.compress((self.z > zmin) & (self.z < zmax),self.sumsfr)
        self.mass = N.compress((self.z > zmin) & (self.z < zmax),self.mass)
        self.sfrmass = N.compress((self.z > zmin) & (self.z < zmax),self.sfrmass)
        self.z = N.compress((self.z > zmin) & (self.z < zmax),self.z)
        print len(self.z), " galaxies after redshift cut"

    def limitrichness(self):
        self.ngal = N.zeros(len(self.id),'f')
        for i in range(len(self.id)):
            #print "length c.member[i], g.sfr = ",len(c.member[i]),len(g.sfr)
            self.ngal[i] = len(N.compress((c.member[i] > 0) & (g.Mabs < mabscut) & (g.memb > 0), g.sfr))
        self.id = N.compress(self.ngal > ngalmin,self.id)
        self.ra = N.compress(self.ngal > ngalmin,self.ra)
        self.dec = N.compress(self.ngal > ngalmin,self.dec)
        self.rvir = N.compress(self.ngal > ngalmin,self.rvir)
        self.sigma = N.compress(self.ngal > ngalmin,self.sigma)
        self.sumsfr = N.compress(self.ngal > ngalmin,self.sumsfr)
        self.mass = N.compress(self.ngal > ngalmin,self.mass)
        self.sfrmass = N.compress(self.ngal > ngalmin,self.sfrmass)
        self.z = N.compress(self.ngal > ngalmin,self.z)
        print len(self.z), " galaxies after richness cut"

    def limitsigma(self):
        self.id = N.compress((self.sigma > sigmamin) & (self.sigma < sigmamax),self.id)
        self.ra = N.compress((self.sigma > sigmamin) & (self.sigma < sigmamax),self.ra)
        self.dec = N.compress((self.sigma > sigmamin) & (self.sigma < sigmamax),self.dec)
        self.rvir = N.compress((self.sigma > sigmamin) & (self.sigma < sigmamax),self.rvir)
        self.sumsfr = N.compress((self.sigma > sigmamin) & (self.sigma < sigmamax),self.sumsfr)
        self.mass = N.compress((self.sigma > sigmamin) & (self.sigma < sigmamax),self.mass)
        self.sfrmass = N.compress((self.sigma > sigmamin) & (self.sigma < sigmamax),self.sfrmass)
        self.z = N.compress((self.sigma > sigmamin) & (self.sigma < sigmamax),self.z)
        self.sigma = N.compress((self.sigma > sigmamin) & (self.sigma < sigmamax),self.sigma)
        print len(self.z), " clusters after sigma cut"

    def limitsubstructure(self):
        self.id = N.compress(self.sub < 1,self.id)
        self.ra = N.compress(self.sub  < 1,self.ra)
        self.dec = N.compress(self.sub < 1,self.dec)
        self.rvir = N.compress(self.sub < 1,self.rvir)
        self.sumsfr = N.compress(self.sub < 1,self.sumsfr)
        self.mass = N.compress(self.sub < 1,self.mass)
        self.sfrmass = N.compress(self.sub < 1,self.sfrmass)
        self.z = N.compress(self.sub < 1,self.z)
        self.sigma = N.compress(self.sub < 1,self.sigma)
        self.sub = N.compress(self.sub < 1,self.sub)
        print len(self.z), " clusters after substructure cut"

    def calcmembers(self):
        self.member=N.zeros([len(self.sigma),len(g.dv)],'f')
        for i in range(len(self.sigma)):
            for j in range(len(g.dv)):
                if abs(g.clusterid[j] - self.id[i]) < 1:
                    self.member[i][j]= 1

    def calcrichness(self):
        self.richness=N.zeros(len(self.sigma),'f')
        self.richness=N.zeros(len(self.sigma),'f')
        self.richness1Mpc=N.zeros(len(self.sigma),'f')
        self.totlum=N.zeros(len(self.sigma),'f')#integrated R-band luminosity
        self.totlum1Mpc=N.zeros(len(self.sigma),'f')#integrated R-band luminosity w/in 1Mpc
        for i in range(len(self.richness)):
            self.richness[i]=len(N.compress((abs(g.clusterid - self.id[i])< 1) &  (g.Mabs < mabscut) & (g.memb > 0), g.sfr))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
            self.richness1Mpc[i]=len(N.compress((abs(g.clusterid - self.id[i])< 1) & (g.Mabs < mabscut) & (g.memb > 0) & (g.dr < 1.), g.sfr))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
            self.totlum[i]=N.sum(N.compress((abs(g.clusterid - self.id[i])< 1) & (g.Mabs < mabscut) & (g.memb > 0), g.rlum))#sum total lum to use a proxy for cluster mass
            self.totlum1Mpc[i]=N.sum(N.compress((abs(g.clusterid - self.id[i])< 1) & (g.Mabs < mabscut) & (g.memb > 0) & (g.dr < 1.), g.rlum))#sum total lum to use a proxy for cluster mass

    def calcsumsfr2(self):
#        self.sumsfr2 = N.zeros(len(self.sigma),'f')
        self.sumsfr2=sumit(nsig,nr,0)
        self.sumsfrbalogh=sumit(6,.5,0)
        self.sumsfrhiz=sumit(3.,0.5,0)
        self.sfr1Mpc=sumitfixedr(3.,1.,0)#dv,r(Mpc),agn?
        #self.sumsfra = sumit(2,1,0)
        #self.sumsfrb = sumit(3,1,0)
        #self.sumsfrc = sumit(6,1,0)
        #self.sumsfrd = sumit(3,0.5,0)
        #self.sumsfre = sumit(3,1,0)
        #self.sumsfrf = sumit(3,nr,0)

        self.sumstellarmass=sumstellmass(nsig,nr)
        self.sumstellarmassbalogh=sumstellmass(6,0.5)
        self.sumstellarmasshiz=sumstellmass(3.,0.5)

        self.sffrac = N.zeros(len(self.sigma),'f')
        self.avesfr = N.zeros(len(self.sigma),'f')
        self.kauffmann= N.zeros(len(self.sigma),'f')#frac of galaxies w/kauffmann stellar mass

        #calculate
        #(1) fraction w/stellar mass measured from kauffmann et al
        #(2) fraction of star-forming galaxies
        #(3) median SFR per cluster

        i=0
        for i in range(len(self.id)):
            a=float(len(N.compress((abs(g.clusterid - self.id[i])< 1) & (g.stellarmass > 0),g.stellarmass)))
            b=float(len(N.compress((abs(g.clusterid - self.id[i])< 1),g.stellarmass)))
            if b == 0:
                print "WARNING - no galaxies w/stellar mass measurements in cluster ",self.id[i]
                self.kauffmann[i]=0
            else:
                self.kauffmann[i]=a/b                

            totalngal=float(len(N.compress((abs(g.clusterid - self.id[i])< 1) & (g.lHa > -999.) & (g.dv < nsig*self.sigma[i]) & (g.dr < nr*self.rvir[i]) & (g.Mabs < mabscut), g.sfr)))
            totalsfgal=float(len(N.compress((abs(g.clusterid - self.id[i])< 1) & (g.ewr > ewmin) & (g.lHa > lHamin) & (g.dv < nsig*self.sigma[i]) & (g.dr < nr*self.rvir[i]) & (g.Mabs < mabscut) & (g.agn < 1), g.sfr)))
            temp=len(N.compress((abs(g.clusterid - self.id[i])< 1) & (g.ewr > ewmin) & (g.lHa > lHamin) & (g.dv < nsig*self.sigma[i]) & (g.dr < nr*self.rvir[i]) & (g.Mabs < mabscut) & (g.agn < 1), g.sfr))
            if (temp > 0):
                self.avesfr[i]=scipy.stats.stats.median(N.compress((abs(g.clusterid - self.id[i])< 1) & (g.ewr > ewmin) & (g.lHa > lHamin) & (g.dv < nsig*self.sigma[i]) & (g.dr < nr*self.rvir[i]) & (g.Mabs < mabscut) & (g.agn < 1), g.sfr))
            else:
                self.avesfr[i] = 0
            if totalngal == 0:
                print "WARNING - no galaxies in cluster ",self.id[i]
                self.sffrac[i]=0
            else:
                self.sffrac[i]=totalsfgal/totalngal

    def calcmass(self):
        self.mass=9.78*(self.sigma/1000.)**3.*1./N.sqrt(.7+.3*(1+self.z)**3)
    def calcr200(self):
        self.r200=N.zeros(len(self.z),'f')
        self.r200=1.41*self.sigma/1000.*1./N.sqrt(.7+.3*(1+self.z)**3)
    def calcsfrmass(self):
        self.sfrmasshiz = N.zeros(len(self.sigma),'f')
        self.sfrmassbalogh = N.zeros(len(self.sigma),'f')
        self.sfrmass=self.sumsfr2/self.mass
        self.sfrtotlum=self.sumsfr2/self.totlum
        self.sfrmasshiz=self.sumsfrhiz/self.mass
        self.sfrmassbalogh=self.sumsfrbalogh/self.mass
        self.sfr1Mpcmass=self.sfr1Mpc/self.mass
        self.sfr1Mpctotlum=self.sfr1Mpc/self.totlum1Mpc
    def assignsub(self):#assign substructure 0 or 1
        sub0=[]
        for line in open("../balogh1/clusters-sub0"):#read in list of clusters with substructure flag = 0
            fields=line.split()
            sub0.append(fields[0])
        #assign substructure flag
        self.sub=N.ones(len(self.id),'f')
        i=0
        for a in self.id:
            for b in sub0:
                if b.find(str(a)) > 0:#does cluster have sub=0
                    self.sub[i]=0.
            #        break
            #else:
            #    self.sub[i]=1.
            #print i, "cluster = ",self.id[i],"substructure = ",self.sub[i]
            i += 1
    #def cdocalcs(self):
    def supercluster(self):
        self.super = N.zeros(len(self.ra),'f')

        for i in range(len(self.super)):
            #print "supercluster ",i,self.ra[i],self.dec[i]
            if (self.ra[i] > superramin) & (self.ra[i] < superramax):
                #print "through ra limits"
                if (self.dec[i] > superdecmin) & (self.dec[i] < superdecmax):
                    #print "through dec limits"
                    if (self.z[i] > superzmin) & (self.z[i] < superzmax):
                        self.super[i] = 1.
                        #print "through z limits!"
class Galaxy:
    def __init__(self):
        self.clusterid = []
        self.clusterrvir = []
        self.clustersigma = []
        self.clusterz = []
        self.ra = []
        self.dec  = []
        self.z = []
        self.Mabs=[]
        self.dr=[]
        self.dv=[]
        self.dLMpc=[]
        self.fHa=[]
        self.lHa=[]
        self.apcor=[]
        self.stellarmass=[]
        self.az=[]
        self.ar=[]
        self.O3Hb=[]
        self.N2Ha=[]
        self.sfr=[]
    def greadfiles(self,clusters):
        for i in range(len(clusters)):
            prefix=clusters[i].split('.')
            j=0
            for line in open(clusters[i]):
                if line.find('#') > -1:
                    continue
                if j == 0:
                    fields=line.split()
                    sum=0
                    if float(fields[5]) < 0:#skip clusters with masses below sigmamin
                        break
                    j=1
                    id=int(fields[0])
                    virial=float(fields[4])
                    sigma=float(fields[5])
                    z=float(fields[3])            
                    continue
                if line.find('*') > -1: #skip any galaxies w/ ****** in one or more fields
                    continue
                fields=line.split()
                self.clusterid.append(id)#keep track of cluster id
                self.clusterrvir.append(virial)#keep track of cluster virial radius
                self.clustersigma.append(sigma)#keep track of cluster sigma
                self.clusterz.append(z)#keep track of cluster sigma
                self.ra.append(float(fields[0]))
                self.dec.append(float(fields[1]))
                self.z.append(float(fields[2]))
                self.Mabs.append(float(fields[3]))
                self.dr.append(float(fields[4]))#projected radial distance in Mpc
                self.dv.append(float(fields[5]))#km/s
                self.dLMpc.append(float(fields[6]))#luminosity distance in Mpc
                self.fHa.append(float(fields[7]))#flux of Halpha in 1E-17 ers/s/cm@
                self.lHa.append(float(fields[8]))#luminosity of Halpha in 1E40 erg/s
                self.apcor.append(float(fields[9]))#aperture correction from flux_R(from image)/flux_R(in fiber)
                self.stellarmass.append(float(fields[10]))#log10(Msun)
                self.az.append(float(fields[11]))#magnitudes of extinction
                self.O3Hb.append(float(fields[12]))#lg(OIII/Hbeta) - AGN diagnostic
                self.N2Ha.append(float(fields[13]))#lg(NII/Halpha) - AGN diagnostic

    def convarray(self):
        self.clusterid = N.array(self.clusterid,'i')
        self.clusterrvir = N.array(self.clusterrvir,'f')
        self.clustersigma = N.array(self.clustersigma,'f')
        self.clusterz = N.array(self.clusterz,'f')
        self.ra = N.array(self.ra,'f')
        self.dec  = N.array(self.dec,'f')
        self.z = N.array(self.z,'f')
        self.Mabs= N.array(self.Mabs,'f')
        self.dr=N.array(self.dr,'f')
        self.dv=N.array(self.dv,'f')
        self.dLMpc=N.array(self.dLMpc,'f')
        self.fHa=N.array(self.fHa,'f')
        self.lHa=N.array(self.lHa,'f')
        self.apcor=N.array(self.apcor,'f')
        self.stellarmass=N.array(self.stellarmass,'f')
        self.az=N.array(self.az,'f')
        self.O3Hb=N.array(self.O3Hb,'f')
        self.N2Ha=N.array(self.N2Ha,'f')
    def limitmemb(self):
        self.clusterid =   N.compress(self.memb > 0,self.clusterid)
        self.clusterrvir = N.compress(self.memb > 0,self.clusterrvir)
        self.clustersigma =N.compress(self.memb > 0,self.clustersigma)
        self.clusterz =    N.compress(self.memb > 0,self.clusterz)
        self.ra =          N.compress(self.memb > 0,self.ra)
        self.dec  =        N.compress(self.memb > 0,self.dec)
        self.z =           N.compress(self.memb > 0,self.z)
        self.Mabs=         N.compress(self.memb > 0,self.Mabs)
        self.dr=           N.compress(self.memb > 0,self.dr)
        self.dv=           N.compress(self.memb > 0,self.dv)
        self.dLMpc=        N.compress(self.memb > 0,self.dLMpc)
        self.fHa=          N.compress(self.memb > 0,self.fHa)
        self.lHa=          N.compress(self.memb > 0,self.lHa)
        self.apcor=        N.compress(self.memb > 0,self.apcor)
        self.stellarmass=  N.compress(self.memb > 0,self.stellarmass)
        self.az=           N.compress(self.memb > 0,self.az)
        self.O3Hb=         N.compress(self.memb > 0,self.O3Hb)
        self.N2Ha=         N.compress(self.memb > 0,self.N2Ha)
    def calcextinction(self):
        self.ar=N.zeros(len(self.az),'f')
        for i in range(len(self.az)):
            if self.az[i] > (-998.):
                self.ar[i]=10**((self.az[i]+0.23)/2.5)
            else:
                self.ar[i]=2.5 #assume 1 mag extinction
    def convlHa2sfr(self):
        self.sfr = N.zeros(len(self.lHa),'f')
        self.sfr = 7.9E-2*self.lHa #lHa in units of 1E40 ergs/s
        self.sfr = self.sfr*self.ar #correct for extinction
        #self.sfr = self.sfr*2.5 #correct for extinction
    def calcmyapcor(self):#calculate lHa/4pidL^2
        self.myapcor = N.zeros(len(self.lHa),'f')
        self.myapcor = self.lHa/(4*3.1415*(self.dLMpc**2)*self.fHa*9.548*10.**(-9))
 
    def calcforplots(self):
        self.sfrperstellarmass=N.zeros(len(self.lHa),'f')
        self.fracrvir=N.zeros(len(self.lHa),'f')
        self.fracrvir=self.dr/self.clusterrvir
        for i in range(len(self.sfr)):
            if  (self.stellarmass[i] > (-900.)) & (self.sfr[i] > 0) :
                    self.stellarmass[i]=10**(self.stellarmass[i])
                    self.sfrperstellarmass[i]=float(self.sfr[i])/float(self.stellarmass[i])
            elif self.stellarmass[i] > -998.:
                self.stellarmass[i]=10**(self.stellarmass[i])
                self.sfrperstellarmass[i]=-999.
            else:
                self.sfrperstellarmass[i]=-999.
    def calcew(self):
        self.lR=N.zeros(len(self.lHa),'f')#rband luminosity
        self.ew=N.zeros(len(self.lHa),'f')
        self.errew=N.zeros(len(self.lHa),'f')
        self.agn1=N.zeros(len(self.lHa),'f')#agn flags from Miller et al 2003
        self.agn2=N.zeros(len(self.lHa),'f')
        self.agn3=N.zeros(len(self.lHa),'f')
        file = open('mycat','r')
        i=0
        for line in file:
            fields=line.split()
            #print "mycat line = ",i,len(self.ew)
            self.ew[i]=float(fields[3])
            self.errew[i]=float(fields[4])
            self.agn1[i]=float(fields[5])
            self.agn2[i]=float(fields[6])
            self.agn3[i]=float(fields[7])
            i += 1
    def calcagn(self): #set agn flag to 0 or 1 according to line ratios
        self.agn=N.zeros(len(self.lHa),'f')
        for i in range(len(self.agn)):
            if (self.agn1[i] < 1) | ((self.agn1[i] > 4) & (self.agn[1] < 9)):
                self.agn[i]=0
            else:
                self.agn[i]=1
    def calcrlum(self):
        self.rlum=N.zeros(len(self.Mabs),'f')
        #self.rlum=10.**(-0.4*(self.Mabs-4.77))#R-band Lum in units of L-solar
        self.rlum=10.**(-0.4*(self.Mabs+24))#R-band Lum in units of L-solar
    def calcewr(self):#rest-frame EW
        self.ewr=N.zeros(len(self.ew),'f')
        self.ewr=self.ew/(1.+self.z)
    def calcmemb(self):#assign membership flag based on position in dv/sigma-dr/rvir plane
        self.memb=N.zeros(len(self.dv),'f')
        for i in range(len(self.memb)):
            dv=self.dv[i]/self.clustersigma[i]
            dr=self.dr[i]/self.clusterrvir[i]
            y1=dr-3.
            y2=3.-dr
            #if ((dv > y1) & (dv < y2)):
            if ((dv < nsig) & (dr < nr)):
                self.memb[i]=1
    def limitmemb(self):
        self.clusterid =   N.compress(self.memb > 0,self.clusterid)
        self.clusterrvir = N.compress(self.memb > 0,self.clusterrvir)
        self.clustersigma =N.compress(self.memb > 0,self.clustersigma)
        self.clusterz =    N.compress(self.memb > 0,self.clusterz)
        self.ra =          N.compress(self.memb > 0,self.ra)
        self.dec  =        N.compress(self.memb > 0,self.dec)
        self.z =           N.compress(self.memb > 0,self.z)
        self.Mabs=         N.compress(self.memb > 0,self.Mabs)
        self.dr=           N.compress(self.memb > 0,self.dr)
        self.dv=           N.compress(self.memb > 0,self.dv)
        self.dLMpc=        N.compress(self.memb > 0,self.dLMpc)
        self.fHa=          N.compress(self.memb > 0,self.fHa)
        self.lHa=          N.compress(self.memb > 0,self.lHa)
        self.apcor=        N.compress(self.memb > 0,self.apcor)
        self.stellarmass=  N.compress(self.memb > 0,self.stellarmass)
        self.az=           N.compress(self.memb > 0,self.az)
        self.O3Hb=         N.compress(self.memb > 0,self.O3Hb)
        self.N2Ha=         N.compress(self.memb > 0,self.N2Ha)
        self.lR   =        N.compress(self.memb > 0,self.lR)
        self.ew=           N.compress(self.memb > 0,self.ew)
        self.errew=        N.compress(self.memb > 0,self.errew)
        self.agn1=         N.compress(self.memb > 0,self.agn1)
        self.agn2=         N.compress(self.memb > 0,self.agn2)
        self.agn3=         N.compress(self.memb > 0,self.agn3)

    def gdocalcs(self):
        self.convarray()
        self.calcew()
        self.calcmemb()
        #self.limitmemb()#limit to member galaxies only

        self.calcagn()
        self.calcrlum()
        self.calcewr()
        self.calcextinction()
        self.convlHa2sfr()
        self.calcmyapcor()
        self.calcforplots()
        self.sfr=N.clip(self.sfr,0,100000)
c=Cluster()
g=Galaxy()
def gotoit():
    #c=Cluster()
    #g=Galaxy()
    clusters = glob.glob("*_DR2*")

    c.creadfiles(clusters)
    g.greadfiles(clusters)
    g.gdocalcs()
    #c.cdocalcs()
    c.convarray()
    print "find cluster members"
    c.calcmembers()

    print "starting limitz"
    c.limitz()
    print "starting limitrichness"
    c.limitrichness()
    print "starting limitsigma"
    c.limitsigma()
    print "starting assignsub"
    c.assignsub() #assign substructure flag
    print "starting supercluster"
    c.supercluster() #assign substructure flag
    #print "starting limitsubtructure"
    #c.limitsubstructure() #assign substructure flag
    print "starting calcmass"
    c.calcmass()
    print "starting calcr200"
    c.calcr200()
    print "starting calcsumsfr2"
    c.calcsumsfr2()
    print "starting calcrichness"
    c.calcrichness()
    print "starting calcsfrmass"
    c.calcsfrmass()
    #print "starting dostats"
    #dostats()

def updateplots():
    print "starting plots"

    allgals=0 #plots showing all galaxies in all clusters
    stellmass=1 #plots w/stellar mass
    sfrm=1 #plots involving sum sfr/Mcl
    hiz=1 #comparison w/hi-z clusters
    if plotall == 0 or plotall == 2:
        psplotinit("richnessvsz.ps")
        plotrichnessz()
        psplotinit("r200rv.ps")
        plotr200rv() #plot R200 versus Rvir

        psplotinit("cumulativesfrlowz.ps")
        plotcumulativesfrlowz()
        
        psplotinit("cumulativesfrewlowz.ps")
        plotcumulativesfrewlowz()
        
        psplotinit("mclvsz.ps")#plot cluster mass vs z 
        plotmclz()
        
        psplotinit("sigmaz.ps")    #plot cluster sigma vs z 
        plotsigmaz()
        
        psplotinit("magcutz.ps")
        plotmagcutz()
        
        psplotinit("richnesssigmaz.ps")
        plotrichnesssigmaz()
        
        psplotinit("richnessr.ps")
        plotrichnessr()
        
        #psplotinit("azvsmabs.ps")
        #plotazmabs()

        #psplotinit("stellarmassvsmabs.ps")
        #plotstellarmassmabs()

        #plot sumsfr versus cluster mass
        psplotinit("sfrvsmcl.ps")
        plotsfrmcl()

        psplotinit("sfrtotlumsigma.ps")
        plotsfrtotlumsigma()

        psplotinit("sfr1Mpctotlumsigma.ps")
        plotsfr1Mpctotlumsigma()

        if (allgals > 0):
            psplotinit("dvr.ps")
            plotdvvsr()
            
            psplotinit("sfrmabs.ps")
            plotsfrmabs()

            psplotinit("ewmabs.ps")
            plotewmabs()

            #plot L-Ha/LR versus r(projected)/rvir
            #psplotinit("ewvsfracrvir.ps")
            #plotewvsfracrvir()
    
            #psplotinit("ewvsfrperstellarmass.ps")
            #plotewvsfrperstellarmass()
    
            #psplotinit("sumsfrvsumsfr2.ps")
            #plotsumsfrvsumsfr2()

        if (stellmass > 0):
            psplotinit("stellarmassz.ps")
            plotstellarmassz()
        
            psplotinit("sfrvsstellarmass.ps")    #plot sumsfr versus stellar mass
            plotsfrstellarmass()
        
            psplotinit("sfrperstellarmassvsmcl.ps")#plot sumsfr/stellarmass versus cluster mass
            plotsfrpstellarmassmcl()
        
            psplotinit("stellarmassvsmcl.ps")#plot sum stellarmass versus cluster mass
            plotstellarmassmcl()
            
            psplotinit("sfrperstellarmassvsz.ps")#plot sumsfr/stellarmass versus z
            plotsfrpstellarmassz()

            psplotinit("sfrstellarmassvfracrvir.ps")#plot sfr per stellar mass versus r(projected)/rvir
            plotsfrpstellarmassfracrvir()
        
        if (sfrm > 0):
            psplotinit("sfrmasssigma.ps")#plot sumsfr/Mcl versus sigma
            plotsfrpmclsigma()
    
            psplotinit("sfrpermclvsmcl.ps")#plot sumsfr/Mcl versus cluster mass
            plotsfrpmclmcl()
    
            psplotinit("sfrzallr.ps")#plot sfr/Mcl vs z for 4 radial cuts
            plotsfrzallr()
            
            psplotinit("sfrzalls.ps")#plot sfr/Mcl vs z for 4 radial cuts
            plotsfrzalls()
        
            psplotinit("sfrpermclvsz.ps")#plot sfr per cluster mass vs z
            plotsfrpmclz()

            psplotinit("sffracsfrmcl.ps")
            plotsffracsfrmcl()

            psplotinit("sfrmclmcl08.ps")#sumsfr/mcl vs mcl for z~.08 cl only
            plotsfrpmclmcl08()
        
            psplotinit("sffracmcl08.ps")
            plotsffracmcl08()
    
        if (hiz > 0):
    
            psplotinit("hizz.ps")#plot sfr/Mcl versus z, including hiz clusters
            plothizz()
    
            psplotinit("hizsigma.ps")#plot sfr/Mcl versus sigma, including hiz clusters
            plothizsigma()

            psplotinit("hizmcl.ps")
            plothizmcl()

            psplotinit("sffracz.ps")
            plotsffracz()
            
            psplotinit("avesfrz.ps")
            plotavesfrz()

            #ppgplot.pgbeg("histmabs.ps/vcps",1,1)
            #ppgplot.pgpap(8.,1.25)
            #ppgplot.pgsch(1.7) #font size
            #ppgplot.pgslw(4)  #line width
            #plothistmabs()


            ##############################
            #plot sumsfr/mcl versus sigma w/hi-z clusters
            psplotinit("sffracsigma.ps")
            plotsffracsigma()
            
            psplotinit("sffracmcl.ps")
            plotsffracmcl()

def updateall():
    if plotall == 1 or plotall == 2:
        ppgplot.pgbeg("all.ps/vcps",3,3)
        ppgplot.pgsch(2.) #font size
        ppgplot.pgslw(4)  #line width
        
        plotmclz()
        
        ppgplot.pgpage
        plotmagcutz()

        ppgplot.pgpage
        plotrichnesssigmaz()
        
        ppgplot.pgpage
        plotrichnessr()

        
        ppgplot.pgpage
        plotrichnessz()
        
        ppgplot.pgpage
        plottotallumz()
        
        ppgplot.pgpage
        plottotallumrichness()
        
        ppgplot.pgpage
        plottotallumprichnessz()
        
        ppgplot.pgpage
        plotrichnessmcl()

    
        #ppgplot.pgpage
        #plotgalsfrz()
    
        #ppgplot.pgpage
        #plotgalfHaz()


        #ppgplot.pgpage
        #plotstellarmassmabs()
        ##############################
        #plot sumsfr versus stellar mass
        ppgplot.pgpage
        plotsfrstellarmass()
        ##############################
        #plot sumsfr versus cluster mass
        ppgplot.pgpage
        plotsfrzallr()
        #plot sumsfr versus cluster mass
        ppgplot.pgpage
        plotsfrmcl()
        ppgplot.pgpage
        plotsfrzalls()
        ppgplot.pgpage
        plotsfrmclzall()
        ##############################
        #plot sum stellarmass versus cluster mass
        ppgplot.pgpage
        plotstellarmassmcl()

        ppgplot.pgpage
        plotsfrz()
        
        ppgplot.pgpage
        plotsfrz08()
        
        ppgplot.pgpage
        plotsfrpmclmcl08()

        ##############################
        #plot sumsfr versus cluster mass, for vel cuts of 1.5, 3, 6 sigma

        #ppgplot.pgpage
        #plot3sigmas()
        ##############################
        #plot sumsfr versus cluster mass, for rad cuts of .5, 1, 2 rvir, +/-3sigma
        #ppgplot.pgpage
        #plot3r()
    
        #ppgplot.pgpage
        #plotdvvsr()
        #ppgplot.pgpage
        #plotsfrmabs()

        #ppgplot.pgpage
        #plotewmabs()

        #ppgplot.pgpage
        #plotewsfr()

        ##############################
        #plot sumsfr/stellarmass versus cluster mass
        ppgplot.pgpage
        plotsfrpstellarmassmcl()
        ##############################
        #plot sumsfr/Mcl versus cluster mass
        ppgplot.pgpage
        plotsfrpmclmcl()
        ##############################
        #plot sumsfr/stellarmass versus z
        ppgplot.pgpage
        plotsfrpstellarmassz()
        ##############################
        #plot sfr per cluster mass vs z
        ppgplot.pgpage
        plotsfrpmclz()
        ##############################
        #plot sumsfr/mcl versus z w/hi-z clusters
        ppgplot.pgpage
        plothizz()
        ##############################
        #plot sumsfr/mcl versus sigma w/hi-z clusters
        ppgplot.pgpage
        plothizsigma()

        ppgplot.pgpage
        plothizmcl()
        
        ppgplot.pgpage
        plotsfrtotlumsigma()
        ppgplot.pgpage
        plotsfr1Mpctotlumsigma()
        ppgplot.pgpage
        plotsfr1Mpctotlumz()
        

        ##############################
        #plot sumsfr/mcl versus sigma w/hi-z clusters
        ppgplot.pgpage
        plotsffracz()

        ##############################
        #plot sumsfr/mcl versus sigma w/hi-z clusters
        ppgplot.pgpage
        plotsffracsigma()

        ppgplot.pgpage
        plotsffracrichness()

        ppgplot.pgpage
        plotavesfrz()

        ppgplot.pgpage
        plotavesfrmcl()
        
        ppgplot.pgpage
        plotsffracmcl()
        
        ##############################
        #plot sumsfr/mcl versus sigma w/hi-z clusters

        ppgplot.pgpage
        plotsffracsfrmcl()

        #ppgplot.pgpage
        #plothistmabs()

        ppgplot.pgpage
        plotcumulativesfrlowz()
        
        ppgplot.pgpage
        plotcumulativesfrewlowz()
        
        ppgplot.pgpage
        plotr200rv()
        
        ppgplot.pgbeg("other-plots.ps/vcps",3,3)
        ppgplot.pgsch(2.) #font size
        ppgplot.pgslw(4)  #line width
        ppgplot.pgpage
        plotkauffmannz()
        
        ppgplot.pgpage
        plotcompapcors()
        ppgplot.pgpage
        plothiststellarmass()
        
        ppgplot.pgpage
        plothistnom()#distribution of Mabs for galaxies w/no M*
    
        #ppgplot.pgpage
        #plothistabsmag()
        ppgplot.pgpage
        plotazmabs()
        
        ppgplot.pgend()

#Main
gotoit()
plotall=2
updateplots()#update *.ps plot files
#updateall()#update all.ps

endtime=time.clock()
print "end time = ",endtime
print "elapsed time = ",endtime-starttime
