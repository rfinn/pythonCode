#!/scisoft/bin/python
"""useage
sdssinter.py mode completenesscorr, where mode
= 0 for updating all
= 1 for updating plots only

and completenesscorr
=0 for no completeness correction for fiber sampling
=1 for yes
"""
import sys, glob
import Numeric as N
#import scipy
import pylab #use matplotlib instead of scipy
from math import *
import ppgplot
import time, os
from mystuff import *
import matplotlib
from matplotlib import rc
import mystuff as my
#matplotlib.use('Agg')
#from pylab import savefig

from matplotlib import rc
rc('font',family='serif', style='normal', variant='normal',weight='bold', stretch='normal', size='large')
#rc('text', usetex=True)
import sqlcl
logfile='log'+str(sys.argv[1])+str(sys.argv[2])
log=open(logfile,'w')

#import matplotlib.matlab as mat
#import ediscssfr
#from sdssplots import *
starttime=time.clock()
print "start time = ",starttime
string= "start time = "+str(starttime)
log.write(string)
mode=int(sys.argv[1])#
print "mode = ",mode
getcompl=0. #calculate completeness: 0=no, 1=yes
completenesscorr=int(sys.argv[2])#apply completeness correction for fiber samplin:, 0=no, 1=yes

print "completeness correction = ",completenesscorr
string = "completeness correction = "+str(completenesscorr)
log.write(string)

totrlumlabel="Total L\dR\u (10\u9\d L\d\(2281)\u)"
sfrsigma3label="\gSSFR/(\gs/1000)\u3\d (M\d\(2281)\u yr\u-1\d/(km s\u-1\d)\u3\d)"
omega0=0.3
omegaL=0.7
zminp=.04#min z for plots
zmaxp=.1#max z for plots
h100=.7
H0=100.*h100
c=3.e5
superramin=12.*15.
superramax=14.5*15.
superdecmin=-10.
superdecmax=10.
superzmin=0.075
superzmax=0.085
richmin=0.
richmax=450.
sfrmclminp=-.5
sfrmclmaxp=50.
ngalmin=0.
cdir=str(os.getcwd())
print "Current directory = ",cdir
if cdir.find('balogh') > -1:
    ngalmin=2.#set min ngal to 30 for clusters, 1 for field
mabsorig=-20.38#store original
#mabscut=-20.4 #corresponds roughly to r=17.5 at z=0.17
mabscut=-20.38 #corresponds roughly to r=17.5 at z=0.17
sigmamin=400. #vel cut that corresponds to M200=2x10^14
fHamin=35.
lHamin=2.14#min lHa detected corresponding to fHa=35X10^-17erg/s/cm^2 at z=0.15 
zmin=0.05 #min redshift of sample
#zmin=0.0 #min redshift of sample
zmax=0.09 #min redshift of sample
zbin=N.arange(0.02,0.1,.04,'f')

fHamin=20.
d=dL(zmax,h100)
print "dL (zmax) = ",d
d=d*3.09*10**24.#convert from Mpc to cm
print "dL (zmax) in cm = ",d
lHamin=fHamin*10**(-17.)*4.*3.1415*(d)**2#min lHa detected corresponding to fHa=35X10^-17erg/s/cm^2 at z=0.15
lHamin=lHamin/(10**40.) #convert to units of 10^40 ergs/s to match catalogs
print "min Halpha luminosity (10^40 erg/s)= ",lHamin
sigmamin=500.
sigmamax=1000.
sigmaminp=sigmamin-100.
sigmamaxp=sigmamax+100.
totrlumminp=50.
totrlummaxp=100.

zminpb=.05
zmaxpb=.1
sigmaminpb=200.
sigmamaxpb=1500.
totrlumminpb=50.
totrlummaxpb=800.
Nmin=35 #minimun # so that richness corresponds to M200=2x10^14
ewmin=4.#2A min equivalent width
nsig=2.
nr=1.
z08min=0.07#min z cut for .08 superstructure
z08max=0.09
#g=Galaxy()
#c=Cluster()
mylines=N.arange(-20.,20.,.4)
mylineswidth=3
def psplotinit(output):
    file=output+"/vcps"
    ppgplot.pgbeg(file,1,1)
    ppgplot.pgpap(8.,1.)
    ppgplot.pgsch(1.7) #font size
    ppgplot.pgslw(7)  #line width

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
        #xbin[i]=scipy.stats.stats.median(x[nmin:nmax])
        #ybin[i]=scipy.stats.stats.median(y[nmin:nmax])
        xbin[i]=median(x[nmin:nmax])
        ybin[i]=median(y[nmin:nmax])
        #xbin[i]=N.average(x[nmin:nmax])
        #ybin[i]=N.average(y[nmin:nmax])
        #ybinerr[i]=scipy.stats.std(y[nmin:nmax])
    return xbin, ybin#, ybinerr

def biniterr(x,y,n):#bin arrays x, y into n bins, returning xbin,ybin
    nx=len(x)
    #x=N.array(x,'f')
    #y=N.array(y,'f')
    y=N.take(y,N.argsort(x))
    x=N.take(x,N.argsort(x))
    xbin=N.zeros(n,'f')
    ybin=N.zeros(n,'f')
    ybinerr=N.zeros(n,'f')
    for i in range(n):
        nmin=i*int(float(nx)/float(n))
        nmax=(i+1)*int(float(nx)/float(n))
        #xbin[i]=scipy.stats.stats.median(x[nmin:nmax])
        #ybin[i]=scipy.stats.stats.median(y[nmin:nmax])
        #xbin[i]=pylab.median(x[nmin:nmax])
        #ybin[i]=pylab.median(y[nmin:nmax])
        xbin[i]=N.average(x[nmin:nmax])
        ybin[i]=N.average(y[nmin:nmax])
        #ybinerr[i]=scipy.stats.std(y[nmin:nmax])/N.sqrt(1.*(nmax-nmin))
        ybinerr[i]=pylab.std(y[nmin:nmax])/N.sqrt(1.*(nmax-nmin))
        #ybinerr[i]=pylab.std(y[nmin:nmax])
    return xbin, ybin, ybinerr

def drawbinned(x,y,nbin):
    xbin,ybin=binit(x,y,nbin)
    #ppgplot.pgsci(2)
    ppgplot.pgline(xbin,ybin)

def drawbinnedsub0(x,y):
    ppgplot.pgsls(4)
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)

    drawbinned(xs,ys,5)
    ppgplot.pgsls(1)
def drawptsub0(x,y):
    defsci=ppgplot.pgqci()
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    #ppgplot.pgsci(4)
    ppgplot.pgpt(xs,ys,17)
    #ppgplot.pgsci(defsci)
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
    for i in range(len(y)):
        membids=list(c.sfmemberids[i])
        for j in membids:
            j=int(j)
            if (g.dv[j] < sig*c.sigma[i]):
                if (g.dr[j] < r*c.r200[i]):
                    if (g.Mabs[j] < mabscut):
                        dagn=float(g.agn[j])-float(agn)
                        if (agn < 0.1):
                            y[i]=y[i]+g.sfr[j]
                                    
    return y

def sumitold(sig,r,agn):
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    if agn < 1:
        for i in range(len(y)):
            y[i]=N.sum(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.ew > ewmin)  & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]) & (g.Mabs < mabscut) & (g.agn < 1) & (g.lHa > lHamin), g.sfr))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
    if agn > 0:#include agn in summed sfr
        for i in range(len(y)):
            y[i]=N.sum(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.ewr > ewmin)  & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]) & (g.Mabs < mabscut) & (g.lHa > lHamin), g.sfr))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr

    return y


def sumitmemb(sig,r,agn):
    y=N.zeros(len(c.z),'f')
    yave=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    for i in range(len(y)):
        ntot=0.
        membids=list(c.sfmemberids[i])
        for j in membids:
            j=int(j)
            if (g.memb[j] > 0):
                if (g.dv[j] < sig*c.sigma[i]):
                    if (g.dr[j] < r*c.r200[i]):
                        if (g.Mabs[j] < mabscut):
                            ntot=ntot+1.
                            dagn=float(g.agn[j])-float(agn)
                            if (dagn < 0.1):
                                y[i]=y[i]+g.sfr[j]
        try:
            yave[i]=float(y[i])/float(ntot)
        except ZeroDivisionError:
            yave[i]=0.
    return y,yave

def sumitmembold(sig,r,agn):
    y=N.zeros(len(c.z),'f')
    if agn < 1:
        for i in range(len(y)):
            ytemp=N.compress((abs(g.clusterid - c.id[i])< 1) & (g.ew > ewmin)   & (g.Mabs < mabscut) & (g.agn < 1) & (g.lHa > lHamin) & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]) & (g.memb > 0), g.sfr)#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
            y[i]=N.sum(ytemp)
            
    if agn > 0:#include agn in summed sfr
        for i in range(len(y)):
            ytemp=N.compress((abs(g.clusterid - c.id[i])< 1) & (g.ewr > ewmin)  & (g.Mabs < mabscut) & (g.lHa > lHamin) & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]) & (g.memb > 0), g.sfr)#keep and then sum g.sfr according to selection criteria on cluster, sfr, st
            y[i]=N.sum(ytemp)
            

    return y

def sumitbin(sig,r,agn):
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    for i in range(len(y)):
        membids=list(c.sfmemberids[i])
        for j in membids:
            j=int(j)
            if (g.dv[j] < sig*c.sigmabin[i]):
                if (g.dr[j] < r*c.r200bin[i]):
                    if (g.Mabs[j] < mabscut):
                        dagn=float(g.agn[j])-float(agn)
                        if (dagn < 0.1):
                            y[i]=y[i]+g.sfr[j]
    return y

def sumitbinold(sig,r,agn):#uses binned sig and r 
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    if agn < 1:
        for i in range(len(y)):
            ytemp=N.compress((abs(g.clusterid - c.id[i])< 1) & (g.ew > ewmin)   & (g.Mabs < mabscut) & (g.agn < 1) & (g.lHa > lHamin) & (g.dv < sig*c.sigmabin[i]) & (g.dr < r*c.r200bin[i]), g.sfr)#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
            y[i]=N.sum(ytemp)
            
    if agn > 0:#include agn in summed sfr
        for i in range(len(y)):
            ytemp=N.compress((abs(g.clusterid - c.id[i])< 1) & (g.ewr > ewmin)  & (g.Mabs < mabscut) & (g.lHa > lHamin) & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]), g.sfr)#keep and then sum g.sfr according to selection criteria on cluster, sfr, st
            y[i]=N.sum(ytemp)
            
    return y

def ngalmemb(sig,r):
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    for i in range(len(y)):
        membids=list(c.memberids[i])
        for j in membids:
            j=int(j)
            #if (g.memb[j] > 0):
	    if (g.dv[j] < sig*c.sigma[i]):
		if (g.dr[j] < r*c.r200[i]):
		    if (g.Mabs[j] < mabscut):
			y[i]=y[i]+1.
    return y

def ngalmembold(sig,r):
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    for i in range(len(y)):
        y[i]=float(len(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.Mabs < mabscut) & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]) & (g.memb > 0), g.sfr)))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
    return y

def nsfgalmemb(sig,r):
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    for i in range(len(y)):
        membids=list(c.sfmemberids[i])
        for j in membids:
            j=int(j)
            if (g.memb[j] > 0):
                if (g.dv[j] < sig*c.sigma[i]):
                    if (g.dr[j] < r*c.r200[i]):
                        if (g.Mabs[j] < mabscut):
                            y[i]=y[i]+1.
    return y


def nsfgalmembold(sig,r):#number of star-forming galaxies
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    for i in range(len(y)):
        y[i]=float(len(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.Mabs < mabscut) &  (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]) & (g.memb > 0) & (g.ew > ewmin) & (g.agn < 1) & (g.lHa > lHamin), g.sfr)))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
    return y

def sffracmemb(sig,r,agn):
    y=N.zeros(len(c.z),'f')
    sf=N.zeros(len(c.z),'f')
    tot=N.zeros(len(c.z),'f')
    for i in range(len(y)):
        membids=list(c.memberids[i])
        for j in membids:
            j=int(j)
            if (g.memb[j] > 0):
                if (g.dv[j] < sig*c.sigma[i]):
                    if (g.dr[j] < r*c.r200[i]):
                        if (g.Mabs[j] < mabscut):
                            if (g.ew[j] > -500.):#make sure EW data available
                                tot[i]=tot[i]+1.
                                if (g.ew[j] > ewmin):
                                    if (g.lHa[j] > lHamin):
                                        dagn=float(g.agn[j])-float(agn)
                                        if (dagn < 0.1):
                                            sf[i]=sf[i]+1.

        try:
            y[i]=(sf[i])/(tot[i])
        except ZeroDivisionError:
            y[i]=0.
    return y,sf,tot

def stellmassfracmemb(sig,r):#fraction of member w/stellar mass estimates
    y=N.zeros(len(c.z),'f')
    sf=N.zeros(len(c.z),'f')
    tot=N.zeros(len(c.z),'f')
    for i in range(len(y)):
        membids=list(c.memberids[i])
        for j in membids:
            j=int(j)
            if (g.memb[j] > 0):
                if (g.dv[j] < sig*c.sigma[i]):
                    if (g.dr[j] < r*c.r200[i]):
                        if (g.Mabs[j] < mabscut):
                            tot[i]=tot[i]+1.
                            if (g.stellarmass[j] > 0):
                                sf[i]=sf[i]+1.

        try:
            y[i]=(sf[i])/(tot[i])
        except ZeroDivisionError:
            y[i]=0.
    return y,sf,tot

def sffracmembold(sig,r,agn):
    y=N.zeros(len(c.z),'f')
    sf=N.zeros(len(c.z),'f')
    tot=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut

    if agn < 1:
        for i in range(len(y)):
            #if ((g.ew[i] > 4) & (g.lHa[i] > lHamin) & (g.agn < 1)):
            sf[i]=float(len(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.ew > ewmin) & (g.Mabs < mabscut) & (g.agn < 1) & (g.lHa > lHamin) & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]), g.sfr)))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
            tot[i]=float(len(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.Mabs < mabscut)  & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]), g.sfr)))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
            if tot[i] == 0:
                y[i]=0
            else:
                y[i]=float(sf[i])/float(tot[i])
    if agn > 0:#include agn in summed sfr
        for i in range(len(y)):
            sf[i]=float(len(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.ew > ewmin)   & (g.Mabs < mabscut) & (g.lHa > lHamin) & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]) & (g.memb > 0), g.sfr)))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
            tot[i]=float(len(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.Mabs < mabscut) & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]) & (g.memb > 0), g.sfr)))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
            if tot[i] == 0:
                y[i]=0
            else:
                y[i]=float(sf[i])/float(tot[i])

    return y,sf,tot

def sumitDn(sig,r,agn):
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    if agn < 1:
        for i in range(len(y)):
            y[i]=N.sum(N.compress((abs(g.clusterid - c.id[i])< 1)   & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]) & (g.Mabs < mabscut) & (g.agn < 1) & (g.Dn > -100.), g.Dn))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
    if agn > 0:#include agn in summed sfr
        for i in range(len(y)):
            y[i]=N.sum(N.compress((abs(g.clusterid - c.id[i])< 1)   & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]) & (g.Mabs < mabscut) & (g.Dn > -100.), g.Dn))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr

    return y


def sumitmembDn(sig,r,agn):
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    print len(g.Dn),len(g.dr)
    if agn < 1:
        for i in range(len(y)):
            y[i]=N.sum(N.compress((abs(g.clusterid - c.id[i])< 1)  & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]) & (g.Mabs < mabscut) & (g.agn < 1) & (g.Dn > -100.) & (g.memb > 0), g.Dn))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
    if agn > 0:#include agn in summed sfr
        for i in range(len(y)):
            y[i]=N.sum(N.compress((abs(g.clusterid - c.id[i])< 1)  & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]) & (g.Mabs < mabscut) & (g.Dn > -100.) & (g.memb > 0), g.Dn))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr

    return y

def sumitfixedr(sig,r,agn):
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    if agn < 1:
        for i in range(len(y)):
            y[i]=N.sum(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.ew > ewmin)  & (g.dv < sig*c.sigma[i]) & (g.dr < r) & (g.Mabs < mabscut) & (g.agn < 1) & (g.lHa > lHamin), g.sfr))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
    if agn > 0:#include agn in summed sfr
        for i in range(len(y)):
            y[i]=N.sum(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.ewr > ewmin)  & (g.dv < sig*c.sigma[i]) & (g.dr < r) & (g.Mabs < mabscut) & (g.lHa > lHamin), g.sfr))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr

    return y

def sumitfixedrv(sig,r,agn):
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    if agn < 1:
        for i in range(len(y)):
            y[i]=N.sum(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.ew > ewmin)  & (g.dv < sig) & (g.dr < r) & (g.Mabs < mabscut) & (g.agn < 1) & (g.lHa > lHamin), g.sfr))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
    if agn > 0:#include agn in summed sfr
        for i in range(len(y)):
            y[i]=N.sum(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.ewr > ewmin)  & (g.dv < sig) & (g.dr < r) & (g.Mabs < mabscut) & (g.lHa > lHamin), g.sfr))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr

    return y

def sumstellmass(sig,r):
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    for i in range(len(y)):
        membids=list(c.memberids[i])
        for j in membids:
            j=int(j)
            if (g.dv[j] < sig*c.sigma[i]):
                if (g.dr[j] < r*c.r200[i]):
                    if (g.Mabs[j] < mabscut):
                        if (g.stellarmass[j] > 0):
                            y[i]=y[i]+g.stellarmass[j]
                                    
    return y

def sumstellmassold(sig,r):
    y=N.zeros(len(c.z),'f')
    for i in range(len(y)):
        y[i]=N.sum(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.stellarmass > 0) & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]) & (g.Mabs < mabscut), g.stellarmass))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
    return y

def sumstellmassbin(sig,r):
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    for i in range(len(y)):
        membids=list(c.memberids[i])
        for j in membids:
            j=int(j)
            if (g.dv[j] < sig*c.sigmabin[i]):
                if (g.dr[j] < r*c.r200bin[i]):
                    if (g.Mabs[j] < mabscut):
                        if (g.stellarmass[j] > 0):
                            y[i]=y[i]+g.stellarmass[j]
                                    
    return y

def sumstellmassbinold(sig,r):
    y=N.zeros(len(c.z),'f')
    for i in range(len(y)):
        y[i]=N.sum(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.stellarmass > 0) & (g.dv < sig*c.sigmabin[i]) & (g.dr < r*c.r200bin[i]) & (g.Mabs < mabscut), g.stellarmass))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
    return y

def sumstellmassmemb(sig,r):
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    for i in range(len(y)):
        membids=list(c.memberids[i])
        for j in membids:
            j=int(j)
            if (g.memb[i] > 0):
                if (g.Mabs[j] < mabscut):
                    if (g.stellarmass[j] > 0):
                        y[i]=y[i]+g.stellarmass[j]
                                    
    return y


def sumstellmassmembold(sig,r):
    y=N.zeros(len(c.z),'f')
    for i in range(len(y)):
        y[i]=N.sum(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.stellarmass > 0) & (g.memb > 0) & (g.Mabs < mabscut), g.stellarmass))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
    return y

def getrichness(v,r,agn):
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    for i in range(len(y)):
        membids=list(c.sfmemberids[i])
        for j in membids:
            j=int(j)
            if (g.dv[j] < sig*c.sigma[i]):
                if (g.dr[j] < r*c.r200[i]):
                    if (g.Mabs[j] < mabscut):
                        dagn=float(g.agn[j])-float(agn)
                        if (dagn < 0.1):
                            y[i]=y[i]+g.sfr[j]
                                    
    return y


def richness(sig,r,agn,mabs):
    y=N.zeros(len(c.z),'f')
    for i in range(len(y)):
        membids=list(c.memberids[i])
        for j in membids:
            j=int(j)
            if (g.dv[j] < sig*c.sigma[i]):
                if (g.dr[j] < r*c.r200[i]):
                    if (g.Mabs[j] < mabs):
                        dagn=float(g.agn[j])-float(agn)
                        if (dagn < 0.1):
                            y[i]=y[i]+1.
                                    
    return y

def richnessmemb(sig,r,agn,mabs):
    y=N.zeros(len(c.z),'f')
    for i in range(len(y)):
        membids=list(c.memberids[i])
        for j in membids:
            j=int(j)
            if (g.memb[j] > 0):
                if (g.dv[j] < sig*c.sigma[i]):
                    if (g.dr[j] < r*c.r200[i]):
                        if (g.Mabs[j] < mabs):
                            dagn=float(g.agn[j])-float(agn)
                            if (dagn < 0.1):
                                y[i]=y[i]+1.
                                    
    return y

def countit(sig,r,agn,mabs):
    y=N.zeros(len(c.z),'f')
    for i in range(len(y)):
        membids=list(c.memberids[i])
        for j in membids:
            j=int(j)
            if (g.dv[j] < sig*c.sigma[i]):
                if (g.dr[j] < r*c.r200[i]):
                    if (g.Mabs[j] < mabs):
                        dagn=float(g.agn[j])-float(agn)
                        if (dagn < 0.1):
                            y[i]=y[i]+1.
                                    
    return y

def countitmemb(sig,r,agn,mabs):
    y=N.zeros(len(c.z),'f')
    for i in range(len(y)):
        membids=list(c.memberids[i])
        for j in membids:
            j=int(j)
            if (g.memb[j] > 0):
                if (g.dv[j] < sig*c.sigma[i]):
                    if (g.dr[j] < r*c.r200[i]):
                        if (g.Mabs[j] < mabs):
                            dagn=float(g.agn[j])-float(agn)
                            if (dagn < 0.1):
                                y[i]=y[i]+1.
                                    
    return y

def qsumit(yin,sig,r,agn,mabs):
    y=N.zeros(len(c.z),'f')
    for i in range(len(y)):
        membids=list(c.memberids[i])
        for j in membids:
            j=int(j)
            if (g.dv[j] < sig[i]):
                if (g.dr[j] < r[i]):
                    if (g.Mabs[j] < mabs):
                        dagn=float(g.agn[j])-float(agn)
                        if (dagn < 0.1):
                            if (g.memb[j] > 0):
                                y[i]=y[i]+yin[j]
                                    
    return y

def qsumitmemb(yin,sig,r,agn,mabs):
    y=N.zeros(len(c.z),'f')
    for i in range(len(y)):
        membids=list(c.memberids[i])
        for j in membids:
            j=int(j)
            if (g.memb[j] > 0):
                if (g.dv[j] < sig*c.sigma[i]):
                    if (g.dr[j] < r*c.r200[i]):
                        if (g.Mabs[j] < mabs):
                            dagn=float(g.agn[j])-float(agn)
                            if (dagn < 0.1):
                                y[i]=y[i]+yin[j]
                                    
    return y

def countmemb(sig,r,agn,mabs):
    y=N.zeros(len(c.z),'f')
    for i in range(len(y)):
        membids=list(c.memberids[i])
        for j in membids:
            j=int(j)
            if (g.memb[j] > 0):
                if (g.dv[j] < sig*c.sigma[i]):
                    if (g.dr[j] < r*c.r200[i]):
                        if (g.Mabs[j] < mabs):
                            dagn=float(g.agn[j])-float(agn)
                            if (dagn < 0.1):
                                y[i]=y[i]+1
                                    
    return y


def richnessold(sig,r,agn,mabs):
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    if agn < 1:
        for i in range(len(y)):
            y[i]=len(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]) & (g.Mabs < mabs) & (g.agn < 1), g.sfr))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
    if agn > 0:#include agn in summed sfr
        for i in range(len(y)):
            y[i]=len(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]) & (g.Mabs < mabs), g.sfr))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
    return y
def richnessmembold(sig,r,agn,mabs):
    y=N.zeros(len(c.z),'f')
    #print "mabscut = ",mabscut
    if agn < 1:
        for i in range(len(y)):
            y[i]=len(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]) & (g.Mabs < mabs) & (g.memb > 0) & (g.agn < 1), g.sfr))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr
    if agn > 0:#include agn in summed sfr
        for i in range(len(y)):
            y[i]=len(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.dv < sig*c.sigma[i]) & (g.dr < r*c.r200[i]) & (g.Mabs < mabs) & (g.memb > 0), g.sfr))#keep and then sum g.sfr according to selection criteria on cluster, sfr, stellarmass, dv, and dr

    return y

#def dospear(x,y):
#    (a,b)=scipy.stats.spearmanr(x,y)
#    print "rank correl = %6.3f %8.f" % (a,b)
#    #(a,b)=scipy.stats.kendalltau(x,y)
#q    #print "kendall tau = %6.3f %6.5f" % (a,b)

def dostats():
    print "total number of clusters in final sample = ",len(c.z)
    print "fraction per cluster with Kauffmann stellar mass"
    #med=scipy.median(c.kauffmann)
    med=pylab.median(c.kauffmann)
    ave=N.average(c.kauffmann)
    #std=scipy.stats.std(c.kauffmann)
    std=pylab.std(c.kauffmann)
    print "ave = %6.3f %6.3f, med = %6.3f" % (ave,std,med)

    #agn contamination

    #s1=sumit(3.,2.,0) #total sfr w/out agn
    #s2=sumit(3.,2.,1) #total sfr w/agn
    #agnfrac=(s2-s1)/(s1)
    #contamination=N.average(agnfrac)
    ##std=scipy.stats.std(agnfrac)
    ##med=scipy.median(agnfrac)
    #std=pylab.std(agnfrac)
    #med=pylab.median(agnfrac)
    #n1=N.sum(g.agn)
    #n2=len(g.agn)
    #frac=float(n1)/float(n2)
    #print "AGN contamination by number w/in 3sig, 2R = ",frac
    #print "AGN contamination in terms of total SFR w/in 3sig,2R= ",contamination,std,med

    #s1=sumit(3.,0.5,0) #total sfr w/out agn
    #s2=sumit(3.,0.5,1) #total sfr w/agn
    #agnfrac=(s2-s1)/(s1)
    #contamination=N.average(agnfrac)
    ##std=scipy.stats.std(agnfrac)
    ##med=scipy.median(agnfrac)
    #std=pylab.std(agnfrac)
    #med=pylab.median(agnfrac)
    #print "AGN contamination in terms of total SFR w/in 3sig,0.5R= ",contamination,std,med
    #y=N.take(c.sumsfr2,N.argsort(c.sumsfr2))
    #x=N.take(c.id,N.argsort(c.sumsfr2))
    #z=N.take(c.z,N.argsort(c.sumsfr2))
#    print "Clusters w/lowest tot SFR"
#    for i in range(9):
#        print x[i],z[i],y[i]
    #temp=N.array(len(c.sumsfrhiz),'f')
    #temp=c.sfr05/c.sumsfrhiz
    #med=median(temp)
    #print "0.5Rvir, 6/3 sigma, med = ",med
    #med=N.average(temp)
    #std=pylab.std(temp)
    #print "0.5Rvir, 6/3 sigma, ave ",med, "+/-",std

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

def plotmclz():
    x=c.z
    y=N.log10(c.mass)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,.2,2.1,0,20)
    ppgplot.pglab("z","M\dcl \u (10\u14 \d M\d\(2281)\u)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)#blue
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
    ppgplot.pgsci(1)
    ppgplot.pgpt(xs,ys,17)
    
def plotmclsigma():
    x=c.sigma
    y=N.log10(c.mass)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    xmax=sigmamax+200.
    xmin=sigmamin-200.
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,.2,2.2,0,20)
    ppgplot.pglab("\gs (km/s)","M\dcl \u (10\u14 \d M\d\(2281)\u)","")
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
    ymax=sigmamax+100.
    ymin=sigmamin-150.

    ppgplot.pgenv(zminp,zmaxp,ymin,ymax,0,0)
    ppgplot.pglab("z","\gs (km/s)","")
    ppgplot.pgpt(x,y,17)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsls(4)#
    #xs=N.compress(c.sub < 1,x)
    #ys=N.compress(c.sub < 1,y)
    #drawbinned(xs,ys,5)
    #ppgplot.pgsci(1)
    #ppgplot.pgpt(xs,ys,17)

def plottotlumsigma():
    x=c.sigma
    y=c.totlum
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,0)
    ppgplot.pglab("\gs (km/s)","Total R Lum","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)#blue
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
    ppgplot.pgsci(1)
    ppgplot.pgpt(xs,ys,17)
def plottotrlumsigma():
    x=c.sigma
    y=c.totrlum
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    xmin=sigmaminp
    xmax=sigmamaxp
    ymin=totrlumminp
    ymax=totrlummaxp
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,0)
    ppgplot.pglab("\gs (km/s)","Total L\dR\u (10\u9\d L\d\(2281))\u)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    #ppgplot.pgsci(4)#blue
    ppgplot.pgsls(4)#blue
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
    ppgplot.pgsci(1)
    ppgplot.pgpt(xs,ys,17)

def plotsigmatotrlum():
    y=N.log10(c.sigma)
    print "min totrlum = ",min(c.totrlum)
    print "max totrlum = ",max(c.totrlum)
    x=N.log10(c.totrlum/100.)#convert to units of 10^11
    ymax=max(x)
    ymin=min(x)
    xmax=max(y)
    xmin=min(y)
    ymin=N.log10(sigmaminp)
    ymax=N.log10(sigmamaxp)
    #xmin=N.log10(totrlumminp)
    #xmax=N.log10(totrlummaxp)
    xmin=-1.
    xmax=2.

    #xmin=0.
    #xmax=5.
    #ymin=0.
    #ymax=5.
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,30)
    ppgplot.pglab("L\dR\u (10\u11\d L\d\(2281)\u)","\gs (km/s)","")
    ppgplot.pgpt(x,y,21)
    ppgplot.pgsci(2)
    #drawbinned(x,y,5)
    #ppgplot.pgsci(4)#blue
    ppgplot.pgsls(4)#dotted
    xs=N.compress(c.gauss > 0,x)
    ys=N.compress(c.gauss > 0,y)
    #drawbinned(xs,ys,5)
    ppgplot.pgsci(1)
    ppgplot.pgpt(xs,ys,17)
    #x=N.arange(100.,1500.,100.)
    x=N.arange(-2.,5.,.1)
    #y=.32*x+2.102
    y=.5*x+1.75+1.
    #y=.2*x+2.42
    #x=x-2.
    ysig = 0.5*(N.log10(c.totrlum/100.))+2.75
    ppgplot.pgsls(4)
    ppgplot.pgsci(4)#blue
    ppgplot.pgline(x,y)
    ppgplot.pgsci(2)
    #ppgplot.pgpt((N.log10(c.totrlum/100.)),ysig,12)
def plotsigmarichness():
    y=c.sigma
    x=c.richness
    ymax=max(x)
    ymin=min(x)
    xmax=max(y)
    xmin=min(y)
    ymin=sigmaminp
    ymax=sigmamaxp
    xmin=richmin
    xmax=richmax
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,0)
    ppgplot.pglab("Richess","\gs (km/s)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    #ppgplot.pgsci(4)#blue
    ppgplot.pgsls(4)#blue
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
    ppgplot.pgsci(1)
    ppgplot.pgpt(xs,ys,17)

def plotrichnesstotrlum():
    x=c.totrlum
    y=c.richness
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ymax=richmax
    ymin=richmin
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,0)
    ppgplot.pglab("Tot R Lum (Miller)","Richness","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsls(4)#dotted
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
    ppgplot.pgsci(1)
    ppgplot.pgpt(xs,ys,17)
    
def plottotlumtotrlum():#compare my tot lum w/Chris Miller's tot lum
    x=c.totrlum
    y=c.totlum
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,0)
    ppgplot.pglab("Total R Lum (Miller)","Total R Lum","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)#blue
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
    ppgplot.pgsci(1)
    ppgplot.pgpt(xs,ys,17)
def plottotlum1Mpcsigma():
    x=c.sigma
    y=c.totlum1Mpc
    xmax=max(x) 
    xmin=min(x) 
    ymax=max(y) 
    ymin=min(y) 
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,0)
    ppgplot.pglab("\gs (km/s)","Total R Lum (dr<1Mpc)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)#blue
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
    ppgplot.pgsci(1)

def plotsigmarvir():
    x=c.rvir
    y=c.sigma
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ymin=sigmamin-200.
    ymax=sigmamax+200.
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(1.,4.2,ymin,ymax,0,0)
    ppgplot.pglab("Rvir (Mpc)","\gs (km/s)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)#blue
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
    ppgplot.pgsci(1)
    ppgplot.pgpt(xs,ys,17)
def plottotlumrvir():
    x=c.rvir
    y=c.totlum
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,0)
    ppgplot.pglab("Rvir (Mpc)","Total R Lum","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)#blue
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)

    ppgplot.pgsci(1)
def plottotrlumrvir():
    x=c.rvir
    y=c.totrlum
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    xmin=1.
    xmax=4.
    ymin=totrlumminp
    ymax=totrlummaxp
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,0)
    ppgplot.pglab("Rvir (Mpc)","Total L\dR\u (10\u9\d L\d\(2281)\u)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsls(4)#dotted
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)

    ppgplot.pgsci(1)

def plottotrlumz():
    x=c.z
    y=c.totrlum
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    xmin=zminp
    xmax=zmaxp
    ymin=totrlumminp
    ymax=totrlummaxp
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,0)
    ppgplot.pglab("z","Total L\dR\u (10\u9\d L\d\(2281))\u)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsls(4)#dotted
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
    ppgplot.pgsci(1)
    ppgplot.pgpt(xs,ys,17)
def plotsubz():
    x=c.z
    y=c.sub
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,-.1,4.1,0,0)
    ppgplot.pglab("z","Substructure (Miller SCF)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    #drawbinned(x,y,5)
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
    ppgplot.pgenv(zminp,zmaxp,0.,250,0)
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
    dv=3
    dr=1
    x=c.z
    y=richnessmemb(dv,dr,0,mabs)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,0.,150,0)
    ppgplot.pglab("z","Richness(M\dR\u)","")
    #ppgplot.pgpt(x,y,3)
    drawbinned(x,y,5)
    ppgplot.pgsci(2)
    mabs=mabs+0.5
    y=richnessmemb(dv,dr,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(3)
    mabs=mabs+0.5
    y=richnessmemb(dv,dr,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)
    mabs=mabs+0.5
    y=richnessmemb(dv,dr,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(5)
    mabs=mabs+0.5
    y=richnessmemb(dv,dr,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(6)
    mabs=mabs+0.5
    y=richnessmemb(dv,dr,0,mabs)
    drawbinned(x,y,5)
    ppgplot.pgsci(7)
    mabs=mabs+0.5
    y=richnessmemb(dv,dr,0,mabs)
    #drawbinned(x,y,5)
    ppgplot.pgsci(8)
    mabs=mabs+0.5
    y=richnessmemb(dv,dr,0,mabs)
    #drawbinned(x,y,5)
    ppgplot.pgsci(1)
    mabs=mabscut
    y=richnessmemb(dv,dr,0,mabs)
    ppgplot.pgsls(4)
    ppgplot.pgslw(4)
    drawbinned(x,y,5)
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
    (y,t)=sumitmemb(3,nr,0)
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

def plotDnmcl():
    y=sumitmembDn(3,3,0)/c.richness
    x=N.log10(c.mass)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,2,ymin,ymax,0,10)
    ppgplot.pglab("Mcl","\gSD\dn\u","")
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
    ppgplot.pgenv(zminp,zmaxp,0,200,0)
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

def plotsfrmclsigmaalls(): #plot all cuts on SFR versus z
    #x=(c.z)
    nbin=5
    y=sumit(1.,nr,0)
    y=y/c.mass
    xin=c.sigma
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(sigmaminp,sigmamaxp,0,20,0)
    ppgplot.pglab("\gs (km/s)","\gSSFR/M\dcl \u (M\d\(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(2)
    y=(sumit(2.,nr,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,20)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(4)
    y=(sumit(3.,nr,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,11)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(3)
    y=(sumit(6.,nr,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,5)
    ppgplot.pgline(x,y)

    ppgplot.pgsci(1)

def plotsfrmclsigmaallr(): #plot all cuts on SFR versus z
    #x=(c.z)
    nbin=5
    sig=3.
    xin=c.sigma
    #y=N.log10(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(sigmaminp,sigmamaxp,0,50,0)
    ppgplot.pglab("\gs (km/s)","\gSSFR/M\dcl \u (M\d\(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","")
    y=sumit(sig,0.5,0)
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    ppgplot.pgpt(x,y,3)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(2)
    y=(sumit(sig,1.,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,20)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(4)
    y=(sumit(sig,2.,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,11)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(3)
    y=(sumit(sig,3.,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,5)
    ppgplot.pgline(x,y)

    ppgplot.pgsci(1)

def plotsfrmclzalls(): #plot all cuts on SFR versus z
    #x=(c.z)
    nbin=5
    y=sumit(1.,nr,0)
    y=y/c.mass
    xin=c.z
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,0,20,0)
    ppgplot.pglab("z","\gSSFR/M\dcl \u (M\d\(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(2)
    y=(sumit(2.,nr,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,20)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(4)
    y=(sumit(3.,nr,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,11)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(3)
    y=(sumit(6.,nr,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,5)
    ppgplot.pgline(x,y)

    ppgplot.pgsci(1)

def plotsfrmclzallr(): #plot all cuts on SFR versus z
    #x=(c.z)
    nbin=5
    sig=3.
    xin=c.z
    #y=N.log10(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,0,50,0)
    ppgplot.pglab("z","\gSSFR/M\dcl \u (M\d\(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","")
    y=sumit(sig,0.5,0)
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    ppgplot.pgpt(x,y,3)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(2)
    y=(sumit(sig,1.,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,20)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(4)
    y=(sumit(sig,2.,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,11)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(3)
    y=(sumit(sig,3.,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,5)
    ppgplot.pgline(x,y)

    ppgplot.pgsci(1)

def plotsfrmcltotrlumalls(): #plot all cuts on SFR versus z
    #x=(c.z)
    nbin=5
    y=sumit(1.,nr,0)
    y=y/c.mass
    xin=c.totrlum
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    xmin=min(xin)
    xmax=max(xin)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,0,20,0)
    ppgplot.pglab("Total R Lum (Miller)","\gSSFR/M\dcl \u (M\d\(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(2)
    y=(sumit(2.,nr,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,20)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(4)
    y=(sumit(3.,nr,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,11)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(3)
    y=(sumit(6.,nr,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,5)
    ppgplot.pgline(x,y)

    ppgplot.pgsci(1)

def plotsfrmcltotrlumallr(): #plot all cuts on SFR versus z
    #x=(c.z)
    nbin=5
    sig=3.
    y=sumit(sig,0.5,0)
    y=y/c.mass
    xin=c.totrlum
    x,y=binit(xin,y,nbin)
    xmin=min(xin)
    xmax=max(xin)
    #y=N.log10(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,0,50,0)
    ppgplot.pglab("Total R Lum (Miller)","\gSSFR/M\dcl \u (M\d\(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(2)
    y=(sumit(sig,1.,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,20)
    ppgplot.pgline(x,y)

    ppgplot.pgsls(4)
    y=(sumit(2.,1.,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,20)
    ppgplot.pgline(x,y)

    y=(sumit(1.,1.,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,20)
    ppgplot.pgline(x,y)

    ppgplot.pgsls(1)
    ppgplot.pgsci(4)
    y=(sumit(sig,2.,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,11)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(3)
    y=(sumit(sig,3.,0))
    y=y/c.mass
    x,y=binit(xin,y,nbin)
    #y=N.log10(y)
    ppgplot.pgpt(x,y,5)
    ppgplot.pgline(x,y)

    ppgplot.pgsci(1)

def plotsfrz08():
    #y=(N.compress((c.z > z08min) & (c.z < z08max),c.sumsfr2))
    #x=N.compress((c.z > z08min) & (c.z < z08max),c.z)
    y=(N.compress((c.super > 0),c.sumsfr2))
    x=N.compress((c.super > 0),c.z)
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
    (y,t)=sumitmemb(3,3,0)
    x=N.log10(c.mass)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,1.6,-5,200,10)
    ppgplot.pglab("M\dcl \u (10\u14 \d M\d\(2281)\u)","\gSSFR (M\d\(2281)\u yr\u-1\d)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)#blue
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
    ppgplot.pgsci(1)

def plotsfrtotrlum():
    #y=c.sumsfr2
    (y,t)=sumitmemb(3,3,0)
    x=c.totrlum
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,-5,200,10)
    ppgplot.pglab("Tot R Lum (10\u9 \d L\d\(2281)\u)","\gSSFR (M\d\(2281)\u yr\u-1\d)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)#blue
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
    ppgplot.pgsci(1)

def plotavesfrmcl():
    y=c.avesfr
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
    ppgplot.pgenv(zminp,zmaxp,richmin,richmax,0)
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
    xmin=sigmamin-200.
    xmax=sigmamax+200.
    #ymin=richmin
    #ymax=richmax
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,richmin,richmax,0)
    ppgplot.pglab("\gs (km/s)","Richness","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)#blue
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)

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

def plothistsffrac():
    ppgplot.pgbeg("sffrachist.ps/vcps",1,1)  #color port.
    ppgplot.pgpap(8.,1.25)
    ppgplot.pgpage
    ppgplot.pgsch(1.2) #font size
    ppgplot.pgslw(4)   #line width

    # 1st panel with symbols w/ stddev errorbars


    x1=.1
    x2=.95
    y1=.15
    y2=.85
    ppgplot.pgsvp(x1,x2,y1,y2)  #sets viewport
    xmin=-.1
    xmax=1.1
    ymin=-1.
    ymax=40.
    ybig=10.
    ysub=5
    if cdir.find('fieldDR3') > -1:
        ymax=300.#set min ngal to 30 for clusters, 1 for field
        ybig=50.
        ysub=5
    if cdir.find('fieldDR4') > -1:
        ymax=300.#set min ngal to 30 for clusters, 1 for field
        ybig=50.
        ysub=5
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('bcnst',.5,5,'bcvnst',ybig,ysub)  #tickmarks and labeling
    ppgplot.pgmtxt('b',2.5,0.5,0.5,"Fraction of SF Galaxies")    #xlabel
    ppgplot.pgmtxt('l',2.1,0.5,0.5,"Number")

    x=c.sffrac05
    dmax=1.
    dmin=0.
    nx=len(x)
    ppgplot.pghist(nx,x,dmin,dmax,10,5)

    ppgplot.pgsci(2)
    #ppgplot.pgsls(2)
    x=c.sffrac1
    nx=len(x)
    ppgplot.pghist(nx,x,dmin,dmax,10,5)

    ppgplot.pgsci(4)
    #ppgplot.pgsls(4)
    x=c.sffrac2
    nx=len(x)
    ppgplot.pghist(nx,x,dmin,dmax,10,5)

    ppgplot.pgsci(1)
    ppgplot.pgsls(1)

    ppgplot.pgend()
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
    #y=scipy.stats.histogram2(x,bins)
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
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(-24.5,-17.5,-0.03,1.03,0)
    ppgplot.pglab("M\dR\u","Cumulative SFR","")

    sfr=N.compress((g.clusterz < .075) & (g.Mabs < 0),g.sfr)
    totsfr=N.sum(sfr)
    mabs=N.compress((g.clusterz < .075) & (g.Mabs < 0),g.Mabs)
    bins=N.arange(min(mabs),max(mabs),0.1,'f')
    nbins=len(bins)
    sfrpart=N.zeros(len(bins),'f')
    for i in range(nbins):
        sfrpart[i]=N.sum(N.compress(mabs < bins[i],sfr))/totsfr
    bins=N.array(bins,'f')
    sfrpart=N.array(sfrpart,'f')
    ppgplot.pgsci(2)
    ppgplot.pgline(bins,sfrpart)

    sfr=N.compress((g.clusterz > .075) & (g.Mabs < 0),g.sfr)
    totsfr=N.sum(sfr)
    mabs=N.compress((g.clusterz > .075) & (g.Mabs < 0),g.Mabs)
    bins=N.arange(min(mabs),max(mabs),0.1,'f')
    nbins=len(bins)
    sfrpart=N.zeros(len(bins),'f')
    for i in range(nbins):
        sfrpart[i]=N.sum(N.compress(mabs < bins[i],sfr))/totsfr
    bins=N.array(bins,'f')
    sfrpart=N.array(sfrpart,'f')
    ppgplot.pgsci(1)
    ppgplot.pgsls(2)
    ppgplot.pgline(bins,sfrpart)

    #print bins
    #print sfrpart
    #ppgplot.pgpt(bins,sfrpart,3)
    #print nx, bins, nfrac
    ppgplot.pgsci(1)
    x=N.array([mabsorig,mabsorig])
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
    xmin=0.
    xmax=3.
    #ymin=-3000.
    #ymax=3000.
    ymin=-2.5
    ymax=2.5

    outfile=open('dvr-phys.dat','w')
    for i in range(len(g.dr)):
        outfile.write("%8.3f %8.3f \n" % (g.dr[i],g.dv[i]))
    outfile.close()

    sigmamin=500.
    t=g.dv/g.clustersigma
    y=N.compress(g.clustersigma > sigmamin,t)
    t=(g.dr)/g.clusterrvir
    x=N.compress(g.clustersigma > sigmamin,t)
    


    #x=N.compress(g.galsub < 1,g.dr)/N.compress(g.galsub < 1,g.clusterr200)
    #x=N.compress(g.galsub < 1,g.dr)/N.compress(g.galsub < 1,g.clusterrvir)
    #x=N.compress(g.galsub < 1,g.dr)
    outfile=open('dvr.dat','w')
    for i in range(len(x)):
        outfile.write("%8.3f %8.3f \n" % (x[i],y[i]))
    outfile.close()
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    #ppgplot.pgenv(0,3.1,-6.5,6.5,0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,0)
    ppgplot.pglab("R/R\d200 \u","\gDv/\gs","")
    #ppgplot.pglab("R","\gD v","")
    ppgplot.pgsch(0.9)
    ppgplot.pgpt(x,y,17)
    #ncontbin=10.
    #xcontour=N.array(xmin,xmax,(xmax-xmin)/ncontbin)
    #ycontour=N.array(ymin,ymax,(ymax-ymin)/ncontbin)
    #A=N.zeros(((ncontbin),(ncontbin)),'f')
    #xbinnumb=N.array(len(x),'f')
    #ybinnumb=N.array(len(x),'f')
    #x1=N.compress((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax),x)
    #y1=N.compress((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax),y) 
    #x=x1
    #y=y1
    #xbinnumb=((x-xmin)*ncontbin/(xmax-xmin))#calculate x  bin number for each point 
    #ybinnumb=((y-ymin)*ncontbin/(ymax-ymin))#calculate x  bin number for each point 
    #binID = zip(xbinnumb,ybinnumb)
    #for (i1,j1) in binID:
    #    i=int(i1)
    #    j=int(j1)
    #    A[j,i]=A[j,i]+1#have to switch i,j in A[] to make contours look right
    #pylab.figure()
    #pylab.contour(xcontour,ycontour,A)
    #pylab.xlabel('dr')
    #pylab.ylabel('dv')
    #savefig('dvvsrpylab.ps')
    #IDIM=int(ncontbin)
    #JDIM=int(ncontbin)
    #I1=0
    #I2=IDIM-1
    #J1=0
    #J2=JDIM-1
    #print IDIM,I1,I2,JDIM, J1,J2
    #TR=N.array([xmin,(0.5*(xmax-xmin)/ncontbin),0.,ymin,0.,(0.5*(ymax-ymin)/ncontbin)],'f')
    #TR=N.array([xmin,(1.*(xmax-xmin)/ncontbin),0.,ymin,0.,(1.*(ymax-ymin)/ncontbin)],'f')
    #amin=min(min(A))
    #amax=max(max(A))
    #ncont=5.
    #contstep=(amax-amin)/(ncont-1.)
    #print amin,amax,contstep
    #C=N.arange(amin,amax,contstep,'f')
    #NC=len(C)
    #ppgplot.pgsci(2)
    #ppgplot.pgsls(1)
    #ppgplot.pgcont_s(A, NC, C,xmin,ymin,xmax,ymax)
    #sfr = N.compress(g.galsub < 1,g.sfr)
    #y=N.compress(sfr > 3, y)
    #x=N.compress(sfr > 3, x)
    #ppgplot.pgsci(3)
    #ppgplot.pgpt(x,y,1)
    ppgplot.pgsci(2)
    ppgplot.pgslw(8)
    x=N.array([0,1.],'f')
    y=N.array([2,2],'f')
    #ppgplot.pgline(x,y)
    y=-1.*y
    #ppgplot.pgline(x,y)
    x=N.array([1.,3.],'f')
    y=N.array([1.,1.],'f')
    #ppgplot.pgline(x,y)
    y=-1.*y
    #ppgplot.pgline(x,y)
    ppgplot.pgsci(1)
    ppgplot.pgslw(4)

def plotdvvsrphys():
    xmin=-0.05
    xmax=2.
    ymin=-2000.
    ymax=2000.
    #ymin=-3.
    #ymax=3.
    sigmamin=400.
    y=N.compress(g.clustersigma > sigmamin,g.dv)
    t=(g.dr)/g.clusterr200
    x=N.compress(g.clustersigma > sigmamin,t)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,0)
    ppgplot.pglab("R/R\d200\u","\gDv","")
    ppgplot.pgsch(.9)
    ppgplot.pgpt(x,y,17)
    ppgplot.pgsci(2)
    ppgplot.pgslw(8)
    x=N.array([0,1.],'f')
    y=N.array([2,2],'f')
    #ppgplot.pgline(x,y)
    y=-1.*y
    #ppgplot.pgline(x,y)
    x=N.array([1.,3.],'f')
    y=N.array([1.,1.],'f')
    #ppgplot.pgline(x,y)
    y=-1.*y
    #ppgplot.pgline(x,y)
    ppgplot.pgsci(1)
    ppgplot.pgslw(4)

    
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
    drawbinned(x,y,5)
    ppgplot.pgsci(4)#blue
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
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
    ppgplot.pgenv(0,xmax,-.02,2.2,0,10)
    ppgplot.pglab("M\dcl \u (10\u14 \d M\d\(2281)\u)","\gSSFR/\gSM\d* \u(M\d\(2281)\u yr\u-1\d/10\u10 \d M\d\(2281)\u)","")
    ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    xs=N.compress(c.super > 0,x)
    ys=N.compress(c.super > 0,y)
    ppgplot.pgpt(xs,ys,3)

    ppgplot.pgsci(4)#blue
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)

    ppgplot.pgsci(1)


def plotsfrpmclsigma():
    x=c.sigma
    y=(c.sfrmass)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(300,1500,0,1000,0)#match hiz axes
    ppgplot.pglab("\gs (km s\u-1\d)","\gSSFR/M\dcl \u (M\d\(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)#blue
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)

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
    ppgplot.pgsci(4)
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
    ppgplot.pgsci(1)

def plotsfrtotrlumsigma():
    x=c.sigma
    y=N.log10(c.sumsfr/c.totrlum)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(300,1500,-1,2,0,20)#match hiz axes
    ppgplot.pglab("\gs (km s\u-1\d)","\gSSFR/L\dR tot\u(M\d\(2281)\u yr\u-1\d/10\u9\dL\d\(2281)\u )(Miller)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
    ppgplot.pgsci(1)

def plotsfrtotlumrvir():
    x=c.rvir
    y=(c.sfrtotlum)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0)#match hiz axes
    ppgplot.pglab("\gs (km s\u-1\d)","\gSSFR/L\dR tot\u(M\d\(2281)\u yr\u-1\d/)","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
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
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
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
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
    ppgplot.pgsci(1)

def plotsfrpmclmcl():
    #y=N.log10(c.sfrmass)
    y=(c.sfrmass)
    x=N.log10(c.mass)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,2,sfrmclminp,sfrmclmaxp,0,10)#match hiz axes
    ppgplot.pglab("M\dcl \u (10\u14 \d M\d\(2281) \u)","\gSSFR/M\dcl \u (M\d\(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","")
    ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(2)
    #y=N.log10(c.sfrmassbalogh)
    #ppgplot.pgpt(x,y,3)
    #ppgplot.pgsci(4)
    #y=N.log10(c.sfrmasshiz)
    #ppgplot.pgpt(x,y,3)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)#blue
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)

    ppgplot.pgsci(1)


def plotsfrpmclmcl08():
    #y=N.log10(c.sfrmass)
    y=(N.compress((c.super > 0),c.sfrmass))

    x=N.log10(N.compress((c.super > 0),c.mass))
    sub=(N.compress((c.super > 0),c.sub))
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,2,sfrmclminp,sfrmclmaxp,0,10)#match hiz axes
    ppgplot.pglab("M\dcl \u (10\u14 \d M\d\(2281) \u)","\gSSFR/M\dcl \u (M\d\(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","")
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(2)
    #y=N.log10(c.sfrmassbalogh)
    #ppgplot.pgpt(x,y,3)
    #ppgplot.pgsci(4)
    #y=N.log10(c.sfrmasshiz)
    #ppgplot.pgpt(x,y,3)
    drawbinned(x,y,3)
    ppgplot.pgsci(4)#blue
    xs=N.compress(sub < 1,x)
    ys=N.compress(sub < 1,y)
    drawbinned(xs,ys,5)

    x=N.log10(kodamamasshiz)
    y=sfrmhiz
    ppgplot.pgsci(4)
    ppgplot.pgpt(x,y,13)
    ppgplot.pgsci(1)

def plotsffracmcl08():
    #y=N.log10(c.sfrmass)
    #y=(N.compress((c.z > z08min) & (c.z < z08max),c.sffrac))
    #x=N.log10(N.compress((c.z > z08min) & (c.z < z08max),c.mass))
    y=(N.compress((c.super > 0),c.sffrac))
    x=N.log10(N.compress((c.super > 0),c.mass))
    sub=(N.compress((c.super > 0),c.sub))
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
    ppgplot.pgsci(4)#blue
    xs=N.compress(sub < 1,x)
    ys=N.compress(sub < 1,y)
    drawbinned(xs,ys,5)

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
    ppgplot.pgenv(zminp,zmaxp,-.02,2.2,0)
    ppgplot.pglab("z","\gSSFR/\gSM\d* \u(M\d\(2281)\u yr\u-1\d/10\u10 \d M\d\(2281)\u)","")
    ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    xs=N.compress(c.super > 0,x)
    ys=N.compress(c.super > 0,y)
    ppgplot.pgpt(xs,ys,3)

    ppgplot.pgsci(4)
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)

    ppgplot.pgsci(1)

def plotsfrpmclz():
    #y=(c.sfrmass)
    (y,t)=sumitmemb(6,3,0)
    y=y/c.mass
    x=c.z
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,sfrmclminp,sfrmclmaxp,0,0)
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
    (x,t)=sumitmemb(6,3,0)
    x=x/c.mass
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(-2.,2.,sfrmclminp,sfrmclmaxp,0)
    ppgplot.pglab("\gS SFR/M\d cl \u (M\d \(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u )","SF FRAC","")
    ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(2)
    #y=N.log10(c.sfrmassbalogh)
    #ppgplot.pgpt(x,y,3)
    #ppgplot.pgsci(4)
    #y=N.log10(c.sfrmasshiz)
    #ppgplot.pgpt(x,y,3)
    drawbinned(x,y,5)
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
    ppgplot.pgsci(1)
def plotavesfrz(): #butcher oemler plot
    y=(c.avesfr)
    x=c.z
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(zminp,zmaxp,-0.2,10,0)
    ppgplot.pglab("z","Median SFR of SF Galaxies","")
    ppgplot.pgpt(x,y,20)
    ppgplot.pgsci(2)
    drawbinned(x,y,5)
    ppgplot.pgsci(4)
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
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
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
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
    ppgplot.pgsci(4)
    xs=N.compress(c.sub < 1,x)
    ys=N.compress(c.sub < 1,y)
    drawbinned(xs,ys,5)
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
    ppgplot.pgsci(4)#blue
    y=N.compress(c.sub < 1,c.sffrac)
    x=N.compress(c.sub < 1,c.sigma)
    drawbinned(x,y,5)

    ppgplot.pgsci(1)

def plotsffracmcl(): #butcher oemler plot
    y=(c.sffrac)
    x=N.log10(c.mass)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
def plotngalsigmapylab():
    x=(c.sigma)
    y05=c.ngal05#sumitmemb(3.,0.5,0)
    y1=c.ngal1#sumitmemb(3.,1.,0)
    y2=c.ngal2#sumitmemb(3.,2.,0)
    fx=(f.sigma)
    fy05=f.ngal05#sumitmemb(3.,0.5,0)
    fy1=f.ngal1#sumitmemb(3.,1.,0)
    fy2=f.ngal2#sumitmemb(3.,2.,0)
    sx=(mil.sigma)
    sy05=mil.ongal05#sumitmemb(3.,0.5,0)
    sy1=mil.ongal1#sumitmemb(3.,1.,0)
    sy2=mil.ongal2#sumitmemb(3.,2.,0)

    pylabsubplot1(x,y05,y1,y2,fx,fy05,fy1,fy2)

    pylab.subplot(312)
    pylab.ylabel(r'$\rm{N_{gal}}$',fontsize=24,fontweight='bold')
    pylab.savefig("ngalmembsigma3.eps")

    print "Ngal Deprojection Check"
    pylab.cla()
    pylab.clf()
    #pylab.plot(x,y1,'k.')
    pylabsubplot11(x,y05,y1,y2,fx,fy05,fy1,fy2,sx,sy05,sy1,sy2)

    pylab.ylabel(r'$\rm{N_{gal}}$',fontsize=24,fontweight='bold')
    pylab.savefig("ngalmembsigma3r200.eps")



def pylabsubplot2(x,y05,y1,y2,fx,fy05,fy1,fy2):#ratio of y/x^3 versus x
    nbinscale=3
    xl=N.arange(100.,5000.,100.)

    pylab.cla()
    pylab.clf()
    xminsig=2.6
    xmaxsig=3.1
    symheight=3.
    sigmin=400.
    sigmax=1200.
    nbin=7

    ymin=0.
    ymax=2.5
    xmin=xminsig
    xmax=xmaxsig

    pylab.subplots_adjust(left=0.125, right=.9,bottom=.1,top=0.9,wspace=.3,hspace=.1)
    p1=pylab.subplot(311)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y05)
    (fxbin,fybin,fybinerr)=binitbins(sigmamin,sigmamax,nbin,fx,fy05)
    ybin=ybin/(xbin/1000.)**3
    ybinerr=ybinerr/(xbin/1000.)**3
    fybin=fybin/(fxbin/1000.)**3
    fybinerr=fybinerr/(fxbin/1000.)**3
    scale=ybin[nbinscale]/fybin[nbinscale]
    fybin=fybin*scale
    fybinnerr=fybinerr*scale
    clabel="Control*%3.1f"%(scale)
    pylab.plot(xbin,ybin,'ko',label="Cluster",markersize=6.)
    pylab.plot(fxbin,fybin,'wo',label=clabel,markersize=6.)
    pylab.hold(True)
    t=pylab.legend(loc='lower left')
    pylab.errorbar(fxbin,fybin,fybinerr,fmt=None,ecolor='k')
    pylab.errorbar(xbin,ybin,ybinerr,fmt='ko')
    c1=ybin[nbinscale]#/(xbin[nbinscale]/1000.)**3
    #yl=c1*(xl/1000.)**3 
    yl=c1*N.ones(len(xl),'f')
    pylab.plot(xl,yl,'k-')

    ymin=.5*min(ybin)
    ymax=1.5*max(ybin)
    pylab.axis([400.,1300.,ymin,ymax])
    ax=pylab.gca()
    ax.set_yscale("log")
    ax.set_xscale("log")
    pylab.setp(ax,xticklabels=[])

    pylab.text(.85,.85,r'$0.5\times R_{200}$',horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=20)
    pylab.subplot(312)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y1)
    (fxbin,fybin,fybinerr)=binitbins(sigmamin,sigmamax,nbin,fx,fy1)
    ybin=ybin/(xbin/1000.)**3.
    ybinerr=ybinerr/(xbin/1000.)**3.
    fybin=fybin/(fxbin/1000.)**3.
    fybinerr=fybinerr/(fxbin/1000.)**3.
    scale=ybin[nbinscale]/fybin[nbinscale]
    fybin=fybin*scale
    fybinnerr=fybinerr*scale
    clabel="Control*%3.1f"%(scale)
    pylab.plot(xbin,ybin,'ko',label="Cluster",markersize=6.)
    pylab.plot(fxbin,fybin,'wo',label=clabel,markersize=6.)
    pylab.hold(True)
    pylab.legend(loc='lower left')
    pylab.errorbar(fxbin,fybin,fybinerr,fmt=None,ecolor='k')
    pylab.errorbar(xbin,ybin,ybinerr,fmt='ko')
    c1=ybin[nbinscale]#/(xbin[nbinscale]/1000.)**3
    #yl=c1*(xl/1000.)**3 
    yl=c1*N.ones(len(xl),'f')
    pylab.plot(xl,yl,'k-')

    ymin=.5*min(ybin)
    ymax=1.5*max(ybin)
    pylab.axis([400.,1300.,ymin,ymax])

    ax=pylab.gca()
    ax.set_yscale("log")
    ax.set_xscale("log")
    pylab.setp(ax,xticklabels=[])
    pylab.text(.85,.85,r'$1.0\times R_{200}$',horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=20)

    pylab.subplot(313)

    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y2)
    (fxbin,fybin,fybinerr)=binitbins(sigmamin,sigmamax,nbin,fx,fy2)
    ybin=ybin/(xbin/1000.)**3
    ybinerr=ybinerr/(xbin/1000.)**3
    fybin=fybin/(fxbin/1000.)**3
    fybinerr=fybinerr/(fxbin/1000.)**3
    scale=ybin[nbinscale]/fybin[nbinscale]
    fybin=fybin*scale
    fybinnerr=fybinerr*scale
    clabel="Control*%3.1f"%(scale)
    pylab.plot(xbin,ybin,'ko',label="Cluster",markersize=6.)
    pylab.plot(fxbin,fybin,'wo',label=clabel,markersize=6.)
    pylab.hold(True)
    pylab.legend(loc='lower left')
    pylab.errorbar(fxbin,fybin,fybinerr,fmt=None,ecolor='k')
    pylab.errorbar(xbin,ybin,ybinerr,fmt='ko')
    c1=ybin[nbinscale]#/(xbin[nbinscale]/1000.)**3
    yl=c1*N.ones(len(xl),'f')
    pylab.plot(xl,yl,'k-')

    ymin=.5*min(ybin)
    ymax=1.5*max(ybin)
    pylab.axis([400.,1300.,ymin,ymax])

    ax=pylab.gca()
    ax.set_yscale("log")
    ax.set_xscale("log",)
    pylab.text(.85,.85,r'$2.0\times R_{200}$',horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=20)
    pylab.xlabel(r'$\sigma \ \rm{(km/s)}$',fontsize=24)
    pylab.setp(ax,xticklabels=['1000.'])

def plotstellarsigma3sigmapylab():
    x=(c.sigma)
    y05=c.stellmass05/(10.**11)
    y1=c.stellmass1/(10.**11)
    y2=c.stellmass2/(10.**11)
    fx=(f.sigma)
    fy05=f.stellmass05/(10.**11)
    fy1=f.stellmass1/(10.**11)
    fy2=f.stellmass2/(10.**11)

    pylabsubplot2(x,y05,y1,y2,fx,fy05,fy1,fy2)
    pylab.subplot(312)
    pylab.ylabel(r'$\rm{\Sigma \ M_* \  (10^{11} M_\odot)}$',fontsize=24)
    pylab.savefig('stellmassmembsigma3sigma.eps')

def plotsfrsigma3sigmapylab():

    x=(c.sigma)
    y05=c.sfr05#sumitmemb(3.,0.5,0)
    y1=c.sfr1#sumitmemb(3.,1.,0)
    y2=c.sfr2#sumitmemb(3.,2.,0)
    fx=(f.sigma)
    fy05=f.sfr05#sumitmemb(3.,0.5,0)
    fy1=f.sfr1#sumitmemb(3.,1.,0)
    fy2=f.sfr2#sumitmemb(3.,2.,0)
    pylabsubplot2(x,y05,y1,y2,fx,fy05,fy1,fy2)

    pylab.subplot(312)
    pylab.ylabel(r'$\rm{\Sigma SFR \  (M_\odot/yr)}$',fontsize=24)
    pylab.savefig('sfrmembsigma3sigma.eps')

def plotngalsigma3sigmapylab():
    x=(c.sigma)
    y05=c.ngal05#sumitmemb(3.,0.5,0)
    y1=c.ngal1#sumitmemb(3.,1.,0)
    y2=c.ngal2#sumitmemb(3.,2.,0)
    fx=(f.sigma)
    fy05=f.ngal05#sumitmemb(3.,0.5,0)
    fy1=f.ngal1#sumitmemb(3.,1.,0)
    fy2=f.ngal2#sumitmemb(3.,2.,0)
    pylabsubplot2(x,y05,y1,y2,fx,fy05,fy1,fy2)

    pylab.subplot(312)
    pylab.ylabel(r'$\rm{N_{gal}}$',fontsize=24,fontweight='bold')
    pylab.savefig('ngalmembsigma3sigma.eps')



def plotsfr(lbin):
    xminsig=2.6
    xmaxsig=3.1
    symheight=3.
    sigmin=400.
    sigmax=1200.
    nbin=7

    ylabel="\gSSFR (M\d\(2281)\u yr\u-1\d)"
    #ymin=0.
    #ymax=105.
    #xmin=400.
    #xmax=1200.

    #log space
    ymin=0.
    ymax=2.5
    xmin=xminsig
    xmax=xmaxsig
    #ymin=-2.
    #ymax=4.
    #xmin=0.
    #xmax=5.


    psplotinit("sfrmembsigma3.ps")
    x=(c.sigma)
    y05=c.sfr05#sumitmemb(3.,0.5,0)
    y1=c.sfr1#sumitmemb(3.,1.,0)
    y2=c.sfr2#sumitmemb(3.,2.,0)
    y3=c.sfr3#sumitmemb(3.,3.,0)
    fx=(f.sigma)
    fy05=f.sfr05*10.#sumitmemb(3.,0.5,0)
    fy1=f.sfr1*10.#sumitmemb(3.,1.,0)
    fy2=f.sfr2*10.#sumitmemb(3.,2.,0)
    fy3=f.sfr3*10.#sumitmemb(3.,3.,0)

    ppgplot.pgbox("",0.0,0,"",0.0,1)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,30)
    ppgplot.pglab("\gs (km/s)",ylabel,"")
    ppgplot.pgsci(1)
    
    print "SFR VS SIGMA"
    print "0.5 Rv"
    dospear(x,y05)
    #print "Regression in log-log space"
    #dolinregress(N.log10(x),N.log10(y05))
    
    print "1.0 Rv"
    dospear(x,y1)
    #print "Regression in log-log space"
    #xt=N.compress(y1>0,x)
    #yt=N.compress(y1>0,y1)
    #print "cut out zeros before running linear regression on field sample"
    #dolinregress(N.log10(xt),N.log10(yt))

    #print "Regression in log-log space"
    #dolinregress(N.log10(x),N.log10(y1))
    print "2.0 Rv"
    dospear(x,y2)
    
    
    #(xbin,ybin,ybinerr)=biniterr(x,y05,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y05)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgsch(symheight)
    ppgplot.pgpt(xbin,ybin,18)

    #(xbin,ybin,ybinerr)=biniterr(fx,fy05,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,fx,fy05)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,12)

    print "1.0 Rv, binned data"
    #(xbin,ybin,ybinerr)=biniterr(x,y1,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y1)
    xbin=N.log10(xbin)
    #ybin=ybin+1
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,16)
    dolinregress(xbin,ybin)

    print "Regression in log-log space using matplotlib."
    print xbin,ybin
    (a,b)=pylab.polyfit(xbin,ybin,1)
    matout = "slope = %5.2f, intercept = %5.2f " %(a,b)
    print matout

    print "1.0 Rv, Control binned data"
    #(xbin,ybin,ybinerr)=biniterr(fx,fy1,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,fx,fy1)
    xbin=N.log10(xbin)
    #ybin=ybin+1
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,6)
    dolinregress(xbin,ybin)

    print "Regression in log-log space using matplotlib."
    print xbin,ybin
    (a,b)=pylab.polyfit(xbin,ybin,1)
    matout = "slope = %5.2f, intercept = %5.2f " %(a,b)
    print matout


    #(xbin,ybin,ybinerr)=biniterr(x,y2,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y2)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,17)

    #(xbin,ybin,ybinerr)=biniterr(fx,fy2,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,fx,fy2)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,22)


    #draw lines showing sfr proportional to mass ~ sigma^3

    ppgplot.pgsls(4)
    ppgplot.pgslw(mylineswidth)
    x=N.arange(0.,5.,1.)
    lines=[-6.81,-7.25,-7.68]
    if cdir.find('field') > -1:
        lines=[-6.7,-7.2,-7.8]
    lines=mylines
    for y0 in lines:  
        y=3*x +y0 
        ppgplot.pgline(x,y)

    ppgplot.pgsls(1)
    ppgplot.pgend()


    #comparison with high-z clusters
    ylabel="\gSSFR (M\d\(2281)\u yr\u-1\d)"

    #log space
    ymin=-0.
    ymax=3.3
    #xmin=2.5
    #xmax=3.2
    #ymin=-2.
    #ymax=4.
    #xmin=0.
    #xmax=5.

    
    psplotinit("sfrmembsigma3hiz.ps")
    x=(c.sigma)
    y05=c.sfr05hiz#sumitmemb(3.,0.5,0)
    ppgplot.pgbox("",0.0,0,"",0.0,1)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,30)
    ppgplot.pglab("\gs (km/s)",ylabel,"")
    ppgplot.pgsci(1)

    #(xbin,ybin,ybinerr)=biniterr(x,y05,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y05)

    #print "sigmbin = ",c.sigmabin
    #print "xbin = ",xbin
    xbin=N.log10(xbin)
    #ybin=ybin+1
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgsch(symheight)
    ppgplot.pgpt(xbin,ybin,18)

    #plot hiz clusters
    ppgplot.pgsci(2)
    x=N.log10(hiz.sigma)
    errorylog(x,hiz.sumsfr,hiz.sumsfrerr)
    y=N.log10(hiz.sumsfr)
    pgpnts(x,y,hiz.symbols)
    #for i in range(len(x)):
    #    print "hey baby",i, hiz.sigma[i],hiz.sumsfr[i],hiz.sumsfrerr[i]
    xlabel = 2.57
    ylabel = 3.1
    ystep = .11
    dy=.05
    dx=.02
    ppgplot.pgsch(1.7)
    ppgplot.pgsci(1)
    hizlabel(xlabel,ylabel,dx,dy,ystep,hiz.ref,hiz.refsymb,hiz.colors)

    #draw lines showing sfr proportional to mass ~ sigma^3
    ppgplot.pgsci(1)
    ppgplot.pgsls(4)
    x=N.arange(0.,5.,1.)
    y0=-6.95
    y=3*x +y0 
    ppgplot.pgline(x,y)

    y0=-7.4
    y=3*x +y0 

    #y=(x/400)**3 + y0
    ppgplot.pgline(x,y)

    y0=-7.9
    y=3*x +y0 
    #y=(x/400)**3 + y0
    ppgplot.pgline(x,y)
    ppgplot.pgsls(1)
    ppgplot.pgend()


    ylabel="\gSSFR/(\gs/1000)\u3\d (M\d\(2281)\u yr\u-1\d/(km/s)\u3\d)"
    #ymin=0.
    #ymax=105.
    #xmin=400.
    #xmax=1200.

    #log space
    ymin=-0.2
    ymax=3.
    #xmin=2.5
    #xmax=3.2
    #ymin=-2.
    #ymax=4.
    #xmin=0.
    #xmax=5.


    psplotinit("sfrsigma3sigma.ps")#SFR/sigma^3 vs sigma
    x=(c.sigma)
    y05=c.sfr05
    y1=c.sfr1
    y2=c.sfr2

    fx=(f.sigma)
    fy05=f.sfr05*10.
    fy1=f.sfr1*10.
    fy2=f.sfr2*10.

    ppgplot.pgbox("",0.0,0,"",0.0,1)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,30)
    ppgplot.pglab("\gs (km/s)",ylabel,"")
    ppgplot.pgsci(1)


    print "SFR/SIGMA^3 VS SIGMA"
    print "0.5 Rv"
    #(xbin,ybin,ybinerr)=biniterr(x,y05,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y05)
    ybin=ybin/((xbin/1000.)**3)
    ybinerr=ybinerr/((xbin/1000.)**3)
    xs=x/1000.
    ys=y05/((x/1000)**3)
    dospear(xs,ys)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgsch(symheight)
    ppgplot.pgpt(xbin,ybin,18)

    ppgplot.pgsls(4)
    ppgplot.pgslw(mylineswidth)
    xl=N.arange(0.,5.,1.)
    yl0=N.average(ybin)
    yl=0*xl +yl0 
    ppgplot.pgline(xl,yl)
    ppgplot.pgsls(1)

    print "0.5 Rv - Control data"
    #(xbin,ybin,ybinerr)=biniterr(fx,fy05,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,fx,fy05)
    ybin=ybin/((xbin/1000.)**3)
    ybinerr=ybinerr/((xbin/1000.)**3)
    xs=x/1000.
    ys=y05/((x/1000)**3)
    dospear(xs,ys)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,12)

    ppgplot.pgsls(4)
    ppgplot.pgslw(mylineswidth)
    xl=N.arange(0.,5.,1.)
    yl0=N.average(ybin)
    yl=0*xl +yl0 
    ppgplot.pgline(xl,yl)
    ppgplot.pgsls(1)

    #(xbin,ybin,ybinerr)=biniterr(x,y1,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y1)
    ybin=ybin/((xbin/1000.)**3)
    ybinerr=ybinerr/((xbin/1000.)**3)
    print "1.0 Rv"
    xs=x/1000.
    ys=y1/((x/1000)**3)
    dospear(xs,ys)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,16)
    ppgplot.pgsls(4)
    xl=N.arange(0.,5.,1.)
    yl0=N.average(ybin)
    yl=0*xl + yl0 
    ppgplot.pgline(xl,yl)
    ppgplot.pgsls(1)

    #(xbin,ybin,ybinerr)=biniterr(fx,fy1,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,fx,fy1)
    ybin=ybin/((xbin/1000.)**3)
    ybinerr=ybinerr/((xbin/1000.)**3)
    print "1.0 Rv - CONTROL DATA"
    xs=x/1000.
    ys=y1/((x/1000)**3)
    dospear(xs,ys)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,6)
    ppgplot.pgsls(4)
    xl=N.arange(0.,5.,1.)
    yl0=N.average(ybin)
    yl=0*xl + yl0 
    ppgplot.pgline(xl,yl)
    ppgplot.pgsls(1)

    #(xbin,ybin,ybinerr)=biniterr(x,y2,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y2)
    ybin=ybin/((xbin/1000.)**3)
    ybinerr=ybinerr/((xbin/1000.)**3)
    print "2.0 Rv"
    xs=x/1000.
    ys=y2/((x/1000)**3)
    dospear(xs,ys)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,17)

    ppgplot.pgsls(4)
    ppgplot.pgslw(mylineswidth)
    xl=N.arange(0.,5.,1.)
    yl0=N.average(ybin)
    yl=0*xl + yl0 
    ppgplot.pgline(xl,yl)
    ppgplot.pgsls(1)

    #(xbin,ybin,ybinerr)=biniterr(fx,fy2,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,fx,fy2)
    ybin=ybin/((xbin/1000.)**3)
    ybinerr=ybinerr/((xbin/1000.)**3)
    print "2.0 Rv - CONTROL DATA"
    xs=x/1000.
    ys=y2/((x/1000)**3)
    dospear(xs,ys)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,22)

    ppgplot.pgsls(4)
    ppgplot.pgslw(mylineswidth)
    xl=N.arange(0.,5.,1.)
    yl0=N.average(ybin)
    yl=0*xl + yl0 
    ppgplot.pgline(xl,yl)
    ppgplot.pgsls(1)

    ppgplot.pgsls(1)
    ppgplot.pgend()




    #comparison with hiz clusters
    ylabel="\gSSFR/(\gs/1000)\u3\d (M\d\(2281)\u yr\u-1\d/(km/s)\u3\d)"
    #log space
    ymin=.9
    ymax=3.3
    #xmin=2.5
    #xmax=3.2
    #ymin=-2.
    #ymax=4.
    #xmin=0.
    #xmax=5.
    psplotinit("sfrsigma3sigmahiz.ps")#SFR/sigma^3 vs sigma
    x=(c.sigma)
    y05=c.sfr05hiz#sumitmemb(3.,0.5,0)
    ppgplot.pgbox("",0.0,0,"",0.0,1)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,30)
    ppgplot.pglab("\gs (km/s)",ylabel,"")
    ppgplot.pgsci(1)

    #nbin=5
    #(xbin,ybin,ybinerr)=biniterr(x,y05,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y05)
    #ybin=ybin+1
    ybin=ybin/((xbin/1000.)**3)
    ybinerr=ybinerr/((xbin/1000.)**3)
    for i in range(len(xbin)):
        print "SFR/sigma^3, C4",i, xbin[i],ybin[i],ybinerr[i]

    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgsch(symheight)
    ppgplot.pgpt(xbin,ybin,18)

    #plot hiz clusters
    ppgplot.pgsci(2)
    x=N.log10(hiz.sigma)
    y=hiz.sumsfr/((hiz.sigma/1000.)**3*hiz.cosm)
    erry=hiz.sumsfrerr/((hiz.sigma/1000.)**3*hiz.cosm)
    errorylog(x,y,erry)
    y=N.log10(y)
    pgpnts(x,y,hiz.symbols)

    xlabel = 2.95
    ylabel = 3.15
    ystep = .09
    dy=.03
    dx=.02
    ppgplot.pgsci(1)
    hizlabel(xlabel,ylabel,dx,dy,ystep,hiz.ref,hiz.refsymb,hiz.colors)


    #draw lines showing sfr proportional to mass ~ sigma^3

    #draw lines showing sfr proportional to mass ~ sigma^3
    ppgplot.pgsls(4)
    x=N.arange(0.,5.,1.)
    y0=2.07
    y=0*x +y0 
    ppgplot.pgline(x,y)

    y0=1.78
    y=0*x +y0 

    #y=(x/400)**3 + y0
    ppgplot.pgline(x,y)

    y0=1.4
    y=0*x +y0 
    #y=(x/400)**3 + y0
    ppgplot.pgline(x,y)

    ppgplot.pgsls(1)
    ppgplot.pgend()


    #comparison with hiz clusters
    ylabel="\gSSFR/(\gs/1000)\u3\d (M\d\(2281)\u yr\u-1\d/(km/s)\u3\d)"
    #log space
    ymin=.9
    ymax=3.3
    xmin=-.01
    xmax=1
    #ymin=-2.
    #ymax=4.
    #xmin=0.
    #xmax=5.
    psplotinit("sfrsigma3zhiz.ps")#SFR/sigma^3 vs sigma

    y05=c.sfr05hiz#sumitmemb(3.,0.5,0)
    ppgplot.pgbox("",0.0,0,"",0.0,1)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,20)
    ppgplot.pglab("z",ylabel,"")
    ppgplot.pgsci(1)

    nbin=1
    x=c.sigma
    x1=c.z
    (xbin,ybin,ybinerr)=biniterr(x,y05,nbin)
    (xbin,zbin,zbinerr)=biniterr(x,x1,nbin)
    #ybin=ybin+1
    ybin=ybin/((xbin/1000.)**3)
    ybinerr=ybinerr/((xbin/1000.)**3)
    print "C4 clusters: z, SFR/sigma^3"
    print zbin,ybin

    #xbin=N.log10(xbin)
    errorylog(zbin,ybin,ybinerr)
    errorx(zbin,ybin,zbinerr)
    ybin=N.log10(ybin)
    ppgplot.pgsch(symheight)
    ppgplot.pgpt(zbin,ybin,18)

    #plot hiz clusters
    ppgplot.pgsci(2)
    #x=N.log10(hiz.sigma)
    x=(hiz.z)
    y=hiz.sumsfr/((hiz.sigma/1000.)**3*hiz.cosm)
    print "High z clusters: z, SFR/sigma^3"
    for i in range(len(x)):
        print i, x[i],y[i]
    erry=hiz.sumsfrerr/((hiz.sigma/1000.)**3*hiz.cosm)
    xave=(N.average(x[0:4]))
    yave=(N.average(y[0:4]))
    yerr=pylab.std((y[0:4]))
    xerr=pylab.std((x[0:4]))

    errorylog(x,y,erry)
    y=N.log10(y)
    pgpnts(x,y,hiz.symbols)

    xlabel = .6
    ylabel = 1.8
    ystep = .09
    dy=.03
    dx=.02
    ppgplot.pgsci(1)
    hizlabel(xlabel,ylabel,dx,dy,ystep,hiz.ref,hiz.refsymb,hiz.colors)
    #plot average for highz
    xave=N.array([xave],'f')
    yave=N.array([yave],'f')
    xerr=N.array([xerr],'f')
    yerr=N.array([yerr],'f')
    ppgplot.pgsci(1)
    errorylog(xave,yave,yerr)
    yave=N.log10(yave)
    errorx(xave,yave,xerr)
    ppgplot.pgsch(5)
    ppgplot.pgpt(xave,yave,17)
    ppgplot.pgsch(1.7)
    ppgplot.pgsci(1)

    ppgplot.pgsls(4)
    x=N.arange(0.,5.,.1)

    y0=19.
    y=7.*(1+x)+y0 
    #y=N.log10(y)
    #ppgplot.pgline(x,y)

    y=(1+x)**7*y0
    y=N.log10(y)
    ppgplot.pgline(x,y)

    y=(1+x)**5*y0
    y=N.log10(y)
    ppgplot.pgline(x,y)

    y=(1+x)**6*y0
    y=N.log10(y)
    ppgplot.pgline(x,y)



    ppgplot.pgend()

    normal=11.
    ymax=2.7
    ymin=0.


    ylabel="\gS M\d*\u (M\d\(2281)\u yr\u-1\d / 10\u"+str(int(normal))+"\d)"
    xmin=xminsig
    xmax=xmaxsig
    
    psplotinit("stellarmembsigma3.ps")
    x=(c.sigma)
    y05=c.stellmass05#sumstellmass(3.,0.5)
    y1=c.stellmass1
    y2=c.stellmass2
    fx=(f.sigma)
    fy05=f.stellmass05*10.
    fy1=f.stellmass1*10.
    fy2=f.stellmass2*10.
    ppgplot.pgbox("",0.0,0,"",0.0,1)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,30)
    ppgplot.pglab("\gs (km/s)",ylabel,"")
    ppgplot.pgsci(1)
    
    print "STELLAR MASS VS SIGMA"
    print "0.5 Rv"
    dospear(x,y05)
    print "0.5 Rv - CONTROL DATA"
    dospear(fx,fy05)
    print "1.0 Rv, unbinned"
    dospear(x,y1)
    print "1.0 Rv, unbinned - CONTROL DATA"
    dospear(fx,fy1)
            

    print "2.0 Rv"
    dospear(x,y2)
    print "2.0 Rv - CONTROL DATA"
    dospear(fx,fy2)
    
    nbin=7
    #(xbin,ybin,ybinerr)=biniterr(x,y05,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y05)
    ybin=(ybin)/(10.**normal)
    ybinerr=(ybinerr)/(10.**normal)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,18)

    #(xbin,ybin,ybinerr)=biniterr(fx,fy05,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,fx,fy05)
    ybin=(ybin)/(10.**normal)
    ybinerr=(ybinerr)/(10.**normal)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,12)


    #(xbin,ybin,ybinerr)=biniterr(x,y1,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y1)
    ybin=(ybin)/(10.**normal)
    ybinerr=(ybinerr)/(10.**normal)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,16)
    print "1.0 Rv, binned data"
    dolinregress(xbin,ybin)

    #(xbin,ybin,ybinerr)=biniterr(fx,fy1,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,fx,fy1)
    ybin=(ybin)/(10.**normal)
    ybinerr=(ybinerr)/(10.**normal)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,6)
    print "1.0 Rv, binned data - CONTROL DATA"
    dolinregress(xbin,ybin)


    #(xbin,ybin,ybinerr)=biniterr(x,y2,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y2)
    ybin=(ybin)/(10.**normal)
    ybinerr=(ybinerr)/(10.**normal)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,17)

    #(xbin,ybin,ybinerr)=biniterr(fx,fy2,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,fx,fy2)
    ybin=(ybin)/(10.**normal)
    ybinerr=(ybinerr)/(10.**normal)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,22)

    #draw lines showing stellar mass proportional to mass ~ sigma^3
    ppgplot.pgsls(4)
    ppgplot.pgslw(mylineswidth)
    x=N.arange(0.,5.,1.)
    lines=[4.14,3.87,3.6]
    if cdir.find('field') > -1:
        lines=[4.15,3.59,3.05]
    lines=mylines
    for y0 in lines:  
        y=3*x +y0 -11.
        ppgplot.pgline(x,y)
        #ppgplot.pgsci(2)
        #y=2.77*x +y0 -11.+.62
        #ppgplot.pgline(x,y)
        #ppgplot.pgsci(1)

    ppgplot.pgsls(1)
    ppgplot.pgend()


    normal=11.
    #ymin=11.5-normal
    #ymax=13.-normal
    ymax=2.7
    ymin=-.15


    ylabel="N\dgal\u "
    xmin=xminsig
    xmax=xmaxsig
    ymin=0.5
    ymax=3.
    psplotinit("ngalmembsigma3.ps")
    x=(c.sigma)
    y05=c.ngal05#ngalmemb
    y1=c.ngal1
    y2=c.ngal2
    fx=(f.sigma)
    fy05=f.ngal05*10.
    fy1=f.ngal1*10.
    fy2=f.ngal2*10.
    ppgplot.pgbox("",0.0,0,"",0.0,1)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,30)
    ppgplot.pglab("\gs (km/s)",ylabel,"")
    ppgplot.pgsci(1)
    
    
    nbin=7
    #(xbin,ybin,ybinerr)=biniterr(x,y05,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y05)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgsch(symheight)
    ppgplot.pgpt(xbin,ybin,18)

    #(xbin,ybin,ybinerr)=biniterr(fx,fy05,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,fx,fy05)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,12)

    #(xbin,ybin,ybinerr)=biniterr(x,y1,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y1)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,16)

    #(xbin,ybin,ybinerr)=biniterr(fx,fy1,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,fx,fy1)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,6)

    #(xbin,ybin,ybinerr)=biniterr(x,y2,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y2)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,17)

    #(xbin,ybin,ybinerr)=biniterr(fx,fy2,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,fx,fy2)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,22)

    #draw lines showing stellar mass proportional to ngal ~ sigma^3
    ppgplot.pgsls(4)
    ppgplot.pgslw(mylineswidth)
    x=N.arange(0.,5.,1.)
    #lines=[4.4,4.14,3.87,3.6]
    lines=mylines
    #if cdir.find('field') > -1:
    #    lines=[4.15,3.59,3.05]
    for y0 in lines:  
        y=3*x +y0 -11.
        ppgplot.pgline(x,y)
    ppgplot.pgsls(1)
    ppgplot.pgend()



    #ymin=0.
    #ymax=105.
    #xmin=400.
    #xmax=1200.

    #log space
    #ymin=0.
    #ymax=2.2
    #xmin=2.5
    #xmax=3.2
    normal=11.
    ymin=11.5-normal
    ymax=13.9-normal


    ylabel="\gS M\d*\u/(\gs/1000)\u3\d (M\d\(2281)\u /10\u"+str(int(normal))+"\d/(km/s)\u3\d)"
    xmin=xminsig
    xmax=xmaxsig
    
    
    psplotinit("stellarsigma3sigma.ps")
    x=(c.sigma)
    y05=c.stellmass05#sumstellmass(3.,0.5)
    y1=c.stellmass1#sumstellmass(3.,1.)
    y2=c.stellmass2#sumstellmass(3.,2.)

    fx=(f.sigma)
    fy05=f.stellmass05*10.#sumstellmass(3.,0.5)
    fy1=f.stellmass1*10.#sumstellmass(3.,1.)
    fy2=f.stellmass2*10.#sumstellmass(3.,2.)
    ppgplot.pgbox("",0.0,0,"",0.0,1)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,30)
    ppgplot.pglab("\gs (km/s)",ylabel,"")
    ppgplot.pgsci(1)

    #(xbin,ybin,ybinerr)=biniterr(x,y05,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y05)
    ybin=(ybin)/(10.**normal)
    ybinerr=(ybinerr)/(10.**normal)
    ybin=ybin/((xbin/1000.)**3)
    ybinerr=ybinerr/((xbin/1000.)**3)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgsch(symheight)
    ppgplot.pgpt(xbin,ybin,18)
    ppgplot.pgsls(4)
    ppgplot.pgslw(mylineswidth)
    xl=N.arange(0.,5.,1.)
    yl0=N.average(ybin)
    yl=0*xl + yl0 
    ppgplot.pgline(xl,yl)
    ppgplot.pgsls(1)

    #(xbin,ybin,ybinerr)=biniterr(fx,fy05,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,fx,fy05)
    ybin=(ybin)/(10.**normal)
    ybinerr=(ybinerr)/(10.**normal)
    ybin=ybin/((xbin/1000.)**3)
    ybinerr=ybinerr/((xbin/1000.)**3)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,12)
    ppgplot.pgsls(4)
    ppgplot.pgslw(mylineswidth)
    xl=N.arange(0.,5.,1.)
    yl0=N.average(ybin)
    yl=0*xl + yl0 
    ppgplot.pgline(xl,yl)
    ppgplot.pgsls(1)

    #(xbin,ybin,ybinerr)=biniterr(x,y1,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y1)
    ybin=(ybin)/(10.**normal)
    ybinerr=(ybinerr)/(10.**normal)
    ybin=ybin/((xbin/1000.)**3)
    ybinerr=ybinerr/((xbin/1000.)**3)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,16)
    ppgplot.pgsls(4)
    ppgplot.pgslw(mylineswidth)
    xl=N.arange(0.,5.,1.)
    yl0=N.average(ybin)
    yl=0*xl + yl0 
    ppgplot.pgline(xl,yl)
    ppgplot.pgsls(1)

    #(xbin,ybin,ybinerr)=biniterr(fx,fy1,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,fx,fy1)
    ybin=(ybin)/(10.**normal)
    ybinerr=(ybinerr)/(10.**normal)
    ybin=ybin/((xbin/1000.)**3)
    ybinerr=ybinerr/((xbin/1000.)**3)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,6)
    ppgplot.pgsls(4)
    ppgplot.pgslw(mylineswidth)
    xl=N.arange(0.,5.,1.)
    yl0=N.average(ybin)
    yl=0*xl + yl0 
    ppgplot.pgline(xl,yl)
    ppgplot.pgsls(1)

    #(xbin,ybin,ybinerr)=biniterr(x,y2,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,x,y2)
    ybin=(ybin)/(10.**normal)
    ybinerr=(ybinerr)/(10.**normal)
    ybin=ybin/((xbin/1000.)**3)
    ybinerr=ybinerr/((xbin/1000.)**3)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,17)
    ppgplot.pgsls(4)
    ppgplot.pgslw(mylineswidth)
    xl=N.arange(0.,5.,1.)
    yl0=N.average(ybin)
    yl=0*xl +yl0 
    ppgplot.pgline(xl,yl)
    ppgplot.pgsls(1)

    #(xbin,ybin,ybinerr)=biniterr(fx,fy2,nbin)
    (xbin,ybin,ybinerr)=binitbins(sigmamin,sigmamax,nbin,fx,fy2)
    ybin=(ybin)/(10.**normal)
    ybinerr=(ybinerr)/(10.**normal)
    ybin=ybin/((xbin/1000.)**3)
    ybinerr=ybinerr/((xbin/1000.)**3)
    xbin=N.log10(xbin)
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,22)
    ppgplot.pgsls(4)
    ppgplot.pgslw(mylineswidth)
    xl=N.arange(0.,5.,1.)
    yl0=N.average(ybin)
    yl=0*xl +yl0 
    ppgplot.pgline(xl,yl)
    ppgplot.pgsls(1)

    ppgplot.pgsls(1)
    ppgplot.pgend()


def getlitsfr(zhiz,sigmahiz):
    mabscut=-15.
    #for reading ediscs files
    #output.write("%8.1f %8.1f %6.2f %6.1f %6.1f %6.2f %6.2f  \n" % (dra,ddec,Mr[i],ew[i],errew[i],sfr[i],errsfr[i]))
    sumsfrhiz=N.ones(len(zhiz),'f')
    errsumsfrhiz=N.ones(len(zhiz),'f')
    #cl1040,cl1054,cl1216,clj0023,A2390,AC114,A1689,cl0024
    #                             BM00  C01   B02   K04
    infiles=['/Users/rfinn/clusters/final/cl1040sdss.dat','/Users/rfinn/clusters/final/cl1054sdss.dat','/Users/rfinn/clusters/final/cl1216sdss.dat','/Users/rfinn/clusters/final/clj0023sdss.dat','/Users/rfinn/clusters/final/lit_data/balogh00table2.dat','/Users/rfinn/clusters/final/lit_data/couchdatafile1.txt','/Users/rfinn/clusters/final/lit_data/balogh02table1.dat','/Users/rfinn/SDSS/cl0024/cl0024_photzmem_finn.dat']
    for i in range(len(sumsfrhiz)):
        allsfr=[]
        sumsfr=0.
        errsumsfr=0.
        z=zhiz[i]
        sigma=sigmahiz[i]
        r2=r200(sigma,z,h100)#r200 in Mpc
        angd=DA(z,h100)#angular diameter distance, kpc/arcsec
        r2arcsec=r2*1000./angd#r200 in arcsec
        r2arcsec=0.5*r2arcsec#limit to 0.5R200
        DL=dL(z,h100)
        DLcm=DL*1.e6*3.09e18
        distmod=5*log10(DL*1.e6/10)
        #print "dL Mpc, dL cm, distmod = ",DL,DLcm,distmod
        if i == 3:#CLJ0023, kcorr J_Rc, band=8 in mystuff.py
            k=kcorr(z,8)
            kc=N.average(k[1:3])#take average of Sbc/Scd gal types
        if i == 4:#A2390, kcorr R_Rc, band=4 in mystuff.py
            k=kcorr(z,4)
            kc=N.average(k[1:3])

        if i == 5:#Couch et al
            dv=6.*sigma/3.e5#apply +/- 3sigma cut in velocity
            dz=dv*(1+z)
            zmin=z-dz
            zmax=z+dz
            mImax=3.21+5*log10(DL)#see notebook for calculation of limiting I mag
            cra=15*(22+58./60.+52.3/3600.)#RA center from NED,convert to deg
            cdec=-1*(34+46/60. +55./3600.)#DEC center from NED
            k=kcorr(z,6)#kcorr I_Rc, band = 6 in mystuff.py
            kc=N.average(k[1:3])

        if i == 6:
            dv=6.*sigma/3.e5
            dz=dv*(1+z)
            zmin=z-dz
            zmax=z+dz
            cra=15*(13+11./60.+34.2/3600.)#RA center from NED,convert to deg
            cdec=-1*(1+21/60. +56./3600.)#DEC center from NED    
            mImax=3.21+5*log10(DL)#see notebook for calculation of limiting I mag
            k=kcorr(z,6)#kcorr I_Rc, band = 6 in mystuff.py
            kc=N.average(k[1:3])
        if i == 7:#CL0024, kcorr R_Rc, band=4 in mystuff.py
            k=kcorr(z,4)
            kc=N.average(k[1:3])


        #dra=[]
        #ddec=[]
        #Mr=[]
        #ew=[]
        #sfr=[]
        infile=open(infiles[i],'r')
        for line in infile:
            if line.find('#') > -1:
                continue
            f=line.split()
            for j in range(len(f)):
                try:
                    f[j]=float(f[j])
                except:
                    f[j]=f[j]
            if i < 3:#f=line.split()#(dra,ddec,Mr[i],ew[i],errew[i],sfr[i],errsfr[i]))
                x=f[0]
                y=f[1]
                m=f[2]
                sfr=f[5]
                errsfr=f[6]
            if i == 3:#f=line.split()#(dra,ddec,mj[i],ew[i],errew[i],sfr[i],errsfr[i]))
                x=f[0]
                y=f[1]
                m=f[2] -distmod -kc
                sfr=f[5]
                errsfr=f[6]
            if i == 4:#Balogh & Morris
                x=f[5]
                y=f[6]
                m=f[9] - distmod -kc #+ 0.23#observed R
                sfr=0
                #memb=str(f[17]) - only means it has CNOC spectroscopy
                #if memb.find('y') > -1:
                SEflag=float(f[14])
                if (SEflag < 1):
                    if (f[10]/f[11]) > 3:#require 3 sigma
                    #if (float(f[12]) > 50):#require 3 sigma
                        const=0.152#see notebook re conversion, from flux to sfr
                        sfr=f[10]*.001*DL**2*const
                        errsfr=f[11]*.001*DL**2*const
                        #flux in Jy*10^-23erg/s/cm^2/Hz*dv
                        #dv=1.591e13, dla1m=348A, lam=8100A
                        sfr=f[10]*1.591e-10*4.*3.1415*(DLcm**2)*7.9e-42
                        errsfr=f[11]*1.591e-10*4.*3.1415*(DLcm**2)*7.9e-42
                        sfr=0.001*sfr#online cat flux units are mJy, not Jy
                        errsfr=0.001*errsfr
                        sfr=sfr*2.5#correct for dust
                        errsfr=errsfr*2.5#correct for dust
            if i == 5:#Couch et al
                #print line
                rah=float(line[5:7])
                ram=float(line[8:10])
                ras=float(line[11:18])
                decd=float(line[19:22])#ignore sign of declination
                decm=float(line[23:25])
                decs=float(line[26:32])

                ra=15.*(rah+ram/60. + ras/3600.)
                dec=-1.*(-1*decd+decm/60. + decs/3600.)
                x=(cra - ra)*3600.
                y=(cdec - dec)*3600.
                #d=sqrt((ra_c - ra)**2+(dec_c - dec)**2)*3600.
                m=float(line[33:39]) - distmod -kc#I mag
                #print rah, ram, ras, decd, decm, decs, ra,dec,m
                try:
                    zg=float(line[52:59])
                    if((zg > zmin) & (zg < zmax)):
                        f=float(line[60:67])
                        errf=float(line[68:75])
                        if ((f/errf) > 0):
                            
                            #sfr=f*((DL)**2)*9.48e-9
                            #errsfr=errf*((DL)**2)*9.48e-9
                            sfr=f*1.e-17*4*3.1415*DLcm**2*7.9e-42
                            errsfr=errf*1.e-17*4*3.1415*DLcm**2*7.9e-42

                            sfr=sfr*2.5#correct for dust
                            errsfr=errsfr*2.5
                        else:
                            continue
                except:
                    continue
            if i == 6:#Balogh et al 2002
                rah=float(line[5:7])
                ram=float(line[8:10])
                ras=float(line[11:17])
                decd=float(line[19:20])#ignore sign of declination
                decm=float(line[21:23])
                decs=float(line[24:29])
                ra=15*(rah+ram/60. + ras/3600.)
                dec=-1*(decd+decm/60. + decs/3600.)
                #print line[5:6],line[8:9],line[11:17]


                x=(cra - ra)*3600.
                y=(cdec - dec)*3600.
                #d=sqrt((ra_c - ra)**2+(dec_c - dec)**2)*3600.
                m=float(line[30:35]) - distmod -kc #I mag

                try:
                    zg=float(line[39:44])
                    if((zg > zmin) & (zg < zmax)):
                        f=float(line[45:50])
                        errf=float(line[51:54])
                        #print "made z cut"
                        if (f/errf > 0):
                            #print "made flux cut"
                            #sfr=f*1.*((DL)**2)*9.48e-9
                            #errsfr=errf*1.*((DL)**2)*9.48e-9
                            #print "HEY",f,DLcm
                            sfr=f*1.e-17*4*3.1415*DLcm**2*7.9e-42
                            errsfr=errf*1.e-17*4*3.1415*DLcm**2*7.9e-42
                            
                            sfr=sfr*2.5#correct for dust
                            errsfr=errsfr*2.5

                            #print "sfr = ",sfr
                            sfrb=float(line[65:70])/h100**2
                            errsfrb=float(line[71:76])/h100**2
                            #print "my/balogh = ",sfr/sfrb,sfrb/sfr
                            #sfr=sfrb
                            #errsfr=errsfrb
                except:
                    continue
                
            if i == 7:# ID(TK) ID(OC) ID(TT) X(pix)  Y(pix)  dRA(') dDec(')  Rc(')    B      error    R      error    z'     error    NB     error  (B-R)c   error  (B-z')c  error  (R-z')c  error  z'-NB912 error  EW(Ha)   error emitter SFR    error  z_phot  error+  error-  mem1  z_spec   qz  mem2
                x=f[5]*60.#convert from ' to "
                y=f[6]*60.#convert from ' to "
                m=f[10] - distmod -kc
                sfr=0
                if f[26] > 0:#check if galaxy is emitter
                    sfr=f[27]*2.5 #1 mag extinc
                    errsfr=f[28]*2.5 #1 mag extinc
            d=N.sqrt(x**2 + y**2)
            if (d <= r2arcsec):
                #print "within 0.5 R200!",m
                if m < mabscut:
                    sumsfr+=sfr
                    errsumsfr+=errsfr**2
                    allsfr.append(sfr)
                    #print "got one!",i,sfr,sumsfr
        infile.close()
        sumsfrhiz[i]=sumsfr
        errsumsfrhiz[i]=errsumsfr
        allsfr=N.array(allsfr,'f')
        #print "lit sfr", i, sumsfrhiz[i],errsumsfrhiz[i]
        if sumsfrhiz[i] == 0:
            sumsfrhiz[i]=1.
    return sumsfrhiz,errsumsfrhiz

            


            
def plothiz():
    sigmaminphiz=300.
    zhiz=N.array([.704,.748,.794,.845,.228,.32,.183,0.39],'f')
    sigmahiz=N.array([418,504,1018,415,1023,1390,1274,561],'f')
    #sumsfrhiz=N.array([85.3,45.8,360.7,78.,249.3,44.4,82.6,(253.)],'f')
    #sumsfrhiz=N.array([30.9,13.1,156.9,38.2,79.9,21.8,40.5,124.],'f')
    #applying Mr < -20.38
    sumsfrhiz=N.array([30.95,5.76,129.5,38.2,79.9,21.8,40.5,124.],'f')
    sumsfrhiz=sumsfrhiz/(h100)**2
    (sumsfrhiz,errsumsfrhiz)=getlitsfr(zhiz,sigmahiz)
    sfrcor=N.array([1.,1,1,1,1,2.8,2.8,1.],'f')#aperture correction
    volcor=N.array([1.,(1./.94),(1./.86),1,1,(1./.99),(1./.74),1.],'f')#aperture correction
    #sumsfrhiz=sumsfrhiz*sfrcor*volcor

    masshiz=12*(sigmahiz/1000.)**3.*1./N.sqrt(.7+.3*(1+zhiz)**3)/h100#x10^14 Mo
    kodamamasshiz=masshiz
    mratio=masshiz/kodamamasshiz
    sfrmhiz=(sumsfrhiz/masshiz)#correct for AGN contamination
    sumsfrhiz=sumsfrhiz
    clname = ["\ca","\cb","\cc","\cj", "Abell~2390","AC~114","Abell~1689","CL0024.0$+$1652"]
    refnumber = ["1","1","1","2", "3","4","5","6"]
    technique = ["I","I","I","I", "I","S","S","I"]
    symbols = [-4,-4,-4,-3,6,4,7,11]
    refsymb = [-4,-5,-6,-3,6,4,7,11,21,17]
    colors =  [ 2, 2, 2, 2,2,2,2, 2, 1,1]
    #ref = ["Finn et al. 2005","Finn et al. 2004","Balogh & Morris 2000","Couch et al. 2001", "Balogh et al. 2002","Kodama et al. 2004","C4 Clusters, sub>0","C4 Clusters, sub=0"]
    symbols = [-4,-5,-6,-3,6,4,7,11]
    #ref = ["This study","Finn et al. 2004","Balogh & Morris 2000","Couch et al. 2001", "Balogh et al. 2002","Kodama et al. 2004","Alternate M\dcl\u"]
    ref = ["Finn et al. 2005","Finn et al. 2004","Balogh & Morris 2000","Couch et al. 2001", "Balogh et al. 2002","Kodama et al. 2004","Alternate M\dcl\u"]
    ref = ["CL1040-1155","CL1054-1245","CL1216-1201","CLJ0023+0423B","Abell 2390","AC 114", "Abell 1689","CL 0024.0+1652","Alternate M\dcl\u"]


    yplabel="\gS SFR + 1 (M\d \(2281)\u yr\u-1\d)"
    psplotinit("hizsfrz.ps")
    #ymax=410.
    #ymin=-1.
    ymax=2.8
    ymin=-.2
    #ymin=-1.
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,1,ymin,ymax,0,20)
    ppgplot.pglab("z",yplabel,"")
    x1=c.z
    y1=N.log10(c.sfr05+1)
    ppgplot.pgsci(1)
    ppgplot.pgpt(x1,y1,21)
    drawptsub0(x1,y1)
    ppgplot.pgsci(2)
    x2=zhiz
    y2=N.log10(sumsfrhiz+1)
    pgpnts(x2,y2,symbols)
    xlabel = .5
    ylabel = .8
    ystep = .12
    dy=.05
    dx=.02
    ppgplot.pgsci(1)
    hizlabel(xlabel,ylabel,dx,dy,ystep,ref,refsymb,colors)
    ppgplot.pgend()

    psplotinit("hizsfrsigma.ps")
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(sigmaminphiz,sigmamaxp,ymin,ymax,0,20)
    ppgplot.pglab("\gs (km/s)",yplabel,"")
    x1=c.sigma
    ppgplot.pgpt(x1,y1,21)
    drawptsub0(x1,y1)
    ppgplot.pgsci(2)
    x2=sigmahiz
    pgpnts(x2,y2,symbols)
    xlabel = 1080
    ylabel = 3.
    ystep = .1
    dy=.06
    dx=30
    ppgplot.pgsci(1)
    #hizlabel(xlabel,ylabel,dx,dy,ystep,ref,refsymb,colors)
    ppgplot.pgend()

    
    psplotinit("hizz.ps")
    ylabelp="\gS SFR/(\gs/1000)\u3\d + 1 (M\d \(2281)\u yr\u-1\d/(km/s)\u3\d )"
    x1=c.z
    ymax=3.2
    ymin=-.2
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,1,ymin,ymax,0,20)
    ppgplot.pglab("z",ylabelp,"")
    y1=N.log10((c.sfr05)/(c.sigma/1000.)**3. + 1.)
    ppgplot.pgsci(1)
    ppgplot.pgpt(x1,y1,21)
    drawptsub0(x1,y1)

    ppgplot.pgsci(2)
    x2=zhiz
    #y=sfrmhiz
    y2=N.log10(sumsfrhiz/(sigmahiz/1000.)**3 + 1.)
    pgpnts(x2,y2,symbols)
    xlabel = .53
    ylabel = 1.
    ystep = .13
    dy=.04
    dx=.02
    ppgplot.pgsci(1)
    hizlabel(xlabel,ylabel,dx,dy,ystep,ref,refsymb,colors)
    ppgplot.pgend()

    
    psplotinit("hizsigma.ps")
    #y=N.log10(c.sfrmass)
    #y=(c.sfrmass)
    x1=c.sigma
    x2=sigmahiz
    xmax=max(x1)
    xmin=min(x1)
    ppgplot.pgbox("",0.0,0,"",0.0)
    ppgplot.pgenv(sigmaminphiz,sigmamaxp,ymin,ymax,0,20)
    ppgplot.pglab("\gs (km s\u-1 \d)","\gS SFR/(\gs/1000) + 1(M\d\(2281)\uyr\u-1\d/km/s)","")
    #ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(1)
    ppgplot.pgpt(x1,y1,21)
    drawptsub0(x1,y1)
    ppgplot.pgsci(2)
    pgpnts(x2,y2,symbols)
    xlabel = 1000
    ylabel = 1.
    ystep = .13
    dy=.04
    dx=10
    ppgplot.pgsci(1)
    #hizlabel(xlabel,ylabel,dx,dy,ystep,ref,refsymb,colors)
    ppgplot.pgsci(1)

    ppgplot.pgend()
    psplotinit("hizmcl.ps")
    y=sumit(6.,0.5,0)/c.mass
    ymin=min(y)
    print "min C4 within 0.5R = ",ymin
    #y=(N.clip(y,.01,1.E7))
    x=(c.mass)
    x2=masshiz
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(0,30,-1,80,0)
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
    ppgplot.pgend()

def hizlabel(xlabel,ylabel,dx,dy,ystep,ref,refsymb,colors):
    ppgplot.pgsls(1)
    ppgplot.pgsch(1.)
    for i in range(len(ref)):
        putlabelptc(xlabel,ylabel,dx,dy,ref[i],refsymb[i],colors[i])
        ylabel = ylabel - ystep

def plotsfrlx():
    lxmin=-.2
    lxmax=2.
    sfr=c.sfr2#sumit(3.,2.,0.)
    #sfr=sumitfixedrv(2000.,1.5,0.)
    sfr=N.compress(c.xrayflag > 0,sfr)
    lx=N.log10(N.compress(c.xrayflag > 0,c.Lx)+1)
    sub=N.compress(c.xrayflag > 0,c.sub)
    sig=N.compress(c.xrayflag > 0,c.sigma)
    rv=N.compress(c.xrayflag > 0,c.rvir)
    totrlum=N.compress(c.xrayflag > 0,c.totrlum)
    xmin=min(lx)
    xmax=max(lx)
    ymin=min(sfr)
    ymax=max(sfr)+10.
    
    psplotinit("sfrlx.ps")
    ppgplot.pgbox("",0.0,0,"",0.0)
    ppgplot.pgenv(lxmin,lxmax,-2,ymax,0,10)
    ppgplot.pglab("L\dX\u (10\u44\d ergs/s)","\gS SFR (M\d\(2281)\u yr\u-1\d)","")
    ppgplot.pgpt(lx,sfr,20)
    y=N.compress(sub < 1,sfr)
    x=N.compress(sub < 1,lx)
    ppgplot.pgpt(x,y,17)
    print "total number of x-ray clusters = ",N.sum(c.xrayflag)

    psplotinit("sigmalx.ps")
    ppgplot.pgbox("",0.0,0,"",0.0)
    #ppgplot.pgenv(lxmin,lxmax,200,sigmamaxp,0,10)
    ppgplot.pgenv(-2.8,2.8,-10,2000,0,10)
    ppgplot.pglab("L\dX\u (10\u44\d ergs/s)","\gs (km/s)","")
    ppgplot.pgpt(lx,sig,20)
    y=N.compress(sub < 1,sig)
    x=N.compress(sub < 1,lx)
    ppgplot.pgpt(x,y,17)
    siglit=[]
    lxlit=[]
    infile=open("/Users/rfinn/SDSS/xray/xue-wu-2000.dat",'r')
    for line in infile:
        if line.find('#') > -1:
            continue
        fields=line.split()
        siglit.append(float(fields[2]))#sigma in km/s
        lxlit.append(float(fields[4]))#Lx in 10^42 erg/s - same units as Miller
    x=N.log10(N.array(lxlit,'f'))
    y=N.array(siglit,'f')
    ppgplot.pgsci(3)
    #ppgplot.pgpt(x,y,4)
    infile.close()

    ppgplot.pgsci(1)
    infile=open("/Users/rfinn/SDSS/xray/madhavi.dat",'r')
    for line in infile:
        if line.find('#') > -1:
            continue
        fields=line.split()
        siglit.append(float(fields[8]))#sigma in km/s
        lxlit.append(float(fields[10]))#Lx in 10^42 erg/s - same units as Miller
    y=10.**(N.array(siglit,'f'))
    x=N.log10((70./50.)**2*10**((N.array(lxlit,'f'))-44.))
    ppgplot.pgpt(x,y,20)
    infile.close()

    psplotinit("rvirlx.ps")
    ppgplot.pgbox("",0.0,0,"",0.0)
    ppgplot.pgenv(lxmin,lxmax,1.,3.5,0,10)
    ppgplot.pglab("L\dX\u (10\u44\d ergs/s)","Rvir (Mpc)","")
    ppgplot.pgpt(lx,rv,20)
    y=N.compress(sub < 1,rv)
    x=N.compress(sub < 1,lx)
    ppgplot.pgpt(x,y,17)

    psplotinit("totrlumlx.ps")
    ymin=totrlumminp
    ymax=totrlummaxp
    ppgplot.pgbox("",0.0,0,"",0.0)
    ppgplot.pgenv(lxmin,lxmax,ymin,ymax,0,10)
    ppgplot.pglab("L\dX\u (10\u44\d ergs/s)","Total R Lum","")
    ppgplot.pgpt(lx,totrlum,20)
    y=N.compress(sub < 1,totrlum)
    x=N.compress(sub < 1,lx)
    ppgplot.pgpt(x,y,17)

    psplotinit("sfrsigmalx.ps")
    y=sfr/(sig/1000.)**3
    ppgplot.pgbox("",0.0,0,"",0.0)
    ppgplot.pgenv(lxmin,lxmax,-2,200,0,10)
    ppgplot.pglab("L\dX\u (10\u44\d ergs/s)","\gSSFR/\gs\u3\d","")
    ppgplot.pgpt(lx,y,20)
    y=N.compress(sub < 1,y)
    x=N.compress(sub < 1,lx)
    ppgplot.pgpt(x,y,17)
    ppgplot.pgend()

def plotsfrmabs():
    ylabel="SFR (M\d\(2281)\u yr\u-1\d)"
    psplotinit("sfrmabshiz.ps")

    x=N.compress((g.Mabs > -999.) & (g.Mabs < 0.) &  (g.agn < 1) & (g.dr < 0.5*g.clusterrvir) & (g.dv < 3.*g.clustersigma) & (g.ew > 10.),g.Mabs)

    y=N.compress((g.Mabs > -999.) & (g.Mabs < 0.) &  (g.agn < 1) & (g.dr < 0.5*g.clusterrvir) & (g.dv < 3.*g.clustersigma) & (g.ew > 10.),g.sfr)


    
    ymin=-.9
    ymax=1.9
    xmin=-23.
    xmax=-17.5

    ppgplot.pgbox("",0.0,0,"",0.0,1)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,20)
    ppgplot.pglab("M\dR\u",ylabel,"")
    ppgplot.pgsci(1)
    
    xmin=-22.
    xmax=-19.
    nbin=6.
    #(xbin,ybin,ybinerr)=biniterr(x,y,nbin)
    (xbin,ybin,ybinerr)=binitbins(xmin,xmax,nbin,x,y)
    #print "HEY!!!!",xbin,ybin,ybinerr
    #print "sigmbin = ",c.sigmabin
    #print "xbin = ",xbin
    #ybin=ybin+1
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,18)
    ybin1=ybin
    xedi=[]
    yedi=[]
    input=open('/Users/rfinn/clusters/final/cl1040mabssfr.dat','r')
    for line in input:
	t=line.split()
	xedi.append(float(t[0]))
	yedi.append(float(t[1]))
    input.close()
    input=open('/Users/rfinn/clusters/final/cl1054mabssfr.dat','r')
    for line in input:
	t=line.split()
	xedi.append(float(t[0]))
	yedi.append(float(t[1]))
    input.close()
    input=open('/Users/rfinn/clusters/final/cl1216mabssfr.dat','r')
    for line in input:
	t=line.split()
	xedi.append(float(t[0]))
	yedi.append(float(t[1]))
    input.close()
    x=N.array(xedi,'f')
    y1=N.array(yedi,'f')
    #(xbin,ybin,ybinerr)=biniterr(x,y1,nbin)
    (xbin,ybin,ybinerr)=binitbins(xmin,xmax,nbin,x,y1)
    #print "HEY!!!!",xbin,ybin,ybinerr
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,16)
    ppgplot.pgend()
    for i in range(len(ybin)):
	print "Ratio of ediscs/SDSS SFR at mabs = ",xbin[i]," is ",10.**(ybin[i] - ybin1[i]),ybin[i],ybin1[i]
def plotewmabshiz():
    ylabel="EW(H\ga) (\(2078))"
    psplotinit("ewmabshiz.ps")

    x=N.compress((g.Mabs > -999.) & (g.Mabs < 0.) &  (g.agn < 1) & (g.dr < 0.5*g.clusterrvir) & (g.dv < 3.*g.clustersigma) & (g.ew > 10.),g.Mabs)

    y=N.compress((g.Mabs > -999.) & (g.Mabs < 0.) &  (g.agn < 1) & (g.dr < 0.5*g.clusterrvir) & (g.dv < 3.*g.clustersigma) & (g.ew > 10.),g.ew)


    
    ymin=.9
    ymax=2.3
    xmin=-23.
    xmax=-17.5

    ppgplot.pgbox("",0.0,0,"",0.0,1)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,20)
    ppgplot.pglab("M\dR\u",ylabel,"")
    ppgplot.pgsci(1)
    
    xmin=-22.
    xmax=-19.
    nbin=6.
    #(xbin,ybin,ybinerr)=biniterr(x,y,nbin)
    (xbin,ybin,ybinerr)=binitbins(xmin,xmax,nbin,x,y)
    #print "HEY!!!!",xbin,ybin,ybinerr
    #print "sigmbin = ",c.sigmabin
    #print "xbin = ",xbin
    #ybin=ybin+1
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,18)

    xedi=[]
    yedi=[]
    input=open('/Users/rfinn/clusters/final/cl1040mabssfr.dat','r')
    for line in input:
	t=line.split()
	xedi.append(float(t[0]))
	yedi.append(float(t[2]))
    input.close()
    input=open('/Users/rfinn/clusters/final/cl1054mabssfr.dat','r')
    for line in input:
	t=line.split()
	xedi.append(float(t[0]))
	yedi.append(float(t[2]))
    input.close()
    input=open('/Users/rfinn/clusters/final/cl1216mabssfr.dat','r')
    for line in input:
	t=line.split()
	xedi.append(float(t[0]))
	yedi.append(float(t[2]))
    input.close()
    x=N.array(xedi,'f')
    y1=N.array(yedi,'f')
    #(xbin,ybin,ybinerr)=biniterr(x,y1,nbin)
    (xbin,ybin,ybinerr)=binitbins(xmin,xmax,nbin,x,y1)
    #print "HEY!!!!",xbin,ybin,ybinerr
    errorylog(xbin,ybin,ybinerr)
    ybin=N.log10(ybin)
    ppgplot.pgpt(xbin,ybin,16)
    ppgplot.pgend()

class sim:
    def __init__(self):
        self.sigma = []
    def readdatafile(self):
	input=open('/Users/rfinn/clusters/millenium/millen.dat','r')
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


	input=open('/Users/rfinn/clusters/millenium/millen.dat','r')
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
                    break
                #if line.find('*') > -1:#skip any galaxies w/ *** in one or more fields
                #    continue
                #fields=line.split()
                #if float(fields[8]) > (-100.): #get rid of large negative values
                #    sum += float(fields[8])
            #self.sumsfr.append(float(sum))
    def creadfieldfiles(self,clusters):
        print "number of clusters = ",len(clusters)

        for i in range(len(clusters)):
            prefix=clusters[i].split('.')
            t=prefix[0].split('_')
	    cid=int(str(t[2])+str(t[3]))
            j=0
            for line in open(clusters[i]):
                if line.find('#') > -1:
                    continue
                if j == 0:
                    fields=line.split()
                    sum=0
                    if float(fields[5]) < 0:#skip clusters with masses below sigmamin
                        break
                    #self.id.append(int(fields[0]))
		    self.id.append(cid)
                    self.ra.append(float(fields[1]))
                    self.dec.append(float(fields[2]))
                    self.z.append(float(fields[3]))
                    self.rvir.append(float(fields[4]))#Rvirial in Mpc
                    self.sigma.append(float(fields[5]))#velocity dispersion in km/s
                    j=1            
                    break
                #if line.find('*') > -1:#skip any galaxies w/ *** in one or more fields
                #    continue
                #fields=line.split()
                #if float(fields[8]) > (-100.): #get rid of large negative values
                #    sum += float(fields[8])

    def convarray(self):
        self.id = N.array(self.id,'i')
        self.ra = N.array(self.ra,'f')
        self.dec  = N.array(self.dec,'f')
        self.z = N.array(self.z,'f')
        self.rvir = N.array(self.rvir,'f')
        self.rvir=self.rvir/h100*.7
        self.sigma = N.array(self.sigma,'f')
        #self.sumsfr = N.array(self.sumsfr,'f')
        self.mass = N.zeros(len(self.sigma),'f')
        self.sfrmass = N.zeros(len(self.sigma),'f')

    def limitz(self): #limit cluster sample to z > zmin
        print len(self.id),len(self.z)
        #for i in range(len(self.id)):
        #    print i,self.id[i],self.z[i]
        print "zmin, zmax = ",zmin,zmax
        self.id = N.compress((self.z > zmin) & (self.z < zmax),self.id)
        self.ra = N.compress((self.z > zmin) & (self.z < zmax),self.ra)
        self.dec = N.compress((self.z > zmin) & (self.z < zmax),self.dec)
        self.rvir = N.compress((self.z > zmin) & (self.z < zmax),self.rvir)
        self.sigma = N.compress((self.z > zmin) & (self.z < zmax),self.sigma)
        #self.sumsfr = N.compress((self.z > zmin) & (self.z < zmax),self.sumsfr)
        self.mass = N.compress((self.z > zmin) & (self.z < zmax),self.mass)
        self.sfrmass = N.compress((self.z > zmin) & (self.z < zmax),self.sfrmass)
        self.z = N.compress((self.z > zmin) & (self.z < zmax),self.z)
        print len(self.z), " galaxies after redshift cut"

    def recordchisq(self,chisq):
        self.chisq=N.zeros(len(self.z),'f')
        for i in range(len(self.chisq)):
            self.chisq[i]=chisq[i]
    def recordchisqfield(self):
        self.chisq=N.ones(len(self.z),'f')
            
    def limitrichness(self):
        print "minimum number of galaxies = ",ngalmin
        self.ngal = N.zeros(len(self.id),'f')
        for i in range(len(self.id)):
            #self.ngal[i] = len(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.Mabs < mabscut) & (g.memb > 0), g.sfr))
            #for j in c.memberids
            self.ngal[i] = len(N.compress((abs(g.clusterid - c.id[i])< 1) & (g.Mabs < mabscut) & (g.memb > 0), g.sfr))
        self.id = N.compress(self.ngal > ngalmin,self.id)
        self.ra = N.compress(self.ngal > ngalmin,self.ra)
        self.dec = N.compress(self.ngal > ngalmin,self.dec)
        self.rvir = N.compress(self.ngal > ngalmin,self.rvir)
        self.sigma = N.compress(self.ngal > ngalmin,self.sigma)
        #self.sumsfr = N.compress(self.ngal > ngalmin,self.sumsfr)
        self.mass = N.compress(self.ngal > ngalmin,self.mass)
        self.sfrmass = N.compress(self.ngal > ngalmin,self.sfrmass)
        self.z = N.compress(self.ngal > ngalmin,self.z)
        print len(self.z), " clusters after richness cut"

    def limitsigma(self):
        self.id = N.compress((self.sigma > sigmamin) & (self.sigma < sigmamax),self.id)
        self.ra = N.compress((self.sigma > sigmamin) & (self.sigma < sigmamax),self.ra)
        self.dec = N.compress((self.sigma > sigmamin) & (self.sigma < sigmamax),self.dec)
        self.rvir = N.compress((self.sigma > sigmamin) & (self.sigma < sigmamax),self.rvir)
        #self.sumsfr = N.compress((self.sigma > sigmamin) & (self.sigma < sigmamax),self.sumsfr)
        self.mass = N.compress((self.sigma > sigmamin) & (self.sigma < sigmamax),self.mass)
        self.sfrmass = N.compress((self.sigma > sigmamin) & (self.sigma < sigmamax),self.sfrmass)
        self.z = N.compress((self.sigma > sigmamin) & (self.sigma < sigmamax),self.z)
        self.sigma = N.compress((self.sigma > sigmamin) & (self.sigma < sigmamax),self.sigma)

        print len(self.z), " clusters after sigma cut"

    def limitsubstructureold(self):
        submax=1
        isub=1#cut on internal substructure
        self.id = N.compress(self.sub < submax,self.id)
        self.ra = N.compress(self.sub  < submax,self.ra)
        self.dec = N.compress(self.sub < submax,self.dec)
        self.rvir = N.compress(self.sub < submax,self.rvir)
        #self.sumsfr = N.compress(self.sub < submax,self.sumsfr)
        self.mass = N.compress(self.sub < submax,self.mass)
        self.sfrmass = N.compress(self.sub < submax,self.sfrmass)
        self.z = N.compress(self.sub < submax,self.z)
        self.sigma = N.compress(self.sub < submax,self.sigma)
        self.totrlum=N.compress(self.sub < submax,self.totrlum)
        self.rh1500=N.compress(self.sub < submax,self.rh1500)
        self.gauss=N.compress(self.sub < submax,self.gauss)
        self.chisq=N.compress(self.sub < submax,self.chisq)
        self.sub = N.compress(self.sub < submax,self.sub)
        
        #print len(self.z)," clusters after sub flag cut"

        #if (isub > 0):#now cut based on internal substructure
        #    #y=.32*N.log10(self.totrlum)+2.102
        #    #y=.2*N.log10(self.totrlum)+2.42
        #    y=.5*N.log10(self.totrlum/100.)+2.75
        #    sig=N.log10(self.sigma)
        #    for i in range(len(self.sigma)):
        #        if sig[i] > y[i]:
        #            print "inflated sigma for cluster ",self.id[i]," diff = ",sig[i]-y[i]

        #    if cdir.find('balogh1') > -1:
        #        print "hey, in balogh1 directory!"
        #        outfile=open('final-cluster-list-gaussian','w')
        #        for i in range(len(self.sigma)):
        #            if sig[i] > y[i]:
        #                print "inflated sigma for cluster ",self.id[i]
        #            if sig[i] < y[i]:
        #                st=str(c.id[i])+"\n"
        #                #print st
        #                outfile.write(st)
        #        outfile.close()
        #    self.id = N.compress(sig < y,self.id)
        #    self.ra = N.compress(sig < y,self.ra)
        #    self.dec = N.compress(sig < y,self.dec)
        #    self.rvir = N.compress(sig < y,self.rvir)
        #    self.sumsfr = N.compress(sig < y,self.sumsfr)
        #    self.mass = N.compress(sig < y,self.mass)
        #    self.sfrmass = N.compress(sig < y,self.sfrmass)
        #    self.z = N.compress(sig < y,self.z)
        #    self.totrlum=N.compress(sig < y,self.totrlum)
        #    self.rh1500=N.compress(sig < y,self.rh1500)
        #    self.sub = N.compress(sig < y,self.sub)
        #    self.sigma = N.compress(sig < y,self.sigma)
        #    self.gauss = N.compress(sig < y,self.gauss)
        #print len(self.z), " clusters after substructure cut"

    def limitsubstructure(self):
        #changed to limit on reduced chi-square for gaussian fit
        #i.e. cut clusters that are significanlty different from gaussian
        chisqmax=2.0
        self.id = N.compress(self.chisq < chisqmax,self.id)
        self.ra = N.compress(self.chisq  < chisqmax,self.ra)
        self.dec = N.compress(self.chisq < chisqmax,self.dec)
        self.rvir = N.compress(self.chisq < chisqmax,self.rvir)
        #self.sumsfr = N.compress(self.chisq < chisqmax,self.sumsfr)
        self.mass = N.compress(self.chisq < chisqmax,self.mass)
        self.sfrmass = N.compress(self.chisq < chisqmax,self.sfrmass)
        self.z = N.compress(self.chisq < chisqmax,self.z)
        self.sigma = N.compress(self.chisq < chisqmax,self.sigma)
        self.totrlum=N.compress(self.chisq < chisqmax,self.totrlum)
        self.rh1500=N.compress(self.chisq < chisqmax,self.rh1500)
        self.sub = N.compress(self.chisq < chisqmax,self.sub)
        self.gauss = N.compress(self.chisq < chisqmax,self.gauss)
        self.chisq = N.compress(self.chisq < chisqmax,self.chisq)
        print len(self.z)," clusters after gaussian cut"


    def assigngaussold(self):
        #infile=open("/Users/rfinn/SDSS/balogh1/final-cluster-list-gaussian",'r')
        infile=open("/Users/rfinn/SDSS/baloghDR3/gaussian-clusters",'r')
        gaussid=[]
        for line in infile:
            if line.find('#') > -1:
                continue
            fields=line.split()
            gaussid.append(float(fields[0]))
        print "number of gaussian clusters from file = ",len(gaussid)
        self.gauss = N.zeros(len(self.id),'i')
        cindex=N.arange(0,len(self.id),1)
        sortcindex=N.take(cindex,N.argsort(self.id))
        sortcid=N.take(self.id,N.argsort(self.id))
        deltaid=0.1
        #gaussid is already sorted numerically
        for i in range(len(gaussid)):
            (match,flag)=findmatch(gaussid[i],sortcid,deltaid)
            try:
                for j in match:
                    k=sortcindex[int(j)]
                    self.gauss[k]=1.
                #if len(match) > 1:
                #    print "Error in assigngauss - more than 1 match to cluster id"
            except TypeError:
                continue

    def assigngauss(self):
        maxchisq=1.
        self.gauss = N.zeros(len(self.id),'f')
        for i in range(len(self.chisq)):
            if self.chisq[i] <= maxchisq:
                self.gauss[i]=1.
    def assignlumbins(self): #bin by total luminosities
        nbin=7
        binnumb=N.arange(0,nbin,1,'i')
        lbin=N.zeros((nbin),'f')
        sigmabin=N.zeros((nbin),'f')
        #print self.gauss
        #print len(self.gauss),len(self.totrlum)
        totrlum=N.compress(self.gauss > 0.,self.totrlum)
        sigma=N.compress(self.gauss > 0.,self.sigma)
        sortsigma=N.take(sigma,N.argsort(totrlum))
        sortlum=N.take(totrlum,N.argsort(totrlum))#sort luminosities in order to define bins
        #sortgauss=N.take(self.gauss,N.argsort(self.totrlum))

        nlum=len(sortlum)-1
        for i in range(nbin):#define luminosity bins
            nmin=int((i)*1.*nlum/(1.*nbin))
            nmax=int((i+1)*1.*nlum/(1.*nbin))
            #print i,nmin,nmax,len(sortlum)
            lbin[i]=sortlum[int((i+1)*1.*nlum/(1.*nbin))]
            sigmabin[i]=N.average(sortsigma[nmin:nmax])
        lbinnumb=N.zeros(len(self.totrlum),'i') #lumin bin number
        self.sigmabin=N.zeros(len(self.totrlum),'f')#defaults to last bin
        self.r200bin=N.zeros(len(self.totrlum),'f')
        for i in range(len(self.totrlum)):#assign each cluster to lumin bin
            for j in range(nbin):
                if self.totrlum[i] < lbin[j]:
                    lbinnumb[i]=j
                    self.sigmabin[i]=sigmabin[j]
                    break
            #print "hey ",i,self.sigmabin[i],lbinnumb[i]
        self.r200bin=1.73*self.sigmabin/1000.*1./N.sqrt(.7+.3*(1+self.z)**3)/h100

    def calcrichness(self):
        #self.richness=richness(nsig,nr,0,mabscut)
        self.totlum05=qsumit(g.rlum,nsig*self.sigma,0.5*self.r200,0,mabscut)/1.e11
        self.totlum1=qsumit(g.rlum,nsig*self.sigma,1.*self.r200,0,mabscut)/1.e11
        self.totlum2=qsumit(g.rlum,nsig*self.sigma,2.*self.r200,0,mabscut)/1.e11
        
        #dr=N.ones(len(self.r200),'f')
        #self.totlum1Mpc=qsumit(g.rlum,nsig*self.sigma,dr,0,mabscut)
        #self.richness1Mpc=countit(nsig*self.sigma,dr,0,mabscut)

    def calcsfrstellmass(self):#calculate sfr and stellar mass within .5, 1, 2Rv
#        self.sumsfr2 = N.zeros(len(self.sigma),'f')
        (self.sfr2,t)=sumitmemb(nsig,2.,0)
        (self.sfr3,t)=sumitmemb(nsig,3.,0)
        (self.sfr05,t)=sumitmemb(nsig,0.5,0)
        (self.sfr05hiz,t)=sumitmemb(3.,0.5,0)
        (self.sfr1,self.avesfr)=sumitmemb(nsig,1,0)
        self.sfr05bin=sumitbin(nsig,.5,0)
        self.sfr1bin=sumitbin(nsig,1.,0)
        self.sfr2bin=sumitbin(nsig,2.,0)
        self.sfr3bin=sumitbin(nsig,3.,0)
        #self.sfr1Mpc=sumitfixedr(3.,1.,0)#dv,r(Mpc),agn?
        self.stellmass05=sumstellmass(nsig,0.5)
        self.stellmass1=sumstellmass(nsig,1.)
        self.stellmass2=sumstellmass(nsig,2.)
        self.stellmass3=sumstellmass(nsig,3.)

        self.stellmass05bin=sumstellmassbin(nsig,0.5)
        self.stellmass1bin=sumstellmassbin(nsig,1.)
        self.stellmass2bin=sumstellmassbin(nsig,2.)
        self.stellmass3bin=sumstellmassbin(nsig,3.)

        (self.sffrac05,sf,tot)=sffracmemb(nsig,0.5,0)
        (self.sffrac1,sf,tot)=sffracmemb(nsig,1.,0)
        (self.sffrac2,sf,tot)=sffracmemb(nsig,2.,0)
        self.sffrac=self.sffrac1
        (self.kauffmann,sf,tot)= (stellmassfracmemb(nsig,1.))#frac of galaxies w/kauffmann stellar mass
        #print "c.kauffmann = ",self.kauffmann
        if (completenesscorr > 0):
            for i in range(len(self.sfr05)):
                try:
                    self.sfr05[i]=self.sfr05[i]/self.compl05[i]
                except ZeroDivisionError:
                    self.sfr05[i]=self.sfr05[i]
                try:
                    self.sfr05hiz[i]=self.sfr05hiz[i]/self.compl05[i]
                except ZeroDivisionError:
                    self.sfr05hiz[i]=self.sfr05hiz[i]
                try:
                    self.sfr1[i]=self.sfr1[i]/self.compl1[i]
                except ZeroDivisionError:
                    self.sfr1[i]=self.sfr1[i]
                try:
                    self.sfr2[i]=self.sfr2[i]/self.compl2[i]
                except ZeroDivisionError:
                    self.sfr2[i]=self.sfr2[i]
                try:
                    self.stellmass05[i]=self.stellmass05[i]/self.compl05[i]
                except ZeroDivisionError:
                    self.stellmass05[i]=self.stellmass05[i]
                try:
                    self.stellmass1[i]=self.stellmass1[i]/self.compl1[i]
                except ZeroDivisionError:
                    self.stellmass1[i]=self.stellmass1[i]
                try:
                    self.stellmass2[i]=self.stellmass2[i]/self.compl2[i]
                except ZeroDivisionError:
                    self.stellmass2[i]=self.stellmass2[i]

    def calcmass(self):
        self.mass=9.78*(self.sigma/1000.)**3.*1./N.sqrt(.7+.3*(1+self.z)**3)
    def calcr200(self):
        self.r200=N.zeros(len(self.z),'f')
        self.r200=1.73*self.sigma/1000.*1./N.sqrt(.7+.3*(1+self.z)**3)/h100
        #self.r200=1.26*self.r200#adjust density contrast to 100x at z=0
        #self.rvir=self.r200
    def calcsfrmass(self):
        self.sfrmasshiz = N.zeros(len(self.sigma),'f')
        self.sfrmassbalogh = N.zeros(len(self.sigma),'f')
        self.sfrmass=self.sfr1/self.mass
        self.sfrtotlum=self.sfr1/self.totlum1
        self.sfrmasshiz=self.sfr05hiz/self.mass
        self.sfrmassbalogh=self.sfr05hiz/self.mass
        #self.sfr1Mpcmass=self.sfr1Mpc/self.mass
        #self.sfr1Mpctotlum=self.sfr1Mpc/self.totlum1Mpc

    def calcvolume(self):
        self.vol05 = N.zeros(len(self.sigma),'f')
        self.vol1 = N.zeros(len(self.sigma),'f')
        self.vol2 = N.zeros(len(self.sigma),'f')
        self.vol2 = N.zeros(len(self.sigma),'f')
        self.vol05=4./3.*3.1415*(0.5*self.rvir)**3
        self.vol1=4./3.*3.1415*(1.*self.rvir)**3
        self.vol2=4./3.*3.1415*(2.*self.rvir)**3
        self.vol3=4./3.*3.1415*(3.*self.rvir)**3        
    def assignsub(self):#assign substructure 0 or 1
        sub0=[]
        cname=[]
        rlum=[]
        rh=[]
        #0-RA_MEAN 1-DEC_MEAN 2-RA_BCG 3-DEC_BCG 4-Z 5-BIWT1500 6-WBIWT1500 7-WMAG1500 8-RH1500 9-RJK1500 10-RV1500 11-RHO200_1500 12-KEEP 13-SINGLE  14-CLUSTER_ID 15-SUB 16-Lum_r
  
        #for line in open("../balogh1/clusters-sub0"):#read in list of clusters with substructure flag = 0
        #for line in open("../cluster_catalogs/combined_dr2_062404_noheader.dat"):#read in list of clusters with substructure flag = 0
        for line in open("../cluster_catalogs/sdss_c4_dr3_unpublished.dat"):#read in list of clusters with substructure flag = 0
            fields=line.split()
            rh.append(float(fields[8]))
            cname.append(float(fields[14]))
            sub0.append(float(fields[15]))
            rlum.append(float(fields[16]))
        #assign substructure flag
        self.sub=N.ones(len(self.id),'f')
        self.totrlum=N.zeros(len(self.sub),'f')
        self.rh1500=N.zeros(len(self.sub),'f')

        cnameindex=N.arange(0,len(cname),1)
        cnameindexsort=N.take(cnameindex,N.argsort(cname))
        cnamesort=N.take(cname,N.argsort(cname))
        deltaid=0.01
        for i in range(len(self.id)):
            name=self.id[i]
	    #if cdir.find('field') > -1:
	#	name=int(float(name)/10.)#field ids have extra number (1-5 appended to end), so divide by 10
            (match,flag)=findmatch(name,cnamesort,deltaid)#look for cid in cluster list
            #match has indices of cnamesort, cindexsort(match) gives indices
            #of cname
            try:
                for j in match:
                    k=cnameindexsort[int(j)]#index for cname,sub0
                    self.sub[i]=float(sub0[k])
                    self.totrlum[i]=float(rlum[k])/(1.e9)
                    self.rh1500[i]=float(rh[k])
                #if len(match) > 1:
                #    print "Error in getxray - more than 1 match for ",name,len(match),match
            except TypeError:
                print "no match in assignsub for cluster ",name
                continue

        #for i in range(len(self.id)):
        #    a=self.id[i]
        #    for j in range(len(cname)):
        #        b=cname[j]
        #        diff=abs(float(b)-float(a))
        #        if diff < 0.1:#match cluster ids
        #            self.sub[i]=float(sub0[j])
        #            self.totrlum[i]=float(rlum[j])/(10.**9)
        #            self.rh1500[i]=float(rh[j])
        #            #print i,"found a match",a,b,self.totrlum[i]
        #            break
        #self.rvir=self.rh1500
	temp=self.rh1500/self.rvir
	print "Average rh1500/rvir = ",N.average(temp)
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

    def getxray(self):
        self.xrayflag = N.zeros(len(self.ra),'f')
        self.Lx = N.zeros(len(self.ra),'f')
        xrayin=open("/Users/rfinn/SDSS/xray/xray.dat",'r')
        cid=[]
        lx=[]
        for line in xrayin:
            if line.find('#') > -1:
                continue
            fields=line.split()
            cid.append(float(fields[0]))
            lx.append(float(fields[1]))#10^37 W = 10^42 ergs/s
        xrayin.close()
        cindex=N.arange(0,len(self.id),1)
        cindexsort=N.take(cindex,N.argsort(self.id))
        cidsort=N.take(self.id,N.argsort(self.id))
        deltaid=0.1
        for i in range(len(cid)):
            name=cid[i]
            (match,flag)=findmatch(name,cidsort,deltaid)#look for cid in cluster list
            try:
                for j in match:
                    k=cindexsort[int(j)]
                    self.xrayflag[k]=1.
                    self.Lx[k]=lx[i]
                if len(match) > 1:
                    print "Error in getxray - more than 1 match for cid",name
            except TypeError:
                continue

            #for j in range(len(cid)):
            #    diff=abs(self.id[i]-cid[j])
            #    if diff < .9:
            #        #print "matched xray data for ",cid[j],self.id[i]
            #        self.xrayflag[i]=1.
            #        self.Lx[i]=lx[j]
            #    continue
    def readdatafile(self,file):
        output=open(file,'r')
        i=-1
        for line in output:
            if (i == -1):
                fields=line.split()
                ncl=int(fields[1])#number of clusters
                print "number of clusters = ",ncl
                self.id=N.zeros(ncl,'f')
                self.z=N.zeros(ncl,'f')
                self.rvir=N.zeros(ncl,'f')
                self.r200=N.zeros(ncl,'f')
                self.sigma=N.zeros(ncl,'f')
                self.totrlum=N.zeros(ncl,'f')
                self.Lx=N.zeros(ncl,'f')
                self.xrayflag=N.zeros(ncl,'f')
                self.stellmass05=N.zeros(ncl,'d')
                self.stellmass1=N.zeros(ncl,'d')
                self.stellmass2=N.zeros(ncl,'d')
                self.stellmass3=N.zeros(ncl,'d')
                self.sfr05=N.zeros(ncl,'f')
                self.sfr05hiz=N.zeros(ncl,'f')
                self.sfr1=N.zeros(ncl,'f')
                self.sfr2=N.zeros(ncl,'f')
                self.sfr3=N.zeros(ncl,'f')
                self.sffrac05=N.zeros(ncl,'f')
                self.sffrac1=N.zeros(ncl,'f')
                self.sffrac2=N.zeros(ncl,'f')
                self.sub=N.zeros(ncl,'f')
                self.totlum05=N.zeros(ncl,'f')
                self.totlum1=N.zeros(ncl,'f')
                self.totlum2=N.zeros(ncl,'f')
                self.r200bin=N.zeros(ncl,'f')
                self.sigmabin=N.zeros(ncl,'f')
                self.sfr05bin=N.zeros(ncl,'f')
                self.sfr1bin=N.zeros(ncl,'f')
                self.sfr2bin=N.zeros(ncl,'f')
                self.sfr3bin=N.zeros(ncl,'f')
                self.stellmass05bin=N.zeros(ncl,'f')
                self.stellmass1bin=N.zeros(ncl,'f')
                self.stellmass2bin=N.zeros(ncl,'f')
                self.stellmass3bin=N.zeros(ncl,'f')
                self.chisq=N.zeros(ncl,'f')
                self.kauffmann=N.zeros(ncl,'f')
                self.ra=N.zeros(ncl,'f')
                self.dec=N.zeros(ncl,'f')
                self.ngal05=N.zeros(ncl,'f')
                self.ngal1=N.zeros(ncl,'f')
                self.ngal2=N.zeros(ncl,'f')
                self.compl05=N.zeros(ncl,'f')
                self.compl1=N.zeros(ncl,'f')
                self.compl2=N.zeros(ncl,'f')
                i += 1
                continue
            if line.find('#') > -1:#skip header
                continue
            f=line.split()
            #print f
            #f=N.array(f,'e')

            for j in range(len(f)):
                #print i,f[i]
                f[j]=float(f[j])
            (self.id[i],self.z[i],self.rvir[i],self.r200[i],self.sigma[i],self.totrlum[i],self.Lx[i],self.xrayflag[i],self.stellmass05[i],self.stellmass1[i],self.stellmass2[i],self.stellmass3[i],self.sfr05[i],self.sfr1[i],self.sfr2[i],self.sfr3[i],self.sffrac05[i],self.sffrac1[i],self.sffrac2[i],self.sub[i],self.totlum05[i],self.totlum1[i],self.totlum2[i],self.r200bin[i],self.sigmabin[i],self.sfr05bin[i],self.sfr1bin[i],self.sfr2bin[i],self.sfr3bin[i],self.stellmass05bin[i],self.stellmass1bin[i],self.stellmass2bin[i],self.stellmass3bin[i],self.chisq[i],self.sfr05hiz[i],self.kauffmann[i],self.ra[i],self.dec[i],self.ngal05[i],self.ngal1[i],self.ngal2[i],self.compl05[i],self.compl1[i],self.compl2[i])=(f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18],f[19],f[20],f[21],f[22],f[23],f[24],f[25],f[26],f[27],f[28],f[29],f[30],f[31],f[32],f[33],f[34],f[35],f[36],f[37],f[38],f[39],f[40],f[41],f[42],f[43])
            i += 1
        output.close()

           
    def getsdsscompleteness(self):  #get number of photometric sources and number of spectroscopic sources within 2R200 and use ratio to estimate completeness of spectroscopic sampling
        ntry=0
        self.nphot05=N.zeros(len(c.z),'f')#number of objects with sdss photometry w/in 0.5R200
        self.nphot1=N.zeros(len(c.z),'f')
        self.nphot2=N.zeros(len(c.z),'f')
        self.nspec05=N.zeros(len(c.z),'f')#number of objects with sdss spec w/in 0.5R200
        self.nspec1=N.zeros(len(c.z),'f')
        self.nspec2=N.zeros(len(c.z),'f')
        self.compl05=N.zeros(len(c.z),'f')#completeness w/in 0.5R200
        self.compl1=N.zeros(len(c.z),'f')
        self.compl2=N.zeros(len(c.z),'f')
        nphot=N.zeros([4,len(c.z)],'f')
        nspec=N.zeros([4,len(c.z)],'f')
        print "elapsed time = ",time.clock()-starttime
        #t=time.clock()-starttime
        try:
            for i in range(len(c.z)):
                print "calculating completness for cluster ",i,"out of ",len(c.z)
                dA=DA(self.z[i],h100)
                r200arcmin=self.r200[i]*1000./dA/60.
                drsearch=2.*r200arcmin#2xR200 in arcmin for sdss query
                mr=mabscut + 5.*N.log10(dL(self.z[i],h100))+25.
                print "ra, dec, dr, mr = %12.8f %12.8f %8.3f %5.2f" % (self.ra[i],self.dec[i],drsearch,mr)
                query="select n.distance from galaxy g, dbo.fGetNearbyObjEq(%12.8f,%12.8f,%8.3f) n where g.objID = n.objID and g.r< %5.2f order by distance" % (self.ra[i],self.dec[i],drsearch,mr)
                lines=sqlcl.query(query).readlines()
                print "got number+1 phot objects = ",len(lines)
                dr=[]
                if (len(lines) > 1):
                    for line in lines[1:]:
                        f=line.split(',')
                        dr.append(float(f[0]))
                dr=N.array(dr,'f')
                binnumb=(dr/(2.*r200arcmin)*4.)
                for j in binnumb:
                    if ((j - 4.) > 0):
                        j=3.
                    j=int(j)
                    nphot[j][i]=nphot[j][i] + 1.
                        
                query="select n.distance from galaxy g, specobj s, dbo.fGetNearbyObjEq(%12.8f,%12.8f,%8.3f) n where g.objid = s.bestobjid and g.objID = n.objID and g.r< %5.2f order by distance" % (self.ra[i],self.dec[i],drsearch,mr)
                lines=sqlcl.query(query).readlines()
                dr=[]
                print "got number + 1 of spec objects = ",len(lines)
                if (len(lines) > 1.):
                    for line in lines[1:]:
                        f=line.split(',')
                        dr.append(float(f[0]))
                dr=N.array(dr,'f')
                binnumb=(dr/(2.*r200arcmin)*4.)
                for j in binnumb:
                    if ((j - 4.) > 0):
                        j=3.
                    j=int(j)
                    nspec[j][i]=nspec[j][i] + 1.
                try:
                    print "cluster ",i, " completeness w/in 0.5xR200 = ",nspec[0][i]/nphot[0][i]
                except ZeroDivisionError:
                    print "cluster ",i, " completeness w/in 0.5xR200 = NO GALAXIES"
                    print "setting completeness to 1. and checking other radii"
                    for k in range(4):
                        if (nphot[k][i] < 1.):
                            nphot[k][i]=1.
                            nspec[k][i]=1.
        finally:
            nclusters=i
            print "EXITING - FINISHED ",nclusters," clusters"
            output1=open('nphot.dat','w')
            for i in range(nclusters):
                s0=' '
                s=str(nspec[0][i])+s0+str(nphot[0][i])+s0+str(nspec[1][i])+s0+str(nphot[1][i])+s0+str(nspec[2][i])+s0+str(nphot[2][i])+s0+str(nspec[3][i])+s0+str(nphot[3][i])
            output1.close()

            self.nphot05=nphot[0][:]
            self.nphot1=nphot[1][:] + self.nphot05
            self.nphot2=nphot[2][:] + nphot[3][:] + self.nphot1
        
            self.nspec05=nspec[0][:]
            self.nspec1=nspec[1][:] + self.nspec05
            self.nspec2=nspec[2][:] + nspec[3][:] + self.nspec1

            for k in range(i):
                try:
                    self.compl05[k]=self.nspec05[k]/self.nphot05[k]
                except ZeroDivisionError:
                    self.compl05[k]=1.
                try:
                    self.compl1[k]=self.nspec1[k]/self.nphot1[k]
                except ZeroDivisionError:
                    self.compl1[k]=1.
                try:
                    self.compl2[k]=self.nspec2[k]/self.nphot2[k]
                except ZeroDivisionError:
                    self.compl2[k]=1.
            outfile='sdsscompleteness'+str(ntry)+'.dat'
            complout=open(outfile,'w')
            for i in range(nclusters):
                s=" "
                outstring=str(self.nspec05[i])+s+str(self.nphot05[i])+s+str(self.compl05[i])+s+str(self.nspec1[i])+s+str(self.nphot1[i])+s+str(self.compl1[i])+s+str(self.nspec2[i])+s+str(self.nphot2[i])+s+str(self.compl2[i])+" \n" 
                complout.write(outstring)
            complout.close()
            print "Error after completing ",nclusters," clusters - bummer!!!"
            print "Output thus far written to ",outfile

        print "average completeness of fiber sampling w/in 0.5xR200 = %5.2f +/- %5.2f"%(N.average(self.compl05),pylab.std(self.compl05))
        print "average completeness of fiber sampling w/in 1.xR200 = %5.2f +/- %5.2f"%(N.average(self.compl1),pylab.std(self.compl1))
        print "average completeness of fiber sampling w/in 2xR200 = %5.2f +/- %5.2f"%(N.average(self.compl2),pylab.std(self.compl2))
    def readsdsscompleteness(self,file):
        self.nspec05=N.zeros(len(self.z),'f')
        self.nphot05=N.zeros(len(self.z),'f')
        self.compl05=N.zeros(len(self.z),'f')
        self.nspec1= N.zeros(len(self.z),'f')
        self.nphot1= N.zeros(len(self.z),'f')
        self.compl1= N.zeros(len(self.z),'f')
        self.nspec2= N.zeros(len(self.z),'f')
        self.nphot2= N.zeros(len(self.z),'f')
        self.compl2= N.zeros(len(self.z),'f')
        complout=open(file,'r')
        i=0
        for line in complout:
            f=line.split()
            #print line
            #print f
            j=0
            for j in range(len(f)):
                f[j]=float(f[j])
            #print f
            (self.nspec05[i],self.nphot05[i],self.compl05[i],self.nspec1[i],self.nphot1[i],self.compl1[i],self.nspec2[i],self.nphot2[i],self.compl2[i])=(f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8])
            #print f[0],f[1],f[2],self.nspec05[i],self.nphot05[i],self.compl05[i]
            #print line
            i=i+1
        complout.close()
        print "got ",i," clusters from sdsscompleteness.dat"
        print "min of c.compl05 = ",min(self.compl05)
        #print self.compl05

    def getmemberids(self):
        self.memberids=N.zeros(len(self.mass),'f')
        self.memberids=list(self.memberids)
        self.sfmemberids=N.zeros(len(self.mass),'f')
        self.sfmemberids=list(self.sfmemberids)
        deltaid=0.1#tolerance for matching
        gindex=N.arange(0,len(g.clusterid),1)
        gindexsort=N.take(gindex,N.argsort(g.clusterid))
        gclusteridsort=N.take(g.clusterid,N.argsort(g.clusterid))
        for i in range(len(self.mass)):
            self.memberids[i]=[]
            (self.memberids[i],match)=findmatch(self.id[i],gclusteridsort,deltaid)
            sfids=[]
            for j in self.memberids[i]:
		try:
		    k=int(gindexsort[int(j)])
		except:
		    print "Error in getmemberids"
		    print "i, j, len(gindexsort) = ",i,j,len(gindexsort)
                if (g.ew[k] > ewmin):
                    if (g.lHa[k] > lHamin):
                        sfids.append(j)
	    self.sfmemberids[i]=sfids

    def getngal(self):
	self.ngal05=ngalmemb(2.,0.5)
	self.ngal1=ngalmemb(2.,1.)
	self.ngal2=ngalmemb(2.,2.)
	self.ngal05=N.array(self.ngal05,'f')
	self.ngal1=N.array(self.ngal1,'f')
	self.ngal2=N.array(self.ngal2,'f')
	if cdir.find('field') > -1:
	    self.ngal05 = (self.ngal05)/5.
	    self.ngal1 = (self.ngal1)/5.
	    self.ngal2 = (self.ngal2)/5.
        if (completenesscorr > 0):
            for i in range(len(self.ngal05)):
                try:
                    self.ngal05[i]=self.ngal05[i]/self.compl05[i]
                except ZeroDivisionError:
                    self.ngal05[i]=self.ngal05[i]
		    print "Warning:  ZeroDivisionError in compl corr ngal05",i
                try:
                    self.ngal1[i]=self.ngal1[i]/self.compl1[i]
                except ZeroDivisionError:
                    self.ngal1[i]=self.ngal1[i]
		    print "Warning:  ZeroDivisionError in compl corr ngal1",i
                try:
                    self.ngal2[i]=self.ngal2[i]/self.compl2[i]
                except ZeroDivisionError:
                    self.ngal2[i]=self.ngal2[i]
		    print "Warning:  ZeroDivisionError in compl corr ngal2",i

    def normalizefield(self):#divide summed quantities by number of control fields per cluster
	nfield=5.#number of control fields per cluster
	self.stellmass05 = self.stellmass05/nfield
	self.stellmass1 = self.stellmass1/nfield
	self.stellmass2 = self.stellmass2/nfield
	self.stellmass3 = self.stellmass3/nfield
	self.sfr05 = self.sfr05/nfield
	self.sfr1 = self.sfr1/nfield
	self.sfr2 = self.sfr2/nfield
	self.sfr3 = self.sfr3/nfield
	self.stellmass05bin = self.stellmass05bin/nfield
	self.stellmass1bin = self.stellmass1bin/nfield
	self.stellmass2bin = self.stellmass2bin/nfield
	self.stellmass3bin = self.stellmass3bin/nfield
	self.sfr05bin = self.sfr05bin/nfield
	self.sfr1bin = self.sfr1bin/nfield
	self.sfr2bin = self.sfr2bin/nfield
	self.sfr3bin = self.sfr3bin/nfield
	self.totlum05 = self.totlum05/nfield
	self.totlum1 = self.totlum1/nfield
	self.totlum2 = self.totlum2/nfield
	self.sfr05hiz = self.sfr05hiz/nfield



#"""
#    def getmemberids(self):
#        self.memberids=N.zeros(len(self.mass),'f')
#        self.memberids=list(self.memberids)
#        self.sfmemberids=N.zeros(len(self.mass),'f')
#        self.sfmemberids=list(self.sfmemberids)
#        deltaid=0.1#tolerance for matching
#        gindex=N.arange(0,len(g.clusterid),1)
#        gindexsort=N.take(gindex,N.argsort(g.clusterid))
#        gclusteridsort=N.take(g.clusterid,N.argsort(g.clusterid))
#	n05=0.
#	n1=0.
#	n2=0.
#        for i in range(len(self.mass)):
#            self.memberids[i]=[]
#            (self.memberids[i],flag)=findmatch(self.id[i],gclusteridsort,deltaid)
#            sfids=[]
#            for j in self.memberids[i]:
#                k=int(gindexsort[int(j)])
#                if (g.ew[k] > ewmin):
#                    if (g.lHa[k] > lHamin):
#                        sfids.append(j)
#            self.sfmemberids[i]=sfids
#"""


class Galaxy:
    def __init__(self):
        self.clusterid = []
        self.clusterrvir = []
        self.clusterr200 = []
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
        self.ew=[]
        self.agn1=[]
        self.agn2=[]
        self.agn3=[]
        self.agn4=[]
        self.agn5=[]        
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
                    virial=float(fields[4])#/h100
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
                self.agn1.append(float(fields[14]))#Miller AGN flag 0
                self.agn2.append(float(fields[15]))#Miller AGN flag 1
                self.agn3.append(float(fields[16]))#Miller AGN flag 2
                self.agn4.append(float(fields[17]))#Miller AGN flag 2
                self.agn5.append(float(fields[18]))#Miller AGN flag 2
                self.ew.append(float(fields[19]))#rest-frame EW(Halpha)
    def greadfieldfiles(self,clusters):
        for i in range(len(clusters)):
            prefix=clusters[i].split('.')
            t=prefix[0].split('_')
	    cid=int(str(t[2])+str(t[3]))
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
                    virial=float(fields[4])#/h100
                    sigma=float(fields[5])
                    z=float(fields[3])            
                    continue
                if line.find('*') > -1: #skip any galaxies w/ ****** in one or more fields
                    continue
                fields=line.split()
                self.clusterid.append(cid)#keep track of cluster id
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
                self.agn1.append(float(fields[14]))#Miller AGN flag 0
                self.agn2.append(float(fields[15]))#Miller AGN flag 1
                self.agn3.append(float(fields[16]))#Miller AGN flag 2
                self.agn4.append(float(fields[17]))#Miller AGN flag 2
                self.agn5.append(float(fields[18]))#Miller AGN flag 2
                self.ew.append(float(fields[19]))#rest-frame EW(Halpha)

    def convarray(self):
        self.clusterid = N.array(self.clusterid,'i')
        self.clusterrvir = N.array(self.clusterrvir,'f')
        self.clustersigma = N.array(self.clustersigma,'f')
        self.clusterz = N.array(self.clusterz,'f')
        self.clusterr200 = 1.73*self.clustersigma/1000.*1./N.sqrt(.7+.3*(1+self.clusterz)**3)/h100
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
        self.ew=N.array(self.ew,'f')
        self.agn1=N.array(self.agn1,'f')
        self.agn2=N.array(self.agn2,'f')
        self.agn3=N.array(self.agn3,'f')

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
        #self.lR=N.zeros(len(self.lHa),'f')#rband luminosity
        #self.ew=N.zeros(len(self.lHa),'f')
        #self.errew=N.zeros(len(self.lHa),'f')
        self.agn1=N.zeros(len(self.lHa),'f')#agn flags from Miller et al 2003
        self.agn2=N.zeros(len(self.lHa),'f')
        self.agn3=N.zeros(len(self.lHa),'f')
        self.sdsscatindex=N.zeros(len(self.lHa),'i')
        file = open('mycat','r')
        i=0
        for line in file:
            fields=line.split()
            #print "mycat line = ",i,len(self.ew)
            #self.ew[i]=float(fields[3])
            #self.errew[i]=float(fields[4])
            self.agn1[i]=float(fields[5])
            self.agn2[i]=float(fields[6])
            self.agn3[i]=float(fields[7])
            self.sdsscatindex[i]=int(fields[8])
            i += 1
    def getDn(self):
        self.Dn=N.zeros(len(self.dv),'f')
        self.errDn=N.zeros(len(self.dv),'f')
        dn=[]
        errdn=[]
        file=open("/Users/rfinn/SDSS/cluster_catalogs/dr2/sdss_dr2_deriv_4000Abreak.dat",'r')
        i=0
        for line in file:
            if line.find('#') > -1:
                continue
            fields=line.split()
            dn.append(float(fields[3]))
            errdn.append(float(fields[4]))
        file.close()

        dn=N.array(dn,'f')
        errdn=N.array(dn,'f')
        for i in range(len(self.dv)):
            self.Dn[i]=dn[self.sdsscatindex[i]]
            self.errDn[i]=errdn[self.sdsscatindex[i]]
            
    def calcagn(self): #set agn flag to 0 or 1 according to line ratios
        self.agn=N.zeros(len(self.lHa),'f')
        for i in range(len(self.agn)):
            #if (self.agn1[i] < 1) | ((self.agn1[i] > 4) & (self.agn[1] < 9)):
            #    self.agn[i]=0
            #else:
            #    self.agn[i]=1
            if (self.O3Hb[i] > -990.):
                self.agn[i]=1.
    def calcrlum(self):
        self.rlum=N.zeros(len(self.Mabs),'f')
        #self.rlum=10.**(-0.4*(self.Mabs-4.62))#R-band Lum in units of L-solar
        try:
            self.rlum=10.**(-0.4*(self.Mabs-4.62))#R-band Lum in units of L-solar
        except OverflowError:
            #print self.Mabs
            print "len(self.Mabs) = ",len(self.Mabs)
            for i in range(100):
                print self.Mabs[i]
            self.rlum=self.Mabs-4.62#R-band Lum in units of L-solar
                       
        #for i in range(len(self.rlum)):
        #    try:
        #        self.rlum[i]=10.**(-0.4*(self.Mabs[i]-4.62))
        #    except:
        #        self.rlum[i]=1.
        #self.rlum=10.**(-0.4*(self.Mabs+24))#R-band Lum in units of L-solar
                       
    def calcewr(self):#rest-frame EW
        self.ewr=self.ew#N.zeros(len(self.ew),'f')
        #self.ewr=self.ew/(1.+self.z)
    def calcmemb(self):#assign membership flag based on position in dv/sigma-dr/rvir plane
        self.memb=N.zeros(len(self.dv),'f')
        for i in range(len(self.memb)):
            dv=self.dv[i]/self.clustersigma[i]
            dr=self.dr[i]/self.clusterrvir[i]
            #y1=3.-(2./1.2)*dr
            #y2=-1*y1
            #y1=dr-3.
            #y2=3.-dr
            #if (dr < 1.2):
            #    if ((dv < y1) & (dv > y2)):
            #        #if ((dv < nsig) & (dr < nr)):
            #        self.memb[i]=1
            #if (dr >= 1.2):
            #    if (dv < 1.) & (dr < 2):
            #        self.memb[i]=1
            if ((dr < 1.) & (dv < nsig)):
                self.memb[i]=1
                continue
            if (dv < 1):
                self.memb[i]=1
    def limitmemb(self):
        self.clusterid =   N.compress(self.memb > 0,self.clusterid)
        self.clusterrvir = N.compress(self.memb > 0,self.clusterrvir)
        self.clusterr200 = N.compress(self.memb > 0,self.clusterr200)
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
        #self.lR   =        N.compress(self.memb > 0,self.lR)
        self.ew=           N.compress(self.memb > 0,self.ew)
        #self.errew=        N.compress(self.memb > 0,self.errew)
        self.agn1=         N.compress(self.memb > 0,self.agn1)
        self.agn2=         N.compress(self.memb > 0,self.agn2)
        self.agn3=         N.compress(self.memb > 0,self.agn3)

    def gdocalcs(self):
        self.convarray()
        #self.calcew()
        self.calcmemb()
        #self.limitmemb()#limit to member galaxies only

        self.calcagn()
        #self.getDn()
        #self.calcrlum()
        self.calcewr()
        self.calcextinction()
        self.convlHa2sfr()
        self.calcmyapcor()
        self.calcforplots()
        self.sfr=N.clip(self.sfr,0,100000)

    def trimclustersold(self):
        self.trim=N.zeros(len(self.clusterid),'f')
        i=0
        for line in self.clusterid:
            line=str(line)
            for name in c.id:
                if line.find(str(name)) > -1:
                    self.trim[i]=1.
                    continue
            i=i+1
        self.clusterid =   N.compress(self.trim > 0,self.clusterid)
        self.clusterrvir = N.compress(self.trim > 0,self.clusterrvir)
        self.clusterr200 = N.compress(self.trim > 0,self.clusterr200)
        self.clustersigma =N.compress(self.trim > 0,self.clustersigma)
        self.clusterz =    N.compress(self.trim > 0,self.clusterz)
        self.ra =          N.compress(self.trim > 0,self.ra)
        self.dec  =        N.compress(self.trim > 0,self.dec)
        self.z =           N.compress(self.trim > 0,self.z)
        self.dr=           N.compress(self.trim > 0,self.dr)
        self.dv=           N.compress(self.trim > 0,self.dv)
        self.dLMpc=        N.compress(self.trim > 0,self.dLMpc)
        self.fHa=          N.compress(self.trim > 0,self.fHa)
        self.lHa=          N.compress(self.trim > 0,self.lHa)
        self.apcor=        N.compress(self.trim > 0,self.apcor)
        self.stellarmass=  N.compress(self.trim > 0,self.stellarmass)
        self.az=           N.compress(self.trim > 0,self.az)
        self.O3Hb=         N.compress(self.trim > 0,self.O3Hb)
        self.N2Ha=         N.compress(self.trim > 0,self.N2Ha)
        #self.lR   =        N.compress(self.trim > 0,self.lR)
        self.ew=           N.compress(self.trim > 0,self.ew)
        #self.errew=        N.compress(self.trim > 0,self.errew)
        self.ewr=           N.compress(self.trim > 0,self.ewr)
        self.agn1=         N.compress(self.trim > 0,self.agn1)
        self.agn2=         N.compress(self.trim > 0,self.agn2)
        self.agn3=         N.compress(self.trim > 0,self.agn3)
        self.agn=         N.compress(self.trim > 0,self.agn)
        self.memb=         N.compress(self.trim > 0,self.memb)
        self.sfr=         N.compress(self.trim > 0,self.sfr)
        #self.Dn=         N.compress(self.trim > 0,self.Dn)
        self.Mabs=         N.compress(self.trim > 0,self.Mabs)

    def trimclusters(self):
        self.trim=N.zeros(len(self.clusterid),'f')
        i=0
        y=N.arange(0,len(self.clusterid),1)
        gindex=N.take(y,N.argsort(self.clusterid))#sort y according to x rankings
        gclusterids=N.take(self.clusterid,N.argsort(self.clusterid))
        deltaid=.01
        for name in c.id:
            (match,flag)=findmatch(name,gclusterids,deltaid)#
            try:
                for i in match:
                    j=gindex[int(i)]
                    self.trim[j]=1.
            except TypeError:
                continue
        self.clusterid =   N.compress(self.trim > 0,self.clusterid)
        self.clusterrvir = N.compress(self.trim > 0,self.clusterrvir)
        self.clusterr200 = N.compress(self.trim > 0,self.clusterr200)
        self.clustersigma =N.compress(self.trim > 0,self.clustersigma)
        self.clusterz =    N.compress(self.trim > 0,self.clusterz)
        self.ra =          N.compress(self.trim > 0,self.ra)
        self.dec  =        N.compress(self.trim > 0,self.dec)
        self.z =           N.compress(self.trim > 0,self.z)
        self.dr=           N.compress(self.trim > 0,self.dr)
        self.dv=           N.compress(self.trim > 0,self.dv)
        self.dLMpc=        N.compress(self.trim > 0,self.dLMpc)
        self.fHa=          N.compress(self.trim > 0,self.fHa)
        self.lHa=          N.compress(self.trim > 0,self.lHa)
        self.apcor=        N.compress(self.trim > 0,self.apcor)
        self.stellarmass=  N.compress(self.trim > 0,self.stellarmass)
        self.az=           N.compress(self.trim > 0,self.az)
        self.O3Hb=         N.compress(self.trim > 0,self.O3Hb)
        self.N2Ha=         N.compress(self.trim > 0,self.N2Ha)
        #self.lR   =        N.compress(self.trim > 0,self.lR)
        self.ew=           N.compress(self.trim > 0,self.ew)
        #self.errew=        N.compress(self.trim > 0,self.errew)
        self.ewr=           N.compress(self.trim > 0,self.ewr)
        self.agn1=         N.compress(self.trim > 0,self.agn1)
        self.agn2=         N.compress(self.trim > 0,self.agn2)
        self.agn3=         N.compress(self.trim > 0,self.agn3)
        self.agn=         N.compress(self.trim > 0,self.agn)
        self.memb=         N.compress(self.trim > 0,self.memb)
        self.sfr=         N.compress(self.trim > 0,self.sfr)
        #self.Dn=         N.compress(self.trim > 0,self.Dn)
        self.Mabs=         N.compress(self.trim > 0,self.Mabs)

    def trimgalaxies(self): #trim galaxies fainter than mabscut
        self.clusterid =   N.compress(self.Mabs < mabscut,self.clusterid)
        self.clusterrvir = N.compress(self.Mabs < mabscut,self.clusterrvir)
        self.clusterr200 = N.compress(self.Mabs < mabscut,self.clusterr200)
        self.clustersigma =N.compress(self.Mabs < mabscut,self.clustersigma)
        self.clusterz =    N.compress(self.Mabs < mabscut,self.clusterz)
        self.ra =          N.compress(self.Mabs < mabscut,self.ra)
        self.dec  =        N.compress(self.Mabs < mabscut,self.dec)
        self.z =           N.compress(self.Mabs < mabscut,self.z)
        self.dr=           N.compress(self.Mabs < mabscut,self.dr)
        self.dv=           N.compress(self.Mabs < mabscut,self.dv)
        self.dLMpc=        N.compress(self.Mabs < mabscut,self.dLMpc)
        self.fHa=          N.compress(self.Mabs < mabscut,self.fHa)
        self.lHa=          N.compress(self.Mabs < mabscut,self.lHa)
        self.apcor=        N.compress(self.Mabs < mabscut,self.apcor)
        self.stellarmass=  N.compress(self.Mabs < mabscut,self.stellarmass)
        self.az=           N.compress(self.Mabs < mabscut,self.az)
        self.O3Hb=         N.compress(self.Mabs < mabscut,self.O3Hb)
        self.N2Ha=         N.compress(self.Mabs < mabscut,self.N2Ha)
        #self.lR   =        N.compress(self.Mabs < mabscut,self.lR)
        self.ew=           N.compress(self.Mabs < mabscut,self.ew)
        #self.errew=        N.compress(self.Mabs < mabscut,self.errew)
        self.ewr=           N.compress(self.Mabs < mabscut,self.ewr)
        self.agn1=         N.compress(self.Mabs < mabscut,self.agn1)
        self.agn2=         N.compress(self.Mabs < mabscut,self.agn2)
        self.agn3=         N.compress(self.Mabs < mabscut,self.agn3)
        self.agn=         N.compress(self.Mabs < mabscut,self.agn)
        self.memb=         N.compress(self.Mabs < mabscut,self.memb)
        self.sfr=         N.compress(self.Mabs < mabscut,self.sfr)
        #self.Dn=         N.compress(self.Mabs < mabscut,self.Dn)
        self.Mabs=         N.compress(self.Mabs < mabscut,self.Mabs)

    def assigngalsub(self):#assign cluster substructure flag to members
        self.galsub=N.ones(len(self.z),'f')
        for i in range(len(c.id)):
            mids=list(c.memberids[i])
            for j in range(len(mids)):
                k=int(mids[j])
                self.galsub[k]=c.sub[i]
    def readdatafile(self,file):
        output=open(file,'r')
        i=-1
        for line in output:
            if (i == -1):
                fields=line.split()
                ncl=int(fields[1])#number of clusters
                print "number of galaxies = ",ncl
                self.clusterid=N.zeros(ncl,'f')
                self.clusterz=N.zeros(ncl,'f')
                self.clusterrvir=N.zeros(ncl,'f')
                self.clusterr200=N.zeros(ncl,'f')
                self.clustersigma=N.zeros(ncl,'f')
                self.dr=N.zeros(ncl,'f')
                self.dv=N.zeros(ncl,'f')
                self.ew=N.zeros(ncl,'f')
                self.lHa=N.zeros(ncl,'f')
                self.sfr=N.zeros(ncl,'f')
                self.stellarmass=N.zeros(ncl,'f')
                self.agn=N.zeros(ncl,'f')
                self.memb=N.zeros(ncl,'f')
                self.Mabs=N.zeros(ncl,'f')
                self.ra=N.zeros(ncl,'f')
                self.dec=N.zeros(ncl,'f')
                self.sigma10=N.zeros(ncl,'f')
                i += 1
                continue
            if line.find('#') > -1:#skip header
                continue
            f=line.split()
            for j in range(len(f)):
                f[j]=float(f[j])
            #(self.clusterid[i],self.clusterz[i],self.clusterrvir[i],self.clusterr200[i],self.clustersigma[i],self.dr[i],self.dv[i],self.ew[i],self.lHa[i],self.sfr[i],self.stellarmass[i],self.agn[i],self.memb[i],self.Mabs[i],self.ra[i],self.dec[i],self.sigma10[i])=(f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16])
            #print i, len(f)
            (self.clusterid[i],self.clusterz[i],self.clusterrvir[i],self.clusterr200[i],self.clustersigma[i],self.dr[i],self.dv[i],self.ew[i],self.lHa[i],self.sfr[i],self.stellarmass[i],self.agn[i],self.memb[i],self.Mabs[i],self.ra[i],self.dec[i])=(f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15])
            
            i += 1
        output.close()
    def getsigma10(self):
        mabs=mabscut
        self.sigma10=N.zeros(len(self.ew),'f')
        self.sigma10flag=N.zeros(len(self.ew),'f')
        for i in range(len(c.mass)):
            da = DA(c.z[i],h100)/1000.#Mpc/arcsec

            for j in c.memberids[i]:
                d=[]
                if (self.Mabs[j] < mabs):
                    self.sigma10flag[j] = 1.
                    for k in c.memberids[i]:
                        if (self.Mabs[k] < mabs):
                            d.append(N.sqrt((self.ra[j]-self.ra[k])**2+(self.dec[j]-self.dec[k])**2))
                    d=N.array(d,'f')
                    dsort=N.take(d,N.argsort(d))
                    self.sigma10[j]=10./(N.pi)/(dsort[10]*3600.*da)**2



class HizClusters:
    def __init__(self):
        self.z=N.array([.704,.748,.794,.845,.228,.32,.183,0.39],'f')
        self.sigma=N.array([418,504,1018,415,1023,1390,1274,561],'f')
        self.sumsfr=N.array([30.9,13.1,156.9,38.2,79.9,21.8,40.5,124.],'f')
        (self.sumsfr,self.sumsfrerr)=getlitsfr(self.z,self.sigma)
        #self.sumsfrerr=N.array([9.3,3.9,47.1,5.,8.5,14.1,18.9,17.],'f')
        #self.sumsfr=self.sumsfr/(h100)**2
        #self.sumsfrerr=self.sumsfrerr/(h100)**2
        #correct for AGN
        #self.sumsfr=self.sumsfr*.69
        #self.sumsfrerr=self.sumsfrerr*.69

        self.sfrcor=N.array([1.,1,1,1,1.7,2.8,2.8,1.],'f')#aperture correction
        self.volcor=N.array([1.,(1./.94),(1./.86),1,1,(1./.99),(1./.74),1.],'f')#aperture correction
        self.sfrcor=N.array([1.,(1.0),(1.16),1.,1.7,2.83,3.78,1.],'f')#aperture correction
        self.sumsfr=self.sumsfr*self.sfrcor
        self.sumsfrerr=self.sumsfrerr*self.sfrcor
        self.mass=12*(self.sigma/1000.)**3.*1./N.sqrt(.7+.3*(1+self.z)**3)/h100#x10^14 Mo
        self.kodamamasshiz=self.mass
        self.sfrmhiz=(self.sumsfr/self.mass)#correct for AGN contamination
        self.clname = ["\ca","\cb","\cc","\cj", "Abell~2390","AC~114","Abell~1689","CL0024.0$+$1652"]
        self.refnumber = ["1","1","1","2", "3","4","5","6"]
        self.technique = ["I","I","I","I", "I","S","S","I"]
        self.symbols = [-4,-4,-4,-3,6,4,7,11]
        self.refsymb = [-4,-5,-6,-3,6,4,7,11,18]
        self.colors =  [ 2, 2, 2, 2,2,2,2, 2, 1]
        self.symbols = [-4,-5,-6,-3,6,4,7,11,18]
        self.ref = ["Finn et al. 2005","Finn et al. 2004","Balogh & Morris 2000","Couch et al. 2001", "Balogh et al. 2002","Kodama et al. 2004","Alternate M\dcl\u"]
        self.ref = ["CL1040-1155","CL1054-1245","CL1216-1201","CLJ0023+0423B","Abell 2390","AC 114", "Abell 1689","CL 0024.0+1652","C4"]
        self.cosm = 1./N.sqrt(omegaL + omega0*(1.+self.z)**3)
            
def gotoit():
    #c=Cluster()
    #g=Galaxy()
    clusters = glob.glob("*_DR*")
    #if cdir.find('field') > -1:
#	clusters=glob.glob("*DR*_1.dat")#read in only first of each field file to establish "cluster" list
    print "reading in cluster files to get cluster parameters"
    c.creadfiles(clusters)        
#    if cdir.find('balogh') > -1:
#	c.creadfiles(clusters)        
#    if cdir.find('field') > -1:
#	c.creadfieldfiles(clusters)        
    clusters = glob.glob("*_DR*")#read in galaxies from all 5 field files
    print "reading in cluster files to get galaxies parameters"
    g.greadfiles(clusters)        
#    if cdir.find('balogh') > -1:
#	g.greadfiles(clusters)        
#    if cdir.find('field') > -1:
#	g.greadfieldfiles(clusters)        

    print "starting gdocalcs()"
    g.gdocalcs()
    c.convarray()
    print "starting limitz"
    c.limitz()
    print "starting limitsigma"
    c.limitsigma()
    #print "starting limitrichness"
    #c.limitrichness()
    print "starting assignsub"
    c.assignsub() #assign substructure flag

    print "plotting velocity histograms"
    if cdir.find('balogh') > -1:
        plotvelhist()
    else:
        c.recordchisqfield()
    #c.recordchisqfield()
        
    #print "starting supercluster"
    #c.supercluster() #assign substructure flag
    #print "starting limitsubtructure"

    print "assigning gaussian flag"
    c.assigngauss()
    print "plotting sigma vs totrlum"
    psplotinit("sigmatotrlum.ps")    #plot cluster sigma vs z 
    plotsigmatotrlum()
    ppgplot.pgend()
    if cdir.find('balogh') > -1:
        c.limitsubstructure() #cuts sample based on sub flag
    #c.limitsubstructure() #cuts sample based on sub flag
    print "number of clusters after substructure cut = ",len(c.z)
    print "number of clusters w/gaussian velocities = ",sum(c.gauss)
    print "assigning luminosity bins"
    c.assignlumbins()
    #print "Sorting cluster id array"
    #cindex=N.arange(0,len(c.id),1)
    #cindexsort=N.take(cindex,N.argsort(c.id))
    #cidsort=N.take(c.id,N.argsort(c.id))
    #deltaid=0.1
    print "matching x-ray clusters"
    c.getxray()
    #c.cdocalcs()
    print "trimming galaxy arrays"
    g.trimclusters()#keep galaxies in surviving clusters
    g.calcrlum()
    #print "trimming galaxy arrays for M > mabscut"
    #g.trimgalaxies()#keep galaxies in surviving clusters
    print "getting member ids for each cluster"
    c.getmemberids()#get indices in galaxy array of members for each cluster
    print "assigning galaxies sub flag"
    g.assigngalsub()#assign cluster sub flag to individual galaxies
    print "starting calcmass"
    c.calcmass()
    print "starting calcr200"
    c.calcr200()
    t=c.r200/c.rvir
    #print "R200/Rvir = ",t
    print "ave R200/rvir = ",N.average(t)
    if completenesscorr > 0:
	if (getcompl > 0):
	    print "starting getsdsscompleteness"
	    c.getsdsscompleteness()
	if (getcompl < 1):
	    c.readsdsscompleteness()
	print "Average completeness within 0.5xR200 = %6.3f +/- %6.3f" % (N.average(c.compl05),pylab.std(c.compl05))
	print "Average completeness within 1.0xR200 = %6.3f +/- %6.3f" % (N.average(c.compl1),pylab.std(c.compl1))
	print "Average completeness within 2.0xR200 = %6.3f +/- %6.3f" % (N.average(c.compl2),pylab.std(c.compl2))

    print "starting calcsfrstellmass"
    c.calcsfrstellmass()#calculate sfr05,sfr1,sfr2,stell05,etc
    print "starting calcrichness"
    c.calcrichness()
    print "starting calcsfrmass"
    #c.calcsfrmass()
    print "starting calcvolume"
    c.calcvolume()
    print "getting ngal for each cluster"
    c.getngal()#get indices in galaxy array of members for each cluster
    if cdir.find('field') > -1:
	c.normalizefield()        

    #g.getsigma10()
    #print "starting dostats"
    #dostats()
    #output=open("final-cluster-list",'w')
    #for name in c.id:
    #    output.write("%s \n" % (name))
def plotvelhist():#plot histogram of velocities with gaussian
    pylab.subplots_adjust(left=0.125, right=.9,bottom=.1,top=0.9,wspace=.3,hspace=.3)
    psplotinit("velhists.ps")
    #bins=N.arange(-6050.,6050.,100)
    bins=N.arange(-6.,6.,.1)
    nsig=2.#number to dv/sigma
    xmin=min(bins)
    xmax=max(bins)
    ymin=-1
    ymax=40
    allv=[]
    chisq=N.zeros(len(c.id),'f')
    for i in range(len(c.id)):
        dvel = N.compress((abs(g.clusterid - c.id[i]) < 1) & (g.dr < 3.), g.dv)
        sigma=c.sigma[i]
        dv2=dvel/sigma
        x1 = dv2.tolist()
        allv = allv + x1
        #print i,c.id[i],len(dvel)
        ppgplot.pgbox("",0.0,0,"",0.0)
        ppgplot.pgenv(min(bins),max(bins),ymin,ymax,0,0)
        z="%5.3f" % (c.z[i])
        r="%5.2f" % (c.rvir[i])
        s="%4.0f" % (c.sigma[i])
        sub="%1.0f" % (c.sub[i])
        rlum="%5.2f" % (c.totrlum[i])
        #print z,r
        #if abs(c.id[i] - 2110.) < 1:
        #    print "hey baby ",c.id[i],c.totrlum[i],c.totlum[i]
        cluster=str(c.id[i])+",N="+str(len(dvel))+",z="+str(z)+",Rv="+str(r)
        cluster2="\gs="+str(s)+",sub="+str(sub)+", Tot L\dR\u = "+str(rlum)
        ppgplot.pglab("\gDv/\gs","Ngal",cluster)
        ppgplot.pgtext(min(bins),42.,cluster2)
        #chisq[i]=velhistsub2(bins,dvel,sigma)
        chisq[i]=velhistsub2(bins,dv2,sigma,i)
        ppgplot.pgsci(4)
        ppgplot.pgsls(4)
        y=N.array([ymin,ymax])
        x=2.*N.ones(len(y),'f')
        ppgplot.pgline(x,y)
        x=-1.*x
        ppgplot.pgline(x,y)
        ppgplot.pgsci(1)
        ppgplot.pgsls(1)

        ppgplot.pgpage
    ppgplot.pgend()
    psplotinit("velhistsall.ps")
    bins=N.arange(-6.025,6.025,.05)
    xmin=min(bins)
    xmax=max(bins)
    ymin=-5
    ymax=500.
    ppgplot.pgbox("",0.0,0,"",0.0)
    ppgplot.pgenv(min(bins),max(bins),ymin,ymax,0,0)
    ppgplot.pglab("\gD v/\gs","Ngal","")
    #print "yikes ",allv
    allv=N.array(allv,'f')
    velhistsub3(bins,allv,1.)
    ppgplot.pgsci(4)
    ppgplot.pgsls(4)
    y=N.array([ymin,ymax])
    x=2.*N.ones(len(y),'f')
    ppgplot.pgline(x,y)
    x=-1.*x
    ppgplot.pgline(x,y)
    ppgplot.pgsci(1)
    ppgplot.pgsls(1)
    ppgplot.pgend()
    c.recordchisq(chisq)
def velhistsub(bins,x,sigma):
    xbinnumb=((x-min(bins))*len(bins)/(max(bins)-min(bins)))#calculate x  bin number for each point 
    y=N.zeros(len(bins),'f')
    for i in range(len(y)):
        for numb in xbinnumb:
            if int(numb) == i:
                y[i]=y[i]+1
    drawhist(bins,y)
    yfit=max(y)*N.exp(-0.5*(bins/sigma)**2)
    ppgplot.pgsci(2)
    ppgplot.pgline(bins,yfit)
    ppgplot.pgsci(1)

def velhistsub2(bins,x,sigma,k):
    xbinnumb=((x-min(bins))*len(bins)/(max(bins)-min(bins)))#calculate x  bin number for each point 
    y=N.zeros(len(bins),'f')
    sigma=1.#velocities already normalized by sigma
    for i in range(len(y)):
        for numb in xbinnumb:
            if int(numb) == i:
                y[i]=y[i]+1
    drawhist(bins,y)
    histsum=0.
    #dv=2.*sigma
    for i in range(len(bins)):
        if (bins[i] >= -2.) & (bins[i] <= 2.):

            histsum = histsum + y[i]
            #print "dude",i, bins[i], y[i], histsum

    #print "histsum = ",histsum,histsum/0.95,histsum/0.95*(bins[1]-bins[0]),histsum/N.sqrt(2.*N.pi*sigma**2)/(bins[1]-bins[0]),max(y)
    ymaxfit = histsum/0.95*(bins[1]-bins[0])#integral of gaussian is N.sqrt(N.pi)
    #area w/in +/-2 sigma = 95% of total area
    #so fit gaussian so that area within 2 sig = number of gal w/in 2 sig
    #yfit=max(y)*N.exp(-0.5*(bins/sigma)**2)
    yfit=ymaxfit/N.sqrt(2.*N.pi*sigma**2)*N.exp(-0.5*(bins/sigma)**2)
    ppgplot.pgsci(2)
    ppgplot.pgline(bins,yfit)
    ppgplot.pgsci(1)
    #calculate difference b/w observed and gaussian histograms
    #dv=abs(bins[1]-bins[0])
    deltav=abs(bins[1]-bins[0])
    gaussbins=bins+deltav#puts bins at center of each bin instead of left side
    ygauss=ymaxfit/N.sqrt(2.*N.pi*sigma**2)*N.exp(-0.5*(gaussbins/sigma)**2)
    #sum=N.sqrt(N.sum((ygauss-y)**2))/(1.*len(bins))
    #the eqn below blows up when there is discrepancy at large delta v
    #alternative is equally weighted errors
    #sum=N.sum((ygauss-y)**2/ygauss)/histsum
    #using equally weighted errors
    ysum=0.
    gysum=0.
    sum=0.
    nbin=0.
    diff=0.
    for i in range(len(bins)):
        if (bins[i] >= -2.):
            if (bins[i] <= 2.):
                if (ygauss[i] < 0):
                    #error=ygauss[i]+2
                    error=1.
                else:
                    error=ygauss[i]
                sum = sum + (ygauss[i]-y[i])**2/(error)
                #try:
                #    sum=sum + (ygauss[i]-y[i])**2/(ygauss[i])
                #except ZeroDivisionError:
                #    sum=sum + (ygauss[i]-y[i])**2
                diff=diff+abs(ygauss[i]-y[i])
                ysum=ysum+y[i]
                gysum=gysum+ygauss[i]
                nbin=nbin + 1.
                #print "%2d bin/sigma = %5.2f ygauss[i] = %4.2f yobs[i] = %4.2f (y-yobs)=%4.2f diff=%4.2f"%(i,(bins[i]/sigma),ygauss[i],y[i],(abs(ygauss[i]-y[i])),diff)
    #sum=sum/(ysum-2)#gary schmidt says gaussian is N-2 degrees of freedom
    sum=sum/(nbin-3)#dennis says use number of bins instead of number of galaxies, and gaussian has N - 3 degrees of freedom b/c fit mean, std dev, and normalization
    ppgplot.pgsci(2)
    defch=ppgplot.pgqch()
    ppgplot.pgsch(1.)
    xlabel=-5.
    s="\gx\u2\d\d\gn\u = %4.2f" %(sum)
    ppgplot.pgtext(xlabel,35.,str(s))
    s="\gSabs(yg-yobs) = %4.1f" %(diff)
    ppgplot.pgtext(xlabel,33.,str(s))
    s="Ntot = %4.1f" %(ysum)
    ppgplot.pgtext(xlabel,31.,str(s))
    s="Nexp = %4.1f" %(gysum)
    ppgplot.pgtext(xlabel,29.,str(s))
    s="\gSdiff/Nexp = %5.2f" %(diff/gysum)
    ppgplot.pgtext(xlabel,27.,str(s))
    ppgplot.pgsch(defch)

    if (k < 20):
	pylab.subplot(5,4,(k+1))

	pylab.plot(bins,yfit,'r',linewidth=1.5)
	pylab.hold(True)
	drawhistpylab(bins,y)
	#pylab.ylabel('Ngal')
	pylab.axis([-3.,3.,-1.,18.])
	pylab.yticks(N.arange(0.,20.,4.))
	pylab.hold(False)
	#pylab.show()
	s="%3.1f" %(sum)
	pylab.text(-2.,10.,r"$\chi^2_\nu  = $",fontsize=14) 
	pylab.text(0.,12.,s,fontsize=12) 
    if (k == 20):
	pylab.text(-12.,-17.,r"$\Delta v/\sigma$", horizontalalignment='center',fontsize=24)
	pylab.text(-30.,55.,r"$N_{gal}$", verticalalignment='center',rotation='vertical',fontsize=24)
	
	pylab.savefig('velhistpy.eps')

    return sum

def velhistsub3(bins,x,sigma):#same as velhist2 but no text on plot
    xbinnumb=((x-min(bins))*len(bins)/(max(bins)-min(bins)))#calculate x  bin number for each point 
    y=N.zeros(len(bins),'f')
    sigma=1.#velocities already normalized by sigma
    for i in range(len(y)):
        for numb in xbinnumb:
            if int(numb) == i:
                y[i]=y[i]+1
    drawhist(bins,y)
    histsum=0.
    #dv=2.*sigma
    for i in range(len(bins)):
        if (bins[i] >= -2.) & (bins[i] <= 2.):

            histsum = histsum + y[i]
            #print "dude",i, bins[i], y[i], histsum

    #print "histsum = ",histsum,histsum/0.95,histsum/0.95*(bins[1]-bins[0]),histsum/N.sqrt(2.*N.pi*sigma**2)/(bins[1]-bins[0]),max(y)
    ymaxfit = histsum/0.95*(bins[1]-bins[0])#integral of gaussian is N.sqrt(N.pi)
    #area w/in +/-2 sigma = 95% of total area
    #so fit gaussian so that area within 2 sig = number of gal w/in 2 sig
    #yfit=max(y)*N.exp(-0.5*(bins/sigma)**2)
    yfit=ymaxfit/N.sqrt(2.*N.pi*sigma**2)*N.exp(-0.5*(bins/sigma)**2)
    ppgplot.pgsci(2)
    ppgplot.pgline(bins,yfit)
    ppgplot.pgsci(1)
    #calculate difference b/w observed and gaussian histograms
    #dv=abs(bins[1]-bins[0])
    deltav=abs(bins[1]-bins[0])
    gaussbins=bins+deltav#puts bins at center of each bin instead of left side
    ygauss=ymaxfit/N.sqrt(2.*N.pi*sigma**2)*N.exp(-0.5*(gaussbins/sigma)**2)
    #sum=N.sqrt(N.sum((ygauss-y)**2))/(1.*len(bins))
    #the eqn below blows up when there is discrepancy at large delta v
    #alternative is equally weighted errors
    #sum=N.sum((ygauss-y)**2/ygauss)/histsum
    #using equally weighted errors
    ysum=0.
    gysum=0.
    sum=0.
    nbin=0.
    diff=0.
    for i in range(len(bins)):
        if (bins[i] >= -2.):
            if (bins[i] <= 2.):
                if (ygauss[i] < 0):
                    #error=ygauss[i]+2
                    error=1.
                else:
                    error=ygauss[i]
                sum = sum + (ygauss[i]-y[i])**2/(error)
                #try:
                #    sum=sum + (ygauss[i]-y[i])**2/(ygauss[i])
                #except ZeroDivisionError:
                #    sum=sum + (ygauss[i]-y[i])**2
                diff=diff+abs(ygauss[i]-y[i])
                ysum=ysum+y[i]
                gysum=gysum+ygauss[i]
                nbin=nbin + 1.
    return sum


def writefiles():
    """
    write out cluster and galaxy files
    cluster file  contains:
    id z rvir r200 sigma totrlum lx stellmassr05 stellmassr1 stellmassr2 sfr05 sfr1 sfr2 sffrac05 sffrac1 sffrac2

    galaxy
    clusterid clusterz clusterrvir clusterr200 dr dv ew lHa agnflag sfflag memb

    """
    output=open("mygalaxies.cat",'w')
    header="# "+str(len(g.z))+ " galaxies \n"
    output.write(header)
    header="#cid cz crvir cr200 csigma dr dv ew lHa stellmass agnflag sfflag memb Mabs"
    output.write("%s \n" % (header))
    for i in range(len(g.z)):
        #             cid   cz    crvir cr200 csigma dr   dv    ew    lHa   sfr   stelm agn   memb
        #output.write('%5.0f %5.4f %5.3f %5.3f %6.1f %7.3f %7.1f %6.2f %6.2f %6.2f %3.2e %2.0f %2.0f %5.2f %12.8f %12.8f %8.2f\n' %(g.clusterid[i],g.clusterz[i],g.clusterrvir[i],g.clusterr200[i],g.clustersigma[i],g.dr[i],g.dv[i],g.ew[i],g.lHa[i],g.sfr[i],g.stellarmass[i],g.agn[i],g.memb[i],g.Mabs[i],g.ra[i],g.dec[i],g.sigma10[i]))
        output.write('%5.0f %5.4f %5.3f %5.3f %6.1f %7.3f %7.1f %6.2f %6.2f %6.2f %3.2e %2.0f %2.0f %5.2f %12.8f %12.8f\n' %(g.clusterid[i],g.clusterz[i],g.clusterrvir[i],g.clusterr200[i],g.clustersigma[i],g.dr[i],g.dv[i],g.ew[i],g.lHa[i],g.sfr[i],g.stellarmass[i],g.agn[i],g.memb[i],g.Mabs[i],g.ra[i],g.dec[i]))
    output.close()

    output=open("myclusters.cat",'w')
    header="# "+str(len(c.z))+ " clusters \n"
    output.write(header)
    header="#id z rvir r200 sigma totrlum lx stellmass05 stellmass1 stellmass2 stellmass3 sfr05 sfr1 sfr2 sfr3 sffrac05 sffrac1 sffrac2 sub r200bin sigmabin sfr05bin sfr1bin sfr2bin sfr3bin stellmass05bin stellmass1bin stellmass2bin stellmass3bin chisq(diff from gauss) sfr05hiz kauffmann ra dec"
    output.write("%s \n" % (header))
    for i in range(len(c.z)):
        #print i,c.id[i],c.z[i],c.rvir[i],c.r200[i],c.sigma[i],c.Lx[i],c.xrayflag[i],c.stellmass05[i],c.stellmass1[i],c.stellmass2[i],c.stellmass3[i],c.sfr05[i],c.sfr1[i],c.sfr2[i],c.sfr3[i],c.sffrac05[i],c.sffrac1[i],c.sffrac2[i]
        #             id    z     rvir  r200  sigma totrl lx    xflg  ste05 stel1 stel2 stel3 sfr05 sfr1  sfr2  sfr3  sff05 sff1  sff2  sub   totl
        output.write('%5.0f %5.4f %5.3f %5.3f %6.1f %6.1f %6.3f %2.0f %3.2e %3.2e %3.2e %3.2e %5.1f %5.1f %5.1f %5.1f %4.2f %4.2f %4.2f %2.0f %4.3e %4.3e %4.3e %5.3f %5.3f %5.1f %5.1f %5.1f %5.1f %3.2e %3.2e %3.2e %3.2e %5.2f %3.2e %5.3f %12.8f %12.8f\n' %(c.id[i],c.z[i],c.rvir[i],c.r200[i],c.sigma[i],c.totrlum[i],c.Lx[i],c.xrayflag[i],c.stellmass05[i],c.stellmass1[i],c.stellmass2[i],c.stellmass3[i],c.sfr05[i],c.sfr1[i],c.sfr2[i],c.sfr3[i],c.sffrac05[i],c.sffrac1[i],c.sffrac2[i],c.sub[i],c.totlum05[i],c.totlum1[i],c.totlum2[i],c.r200bin[i],c.sigmabin[i],c.sfr05bin[i],c.sfr1bin[i],c.sfr2bin[i],c.sfr3bin[i],c.stellmass05bin[i],c.stellmass1bin[i],c.stellmass2bin[i],c.stellmass3bin[i],c.chisq[i],c.sfr05hiz[i],c.kauffmann[i],c.ra[i],c.dec[i]))
    output.close()
#Main
hiz=HizClusters()

print "hey, mode = ",mode
if (mode == 0):
    print "looks like mode < 1"
    gotoit()
    writefiles()
print "hey, mode = ",mode
if (mode == 1):
    print "looks like mode > 0"
    c=Cluster()
    c.readdatafile("/Users/rfinn/SDSS/fieldDR4/myclusters.cat")
    print "reading in control galaxy data file"
    #g.readdatafile("/Users/rfinn/SDSS/fieldDR4/mygalaxies.cat")

    f=c

    c=Cluster()
    print "reading in cluster data file"
    c.readdatafile("/Users/rfinn/SDSS/baloghDR4/myclusters.cat")
    print "reading in galaxy data file"
    #g.readdatafile("/Users/rfinn/SDSS/baloghDR4/mygalaxies.cat")
    print "got ",len(c.z)," clusters"
mil=sim()
mil.readdatafile()    
#plotsfrlx()
endtime=time.clock()
print "elapsed time = ",endtime-starttime
#plotbinneddata()
#plotewsigma10()
#plotsfrmabs()
#plotsfr(0)


plotsfrsigmapylab()
plotngalsigmapylab()
plotstellarsigmapylab()
plotsfrstellarsigmapylab()

#plotsfrsigma3sigmapylab()
#plotngalsigma3sigmapylab()
#plotstellarsigma3sigmapylab()



#plotsfrsigmaratiopylab()
endtime=time.clock()
print "end time = ",endtime
print "elapsed time = ",endtime-starttime
log.close()
