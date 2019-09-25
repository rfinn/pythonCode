#!/usr/bin/env python
import sys, os
import Numeric as N
#import scipy
from math import *
from mystuff import *
import ppgplot
import random

deflw=3
ewlabel = "EW H\ga+[NII] (\(2078))"
ewo2label = "EW [OII] (\(2078))"
sfrlabel = "SFR (h\d100\u\u-2\d M\d\(2281)\u yr\u-1\d)"
sfrlabel70 = "SFR (M\d\(2281)\u yr\u-1\d)"#label for ms1054, h100=.7
sumsfrlabel = "\gSSFR (h\d100\u\u-2\d M\d\(2281)\u yr\u-1\d)"
sfrmlabel = "\gSSFR/M\dcl\u (h\d100\u\u-3\d M\d\(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u)"
mcllabel = "M\dcl\u (h\d100\u\u-1\d 10\u14\d M\d\(2281)\u)"
omega0=0.3
omegaL=0.7
h100=.7
dlambdaJ=2500 #bandwidth of J filter
nstar=0
pscale=0.18
sfrmin=.5
sfrcomp=.5#sfr cut for comparision w/literature
#parameters for comparison w/literature
f200a=.5 #fraction of R200 to include in integrated SFR and femiss
f200b=2. #fraction of R200 to include in integrated SFR and femiss
lhamin=sfrmin*1.27#lhamin in units of 1d40 ergs/s
snmin=3. #min s/n cut for continuum-subtracted flux
ewmin=10.
ewcomp=20#ew cut for comparison w/low-z studies
dm=0.
ewburst=40 #ew cut for defining starburst
sfrcorr=2.5*0.77#sfr correction factors for dust and NII

def matchname(x1,x2):#match pair  x1,y1 with names in array x2, y2
	matchflag=0
	imatch=0
	for i in range(len(x2)):
		if x1.find(x2[i]) > -1:
		    imatch=i
		    matchflag=1
	return imatch, matchflag



def binitsum(x,y,xmin,xmax,magstep):#bin arrays x, y into n bins, returning xbin,ybin
    nx=len(x)
    #x=N.array(x,'f')
    #y=N.array(y,'f')
    nbin=int((xmax-xmin)/magstep)
    xbin=N.zeros(nbin,'f')
    ybin=N.zeros(nbin,'f')
    #ybinerr=N.zeros(n,'f')
    for i in range(nbin):
        nxmin=xmin+i*magstep
        nxmax=xmin+(i+1)*magstep
        #xbin[i]=(nxmin + nxmax)/2.
        xbin[i]=nxmin
        ybin[i]=N.sum(N.compress((x > nxmin) & (x < nxmax),y))#sum y
    return xbin, ybin
def binithist(x,xmin,xmax,magstep):#bin arrays x, y into n bins, returning xbin,ybin
    nx=len(x)
    nbin=int((xmax-xmin)/magstep)
    xbin=N.zeros(nbin,'f')
    ybin=N.zeros(nbin,'f')
    for i in range(nbin):
        nxmin=xmin+i*magstep
        nxmax=xmin+(i+1)*magstep
        #xbin[i]=(nxmin + nxmax)/2.
        xbin[i]=nxmin
        ybin[i]=float(len(N.compress((x > nxmin) & (x < nxmax),x)))#sum y
    return xbin, ybin

def checkposition(x,y):
    if (ncl == 3):#CLJ0023
        pos = checkfilter1(x,y)
    if (ncl == 2):#CL1216
        pos = checkfilter4(x,y)
    if (ncl == 1):#CL1054-12
        pos =checkfilter6(x,y)
    if (ncl == 0):#CL1040
        pos = checkfilter8(x,y)
    if (ncl == 4):#RXJ0152
        pos = checkfilter2(x,y)
    if (ncl == 5):#MS1054.4-0321
        pos = checkfilterms1054(x,y)
    if (ncl == 6):#MS1054.4-0321
        pos = checkfilterHDFN1(x,y)
    if (x < 1):
        pos = 0
    if (y < 1):
        pos = 0
    if (x > c.xmax[ncl]):
        pos = 0
    if (y > c.ymax[ncl]):
        pos = 0
    return pos

def checkfilter4(x,y):#CL1216
    pos = 1
    #check for bottom right corner
    if ((x > 439) & (y < 431)):
        d=sqrt((x-326)**2 + (y-468)**2)
        if(d > 500) :
            pos = 0
    #check for top left corner
    if ((x < 200) & (y > 590)):
        d=sqrt((x-452)**2 + (y-316)**2)
        if (d > 500):
            pos = 0
    #check for bottom left corner
    if ((x < 300) & (y < 300)):
        d=sqrt((x-400)**2 + (y-486)**2)
        if (d > 500):
            pos = 0
    #check for top right corner
    if ((x > 500) & (y > 500)):
        d=sqrt((x-366)**2 + (y-382)**2)
        if (d > 500):
            pos = 0
    #xystars = [(489.,757.),(75.53,454.1),(492.,386.),(170.,454.1),(181.,336.9),(347.6,615.9),(746.4,496.6),(520.5,499.),(703.1,573.8),(580.1,592.7),(347.6,615.9),(93.,760.),(444.14,139.14)]
    xystars=[]
    xystars=[(0.,0.)]
    for (xs,ys) in xystars:
        d = sqrt((xs-x)**2 + (ys-y)**2)
        #if d < 4:
        #    pos = 0
    return pos
def checkfilter6(x,y):
    pos = 1
    #check for top right and left corners
    if (y > 500):
	d=sqrt((x-430)**2 + (y-423)**2)
        if(d > 495):
            pos = 0
    #check for bottom left corner
    if ((y < 250)):
	d=sqrt((x-465)**2 + (y-518)**2)
        if (d > 500):
            pos = 0
    #check for bottom right corner
    if ((x > 514) & (y < 384)):
	d=sqrt((x-514)**2 + (y-384)**2)
        if(d > 400):
            pos = 0
    #bottom edge
    if (y < 10):
        pos = 0
    #reject stars
    #xystars = [(297.,84.),(827.,653.),(704.,286.),(715.,68.)]#ghost images
    xystars=[]
    xystars=[(0.,0.)]
    for (xs,ys) in xystars:
      d = sqrt((xs-x)**2 + (ys-y)**2)
      if d < 4:
          pos = 1



    return pos

def checkfilter8(x,y):
    #check for right side
    pos = 1
    if (y > 792.):
        pos = 0
    if ((x > 334)):
	d=sqrt((x-334)**2 + (y-464)**2)
        if(d > 440):
            pos = 0
    #check for larger inscribed circle
    d=sqrt((x-400)**2 + (y-400)**2)
    if (d > 430.):
        pos = 0
    #reject foreground galaxies that have ap problems & stars

    xystars = [(362.,631.)]#ghost image
    xystars=[(0.,0.)]
    for (xs,ys) in xystars:
        d = sqrt((xs-x)**2 + (ys-y)**2)
        if d < 4:
            pos = 0
    return pos

def checkfilter1(x,y):
    pos = 1
    d=sqrt((x-485)**2 + (y-460)**2)
    if (d > 480.):
	pos=0
    if (y > 810):
        pos = 0
    if (x > 800):
        pos = 0
    #reject stars
    xystars = [(389.27,444.69),(634.71,281.38),(415.,261.)]
    for (xs,ys) in xystars:
        d = sqrt((xs-x)**2 + (ys-y)**2)
        if d < 4:
            pos = 0
    return pos
def checkfilter2(x,y):
    pos = 1
    if (x > 450):
        d=sqrt((x-391)**2 + (y-470)**2)
        if ( d > 510):
            pos = 0

    if (x < 450) & (y < 450):
        d=sqrt((x-485)**2 + (y-528)**2)
        if ( d > 510):
            pos = 0

    if (y > 516):
        d=sqrt((x-481)**2 + (y-418)**2)
        if ( d > 510):
            pos = 0

    if (y < 20):
        pos = 0
    
    xystars = [(47.,541.),(697.,349.),(813.816,263.114),(840.,498.6)]
    for (xs,ys) in xystars:
      d = sqrt((xs-x)**2 + (ys-y)**2)
      if d < 4:
          pos = 0
    return pos

def checkfilterms1054(x,y):
    pos = 1
    d = sqrt((437.-x)**2 + (475.-y)**2)
    if (d > 400):
        pos = 0
    #xystars = [(666.,308.),(814.,458.),(420.,151.),(341.,756.)]
    #for (xs,ys) in xystars:
    #  d = sqrt((xs-x)**2 + (ys-y)**2)
    #  if d < 4:
    #      pos = 0
    return pos

def checkfilterHDFN1(x,y):
    pos = 1
    if (y > 880.):
        pos = 0
    d = sqrt((457.-x)**2 + (484.-y)**2)
    if (d > 440):
        pos = 0
    xystars = [(502.,306.),(458.,306.),(202.,382.),(240.,546.),(194.,346.),(286.,650.)]
    for (xs,ys) in xystars:
      d = sqrt((xs-x)**2 + (ys-y)**2)
      if d < 4:
          pos = 0
    return pos

def getdL(z):
  H0=h100*100
  c=3.*10**5
  nstep=1000
  dz=z/(nstep-1)
  dL=0
  i=0
  for i in range(nstep):
    zi=(i-1)*dz
    Ez=sqrt(omega0*(1+zi)**3+omegaL)
    #print "H0 z dz Ez omega0 omegaL \n"
    dL=dL+c/H0*(1+z)*dz/Ez
  dL=dL*10**6.*3.09*10**18
  return dL

def psplotinit(output):
    file=output+"/vcps"
    ppgplot.pgbeg(file,1,1)
    ppgplot.pgpap(8.,1.)
    ppgplot.pgsch(1.7) #font size
    ppgplot.pgslw(6)  #line width

def plothistsfr():
    DATAMIN = -4.
    DATAMAX = 15.
    NBIN = int((DATAMAX-DATAMIN)*2.)    
    #print "ngal = ",len(g0.sfr)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,0,45,0)
    ppgplot.pglab("SFR (h\d100\u\u-2\d M\d\(2281)\u yr\u-1 \d)","Number of Galaxies","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)
    x=N.compress((g0.final > 0),g0.sfrc)
    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    xlabel = 6.5
    ylabel = 38.
    ystep = 3.
    dy=.4
    dxl=3
    dxr=.5
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgline(xlin,ylin)

    ppgplot.pgslw(5)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    #ppgplot.pgsci(3)
    #x=N.compress((abs(g1.ew) > ewmin),g1.sfr)
    x=N.compress((g1.final > 0),g1.sfrc)
    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)

    ylabel = ylabel - ystep

    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgline(xlin,ylin)
    ppgplot.pgsls(1)
    ppgplot.pgslw(deflw)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    
    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    x=N.compress((g2.final > 0),g2.sfrc)
    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgslw(2)  #line width
    ppgplot.pgline(xlin,ylin)


    #print "Number in g2.ratios = ",len(g2.ratio)
    #for ratio in g2.ratio:
    #    print ratio
    #drawbinned(x,y,5)
    ppgplot.pgsci(1)
    #ppgplot.pgend()

def plothistew(ew):
    DATAMIN = -60.
    DATAMAX = 300.
    NBIN = int((DATAMAX-DATAMIN)/2.5)    
    #print "ngal = ",len(g0.sfr)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,0,60,0)
    ppgplot.pglab("EW H\ga+[NII] (\(2078))","Number of Galaxies","")#angstrom = 2078
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    x=ew
    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    y=N.arange(-11.,101.,2.)
    x=ewmin*N.ones(len(y),'f')
    ppgplot.pgsls(4)
    ppgplot.pgline(x,y)
    x=-1.*x
    ppgplot.pgline(x,y)
    ppgplot.pgsls(1)
    x=0.*x
    ppgplot.pgline(x,y)
def plothistewall():
    DATAMIN = -60
    DATAMAX = 300.
    NBIN = int((DATAMAX-DATAMIN)/10.)    
    #print "ngal = ",len(g0.sfr)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,0,45,0)
    ppgplot.pglab("EW H\ga+[NII] (\(2078))","Number of Galaxies","")#angstrom = 2078
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=g0.ew
    x=N.compress(g0.final > 0,g0.ewc)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)
    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    xlabel = 120.
    ylabel = 38.
    ystep = 4.
    ppgplot.pgslw(deflw)
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    ppgplot.pgslw(4)
    #ppgplot.pgsch(1)
    xlin = N.array([xlabel-50.,xlabel-5.],'f')
    ylin = N.array([ylabel+1.,ylabel+1.],'f')
    ppgplot.pgline(xlin,ylin)

    ppgplot.pgslw(3)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    #ppgplot.pgsci(3)
    #x=N.compress((abs(g1.ew) > ewmin),g1.sfr)
    #x=g1.ew
    x=N.compress(g1.final > 0,g1.ewc)
    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    xlin = N.array([xlabel-50.,xlabel-5.],'f')
    ylin = N.array([ylabel+1.,ylabel+1.],'f')
    ppgplot.pgline(xlin,ylin)

    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    #x=g2.ew
    x=N.compress(g2.final > 0,g2.ewc)
    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(3)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-50.,xlabel-5.],'f')
    ylin = N.array([ylabel+1.,ylabel+1.],'f')
    ppgplot.pgslw(2)  #line width
    ppgplot.pgline(xlin,ylin)

    #print "Number in g2.ratios = ",len(g2.ratio)
    #for ratio in g2.ratio:
    #    print ratio
    #drawbinned(x,y,5)

    ppgplot.pgsci(1)
    ppgplot.pgslw(1)
    ppgplot.pgsls(4)

    xlin = N.array([ewmin,ewmin],'f')
    ylin = N.array([-20.,200.],'f')
    ppgplot.pgline(xlin,ylin)
    xlin = -1.*xlin
    ppgplot.pgline(xlin,ylin)

    
    #ppgplot.pgend()

def plotnkall():
    DATAMIN = 15.
    DATAMAX = 24.
    NBIN = int((DATAMAX-DATAMIN)*2.)    
    #print "ngal = ",len(g0.sfr)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,0,15,0)
    ppgplot.pglab("K\ds\u (r < 2'')","Number of Galaxies","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)
    x=N.compress((g0.final > 0) & (g0.finalsf > 0),g0.vltk2)
    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    xlabel = 21.2
    ylabel = 13.
    ystep = 1.5
    dy=.4
    dxl=1
    dxr=.2
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgline(xlin,ylin)

    ppgplot.pgslw(3)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    #ppgplot.pgsci(3)
    #x=N.compress((g1.final > 0),g1.sfr)
    x=N.compress((g1.final > 0) & (g1.finalsf > 0),g1.vltk2)
    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)

    ylabel = ylabel - ystep

    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgline(xlin,ylin)
    ppgplot.pgsls(1)
    ppgplot.pgslw(deflw)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    
    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    x=N.compress((g2.final > 0) & (g2.finalsf > 0),g2.vltk2)

    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgslw(2)  #line width
    ppgplot.pgline(xlin,ylin)


    #print "Number in g2.ratios = ",len(g2.ratio)
    #for ratio in g2.ratio:
    #    print ratio
    #drawbinned(x,y,5)
    ppgplot.pgsci(1)
    #ppgplot.pgend()

def plotnkstarburst():
    DATAMIN = -24.
    DATAMAX = -15.
    ewstarburst = 40.
    mvcut=-20.52
    NBIN = int((DATAMAX-DATAMIN)*2.)    
    #print "ngal = ",len(g0.sfr)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,0,9,0)
    ppgplot.pglab("M\dV\u + 5log(h\d100\u)","Number of Galaxies","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)
    x=N.compress((g0.final > 0) & (g0.finalsf > 0) & (g0.ewc > ewstarburst) & (g0.dkpc < (0.5*c.r200Mpc[0]*1000.)),g0.Mv)
    x0 = x.tolist()
    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    xlabel = -22.5
    ylabel = 9.
    ystep = .6
    dy=.2
    dxl=1
    dxr=.2
    ppgplot.pgslw(deflw)  #line width
    #ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgline(xlin,ylin)

    ppgplot.pgslw(3)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    #ppgplot.pgsci(3)
    #x=N.compress((g1.final > 0),g1.sfr)
    x=N.compress((g1.final > 0) & (g1.finalsf > 0) & (g1.ewc > ewstarburst) & (g1.dkpc < (0.5*c.r200Mpc[1]*1000.)),g1.Mv)
    x1 = x.tolist()
    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)

    ylabel = ylabel - ystep

    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    #ppgplot.pgline(xlin,ylin)
    ppgplot.pgsls(1)
    ppgplot.pgslw(deflw)
    #ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    
    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    x=N.compress((g2.final > 0) & (g2.finalsf > 0) & (g2.ewc > ewstarburst) & (g2.dkpc < (0.5*c.r200Mpc[2]*1000.)),g2.Mv)
    x2 = x.tolist()
    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    xall = x0+x1+x2
    xall = N.array(xall,'f')
    x = xall
    ppgplot.pgslw(6)
    #ppgplot.pgsci(2)
    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    #ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgslw(2)  #line width
    #ppgplot.pgline(xlin,ylin)
    fr=.5

    x=N.compress((g0.final > 0) & (g0.finalsf > 0) & (g0.ewc > ewstarburst) & (g0.Mv < mvcut) & (g0.dkpc < (fr*c.r200Mpc[0]*1000)),g0.Mv)
    index = N.arange(0,len(g0.ewc),1,'i')
    xname=(N.compress((g0.specmemb > 0) & (g0.final > 0) & (g0.finalsf > 0) & (g0.ewc > ewstarburst) & (g0.Mv < mvcut) & (g0.dkpc < (fr*c.r200Mpc[0]*1000)),index))
    print "CL 1040 EW>40 spec members = ",len(xname)
    x0 = float(len(x))
    x=N.compress((g1.final > 0) & (g1.finalsf > 0) & (g1.ewc > ewstarburst) & (g1.Mv < mvcut) & (g1.dkpc < (fr*c.r200Mpc[1]*1000)),g1.Mv)
    index = N.arange(0,len(g1.ewc),1,'i')
    x1name=(N.compress((g1.specmemb > 0) & (g1.final > 0) & (g1.finalsf > 0) & (g1.ewc > ewstarburst) & (g1.Mv < mvcut) & (g1.dkpc < (fr*c.r200Mpc[1]*1000)),index))
    print "CL 1054 EW>40 spec members = ",len(x1name)
    x1 = float(len(x))
    x=N.compress((g2.final > 0) & (g2.finalsf > 0) & (g2.ewc > ewstarburst) & (g2.Mv < mvcut) & (g2.dkpc < (fr*c.r200Mpc[2]*1000)),g2.Mv)
    index = N.arange(0,len(g2.ewc),1,'f')
    x2name=(N.compress((g2.specmemb > 0) & (g2.final > 0) & (g2.finalsf > 0) & (g2.ewc > ewstarburst) & (g2.Mv < mvcut) & (g2.dkpc < (fr*c.r200Mpc[2]*1000)),index))
    print "CL 1216 EW>40 spec members = ",len(x2name)
    for i in x2name:
        i=int(i)
        print g2.name[i],g2.ewc[i],g2.ewo2[i],g2.sfrc[i]
    x2 = float(len(x))
    nstarburst = x0+x1+x2

    n0=float(len(N.compress((g0.final > 0) & (g0.Mv < mvcut) & (g0.dkpc < (fr*c.r200Mpc[0]*1000)),g0.Mv)))
    n1=float(len(N.compress((g1.final > 0) & (g1.Mv < mvcut) & (g1.dkpc < (fr*c.r200Mpc[1]*1000)),g1.Mv)))
    n2=float(len(N.compress((g2.final > 0) & (g2.Mv < mvcut) & (g2.dkpc < (fr*c.r200Mpc[2]*1000)),g2.Mv)))
    nall = float(n0+n1+n2)
    (a,b)=ratioerror(x0,n0)
    print "CL1040: Fraction of EW>40, Mv>",mvcut," = ",a,"+/-",b,"(",x0,"/",n0,")"

    (a,b)=ratioerror(x1,n1)
    print "CL1054: Fraction of EW>40, Mv>",mvcut," = ",a,"+/-",b,"(",x1,"/",n1,")"
    (a,b)=ratioerror(x2,n2)
    print "CL1216: Fraction of EW>40, Mv>",mvcut," = ",a,"+/-",b,"(",x2,"/",n2,")"
    (a,b)=ratioerror(nstarburst,nall)
    print "All: Fraction of EW>40, Mv>",mvcut," = ",a,"+/-",b,"(",nstarburst,"/",nall,")"

    #ewstarburst = 100.
    x = xall
    ppgplot.pgslw(2)
    ppgplot.pgsls(2)
    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    #print "Number in g2.ratios = ",len(g2.ratio)
    #for ratio in g2.ratio:
    #    print ratio
    #drawbinned(x,y,5)
    ppgplot.pgsci(1)
    xlin = N.array([-20.27,-20.27],'f')
    ylin = N.array([-10.,100],'f')
    ppgplot.pgsls(4)  #line width
    ppgplot.pgline(xlin,ylin)
    ppgplot.pgsls(1)  #line width
    #ppgplot.pgend()

def plotsfrkall():
    DATAMIN = 15.
    DATAMAX = 24.
    #print "ngal = ",len(g0.sfr)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-.5,15,0)
    ppgplot.pglab("K\ds\u (r < 2'')","SFR (h\d100\u\u-2\d M\d\(2281) \u yr\u-1\d)","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)
    x=N.compress((g0.finalsf > 0),g0.vltk2)
    #x=N.compress((g0.finalsf > 0),g0.mag6j)
    y=N.compress((g0.finalsf > 0),g0.sfrc)
    erry = N.compress((g0.finalsf > 0),g0.errsfrc)
    ppgplot.pgpt(x,y,3)
    errory(x,y,erry)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)
    xlabel = 21.2
    ylabel = 13.
    ystep = 1.3
    dy=.4
    dxl=1
    dxr=.4
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ppgplot.pgslw(3)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    x=N.compress((g1.finalsf > 0),g1.vltk2)
    #x=N.compress((g1.finalsf > 0),g1.mag6j)
    y=N.compress((g1.finalsf > 0),g1.sfrc)
    erry = N.compress((g1.finalsf > 0),g1.errsfrc)
    errory(x,y,erry)
    ppgplot.pgpt(x,y,17)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)

    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)

    ylabel = ylabel - ystep

    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgpt(xlin,ylin,17)
    ppgplot.pgsls(1)
    ppgplot.pgslw(deflw)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    
    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    x=N.compress((g2.finalsf > 0),g2.vltk2)
    #x=N.compress((g2.finalsf > 0),g2.mag6j)
    y=N.compress((g2.finalsf > 0),g2.sfrc)
    erry = N.compress((g2.finalsf > 0),g2.errsfrc)
    ppgplot.pgpt(x,y,7)
    errory(x,y,erry)


    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(2)  #line width
    ppgplot.pgpt(xlin,ylin,7)


    #print "Number in g2.ratios = ",len(g2.ratio)
    #for ratio in g2.ratio:
    #    print ratio
    #drawbinned(x,y,5)

    ZP = 24.94
    sfrconv = 2.
    mag = N.arange(14.,24,.05)
    #set JK=(23-mag)*1.2/5
    JK = 1.2*mag/mag
    #JK = 0.
    sfrlimit = 10**((ZP-mag-JK)/2.5)*ewmin/dlambdaJ*sfrconv
    ppgplot.pgsls(2)
    ppgplot.pgline(mag,sfrlimit)
    ppgplot.pgsci(1)
    #ppgplot.pgend()

def plotratioj(magj,fluxj,errfluxj,fluxn,errfluxn):
    DATAMIN = 17.
    DATAMAX = 25.
    ratio=fluxn/fluxj
    errratio=errfluxn/fluxj
    #print "ngal = ",len(g0.sfr)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    xmin=min(magj)
    xmax=max(magj)
    ymin=min(ratio)
    ymax=max(ratio)
    #ppgplot.pgenv(DATAMIN,DATAMAX,.0,.15,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0)
    ppgplot.pglab("J (mag-iso)","F\dNB\u/F\dJ\u","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    x=magj
    y=ratio
    erry = errratio
    ppgplot.pgpt(x,y,3)
    errory(x,y,erry)

def plotsfrj(magj,sfr,errsfr):
    DATAMIN = 17.
    DATAMAX = 25.
    #print "ngal = ",len(g0.sfr)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-.5,15.,0)
    ppgplot.pglab("J (mag-iso)","SFR (h\d100\u\u-2\d M\d\(2281) \u yr\u-1\d)","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    x=magj
    y=sfr
    erry = errsfr
    ppgplot.pgpt(x,y,3)
    errory(x,y,erry)
def plotsnsfr(sn,sfr,errsfr):
    DATAMIN = -1.
    DATAMAX = 6.
    sn=abs(sn)
    #print "ngal = ",len(g0.sfr)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(-15.,30.,DATAMIN,DATAMAX,0)
    if ncl == 5:
        ppgplot.pglab("SFR (h\d100\u\u-2\d M\d\(2281) \u yr\u-1\d)","SNR","MS1054-0321 (z=0.832)")
    if ncl == 6:
        ppgplot.pglab("SFR (h\d100\u\u-2\d M\d\(2281) \u yr\u-1\d)","SNR","HDFN1 (z=0.804)")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    x=sfr
    y=sn
    errx = errsfr
    ppgplot.pgpt(x,y,3)
    errorx(x,y,errx)
    y=N.arange((DATAMIN-2.),(DATAMAX+2.),1.,'f')
    x=2.*N.ones(len(y),'f')
    ppgplot.pgsls(4)
    ppgplot.pgslw(5)
    ppgplot.pgline(x,y)

    x=5.*N.ones(len(y),'f')
    ppgplot.pgsls(4)
    ppgplot.pgslw(5)
    #ppgplot.pgline(x,y)

    x=N.zeros(len(y),'f')
    ppgplot.pgsls(1)
    ppgplot.pgslw(5)
    ppgplot.pgline(x,y)

    x=N.arange(-50.,100.,1.,'f')
    y=N.zeros(len(x),'f')
    ppgplot.pgsls(1)
    ppgplot.pgslw(5)
    ppgplot.pgline(x,y)
    
    y=3*N.ones(len(x),'f')
    ppgplot.pgsls(4)
    ppgplot.pgslw(5)
    ppgplot.pgline(x,y)
    y=y/1.8
    ppgplot.pgsls(3)
    ppgplot.pgslw(5)
    #ppgplot.pgline(x,y)

def plotsfrjall():
    DATAMIN = 17.
    DATAMAX = 25.
    #print "ngal = ",len(g0.sfr)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-.5,15.,0)
    ppgplot.pglab("J (mag-iso)","SFR (h\d100\u\u-2\d M\d\(2281) \u yr\u-1\d)","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)
    #x=N.compress((g0.finalsf > 0),g0.vltk2)
    x=N.compress((g0.finalsf > 0),g0.magj)
    y=N.compress((g0.finalsf > 0),g0.sfrc)
    erry = N.compress((g0.finalsf > 0),g0.errsfrc)
    ppgplot.pgpt(x,y,3)
    errory(x,y,erry)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)
    xlabel = 21.8
    ylabel = 13.
    ystep = 1.3
    dy=.4
    dxl=1
    dxr=.2
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ppgplot.pgslw(3)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    #x=N.compress((g1.finalsf > 0),g1.vltk2)
    x=N.compress((g1.finalsf > 0),g1.magj)
    y=N.compress((g1.finalsf > 0),g1.sfrc)
    erry = N.compress((g1.finalsf > 0),g1.errsfrc)
    errory(x,y,erry)
    ppgplot.pgpt(x,y,17)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)

    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)

    ylabel = ylabel - ystep

    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgpt(xlin,ylin,17)
    ppgplot.pgsls(1)
    ppgplot.pgslw(deflw)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    
    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    #x=N.compress((g2.finalsf > 0),g2.vltk2)
    x=N.compress((g2.finalsf > 0),g2.magj)
    y=N.compress((g2.finalsf > 0),g2.sfrc)
    erry = N.compress((g2.finalsf > 0),g2.errsfrc)
    ppgplot.pgpt(x,y,7)
    errory(x,y,erry)


    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(2)  #line width
    ppgplot.pgpt(xlin,ylin,7)


    #print "Number in g2.ratios = ",len(g2.ratio)
    #for ratio in g2.ratio:
    #    print ratio
    #drawbinned(x,y,5)

    ZP = 24.92
    sfrconv = 2.
    mag = N.arange(14.,28,.05)
    #set JK=(23-mag)*1.2/5
    #JK = 1.2*mag/mag
    JK = 0.
    sfrlimit = 10**((ZP-mag-JK)/2.5)*ewmin/dlambdaJ*sfrconv
    ppgplot.pgsls(4)
    ppgplot.pgline(mag,sfrlimit)


    x0=N.compress((g0.finalsf > 0),g0.magj)
    a0=N.compress((g0.finalsf > 0),g0.isoarea)
    x1=N.compress((g1.finalsf > 0),g1.magj)
    a1=N.compress((g1.finalsf > 0),g1.isoarea)
    x2=N.compress((g2.finalsf > 0),g2.magj)
    a2=N.compress((g2.finalsf > 0),g2.isoarea)

    mag=x0.tolist() + x1.tolist() + x2.tolist()
    area=a0.tolist() + a1.tolist() + a2.tolist()
    bins=N.arange(16,25,.1)
    (magbin,areabin)=binitmin(mag,area,15)
    errfj = N.sqrt(areabin)*c.jnoisea[0]*(1.+c.jnoiseb[0]*N.sqrt(areabin))#use CL1216 noise
    errfj = N.sqrt(10**((24.92-bins)/2.5))


    sfrlimit=.24*errfj

    ppgplot.pgsls(2)
    #ppgplot.pgline(mag,sfrlimit)
    ppgplot.pgline(bins,sfrlimit)
    ppgplot.pgsci(1)
    errfj = 12.*N.sqrt(10**((24.92-bins)/2.5))#fit min iso-area versus mag
    errfj = N.sqrt(errfj)*c.jnoisea[2]*(1.+c.jnoiseb[2]*N.sqrt(errfj))#use CL1216 noise
    sfrlimit=1.*errfj
    ppgplot.pgsci(2)
    #ppgplot.pgline(bins,sfrlimit)
    ppgplot.pgsci(1)
    #ppgplot.pgend()
def plotmagisoarea():
    DATAMIN = 17.
    DATAMAX = 25.
    #print "ngal = ",len(g0.sfr)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-.5,1000.,0)
    ppgplot.pglab("J (mag-iso)","Area","")
    ppgplot.pgsls(1)#dotted

    x0=N.compress((g0.finalsf > 0),g0.magj)
    a0=N.compress((g0.finalsf > 0),g0.isoarea)
    x1=N.compress((g1.finalsf > 0),g1.magj)
    a1=N.compress((g1.finalsf > 0),g1.isoarea)
    x2=N.compress((g2.finalsf > 0),g2.magj)
    a2=N.compress((g2.finalsf > 0),g2.isoarea)

    mag=x0.tolist() + x1.tolist() + x2.tolist()
    area=a0.tolist() + a1.tolist() + a2.tolist()

    y=N.array(area,'f')
    x=N.array(mag,'f')
    ppgplot.pgpt(x,y,3)
    bins=N.arange(16,25,.1)
    (magbin,areabin)=binitmin(mag,area,15)
    errfj = N.sqrt(areabin)*c.jnoisea[0]*(1.+c.jnoiseb[0]*N.sqrt(areabin))#use CL1216 noise
    errfj = 12.*N.sqrt(10**((24.92-bins)/2.5))
    ppgplot.pgsls(2)
    #ppgplot.pgline(mag,sfrlimit)
    ppgplot.pgline(bins,errfj)


def plotsfrjallnsf():
    DATAMIN = 17.
    DATAMAX = 25.
    #print "ngal = ",len(g0.sfr)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-.5,20.,0)
    ppgplot.pglab("J (mag-iso)","SFR (h\d100\u\u-2\d M\d\(2281) \u yr\u-1\d)","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)
    #x=N.compress((g0.finalsf > 0),g0.vltk2)
    x=N.compress((g0.finalsf > 0),g0.magj)
    y=N.compress((g0.finalsf > 0),g0.sfrc)
    erry = N.compress((g0.finalsf > 0),g0.errsfrc)
    ppgplot.pgpt(x,y,3)
    errory(x,y,erry)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)
    xlabel = 22.5
    ylabel = 13.
    ystep = 1.3
    dy=.4
    dxl=1
    dxr=.4
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ppgplot.pgslw(3)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    #x=N.compress((g1.finalsf > 0),g1.vltk2)
    x=N.compress((g1.finalsf > 0),g1.magj)
    y=N.compress((g1.finalsf > 0),g1.sfrc)
    erry = N.compress((g1.finalsf > 0),g1.errsfrc)
    errory(x,y,erry)
    ppgplot.pgpt(x,y,17)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)

    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)

    ylabel = ylabel - ystep

    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgpt(xlin,ylin,17)
    ppgplot.pgsls(1)
    ppgplot.pgslw(deflw)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    
    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    #x=N.compress((g2.finalsf > 0),g2.vltk2)
    x=N.compress((g2.finalsf > 0),g2.magj)
    y=N.compress((g2.finalsf > 0),g2.sfrc)
    erry = N.compress((g2.finalsf > 0),g2.errsfrc)
    ppgplot.pgpt(x,y,7)
    errory(x,y,erry)


    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(2)  #line width
    ppgplot.pgpt(xlin,ylin,7)


    #print "Number in g2.ratios = ",len(g2.ratio)
    #for ratio in g2.ratio:
    #    print ratio
    #drawbinned(x,y,5)

    ZP = 25.65
    sfrconv = 2.
    mag = N.arange(14.,28,.05)
    #set JK=(23-mag)*1.2/5
    #JK = 1.2*mag/mag
    JK = 0.
    sfrlimit = 10**((ZP-mag-JK)/2.5)*ewmin/dlambdaJ*sfrconv
    ppgplot.pgsls(2)
    ppgplot.pgline(mag,sfrlimit)
    ppgplot.pgsci(1)
    #ppgplot.pgend()
    mag = N.arange(14.,25,.05)
    y = N.ones(len(mag),'f')
    ppgplot.pgsls(4)
    ppgplot.pgline(mag,y)

def plotsfrewall():
    DATAMIN = -20.
    DATAMAX = 200.
    #print "ngal = ",len(g0.sfr)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-.5,20.,0)
    ppgplot.pglab(ewlabel,"SFR (h\d100\u\u-2\d M\d\(2281) \u yr\u-1\d)","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)
    #x=N.compress((g0.finalsf > 0),g0.vltk2)
    x=N.compress((g0.finalsf > 0),g0.ewc)
    y=N.compress((g0.finalsf > 0),g0.sfrc)
    erry = N.compress((g0.finalsf > 0),g0.errsfrc)
    ppgplot.pgpt(x,y,3)
    errory(x,y,erry)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)
    xlabel = 120
    ylabel = 18.
    ystep = 1.3
    dy=.4
    dxl=5
    dxr=5
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ppgplot.pgslw(3)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    x=N.compress((g1.finalsf > 0),g1.ewc)
    y=N.compress((g1.finalsf > 0),g1.sfrc)
    erry = N.compress((g1.finalsf > 0),g1.errsfrc)
    errory(x,y,erry)
    ppgplot.pgpt(x,y,17)

    ylabel = ylabel - ystep

    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgpt(xlin,ylin,17)
    ppgplot.pgsls(1)
    ppgplot.pgslw(deflw)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    
    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    x=N.compress((g2.finalsf > 0),g2.ewc)
    y=N.compress((g2.finalsf > 0),g2.sfrc)
    erry = N.compress((g2.finalsf > 0),g2.errsfrc)
    ppgplot.pgpt(x,y,7)
    errory(x,y,erry)


    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(2)  #line width
    ppgplot.pgpt(xlin,ylin,7)



def plotewkall():
    DATAMIN = 15.
    DATAMAX = 24.
    #print "ngal = ",len(g0.sfr)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-10,280,0)
    ppgplot.pglab("K\ds\u (r < 2'')","EW H\ga+[NII] (\(2078))","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)
    x=N.compress((g0.finalsf > 0),g0.vltk2)
    y=N.compress((g0.finalsf > 0),g0.ewc)
    erry = N.compress((g0.finalsf > 0),g0.errewc)
    errory(x,y,erry)
    ppgplot.pgpt(x,y,3)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)
    xlabel = 16.2
    ylabel = 250.
    ystep = 22.
    dy=5
    dxl=1
    dxr=.4
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ppgplot.pgslw(3)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    x=N.compress((g1.finalsf > 0),g1.vltk2)
    y=N.compress((g1.finalsf > 0),g1.ewc)
    erry = N.compress((g1.finalsf > 0),g1.errewc)
    ppgplot.pgpt(x,y,17)
    ylabel = ylabel - ystep
    errory(x,y,erry)

    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgpt(xlin,ylin,17)
    ppgplot.pgsls(1)
    ppgplot.pgslw(deflw)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    
    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    x=N.compress((g2.finalsf > 0),g2.vltk2)
    y=N.compress((g2.finalsf > 0),g2.ewc)
    erry = N.compress((g2.finalsf > 0),g2.errewc)
    ppgplot.pgpt(x,y,7)
    errory(x,y,erry)


    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(2)  #line width
    ppgplot.pgpt(xlin,ylin,7)

    ZP = 24.94
    sfrconv = 2.
    sfrmin = .2
    mag = N.arange(14.,24,.05)
    y = N.ones(len(mag),'f')
    y = ewmin*y
    ppgplot.pgsls(4)
    ppgplot.pgline(mag,y)
    y = y + 30.
    ppgplot.pgline(mag,y)
    y = y + 60.
    ppgplot.pgline(mag,y)
    JK = 1.2*mag/mag
    ewlimit = sfrmin/sfrconv/10**((ZP-mag-JK)/2.5)*dlambdaJ
    ppgplot.pgsls(2)
    ppgplot.pgline(mag,ewlimit)
    ppgplot.pgsci(1)
    #ppgplot.pgend()

def plotewj(magj,ew,errew):
    DATAMIN = 17.
    DATAMAX = 25.
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-100,400,0)
    ppgplot.pglab("J (mag-iso)","EW H\ga+[NII] (\(2078))","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    x=magj
    y=ew
    erry = errew
    errory(x,y,erry)
    ppgplot.pgpt(x,y,3)

    x=N.arange(0.,30.,3.)
    y=N.zeros(len(x),'f')
    ppgplot.pgsls(4)
    ppgplot.pgline(x,y)
def plotewjall():
    DATAMIN = 17.
    DATAMAX = 25.
    #print "ngal = ",len(g0.sfr)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-10,280,0)
    ppgplot.pglab("J (mag-iso)","EW H\ga+[NII] (\(2078))","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)
    x=N.compress((g0.finalsf > 0),g0.magj)
    y=N.compress((g0.finalsf > 0),g0.ewc)
    erry = N.compress((g0.finalsf > 0),g0.errewc)
    errory(x,y,erry)
    ppgplot.pgpt(x,y,3)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)
    xlabel = 18.
    ylabel = 250.
    ystep = 22.
    dy=5
    dxl=1
    dxr=.4
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ppgplot.pgslw(3)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    x=N.compress((g1.finalsf > 0),g1.magj)
    y=N.compress((g1.finalsf > 0),g1.ewc)
    erry = N.compress((g1.finalsf > 0),g1.errewc)
    ppgplot.pgpt(x,y,17)
    ylabel = ylabel - ystep
    errory(x,y,erry)

    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgpt(xlin,ylin,17)
    ppgplot.pgsls(1)
    ppgplot.pgslw(deflw)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    
    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    x=N.compress((g2.finalsf > 0),g2.magj)
    y=N.compress((g2.finalsf > 0),g2.ewc)
    erry = N.compress((g2.finalsf > 0),g2.errewc)
    ppgplot.pgpt(x,y,7)
    errory(x,y,erry)


    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(2)  #line width
    ppgplot.pgpt(xlin,ylin,7)

    ZP = 24.9
    sfrconv = 2.
    sfrmin = .2
    mag = N.arange(14.,27,.05)
    y = N.ones(len(mag),'f')
    y = ewmin*y
    ppgplot.pgsls(4)
    ppgplot.pgline(mag,y)
    y = 40*N.ones(len(mag),'f')
    ppgplot.pgline(mag,y)
    JK = 0
    #ewlimit = ewmin+sfrmin/sfrconv/10**((ZP-mag-JK)/2.5)*dlambdaJ
    ewlimit = sfrmin/sfrconv/10**((ZP-mag-JK)/2.5)*dlambdaJ
    ppgplot.pgsls(2)
    #ppgplot.pgline(mag,ewlimit)
    ppgplot.pgsci(1)
    #ppgplot.pgend()

    bins=N.arange(16,25,.1)
    errfj = 12.*N.sqrt(10**((24.92-bins)/2.5))#fit min iso-area versus mag
    errfj = N.sqrt(errfj)*c.jnoisea[2]*(1.+c.jnoiseb[2]*N.sqrt(errfj))#use CL1216 noise
    ewlimit=errfj/sfrconv/10**((ZP-bins-JK)/2.5)*dlambdaJ/2.
    #sfrlimit=.24*errfj
    ppgplot.pgsci(1)
    ppgplot.pgsls(2)
    #ppgplot.pgline(mag,sfrlimit)
    ppgplot.pgline(bins,ewlimit)
    ppgplot.pgsci(1)


def plotcompmmtedij():
    DATAMIN = 17.
    DATAMAX = 25.
    psplotinit("compmmtedij.ps")
    #print "ngal = ",len(g0.sfr)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-3,3,0)
    ppgplot.pglab("J\dEDisCS\u (r < 2'')","J\dMMT - J\dEDisCS\u (r < 2'')","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    compsub(g0.mag6j,g0.errmag6j,g0.vltj2,g0.errvltj2,g0.memb)
    compsub(g1.mag6j,g1.errmag6j,g1.vltj2,g1.errvltj2,g1.memb)
    compsub(g2.mag6j,g2.errmag6j,g2.vltj2,g2.errvltj2,g2.memb)

    ppgplot.pgslw(3)
    x=N.arange((DATAMIN-2),(DATAMAX+2))
    y=N.zeros(len(x),'f')
    ppgplot.pgline(x,y)
    ppgplot.pgsls(2)

    y=y+.15
    ppgplot.pgline(x,y)
    y=-1.*y
    ppgplot.pgline(x,y)
    ppgplot.pgsls(1)
    ppgplot.pgend()
    #os.system("cp compmmtedij.ps /Users/rfinn/clusters/papers/paper2/.")
    #os.system("cp compmmtedij.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")
def compsub(mag6j,errmag6j,vltj2,errvltj2,final):
    mag6j=N.compress(final > 0,mag6j)
    errmag6j=N.compress(final > 0,errmag6j)
    vltj2=N.compress(final > 0,vltj2)
    errvltj2=N.compress(final > 0,errvltj2)
    x=mag6j    
    x2=vltj2
    y=x-x2
    print "Average diff in J<23 mags = ",N.average(N.compress(vltj2<23,y))
    print "std in J mags = ",pylab.std(N.compress(vltj2<23,y))
    print "Average diff in J<21.5 mags = ",N.average(N.compress(vltj2<20.5,y))
    print "std in J mags = ",pylab.std(N.compress(vltj2<20.5,y))
    print "Median diff in J mags = ",pylab.median(N.compress(vltj2<23,y))
    ppgplot.pgpt(x2,y,17)
    err=errmag6j
    err2=errvltj2
    erry = N.sqrt(err**2 + err2**2)
    errory(x2,y,erry)


def plothistj():
    minmemb = 0.
    DATAMIN = 17.
    DATAMAX = 26.
    NBIN = int((DATAMAX-DATAMIN)*2.)    
    #print "ngal = ",len(g0.sfr)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,0,35,0)
    ppgplot.pglab("J (r<2'')","Number of Galaxies","")
    #ppgplot.pgsls(2)#dotted
    ppgplot.pgslw(4)  #line width
    ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)
    x = g0.mag6j
    memb=g0.memb
    x=N.compress(memb > minmemb, x)
    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ppgplot.pgslw(1)#dotted
    x = ntt0.magj
    memb=ntt0.memb
    x=N.compress(memb > minmemb, x)

    #print "number of NTT galaxies = ",len(x)
    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ppgplot.pgsls(1)#solid
    ppgplot.pgslw(4)
    #ppgplot.pgsls(4)#dot-dash-dot-dash
    ppgplot.pgsci(3)
    #x=N.compress((abs(g1.ew) > ewmin),g1.sfr)
    x = g1.mag6j
    memb=g1.memb
    x=N.compress(memb > minmemb, x)

    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ppgplot.pgslw(1)#dotted
    x = ntt1.magj
    memb=ntt1.memb
    x=N.compress(memb > minmemb, x)

    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ppgplot.pgsls(1)#dot-dash-dot-dash
    ppgplot.pgsci(2)
    ppgplot.pgslw(4)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    x = g2.mag6j
    memb=g2.memb
    x=N.compress(memb > minmemb, x)
    
    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ppgplot.pgslw(1)#dotted
    x = ntt2.magj
    memb=ntt2.memb
    x=N.compress(memb > minmemb, x)

    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ppgplot.pgsls(1)#solid

    #print "Number in g2.ratios = ",len(g2.ratio)
    #for ratio in g2.ratio:
    #    print ratio
    #drawbinned(x,y,5)
    ppgplot.pgsci(1)
    #ppgplot.pgend()

def plotnbjvj(x,y,ratio):#NB/J vs J
    ppgplot.pgslw(4)  #line width
    #y=g0.fn/g0.fj
    #x=g0.magj
    xmin=15
    xmax=26
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(xmin,xmax,0.,.22,0)
    ppgplot.pglab("J (magiso)","F\dn\u/F\dJ\u","")
    ppgplot.pgsci(4)
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsci(1)
    ppgplot.pgsls(4)
    xl = N.arange(xmin-2.,xmax+2.,1,'f')
    yl = N.ones(len(xl),'f')
    yl = ratio*yl
    #print xl
    #print yl

    ppgplot.pgline(xl,yl)
    ppgplot.pgsls(1)
def plotsumsfrmag(magmmtin,sfrmmtin,membin,magnttin,membinntt,signifin):
    #print signifin
    ppgplot.pgslw(6)
    magmmt = N.compress((membin > 0),magmmtin)
    sfrmmt = N.compress((membin > 0),sfrmmtin)
    signif = N.compress((membin > 0),signifin)
    magntt = N.compress((membinntt > 0),magnttin)

    ppgplot.pgslw(6)
    minmag = 17.
    maxmag = 26.
    magstep = 1.
    mlim = 23. #max mag for calculating completeness
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(minmag,maxmag,-0.1,38.,0)
    ppgplot.pglab("J (2'')","\gSSFR/\gDmag","")
    #sum sfr for all significant detections
    #for i in range(len(magmmt)):
    #    print magmmt[i],sfrmmt[i],signif[i]

    x = N.compress(signif > 0,magmmt)
    y = N.compress(signif > 0,sfrmmt)
    (xbin,ybin)=binitsum(x,y,minmag,maxmag,magstep)#get total sfr in each mag bin
    drawhist(xbin,ybin)
    #completeness
    #x = g0.mag6j
    x = magmmt
    (xbin,y) = binithist(x,minmag,maxmag,magstep)#get total # of gal in each mag bin
    #x = ntt0.magj
    x = magntt
    (xbin,yntt) = binithist(x,minmag,maxmag,magstep)#get total # of gal in each mag bin for ediscs sample
    ybincorr = N.zeros(len(ybin),'f')
    for i in range(len(ybin)):
        if (y[i] > 0):
            ybincorr[i] = ybin[i]*yntt[i]/y[i]
    ppgplot.pgsls(4)
    #ppgplot.pgsci(2)
    drawhist(xbin,ybincorr)
    #for i in range(len(xbin)):
    #    print i, xbin[i], ybin[i],y[i],yntt[i],ybincorr[i]
    a = float(sum(N.compress(xbin < mlim,y)))
    b = float(sum(N.compress(xbin < mlim,yntt)))
    print "Total Ngal detected = ",a
    print "Total Ngal detected by ediscs= ",b
    #print "Total Ngal w/completeness correction = ",(ybincorr)
    #print "Completeness within m < 23 = ",a/b
    print "Total SFR detected w/m<",mlim," = ",N.sum(N.compress(xbin < mlim,ybin))
    print "Total w/completeness correction = ",N.sum(N.compress(xbin < mlim,ybincorr))
    print "Completeness of SFR within m < 23 = ",N.sum(N.compress(xbin < mlim,ybin))/N.sum(N.compress(xbin < mlim,ybincorr))
    ppgplot.pgsls(1)
    ppgplot.pgsci(1)
    compl = N.sum(N.compress(xbin < mlim,ybin))/N.sum(N.compress(xbin < mlim,ybincorr))
    return compl

def plotgal(xg,yg,final,finalsf,ncl):
    ppgplot.pgslw(6)
    ppgplot.pgsls(1)

    xmin = (-1.*c.pscale[ncl]*(c.xc[ncl]-1))
    xmax = (c.pscale[ncl]*(c.xmax[ncl]-c.xc[ncl]))    
    ymin = (-1.*c.pscale[ncl]*(c.yc[ncl]-1))
    ymax = (c.pscale[ncl]*(c.ymax[ncl]-c.yc[ncl]))
    ppgplot.pgbox("",0.0,0,"L",0.0,0)

    dx=5.
    ppgplot.pgenv(xmin-dx,xmax+dx,ymin-dx,ymax+dx,0)
    ppgplot.pglab("\gD DEC (\")","\gD RA (\")","")
    ppgplot.pgtext(-4,-4,"X")
    r = (0.5*c.r200pix[ncl]*c.pscale[ncl])
    ppgplot.pgslw(1)
    ppgplot.pgsls(2)
    ppgplot.pgsfs(2)
    ppgplot.pgcirc(0,0,r)
    #print "cluster ",ncl," r200: ",c.r200pix[ncl],c.r200Mpc[ncl], " Mpc"

    ppgplot.pgslw(3)
    ppgplot.pgsls(1)

    x = (xg - c.xc[ncl])*c.pscale[ncl]
    y = (yg - c.yc[ncl])*c.pscale[ncl]
    #ppgplot.pgpt(x,y,1)
    x = (N.compress((final > 0) & (finalsf < 1), xg) - c.xc[ncl])*c.pscale[ncl]
    y = (N.compress((final > 0) & (finalsf < 1), yg) - c.yc[ncl])*c.pscale[ncl]
    #ppgplot.pgpt(x,y,4)
    ppgplot.pgpt(x,y,22)
    x = (N.compress((final > 0) & (finalsf > 0), xg) - c.xc[ncl])*c.pscale[ncl]
    y = (N.compress((final > 0) & (finalsf > 0), yg) - c.yc[ncl])*c.pscale[ncl]
    ppgplot.pgpt(x,y,18)


def plotsfrdall():

    DATAMIN = -15.
    DATAMAX = 600.

    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-.5,18,0)
    ppgplot.pglab("d (kpc)","SFR (h\d100\u\u-2\d M\d\(2281) \u yr\u-1\d)","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)
    x = N.compress((g0.finalsf > 0),g0.dkpc)
    y = N.compress((g0.finalsf > 0),g0.sfrc)
    erry = N.compress((g0.finalsf > 0),g0.errsfrc)
    ppgplot.pgpt(x,y,3)
    xs = []
    ys = []
    x1 = x.tolist()
    y1 = y.tolist()
    #for i in range(len(x)):
    #    print x[i],x1[i]
    #drawbinned(x,y,3)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)
    xlabel = 40.
    ylabel = 16.2
    ystep = 1.2
    dy=.4
    dxl=10
    dxr=10.
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ppgplot.pgslw(3)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    x=N.compress((g1.finalsf > 0),g1.dkpc)
    y=N.compress((g1.finalsf > 0),g1.sfrc)
    erry = N.compress((g1.finalsf > 0),g1.errsfrc)
    ppgplot.pgpt(x,y,17)
    x2 = x.tolist()
    y2 = y.tolist()
    #drawbinned(x,y,5)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)

    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)

    ylabel = ylabel - ystep

    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgpt(xlin,ylin,17)
    ppgplot.pgsls(1)
    ppgplot.pgslw(deflw)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    
    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    x=N.compress((g2.finalsf > 0),g2.dkpc)
    y=N.compress((g2.finalsf > 0),g2.sfrc)
    erry = N.compress((g2.finalsf > 0),g2.errsfrc)
    ppgplot.pgpt(x,y,7)
    x3 = x.tolist()
    y3 = y.tolist()
    xs = x1+x2+x3
    ys = y1+y2+y3
    xs = N.array(xs,'f')
    ys = N.array(ys,'f')

    #for i in range(len(xs)):
    #    print i, xs[i],ys[i]
    #print "total should be ",(N.sum(g0.finalsf)+N.sum(g1.finalsf) + N.sum(g2.finalsf)) 
    ppgplot.pgslw(6)
    drawbinned(xs,ys,5)
    print "SFR vs d"
    dospear(xs,ys)

    ppgplot.pgslw(deflw)
    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(2)  #line width
    ppgplot.pgpt(xlin,ylin,7)

    y = N.arange(-20.,50.,1.)
    x = N.ones(len(y),'f')
    x = 305.*x
    ppgplot.pgsls(2)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(1)



def plotewdall():

    DATAMIN = -15.
    DATAMAX = 600.

    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-10.,280.,0)
    #ppgplot.pgenv(DATAMIN,DATAMAX,-10.,28.,0)
    ppgplot.pglab("d (kpc)","EW H\ga+[NII] (\(2078))","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    x = N.compress((g0.finalsf > 0),g0.dkpc)
    y = N.compress((g0.finalsf > 0),g0.ewc)
    erry = N.compress((g0.finalsf > 0),g0.errewc)
    ppgplot.pgpt(x,y,3)
    x1 = x.tolist()
    y1 = y.tolist()

    #drawbinned(x,y,3)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)
    xlabel = 40.
    ylabel = 250.
    ystep = 20.
    dy=4
    dxl=10
    dxr=10.
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ppgplot.pgslw(3)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    x=N.compress((g1.finalsf > 0),g1.dkpc)
    y=N.compress((g1.finalsf > 0),g1.ewc)
    erry = N.compress((g1.finalsf > 0),g1.errewc)
    ppgplot.pgpt(x,y,17)
    x2 = x.tolist()
    y2 = y.tolist()

    #drawbinned(x,y,5)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)

    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)

    ylabel = ylabel - ystep

    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgpt(xlin,ylin,17)
    ppgplot.pgsls(1)
    ppgplot.pgslw(deflw)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    
    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    x=N.compress((g2.finalsf > 0),g2.dkpc)
    y=N.compress((g2.finalsf > 0),g2.ewc)
    erry = N.compress((g2.finalsf > 0),g2.errewc)
    ppgplot.pgpt(x,y,7)
    x3 = x.tolist()
    y3 = y.tolist()
    xs = x1 + x2+ x3
    ys = y1 + y2+ y3
    xs = N.array(xs,'f')
    ys = N.array(ys,'f')
    ppgplot.pgslw(6)
    drawbinned(xs,ys,5)
    print "EW vs d"
    dospear(xs,ys)

    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(2)  #line width
    ppgplot.pgpt(xlin,ylin,7)

    y = N.arange(-50.,500.,1.)
    x = N.ones(len(y),'f')
    x = 305.*x
    ppgplot.pgsls(2)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(1)


def plotsffracd():
    dmin = -15.
    dmax = 600.
    dstep = 50
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(dmin,dmax,0.,1.,0)
    ppgplot.pglab("d (kpc)","Star-Forming Fraction","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)

    x = N.compress((g0.final > 0),g0.dkpc)
    y1 = N.compress((g0.final > 0),g0.finalsf)
    y2 = N.compress((g0.final > 0),g0.final)
    (xbin,y1sum)=binitsumequal(x,y1,5)#get total sfr in each mag bin
    (xbin,y2sum)=binitsumequal(x,y2,5)#get total sfr in each mag bin
    sffrac = (y1sum)/(y2sum)
    sffrac = N.array(sffrac,'f')
    ppgplot.pgpt(xbin,sffrac,3)
    xs=[]
    y1s=[]
    y2s=[]
    xs=x.tolist()
    y1s=y1.tolist()
    y2s=y2.tolist()


    x = N.compress((g1.final > 0),g1.dkpc)
    y1 = N.compress((g1.final > 0),g1.finalsf)
    y2 = N.compress((g1.final > 0),g1.final)
    (xbin,y1sum)=binitsumequal(x,y1,5)#get total sfr in each mag bin
    (xbin,y2sum)=binitsumequal(x,y2,5)#get total sfr in each mag bin
    sffrac = (y1sum)/(y2sum)
    sffrac = N.array(sffrac,'f')
    ppgplot.pgpt(xbin,sffrac,17)
    xs=xs+x.tolist()
    y1s=y1s+y1.tolist()
    y2s=y2s+y2.tolist()


    x = N.compress((g2.final > 0),g2.dkpc)
    y1 = N.compress((g2.final > 0),g2.finalsf)
    y2 = N.compress((g2.final > 0),g2.final)
    (xbin,y1sum)=binitsumequal(x,y1,5)#get total sfr in each mag bin
    (xbin,y2sum)=binitsumequal(x,y2,5)#get total sfr in each mag bin
    sffrac = (y1sum)/(y2sum)
    sffrac = N.array(sffrac,'f')
    ppgplot.pgpt(xbin,sffrac,7)

    xs=xs+x.tolist()
    y1s=y1s+y1.tolist()
    y2s=y2s+y2.tolist()
    x=N.array(xs,'f')
    y1=N.array(y1s,'f')
    y2=N.array(y2s,'f')
    (xbin,y1sum)=binitsumequal(x,y1,5)#get total sfr in each mag bin
    (xbin,y2sum)=binitsumequal(x,y2,5)#get total sfr in each mag bin
    sffrac = (y1sum)/(y2sum)
    sffrac = N.array(sffrac,'f')
    ppgplot.pgsls(1)
    ppgplot.pgline(xbin,sffrac)
    #for i in range(len(xbin)):
    #    print "hey", i, xbin[i], sffrac[i]
    xlabel = 40.
    ylabel = .9
    ystep = 0.07
    dy=.015
    dxl=10
    dxr=15.
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,17)

    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,7)
def plotsfrdr200all():

    DATAMIN = -.05
    DATAMAX = 1.

    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-.5,19,0)
    ppgplot.pglab("d/R\d200\u","SFR (h\d100\u\u-2\d M\d\(2281) \u yr\u-1\d)","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)
    x = N.compress((g0.finalsf > 0),g0.dr200)
    y = N.compress((g0.finalsf > 0),g0.sfrc)
    erry = N.compress((g0.finalsf > 0),g0.errsfrc)
    ppgplot.pgpt(x,y,3)
    xs = []
    ys = []
    x1 = x.tolist()
    y1 = y.tolist()
    #for i in range(len(x)):
    #    print x[i],x1[i]
    #drawbinned(x,y,3)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)
    xlabel = 0.57
    ylabel = 16.7
    ystep = 1.2
    dy=.4
    dxl=10
    dxr=.05
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ppgplot.pgslw(3)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    x=N.compress((g1.finalsf > 0),g1.dr200)
    y=N.compress((g1.finalsf > 0),g1.sfrc)
    erry = N.compress((g1.finalsf > 0),g1.errsfrc)
    ppgplot.pgpt(x,y,17)
    x2 = x.tolist()
    y2 = y.tolist()
    #drawbinned(x,y,5)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)

    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)

    ylabel = ylabel - ystep

    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgpt(xlin,ylin,17)
    ppgplot.pgsls(1)
    ppgplot.pgslw(deflw)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    
    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    x=N.compress((g2.finalsf > 0),g2.dr200)
    y=N.compress((g2.finalsf > 0),g2.sfrc)
    erry = N.compress((g2.finalsf > 0),g2.errsfrc)
    ppgplot.pgpt(x,y,7)
    x3 = x.tolist()
    y3 = y.tolist()
    xs = x1+x2+x3
    ys = y1+y2+y3
    xs = N.array(xs,'f')
    ys = N.array(ys,'f')

    #for i in range(len(xs)):
    #    print i, xs[i],ys[i]
    #print "total should be ",(N.sum(g0.finalsf)+N.sum(g1.finalsf) + N.sum(g2.finalsf)) 
    ppgplot.pgslw(6)
    drawbinned(xs,ys,3)
    print "SFR vs d/R200"
    dospear(xs,ys)

    ppgplot.pgslw(deflw)
    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(2)  #line width
    ppgplot.pgpt(xlin,ylin,7)

    y = N.arange(-20.,50.,1.)
    x = N.ones(len(y),'f')
    x = 305.*x
    ppgplot.pgsls(2)
    #ppgplot.pgline(x,y)
    ppgplot.pgsci(1)

def plotewdr200all():

    DATAMIN = -.05
    DATAMAX = 1.

    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-10.,280.,0)
    #ppgplot.pgenv(DATAMIN,DATAMAX,-10.,28.,0)
    ppgplot.pglab("d/R\d200\u","EW H\ga+[NII] (\(2078))","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    x = N.compress((g0.finalsf > 0),g0.dr200)
    y = N.compress((g0.finalsf > 0),g0.ewc)
    erry = N.compress((g0.finalsf > 0),g0.errewc)
    ppgplot.pgpt(x,y,3)
    ppgplot.pgslw(3)
    ppgplot.pgsls(1)

    #drawbinned(x,y,3)
    x1 = x.tolist()
    y1 = y.tolist()

    #drawbinned(x,y,3)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)
    xlabel = 0.6
    ylabel = 250.
    ystep = 20.
    dy=4
    dxl=10
    dxr=.05
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)
    ppgplot.pgslw(3)
    ppgplot.pgsls(1)

    #drawbinned(x,y,3)

    ppgplot.pgslw(3)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    x=N.compress((g1.finalsf > 0),g1.dr200)
    y=N.compress((g1.finalsf > 0),g1.ewc)
    erry = N.compress((g1.finalsf > 0),g1.errewc)
    ppgplot.pgpt(x,y,17)
    ppgplot.pgslw(2)
    ppgplot.pgsls(2)

    #drawbinned(x,y,3)
    x2 = x.tolist()
    y2 = y.tolist()

    #drawbinned(x,y,5)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)

    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)

    ylabel = ylabel - ystep

    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgpt(xlin,ylin,17)
    ppgplot.pgsls(1)
    ppgplot.pgslw(deflw)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    
    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    x=N.compress((g2.finalsf > 0),g2.dr200)
    y=N.compress((g2.finalsf > 0),g2.ewc)
    erry = N.compress((g2.finalsf > 0),g2.errewc)
    ppgplot.pgpt(x,y,7)
    ppgplot.pgslw(1)
    ppgplot.pgsls(1)

    #drawbinned(x,y,3)

    x3 = x.tolist()
    y3 = y.tolist()
    xs = x1 + x2+ x3
    ys = y1 + y2+ y3
    xs = N.array(xs,'f')
    ys = N.array(ys,'f')
    ppgplot.pgslw(6)
    drawbinned(xs,ys,3)
    print "EW vs d/R200"
    dospear(xs,ys)

    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(2)  #line width
    ppgplot.pgpt(xlin,ylin,7)

    y = N.arange(-50.,500.,1.)
    x = N.ones(len(y),'f')
    x = 305.*x
    ppgplot.pgsls(2)
    #ppgplot.pgline(x,y)
    ppgplot.pgsci(1)

def plotsfrdNNall():

    DATAMIN = -.05
    DATAMAX = 50.

    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-.5,18,0)
    ppgplot.pglab("d\dNN\u (arcsec)","SFR (h\d100\u\u-2\d M\d\(2281) \u yr\u-1\d)","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)
    #x = N.compress((g0.finalsf > 0),g0.dNN)
    x = N.compress((g0.finalsf > 0),g0.dNN)
    y = N.compress((g0.finalsf > 0),g0.sfr)
    erry = N.compress((g0.finalsf > 0),g0.errsfrc)

    ppgplot.pgpt(x,y,3)
    xs = []
    ys = []
    x1 = x.tolist()
    y1 = y.tolist()
    #for i in range(len(x)):
    #    print x[i],x1[i]
    #drawbinned(x,y,3)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)
    xlabel = 0.
    ylabel = 16.2
    ystep = 1.2
    dy=.4
    dxl=10
    dxr=.10
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ppgplot.pgslw(3)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    x=N.compress((g1.finalsf > 0),g1.dNN)
    y=N.compress((g1.finalsf > 0),g1.sfrc)
    erry = N.compress((g1.finalsf > 0),g1.errsfrc)
    ppgplot.pgpt(x,y,17)
    x2 = x.tolist()
    y2 = y.tolist()
    #drawbinned(x,y,5)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)

    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)

    ylabel = ylabel - ystep

    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgpt(xlin,ylin,17)
    ppgplot.pgsls(1)
    ppgplot.pgslw(deflw)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    
    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    x=N.compress((g2.finalsf > 0),g2.dNN)
    y=N.compress((g2.finalsf > 0),g2.sfrc)
    erry = N.compress((g2.finalsf > 0),g2.errsfrc)
    ppgplot.pgpt(x,y,7)
    x3 = x.tolist()
    y3 = y.tolist()
    xs = x1+x2+x3
    ys = y1+y2+y3
    xs = N.array(xs,'f')
    ys = N.array(ys,'f')

    #for i in range(len(xs)):
    #    print i, xs[i],ys[i]
    #print "total should be ",(N.sum(g0.finalsf)+N.sum(g1.finalsf) + N.sum(g2.finalsf)) 
    ppgplot.pgslw(6)
    drawbinned(xs,ys,3)
    print "SFR vs dNN"
    dospear(xs,ys)

    ppgplot.pgslw(deflw)
    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(2)  #line width
    ppgplot.pgpt(xlin,ylin,7)

    y = N.arange(-20.,50.,1.)
    x = N.ones(len(y),'f')
    x = 305.*x
    ppgplot.pgsls(2)
    #ppgplot.pgline(x,y)
    ppgplot.pgsci(1)

def plotewdNNall():

    DATAMIN = -.05
    DATAMAX = 60.

    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-10.,280.,0)
    #ppgplot.pgenv(DATAMIN,DATAMAX,-10.,28.,0)
    ppgplot.pglab("d/R\d200\u","EW H\ga+[NII] (\(2078))","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    x = N.compress((g0.final > 0),g0.dNN)
    y = N.compress((g0.final > 0),g0.ewc)
    erry = N.compress((g0.final > 0),g0.errewc)
    ppgplot.pgpt(x,y,3)
    x1 = x.tolist()
    y1 = y.tolist()

    #drawbinned(x,y,3)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)
    xlabel = 0.
    ylabel = 250.
    ystep = 20.
    dy=4
    dxl=10
    dxr=.10
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ppgplot.pgslw(3)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    x=N.compress((g1.final > 0),g1.dNN)
    y=N.compress((g1.final > 0),g1.ewc)
    erry = N.compress((g1.final > 0),g1.errewc)
    ppgplot.pgpt(x,y,17)
    x2 = x.tolist()
    y2 = y.tolist()

    #drawbinned(x,y,5)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)

    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)

    ylabel = ylabel - ystep

    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgpt(xlin,ylin,17)
    ppgplot.pgsls(1)
    ppgplot.pgslw(deflw)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    
    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    x=N.compress((g2.final > 0),g2.dNN)
    y=N.compress((g2.final > 0),g2.ewc)
    erry = N.compress((g2.final > 0),g2.errewc)
    ppgplot.pgpt(x,y,7)
    x3 = x.tolist()
    y3 = y.tolist()
    xs = x1 + x2+ x3
    ys = y1 + y2+ y3
    xs = N.array(xs,'f')
    ys = N.array(ys,'f')
    ppgplot.pgslw(6)
    drawbinned(xs,ys,5)
    print "EW vs dNN"
    dospear(xs,ys)

    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(2)  #line width
    ppgplot.pgpt(xlin,ylin,7)

    y = N.arange(-50.,500.,1.)
    x = N.ones(len(y),'f')
    x = 305.*x
    ppgplot.pgsls(2)
    #ppgplot.pgline(x,y)
    ppgplot.pgsci(1)

def plotsfrsigma10all():

    DATAMIN = -.05
    DATAMAX = 50.

    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-.5,18,0)
    ppgplot.pglab("\gSd\5\u (Mpc\u-2\d)","SFR (h\d100\u\u-2\d M\d\(2281) \u yr\u-1\d)","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)
    #x = N.compress((g0.finalsf > 0),g0.dNN)
    x = N.compress((g0.finalsf > 0),g0.sigma10)
    y = N.compress((g0.finalsf > 0),g0.sfr)
    erry = N.compress((g0.finalsf > 0),g0.errsfrc)

    ppgplot.pgpt(x,y,3)
    xs = []
    ys = []
    x1 = x.tolist()
    y1 = y.tolist()
    #for i in range(len(x)):
    #    print x[i],x1[i]
    #drawbinned(x,y,3)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)
    xlabel = 0.
    ylabel = 16.2
    ystep = 1.2
    dy=.4
    dxl=10
    dxr=.10
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ppgplot.pgslw(3)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    x=N.compress((g1.finalsf > 0),g1.sigma10)
    y=N.compress((g1.finalsf > 0),g1.sfrc)
    erry = N.compress((g1.finalsf > 0),g1.errsfrc)
    ppgplot.pgpt(x,y,17)
    x2 = x.tolist()
    y2 = y.tolist()
    #drawbinned(x,y,5)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)

    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)

    ylabel = ylabel - ystep

    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgpt(xlin,ylin,17)
    ppgplot.pgsls(1)
    ppgplot.pgslw(deflw)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    
    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    x=N.compress((g2.finalsf > 0),g2.sigma10)
    y=N.compress((g2.finalsf > 0),g2.sfrc)
    erry = N.compress((g2.finalsf > 0),g2.errsfrc)
    ppgplot.pgpt(x,y,7)
    x3 = x.tolist()
    y3 = y.tolist()
    xs = x1+x2+x3
    ys = y1+y2+y3
    xs = N.array(xs,'f')
    ys = N.array(ys,'f')

    #for i in range(len(xs)):
    #    print i, xs[i],ys[i]
    #print "total should be ",(N.sum(g0.finalsf)+N.sum(g1.finalsf) + N.sum(g2.finalsf)) 
    ppgplot.pgslw(6)
    drawbinned(xs,ys,3)
    print "SFR vs dNN"
    dospear(xs,ys)

    ppgplot.pgslw(deflw)
    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(2)  #line width
    ppgplot.pgpt(xlin,ylin,7)

    y = N.arange(-20.,50.,1.)
    x = N.ones(len(y),'f')
    x = 305.*x
    ppgplot.pgsls(2)
    #ppgplot.pgline(x,y)
    ppgplot.pgsci(1)



def plotewsigma10all():

    DATAMIN = 1.5
    DATAMAX = 4.

    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-20.,200.,0,10)
    #ppgplot.pgenv(DATAMIN,DATAMAX,-10.,28.,0)
    ppgplot.pglab("\gS\d5\u (Mpc\u-2\d)","EW H\ga+[NII] (\(2078))","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #for i in range(len(g0.sigma10)):
    #    print g0.name[i],g0.sigma10[i],g0.sigma10flag[i],g0.vlti[i]
    a=g0.final
    x = N.log10(N.compress((a > 0) & (g0.sigma10flag < 1.),g0.sigma10))
    y = N.compress((a > 0) & (g0.sigma10flag < 1.),g0.ewc)
    erry = N.compress((a > 0) & (g0.sigma10flag < 1.),g0.errewc)
    ppgplot.pgpt(x,y,3)
    errory(x,y,erry)
    x1 = x.tolist()
    y1 = y.tolist()

    #drawbinned(x,y,3)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)
    xlabel = 1.8
    ylabel = 170.
    ystep = 10.
    dy=3.6
    dxl=10
    dxr=.10
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ppgplot.pgslw(3)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    a=g1.final
    x=N.log10(N.compress((a > 0) & (g1.sigma10flag < 1.),g1.sigma10))
    y=N.compress((a > 0) & (g1.sigma10flag < 1.),g1.ewc)
    erry = N.compress((a > 0) & (g1.sigma10flag < 1.),g1.errewc)
    ppgplot.pgpt(x,y,17)
    errory(x,y,erry)
    x2 = x.tolist()
    y2 = y.tolist()

    #drawbinned(x,y,5)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)

    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)

    ylabel = ylabel - ystep

    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgpt(xlin,ylin,17)
    ppgplot.pgsls(1)
    ppgplot.pgslw(deflw)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    
    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    a=g2.final
    x=N.log10(N.compress((a > 0) & (g2.sigma10flag < 1.),g2.sigma10))
    y=N.compress((a > 0) & (g2.sigma10flag < 1.),g2.ewc)
    erry = N.compress((a > 0) & (g2.sigma10flag < 1.),g2.errewc)
    ppgplot.pgpt(x,y,7)
    errory(x,y,erry)
    x3 = x.tolist()
    y3 = y.tolist()
    xs = x1 + x2+ x3
    ys = y1 + y2+ y3
    xs = N.array(xs,'f')
    ys = N.array(ys,'f')
    ppgplot.pgslw(6)
    drawbinned(xs,ys,3)
    x=N.arange(0.,10.,1.)
    y=10*N.ones(len(x),'f')
    ppgplot.pgsls(4)
    ppgplot.pgline(x,y)
    y=-1.*y
    ppgplot.pgline(x,y)
    ppgplot.pgsls(1)
    print "EW vs sigma10"
    dospear(xs,ys)

    #ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(2)  #line width
    ppgplot.pgpt(xlin,ylin,7)

    y = N.arange(-50.,500.,1.)
    x = N.ones(len(y),'f')
    x = 305.*x
    ppgplot.pgsls(2)
    #ppgplot.pgline(x,y)
    ppgplot.pgsci(1)



def plotsffracdr200():
    dmin = -.05
    dmax = 1.
    dstep = 50
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(dmin,dmax,0.,1.,0)
    ppgplot.pglab("d/R\d200\u","Star-Forming Fraction","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)

    x = N.compress((g0.final > 0),g0.dr200)
    y1 = N.compress((g0.final > 0),g0.finalsf)
    y2 = N.compress((g0.final > 0),g0.final)
    (xbin,y1sum)=binitsumequal(x,y1,3)#get total sfr in each mag bin
    (xbin,y2sum)=binitsumequal(x,y2,3)#get total sfr in each mag bin
    sffrac = (y1sum)/(y2sum)
    sffrac = N.array(sffrac,'f')
    ppgplot.pgpt(xbin,sffrac,3)
    ppgplot.pgslw(3)
    ppgplot.pgsls(1)
    ppgplot.pgline(xbin,sffrac)

    xs=[]
    y1s=[]
    y2s=[]
    xs=x.tolist()
    y1s=y1.tolist()
    y2s=y2.tolist()


    x = N.compress((g1.final > 0),g1.dr200)
    y1 = N.compress((g1.final > 0),g1.finalsf)
    y2 = N.compress((g1.final > 0),g1.final)
    (xbin,y1sum)=binitsumequal(x,y1,3)#get total sfr in each mag bin
    (xbin,y2sum)=binitsumequal(x,y2,3)#get total sfr in each mag bin
    sffrac = (y1sum)/(y2sum)
    sffrac = N.array(sffrac,'f')
    ppgplot.pgpt(xbin,sffrac,17)
    ppgplot.pgslw(2)
    ppgplot.pgsls(2)
    ppgplot.pgline(xbin,sffrac)

    xs=xs+x.tolist()
    y1s=y1s+y1.tolist()
    y2s=y2s+y2.tolist()


    x = N.compress((g2.final > 0),g2.dr200)
    y1 = N.compress((g2.final > 0),g2.finalsf)
    y2 = N.compress((g2.final > 0),g2.final)
    (xbin,y1sum)=binitsumequal(x,y1,3)#get total sfr in each mag bin
    (xbin,y2sum)=binitsumequal(x,y2,3)#get total sfr in each mag bin
    sffrac = (y1sum)/(y2sum)
    sffrac = N.array(sffrac,'f')
    ppgplot.pgpt(xbin,sffrac,7)
    ppgplot.pgslw(1)
    ppgplot.pgsls(1)
    ppgplot.pgline(xbin,sffrac)

    xs=xs+x.tolist()
    y1s=y1s+y1.tolist()
    y2s=y2s+y2.tolist()
    x=N.array(xs,'f')
    y1=N.array(y1s,'f')
    y2=N.array(y2s,'f')
    (xbin,y1sum)=binitsumequal(x,y1,3)#get total sfr in each mag bin
    (xbin,y2sum)=binitsumequal(x,y2,3)#get total sfr in each mag bin
    sffrac = (y1sum)/(y2sum)
    sffrac = N.array(sffrac,'f')
    ppgplot.pgsls(1)
    #ppgplot.pgline(xbin,sffrac)
    #for i in range(len(xbin)):
    #    print "hey", i, xbin[i], sffrac[i]
    xlabel = 0.1
    ylabel = .9
    ystep = 0.07
    dy=.015
    dxl=10
    dxr=.05
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,17)

    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,7)
def plotsffracsigma10():
    dmin = 1.5
    dmax = 3.5
    dstep = 50
    nbin=3
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(dmin,dmax,0.,1.,0,10)
    ppgplot.pglab("\gS\d5\u (Mpc\u-2\d)","Star-Forming Fraction","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)


    x = N.log10(N.compress((g0.final > 0) & (g0.sigma10flag < 1.),g0.sigma10))
    y1 = N.compress((g0.final > 0) & (g0.sigma10flag < 1.),g0.finalsf)
    y2 = N.compress((g0.final > 0) & (g0.sigma10flag < 1.),g0.final)
    (xbin,y1sum)=binitsumequal(x,y1,nbin)#get total sfr in each mag bin
    (xbin,y2sum)=binitsumequal(x,y2,nbin)#get total sfr in each mag bin
    (sffrac,sffracerr) = ratioerror(y1sum,y2sum)
    sffrac = N.array(sffrac,'f')
    ppgplot.pgpt(xbin,sffrac,3)
    ppgplot.pgslw(3)
    ppgplot.pgsls(1)
    ppgplot.pgline(xbin,sffrac)
    
    #errory(xbin,sffrac,sffracerr)
    xs=[]
    y1s=[]
    y2s=[]
    xs=x.tolist()
    y1s=y1.tolist()
    y2s=y2.tolist()


    x = N.log10(N.compress((g1.final > 0) & (g1.sigma10flag < 1.),g1.sigma10))
    y1 = N.compress((g1.final > 0) & (g1.sigma10flag < 1.),g1.finalsf)
    y2 = N.compress((g1.final > 0) & (g1.sigma10flag < 1.),g1.final)
    (xbin,y1sum)=binitsumequal(x,y1,nbin)#get total sfr in each mag bin
    (xbin,y2sum)=binitsumequal(x,y2,nbin)#get total sfr in each mag bin
    sffrac = (y1sum)/(y2sum)
    sffrac = N.array(sffrac,'f')
    ppgplot.pgpt(xbin,sffrac,17)
    ppgplot.pgslw(2)
    ppgplot.pgsls(2)
    ppgplot.pgline(xbin,sffrac)
    
    xs=xs+x.tolist()
    y1s=y1s+y1.tolist()
    y2s=y2s+y2.tolist()


    x = N.log10(N.compress((g2.final > 0) & (g2.sigma10flag < 1.),g2.sigma10))
    y1 = N.compress((g2.final > 0) & (g2.sigma10flag < 1.),g2.finalsf)
    y2 = N.compress((g2.final > 0) & (g2.sigma10flag < 1.),g2.final)
    (xbin,y1sum)=binitsumequal(x,y1,nbin)#get total sfr in each mag bin
    (xbin,y2sum)=binitsumequal(x,y2,nbin)#get total sfr in each mag bin
    sffrac = (y1sum)/(y2sum)
    sffrac = N.array(sffrac,'f')
    ppgplot.pgpt(xbin,sffrac,7)
    ppgplot.pgslw(1)
    ppgplot.pgsls(1)
    ppgplot.pgline(xbin,sffrac)

    xs=xs+x.tolist()
    y1s=y1s+y1.tolist()
    y2s=y2s+y2.tolist()
    x=N.array(xs,'f')
    y1=N.array(y1s,'f')
    y2=N.array(y2s,'f')
    (xbin,y1sum)=binitsumequal(x,y1,nbin)#get total sfr in each mag bin
    (xbin,y2sum)=binitsumequal(x,y2,nbin)#get total sfr in each mag bin
    sffrac = (y1sum)/(y2sum)
    sffrac = N.array(sffrac,'f')
    ppgplot.pgsls(1)
    #ppgplot.pgline(xbin,sffrac)
    #for i in range(len(xbin)):
    #    print "hey", i, xbin[i], sffrac[i]
    xlabel = 1.7
    ylabel = .9
    ystep = 0.07
    dy=.015
    dxl=10
    dxr=.05
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,17)

    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,7)

def plotewjms(magj,ew,errew,finalsf,specspir):
    xmin = 16.
    xmax = 23.
    ymin=-50.
    ymax=250.
    z=c.z[ncl]
    DL=dL(z,h100)
    distmod=5*log10(DL*1.e6/10)

    k=kcorr(z,8)#CLJ0023, kcorr J_Rc, band=8 in mystuff.py
    kc=N.average(k[1:3])#take average of Sbc/Scd gal types

    #print "ngal = ",len(g0.sfr)

    ppgplot.pgbeg("ms1054ewj.ps/vcps",1,1)  #color port.
    ppgplot.pgpap(8.,1.25)
    ppgplot.pgpage
    ppgplot.pgsch(1.2) #font size
    ppgplot.pgslw(3)   #line width
    x1=.2
    x2=.45
    x3=.6
    x4=.95
    y1=.2
    y2=.425
    y3=.575
    y4=.85
    xlabel=14.1-14.
    ylabel=1.15
    ppgplot.pgsvp(x1,x4,y1,y4)  #sets viewport
    ppgplot.pgsch(1.3)
    ppgplot.pgslw(5)   #line width
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('bnst',2.,2,'bcvnst',100.,5)  #tickmarks and labeling
    ppgplot.pgmtxt('b',2.5,0.5,0.5,"J (mag-iso)")    #xlabel
    ppgplot.pgmtxt('l',3.1,0.5,0.5,"EW H\ga+[NII] (\(2078))")

    ppgplot.pgsch(3.)    
    #ppgplot.pgsls(1)#dotted
    #ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)
    ppgplot.pgpt(magj,ew,1)
    ppgplot.pgsch(2.)    
    x=N.compress((finalsf > 0),magj)
    #ppgplot.pgpt(MR,ew,1)
    #x=N.compress((finalsf > 0),MR)
    y=N.compress((finalsf > 0),ew)
    erry = N.compress((finalsf > 0),errew)
    errory(x,y,erry)
    ppgplot.pgpt(x,y,17)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)

    x=N.compress(((specspir > 0) & (finalsf > 0)),magj)
    y=N.compress(((specspir > 0) & (finalsf > 0)),ew)
    ppgplot.pgsci(4)
    ppgplot.pgpt(x,y,18)

    x=N.compress(((specspir > 0) & (finalsf < 1)),magj)
    y=N.compress(((specspir > 0) & (finalsf < 1)),ew)
    erry=N.compress(((specspir > 0) & (finalsf < 1)),errew)
    ppgplot.pgsci(8)
    errory(x,y,erry)
    ppgplot.pgpt(x,y,18)
    ppgplot.pgsch(1.2)
    ppgplot.pgsci(1)
    ZP = 22.04
    sfrconv = c.sfrconv[ncl]
    sfrmin = 3.
    mag = N.arange(14.,27,.05)
    y = N.ones(len(mag),'f')
    y = ewmin*y
    ppgplot.pgsls(4)
    ppgplot.pgline(mag,y)
    y = 40*N.ones(len(mag),'f')
    ppgplot.pgline(mag,y)
    JK = 0
    ppgplot.pgsls(2)
    #ppgplot.pgline(mag,ewlimit)
    ppgplot.pgsci(1)
    #ppgplot.pgend()

    bins=N.arange(16,25,.1)
    errfj = 92.*N.sqrt(10**((22.-bins)/2.5))#fit min iso-area versus mag
    errfj = N.sqrt(errfj)*c.jnoisea[ncl]*(1.+c.jnoiseb[ncl]*N.sqrt(errfj))#use CL1216 noise
    ewlimit=2.4*errfj/sfrconv/h100/10**((ZP-bins-JK)/2.5)*dlambdaJ
    #sfrlimit=.24*errfj
    ppgplot.pgsci(1)
    ppgplot.pgsls(2)
    #ppgplot.pgline(mag,sfrlimit)
    ppgplot.pgline(bins,ewlimit)
    #ewlimit=-1*ewlimit
    #ppgplot.pgline(bins,ewlimit)
    ppgplot.pgsci(1)
    ppgplot.pgsls(1)   #line width
    ppgplot.pgsch(1.3)
    ppgplot.pgslw(5)   #line width
    
    xmin=xmin -distmod -kc
    xmax=xmax -distmod -kc
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('cmst',2.,2,'',100.,5)  #tickmarks and labeling
    #ppgplot.pgmtxt('t',2.4,0.5,0.5,"M\dR\u + 5 log h\d100\u")    #xlabel
    ppgplot.pgmtxt('t',2.4,0.5,0.5,"M\dR\u")    #xlabel
    ppgplot.pgend()

    #os.system("cp ms1054ewj.ps /Users/rfinn/clusters/papers/ms1054/.")

def plotsfrjms(magj,sfr,errsfr,finalsf,specspir):
    xmin = 16.
    xmax = 23.
    ymin=-20.
    ymax=140.
    z=c.z[ncl]
    DL=dL(z,h100)
    distmod=5*log10(DL*1.e6/10)

    k=kcorr(z,8)#CLJ0023, kcorr J_Rc, band=8 in mystuff.py
    kc=N.average(k[1:3])#take average of Sbc/Scd gal types

    #print "ngal = ",len(g0.sfr)

    ppgplot.pgbeg("ms1054sfrj.ps/vcps",1,1)  #color port.
    ppgplot.pgpap(8.,1.25)
    ppgplot.pgpage
    ppgplot.pgsch(1.2) #font size
    ppgplot.pgslw(3)   #line width
    x1=.2
    x2=.45
    x3=.6
    x4=.95
    y1=.2
    y2=.425
    y3=.575
    y4=.85
    xlabel=14.1-14.
    ylabel=1.15
    ppgplot.pgsvp(x1,x4,y1,y4)  #sets viewport
    ppgplot.pgsch(1.3)
    ppgplot.pgslw(5)   #line width
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('bnst',2.,2,'bcvnst',20.,2)  #tickmarks and labeling
    ppgplot.pgmtxt('b',2.5,0.5,0.5,"J (mag-iso)")    #xlabel
    ppgplot.pgmtxt('l',3.1,0.5,0.5,sfrlabel70)

    ppgplot.pgsch(3.)
    #ppgplot.pgsls(1)#dotted
    #ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)
    ppgplot.pgpt(magj,sfr,1)
    ppgplot.pgsch(2.)
    x=N.compress((finalsf > 0),magj)
    #ppgplot.pgpt(MR,ew,1)
    #x=N.compress((finalsf > 0),MR)
    y=N.compress((finalsf > 0),sfr)
    erry = N.compress((finalsf > 0),errsfr)
    errory(x,y,erry)
    ppgplot.pgpt(x,y,17)
    #ppgplot.pgerry(len(x),x,y+erry,y-erry,.1)
    #ppgplot.pgerr1(2,x,y,erry,.1)

    x=N.compress(((specspir > 0) & (finalsf > 0)),magj)
    y=N.compress(((specspir > 0) & (finalsf > 0)),sfr)
    ppgplot.pgsci(4)
    ppgplot.pgpt(x,y,18)

    x=N.compress(((specspir > 0) & (finalsf < 1)),magj)
    y=N.compress(((specspir > 0) & (finalsf < 1)),sfr)
    erry=N.compress(((specspir > 0) & (finalsf < 1)),errsfr)
    ppgplot.pgsci(8)
    errory(x,y,erry)
    ppgplot.pgpt(x,y,18)

    ppgplot.pgsch(1.2) #font size
    ppgplot.pgsci(1)
    ZP = 22.
    sfrconv = c.sfrconv[ncl]
    mag = N.arange(14.,28,.05)
    #set JK=(23-mag)*1.2/5
    #JK = 1.2*mag/mag
    sfrlimit = 10**((ZP-mag)/2.5)*ewmin/dlambdaJ*sfrconv/h100**2
    ppgplot.pgsls(4)
    ppgplot.pgline(mag,sfrlimit)

    area=g5.isoarea
    bins=N.arange(16,25,.1)
    (magbin,areabin)=binitmin(magj,area,15)
    errfj = 2.*N.sqrt(areabin)*c.jnoisea[ncl]*(1.+c.jnoiseb[ncl]*N.sqrt(areabin))#use CL1216 noise
    errfj = N.sqrt(10**((ZP-bins)/2.5))


    #sfrlimit=0.08*c.sfrconv[ncl]*errfj
    sfrlimit=0.055*c.sfrconv[ncl]*errfj

    ppgplot.pgsls(2)
    #ppgplot.pgline(mag,sfrlimit)
    ppgplot.pgline(bins,sfrlimit)
    #sfrlimit=-1*sfrlimit
    #ppgplot.pgline(bins,sfrlimit)
    ppgplot.pgsci(1)
    ppgplot.pgsls(1)
    y=0.*bins
    #ppgplot.pgline(bins,y)
    ppgplot.pgsci(1)
    ppgplot.pgsls(1)   #line width
    ppgplot.pgsch(1.3)
    ppgplot.pgslw(5)   #line width
    
    xmin=xmin -distmod -kc
    xmax=xmax -distmod -kc
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('cmst',2.,2,'',100.,5)  #tickmarks and labeling
    #ppgplot.pgmtxt('t',2.4,0.5,0.5,"M\dR\u + 5 log h\d100\u")    #xlabel
    ppgplot.pgmtxt('t',2.4,0.5,0.5,"M\dR\u")    #xlabel
    ppgplot.pgend()

    os.system("cp ms1054sfrj.ps /Users/rfinn/clusters/papers/ms1054/.")
def plotnsfsim(): #simulation to estimate errors for NSF
    ngal=10000
    nobs=1000
    boost=1.20
    spreadboost=.1
    psplotinit("nsfsim.ps")
    DATAMIN = 18.
    DATAMAX = 23.
    bins=N.arange(DATAMIN,DATAMAX+1,1.,'f')
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-.5,10,0)
    ppgplot.pglab("J (mag-iso)","SFR (h\d100\u\u-2\d M\d\(2281) \u yr\u-1\d)","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    x0=N.compress((g0.finalsf > 0),g0.magj)
    y0=N.compress((g0.finalsf > 0),g0.sfrc)
    erry0 = N.compress((g0.finalsf > 0),g0.errsfrc)
    x1=N.compress((g1.finalsf > 0),g1.magj)
    y1=N.compress((g1.finalsf > 0),g1.sfrc)
    erry1 = N.compress((g1.finalsf > 0),g1.errsfrc)
    x2=N.compress((g2.finalsf > 0),g2.magj)
    y2=N.compress((g2.finalsf > 0),g2.sfrc)
    erry2 = N.compress((g2.finalsf > 0),g2.errsfrc)

    mag = x0.tolist() + x1.tolist() + x2.tolist()
    sfr = y0.tolist() + y1.tolist() + y2.tolist()

    mag=N.array(mag,'f')
    sfr=N.array(sfr,'f')
    mag=N.compress(sfr > 1,mag)
    sfr=N.compress(sfr > 1,sfr)

    #generate random sample w/1000 galaxies
    magsim=N.zeros(ngal,'f')
    sfrsim=N.zeros(ngal,'f')
    magsimfield=N.zeros(ngal,'f')
    sfrsimfield=N.zeros(ngal,'f')


    gal=N.arange(0,len(mag),1,'i')
    for i in range(ngal):
        c=int(random.choice(gal))
        magsim[i]=mag[c]
        sfrsim[i]=sfr[c]
        c=int(random.choice(gal))
        magsimfield[i]=mag[c]
        sfrsimfield[i]=random.gauss(boost,spreadboost)*sfr[c]
    avecluster=N.average(sfrsim)
    erravecluster=pylab.std(sfrsim)/N.sqrt(1.*ngal)
    avefield=N.average(sfrsimfield)
    erravefield=pylab.std(sfrsimfield)/N.sqrt(1.*ngal)
    print "NSF simulation"
    print "cluster, field ave sfr = %5.3f +/- %5.3f, %5.3f +/- %5.3f" %(avecluster,erravecluster,avefield,erravefield)
    print "cluster+3sigma, field-3sigma = %5.3f , %5.3f " %(avecluster+3*erravecluster,avefield-3*erravefield)
    #(med,errmin,errmax)=bootstrap(sfrsimfield)
    #print "bootstrap cluster med %5.3f + %5.3f - %5.3f" % (float(med),float(errmin),float(errmax))
    #(med,errmin,errmax)=bootstrap(sfrsim)
    #print "bootstrap field med %5.3f + %5.3f - %5.3f" % (float(med),float(errmin),float(errmax))


    ybin=N.zeros((len(bins)-1),'f')
    magbin=N.zeros((len(bins)-1),'f')
    errybin=N.zeros((len(bins)-1),'f')
    ybinfield=N.zeros((len(bins)-1),'f')
    errybinfield=N.zeros((len(bins)-1),'f')
    for i in range(len(bins)-1):
        magbin[i]=(bins[i]+bins[i+1])/2.
        #print i,len(bins),bins[i],bins[i+1],magbin[i]
        sfrbin=N.compress((magsim > bins[i]) & (magsim < bins[i+1]),sfrsim)
        ybin[i] = N.average(sfrbin)
        errybin[i]= pylab.std(sfrbin)/N.sqrt(1.*len(sfrbin)*nobs/ngal) 
        sfrbin=N.compress((magsimfield > bins[i]) & (magsimfield < bins[i+1]),sfrsimfield)
        ybinfield[i] = N.average(sfrbin)
        errybinfield[i]= pylab.std(sfrbin)/N.sqrt(1.*len(sfrbin)*nobs/ngal)  

    ppgplot.pgpt(magbin,ybin,-3)
    errory(magbin,ybin,1.*errybin)

    ppgplot.pgpt(magbin,ybinfield,6)
    errory(magbin,ybinfield,1.*errybinfield)

    ppgplot.pgend()

def plotnsfsim2(): #simulation of factor of 2 increase in EW for 20% of field population to estimate errors for NSF
    ngal=10000
    nobs=1000
    boost=2.
    spreadboost=.1
    mergerfrac = .22
    psplotinit("nsfsim2.ps")
    DATAMIN = 18.
    DATAMAX = 23.
    bins=N.arange(DATAMIN,DATAMAX+1,1.,'f')
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-5,150,0)
    ppgplot.pglab("J (mag-iso)","EW H\ga+[NII] (\(2078))","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    x0=N.compress((g0.finalsf > 0),g0.magj)
    ew0=N.compress((g0.finalsf > 0),g0.ewc)
    erry0 = N.compress((g0.finalsf > 0),g0.errewc)
    y0=N.compress((g0.finalsf > 0),g0.sfrc)
    x1=N.compress((g1.finalsf > 0),g1.magj)
    ew1=N.compress((g1.finalsf > 0),g1.ewc)
    erry1 = N.compress((g1.finalsf > 0),g1.errewc)
    y1=N.compress((g1.finalsf > 0),g1.sfrc)

    x2=N.compress((g2.finalsf > 0),g2.magj)
    ew2=N.compress((g2.finalsf > 0),g2.ewc)
    erry2 = N.compress((g2.finalsf > 0),g2.errewc)
    y2=N.compress((g2.finalsf > 0),g2.sfrc)
    mag = x0.tolist() + x1.tolist() + x2.tolist()
    sfr = y0.tolist() + y1.tolist() + y2.tolist()
    ew = ew0.tolist() + ew1.tolist() + ew2.tolist()
    mag=N.array(mag,'f')
    sfr=N.array(sfr,'f')
    ew = N.array(ew,'f')
    mag=N.compress(sfr > 1,mag)
    ew=N.compress(sfr > 1,ew)

    #generate random sample w/1000 galaxies
    magsim=N.zeros(ngal,'f')
    sfrsim=N.zeros(ngal,'f')
    magsimfield=N.zeros(ngal,'f')
    sfrsimfield=N.zeros(ngal,'f')


    gal=N.arange(0,len(mag),1,'i')
    for i in range(ngal):
        c=int(random.choice(gal))
        magsim[i]=mag[c]
        sfrsim[i]=ew[c]
        c=int(random.choice(gal))
        magsimfield[i]=mag[c]
        if random.uniform(0,1) < mergerfrac:
            sfrsimfield[i]=random.gauss(boost,spreadboost)*ew[c]
        else:
            sfrsimfield[i]=ew[c]
    avecluster=N.average(sfrsim)
    erravecluster=pylab.std(sfrsim)/N.sqrt(1.*ngal)
    avefield=N.average(sfrsimfield)
    erravefield=pylab.std(sfrsimfield)/N.sqrt(1.*ngal)
    print "NSF simulation 2"
    print "cluster, field ave ew = %5.3f +/- %5.3f, %5.3f +/- %5.3f" %(avecluster,erravecluster,avefield,erravefield)
    print "cluster+3sigma, field-3sigma = %5.3f , %5.3f " %(avecluster+3*erravecluster,avefield-3*erravefield)
    ybin=N.zeros((len(bins)-1),'f')
    magbin=N.zeros((len(bins)-1),'f')
    errybin=N.zeros((len(bins)-1),'f')
    ybinfield=N.zeros((len(bins)-1),'f')
    errybinfield=N.zeros((len(bins)-1),'f')
    for i in range(len(bins)-1):
        magbin[i]=(bins[i]+bins[i+1])/2.
        print i,len(bins),bins[i],bins[i+1],magbin[i]
        sfrbin=N.compress((magsim > bins[i]) & (magsim < bins[i+1]),sfrsim)
        ybin[i] = N.average(sfrbin)
        errybin[i]= pylab.std(sfrbin)/N.sqrt(1.*len(sfrbin)*nobs/ngal) 
        sfrbin=N.compress((magsimfield > bins[i]) & (magsimfield < bins[i+1]),sfrsimfield)
        ybinfield[i] = N.average(sfrbin)
        errybinfield[i]= pylab.std(sfrbin)/N.sqrt(1.*len(sfrbin)*nobs/ngal)  

    ppgplot.pgpt(magbin,ybin,-3)
    errory(magbin,ybin,1.*errybin)

    ppgplot.pgpt(magbin,ybinfield,6)
    errory(magbin,ybinfield,1.*errybinfield)

    ppgplot.pgend()

def plotnsfsim3(): #simulation to estimate errors for NSF, morphology, strangulation
    ngal=10000
    nobs=1000
    boost=1.2
    spreadboost=.1
    mergerfrac = .22
    psplotinit("nsfsim3.ps")
    DATAMIN = -0.05
    DATAMAX = 1.05
    bins=N.array([0.,0.3,0.6,1.0],'f')
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-.2,8,0)
    #ppgplot.pglab("B/T","EW H\ga+[NII] (\(2078))","")
    ppgplot.pglab("B/T","SFR (h\d100\u\u-2\d M\d\(2281) \u yr\u-1\d)","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    x0=N.compress((g0.finalsf > 0),g0.magj)
    ew0=N.compress((g0.finalsf > 0),g1.ewc)
    morph0=N.compress((g0.finalsf > 0),g0.BTiband)
    erry0 = N.compress((g0.finalsf > 0),g0.errewc)
    y0=N.compress((g0.finalsf > 0),g0.sfrc)

    x1=N.compress((g1.finalsf > 0),g1.magj)
    ew1=N.compress((g1.finalsf > 0),g1.ewc)
    morph1=N.compress((g1.finalsf > 0),g1.BTiband)
    erry1 = N.compress((g1.finalsf > 0),g1.errewc)
    y1=N.compress((g1.finalsf > 0),g1.sfrc)

    x2=N.compress((g2.finalsf > 0),g2.magj)
    ew2=N.compress((g2.finalsf > 0),g2.ewc)
    morph2=N.compress((g2.finalsf > 0),g2.BTiband)
    erry2 = N.compress((g2.finalsf > 0),g2.errewc)
    y2=N.compress((g2.finalsf > 0),g2.sfrc)
    mag = x0.tolist() + x1.tolist() + x2.tolist()
    sfr = y0.tolist() + y1.tolist() + y2.tolist()
    ew = ew0.tolist() + ew1.tolist() + ew2.tolist()
    morph = morph0.tolist() + morph1.tolist() + morph2.tolist()
    mag=N.array(mag,'f')
    sfr=N.array(sfr,'f')
    ew = N.array(ew,'f')
    morph = N.array(morph,'f')
    mag=N.compress(sfr > 1,mag)
    ew=N.compress(sfr > 1,ew)
    morph=N.compress(sfr > 1,morph)
    sfr=N.compress(sfr > 1,sfr)

    #ppgplot.pgpt(morph,ew,3)
    #drawbinned(morph,ew,3)
    #generate random sample w/1000 galaxies
    xsim=N.zeros(ngal,'f')
    ysim=N.zeros(ngal,'f')
    xsimfield=N.zeros(ngal,'f')
    ysimfield=N.zeros(ngal,'f')


    gal=N.arange(0,len(mag),1,'i')
    for i in range(ngal):
        c=int(random.choice(gal))
        xsim[i]=morph[c]
        ysim[i]=sfr[c]
        c=int(random.choice(gal))
        xsimfield[i]=morph[c]
        ysimfield[i]=random.gauss(boost,spreadboost)*sfr[c]
    avecluster=N.average(ysim)
    erravecluster=pylab.std(ysim)/N.sqrt(1.*nobs)
    avefield=N.average(ysimfield)
    erravefield=pylab.std(ysimfield)/N.sqrt(1.*nobs)
    print "NSF simulation 3"
    print "cluster, field ave ew = %5.3f +/- %5.3f, %5.3f +/- %5.3f" %(avecluster,erravecluster,avefield,erravefield)
    print "cluster+3sigma, field-3sigma = %5.3f , %5.3f " %(avecluster+3*erravecluster,avefield-3*erravefield)
    ybin=N.zeros((len(bins)-1),'f')
    magbin=N.zeros((len(bins)-1),'f')
    errybin=N.zeros((len(bins)-1),'f')
    errypbin=N.zeros((len(bins)-1),'f')
    errymbin=N.zeros((len(bins)-1),'f')
    ybinfield=N.zeros((len(bins)-1),'f')
    errybinfield=N.zeros((len(bins)-1),'f')
    errypbinfield=N.zeros((len(bins)-1),'f')
    errymbinfield=N.zeros((len(bins)-1),'f')
    for i in range(len(bins)-1):
        magbin[i]=(bins[i]+bins[i+1])/2.
        #print i,len(bins),bins[i],bins[i+1],magbin[i]
        sfrbin=N.compress((xsim > bins[i]) & (xsim < bins[i+1]),ysim)
        #(med,errm,errp)=bootstrap(sfrbin)
        #ybin[i]=float(med)
        #errymbin[i]=float(errm)
        #errypbin[i]=float(errp)

        ybin[i] = N.average(sfrbin)
        errybin[i]= pylab.std(sfrbin)/N.sqrt(1.*len(sfrbin)*nobs/ngal) 
        sfrbin=N.compress((xsimfield >= bins[i]) & (xsimfield < bins[i+1]),ysimfield)
        #(med,errm,errp)=bootstrap(sfrbin)
        ybinfield[i] = N.average(sfrbin)
        errybinfield[i]= pylab.std(sfrbin)/N.sqrt(1.*len(sfrbin)*nobs/ngal)  
        #ybinfield[i]=float(med)
        #errymbinfield[i]=float(errm)
        #errypbinfield[i]=float(errp)
    ppgplot.pgpt(magbin,ybin,-3)
    errory(magbin,ybin,1.*errybin)
    #errorym(magbin,ybin,1.*errymbin)
    #erroryp(magbin,ybin,1.*errypbin)

    ppgplot.pgpt(magbin,ybinfield,6)
    errory(magbin,ybinfield,1.*errybinfield)
    #erroryp(magbin,ybinfield,1.*errypbinfield)
    #errorym(magbin,ybinfield,1.*errymbinfield)
    
    ppgplot.pgend()
def plotnsfsim4(): #simulation to estimate errors for NSF, morphology, strangulation
    ngal=10000
    nobs=500
    boost=2.
    spreadboost=.1
    mergerfrac = .22
    psplotinit("nsfsim4.ps")
    DATAMIN = -0.05
    DATAMAX = 1.05
    bins=N.array([0.,0.3,0.6,1.0],'f')
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,-5,80,0)
    ppgplot.pglab("B/T","EW H\ga+[NII] (\(2078))","")
    #ppgplot.pglab("B/T","SFR (h\d100\u\u-2\d M\d\(2281) \u yr\u-1\d)","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    x0=N.compress((g0.finalsf > 0),g0.magj)
    ew0=N.compress((g0.finalsf > 0),g1.ewc)
    morph0=N.compress((g0.finalsf > 0),g0.BTiband)
    erry0 = N.compress((g0.finalsf > 0),g0.errewc)
    y0=N.compress((g0.finalsf > 0),g0.sfrc)

    x1=N.compress((g1.finalsf > 0),g1.magj)
    ew1=N.compress((g1.finalsf > 0),g1.ewc)
    morph1=N.compress((g1.finalsf > 0),g1.BTiband)
    erry1 = N.compress((g1.finalsf > 0),g1.errewc)
    y1=N.compress((g1.finalsf > 0),g1.sfrc)

    x2=N.compress((g2.finalsf > 0),g2.magj)
    ew2=N.compress((g2.finalsf > 0),g2.ewc)
    morph2=N.compress((g2.finalsf > 0),g2.BTiband)
    erry2 = N.compress((g2.finalsf > 0),g2.errewc)
    y2=N.compress((g2.finalsf > 0),g2.sfrc)
    mag = x0.tolist() + x1.tolist() + x2.tolist()
    sfr = y0.tolist() + y1.tolist() + y2.tolist()
    ew = ew0.tolist() + ew1.tolist() + ew2.tolist()
    morph = morph0.tolist() + morph1.tolist() + morph2.tolist()
    mag=N.array(mag,'f')
    sfr=N.array(sfr,'f')
    ew = N.array(ew,'f')
    morph = N.array(morph,'f')
    mag=N.compress(sfr > 1,mag)
    ew=N.compress(sfr > 1,ew)
    morph=N.compress(sfr > 1,morph)
    sfr=N.compress(sfr > 1,sfr)

    #ppgplot.pgpt(morph,ew,3)
    #drawbinned(morph,ew,3)
    #generate random sample w/1000 galaxies
    xsim=N.zeros(ngal,'f')
    ysim=N.zeros(ngal,'f')
    xsimfield=N.zeros(ngal,'f')
    ysimfield=N.zeros(ngal,'f')


    gal=N.arange(0,len(mag),1,'i')
    for i in range(ngal):
        c=int(random.choice(gal))
        xsim[i]=morph[c]
        ysim[i]=ew[c]
        c=int(random.choice(gal))
        xsimfield[i]=morph[c]
        if random.uniform(0,1) < mergerfrac:
            ysimfield[i]=random.gauss(boost,spreadboost)*ew[c]
        else:
            ysimfield[i]=ew[c]
    avecluster=N.average(ysim)
    erravecluster=pylab.std(ysim)/N.sqrt(1.*nobs)
    avefield=N.average(ysimfield)
    erravefield=pylab.std(ysimfield)/N.sqrt(1.*nobs)
    print "NSF simulation 3"
    print "cluster, field ave ew = %5.3f +/- %5.3f, %5.3f +/- %5.3f" %(avecluster,erravecluster,avefield,erravefield)
    print "cluster+3sigma, field-3sigma = %5.3f , %5.3f " %(avecluster+3*erravecluster,avefield-3*erravefield)
    ybin=N.zeros((len(bins)-1),'f')
    magbin=N.zeros((len(bins)-1),'f')
    errybin=N.zeros((len(bins)-1),'f')
    errypbin=N.zeros((len(bins)-1),'f')
    errymbin=N.zeros((len(bins)-1),'f')
    ybinfield=N.zeros((len(bins)-1),'f')
    errybinfield=N.zeros((len(bins)-1),'f')
    errypbinfield=N.zeros((len(bins)-1),'f')
    errymbinfield=N.zeros((len(bins)-1),'f')
    for i in range(len(bins)-1):
        magbin[i]=(bins[i]+bins[i+1])/2.
        #print i,len(bins),bins[i],bins[i+1],magbin[i]
        sfrbin=N.compress((xsim > bins[i]) & (xsim < bins[i+1]),ysim)
        ybin[i] = N.average(sfrbin)
        errybin[i]= pylab.std(sfrbin)/N.sqrt(1.*len(sfrbin)*nobs/ngal) 
        sfrbin=N.compress((xsimfield >= bins[i]) & (xsimfield < bins[i+1]),ysimfield)
        ybinfield[i] = N.average(sfrbin)
        errybinfield[i]= pylab.std(sfrbin)/N.sqrt(1.*len(sfrbin)*nobs/ngal)  
    ppgplot.pgpt(magbin,ybin,-3)
    errory(magbin,ybin,1.*errybin)

    ppgplot.pgpt(magbin,ybinfield,6)
    errory(magbin,ybinfield,1.*errybinfield)
    
    ppgplot.pgend()

def plotsfrlit():

    #print "sfrlit"
    sumsfra = g0.totalsfrr200#/cl1040compl
    sumsfrb = g1.totalsfrr200#/cl1054compl
    sumsfrc = g2.totalsfrr200#/cl1216compl 
    errsumsfra = g0.totalsfrr200err#/cl1040compl
    errsumsfrb = g1.totalsfrr200err#/cl1054compl
    errsumsfrc = g2.totalsfrr200err#/cl1216compl 
    #errsumsfra = 0.3*sumsfra
    #errsumsfrb = 0.3*sumsfrb
    #errsumsfrc = 0.3*sumsfrc

    #print sumsfra,g0.totalsfr,cl1040compl
    #print sumsfrb,g1.totalsfr,cl1054compl
    #print sumsfrc,g2.totalsfr,cl1216compl
    #with RXJ0152
    #sumsfrhiz=N.array([sumsfra,sumsfrb,sumsfrc,28.6,40.1,47,7.7,10.7,(253.*(70./100.)**2)],'f')
    #sfrcor=N.array([1.,1,1,1,1,1,2.8,2.8,1.],'f')#aperture correction
    #volcor=N.array([1.,(1./.94),(1./.86),1,1,1,(1./.99),(1./.74),1.],'f')#aperture correction
    #zhiz=N.array([.704,.748,.794,.845,0.833,.228,.32,.183,0.39],'f')
    #symbols = [-4,-4,-4,-3,-16,6,4,7,11]
    #ref = ["This study","Finn et al. 2004","Finn et al. 2005","Balogh & Morris 2000","Couch et al. 2001", "Balogh et al. 2002","Kodama et al. 2004"]
    #refsymb = [-4,-3,-16,6, 4,7,11]
    #sigmahiz=N.array([418,504,1018,415,1000,1023,1390,1274,560],'f')
    #kodamamasshiz=N.array([.8,1.4,11.4,(2.3*.7),(11.3*.7),(13.6*.7),(.7*7.3),(8.1*.7),(5.7*.7)],'f')#lensing masses

    #w/out RXJ0152
    sumsfrhiz=N.array([sumsfra,sumsfrb,sumsfrc,38.2/h100**2,46.98/h100**2,7.7/h100**2,10.7/h100**2,(253.*(70./100./h100)**2)],'f')
    errsumsfrhiz=N.array([errsumsfra,errsumsfrb,errsumsfrc,5,5,5,5,17.],'f')
    sfrcor=N.array([1.,1,1,1,1,(2.8*1),(2.8*1),1.],'f')#aperture correction
    volcor=N.array([1.,(1./.94),(1./.86),1,1,(1./.99),(1./.74),1.],'f')#volume correction

    ewcor=N.array([1.,1.,1.,1.,1.7,1.,1.,1.],'f')#ew correction
    fracvol=("1.00","0.94","0.86","1.00","1.00","0.99","0.74","1.00")#aperture correction
    zhiz=N.array([.704,.748,.794,.845,.228,.32,.183,0.395],'f')
    refsymb = [-4,-5,-6,-3,6, 4,7,11,3]
    sigmahiz=N.array([418,504,1018,415,1023,1390,1274,561],'f')
    errpsigmahiz=N.array([55,113,73,102,102,139,127,95],'f')
    errmsigmahiz=N.array([46,65,77,63,102,139,127,83],'f')
    errsigmahiz=0.5*(errpsigmahiz+errmsigmahiz)
    da=N.array([.704,.748,.794,.845,.228,.32,.183,0.39],'f')
    r200arcmin=N.array([c.r200pix[0]*c.pscale[0]/60.,c.r200pix[1]*c.pscale[1]/60.,c.r200pix[2]*c.pscale[2]/60.,c.r200pix[3]*c.pscale[3]/60.,10.29,9.54,16.99,99.],'f')
    r200Mpc=N.array([c.r200Mpc[0],c.r200Mpc[1],c.r200Mpc[2],c.r200Mpc[3],1.58,1.87,2.20,99.],'f')
    surveyrarcmin = N.array([1.19,1.31,1.15,1.29,8.00,4.35,4.35,15.],'f')
    surveyr200 =    N.array([0.68,0.41,0.32,0.96,0.78,0.46,0.27,1.00],'f')
    for i in range(len(r200Mpc)):
        r200Mpc[i]=r200(sigmahiz[i],zhiz[i],1)
        r200arcmin[i]=r200Mpc[i]*1000./DA(zhiz[i],1.)/60.
        surveyr200[i]=surveyrarcmin[i]/r200arcmin[i]
    clname = ["\ca","\cb","\cc","\cj", "Abell~2390","AC~114","Abell~1689","CL0024.0$+$1652"]
    refnumber = ["1","1","1","2", "3","4","5","6"]
    technique = ["I","I","I","I", "I","S","S","I"]
    symbols = [-4,-5,-6,-3,6,4,7,11]
    #ref = ["This study","Finn et al. 2004","Balogh & Morris 2000","Couch et al. 2001", "Balogh et al. 2002","Kodama et al. 2004","Alternate M\dcl\u"]
    ref = ["Finn et al. 2005","Finn et al. 2004","Balogh & Morris 2000","Couch et al. 2001", "Balogh et al. 2002","Kodama et al. 2004","Alternate M\dcl\u"]
    ref = ["CL1040-1155","CL1054-1245","CL1216-1201","CLJ0023+0423B","Abell 2390","AC 114", "Abell 1689","CL 0024.0+1652","Alternate M\dcl\u"]
    
    kodamamasshiz=N.array([-1.,-1.,-1.,(2.3/.7),(13.6/.7),(7.3/.7),(8.1/.7),(5.7/.7)],'f')#lensing masses

    sumsfrhiz=sumsfrhiz*sfrcor*volcor*ewcor
    totalcor=sfrcor*volcor*ewcor
    errsumsfrhiz=errsumsfrhiz*sfrcor*volcor*ewcor
    #sfrmhiz=N.log10(N.array([70,20,26,58,5.8,4,.3,.8],'f'))
    masshiz=kodamamasshiz
    specmass=[0,1,2]
    #masshiz=9.78*(sigmahiz/1000.)**3.*1./N.sqrt(.7+.3*(1+zhiz)**3.)
    masshiz=12.*(sigmahiz/1000.)**3.*1./N.sqrt(.7+.3*(1+zhiz)**3.)
    #for i in specmass:
    #    masshiz[i]=12*(sigmahiz[i]/1000.)**3.*1./N.sqrt(.7+.3*(1+zhiz[i]))
    #mratio=masshiz/kodamamasshiz
    #sfrmhiz=(sumsfrhiz/kodamamasshiz)*(.64)#correct for AGN contamination
    #sfrmhiz=(sumsfrhiz/kodamamasshiz)
    sfrmhiz=(sumsfrhiz/masshiz)
    errsfrmhiz=(errsumsfrhiz/masshiz)
    sfrmkodama=(sumsfrhiz/kodamamasshiz)
    errsfrmkodama=(errsumsfrhiz/kodamamasshiz)
    x1=sfrmhiz[4:]
    x2=sfrmkodama[4:]
    loz1=x1.tolist()+x2.tolist()
    #loz2=N.average(sfrmhiz[4:]/kodamamasshiz[4:])
    x1=sfrmhiz[0:3]
    x2=sfrmkodama[3]
    hiz1=x1.tolist()
    hiz1.append(x2)
    sfrmall=hiz1+loz1
    sfrmall=N.array(sfrmall,'f')
    #print "length of sfrmall = ",len(sfrmall)
    zall=     N.zeros(len(sfrmall),'f')
    symbolall=N.zeros(len(sfrmall),'i')
    massall=  N.zeros(len(sfrmall),'f')
    sfrall=   N.zeros(len(sfrmall),'f')
    zall=     []
    symbolall=[]
    massall=  []
    sfrall=   []
    for i in range(len(zhiz)):
        zall.append(float(zhiz[i]))
        symbolall.append(symbols[i])
        massall.append(float(masshiz[i]))
        sfrall.append(float(sumsfrhiz[i]))
    #for i in range(len(kodamamasshiz)):
    #    if (kodamamasshiz[i] > 0):
    #        zall.append(float(zhiz[i]))
    #        symbolall.append(3)
    #        massall.append(float(kodamamasshiz[i]))
    #        sfrall.append(float(sumsfrhiz[i]))
    zall=     N.array(zall,'f')
    symbolall=N.array(symbolall,'i')
    massall=  N.array(massall,'f')
    sfrall=   N.array(sfrall,'f')

    sfrmall=sfrall/massall
    loz1=N.array(loz1,'f')
    hiz1=N.array(hiz1,'f')
    lozave=N.average(loz1)

    hizave=N.average(hiz1)

    lozsfrm=N.average(N.compress(zall < 0.5,sfrmall))
    lozmass=N.average(N.compress(zall < 0.5,massall))
    loz=N.average(N.compress(zall < 0.5,zall))
    hizsfrm=N.average(N.compress(zall > 0.5,sfrmall))
    hizmass=N.average(N.compress(zall > 0.5,massall))
    hiz=N.average(N.compress(zall > 0.5,zall))

    lozsfrm=pylab.median(N.compress(zall < 0.5,sfrmall))
    lozmass=pylab.median(N.compress(zall < 0.5,massall))
    hizsfrm=pylab.median(N.compress(zall > 0.5,sfrmall))
    hizmass=pylab.median(N.compress(zall > 0.5,massall))

    aveincrease = (hizsfrm/lozsfrm)
    print "Ave increase in SFR/Mcl from loz to hiz = ",aveincrease,"from ",lozsfrm," to ",hizsfrm
    print "for average z of ",loz, " to ",hiz

    print "loz sfrm values = ",N.compress(zall < 0.5,sfrmall)
    m=(hizsfrm-lozsfrm)/(hiz-loz)
    b=hizsfrm - m*(hiz)
    
    #predictedsfrm=m*zhiz+b
    #dev=sfrmhiz-predictedsfrm
    predictedsfrm=m*zall+b
    dev=sfrmall-predictedsfrm


    m=(hizsfrm-lozsfrm)/(hizmass-lozmass)
    b=hizsfrm - m*(hizmass)
    #predictedsfrmz=m*masshiz+b
    #devz=sfrmhiz-predictedsfrmz
    predictedsfrmz=m*massall+b
    devz=sfrmall-predictedsfrmz

    psplotinit("sfrz.ps")
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(0.,0.9,0.,300.,0)
    ppgplot.pglab("z",sumsfrlabel,"")
    pgpnts(zhiz,sumsfrhiz,symbols)
    errory(zhiz,sumsfrhiz,errsumsfrhiz)

    xlabel = .1
    ylabel = 280
    ystep = 15
    dy=3
    dx=.05

    ppgplot.pgsls(1)
    ppgplot.pgsch(1.)
    ppgplot.pgslw(deflw)
    for i in range(len(ref)):
        putlabelpt(xlabel,ylabel,dx,dy,ref[i],refsymb[i])
        ylabel = ylabel - ystep

    ppgplot.pgend()
    os.system("cp sfrz.ps /Users/rfinn/clusters/papers/paper2/.")
    os.system("cp sfrz.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

    psplotinit("sfrmcl.ps")
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(-.6,1.8,0.,280.,0,10)
    ppgplot.pglab(mcllabel,sumsfrlabel,"")
    #ppgplot.pgpt(N.log10(masshiz),sumsfrhiz,3)
    pgpnts(N.log10(masshiz),sumsfrhiz,symbols)
    errory(N.log10(masshiz),sumsfrhiz,errsumsfrhiz)
    x=N.log10(N.compress(kodamamasshiz > 0,kodamamasshiz))
    y=(N.compress(kodamamasshiz > 0,sumsfrhiz))
    x1=N.log10(N.compress(kodamamasshiz > 0,masshiz))
    y1=(N.compress(kodamamasshiz > 0,sumsfrhiz))
    ppgplot.pgpt(x,y,3)
    ppgplot.pgsls(1)
    ppgplot.pgslw(1)
    for i in range(len(x1)):
        #x = N.array([N.log10(masshiz[i]),N.log10(kodamamasshiz[i])],'f')
        #y = N.array([sumsfrhiz[i],sumsfrhiz[i]],'f')
        x2 = N.array([x1[i],x[i]],'f')
        y2 = N.array([y1[i],y[i]],'f')
        ppgplot.pgline(x2,y2)

    for i in range(len(masshiz)):
        if ewcor[i] > 1.1:
            y1=sumsfrhiz[i]
            y2=sumsfrhiz[i]*ewcor[i]
            x1 = N.log10(masshiz[i])
            x = N.array([N.log10(masshiz[i]),N.log10(masshiz[i])],'f')
            y = N.array([y1,y2],'f')
            #ppgplot.pgline(x,y)
            y2=N.array([y2],'f')
            x1=N.array([x1],'f')
            #ppgplot.pgpt(x1,y2,3)
    xlabel = -.4
    ylabel = 250
    ystep = 12
    dy=4
    dx=.05
    ppgplot.pgsch(1.)
    ppgplot.pgslw(deflw)
    for i in range(len(ref)):
        putlabelpt(xlabel,ylabel,dx,dy,ref[i],refsymb[i])
        ylabel = ylabel - ystep

    ppgplot.pgend()
    os.system("cp sfrmcl.ps /Users/rfinn/clusters/papers/paper2/.")
    os.system("cp sfrmcl.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

    psplotinit("sfrmclz.ps")
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(0,1,-.8,2.2,0,20)
    ppgplot.pglab("z",sfrmlabel,"")
    y1 =  N.log10(sumsfrhiz/masshiz)
    pgpnts(zhiz,y1,symbols)
    errorylog(zhiz,sfrmhiz,errsfrmhiz)
    #x=(N.compress(kodamamasshiz > 0,kodamamasshiz))
    #y=(N.compress(kodamamasshiz > 0,sumsfrhiz))
    #z=(N.compress(kodamamasshiz > 0,zhiz))
    #y3 = (N.compress(kodamamasshiz > 0,y1))
    #y2 =  N.log10(y/x)
    #ppgplot.pgpt(z,y2,3)
    #x=zhiz.tolist()+z.tolist()
    #y=y1.tolist() + y2.tolist()
    ppgplot.pgsls(2)
    ppgplot.pgslw(1)
    #drawbinned(x,y,2)

    y = N.log10(N.array([hizsfrm,lozsfrm],'f'))
    x = N.array([hiz,loz],'f')
    ppgplot.pgline(x,y)
    ppgplot.pgsls(1)
    ppgplot.pgslw(1)

    #for i in range(len(z)):
    #    #y = N.array([y3[i],y2[i]],'f')
    #    #x = N.array([z[i],z[i]],'f')
    #    y = N.array([y3[i],y2[i]],'f')
    #    x = N.array([z[i],z[i]],'f')
    #    ppgplot.pgline(x,y)

    xlabel = .55
    ylabel = 0.5
    ystep = 0.15
    dy=.05
    dx=.02

    ppgplot.pgsls(1)
    ppgplot.pgsch(1.)
    ppgplot.pgslw(deflw)
    for i in range(len(ref)-1):
        putlabelpt(xlabel,ylabel,dx,dy,ref[i],refsymb[i])
        ylabel = ylabel - ystep
    ppgplot.pgend()
    os.system("cp sfrmclz.ps /Users/rfinn/clusters/papers/paper2/.")
    os.system("cp sfrmclz.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

    psplotinit("residualz.ps")
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(0,1,-100,100.,0,00)
    ppgplot.pglab("z","\gD"+sfrmlabel,"")
    #ppgplot.pgpt(N.log10(masshiz),dev,3)
    pgpnts(zall,devz,symbolall)
    #errory(N.log10(masshiz),sumsfrhiz,errsumsfrhiz)
    ppgplot.pgsls(1)
    ppgplot.pgslw(1)
    #for i in range(len(x1)):
    #    #x = N.array([N.log10(masshiz[i]),N.log10(kodamamasshiz[i])],'f')
    #    #y = N.array([sumsfrhiz[i],sumsfrhiz[i]],'f')
    #    x2 = N.array([x1[i],x[i]],'f')
    #    y2 = N.array([y1[i],y[i]],'f')
    #    ppgplot.pgline(x2,y2)
    x=N.arange(-5,5,.1)
    y=N.zeros(len(x),'f')
    ppgplot.pgsls(2)
    ppgplot.pgline(x,y)
    ppgplot.pgsls(1)

    xlabel = .5
    ylabel = -40
    ystep = 8
    dy=4
    dx=.05
    ppgplot.pgsch(1.)
    ppgplot.pgslw(deflw)
    for i in range(len(ref)-1):
        putlabelpt(xlabel,ylabel,dx,dy,ref[i],refsymb[i])
        ylabel = ylabel - ystep

    ppgplot.pgend()
    os.system("cp residualz.ps /Users/rfinn/clusters/papers/paper2/.")
    os.system("cp residualz.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

    psplotinit("residualmcl.ps")
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(-.6,1.8,-70,70.,0,10)
    ppgplot.pglab(mcllabel,"\gD"+sfrmlabel,"")
    #ppgplot.pgpt(N.log10(masshiz),dev,3)
    pgpnts(N.log10(massall),dev,symbolall)
    ppgplot.pgsls(2)
    ppgplot.pgslw(1)
    drawbinned(N.log10(massall),dev,2)
    print "deviation vs. Mcl"
    dospear(massall,dev)
    ppgplot.pgsls(1)
    ppgplot.pgslw(1)

    x=N.arange(-5,5,.1)
    y=N.zeros(len(x),'f')
    ppgplot.pgsls(4)
    ppgplot.pgline(x,y)
    ppgplot.pgsls(1)
    xlabel = .7
    ylabel = 55
    ystep = 5
    dy=2
    dx=.05
    ppgplot.pgsch(1.)
    ppgplot.pgslw(deflw)
    for i in range(len(ref)-1):
        putlabelpt(xlabel,ylabel,dx,dy,ref[i],refsymb[i])
        ylabel = ylabel - ystep

    ppgplot.pgend()
    os.system("cp residualmcl.ps /Users/rfinn/clusters/papers/paper2/.")
    os.system("cp residualmcl.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

    psplotinit("sfrmclmcl.ps")
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(-.5,1.8,-.8,2.2,0,30)
    ppgplot.pglab(mcllabel,sfrmlabel,"")
    y1 =  N.log10(sumsfrhiz/masshiz)
    x1 = N.log10(masshiz)
    pgpnts(x1,y1,symbols)
    errorylog(x1,sfrmhiz,errsfrmhiz)
    km=(N.compress(kodamamasshiz > 0,kodamamasshiz))
    ksfr=(N.compress(kodamamasshiz > 0,sumsfrhiz))
    kz=(N.compress(kodamamasshiz > 0,zhiz))
    #x1=N.log10(N.compress(kodamamasshiz > 0,masshiz))
    #y2=N.log10(N.compress(kodamamasshiz > 0,sumsfrhiz))
    mysfrm = (N.compress(kodamamasshiz > 0,y1))
    mym = (N.compress(kodamamasshiz > 0,x1))
    y2 =  N.log10(ksfr/km)
    x2 = N.log10(km)
    #ppgplot.pgpt(x2,y2,3)

    ppgplot.pgsls(1)
    ppgplot.pgslw(1)
    for i in range(len(km)):
        y = N.array([mysfrm[i],y2[i]],'f')
        x = N.array([mym[i],x2[i]],'f')
        #ppgplot.pgline(x,y)

    xlabel = -.3
    ylabel = .5
    ystep = 0.15
    dy=.05
    dx=.05

    ppgplot.pgsls(1)
    ppgplot.pgsch(1.)
    ppgplot.pgslw(deflw)
    for i in range(len(ref)-1):
        putlabelpt(xlabel,ylabel,dx,dy,ref[i],refsymb[i])
        ylabel = ylabel - ystep
    ppgplot.pgend()
    os.system("cp sfrmclmcl.ps /Users/rfinn/clusters/papers/paper2/.")
    os.system("cp sfrmclmcl.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

    output=open("sfrmassapj.tex",'w')

    output.write("\\begin{deluxetable}{lccccccccccc} \n")
    output.write("\\tablecaption{Integrated \ha \ SFRs of Galaxy Clusters \\label{cljlitsfrs}}\n")
    output.write("\\tablehead{\colhead{Name}  	&\colhead{z}  	& \colhead{$\sigma$} & \multicolumn{2}{c}{$R_{200}$} & \multicolumn{2}{c}{Survey Radius} & \colhead{SFR Cor} & \colhead{$\Sigma$~SFR}& \colhead{$\Sigma$~SFR/M$_{cl}$} & \colhead{Tech} & \colhead{Ref.} \\\\ \colhead{(1)} & (2) & (3) & (4) & (5) & (6) & (7) & (8) & (9) & (10) & (11) &(12)} \n")
    output.write("\startdata \n")
    for i in range(len(zhiz)):
        output.write("%s & %5.3f & %4.0f $\pm$ %3.0f & %5.2f & %5.2f & %5.2f & %5.2f & %5.2f & %5.1f $\pm$ %5.2f & %5.2f $\pm$ %5.2f & %s & %s \\\\ \n" % (clname[i],zhiz[i],sigmahiz[i],errsigmahiz[i],r200arcmin[i],r200Mpc[i],surveyrarcmin[i],surveyr200[i],totalcor[i],sumsfrhiz[i],errsumsfrhiz[i],sfrmhiz[i],errsfrmhiz[i],technique[i],refnumber[i]))
    output.write("\\tablerefs{(1) This work; (2) \citeauthor{finn04a} \citeyear{finn04a}; (3) \citeauthor{balogh00} \citeyear{balogh00}; (4) \citeauthor{couch01} \citeyear{couch01}; (5) \citeauthor{balogh02} \citeyear{balogh02}; (6) \citeauthor{kodama04} \\citeyear{kodama04}.} \n")
    output.write("\enddata")
    output.write("\\tablecomments{Columns: (1) Cluster name. (2) Redshift. (3) Velocity dispersion in $\\rm km\ s^{-1}$. Velocity dispersions for AC~114 and Abell~1689 are calculated from $L_X$ using best-fit $L_X - \\sigma$ relation of \citet{mahdavi01} because measured dispersions are inflated by substructure. (4) $R_{200}$ in arcmin. (5) $R_{200}$ in \h \ Mpc. (6) Survey radius in arcmin. (7) Survey radius in units of $R_{200}$. (8) Total correction applied to integrated SFR to account for incomplete sampling within $0.5 \\times R_{200}$, aperture corrections, and different EW limits.  (9) Integrated SFR in \smy. (10) Integrated SFR per cluster mass, in units of $h_{100}^{-3}\ \\rm M_\odot \ yr^{-1}\ / \ 10^{14} \  M_\odot$. (11) Observing Technique: I = narrowband imaging, S = spectroscopic survey. (12) References.} \n")
#
    output.write("\end{deluxetable} \n")
    output.close()
    os.system("cp sfrmassapj.tex /Users/rfinn/clusters/papers/paper2/.")
    os.system("cp sfrmassapj.tex /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

def plotspecewz(z,ew,errew,ncl):
    ppgplot.pgslw(6)
    #psplotinit("sfrmcl.ps")
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(.52,1.,-18.,90.,0,0)
    ppgplot.pglab("z",ewlabel,"")
    ppgplot.pgpt(z,ew,7)
    #error_y z ew errew
    
    
    x = N.arange(-1,11,2)
    y1 = N.ones(len(x),'f')
    y = 0*y1
    ppgplot.pgline(x,y)
    
    y = ewmin*y1
    ppgplot.pgsls(4)
    ppgplot.pgline(x,y)
    y = -1.*y
    ppgplot.pgline(x,y)
    
    
    #ctype red
    ppgplot.pgslw(1)
    ppgplot.pgsls(1)
    y = N.arange(-501,501,2)
    x = 0.98*c.z[ncl]*N.ones(len(y),'f')
    ppgplot.pgline(x,y)

    y = N.arange(-501,501,2)
    x = 1.02*c.z[ncl]*N.ones(len(y),'f')
    ppgplot.pgline(x,y)

    ppgplot.pgslw(4)
    if (ncl == 0):
        label = "(a) CL1040"
    if (ncl == 1):
        label = "(b) CL1054-12"
    if (ncl == 2):
        label = "(c) CL1216"
    ppgplot.pgtext(.81,77,label)

def plotspecewzall():
    psplotinit("haewzall.ps")
    ppgplot.pgslw(6)
    #psplotinit("sfrmcl.ps")
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(.52,.93,-18.,90.,0,0)
    ppgplot.pglab("z",ewlabel,"")

    z = N.compress((g0.spec > 0) & (g0.final > 0),g0.specz)
    ew = N.compress((g0.spec > 0) & (g0.final > 0),g0.ew) 
    erry = N.compress((g0.spec > 0) & (g0.final > 0),g0.errew)
    ppgplot.pgpt(z,ew,3)

    z = N.compress((g1.spec > 0) & (g1.final > 0),g1.specz)
    ew = N.compress((g1.spec > 0) & (g1.final > 0),g1.ew) 
    errew = N.compress((g1.spec > 0) & (g1.final > 0),g1.errew)
    ppgplot.pgpt(z,ew,17)
    
    z = N.compress((g2.spec > 0) & (g2.final > 0),g2.specz)
    ew = N.compress((g2.spec > 0) & (g2.final > 0),g2.ew) 
    errew = N.compress((g2.spec > 0) & (g2.final > 0),g2.errew)
    ppgplot.pgpt(z,ew,7)

    #error_y z ew errew
    

    
    x = N.arange(-1,11,2)
    y1 = N.ones(len(x),'f')
    y = 0*y1
    ppgplot.pgline(x,y)
    
    y = ewmin*y1
    ppgplot.pgsls(4)
    ppgplot.pgline(x,y)
    y = -1.*y
    ppgplot.pgline(x,y)
        
    
    #ctype red
    ppgplot.pgslw(1)
    ppgplot.pgsls(1)
    y = N.arange(-501,501,2)
    x = 0.98*c.z[0]*N.ones(len(y),'f')
    ppgplot.pgline(x,y)

    y = N.arange(-501,501,2)
    x = 1.02*c.z[0]*N.ones(len(y),'f')
    ppgplot.pgline(x,y)

    x = 0.98*c.z[1]*N.ones(len(y),'f')
    ppgplot.pgline(x,y)

    y = N.arange(-501,501,2)
    x = 1.02*c.z[1]*N.ones(len(y),'f')
    ppgplot.pgline(x,y)

    #x = 0.98*c.z[2]*N.ones(len(y),'f')
    #account for mismatch b/w cluster and filter
    x = 0.98*.802*N.ones(len(y),'f')
    ppgplot.pgline(x,y)

    y = N.arange(-501,501,2)
    #x = 1.02*c.z[2]*N.ones(len(y),'f')
    x = 1.02*.802*N.ones(len(y),'f')
    ppgplot.pgline(x,y)

    xlabel = .56
    ylabel = 80
    ystep = 7
    dy=2
    dxl=.01
    dxr=.01
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,17)

    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,7)

    ppgplot.pgend()

def plotspecewz(z,ew,errew,ncl):
    ppgplot.pgslw(6)
    #psplotinit("sfrmcl.ps")
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(.52,1.,-18.,90.,0,0)
    ppgplot.pglab("z",ewlabel,"")
    ppgplot.pgpt(z,ew,7)
    #error_y z ew errew
    
    
    x = N.arange(-1,11,2)
    y1 = N.ones(len(x),'f')
    y = 0*y1
    ppgplot.pgline(x,y)
    
    y = ewmin*y1
    ppgplot.pgsls(4)
    ppgplot.pgline(x,y)
    y = -1.*y
    ppgplot.pgline(x,y)
    
    
    #ctype red
    ppgplot.pgslw(1)
    ppgplot.pgsls(1)
    y = N.arange(-501,501,2)
    x = 0.98*c.z[ncl]*N.ones(len(y),'f')
    ppgplot.pgline(x,y)

    y = N.arange(-501,501,2)
    x = 1.02*c.z[ncl]*N.ones(len(y),'f')
    ppgplot.pgline(x,y)

    ppgplot.pgslw(4)
    if (ncl == 0):
        label = "(a) CL1040"
    if (ncl == 1):
        label = "(b) CL1054-12"
    if (ncl == 2):
        label = "(c) CL1216"
    ppgplot.pgtext(.81,77,label)

def plotspecewzall2():#don't use photozs
    psplotinit("haewzall2.ps")
    ppgplot.pgslw(6)
    #psplotinit("sfrmcl.ps")
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(0.52,.95,-18.,70.,0,0)
    ppgplot.pglab("z",ewlabel,"")

    x= g0.sn
    xmin = -20
    z = N.compress((g0.spec > 0) & (x > xmin),g0.specz)
    ew = N.compress((g0.spec > 0) & (x > xmin),g0.ew) 
    erry = N.compress((g0.spec > 0) & (x > xmin),g0.errew)
    ppgplot.pgpt(z,ew,3)

    x= g1.sn
    z = N.compress((g1.spec > 0) & (x > xmin),g1.specz)
    ew = N.compress((g1.spec > 0) & (x > xmin),g1.ew) 
    errew = N.compress((g1.spec > 0) & (x > xmin),g1.errew)
    ppgplot.pgpt(z,ew,17)

    x= g2.sn
    z = N.compress((g2.spec > 0) & (x > xmin),g2.specz)
    ew = N.compress((g2.spec > 0) & (x > xmin),g2.ew) 
    errew = N.compress((g2.spec > 0) & (x > xmin),g2.errew)
    ppgplot.pgpt(z,ew,7)

    #error_y z ew errew
    

    
    x = N.arange(-1,11,2)
    y1 = N.ones(len(x),'f')
    y = 0*y1
    ppgplot.pgline(x,y)
    
    y = ewmin*y1
    ppgplot.pgsls(4)
    ppgplot.pgline(x,y)
    y = -1.*y
    ppgplot.pgline(x,y)
        
    
    #ctype red
    ppgplot.pgslw(1)
    ppgplot.pgsls(1)
    y = N.arange(-501,501,2)
    x = 0.98*c.z[0]*N.ones(len(y),'f')
    ppgplot.pgline(x,y)

    y = N.arange(-501,501,2)
    x = 1.02*c.z[0]*N.ones(len(y),'f')
    ppgplot.pgline(x,y)

    x = 0.98*c.z[1]*N.ones(len(y),'f')
    ppgplot.pgline(x,y)

    y = N.arange(-501,501,2)
    x = 1.02*c.z[1]*N.ones(len(y),'f')
    ppgplot.pgline(x,y)

    x = 0.98*c.z[2]*N.ones(len(y),'f')
    ppgplot.pgline(x,y)

    y = N.arange(-501,501,2)
    x = 1.02*c.z[2]*N.ones(len(y),'f')
    ppgplot.pgline(x,y)

    xlabel = .56
    ylabel = 60
    ystep = 7
    dy=2
    dxl=.01
    dxr=.01
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,17)

    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,7)

    ppgplot.pgend()

def bootstrap(x):
    bootin = open("bootstrap.in",'w')
    for xin in x:
        s = str(xin)+"\n"
        bootin.write(s)
    bootin.close()
    os.system("bootstrap \n")
    for line in open("bootstrap.out"):
        (med,errmin,errmax) = line.split()
    return med,errmin,errmax
def bootfieldsub(z,x,ptype):
    (med,errmin,errmax) = bootstrap(x)
    med=float(med)
    errmin=float(errmin)
    errmax=float(errmax)
    y = N.array([med],'f')
    ppgplot.pgpt(z,y,ptype)
    #print z,y
    #drawerr
    sch = ppgplot.pgqch()
    lw = ppgplot.pgqlw()
    ppgplot.pgsch(.5)
    ppgplot.pgslw(1)
    ppgplot.pgsah(2,180.,1)
    ppgplot.pgarro(z,med,z,(med+errmax))
    ppgplot.pgarro(z,med,z,(med-errmin))

    ppgplot.pgsch(sch)
    ppgplot.pgslw(lw)
    

def bootfield():
    psplotinit("bootfield.ps")
    ppgplot.pgslw(6)
    #psplotinit("sfrmcl.ps")
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(.6,.9,0.,11.,0,0)
    ppgplot.pglab("z",sfrlabel,"")
    ncl = 0
    sfrmin=2.
    z = N.array([c.z[ncl]],'f')
    x = N.compress((g0.finalsf > 0) & (g0.sfr > sfrmin),g0.sfr)
    bootfieldsub(z,x,3)
    ncl = 1
    z = N.array([c.z[ncl]],'f')
    x = N.compress((g1.finalsf > 0) & (g1.sfr > sfrmin),g1.sfr)
    bootfieldsub(z,x,17)
    ncl = 2
    z = N.array([c.z[ncl]],'f')
    x = N.compress((g2.finalsf > 0) & (g2.sfr > sfrmin),g2.sfr)
    bootfieldsub(z,x,7)

    ncl = 3
    z = N.array([c.z[ncl]],'f')
    x = N.compress((g3.ew > ewmin) & (g3.sn > 3) & (g3.sfr > sfrmin),g3.sfr)
    bootfieldsub(z,x,0)

    #field
    z = []
    sfr = []
    errsfr = []
    for line in open("/Users/rfinn/clusters/final/tresse_sfr.dat",'r'):
        t = line.split()
        z.append(float(t[0]))
        sfr.append(float(t[1]))
        errsfr.append(float(t[2]))
    z = N.array(z,'f')
    sfr = N.array(sfr,'f')
    sfr = N.compress(sfr > sfrmin,sfr)
    errsfr = N.array(errsfr,'f')
    (med,errmin,errmax) = bootstrap(sfr)
    x = N.arange(.5,1.5,.1)
    y = N.ones(len(x),'f')
    y1 = float(med)*y
    print "field sfr = ",float(med)
    ppgplot.pgsls(1)
    ppgplot.pgslw(2)
    ppgplot.pgline(x,y1)
    ppgplot.pgsls(4)
    ppgplot.pgslw(2)
    y1 = (float(med)+float(errmax))*y
    ppgplot.pgline(x,y1)
    y1 = (float(med)-float(errmin))*y
    ppgplot.pgline(x,y1)

    xlabel = .64
    ylabel = 9.8
    ystep = .6
    dy=.2
    dxl=.010
    dxr=.01
    sch = ppgplot.pgqch()
    ppgplot.pgsch(1.3)
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,3)

    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,17)

    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,7)

    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CLJ0023+0423B")
    xlin = N.array([xlabel-dxr],'f')
    ylin = N.array([ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgpt(xlin,ylin,0)

    ppgplot.pgtext(.84,6,"Field")
    ppgplot.pgend()
    os.system("cp bootfield.ps /Users/rfinn/clusters/papers/paper2/.")
    os.system("cp bootfield.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

def plotcumulativesfrew():
    psplotinit("cumsfrew.ps")
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(0,260,-0.03,1.03,0)
    ppgplot.pglab("EW H\ga+[NII] (\(2078))","Cumulative SFR","")
    bins=N.arange(0,300,5,'f')


    cumsfrewsub(bins,g0.finalsf,g0.sfrc,g0.ewc,1,deflw)
    cumsfrewsub(bins,g1.finalsf,g1.sfrc,g1.ewc,3,deflw)
    cumsfrewsub(bins,g2.finalsf,g2.sfrc,g2.ewc,2,deflw)
    ppgplot.pgslw(deflw)
    ppgplot.pgsci(1)
    x=N.array([10.,10.])
    y=N.array([-50,500])
    ppgplot.pgsls(4)#dashed line
    ppgplot.pgline(x,y)

    x=N.array([50.,50.])
    y=N.array([-50,500])
    ppgplot.pgsls(4)#dashed line
    ppgplot.pgline(x,y)

    ppgplot.pgsls(1)#solid line

    xlabel = 180.
    ylabel = .26
    ystep = .08
    dxl = 40.
    dxr = 5
    dy = .02
    ppgplot.pgslw(4)
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    ppgplot.pgsls(1)
    #ppgplot.pgsch(1)
    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgline(xlin,ylin)

    ylabel = ylabel - ystep
    ppgplot.pgslw(4)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    ppgplot.pgsls(3)
    #ppgplot.pgsch(1)
    #ppgplot.pgslw(2)
    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgline(xlin,ylin)

    ylabel = ylabel - ystep
    ppgplot.pgslw(4)
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    ppgplot.pgsls(2)
    #ppgplot.pgslw(2)
    #ppgplot.pgsch(1)
    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgline(xlin,ylin)


    ppgplot.pgend()
    os.system("cp cumsfrew.ps /Users/rfinn/clusters/papers/paper2/.")
def cumsfrewsub(bins,finalsf,sfrc,ewc,sls,slw):
    nbins=len(bins)
    sfr=N.compress((finalsf > 0),sfrc)
    totsfr=N.sum(sfr)
    mabs=N.compress(finalsf > 0,ewc)
    #bins=N.arange(min(mabs),max(mabs),10,'f')
    sfrpart=N.zeros(len(bins),'f')
    for i in range(nbins):
        sfrpart[i]=N.sum(N.compress(mabs < bins[i],sfr))/totsfr
    bins=N.array(bins,'f')
    sfrpart=N.array(sfrpart,'f')
    ppgplot.pgsci(1)
    ppgplot.pgsls(sls)
    ppgplot.pgslw(slw)
    ppgplot.pgline(bins,sfrpart)

def plotcumulativesfrk():
    psplotinit("cumsfrk.ps")
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(15.5,24,-0.03,1.03,0)
    ppgplot.pglab("K\dS\u (r < 2'')","Cumulative SFR","")
    bins=N.arange(16,24,.5,'f')


    cumsfrewsub(bins,g0.finalsf,g0.sfrc,g0.vltk2,1,deflw)
    cumsfrewsub(bins,g1.finalsf,g1.sfrc,g1.vltk2,3,deflw)
    cumsfrewsub(bins,g2.finalsf,g2.sfrc,g2.vltk2,2,deflw)
    ppgplot.pgslw(deflw)
    ppgplot.pgsci(1)
    x=N.array([10.,10.])
    y=N.array([-50,500])
    ppgplot.pgsls(4)#dashed line
    ppgplot.pgline(x,y)
    ppgplot.pgsls(1)#solid line
    x=N.array([10.,30.])
    y=N.array([.5,.5])
    ppgplot.pgsls(4)#dashed line
    ppgplot.pgline(x,y)

    xlabel = 21.2
    ylabel = .26
    ystep = .08
    dxl = 1.
    dxr = .1
    dy = .02
    ppgplot.pgslw(4)
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    ppgplot.pgsls(1)
    #ppgplot.pgsch(1)
    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgline(xlin,ylin)

    ylabel = ylabel - ystep
    ppgplot.pgslw(4)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    ppgplot.pgsls(3)
    #ppgplot.pgsch(1)
    #ppgplot.pgslw(2)
    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgline(xlin,ylin)

    ylabel = ylabel - ystep
    ppgplot.pgslw(4)
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    ppgplot.pgsls(2)
    #ppgplot.pgslw(2)
    #ppgplot.pgsch(1)
    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgline(xlin,ylin)

    ppgplot.pgend()
    os.system("cp cumsfrk.ps /Users/rfinn/clusters/papers/paper2/.")

def plotcumulativesfrMR():

    minmag = -24.
    maxmag = -17.
    magstep = .5

    psplotinit("mmtediMR.ps")
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(minmag,maxmag,-0.03,30,0)
    ppgplot.pglab("M\dR\u","N\dgal\u","")
    bins=N.arange(minmag,maxmag,magstep,'f')

    magmmt = N.compress((g0.memb > 0),g0.MR)
    sfrmmt = N.compress((g0.memb > 0),g0.sfrc)
    signif = N.compress((g0.memb > 0),g0.signifc)
    magntt = N.compress((ntt0.memb > 0),ntt0.MR)
    y0compl=cumsfrMRsub(bins,magmmt,sfrmmt,signif,magntt)

    magmmt = N.compress((g1.memb > 0),g1.MR)
    sfrmmt = N.compress((g1.memb > 0),g1.sfrc)
    signif = N.compress((g1.memb > 0),g1.signifc)
    magntt = N.compress((ntt1.memb > 0),ntt1.MR)
    y1compl=cumsfrMRsub(bins,magmmt,sfrmmt,signif,magntt)

    magmmt = N.compress((g2.memb > 0),g2.MR)
    sfrmmt = N.compress((g2.memb > 0),g2.sfrc)
    signif = N.compress((g2.memb > 0),g2.signifc)
    magntt = N.compress((ntt2.memb > 0),ntt2.MR)
    y2compl=cumsfrMRsub(bins,magmmt,sfrmmt,signif,magntt)
    
    psplotinit("cumsfrMR.ps")
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(minmag,maxmag,-0.03,1.03,0)
    ppgplot.pglab("M\dR\u","Cumulative SFR","")
    #bins=N.arange(minmag,maxmag,.5,'f')

    sfr=N.compress((g0.finalsf > 0),g0.sfrc)
    mabs=N.compress((g0.finalsf > 0),g0.Mr)
    ew=N.compress((g0.finalsf > 0),g0.ewc)
    output=open('cl1040mabssfr.dat','w')
    for i in range(len(sfr)):
	string = "%6.2f %6.2f %6.2f \n" %(mabs[i],sfr[i],ew[i])
	output.write(string)
    output.close()
    sfr=N.compress((g1.finalsf > 0),g1.sfrc)
    mabs=N.compress((g1.finalsf > 0),g1.Mr)
    ew=N.compress((g1.finalsf > 0),g1.ewc)
    output=open('cl1054mabssfr.dat','w')
    for i in range(len(sfr)):
	string = "%6.2f %6.2f %6.2f \n" %(mabs[i],sfr[i],ew[i])
	output.write(string)
    output.close()
    sfr=N.compress((g2.finalsf > 0),g2.sfrc)
    mabs=N.compress((g2.finalsf > 0),g2.Mr)
    ew=N.compress((g2.finalsf > 0),g2.ewc)
    output=open('cl1216mabssfr.dat','w')
    for i in range(len(sfr)):
	string = "%6.2f %6.2f %6.2f \n" %(mabs[i],sfr[i],ew[i])
	output.write(string)
    output.close()
    cumsfrMRsub2(bins,g0.finalsf,g0.sfrc,g0.MR,1,deflw,y0compl)
    cumsfrMRsub2(bins,g1.finalsf,g1.sfrc,g1.MR,3,deflw,y1compl)
    cumsfrMRsub2(bins,g2.finalsf,g2.sfrc,g2.MR,2,deflw,y2compl)
    ppgplot.pgslw(deflw)
    ppgplot.pgsci(1)
    x=N.array([10.,10.])
    y=N.array([-50,500])
    ppgplot.pgsls(4)#dashed line
    ppgplot.pgline(x,y)
    ppgplot.pgsls(1)#solid line
    x=N.array([-25.,-15.])
    y=N.array([.5,.5])
    ppgplot.pgsls(4)#dashed line
    ppgplot.pgline(x,y)

    xlabel = -22.5
    ylabel = .92
    ystep = .08
    dxl = 1.
    dxr = .1
    dy = .02
    ppgplot.pgslw(4)
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    ppgplot.pgsls(1)
    #ppgplot.pgsch(1)
    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgline(xlin,ylin)

    ylabel = ylabel - ystep
    ppgplot.pgslw(4)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    ppgplot.pgsls(3)
    #ppgplot.pgsch(1)
    #ppgplot.pgslw(2)
    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgline(xlin,ylin)

    ylabel = ylabel - ystep
    ppgplot.pgslw(4)
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    ppgplot.pgsls(2)
    #ppgplot.pgslw(2)
    #ppgplot.pgsch(1)
    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgline(xlin,ylin)

    ppgplot.pgend()
    os.system("cp cumsfrMR.ps /Users/rfinn/clusters/papers/paper2/.")
    os.system("cp cumsfrMR.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

def cumsfrMRsub(bins,magmmt,sfrmmt,signif,magntt):
    x = magmmt
    ybin=N.zeros(len(bins),'f')
    for i in range(len(bins)-1):
        ybin[i]= len(N.compress((x >= bins[i]) & (x < bins[i+1]),x))
    drawhist(bins,ybin)
    ybinmmt=ybin
    x = magntt
    ybin=N.zeros(len(bins),'f')
    for i in range(len(bins)-1):
        ybin[i]= len(N.compress((x >= bins[i]) & (x < bins[i+1]),x))
    ppgplot.pgsls(2)
    drawhist(bins,ybin)
    ybinntt=ybin
    ycompl = N.zeros(len(ybinntt),'f')
    for i in range(len(ybinntt)):
        if (ybinmmt[i] > 0):
            if (ybinntt[i] > 0):
                ycompl[i]=ybinmmt[i]/ybinntt[i]
                if ycompl[i] > 1:
                    ycompl[i] = 1.
    return ycompl

def cumsfrMRsub2(bins,finalsf,sfrc,mag,sls,slw,compl):
    nbins=len(bins)
    #print len(bins),len(finalsf),len(sfrc),len(mag)
    sfr=N.compress((finalsf > 0),sfrc)
    mabs=N.compress((finalsf > 0),mag)
    #print len(sfr),len(mabs)
    sfrpart=N.zeros((len(bins)-1),'f')
    x=N.zeros((len(bins)-1),'f')
    temp=N.zeros((len(bins)-1),'f')
    for i in range(len(bins)-1):
        sfrt=0.
        temp[i]=0.
        for k in range(len(mabs)):
            if (mabs[k] > bins[i]) & (mabs[k] <= bins[i+1]):
                temp[i]=temp[i]+sfr[k]  
        if compl[i] > 0:
            temp[i]=temp[i]/compl[i]
        x[i]=(bins[i]+bins[i+1])/2.

    for i in range(len(temp)):
        for j in range((i+1)):
            sfrpart[i]=N.sum(temp[:j])
        #print i,len(bins),bins[i],bins[i+1],compl[i],temp[i],sfrpart[i]
    sfrtot=sfrpart[len(sfrpart)-1]
    print "sfrtot = ",sfrtot
    sfrpart=sfrpart/sfrtot
    bins=N.array(bins,'f')
    #sfrpart=N.array(sfrpart,'f')
    #for i in range(len(sfrpart)):
    #    print x[i], sfrpart[i]
    ppgplot.pgsci(1)
    ppgplot.pgsls(sls)
    ppgplot.pgslw(slw)
    ppgplot.pgline(x,sfrpart)

def plotfitzp(xj,yj,fj,errfj,ncl):
    #if ncl == 5:
    infile = open("MS1054/2mass.xy.coords", 'r')
    zp=22.09
    minmag=13.
    maxmag=19.

    if ncl == 6:
        infile = open("~/field/final/HDFN1/2mass.xy.coords", 'r')
        zp=21.85
        minmag=13.
        maxmag=17.
    x2m=[]
    y2m=[]
    mj2mass=[]
    errmj2mass=[]
    for line in infile:
        if line.find('#') > -1: #skip lines with '#' in them
            continue
        t=line.split()
        x2m.append(float(t[0]))
        y2m.append(float(t[1]))
        mj2mass.append(float(t[2]))
        errmj2mass.append(float(t[3]))
    infile.close()
    x2m=N.array(x2m,'f')
    y2m=N.array(y2m,'f')
    mj2mass=N.array(mj2mass,'f')
    errmj2mass=N.array(errmj2mass,'f')
    fjmatch=N.zeros(len(x2m),'f')
    errfjmatch=N.zeros(len(x2m),'f')
    ind=N.arange(0,len(xj),1,'i')
    #print "len(fjmatch) = ",len(fjmatch),len(x2m)
    for i in range(len(x2m)):
        d=N.sqrt((x2m[i]-xj)**2+(y2m[i]-yj)**2)
        sequence = zip(d,ind)
        (dmin,j)=min(sequence)
        #print i,j,dmin,len(fjmatch),len(x2m),len(fj)
        fjmatch[i]=fj[j]
        errfjmatch[i]=errfj[j]

    magj=zp-2.5*N.log10(fjmatch)
    errmagj=2.5*N.log10(1. + errfjmatch/fjmatch)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(minmag,maxmag,-1.5,1.5,0)
    #ppgplot.pglab("2.5log(f\dJ\u)","m\dJ\u(2MASS) - m_J","")
    ppgplot.pglab("m\dJ\u(2MASS)","m\dJ\u(2MASS) - m\dJ\u","")

    mdiff=mj2mass-magj
    #ppgplot.pgpt(magj,mj2mass)
    #errory(magj,mj2mass,errmj2mass)
    #errorx(magj,mj2mass,errmagj)

    tzp=N.arange(21.9,22.2,.001)
    zpbest=0.
    minvar=10000.
    for zp in tzp:
        magj=zp-2.5*N.log10(fjmatch)
        errmagj=2.5*N.log10(1. + errfjmatch/fjmatch)
        mdiff=abs(mj2mass-magj)
        var=0.
        n=0
        for i in range(len(mdiff)):
            #print i,mdiff[i],mj2mass[i],var,n
            if (mdiff[i] < 1.):
                if (mj2mass[i] < 17.):
                    var += mdiff[i] 
                    n+=1
                #var=(var/(1.*n))
        #print "zp = ", zp, "var = ",var, "n = ",n
        try:
            var=var/(1.*n)
        except:
            continue
        if var < minvar:
            zpbest=zp
            minvar=var
            #print "zp = ",zp, "var = ",var, "min var = ",minvar,zpbest
    x=N.arange((minmag-2.),(maxmag+2.),1)
    y=0.*x
    magj=zpbest-2.5*N.log10(fjmatch)
    errmagj=2.5*N.log10(1. + errfjmatch/fjmatch)
    mdiff=mj2mass-magj
    ppgplot.pgpt(mj2mass,mdiff)
    errory(mj2mass,mdiff,errmagj)
    errorx(mj2mass,mdiff,errmj2mass)

    ppgplot.pgsls(1)
    ppgplot.pgline(x,y)
    svar="ZP =%6.3f , std=%4.2f" %(zpbest,minvar)
    ppgplot.pgtext(13.5,1.2,svar)
def dostats():
    #membanalysis("CL1040-1155",g0.sn,g0.ew,g0.memb,g0.sfr)
    print "CL1040: memb=1 = ",N.sum(g0.nmemb)," specmemb = ",N.sum(g0.nspecmemb)," photz+emission = ",N.sum(g0.nphotoznemission), "nstarforming in final sample = ",N.sum(g0.finalsf)
    nonmembers("CL1040-1155",g0.signif,g0.signifc,g0.memb,g0.sfr,g0.sfrc,g0.ewc,g0.zroser,g0.zroserlo,g0.zroserhi,g0.zgreg,g0.zgreglo,g0.zgreghi,g0.x,g0.y,g0.starflag,g0.pclust,0,g0.specmemb,g0.final)
    #membanalysis("CL1054-12-1245",g1.sn,g1.ew,g1.memb,g1.sfr)
    print "CL1054: memb=1 = ",N.sum(g1.nmemb)," specmemb = ",N.sum(g1.nspecmemb)," photz+emission = ",N.sum(g1.nphotoznemission), "nstarforming in final sample = ",N.sum(g1.finalsf)
    nonmembers("CL1054-12-1245",g1.signif,g1.signifc,g1.memb,g1.sfr,g1.sfrc,g1.ewc,g1.zroser,g1.zroserlo,g1.zroserhi,g1.zgreg,g1.zgreglo,g1.zgreghi,g1.x,g1.y,g1.starflag,g1.pclust,1,g1.specmemb,g1.final)
    #membanalysis("CL1216-1201",g2.sn,g2.ew,g2.memb,g2.sfr)
    print "CL1216: memb=1 = ",N.sum(g2.nmemb)," specmemb = ",N.sum(g2.nspecmemb)," photz+emission = ",N.sum(g2.nphotoznemission), "nstarforming in final sample = ",N.sum(g2.finalsf)
    nonmembers("CL1216-1201",g2.signif,g2.signifc,g2.memb,g2.sfr,g2.sfrc,g2.ewc,g2.zroser,g2.zroserlo,g2.zroserhi,g2.zgreg,g2.zgreglo,g2.zgreghi,g2.x,g2.y,g2.starflag,g2.pclust,2,g2.specmemb,g2.final)

    #total SFR
    print "CL1040"
    print "TOTAL SFR = ", N.sum(N.compress(g0.finalsf > 0, g0.sfrc)), " or ",g0.totalsfr,"+/-",g0.totalsfrerr
    print "Number of SF galaxies = ", len(N.compress(g0.finalsf > 0, g0.sfrc)),len(N.compress(g0.final > 0, g0.sfr)),"final = 1",len(N.compress(g0.memb > 0,g0.sfr)), "(memb = 1)"
    print "Fraction of SF galaxies = ", g0.sffrac
    print "Number with memb=0 specmemb=1 = ", len(N.compress((g0.memb < 1) & (g0.specmemb > 0), g0.sfr))
    print "CL1054-12"
    print "TOTAL SFR = ", N.sum(N.compress(g1.finalsf > 0, g1.sfrc)), " or ",g1.totalsfr,"+/-",g1.totalsfrerr
    print "Number of SF galaxies = ", len(N.compress(g1.finalsf > 0, g1.sfr)),len(N.compress(g1.final > 0, g1.sfr)),"(final = 1),",len(N.compress(g1.memb > 0,g1.sfr)), "(memb = 1)"
    print "Fraction of SF galaxies = ", g1.sffrac
    print "Number with memb=0 specmemb=1 = ", len(N.compress((g1.memb < 1) & (g1.specmemb > 0), g1.sfr))
    print "CL1216"
    print "TOTAL SFR = ", N.sum(N.compress(g2.finalsf > 0, g2.sfrc)), " or ",g2.totalsfr,"+/-",g2.totalsfrerr
    print "Number of SF galaxies = ", len(N.compress(g2.finalsf > 0, g2.sfr)),len(N.compress(g2.final > 0, g2.sfr)),"(final = 1),",len(N.compress(g2.memb > 0,g2.sfr)), "(memb = 1)"
    print "Fraction of SF galaxies = ", g2.sffrac
    print "Number with memb=0 specmemb=1 = ", len(N.compress((g2.memb < 1) & (g2.specmemb > 0), g2.sfr))
    print "Number with memb=0 signifc=1 = ", len(N.compress((g2.memb < 1) & (g2.signif > 0), g2.sfrc))
    
    print "Median SFR for galaxies with SFR > 2:"
    print "CL1040 =  %5.2f  %5.2f (ave) +/- %5.2f" % (pylab.median(N.compress((g0.finalsf > 0) & (g0.sfr > 2),g0.sfr)), N.average(N.compress((g0.finalsf > 0) & (g0.sfr > 2),g0.sfr)),pylab.std(N.compress((g0.finalsf > 0) & (g0.sfr > 2),g0.sfr)))
    print "CL1054-12 =  %5.2f  %5.2f (ave) +/- %5.2f" % (pylab.median(N.compress((g1.finalsf > 0) & (g1.sfr > 2),g1.sfr)), N.average(N.compress((g1.finalsf > 0) & (g1.sfr > 2),g1.sfr)),pylab.std(N.compress((g1.finalsf > 0) & (g1.sfr > 2),g1.sfr)))
    print "CL1216 =  %5.2f  %5.2f (ave) +/- %5.2f" % (pylab.median(N.compress((g2.finalsf > 0) & (g2.sfr > 2),g2.sfr)), N.average(N.compress((g2.finalsf > 0) & (g2.sfr > 2),g2.sfr)),pylab.std(N.compress((g2.finalsf > 0) & (g2.sfr > 2),g2.sfr)))

    print "Fraction of SFR detected by spect :"
    temp1 = N.compress((g0.specmemb > 0) & (g0.finalsf > 0),g0.sfrc)
    temp2 = N.compress((g0.finalsf > 0),g0.sfrc)
    a = N.sum(N.compress((g0.specmemb > 0) & (g0.finalsf > 0),g0.sfrc))
    b = g0.totalsfr
    print "CL1040 %5.1f %5.1f %5.2f (%5.1f spec /%5.1f finalsf, %5.2f)" %(a,b,a/b,len(temp1),len(temp2),float(len(temp1))/float(len(temp2)))

    a = N.sum(N.compress((g1.specmemb > 0) & (g1.finalsf > 0),g1.sfrc))
    b = g1.totalsfr
    temp1 = N.compress((g1.specmemb > 0) & (g1.finalsf > 0),g1.sfrc)
    temp2 = N.compress((g1.finalsf > 0),g1.sfrc)
    print "CL1054-12 %5.1f %5.1f %5.2f (%5.1f spec /%5.1f finalsf, %5.2f)" %(a,b,a/b,len(temp1),len(temp2),float(len(temp1))/float(len(temp2)))

    a = N.sum(N.compress((g2.specmemb > 0) & (g2.finalsf > 0),g2.sfrc))
    b = g2.totalsfr
    temp1 = N.compress((g2.specmemb > 0) & (g2.finalsf > 0),g2.sfrc)
    temp2 = N.compress((g2.finalsf > 0),g2.sfrc)
    print "CL1216 %5.1f %5.1f %5.2f (%5.1f spec /%5.1f finalsf, %5.2f)" %(a,b,a/b,len(temp1),len(temp2),float(len(temp1))/float(len(temp2)))
    
def nonmembers(clust,signif,signifc,memb,sfr,sfrc,ewc,zroser,zroserlo,zroserhi,zgreg,zgreglo,zgreghi,xpos,ypos,starflag,pclust,ncl,specmemb,final):
    zcl=c.z[ncl]
    print clust," NONMEMBERS:"
    print "Ngal final=1 = ",len(N.compress((final > 0),sfr)),N.sum(final)
    print "Ngal final=0, signif=1 starflag<1= ",len(N.compress((final < 1) & (signif > 0) & (starflag < 1),sfr)),N.sum(N.compress((final < 1) & (signif > 0) & (starflag < 1),sfr)),N.sum(N.compress((final < 1) & (signif > 0) & (starflag < 1),sfrc))
    print "Ngal final=0, signifc=1 starflag<1 = ",len(N.compress((final < 1)  & (starflag < 1) & (signifc > 0),sfr)),N.sum(N.compress((final < 1)  & (starflag < 1)& (signifc > 0),sfr)),N.sum(N.compress((final < 1)  & (starflag < 1) & (signifc > 0),sfrc))
    print "Ngal memb=1 = ",len(N.compress((memb > 0),sfr))

    print "Ngal memb=0, specmemb=1 = ",len(N.compress((memb < 1) & (specmemb > 0),sfr))
    print "Ngal memb=0, specmemb=1, starflag<1 = ",len(N.compress((final > 1) & (memb < 1) & (starflag < 1) & (specmemb > 0),sfr))
    print "Ngal memb=0, specmemb=1, signif=1 = ",len(N.compress((signif > 0) & (memb < 1) & (specmemb > 0),sfr))
    print "Ngal memb=0, specmemb=1, signifc=1 = ",len(N.compress( (signifc > 0) & (memb < 1) & (specmemb > 0),sfr))
    #print "Ngal memb=0, signif=1, signifc=1 = ",N.sum(N.compress((signif > 0) & (signifc > 0) & (memb < 1),sfr)),len(N.compress((signif > 0) & (signifc > 0) & (memb < 1),sfr))
    print "Ngal memb=0,signifc>0, specmemb<1,starflag<1,(zglo < zcl < zghi or Pclust>0.2) = ",len(N.compress((signifc > 0 ) & (memb < 1) & (starflag < 1) & (specmemb < 1) & (((zgreglo < zcl) & (zgreghi > zcl)) | (pclust > 0.2)),sfrc))
    print "Ngal memb=0,signifc>0, starflag<1,specmemb<1  = ",len(N.compress((signifc > 0 ) & (memb < 1) & (starflag < 1) & (specmemb < 1),sfrc))
    print "Ngal final=0,signifc>0  = ",len(N.compress((final < 1) & (signif > 0),sfrc))
    print "Ngal final=0,signif>0  = ",len(N.compress((signif > 0 ) & (final < 1),sfrc))
    print "Ngal memb=0, signif=1, starflag <1 = ",len(N.compress((signif > 0)  & (starflag < 1) & (memb < 1),sfr))
    print "Sum sfr for final=1 EW>50, fraction ",N.sum(N.compress((final > 0) & (signifc > 0) & (ewc > 50),sfr)),N.sum(N.compress((final > 0) & (signifc > 0) & (ewc > 50),sfr))/N.sum(N.compress((final > 0) & (signifc > 0),sfr)),N.sum(N.compress((final > 0) & (signifc > 0) & (ewc < 50),sfr))/N.sum(N.compress((final > 0) & (signifc > 0),sfr))
    sf =   N.compress((signifc > 0) & (memb < 1),sfr)
    sfc =  N.compress((signifc > 0) & (memb < 1),sfrc)
    zr =   N.compress((signifc > 0) & (memb < 1),zroser)
    zrlo = N.compress((signifc > 0) & (memb < 1),zroserlo)
    zrhi = N.compress((signifc > 0) & (memb < 1),zroserhi)
    zg =   N.compress((signifc > 0) & (memb < 1),zgreg)
    zglo = N.compress((signifc > 0) & (memb < 1),zgreglo)
    zghi = N.compress((signifc > 0) & (memb < 1),zgreghi)
    x =    N.compress((signifc > 0) & (memb < 1),xpos)
    y =    N.compress((signifc > 0) & (memb < 1),ypos)
    star = N.compress((signifc > 0) & (memb < 1),starflag)
    pc =   N.compress((signifc > 0) & (memb < 1),pclust)


    #for i in range(len(x)):
    #    print "%6.2f %6.2f %6.1f %6.1f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f" % (sf[i],sfc[i],x[i],y[i],zrlo[i],zr[i],zrhi[i],zglo[i],zg[i],zghi[i],star[i],pc[i])
    #print "total number with memb=1 = ",len(N.compress((memb > 0),sfr))
    #print "total SFR signif=1 = ",N.sum(N.compress((signif > 0),sfr))
    #print "total SFR of memb=0 = ",N.sum(N.compress((signif > 0 ) & (memb < 1),sfr)),len((N.compress((signif > 0 ) & (memb < 1),sfr)))
    #print "total SFR of memb=1 = ", N.sum(N.compress((signif > 0) & (memb > 0),sfr)),len((N.compress((signif > 0) & (memb > 0),sfr)))
    #print "total SFR of memb=0,signif>0,signifc>0 = ",N.sum(N.compress((signif > 0 ) & (signifc > 0 ) & (memb < 1),sfr)),len((N.compress((signif > 0 ) & (memb < 1),sfr)))
    #print "total SFR of memb=0,signif>0,signifc>0, Pclust > 0.2 = ",N.sum(N.compress((signif > 0 ) & (signifc > 0 ) & (memb < 1) & (pclust > .2),sfrc)),len((N.compress((signif > 0 ) & (signifc > 0 ) & (memb < 1) & (pclust > .2),sfrc)))
    #print "total SFR of memb=0,signif>0,signifc>0, zglo < zcl < zghi = ",N.sum(N.compress((signif > 0 ) & (signifc > 0 ) & (memb < 1) & (zgreglo < zcl) & (zgreghi > zcl),sfrc)),len(N.compress((signif > 0 ) & (signifc > 0 ) & (memb < 1) & (zgreglo < zcl) & (zgreghi > zcl),sfrc))
    #print "total SFR of memb=0,signifc>0, (zglo < zcl < zghi or Pclust>0.2) = ",N.sum(N.compress((signifc > 0 ) & (memb < 1) & (((zgreglo < zcl) & (zgreghi > zcl)) | (pclust > 0.2)),sfrc)),len(N.compress((signifc > 0 ) & (memb < 1) & (((zgreglo < zcl) & (zgreghi > zcl)) | (pclust > 0.2)),sfrc))

def membanalysis(cluster,sn,ew,memb,sfr):
    print cluster,": ",len(N.compress((sn > snmin) & (ew > ewmin) & (memb < 1),sfr))

    print "total SFR signif = 1 = ",N.sum(N.compress((sn > snmin) & (ew > ewmin),sfr))
    print "total SFR of memb=0 = ",N.sum(N.compress((sn > snmin) & (ew > ewmin) & (memb < 1),sfr)),len(N.compress((sn > snmin) & (ew > ewmin) & (memb < 1),sfr))
    print "total SFR of memb=1 = ", N.sum(N.compress((sn > snmin) & (ew > ewmin) & (memb > 0),sfr)),len(N.compress((sn > snmin) & (ew > ewmin) & (memb > 0),sfr))
    x = N.compress((sn > snmin) & (ew > ewmin) & (memb < 1),sfr)
    y = N.compress((sn > snmin) & (ew > ewmin) & (memb < 1),ew)
#    for i in range(len(x)):
#        print x[i], y[i]


def writeapjtable(name,x,y,fn,errfn,fj,errfj,ratio,errratio,contsub,errcontsub,ew,errew,sfr,errsfr,final,finalsf,tablename,ncl,magj,errmagj,ediscsra):
    snc=abs(contsub/errcontsub)
    i=0
    j=0
    v=2
    index=N.arange(0,len((name)),1)
    ra=N.zeros(len((name)),'f')
    dec=N.zeros(len((name)),'f')
    #n=str(name)
    dra=(y-c.yc[ncl])*c.pscale[ncl]
    ddec=(x-c.xc[ncl])*c.pscale[ncl]
    sortindex=N.take(index,N.argsort(ediscsra))
    textable=tablename+"sfrtableapj.tex"
    output = open(textable, 'w')
    dattable=tablename+"sfrtableapj.dat"
    output2 = open(dattable, 'w')
    output.write("\\begin{deluxetable}{rrrrrrrcrrc} \n")
    caption = "\\tablecaption{\ha \ Data for "+c.fullclustername[ncl]+" Galaxies \label{"+tablename+"}} \n"
    output.write(caption)
#    output.write("\\tablehead{\colhead{Name} & \colhead{$\delta$RA} & \colhead{$\delta$Dec} & \colhead{x} & \colhead{y} & \colhead{Flux$_n$} & \colhead{Flux$_J$} & \colhead{Ratio} & \colhead{Cont. Sub} & \colhead{\ewr(\ha)} & \colhead{SFR} \\\\  \colhead{(1)} &  \colhead{(2)} &  \colhead{(3)} &  \colhead{(4)} &  \colhead{(5)} &  \colhead{(6)} &  \colhead{(7)} &  \colhead{(8)} &  \colhead{(9)} &  \colhead{(10)} &  \colhead{(11)}} \n")
    output.write("\\tablehead{\colhead{Name} & \colhead{$\delta$RA} & \colhead{$\delta$Dec} & \colhead{Flux$_n$} & \colhead{Flux$_J$} & \colhead{J} & \colhead{Cont. Sub} & \colhead{SNR} & \colhead{\ewr(\ha)} & \colhead{SFR} & \colhead{SF}\\\\  \colhead{(1)} &  \colhead{(2)} &  \colhead{(3)} &  \colhead{(4)} &  \colhead{(5)} &  \colhead{(6)} &  \colhead{(7)} &  \colhead{(8)} &  \colhead{(9)} & \colhead{(10)} & \colhead{(11)}} \n")
    output.write( "\startdata \n")
    
    output2.write("#see Finn et al. 2005, ApJ, Aug 20 issue for details on units (astro-ph/0504578) \n")
    output2.write("#Name                 dRA      dDec     Fluxn      err    FluxJ      err     magJ      err  ContSub      err      SNR EW(Ha)  err   SFR   err SF_flag\n")

    #for i in range(ndata):
    i = 0
    j = 0
    for j in range(len(x)):
        i=int(sortindex[j])
        if (final[i] > 0):
            totname = name[i]
            shortname = totname[6:]
            sf = 0
            if (finalsf[i] > 0):
                sf=1
#            output.write("%s & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f $\pm$ %8.1f & %8.1f $\pm$ %8.1f & %8.4f $\pm$ %8.4f & %8.1f $\pm$ %8.1f & %5.1f $\pm$ %5.1f & %5.1f $\pm$ %5.1f \\\\ \n" % (shortname,dra,ddec,x[i],y[i],fn[i],errfn[i],fj[i],errfj[i],ratio[i],errratio[i],contsub[i],errcontsub[i],ew[i],errew[i],sfr[i],errsfr[i]))
#            output.write("%s & %8.1f & %8.1f & %8.1f $\pm$ %8.1f & %8.1f $\pm$ %8.1f & %8.3f $\pm$ %8.3f & %8.1f $\pm$ %8.1f & %5.1f $\pm$ %5.1f & %5.1f $\pm$ %5.1f \\\\ \n" % (shortname,dra,ddec,fn[i],errfn[i],fj[i],errfj[i],c.defratio[ncl],errratio[i],contsub[i],errcontsub[i],ew[i],errew[i],sfr[i],errsfr[i]))

            output.write("%s & %8.1f & %8.1f & %8.1f $\pm$ %8.1f & %8.1f $\pm$ %8.1f & %8.2f $\pm$ %8.2f & %8.2f $\pm$ %8.2f & %8.1f & %5.1f $\pm$ %5.1f & %5.1f $\pm$ %5.1f & %i \\\\ \n" % (shortname,dra[i],ddec[i],fn[i],errfn[i],fj[i],errfj[i],magj[i],errmagj[i],contsub[i],errcontsub[i],snc[i],ew[i],errew[i],sfr[i],errsfr[i],sf))
            output2.write("%s  %8.1f  %8.1f  %8.1f %8.1f %8.1f %8.1f %8.2f %8.2f %8.2f %8.2f %8.1f %5.1f %5.1f %5.1f %5.1f %i \n" % (shortname,dra[i],ddec[i],fn[i],errfn[i],fj[i],errfj[i],magj[i],errmagj[i],contsub[i],errcontsub[i],snc[i],ew[i],errew[i],sfr[i],errsfr[i],sf))
            #RAF temporarily changed ra and dec to x and y
            #output.write("%s & %8.1f & %8.1f & %8.1f $\pm$ %8.1f & %8.1f $\pm$ %8.1f & %8.2f $\pm$ %8.2f & %8.1f $\pm$ %8.1f & %5.1f $\pm$ %5.1f & %5.1f $\pm$ %5.1f \\\\ \n" % (shortname,y[i],x[i],fn[i],errfn[i],fj[i],errfj[i],magj[i],errmagj[i],contsub[i],errcontsub[i],ew[i],errew[i],sfr[i],errsfr[i]))
            j = j +1 
#            if (j > 75):
#                output.write("\enddata \n")
#                output.write("\end{deluxetable} \n")
#                output.write("\clearpage \n")
#                output.write("\\begin{deluxetable}{rrrrrrrrr} \n")
#                caption = "Table \\ref{" + tablename +"} \ha \ Data for cluster Galaxies \\\n"
#                output.write(caption)
#                #    output.write("\\tablehead{\colhead{Name} & \colhead{$\delta$RA} & \colhead{$\delta$Dec} & \colhead{x} & \colhead{y} & \colhead{Flux$_n$} & \colhead{Flux$_J$} & \colhead{Ratio} & \colhead{Cont. Sub} & \colhead{\ewr(\ha)} & \colhead{SFR} \\\\  \colhead{(1)} &  \colhead{(2)} &  \colhead{(3)} &  \colhead{(4)} &  \colhead{(5)} &  \colhead{(6)} &  \colhead{(7)} &  \colhead{(8)} &  \colhead{(9)} &  \colhead{(10)} &  \colhead{(11)}} \n")
#                output.write("\\tablehead{\colhead{Name} & \colhead{$\delta$RA} & \colhead{$\delta$Dec} & \colhead{Flux$_n$} & \colhead{Flux$_J$} & \colhead{Ratio} & \colhead{Cont. Sub} & \colhead{\ewr(\ha)} & \colhead{SFR} \\\\  \colhead{(1)} &  \colhead{(2)} &  \colhead{(3)} &  \colhead{(4)} &  \colhead{(5)} &  \colhead{(6)} &  \colhead{(7)} &  \colhead{(8)} &  \colhead{(9)}}  \n")
#                output.write( "\startdata \n")
#                j = 0

            
    output.write("\enddata \n")
#    output.write("\\tablecomments{Columns:  (1) Name is EDCSNJ followed by number listed in column.  (2) RA offset from BCG in arcseconds.  (3) DEC offset from BCG in arcseconds.  (4) Image x-position in pixels.  (5) Image y-position in pixels.  (6) Narrow-band flux in ADU/s. (7)  $J$-band flux in ADU/s. (8) Narrow-to-J flux ratio.  (9) Continuum-subtracted flux in ADU/s.  (10) Narrow-band \ewr \ in \AA.  (11) SFR in units of $h_{100}^{-2}\ \rm M_\odot \ yr^{-1}$. }\n")
#    output.write("\\tablecomments{Columns:  (1) Name is EDCSNJ followed by number listed in column.  (2) RA offset from BCG in arcseconds.  (3) DEC offset from BCG in arcseconds.  (4) Narrow-band flux in ADU/s. (5)  $J$-band flux in ADU/s. (6) Ratio used to estimate NB continuum from $J$-band flux.  (7) Continuum-subtracted flux in ADU/s.  (8) Narrow-band \ewr \ in \AA.  (9) SFR in units of $h_{100}^{-2}\ \\rm M_\odot \ yr^{-1}$. }\n")
    output.write("\\tablecomments{Columns:  (1) Name is EDCSNJ followed by number listed in column.  (2) RA offset from BCG in arcseconds.  (3) DEC offset from BCG in arcseconds.  (4) Narrow-band flux in ADU/s. (5)  $J$-band flux in ADU/s. (6) J isophotal magnitude with SExtractor error.  (7) Continuum-subtracted flux in ADU/s.  (8) Signal-to-noise ratio of continuum-subtracted flux.  (9) Narrow-band \ewr \ in \AA.  (10) SFR in units of $h_{100}^{-2}\ \\rm M_\odot \ yr^{-1}$.  (11) Star-forming galaxies that meet minimum continuum-subtract flux and \ewr \ cuts are denoted with 1. }\n")
        
    output.write("\end{deluxetable} \n")
    output.close()
    output2.close()


def writesdsstable(x,y,ew,errew,sfr,errsfr,finalsf,tablename,ncl,MR,nttmatchflag):
    textable=tablename+"sdss.dat"
    output = open(textable, 'w')
    for i in range(len(x)):
        if (finalsf[i] > 0):
            dra=(y[i]-c.yc[ncl])*c.pscale[ncl]
            ddec=(x[i]-c.xc[ncl])*c.pscale[ncl]
            output.write("%8.1f %8.1f %6.2f %6.1f %6.1f %6.2f %6.2f %s %i\n" % (dra,ddec,MR[i],ew[i],errew[i],sfr[i],errsfr[i],nttmatchflag[i],i))
    output.close()
    

class Cluster:
    def __init__(self):
        #print "dude - a cluster!"
        self.id = [1040,1054,1216,0023,0152,1054]


        self.filt=[8,6,4,1,2,1201]#pisces filter positions
        #self.jnoisea= [.210,.165,.095,.100,.100,.105]#empirical fit to sky  noise
        #self.jnoiseb= [.030,.040,.035,.035,.035,.015]#empirical fit to sky  noise
        #self.nbnoisea=[.020,.020,.005,.020,.020,.01]#empirical fit to sky  noise
        #self.nbnoiseb=[.030,.050,.575,.040,.040,.06]#empirical fit to sky  noise
        #dimsum flattening
        self.jnoisea= [.185,.125,.100,.175,.155,.100,.075]#empirical fit to sky  noise
        self.jnoiseb= [.030,.055,.015,.055,.035,.015,.035]#empirical fit to sky  noise
        self.nbnoisea=[.015,.025,.010,.015,.020,.01,.005]#empirical fit to sky  noise
        self.nbnoiseb=[.050,.045,.055,.070,.075,.05,.075]#empirical fit to sky  noise
        self.fullclustername=["CL 1040-1155","CL 1054-1245","CL1216-1201","CL J0023+0423B","CL0152-1357", "MS1054-0321", "HDFN-1"]
        self.tablename=["cl1040sfrtable","cl1054sfrtable","cl1216sfrtable","clj0023sfrtable","cl0152sfrtable","ms1054sfrtable","HDFN1sfrtable"]
        self.errzp=N.array([0.05,0.06,.03,0.05,.05,.05,.05]) #percent error in ZP
        self.sfrconv=N.array([2.4,2.3,2.1,1.19,3.18,(46.29),(7.94*2.5)]) #1 ADU/sec = sfrconv Mo/yr
        self.sfrconv=self.sfrconv/(h100**2) #account for cosmology
        self.fluxzp=N.array([8.55,8.46,9.54,9.21,9.21,9.21,9.21])#flux ZP in 1d-17 erg/s/cm^2
        self.jzp=N.array([24.92-0.06,24.93-0.08,24.93-0.08,24.94-0.06,24.94-0.06,22.04,21.85])
        self.z=N.array([0.704,0.748,0.794,0.845,0.833,0.832,0.804])
        self.pscale=N.array([0.18,0.18,0.18,0.18,0.18,0.5,.5])
        #self.defratio=N.array([0.0723,0.0718,0.0721,0.0488,0.072]) #ratio of filter fluxes
        #self.errdefratio=N.array([0.31,0.31,.031,0.0,0.006]) #error associated w/ratio
        self.ratiofitinter = N.array([.065,.066,.073,.0488,.085,.0648,.094],'f')
        #min to keep peak just shy of 10A
        #self.ratiofitinter = N.array([.062,.06,.072,.0488,.072],'f')
        #min to keep peak just > of -10A
        #self.ratiofitinter = N.array([.072,.07,.082,.0488,.072],'f')
        self.ratiofitinter = N.array([.064,.065,.077,.0488,.085,.065,.094],'f')
        self.ratiofitslope = N.array([-0.005,-.004,-.004,0.,0.,0.,0.],'f')

        self.ratiofitinterhiz = N.array([.0975,.116,.1075,.0488,.085,.065,.094],'f')
        self.ratiofitslopehiz = N.array([-0.03,-.036,-.02,0.,0.,0.,0.],'f')
        self.defratio=self.ratiofitinter + self.ratiofitslope*self.z
        self.errdefratio=N.array([0.005,0.005,.005,0.0,0.006,.006,.006]) #error associated w/ratio
        self.xc=N.array([394.,524.,389.,555.,428.,454.,469.])
        self.yc=N.array([348.,519.,341.,417.,422.,417.,450.])
        self.sig=N.array([418.,504.,1018.,415.,1000.,1170.,1000.])
        self.errsigp=N.array([41.,113.,46.,0.1*415.,0.1*self.sig[4],150.,150.],'f')
        self.errsigm=N.array([41.,65.,46.,0.1*415,0.1*self.sig[4],150.,150.],'f')
        self.sfrmin=N.array([0.5,0.5,0.5,0.24,0.24,.5])
        self.dL=N.array([3005.,3236.,3483.,3760.,3694.,3688.,3536.]) #Mpc/h
        #self.dA=N.array([5.02,5.14,5.25,5.35,5.33]) #kpc/arcsec/h
        self.dA=N.zeros(len(self.z),'f')
        for i in range(len(self.dA)):
            self.dA[i]=DA(self.z[i],h100)
        self.flux0=N.array([8.55,8.46,9.54,9.21,12.8,12.8,12.8])#flux ZP times 1d-17 in ergs/s/cm^2
        self.fluxmin=N.array([0.17,0.17,0.081,0.11,0.11,.11,.11])#3sig narrow-band flux in adu/s
        self.r200Mpc = N.zeros(len(self.z),'f')
        #self.r200Mpc = 1.73*self.sig/1000./N.sqrt(omegaL+omega0*(1+self.z)**3) #Mpc/h
        self.r200Mpc = r200(self.sig,self.z,h100)

        self.r200pix=self.r200Mpc*1000./self.dA/self.pscale #convert to kpc,arcsec, then pix
        self.mass = N.ones(len(self.r200Mpc),'f')

        self.xmax = N.array([791.,909.,822.,861.,973.,934.,938.],'f')
        self.ymax = N.array([821.,872.,769.,861.,973.,910.,901.],'f')
                             
    def calcstats(self):
        f_emission=n_emission/ntot
        f_burst=nburst/ntot
        print "Emission line = n_emission \n"
        print "Starburst = nburst \n"
        
        print "Total # gal   = ntot \n"
        nel=n_emission
        errfel=sqrt(nel+nel/ntot**2)/ntot
        print "Fraction of EL galaxies = %5.2f +/- %5.2f\n",f_emission,errfel
        errfs=sqrt(nburst+nburst/ntot**2)/ntot
        print "Fraction of starburst galaxies = %5.2f +/- %5.2f\n",f_burst,errfs 
        print "number of starforming galaxies = nsfr \n"


class Galaxy:
    def __init__(self):#individual galaxy properties
        #print "dude - a galaxy!"
        self.magj = []
        self.errmagj = []
        self.isoarea = []
        self.fj = []
        self.errfj = []
        self.x = []
        self.y = []
        self.fwhm = []
        self.flag = []
        self.seclass  = []
        self.f1j = []
        self.f2j = []
        self.f3j = []
        self.f4j = []
        self.f5j = []
        self.f6j = []
        self.f7j = []
        self.errf1j = []
        self.errf2j = []
        self.errf3j = []
        self.errf4j = []
        self.errf5j = []
        self.errf6j = []
        self.errf7j = []
        self.mag6j = []
        self.errmag6j = []
        self.magn = []
        self.errmagn = []
        self.fn = []
        self.errfn = []
        self.f1n = []
        self.f2n = []
        self.f3n = []
        self.f4n = []
        self.f5n = []
        self.f6n = []
        self.f7n = []
        self.errf1n = []
        self.errf2n = []
        self.errf3n = []
        self.errf4n = []
        self.errf5n = []
        self.errf6n = []
        self.errf7n = []
        self.mag6n = []
        self.errmag6n = []

        self.clusterid = []
        self.clusterr200 = []
        self.clustersig  = []
        self.clusterz = []
        self.clusterx = []
        self.clustery = []
        
        self.sfr = []
        self.errsfr = []
        self.ew = []
        self.errew = []
    def readjcat(self):
        name="testj.cat"
        j=0
        for line in open(name):
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            t=line.split()
            #print t[0]
            #    number 1 magiso 2 magisoerr 3 magauto 6 magautoerr 7 fluxiso 9 fluxisoerr 10 x 11 y 12 magbest 13 magbesterr 14 elong 15 ellip 16 fwhm 17 flag 18 class 19 fluxap1 20 fluxap2 21 fluxap3 22 fluxap1err 23 fluxap2err 24 fluxap3err 25
            pos = checkposition(float(t[10]),float(t[11]))
            #print pos,"fj = ", float(t[3])
            if pos < 1:
                continue
            #print "continuing...",j
            j=j+1
            self.magj.append(float(t[1]))#magiso
            self.errmagj.append(float(t[2]))
            self.isoarea.append(float(t[7]))
            self.fj.append(float(t[8]))
            self.errfj.append(float(t[9]))
            self.x.append(float(t[10]))
            self.y.append(float(t[11]))
            self.fwhm.append(float(t[16]))
            self.flag.append(float(t[17]))
            self.seclass.append(float(t[18]))
            self.f1j.append(float(t[19]))
            self.f2j.append(float(t[20]))
            self.f3j.append(float(t[21]))
            self.f4j.append(float(t[22]))
            self.f5j.append(float(t[23]))
            self.f6j.append(float(t[24]))
            self.f7j.append(float(t[25]))
            self.errf1j.append(float(t[26]))
            self.errf2j.append(float(t[27]))
            self.errf3j.append(float(t[28]))
            self.errf4j.append(float(t[39]))
            self.errf5j.append(float(t[30]))
            self.errf6j.append(float(t[31]))
            self.errf7j.append(float(t[32]))
            self.mag6j.append(float(t[38])) #mag in 2" aperture
            self.errmag6j.append(float(t[45]))

            self.clusterid.append(ncl)#keep track of cluster id
            self.clusterr200.append(c.r200pix[ncl])#keep track of cluster virial radius
            self.clustersig.append(c.sig[ncl])#keep track of cluster sig
            self.clusterz.append(c.z[ncl])#keep track of cluster sig
            self.clusterx.append(c.xc[ncl])#keep track of cluster sig
            self.clustery.append(c.yc[ncl])#keep track of cluster sig
    def readncat(self):
        #print "hey narrow"
        name="testn.cat"
        for line in open(name):
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            t=line.split()
            #print "blah"
            #    number 1 magiso 2 magisoerr 3 magauto 6 magautoerr 7 fluxiso 9 fluxisoerr 10 x 11 y 12 magbest 13 magbesterr 14 elong 15 ellip 16 fwhm 17 flag 18 class 19 fluxap1 20 fluxap2 21 fluxap3 22 fluxap1err 23 fluxap2err 24 fluxap3err 25
            pos = checkposition(float(t[10]),float(t[11]))
            if pos < 1:
                continue
            self.magn.append(float(t[1]))#magiso
            self.errmagn.append(float(t[2]))
            self.fn.append(float(t[8]))
            self.errfn.append(float(t[9]))
            self.f1n.append(float(t[19]))
            self.f2n.append(float(t[20]))
            self.f3n.append(float(t[21]))
            self.f4n.append(float(t[22]))
            self.f5n.append(float(t[23]))
            self.f6n.append(float(t[24]))
            self.f7n.append(float(t[25]))
            self.errf1n.append(float(t[26]))
            self.errf2n.append(float(t[27]))
            self.errf3n.append(float(t[28]))
            self.errf4n.append(float(t[29]))
            self.errf5n.append(float(t[30]))
            self.errf6n.append(float(t[31]))
            self.errf7n.append(float(t[32]))
            self.mag6n.append(float(t[38])) #mag in 2" aperture
            self.errmag6n.append(float(t[45]))
                
                
    def getediscs(self):
        self.vltv = []
        self.errvltv = []
        self.vltr = []
        self.errvltr = []
        self.vlti = []
        self.errvlti = []
        self.vltj = []
        self.errvltj = []
        self.vltk = []
        self.errvltk = []
        self.vltv2 = []
        self.errvltv2 = []
        self.vltr2 = []
        self.errvltr2 = []
        self.vlti2 = []
        self.errvlti2 = []
        self.vltj2 = []
        self.errvltj2 = []
        self.vltk2 = []
        self.errvltk2 = []
        self.vltv1 = []
        self.errvltv1 = []
        self.vltr1 = []
        self.errvltr1 = []
        self.vlti1 = []
        self.errvlti1 = []
        self.vltj1 = []
        self.errvltj1 = []
        self.vltk1 = []
        self.errvltk1 = []
        self.memb = []
        self.vlumlo = []
        self.vlum = []
        self.vlumhi = []
        self.Rlumlo = []
        self.Rlum = []
        self.Rlumhi = []
        self.name = []
        self.ediscsra = []
        self.BTiband = []
        self.errminBTiband = []
        self.errmaxBTiband = []
        for line in open("vltphot"):
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            t=line.split()
            #print "number of fields in vltphot = ", len(t)
            self.vltv.append(float(t[0]))
            self.errvltv.append(float(t[1]))
            self.vltr.append(float(t[2]))
            self.errvltr.append(float(t[3]))
            self.vlti.append(float(t[4]))
            self.errvlti.append(float(t[5]))
            self.vltj.append(float(t[6]))
            self.errvltj.append(float(t[7]))
            self.vltk.append(float(t[8]))
            self.errvltk.append(float(t[9]))
            self.vltv2.append(float(t[10]))
            self.errvltv2.append(float(t[11]))
            self.vltr2.append(float(t[12]))
            self.errvltr2.append(float(t[13]))
            self.vlti2.append(float(t[14]))
            self.errvlti2.append(float(t[15]))
            self.vltj2.append(float(t[16]))
            self.errvltj2.append(float(t[17]))
            self.vltk2.append(float(t[18]))
            self.errvltk2.append(float(t[19]))
            self.vltv1.append(float(t[20]))
            self.errvltv1.append(float(t[21]))
            self.vltr1.append(float(t[22]))
            self.errvltr1.append(float(t[23]))
            self.vlti1.append(float(t[24]))
            self.errvlti1.append(float(t[25]))
            self.vltj1.append(float(t[26]))
            self.errvltj1.append(float(t[27]))
            self.vltk1.append(float(t[28]))
            self.errvltk1.append(float(t[29]))
            self.memb.append(float(t[30]))
            self.Rlumlo.append(float(t[31]))
            self.Rlum.append(float(t[32]))
            self.Rlumhi.append(float(t[33]))
            self.vlumlo.append(float(t[34]))
            self.vlum.append(float(t[35]))
            self.vlumhi.append(float(t[36]))
            self.BTiband.append(float(t[37]))
            self.errminBTiband.append(float(t[38]))
            self.errmaxBTiband.append(float(t[39]))
        for line in open("names"):
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            t=line.split()
            self.name.append(t[0])
            try:
                self.ediscsra.append(float(t[1]))
            except:
                self.ediscsra.append(-999.)
        self.vltv =     N.array(self.vltv,'f')
        self.errvltv =  N.array(self.errvltv,'f')
        self.vltr    =  N.array(self.vltr,'f')
        self.errvltr =  N.array(self.errvltr,'f')
        self.vlti    =  N.array(self.vlti,'f')
        self.errvlti =  N.array(self.errvlti,'f')
        self.vltj    =  N.array(self.vltj,'f')
        self.errvltj =  N.array(self.errvltj,'f')
        self.vltk =     N.array(self.vltk,'f')
        self.errvltk =  N.array(self.errvltk,'f')
        self.vltv2 =    N.array(self.vltv2,'f')
        self.errvltv2 = N.array(self.errvltv2,'f')
        self.vltr2 =    N.array(self.vltr2,'f')
        self.errvltr2 = N.array(self.errvltr2,'f')
        self.vlti2 =    N.array(self.vlti2,'f')
        self.errvlti2 = N.array(self.errvlti2,'f')
        self.vltj2 =    N.array(self.vltj2,'f')
        self.errvltj2 = N.array(self.errvltj2,'f')
        self.vltk2 =    N.array(self.vltk2,'f')
        self.errvltk2 = N.array(self.errvltk2,'f')
        self.vltv1 =    N.array(self.vltv1,'f')
        self.errvltv1 = N.array(self.errvltv1,'f')
        self.vltr1 =    N.array(self.vltr1,'f')
        self.errvltr1 = N.array(self.errvltr1,'f')
        self.vlti1 =    N.array(self.vlti1,'f')
        self.errvlti1 = N.array(self.errvlti1,'f')
        self.vltj1 =    N.array(self.vltj1,'f')
        self.errvltj1 = N.array(self.errvltj1,'f')
        self.vltk1 =    N.array(self.vltk1,'f')
        self.errvltk1 = N.array(self.errvltk1,'f')
        self.memb =     N.array(self.memb,'i')
        self.ediscsra =     N.array(self.ediscsra,'f')
        self.Rlumlo =   N.array(self.Rlumlo,'f')
        self.Rlum =   N.array(self.Rlum,'f')
        self.Rlumhi =   N.array(self.Rlumhi,'f')
        self.vlumlo =   N.array(self.vlumlo,'f')
        self.vlum =   N.array(self.vlum,'f')
        self.vlumhi =   N.array(self.vlumhi,'f')
        self.BTiband =   N.array(self.BTiband,'f')
        self.errminBTiband =   N.array(self.BTiband,'f')
        self.errmaxBTiband =   N.array(self.BTiband,'f')
        self.sigma10 =   N.ones(len(self.vlumhi),'f') #surface density to 10th nearest neighbor
        self.sigma10flag =   N.ones(len(self.vlumhi),'f') #surface density to 10th nearest neighbor, background corrected
        tempname=[]
        tempsigma=[]
        #tempflag=[]
        for line in open("tsigma.cat"):
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            t=line.split()
            #tempname.append(t[1])
            #tempsigma.append(t[6])
            #tempflag.append(t[10])
            tempname.append(t[0])
            tempsigma.append(t[1])
        for i in range(len(self.sigma10)):
            for j in range(len(tempname)):
                t=tempname[j]
                if t.find(self.name[i]) > -1: 
                    self.sigma10[i]=float(tempsigma[j])
                    self.sigma10flag[i]=0.
                    #print ncl,"sigma10 = ",self.sigma10[i],N.log10(self.sigma10[i]),self.sigma10flag[i]
    def conv2array(self):
        self.magj = N.array(self.magj,'f')
        self.errmagj = N.array(self.errmagj,'f')
        self.isoarea = N.array(self.isoarea,'f')
        self.fj = N.array(self.fj,'f')
        self.errfj = N.array(self.errfj,'f')
        self.errfj = N.sqrt(self.isoarea)*c.jnoisea[ncl]*(1.+c.jnoiseb[ncl]*N.sqrt(self.isoarea))
        try:
            self.errmagj=2.5*N.log10(1. + self.errfj/self.fj)
        except:
            self.errmagj=99.
        self.x = N.array(self.x,'f')
        self.y = N.array(self.y,'f')
        self.fwhm = N.array(self.fwhm,'f')
        self.flag = N.array(self.flag,'f')
        self.seclass  = N.array(self.seclass,'f')
        self.f1j = N.array(self.f1j,'f')
        self.f2j = N.array(self.f2j,'f')
        self.f3j = N.array(self.f3j,'f')
        self.f4j = N.array(self.f4j,'f')
        self.f5j = N.array(self.f5j,'f')
        self.f6j = N.array(self.f6j,'f')
        self.f7j = N.array(self.f7j,'f')
        self.errf1j = N.array(self.errf1j,'f')
        self.errf2j = N.array(self.errf2j,'f')
        self.errf3j = N.array(self.errf3j,'f')
        self.errf4j = N.array(self.errf4j,'f')
        self.errf5j = N.array(self.errf5j,'f')
        self.errf6j = N.array(self.errf6j,'f')
        area=3.1415*(2./.18)**2
        self.errf6j = N.sqrt(area)*c.jnoisea[ncl]*(1.+c.jnoiseb[ncl]*N.sqrt(area))
        self.errf7j = N.array(self.errf7j,'f')
        self.mag6j = N.array(self.mag6j,'f')
        self.errmag6j = N.array(self.errmag6j,'f')
        self.errmag6j=2.5*N.log10(1. + abs(self.errf6j/self.f6j))
        self.magn = N.array(self.magn,'f')
        self.errmagn = N.array(self.errmagn,'f')

        self.fn = N.array(self.fn,'f')
        self.errfn = N.array(self.errfn,'f')
        self.errfn = N.sqrt(self.isoarea)*c.nbnoisea[ncl]*(1.+c.nbnoiseb[ncl]*N.sqrt(self.isoarea))
        self.f1n = N.array(self.f1n,'f')
        self.f2n = N.array(self.f2n,'f')
        self.f3n = N.array(self.f3n,'f')
        self.f4n = N.array(self.f4n,'f')
        self.f5n = N.array(self.f5n,'f')
        self.f6n = N.array(self.f6n,'f')
        self.f7n = N.array(self.f7n,'f')
        self.errf1n = N.array(self.errf1n,'f')
        self.errf2n = N.array(self.errf2n,'f')
        self.errf3n = N.array(self.errf3n,'f')
        self.errf4n = N.array(self.errf4n,'f')
        self.errf5n = N.array(self.errf5n,'f')
        self.errf6n = N.array(self.errf6n,'f')
        self.errf7n = N.array(self.errf7n,'f')
        self.mag6n = N.array(self.mag6n,'f')
        self.errmag6n = N.array(self.errmag6n,'f')

        self.clusterid = N.array(self.clusterid,'i')
        self.clusterr200 = N.array(self.clusterr200,'f')
        self.clustersig  = N.array(self.clustersig,'f')
        self.clusterz = N.array(self.clusterz,'f')
        self.clusterx = N.array(self.clusterx,'f')
        self.clustery = N.array(self.clustery,'f')
        

    def getratiosediscs(self):
        self.ratio = N.zeros(len(self.vltj),'f')
        self.errratio = N.zeros(len(self.vltj),'f')
        self.zroser = N.zeros(len(self.vltj),'f')
        self.zroserlo = N.zeros(len(self.vltj),'f')
        self.zroserhi = N.zeros(len(self.vltj),'f')
        self.zgreg = N.zeros(len(self.vltj),'f')
        self.zgreglo = N.zeros(len(self.vltj),'f')
        self.zgreghi = N.zeros(len(self.vltj),'f')
        self.starflag = N.zeros(len(self.vltj),'f')
        self.pclust = N.zeros(len(self.vltj),'f')
        i = 0
        for line in open("ratios.out"):
            length=len(line)
            t = line.split()
            #print length,i, len(self.ratio),len(line)
            #print line
            if (length < 50):
                #self.ratio[i]=substr(line,0,8)
                self.ratio[i]=float(line[0:8])
                self.errratio[i]=0.005
                self.zroser[i]=0.
                self.zroserlo[i]=0.
                self.zroserhi[i]=0.
                self.zgreg[i]=0.
                self.zgreglo[i]=0.
                self.zgreghi[i]=0.
                self.starflag[i]=0.
                self.pclust[i]=0.

            if (length > 50):
                #self.ratio[i]=substr(line,19,7)
                #self.errratio[i]=substr(line,28,7)
                self.ratio[i]=float(line[19:26])
                #self.errratio[i]=float(line[28:35])
                self.errratio[i]=.005
                self.zroser[i]=float(t[4])
                self.zroserlo[i]=float(t[5])
                self.zroserhi[i]=float(t[6])
                self.zgreg[i]=float(t[7])
                self.zgreglo[i]=float(t[8])
                self.zgreghi[i]=float(t[9])
                self.starflag[i]=float(t[13])
                self.pclust[i]=float(t[11])
            #reset ratios according to fits above
            if (self.zgreg[i] < .001):
                self.ratio[i] = c.defratio[ncl]
            if (self.zgreg[i] < 1.5):
                self.ratio[i] = c.ratiofitinter[ncl]-c.ratiofitslope[ncl]*self.zgreg[i]
            if (self.zgreg[i] > 1.5):
                self.ratio[i] = c.ratiofitinterhiz[ncl]-c.ratiofitslopehiz[ncl]*self.zgreg[i]
            i = i + 1
            

    def getspecmatch(self):
        self.specz = N.zeros(len(self.vltj),'f')
        self.ewo2 = N.zeros(len(self.vltj),'f')
        self.errewo2 = N.zeros(len(self.vltj),'f')
        self.spec = N.zeros(len(self.vltj),'f')
        self.spectype = N.zeros(len(self.vltj),'f')
        self.specmemb = N.zeros(len(self.vltj),'f')
        self.vltmatch = N.zeros(len(self.vltj),'f')
        i = 0
        for line in open("specmatch.out"):
            length=len(line)
            t = line.split()
            self.specz[i] = float(t[0])
            self.ewo2[i] = float(t[1])
            self.errewo2[i] = float(t[2])
            self.spec[i] = float(t[4])

            #self.spectype[i] = float(t[3])
            self.vltmatch[i] = float(t[5])
            #print ncl, i,self.specz[i], self.specmemb[i]
            if (self.spec[i] > 0) & (abs(self.specz[i] - c.z[ncl]) < 0.02):
                self.specmemb[i] = 1.
            i = i + 1
    def getspecmemb(self):#match spec memb for ms1054
        self.specspir = N.zeros(len(self.x),'f')
        spirx=[]
        spiry=[]
        for line in open("MS1054/ms1054_xy_spir_pisces.cat"):
            t = line.split()
            spirx.append(float(t[0]))
            spiry.append(float(t[1]))
        spirx=N.array(spirx,'f')
        spiry=N.array(spiry,'f')
        gxindex=N.arange(0,len(self.x),1)
        gxindexsort=N.take(gxindex,N.argsort(self.x))
        gxsort=N.take(self.x,N.argsort(self.x))
        gyindex=N.arange(0,len(self.y),1)
        gyindexsort=N.take(gyindex,N.argsort(self.y))
        gysort=N.take(self.y,N.argsort(self.y))
        delta=10.#matching radius in pixels
        for i in range(len(spirx)):
            xpos=spirx[i]
            match=findmatch(xpos,gxsort,delta)#look for cid in cluster list
            mindist=100.
            for j in match:
                k=gxindexsort[j]#index in unsorted x array
                d=N.sqrt((spirx[i]-self.x[k])**2+(spiry[i]-self.y[k])**2)
                if d < mindist:
                    kmatch=k
                    mindist=d
            self.specspir[kmatch]=1.

    def getratios(self):
        self.ratio = N.ones(len(self.fj),'f')
        self.errratio = N.zeros(len(self.fj),'f')
        self.ratio = c.defratio[ncl]*self.ratio
            
    def getnttname(self,nttname,nttrlum,nttrlumlo,nttrlumhi):
	self.intt=N.zeros(len(self.name),'i')
	self.inttmatchflag=N.zeros(len(self.name),'i')
	self.rlum=N.zeros(len(self.name),'f')
	self.rlumlo=N.zeros(len(self.name),'f')
	self.rlumhi=N.zeros(len(self.name),'f')
        self.Mr = N.zeros(len(self.Rlum),'f')
        self.errMrlo = N.zeros(len(self.Rlum),'f')
        self.errMrhi = N.zeros(len(self.Rlum),'f')

	for i in range(len(self.name)):
	    (self.intt[i],self.inttmatchflag[i])=matchname(self.name[i],nttname)
	    if (self.inttmatchflag[i] < 1):
		print "No match for ",self.name[i]
	    if (self.inttmatchflag[i] > 0) :
		#print "found a match ",self.name[i],nttname[self.intt[i]]
		self.rlum[i]=nttrlum[int(self.intt[i])]
		self.rlumlo[i]=nttrlumlo[int(self.intt[i])]
		self.rlumhi[i]=nttrlumhi[int(self.intt[i])]
                self.Mr[i] = 4.62 - 2.5*N.log10(self.rlum[i])-25.#sdss r but for now just swapping
                self.errMrlo[i] = self.Mr[i] - (4.62 - 2.5*N.log10(self.rlumlo[i])-25.)
                self.errMrhi[i] = (4.62 - 2.5*N.log10(self.rlumhi[i]) -25.) - self.Mr[i]
		print self.name[i],self.Mr[i],self.rlum[i],nttrlum[int(self.intt[i])]

	    
    def calcstuff(self):#calculate projected radial distance in pixels and arcsec
        self.dr = N.zeros(len(self.vltj),'f')
        self.dr = N.sqrt((self.x - self.clusterx)**2 + (self.y - self.clustery)**2)
        self.dtheta = self.dr*c.pscale[ncl]
        self.dkpc = N.zeros(len(self.dr),'f')
        self.dkpc = self.dtheta*c.dA[ncl] #projected distance in kpc
        self.dr200 = N.zeros(len(self.dr),'f')
        self.dr200 = self.dtheta*c.dA[ncl]/(c.r200Mpc[ncl]*1000) #projected distance in kpc
        self.MR = N.zeros(len(self.Rlum),'f')
        self.errMRlo = N.zeros(len(self.Rlum),'f')
        self.errMRhi = N.zeros(len(self.Rlum),'f')
        self.Mv = N.zeros(len(self.Rlum),'f')
        self.errMvlo = N.zeros(len(self.Rlum),'f')
        self.errMvhi = N.zeros(len(self.Rlum),'f')
        for i in range(len(self.Rlum)):
            if (self.Rlum[i] > 0):
                #self.MR[i] = 4.76 - 2.5*N.log10(self.Rlum[i])-25.
                #self.errMRlo[i] = self.MR[i] - (4.76 - 2.5*N.log10(self.Rlumlo[i])-25.)
                #self.errMRhi[i] = (4.76 - 2.5*N.log10(self.Rlumhi[i]) -25.) - self.MR[i] 
                #self.Mv[i] = 4.82 - 2.5*N.log10(self.vlum[i]) -25.
                #self.errMvlo[i] = self.Mv[i] - (4.82 - 2.5*N.log10(self.vlumlo[i]) -25. )
                #self.errMvhi[i] = (4.82 - 2.5*N.log10(self.vlumhi[i]) - 25. ) - self.Mv[i] 

                self.MR[i] = 4.28 - 2.5*N.log10(self.Rlum[i])-25.
                self.errMRlo[i] = self.MR[i] - (4.28 - 2.5*N.log10(self.Rlumlo[i])-25.)
                self.errMRhi[i] = (4.28 - 2.5*N.log10(self.Rlumhi[i]) -25.) - self.MR[i] 
                self.Mv[i] = 4.82 - 2.5*N.log10(self.vlum[i]) -25.
                self.errMvlo[i] = self.Mv[i] - (4.82 - 2.5*N.log10(self.vlumlo[i]) -25. )
                self.errMvhi[i] = (4.82 - 2.5*N.log10(self.vlumhi[i]) - 25. ) - self.Mv[i] 

    def calcsfr(self):
        sfconv=c.sfrconv[ncl]
        self.contsub = N.zeros(len(self.fj),'f')
        self.errcontsub = N.zeros(len(self.fj),'f')
        self.sfr = N.zeros(len(self.fj),'f')
        self.errsfr = N.zeros(len(self.fj),'f')
        self.errsfrtot = N.zeros(len(self.fj),'f')
        self.sfr1 = N.zeros(len(self.fj),'f')
        self.errsfr1 = N.zeros(len(self.fj),'f')
        self.sfr2 = N.zeros(len(self.fj),'f')
        self.errsfr2 = N.zeros(len(self.fj),'f')
        self.sfr3 = N.zeros(len(self.fj),'f')
        self.errsfr3 = N.zeros(len(self.fj),'f')
        print "length contsub,fn,ratio,fj = ",len(self.contsub),len(self.fn),len(self.ratio),len(self.fj),
        self.contsub = self.fn - self.ratio*self.fj
        self.errcontsub = N.sqrt(self.errfn**2 + (self.ratio*self.errfj)**2)
        #self.errcontsub = N.sqrt(self.isoarea)*c.noisea[ncl]*(1+ c.noiseb[ncl]*N.sqrt(self.isoarea))
        self.sfr = self.contsub*sfconv


        #self.errsfr = N.sqrt(self.errfn**2 + (self.ratio*self.errfj)**2 + (0.3*self.sfr)**2)*sfconv
        self.errsfr = self.errcontsub*sfconv
        self.fullerrsfr = N.sqrt(self.errsfr**2+(self.errratio*self.fj*sfconv)**2)

        self.errsfrtot = N.sqrt(self.errsfr**2+(self.sfr*c.errzp[ncl])**2+(self.errratio*self.fj)**2)
	self.sfr1 =(self.f1n-self.ratio*self.f1j)*sfconv
	self.errsfr1 = N.sqrt(self.errf1n**2 + (self.ratio*self.errf1j)**2)*sfconv
	self.sfr2 =(self.f2n-self.ratio*self.f2j)*sfconv
	self.errsfr2 = N.sqrt(self.errf2n**2 + (self.ratio*self.errf2j)**2)*sfconv
	self.sfr3=(self.f6n-self.ratio*self.f6j)*sfconv
	self.errsfr3 = N.sqrt(self.errf6n**2 + (self.ratio*self.errf6j)**2)*sfconv
        self.sn = N.array(len(self.sfr),'f')
        #self.sn = self.sfr / self.errsfr
        self.sn = self.contsub / self.errcontsub
        #calculate SFR if object at cluster z
        self.contsubc = N.zeros(len(self.fj),'f')
        self.errcontsubc = N.zeros(len(self.fj),'f')
        self.contsubc = self.fn - c.defratio[ncl]*self.fj
        self.errcontsubc = N.sqrt(self.errfn**2 + (c.defratio[ncl]*self.errfj)**2)
        self.sfrc = N.zeros(len(self.fj),'f')
        self.errsfrc = N.zeros(len(self.fj),'f')
        self.sfrc = sfconv*(self.fn - c.defratio[ncl]*self.fj)
        self.errsfrc = self.errcontsubc*sfconv
        self.fullerrsfrc = N.sqrt(self.errsfrc**2+(self.errratio*self.fj*sfconv)**2)
        #self.snc = N.array(len(self.sfr),'f')
        self.snc = self.contsubc / self.errcontsubc
        self.snc = self.sfrc / self.errsfrc

        self.sfrmax = sfconv*(self.fn-(c.defratio[ncl]-.005)*self.fj)

    def calcew(self):
	self.ew = (self.fn/self.fj-self.ratio)*dlambdaJ/(1+self.clusterz)
        #this error is wrong - recheck!!
	self.errew = N.sqrt( (dlambdaJ/self.fj*self.errfn)**2 + (dlambdaJ*self.fn/self.fj**2*self.errfj)**2)/(1+self.clusterz)
        self.ewc = (self.fn/self.fj-c.defratio[ncl])*dlambdaJ/(1+self.clusterz)
        self.errewc = N.sqrt( (dlambdaJ/self.fj*self.errfn)**2 + (dlambdaJ*self.fn/self.fj**2*self.errfj)**2)/(1+self.clusterz)
        self.signif = N.zeros(len(self.ew),'f')
        self.signifc = N.zeros(len(self.ew),'f')
        self.ewmax = (self.fn/self.fj-(c.defratio[ncl]-0.005))*dlambdaJ/(1+self.clusterz)
        for i in range(len(self.signif)):
            if ((self.ew[i] > ewmin) & (self.sn[i] > snmin)):
                self.signif[i] = 1
            if ((self.ewc[i] > ewmin) & (self.snc[i] > snmin)):
                self.signifc[i] = 1

    def finalsample(self):
        self.final = N.zeros(len(self.sfr),'f')#flag indicating galaxy is in final sample
        self.finalsf = N.zeros(len(self.sfr),'f')#flag indicating galaxy is in final sample w/sig sf
        self.nmemb= N.zeros(len(self.sfr),'f')
        self.nspecmemb= N.zeros(len(self.sfr),'f')
        self.nphotoznemission= N.zeros(len(self.sfr),'f')
        zcl=c.z[ncl]
        for i in range(len(self.final)):
            if (self.starflag[i] < 1):
                if (self.memb[i] > 0):
                    self.final[i] = 1.
                    self.nmemb[i]= 1
                    if (self.signifc[i] > 0):
                        self.finalsf[i] = 1.
                    continue
                if (self.memb[i] < 1) & (self.specmemb[i] > 0):
                    self.final[i] = 1.
                    self.nspecmemb[i]=1
                    if (self.signifc[i] > 0):
                        self.finalsf[i] = 1
                    continue
                if ((self.memb[i] < 1) & (self.specmemb[i] < 1) & (self.signifc[i] > 0) & (((self.zgreglo[i] < zcl) & (self.zgreghi[i] > zcl)) | (self.pclust[i] > 0.2))):
                    self.final[i] = 1.
                    self.finalsf[i] = 1.
                    self.nphotoznemission[i]= 1
                if ncl == 0: #pull out some obvious stars...
                    xystars = [(170.,666.3),(362.,714.),(640.,417.)]
                if ncl == 1:
                    xystars = [(879.,598.),(349.,620.)]
                if ncl == 2:
                    xystars = [(444.14,139.14)]
                for (xs,ys) in xystars:
                    #for i in range(len(self.final)):
                    d = sqrt((xs-self.x[i])**2 + (ys-self.y[i])**2)
                    if d < 4:
                        print ncl,"found match to star"
                        self.final[i]=0
                        self.finalsf[i]=0
                        self.nphotoznemission[i]=0
                        self.nspecmemb[i]=0
                        self.nmemb[i]=0

        self.nfinal = N.sum(N.compress(self.dkpc < (0.5*c.r200Mpc[ncl]*1000.),self.final))
        self.nfinalsf = N.sum(N.compress(self.dkpc < (0.5*c.r200Mpc[ncl]*1000.),self.finalsf))
        (self.sffrac,self.errsffrac) = ratioerror(self.nfinalsf,self.nfinal)
        self.nfinalsf40 = float(len(N.compress((self.final > 0) & (self.ewc > 40) & (self.finalsf> 0) & (self.dkpc < (0.5*c.r200Mpc[ncl]*1000.)),self.sfr)))
        (self.sffrac40,self.errsffrac40) = ratioerror(self.nfinalsf40,self.nfinal)
        self.nfinalsf100 = float(len(N.compress((self.final > 0) & (self.ewc > 100) & (self.finalsf > 0) & (self.dkpc < (0.5*c.r200Mpc[ncl]*1000.)),self.sfr)))
        (self.sffrac100,self.errsffrac100) = ratioerror(self.nfinalsf100,self.nfinal)
        self.totalsfr = N.sum(N.compress(self.finalsf > 0, self.sfrc))
        #self.totalsfrerr = 0.3*self.totalsfr
        self.totalsfrerr = N.sqrt(N.sum(N.compress(self.finalsf > 0, self.errsfrc)**2))
        self.totalsfrr200 = N.sum(N.compress((self.finalsf > 0) & (self.dkpc < (0.5*c.r200Mpc[ncl]*1000.)), self.sfrc))
        self.totalsfrr200err = 0.3*self.totalsfrr200
        #for comparison with SDSS clusters
        MRsdss=-20.38+5*N.log10(h100/.7)
        self.totalsfrsdss = N.sum(N.compress((self.finalsf > 0) & (self.dkpc < (0.5*c.r200Mpc[ncl]*1000.)) & (self.MR < MRsdss), self.sfrc))
        self.totalsfrsdsserr = N.sqrt(N.sum(N.compress((self.finalsf > 0) & (self.dkpc < (0.5*c.r200Mpc[ncl]*1000.)) & (self.MR < MRsdss), self.errsfrc)**2))
        #write out noise data
        noisedat = open("contsub-noise-area.dat",'w')
        for i in range(len(self.contsubc)):
            if self.final[i] > 0:
                noisedat.write("%8.4f %8.4f %8.3f \n" % (self.contsubc[i],self.errcontsubc[i],self.isoarea[i]))
                    #print self.contsubc[i],self.errcontsubc[i],self.isoarea[i]
        noisedat.close()
        self.dNN = N.zeros(len(self.Rlum),'f')
        for i in range(len(self.dNN)):
            temp=N.sqrt((self.x[i]-self.x)**2+(self.y[i]-self.y)**2)*c.pscale[ncl]
            temp=N.sort(N.compress(self.final > 0,temp))
            self.dNN[i]=temp[10]#distance to 5th NN, no edge effect accounted for

                                  
class NTT:
    def __init__(self):#individual galaxy properties
        #print "dude - a galaxy!"
        self.magj = []
        self.errmagj = []
        self.x = []
        self.y = []
        self.seclass  = []
        self.memb = []
        self.inmmtfield = []
        self.vlumlo = []
        self.vlum = []
        self.vlumhi = []
        self.Rlumlo = []
        self.Rlum = []
        self.Rlumhi = []
        self.rlumlo = []
        self.rlum = []
        self.rlumhi = []
	self.name=[]
	self.oldname=[]
        for line in open("vlt2mmtjcat.dat"):
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            t=line.split()
            # xmmt ymmt j(2") errj SEclass memb
            pos = checkposition(float(t[0]),float(t[1]))
            self.inmmtfield.append(float(pos))
            self.x.append(float(t[0]))
            self.y.append(float(t[1]))
            self.magj.append(float(t[2]))#magiso
            self.errmagj.append(float(t[3]))
            self.seclass.append(float(t[4]))
            self.memb.append(float(t[5]))
        for line in open("photoz.dat"):
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            t=line.split()
	    self.name.append(t[0])
	    self.oldname.append(t[1])
            self.vlumlo.append(float(t[38]))
            self.vlum.append(float(t[39]))
            self.vlumhi.append(float(t[40]))
            self.Rlumlo.append(float(t[41]))
            self.Rlum.append(float(t[42]))
            self.Rlumhi.append(float(t[43]))
            self.rlumlo.append(float(t[72]))#SDSS r luminosity
            self.rlum.append(float(t[73]))
            self.rlumhi.append(float(t[74]))

        self.inmmtfield = N.array(self.inmmtfield,'f')
        #self.x = N.compress(self.inmmtfield > .5,N.array(self.x,'f'))
        #self.y = N.compress(self.inmmtfield > .5,N.array(self.y,'f'))
        #self.magj = N.compress(self.inmmtfield > .5,N.array(self.magj,'f'))
        #self.errmagj = N.compress(self.inmmtfield > .5,N.array(self.errmagj,'f'))
        #self.seclass = N.compress(self.inmmtfield > .5,N.array(self.seclass,'f'))
        #self.memb = N.compress(self.inmmtfield > .5,N.array(self.memb,'f'))#

        #self.vlumlo = N.compress(self.inmmtfield > .5,N.array(self.vlumlo,'f'))
        #self.vlum = N.compress(self.inmmtfield > .5,N.array(self.vlum,'f'))
        #self.vlumhi = N.compress(self.inmmtfield > .5,N.array(self.vlumhi,'f'))
        #self.Rlumlo = N.compress(self.inmmtfield > .5,N.array(self.Rlumlo,'f'))
        #self.Rlum = N.compress(self.inmmtfield > .5,N.array(self.Rlum,'f'))
        #self.Rlumhi = N.compress(self.inmmtfield > .5,N.array(self.Rlumhi,'f'))
        #self.rlumlo = N.compress(self.inmmtfield > .5,N.array(self.rlumlo,'f'))
        #self.rlum = N.compress(self.inmmtfield > .5,N.array(self.rlum,'f'))
        #self.rlumhi = N.compress(self.inmmtfield > .5,N.array(self.rlumhi,'f'))
        #self.name = N.compress(self.inmmtfield > .5,self.name)

	self.x = N.array(self.x,'f')
        self.y = N.array(self.y,'f')
        self.magj = N.array(self.magj,'f')
        self.errmagj = N.array(self.errmagj,'f')
        self.seclass = N.array(self.seclass,'f')
        self.memb = N.array(self.memb,'f')

        self.vlumlo = N.array(self.vlumlo,'f')
        self.vlum = N.array(self.vlum,'f')
        self.vlumhi = N.array(self.vlumhi,'f')
        self.Rlumlo = N.array(self.Rlumlo,'f')
        self.Rlum = N.array(self.Rlum,'f')
        self.Rlumhi = N.array(self.Rlumhi,'f')
        self.rlumlo = N.array(self.rlumlo,'f')
        self.rlum = N.array(self.rlum,'f')
        self.rlumhi = N.array(self.rlumhi,'f')


        
        self.MR = N.zeros(len(self.Rlum),'f')
        self.errMRlo = N.zeros(len(self.Rlum),'f')
        self.errMRhi = N.zeros(len(self.Rlum),'f')
        self.Mv = N.zeros(len(self.Rlum),'f')
        self.errMvlo = N.zeros(len(self.Rlum),'f')
        self.errMvhi = N.zeros(len(self.Rlum),'f')
        self.Mr = N.zeros(len(self.Rlum),'f')
        self.errMrlo = N.zeros(len(self.Rlum),'f')
        self.errMrhi = N.zeros(len(self.Rlum),'f')
        for i in range(len(self.Rlum)):
            if (self.Rlum[i] > 0):
                #self.MR[i] = 4.76 - 2.5*N.log10(self.Rlum[i])-25.
                self.MR[i] = 4.28 - 2.5*N.log10(self.Rlum[i])-25.
                self.errMRlo[i] = self.MR[i] - (4.28 - 2.5*N.log10(self.Rlumlo[i])-25.)
                self.errMRhi[i] = (4.28 - 2.5*N.log10(self.Rlumhi[i]) -25.) - self.MR[i] 
                self.Mv[i] = 4.82 - 2.5*N.log10(self.vlum[i]) -25.
                self.errMvlo[i] = self.Mv[i] - (4.82 - 2.5*N.log10(self.vlumlo[i]) -25. )
                self.errMvhi[i] = (4.82 - 2.5*N.log10(self.vlumhi[i]) - 25. ) - self.Mv[i] 
                self.Mr[i] = 4.62 - 2.5*N.log10(self.rlum[i])-25.
                self.errMrlo[i] = self.MR[i] - (4.62 - 2.5*N.log10(self.rlumlo[i])-25.)
                self.errMrhi[i] = (4.62 - 2.5*N.log10(self.rlumhi[i]) -25.) - self.Mr[i] 
        

c = Cluster()
ncl2 = 3
#print c.sig
if (ncl2 < 4):
    ncl=0
    os.system("cp /Users/rfinn/clusters/final/CL1040/test*.cat .")
    os.system("cp /Users/rfinn/clusters/final/CL1040/vltphot .")
    os.system("cp /Users/rfinn/clusters/final/CL1040/ratios.out .")
    os.system("cp /Users/rfinn/clusters/final/CL1040/specmatch.out .")
    os.system("cp /Users/rfinn/clusters/final/CL1040/vlt2mmtjcat.dat .")
    os.system("cp /Users/rfinn/clusters/final/CL1040/names .")
    #os.system("cp /Users/rfinn/clusters/ediscs/catalogs/photoz/v23b/zcat.final.cl1040-1155vrijk.v2.3b.dat photoz.dat")
    os.system("cp /Users/rfinn/clusters/ediscs/catalogs/photoz/v241/zcat.final.cl1040-1155vrijk.v2.4.1.dat photoz.dat")
    os.system("cp /Users/rfinn/clusters/ediscs/roser/cl1040-1155.sigma10 tsigma.cat")
    #ncl=0
    print "CL1040"
    g0 = Galaxy()
    g0.readjcat()
    g0.readncat()
    g0.getediscs()
    g0.conv2array()
    g0.getratiosediscs()
    print len(g0.fj), len(g0.fn), len(g0.ratio)
    g0.getspecmatch()
    ntt0=NTT()
    g0.getnttname(ntt0.name,ntt0.rlum,ntt0.rlumlo,ntt0.rlumhi)
    g0.calcstuff()
    g0.calcsfr()
    g0.calcew()
    g0.finalsample()
    a1=N.compress(g0.finalsf > 0,g0.sfrc)
    b1=N.compress(g0.finalsf > 0,g0.sfrmax)
    c1=(b1-a1)/a1
    print "sum sfr (finalsf > 0)= ",N.sum(a1)
    print "sum sfr max = ",N.sum(b1)
    print "average increase = ",N.average(c1)
    print "# of gal w/SFR > 1 = ",len(N.compress((g0.finalsf > 0) & (g0.sfrc > 1),g0.sfrc))
    print "# of gal w/SFR > 1 and EW > 40= ",len(N.compress((g0.finalsf > 0) & (g0.sfrc > 1) & (g0.ewc > 40),g0.sfrc))
    print "names # of gal w/SFR > 1 and EW > 40= ",len(N.compress((g0.finalsf > 0) & (g0.sfrc > 1) & (g0.ewc > 40),g0.sfrc))
    output=open("CL1040/sfrgt1",'w')
    output2=open("CL1040/sfrlt1",'w')
    for i in range(len(g0.sfr)):
        if (g0.finalsf[i] > 0):
            output.write("%5.1f %5.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n" % (g0.x[i],g0.y[i],g0.ew[i],g0.errew[i],g0.sfr[i],g0.errsfr[i],g0.magj[i]))
        if (g0.sn[i] > snmin) & (g0.ew[i] < ewmin):
            output2.write("%5.1f %5.1f %6.1f %6.1f %6.1f %6.1f \n" % (g0.x[i],g0.y[i],g0.ew[i],g0.errew[i],g0.sfr[i],g0.errsfr[i]))
    output.close()
    output2.close()
if (ncl < 4):
    print "CL1054"
    ncl = 1
    os.system("cp contsub-noise-area.dat /Users/rfinn/clusters/final/CL1040/. ")
    os.system("cp /Users/rfinn/clusters/final/CL1054/test*.cat .")
    os.system("cp /Users/rfinn/clusters/final/CL1054/vltphot .")
    os.system("cp /Users/rfinn/clusters/final/CL1054/ratios.out .")
    os.system("cp /Users/rfinn/clusters/final/CL1054/specmatch.out .")
    os.system("cp /Users/rfinn/clusters/final/CL1054/vlt2mmtjcat.dat .")
    os.system("cp /Users/rfinn/clusters/final/CL1054/names .")
    #os.system("cp /Users/rfinn/clusters/ediscs/catalogs/photoz/v23b/zcat.final.cl1054-1245vrijk.v2.3b.dat photoz.dat")
    os.system("cp /Users/rfinn/clusters/ediscs/catalogs/photoz/v241/zcat.final.cl1054-1245vrijk.v2.4.1.dat photoz.dat")
    os.system("cp /Users/rfinn/clusters/ediscs/roser/cl1054-1245.sigma10 tsigma.cat")
    g1=Galaxy() 
    print "ncl = ",ncl
    g1.readjcat()
    g1.readncat()
    g1.getediscs()
    g1.conv2array()
    g1.getratiosediscs()
    g1.getspecmatch()
    ntt1=NTT()
    g1.getnttname(ntt1.name,ntt1.rlum,ntt1.rlumlo,ntt1.rlumhi)
    g1.calcstuff()
    g1.calcsfr()
    g1.calcew()
    #ntt1=NTT()
    g1.finalsample()
    a1=N.compress(g1.finalsf > 0,g1.sfrc)
    b1=N.compress(g1.finalsf > 0,g1.sfrmax)
    c1=(b1-a1)/a1
    print "sum sfr = ",N.sum(a1)
    print "sum sfr max = ",N.sum(b1)
    print "average increase = ",N.average(c1)
    print "# of gal w/SFR > 1 = ",len(N.compress((g1.finalsf > 0) & (g1.sfrc > 1),g1.sfrc))
    print "# of gal w/SFR > 1 and EW > 40= ",len(N.compress((g1.finalsf > 0) & (g1.sfrc > 1) & (g1.ewc > 40),g1.sfrc))
    os.system("cp contsub-noise-area.dat /Users/rfinn/clusters/final/CL1054/. ")
    # if (ncl == 2):
    ncl=2
    os.system("cp /Users/rfinn/clusters/final/CL1216/test*.cat .")
    os.system("cp /Users/rfinn/clusters/final/CL1216/vltphot .")
    os.system("cp /Users/rfinn/clusters/final/CL1216/ratios.out .")
    os.system("cp /Users/rfinn/clusters/final/CL1216/specmatch.out .")
    os.system("cp /Users/rfinn/clusters/final/CL1216/vlt2mmtjcat.dat .")
    os.system("cp /Users/rfinn/clusters/final/CL1216/names .")
    #os.system("cp /Users/rfinn/clusters/ediscs/catalogs/photoz/v23b/zcat.final.cl1216-1201vrijk.v2.3b.dat photoz.dat")
    os.system("cp /Users/rfinn/clusters/ediscs/catalogs/photoz/v241/zcat.final.cl1216-1201vrijk.v2.4.1.dat photoz.dat")
    os.system("cp /Users/rfinn/clusters/ediscs/roser/cl1216-1201.sigma10 tsigma.cat")
    g2=Galaxy()
    g2.readjcat()
    g2.readncat()
    g2.getediscs()
    g2.conv2array()
    g2.getratiosediscs()
    g2.getspecmatch()
    ntt2=NTT()
    g2.getnttname(ntt2.name,ntt2.rlum,ntt2.rlumlo,ntt2.rlumhi)
    g2.calcstuff()
    g2.calcsfr()
    g2.calcew()
    ntt2=NTT()
    g2.finalsample()
    os.system("cp contsub-noise-area.dat /Users/rfinn/clusters/final/CL1216/. ")
    a1=N.compress(g2.finalsf > 0,g2.sfrc)
    b1=N.compress(g2.finalsf > 0,g2.sfrmax)
    c1=(b1-a1)/a1
    print "sum sfr (finalsf > 0)= ",N.sum(a1)
    print "sum sfrmax = ",N.sum(b1)
    print "average increase = ",N.average(c1)
    a2=N.compress((g2.finalsf > 0) & (g2.sfrc > 2),g2.sfrc)
    b2=N.compress((g2.finalsf > 0) & (g2.sfrc > 2),g2.sfrmax)
    c2=(b2-a2)/a2
    print "sum sfr (> 2) = ",N.sum(a2)
    print "sum sfr max = ",N.sum(b2)
    print "average increase = ",N.average(c2)
    a2=N.compress((g2.finalsf > 0) & (g2.sfrc > 2),g2.ewc)
    b2=N.compress((g2.finalsf > 0) & (g2.sfrc > 2),g2.ewmax)
    c2=(b2-a2)/a2
    print "sum ew (sfr > 2) = ",N.sum(a2)
    print "sum ew max = ",N.sum(b2)
    print "average increase in ew = ",N.average(c2)
    print "# CL1216 of gal w/SFR > 1 = ",len(N.compress((g2.finalsf > 0) & (g2.sfrc > 1),g2.sfrc))
    print "# of gal w/SFR > 1 and EW > 40= ",len(N.compress((g2.finalsf > 0) & (g2.sfrc > 1) & (g2.ewc > 40),g2.sfrc))
    print "CL1040 - sf fracs (10, 40, 100A): %6.2f +/- %6.2f (%i/%i) %6.2f +/- %6.2f (%i) %6.3f +/- %6.3f (%i)" % (g0.sffrac,g0.errsffrac,g0.nfinalsf,g0.nfinal, g0.sffrac40,g0.errsffrac40,g0.nfinalsf40,g0.sffrac100,g0.errsffrac,g0.nfinalsf100)
    print "CL1040 - sf fracs (10, 40, 100A): %6.2f +/- %6.2f (%i/%i) %6.2f +/- %6.2f (%i) %6.3f +/- %6.3f (%i)" % (g1.sffrac,g1.errsffrac,g1.nfinalsf,g1.nfinal, g1.sffrac40,g1.errsffrac40,g1.nfinalsf40,g1.sffrac100,g1.errsffrac,g1.nfinalsf100)
    print "CL1040 - sf fracs (10, 40, 100A): %6.2f +/- %6.2f (%i/%i) %6.2f +/- %6.2f (%i) %6.3f +/- %6.3f (%i)" % (g2.sffrac,g2.errsffrac,g2.nfinalsf,g2.nfinal, g2.sffrac40,g2.errsffrac40,g2.nfinalsf40,g2.sffrac100,g2.errsffrac,g2.nfinalsf100)

    output=open("CL1216/sfrgt1",'w')
    output2=open("CL1216/sfrlt1",'w')
    for i in range(len(g2.sfr)):
        if (g2.finalsf[i] > 0):
            output.write("%5.1f %5.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n" % (g2.x[i],g2.y[i],g2.ew[i],g2.errew[i],g2.sfr[i],g2.errsfr[i],g2.magj[i]))
        if (g2.sn[i] > snmin) & (g2.ew[i] < ewmin):
            output2.write("%5.1f %5.1f %6.1f %6.1f %6.1f %6.1f \n" % (g2.x[i],g2.y[i],g2.ew[i],g2.errew[i],g2.sfr[i],g2.errsfr[i]))
    output.close()
    output2.close()
    ncl=3
    os.system("cp /Users/rfinn/clusters/final/CLJ0023/test*.cat .")
    g3=Galaxy()
    g3.readjcat()
    g3.readncat()
    g3.conv2array()
    g3.getratios()
    g3.calcsfr()
    g3.calcew()
    g3.finalsf=N.zeros(len(g3.ewc),'f')
    g3.name=[]
    for i in range(len(g3.finalsf)):
	    g3.name.append(str(i))
	    if (g3.ewc[i] > ewmin):
		    sn=g3.contsub[i]/g3.errcontsub[i]
		    if (sn > snmin):
			    g3.finalsf[i]=1

if (ncl == 4):

    os.system("cp /Users/rfinn/clusters/final/RXJ0152/test*.cat .")
    g4=Galaxy()
    g4.readjcat()
    g4.readncat()
    g4.conv2array()
    g4.getratios()
    g4.calcsfr()
    g4.calcew()
    output=open("rxjfinalxy.cat",'w')
    output2=open("rxjfinalsfxy.cat",'w')
    for i in range(len(g4.x)):
        output.write("%8.3f %8.3f \n" % (g4.x[i],g4.y[i]))
        if (g4.sn[i] > snmin) & (g4.ew[i] > ewmin):
            output2.write("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n" % (g4.x[i],g4.y[i],g4.sfr[i],g4.errsfr[i],g4.ew[i],g4.errew[i]))
    output.close()
    output2.close()
    os.system("cp rxjfinalxy.cat /Users/rfinn/clusters/final/RXJ0152/.")
    output=open("RXJ0152/sfrgt1",'w')
    output2=open("RXJ0152/sfrlt1",'w')
    for i in range(len(g4.sfr)):
        if (g4.sn[i] > snmin) & (g4.ew[i] > ewmin):
            output.write("%5.1f %5.1f %6.1f %6.1f %6.1f %6.1f %5.2f %5.2f\n" % (g4.x[i],g4.y[i],g4.ew[i],g4.errew[i],g4.sfr[i],g4.errsfr[i],g4.sn[i],g4.magj[i]))
            print("%5.1f %5.1f %6.1f %6.1f %6.1f %6.1f %5.2f %5.2f\n" % (g4.x[i],g4.y[i],g4.ewc[i],g4.errewc[i],g4.sfr[i],g4.errsfr[i],g4.sn[i],g4.magj[i]))
        if (g4.sn[i] > snmin) & (g4.ew[i] < ewmin):
            output2.write("%5.1f %5.1f %6.1f %6.1f %6.1f %6.1f \n" % (g4.x[i],g4.y[i],g4.ew[i],g4.errew[i],g4.sfr[i],g4.errsfr[i]))
    output.close()
    output2.close()

#ncl = 5
if (ncl == 5):
    ewmin=10.
    snmin=2.
    #os.system("cp /Volumes/giga/clusters/final/MS1054/test*.cat .")
    os.system("cp MS1054/test*.cat .")
    g5=Galaxy()
    g5.readjcat()
    g5.readncat()
    g5.conv2array()
    g5.getratios()
    g5.calcsfr()
    g5.calcew()
    g5.getspecmemb()
    psplotinit("ms1054ratioj.ps")
    plotratioj(g5.magj,g5.fj,g5.errfj,g5.fn,g5.errfn)
    ppgplot.pgend()
    psplotinit("ms1054snsfr.ps")
    plotsnsfr(g5.sn,g5.sfr,g5.errsfr)
    ppgplot.pgend()
    output=open("MS1054/sfrgt1",'w')
    output2=open("MS1054/sfrlt1",'w')
    g5.finalsf=N.zeros(len(g5.sfr),'f')
    for i in range(len(g5.sfr)):
        if (g5.sn[i] > snmin) & (g5.ewc[i] > ewmin):
            g5.finalsf[i]=1.
            output.write("%5.1f %5.1f %6.1f %6.1f %6.1f %6.1f %5.2f %5.2f\n" % (g5.x[i],g5.y[i],g5.ewc[i],g5.errewc[i],g5.sfr[i],g5.errsfr[i],g5.sn[i],g5.magj[i]))
        if (g5.sn[i] > snmin) & (g5.ewc[i] < ewmin):
            g5.finalsf[i]=0.
            output2.write("%5.1f %5.1f %6.1f %6.1f %6.1f %6.1f \n" % (g5.x[i],g5.y[i],g5.ew[i],g5.errew[i],g5.sfr[i],g5.errsfr[i]))
    output.close()
    output2.close()
    psplotinit("ms1054fitzp.ps")
    plotfitzp(g5.x,g5.y,g5.f4j,g5.errf4j,ncl)
    ppgplot.pgend()
    plotewjms(g5.magj,g5.ew,g5.errew,g5.finalsf,g5.specspir)
    plotsfrjms(g5.magj,g5.sfr,g5.errsfr,g5.finalsf,g5.specspir)
    psplotinit("ms1054histew.ps")
    plothistew(g5.ew)
    ppgplot.pgend()
    g5.final=N.ones(len(g5.sfr),'f')
    g5.finalsf=N.zeros(len(g5.sfr),'f')
    psplotinit("ms1054plotgal.ps")
    plotgal(g5.x,g5.y,g5.final,g5.finalsf,5)
    ppgplot.pgend()
if (ncl == 6):
    os.system("cp ~/field/final/HDFN1/test*.cat .")
    g6=Galaxy()
    g6.readjcat()
    g6.readncat()
    g6.conv2array()
    g6.getratios()
    g6.calcsfr()
    g6.calcew()
    psplotinit("HDFN1ewj.ps")
    plotewj(g6.magj,g6.ew,g6.errew)
    ppgplot.pgend()
    psplotinit("HDFN1sfrj.ps")
    plotsfrj(g6.magj,g6.sfr,g6.errsfr)
    ppgplot.pgend()
    psplotinit("HDFN1snsfr.ps")
    plotsnsfr(g6.sn,g6.sfr,g6.errsfr)
    ppgplot.pgend()
    psplotinit("HDFN1histew.ps")
    plothistew(g6.ew)
    ppgplot.pgend()
    psplotinit("HDFN1ratioj.ps")
    plotratioj(g6.magj,g6.fj,g6.errfj,g6.fn,g6.errfn)
    ppgplot.pgend()
    g6.final=N.ones(len(g6.sfr),'f')
    g6.finalsf=N.zeros(len(g6.sfr),'f')
    output=open("~/field/final/HDFN1/sfrgt1",'w')
    output2=open("~/field/final/HDFN1/sfrlt1",'w')
    for i in range(len(g6.sfr)):
        if (g6.sn[i] > snmin) & (g6.ewc[i] > ewmin):
            g6.finalsf[i]=1.
            output.write("%5.1f %5.1f %6.1f %6.1f %6.1f %6.1f %5.2f %5.2f\n" % (g6.x[i],g6.y[i],g6.ewc[i],g6.errewc[i],g6.sfr[i],g6.errsfr[i],g6.sn[i],g6.magj[i]))
        if (g6.sn[i] > snmin) & (g6.ewc[i] < ewmin):
            g6.finalsf[i]=0.
            output2.write("%5.1f %5.1f %6.1f %6.1f %6.1f %6.1f \n" % (g6.x[i],g6.y[i],g6.ew[i],g6.errew[i],g6.sfr[i],g6.errsfr[i]))
    output.close()
    output2.close()
    psplotinit("HDFN1plotgal.ps")
    plotgal(g6.x,g6.y,g6.final,g6.finalsf,6)
    ppgplot.pgend()
    psplotinit("HDFN1fitzp.ps")
    plotfitzp(g6.x,g6.y,g6.f4j,g6.errf4j,ncl)
    ppgplot.pgend()

psplotinit("cl1216snsfr.ps")
plotsnsfr(g3.sn,g3.sfr,g3.errsfr)
ppgplot.pgend()


ppgplot.pgbeg("all.ps/vcps",2,2)
ppgplot.pgsch(2.) #font size
ppgplot.pgslw(4)  #line width
#ppgplot.pgpage()
plothistsfr()
#ppgplot.pgpage()
x=N.compress((g0.memb>0) & (g0.flag < 1),g0.magj)
y=N.compress((g0.memb>0) & (g0.flag < 1),(g0.fn/g0.fj))
ratio=0.0625
plotnbjvj(x,y,ratio)
#ppgplot.pgpage()
x=N.compress((g1.memb>0) & (g1.flag < 1),g1.magj)
y=N.compress((g1.memb>0) & (g1.flag < 1),(g1.fn/g1.fj))
ratio=0.0685
plotnbjvj(x,y,ratio)

#ppgplot.pgpage()
x=N.compress((g2.memb>0) & (g2.flag < 1),g2.magj)
y=N.compress((g2.memb>0) & (g2.flag < 1),(g2.fn/g2.fj))
ratio=0.0711
plotnbjvj(x,y,ratio)
#ppgplot.pgpage()
plothistj()

plotgal(g0.x,g0.y,g0.final,g0.finalsf,0)
plotgal(g1.x,g1.y,g1.final,g1.finalsf,1)
plotgal(g2.x,g2.y,g2.final,g2.finalsf,2)

plotnkall()

ppgplot.pgend()


#individual plots
psplotinit("cl1040sumsfrmag.ps")
print "CL1040"
magmmt = g0.mag6j
sfrmmt = g0.sfr
membmmt = g0.memb
magntt = ntt0.magj
membntt = ntt0.memb
signif = g0.signif
cl1040compl = plotsumsfrmag(magmmt,sfrmmt,membmmt,magntt,membntt,signif)
xlsum = 21.8
ylsum = 32.
ppgplot.pgtext(xlsum,ylsum,"(a) CL1040")
ppgplot.pgend()

psplotinit("cl1054sumsfrmag.ps")
print "CL1054-12"
magmmt = g1.mag6j
sfrmmt = g1.sfr
membmmt = g1.memb
magntt = ntt1.magj
membntt = ntt1.memb
signif = g1.signif
cl1054compl = plotsumsfrmag(magmmt,sfrmmt,membmmt,magntt,membntt,signif)
ppgplot.pgtext(xlsum-.6,ylsum,"(b) CL1054-12")
ppgplot.pgend()

psplotinit("cl1216sumsfrmag.ps")
print "CL1216"
magmmt = g2.mag6j
sfrmmt = g2.sfr
membmmt = g2.memb
magntt = ntt2.magj
membntt = ntt2.memb
signif = g2.signif
cl1216compl = plotsumsfrmag(magmmt,sfrmmt,membmmt,magntt,membntt,signif)
ppgplot.pgtext(xlsum,ylsum,"(c) CL1216")
ppgplot.pgend()
os.system("cp *sumsfrmag.ps /Users/rfinn/clusters/papers/paper2/.")
os.system("cp *sumsfrmag.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")


psplotinit("ngalsfr.ps")
plothistsfr()
ppgplot.pgend()
os.system("cp ngalsfr.ps /Users/rfinn/clusters/papers/paper2/.")
os.system("cp ngalsfr.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

psplotinit("histewall.ps")
plothistewall()
ppgplot.pgend()
os.system("cp histewall.ps /Users/rfinn/clusters/papers/paper2/.")
os.system("cp histewall.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")


psplotinit("cl1040plotgal.ps")
plotgal(g0.x,g0.y,g0.final,g0.finalsf,0)
psplotinit("cl1054plotgal.ps")
plotgal(g1.x,g1.y,g1.final,g1.finalsf,1)
psplotinit("cl1216plotgal.ps")
plotgal(g2.x,g2.y,g2.final,g2.finalsf,2)
ppgplot.pgend()
os.system("cp *plotgal.ps /Users/rfinn/clusters/papers/paper2/.")
os.system("cp *plotgal.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

psplotinit("nkall.ps")
plotnkall()
ppgplot.pgend()
os.system("cp nkall.ps /Users/rfinn/clusters/papers/paper2/.")
os.system("cp nkall.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

psplotinit("nMvstarburst.ps")
plotnkstarburst()
ppgplot.pgend()
os.system("cp nMvstarburst.ps /Users/rfinn/clusters/papers/paper2/.")
os.system("cp nMvstarburst.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

psplotinit("sfrkall.ps")
plotsfrkall()
ppgplot.pgend()
os.system("cp sfrkall.ps /Users/rfinn/clusters/papers/paper2/.")

psplotinit("sfrjall.ps")
plotsfrjall()
ppgplot.pgend()
os.system("cp sfrjall.ps /Users/rfinn/clusters/papers/paper2/.")
os.system("cp sfrjall.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

psplotinit("magisoarea.ps")
plotmagisoarea()
ppgplot.pgend()
#os.system("cp sfrjall.ps /Users/rfinn/clusters/papers/paper2/.")

psplotinit("sfrewall.ps")
plotsfrewall()
ppgplot.pgend()
os.system("cp sfrewall.ps /Users/rfinn/clusters/papers/paper2/.")

psplotinit("sfrjallnsf.ps")
plotsfrjallnsf()
ppgplot.pgend()

psplotinit("ewkall.ps")
plotewkall()
ppgplot.pgend()
os.system("cp ewkall.ps /Users/rfinn/clusters/papers/paper2/.")
os.system("cp ewkall.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

psplotinit("ewjall.ps")
plotewjall()
ppgplot.pgend()
os.system("cp ewjall.ps /Users/rfinn/clusters/papers/paper2/.")
os.system("cp ewjall.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

psplotinit("sfrdall.ps")
plotsfrdall()
ppgplot.pgend()
os.system("cp sfrdall.ps /Users/rfinn/clusters/papers/paper2/.")
os.system("cp sfrdall.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")
psplotinit("ewdall.ps")
plotewdall()
ppgplot.pgend()
os.system("cp ewdall.ps /Users/rfinn/clusters/papers/paper2/.")
os.system("cp ewdall.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

psplotinit("sffracd.ps")
plotsffracd()
ppgplot.pgend()
os.system("cp sffracd.ps /Users/rfinn/clusters/papers/paper2/.")

psplotinit("sfrdr200all.ps")
plotsfrdr200all()
ppgplot.pgend()
os.system("cp sfrdr200all.ps /Users/rfinn/clusters/papers/paper2/.")
os.system("cp sfrdr200all.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")
psplotinit("ewdr200all.ps")
plotewdr200all()
ppgplot.pgend()
os.system("cp ewdr200all.ps /Users/rfinn/clusters/papers/paper2/.")
os.system("cp ewdr200all.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")
#plot in terms of distance to nearest neighbor
psplotinit("sfrdNNall.ps")
plotsfrdNNall()
ppgplot.pgend()
os.system("cp sfrdNNall.ps /Users/rfinn/clusters/papers/paper2/.")
psplotinit("ewdNNall.ps")
plotewdNNall()
ppgplot.pgend()
os.system("cp ewdNNall.ps /Users/rfinn/clusters/papers/paper2/.")

psplotinit("sfrsigma10all.ps")
plotsfrsigma10all()
ppgplot.pgend()
os.system("cp sfrsigma10all.ps /Users/rfinn/clusters/papers/paper2/.")
psplotinit("ewsigma10all.ps")
plotewsigma10all()
ppgplot.pgend()
os.system("cp ewsigma10all.ps /Users/rfinn/clusters/papers/paper2/.")

psplotinit("sffracdr200.ps")
plotsffracdr200()
ppgplot.pgend()
os.system("cp sffracdr200.ps /Users/rfinn/clusters/papers/paper2/.")
os.system("cp sffracdr200.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

psplotinit("sffracsigma10.ps")
plotsffracsigma10()
ppgplot.pgend()
os.system("cp sffracsigma10.ps /Users/rfinn/clusters/papers/paper2/.")
os.system("cp sffracsigma10.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

plotsfrlit()

psplotinit("cl1040haz.ps")
z = N.compress((g0.specmemb > 0) & (g0.final > 0),g0.specz)
ew = N.compress((g0.specmemb > 0) & (g0.final > 0),g0.ew) 
errew = N.compress((g0.specmemb > 0) & (g0.final > 0),g0.errew)
plotspecewz(z,ew,errew,0)
ppgplot.pgend()
psplotinit("cl1054haz.ps")
z = N.compress((g1.specmemb > 0) & (g1.final > 0),g1.specz)
ew = N.compress((g1.specmemb > 0) & (g1.final > 0),g1.ew) 
errew = N.compress((g1.specmemb > 0) & (g1.final > 0),g1.errew)
plotspecewz(z,ew,errew,1)
ppgplot.pgend()

psplotinit("cl1216haz.ps")
z = N.compress((g2.specmemb > 0) & (g2.final > 0),g2.specz)
ew = N.compress((g2.specmemb > 0) & (g2.final > 0),g2.ew) 
errew = N.compress((g2.specmemb > 0) & (g2.final > 0),g2.errew)
plotspecewz(z,ew,errew,2)
ppgplot.pgend()
os.system("cp *haz.ps /Users/rfinn/clusters/papers/paper2/.")
plotspecewzall()
os.system("cp haewzall.ps /Users/rfinn/clusters/papers/paper2/.")
os.system("cp haewzall.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")
plotspecewzall2()
os.system("cp haewzall2.ps /Users/rfinn/clusters/papers/paper2/.")
os.system("cp haewzall2.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

#compare EW[OII] vs EW[Halpha]
psplotinit("ewhao2finalsf.ps")
temp=g0.finalsf
ew = N.compress((g0.specmemb > 0) & (temp > 0),g0.ewc) 
errew = N.compress((g0.specmemb > 0) & (temp > 0),g0.errewc)
ewo2 = -1.*N.compress((g0.specmemb > 0) & (temp > 0),g0.ewo2) 
errewo2 = N.compress((g0.specmemb > 0) & (temp > 0),g0.errewo2)
print len(g0.specmemb),len(temp),len(g0.name)
#name = N.compress((g0.specmemb > 0) & (temp > 0),g0.name)
#names = N.compress((ew > 40),name)
#print "Spec memb w/EW > 40", names 

ppgplot.pgslw(6)
#psplotinit("sfrmcl.ps")
ppgplot.pgbox("",0.0,0,"L",0.0,0)
ppgplot.pgenv(-0,85.,-8.,50.,0,0)
ppgplot.pglab(ewlabel,ewo2label,"")
ppgplot.pgpt(ew,ewo2,3)
errorx(ew,ewo2,errew)
errory(ew,ewo2,errewo2)

temp=g1.finalsf
ew = N.compress((g1.specmemb > 0) & (temp > 0),g1.ewc) 
errew = N.compress((g1.specmemb > 0) & (temp > 0),g1.errewc)
ewo2 = -1.*N.compress((g1.specmemb > 0) & (temp > 0),g1.ewo2) 
errewo2 = N.compress((g1.specmemb > 0) & (temp > 0),g1.errewo2)
#name = N.compress((g1.specmemb > 0) & (temp > 0),g1.name)
#names = N.compress(ew > 40,name)
#print "Spec memb w/EW > 40",names 

ppgplot.pgpt(ew,ewo2,17)
errorx(ew,ewo2,errew)
errory(ew,ewo2,errewo2)

temp=g2.finalsf
ew = N.compress((g2.specmemb > 0) & (temp > 0),g2.ewc) 
errew = N.compress((g2.specmemb > 0) & (temp > 0),g2.errewc)
ewo2 = -1.*N.compress((g2.specmemb > 0) & (temp > 0),g2.ewo2) 
errewo2 = N.compress((g2.specmemb > 0) & (temp > 0),g2.errewo2)

ppgplot.pgpt(ew,ewo2,7)
errorx(ew,ewo2,errew)
errory(ew,ewo2,errewo2)
x=N.arange(-300,400,10)
y=0.4*x
ppgplot.pgline(x,y)

xlabel = 50.
ylabel = 45.
ystep = 4.
dy=1.
dxl=1
dxr=3.
ppgplot.pgslw(deflw)  #line width
ppgplot.pgtext(xlabel,ylabel,"CL1040")
xlin = N.array([xlabel-dxr],'f')
ylin = N.array([ylabel+dy],'f')
ppgplot.pgslw(4)  #line width
ppgplot.pgpt(xlin,ylin,3)

ylabel = ylabel - ystep
xlin = N.array([xlabel-dxr],'f')
ylin = N.array([ylabel+dy],'f')
ppgplot.pgpt(xlin,ylin,17)
ppgplot.pgsls(1)
ppgplot.pgslw(deflw)
ppgplot.pgtext(xlabel,ylabel,"CL1054-12")

ylabel = ylabel - ystep
ppgplot.pgslw(deflw)  #line width
ppgplot.pgtext(xlabel,ylabel,"CL1216")
xlin = N.array([xlabel-dxr],'f')
ylin = N.array([ylabel+dy],'f')
ppgplot.pgslw(2)  #line width
ppgplot.pgpt(xlin,ylin,7)
ppgplot.pgend()
os.system("cp ewhao2finalsf.ps /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

ew=g0.ewc
name=g0.name

starburstname= N.compress

bootfield()
plotcumulativesfrew()
plotcumulativesfrk()
print "TEST 1, g0.MR[9]=",g0.MR[9],g0.sfrc[9],g0.Mr[9]
plotcumulativesfrMR()
print "TEST 2, g0.MR[2]=",g0.MR[9],g0.sfrc[9],g0.Mr[9]
plotcompmmtedij()
print "TEST 3, g0.MR[2]=",g0.MR[9],g0.sfrc[9],g0.Mr[9]
dostats()
print "TEST 4, g0.MR[2]=",g0.MR[9],g0.sfrc[9],g0.Mr[9]
writeapjtable(g0.name,g0.x,g0.y,g0.fn,g0.errfn,g0.fj,g0.errfj,g0.ratio,g0.errratio,g0.contsubc,g0.errcontsubc,g0.ewc,g0.errewc,g0.sfrc,g0.errsfrc,g0.final,g0.finalsf,"cl1040",0,g0.magj,g0.errmagj,g0.ediscsra)
print "TEST 5, g0.MR[2]=",g0.MR[9],g0.sfrc[9],g0.Mr[9]
writeapjtable(g1.name,g1.x,g1.y,g1.fn,g1.errfn,g1.fj,g1.errfj,g1.ratio,g1.errratio,g1.contsubc,g1.errcontsubc,g1.ewc,g1.errewc,g1.sfrc,g1.errsfrc,g1.final,g1.finalsf,"cl1054",1,g1.magj,g1.errmagj,g1.ediscsra)

writeapjtable(g2.name,g2.x,g2.y,g2.fn,g2.errfn,g2.fj,g2.errfj,g2.ratio,g2.errratio,g2.contsubc,g2.errcontsubc,g2.ewc,g2.errewc,g2.sfrc,g2.errsfrc,g2.final,g2.finalsf,"cl1216",2,g2.magj,g2.errmagj,g2.ediscsra)
print "TEST 6, g0.MR[2]=",g0.MR[9],g0.sfrc[9],g0.Mr[9]
writesdsstable(g0.x,g0.y,g0.ewc,g0.errewc,g0.sfrc,g0.errsfrc,g0.finalsf,"cl1040",0,g0.Mr,g0.name)
print "TEST 7, g0.MR[2]=",g0.MR[9],g0.sfrc[9],g0.Mr[9]
writesdsstable(g1.x,g1.y,g1.ewc,g1.errewc,g1.sfrc,g1.errsfrc,g1.finalsf,"cl1054",1,g1.Mr,g1.name)
writesdsstable(g2.x,g2.y,g2.ewc,g2.errewc,g2.sfrc,g2.errsfrc,g2.finalsf,"cl1216",2,g2.Mr,g2.name)
writesdsstable(g3.x,g3.y,g3.ewc,g3.errewc,g3.sfrc,g3.errsfrc,g3.finalsf,"clj0023",3,g3.magj,g3.name)
print "TEST 8, g0.MR[2]=",g0.MR[9],g0.sfrc[9],g0.Mr[9]
os.system("cp *apj.tex /Users/rfinn/clusters/papers/paper2/.")
os.system("cp *apj.tex /Users/rfinn/clusters/papers/paper2/submit/resubmit/.")

#plotnsfsim()
#plotnsfsim2()
#plotnsfsim3()
#plotnsfsim4()

#calculate limiting SFR and J magnitude
for ncl in range(6):#loop over clusters
    minarea=12.
    if ncl == 5:
        minarea = 10.
    minarea=N.sqrt(minarea)
    print c.fullclustername[ncl]
    errj = (minarea)*c.jnoisea[ncl]*(1.+c.jnoiseb[ncl]*minarea)
    errn = (minarea)*c.nbnoisea[ncl]*(1.+c.nbnoiseb[ncl]*minarea)    
    errcontsubc = N.sqrt(errn**2 + (c.defratio[ncl]*errj)**2)
    print "1 sigma nb flux = ",errn,errn*c.sfrconv[ncl]
    print "1 sigma cont sub flux = ",errcontsubc
    print "3 sigma cont sub flux = ",3*errcontsubc
    print "1 sigma SFR = ",errcontsubc*c.sfrconv[ncl]
    print "2 sigma SFR = ",2*errcontsubc*c.sfrconv[ncl]
    minarea = N.sqrt(3.1415*(2./.18)**2)
    errj = (minarea)*c.jnoisea[ncl]*(1.+c.jnoiseb[ncl]*minarea)
    print "1 sigma J in 2arcsec ap= ",c.jzp[ncl] - 2.5*N.log10(errj) 

for ncl in range(6):#loop over clusters
    if ncl == 0:
        final=g0.final
        isoarea=g0.isoarea
    if ncl == 1:
        final=g1.final
        isoarea=g1.isoarea
    if ncl == 2:
        final=g2.final
        isoarea=g2.isoarea
    if ncl == 3:
        final=N.ones(len(g3.isoarea),'f')
        isoarea=g3.isoarea
    if ncl == 4:
        final=N.ones(len(g4.isoarea),'f')
        isoarea=g4.isoarea
    if ncl == 5:
        final=N.ones(len(g5.isoarea),'f')
        isoarea=g5.isoarea
    minarea=N.average(N.compress(final>0,isoarea))
    minarea=pylab.median(N.compress(final>0,isoarea))
    minarea=N.sqrt(minarea)
    print c.fullclustername[ncl]
    errj = (minarea)*c.jnoisea[ncl]*(1.+c.jnoiseb[ncl]*minarea)
    errn = (minarea)*c.nbnoisea[ncl]*(1.+c.nbnoiseb[ncl]*minarea)    
    errcontsubc = N.sqrt(errn**2 + (c.defratio[ncl]*errj)**2)
    print "ave isoarea  = ",minarea**2.
    print "ave 1 sigma nb flux = ",errn,errn*c.sfrconv[ncl]

    print "ave 1 sigma cont sub flux = ",errcontsubc
    print "ave 3 sigma cont sub flux = ",3*errcontsubc
    print "ave 1 sigma SFR = ",errcontsubc*c.sfrconv[ncl]
    print "ave 3 sigma SFR = ",3*errcontsubc*c.sfrconv[ncl]
    minarea = N.sqrt(3.1415*(2./.18)**2)
    errj = (minarea)*c.jnoisea[ncl]*(1.+c.jnoiseb[ncl]*minarea)
    print "ave 1 sigma J in 2arcsec ap= ",c.jzp[ncl] - 2.5*N.log10(errj) 


print "RELATIVE COMPLETENESS"

sfr=g2.sfrc
finalsf=g2.finalsf

lim1=1.98
lim2=2.49
print "fraction of CL1216 tot SFR below ",lim1," = ",N.sum(N.compress((finalsf > 0) & (sfr > lim1),sfr))/N.sum(N.compress(finalsf > 0,sfr))

print "fraction of CL1216 tot SFR below ",lim2," = ",N.sum(N.compress((finalsf > 0) & (sfr > lim2),sfr))/N.sum(N.compress(finalsf > 0,sfr))

sfr=g1.sfrc
finalsf=g1.finalsf

print "fraction of CL1054 tot SFR below ",lim1," = ",N.sum(N.compress((finalsf > 0) & (sfr > lim1),sfr))/N.sum(N.compress(finalsf > 0,sfr))

print "fraction of CL1054 tot SFR below ",lim2," = ",N.sum(N.compress((finalsf > 0) & (sfr > lim2),sfr))/N.sum(N.compress(finalsf > 0,sfr))

sfr=g0.sfrc
finalsf=g0.finalsf

print "fraction of CL1040 tot SFR below ",lim1," = ",N.sum(N.compress((finalsf > 0) & (sfr > lim1),sfr))/N.sum(N.compress(finalsf > 0,sfr))

print "fraction of CL1040 tot SFR below ",lim2," = ",N.sum(N.compress((finalsf > 0) & (sfr > lim2),sfr))/N.sum(N.compress(finalsf > 0,sfr))


psplotinit("rxj0152ewj.ps")
plotewj(g4.magj,g4.ew,g4.errew)
ppgplot.pgend()
psplotinit("rxj0152sfrj.ps")
plotsfrj(g4.magj,g4.sfr,g4.errsfr)
ppgplot.pgend()
psplotinit("rxj0152snsfr.ps")
plotsnsfr(g4.sn,g4.sfr,g4.errsfr)
ppgplot.pgend()
psplotinit("rxj0152histew.ps")
plothistew(g4.ew)
ppgplot.pgend()
psplotinit("rxj0152ratioj.ps")
plotratioj(g4.magj,g4.fj,g4.errfj,g4.fn,g4.errfn)
ppgplot.pgend()
g4.final=N.ones(len(g4.sfr),'f')
g4.finalsf=N.zeros(len(g4.sfr),'f')
#output=open("RXJ0152/sfrgt1",'w')
#output2=open("RXJ0152/sfrlt1",'w')
#for i in range(len(g4.sfr)):
#    if (g4.sn[i] > snmin) & (g4.ew[i] > ewmin):
#        g4.finalsf[i]=1.
#        output.write("%5.1f %5.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n" % (g4.x[i],g4.y[i],g4.ew[i],g4.errew[i],g4.sfr[i],g4.errsfr[i],g4.magj[i]))
#    if (g4.sn[i] > snmin) & (g4.ew[i] < ewmin):
#        g4.finalsf[i]=0.
#        output2.write("%5.1f %5.1f %6.1f %6.1f %6.1f %6.1f \n" % (g4.x[i],g4.y[i],g4.ew[i],g4.errew[i],g4.sfr[i],g4.errsfr[i]))
#output.close()
#output2.close()
psplotinit("rxj0152plotgal.ps")
plotgal(g4.x,g4.y,g4.final,g4.finalsf,5)
ppgplot.pgend()
psplotinit("rxj0152fitzp.ps")
plotfitzp(g4.x,g4.y,g4.f4j,g4.errf4j,4)
ppgplot.pgend()
print "SDSS comparison"
print "CL1040 total SFR = ",g0.totalsfrsdss, g0.totalsfrsdsserr
print "CL1054 total SFR = ",g1.totalsfrsdss, g1.totalsfrsdsserr
print "CL1216 total SFR = ",g2.totalsfrsdss, g2.totalsfrsdsserr

