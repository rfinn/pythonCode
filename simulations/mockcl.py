#!/usr/bin/env python
#"""
#measure contamination of velocity/radius member selection
#using mock cluster catalogs from Mo, Yang, van den Bosch
#"""
import sys
import sets
import copy
import os
#import Numeric as N
import numarray as N
import scipy
import pylab #use matplotlib instead of scipy pylab.std
import mystuff as my
import ppgplot
import time
starttime=time.clock()
print "start time = ",starttime

schdef=1.2
slwdef=4.
Mbsol=5.47 #absolute magnitude of sun in B (Cox 2000); filter from Bessell 1990
omega0=0.3
omegaL=0.7
h=.7
H0=h*100.#km/s/Mpc

G=4.31e-9#gravitational constant in (km/s)^2 Mpc Msolar^-1
cl=3.e5#speed of light in km/s
Mclmin=100.#min cluster mass in units of 10^12 Msun

bJmax=-16.
lbjmin=N.log10(10.**((Mbsol-bJmax)/2.5))
bjdefault=-18.
DH=cl/H0
#nr=1. #number of virial radii away from cluster 
#nv=3. #number of sigma away from cluster
#read in catalog of halos, keeping those w/M>=Mclmin

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
    dV=N.zeros(len(z),'f')
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


def plotdLdz():
    nv=3.
    nr=1.
    ppgplot.pgbeg("dLdz.ps/vcps",1,1)  #color port.
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
    xmax=1.01
    ymin=.95
    ymax=6.5

    ppgplot.pgsvp(x1,x4,y1,y4)  #sets viewport
    ppgplot.pgslw(slwdef)   #line width
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('bcnst',.2,2,'bcvnst',1.,5)  #tickmarks and labeling
    ppgplot.pgmtxt('b',2.5,0.5,0.5,"z")    #xlabel
    ppgplot.pgmtxt('l',2.6,0.5,0.5,"\gDL/\gDL(z=0)")

    z=N.arange(0.,1.2,.1)
    beta=((1+z)**2 - 1)/((1+z)**2 + 1)
    dV=N.zeros(len(z),'f')
    sigma=[400.,600.,800.,1000.]
    color=[1,2,3,4]
    for j in range(len(sigma)):
        dv=6*sigma[j]/3.e5#dv/c
        for i in range(len(z)):
            dz=dv/(1+z[i])*(1- ((1+z[i])**2 -1)/((1+z[i])**2+1))**(-2)
            z1=z[i]-0.5*dz
            z2=z[i]+0.5*dz
            dV[i]=my.dL(z2,h) - my.dL(z1,h)
            #dA=my.DA(z[i],h)*206264./1000.
            #dV[i]=DH*(1+z[i])*(dA)**2/(my.E(z[i]))/(1-beta[i])**2/DH**3
            #dV[i]=DH*(1+z[i])**2*(dA)**2/(my.E(z[i]))/DH**3#for comparison w/Hogg
            if z[i] < 1:
                print sigma[j],i,z[i],dV[i],dV[i]/dV[0]

        dV=dV/dV[0]
        ppgplot.pgsci(color[j])
        ppgplot.pgline(z,dV) 

    ppgplot.pgend()


def plotdLdz2():
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
    ppgplot.pgsvp(x1,x4,y1,y4)  #sets viewport

    schdef=1.2
    slwdef=4
    ppgplot.pgsch(schdef)
    xmin=0.
    xmax=1.1
    ymin=0.
    ymax=1.2
    z=N.arange(0.,5.,.1)
    beta=((1+z)**2 - 1)/((1+z)**2 + 1)
    dV=N.zeros(len(z),'f')
    for i in range(len(z)):
        #dV[i]=(1./DH)**3*(1+z[i])*(my.DA(z[i],h))**2/H0/(my.E(z[i]))/(1-beta[i])**2
        dA=my.DA(z[i],h)*206264./1000.
        dV[i]=DH*(1+z[i])*(dA)**2/(my.E(z[i]))/(1-beta[i])**2/DH**3
        #dV[i]=DH*(1+z[i])**2*(dA)**2/(my.E(z[i]))/DH**3#for comparison w/Hogg
    #nbin=5
    ppgplot.pgslw(slwdef)   #line width
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('bcnst',.2,2,'bcvnst',.2,2)  #tickmarks and labeling
    ppgplot.pgmtxt('b',2.5,0.5,0.5,"z")    #xlabel
    ppgplot.pgmtxt('l',2.6,0.5,0.5,"(1/DH)\u3\d c dV\dc\u/dv/d\gW")

    ppgplot.pgline(z,dV) 


    ppgplot.pgend()


def plotmagcuts():
    nv=3.
    nr=1.
    ppgplot.pgbeg("compcontmhalo-magcut.ps/vcps",1,1)  #color port.
    ppgplot.pgpap(8.,1.25)
    ppgplot.pgpage
    ppgplot.pgsch(1.2) #font size
    ppgplot.pgslw(3)   #line width

    # 1st panel with symbols w/ stddev errorbars


    x1=.1
    x2=.45
    x3=.6
    x4=.95
    y1=.15
    y2=.425
    y3=.575
    y4=.85
    xlabel=14.1-14.
    ylabel=1.15
    ppgplot.pgsvp(x1,x2,y3,y4)  #sets viewport
    bbJmax=magcut[0]
    g.cutonlbj(bbJmax)
    c.measurecontam(nv,nr,g)
    sub1plotmagcuts(c.mass,c.contam,c.compl)
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    label="M\db\dJ\u\u<"+str(bJmax)
    ppgplot.pgtext(xlabel,ylabel,label) 

    ppgplot.pgsvp(x3,x4,y3,y4)  #sets viewport
    #ppgplot.pgpanl(2,2)
    bbJmax=magcut[1]
    g.cutonlbj(bbJmax)
    c.measurecontam(nv,nr,g)
    sub1plotmagcuts(c.mass,c.contam,c.compl)
    label="M\db\dJ\u\u < "+str(bbJmax)
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    ppgplot.pgtext(xlabel,ylabel,label) 

    ppgplot.pgsvp(x1,x2,y1,y2)  #sets viewport
    #ppgplot.pgpanl(1,1)
    bbJmax=magcut[2]
    g.cutonlbj(bbJmax)
    c.measurecontam(nv,nr,g)
    sub1plotmagcuts(c.mass,c.contam,c.compl)
    label="M\db\dJ\u\u < "+str(bbJmax)
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    ppgplot.pgtext(xlabel,ylabel,label) 

    ppgplot.pgsvp(x3,x4,y1,y2)  #sets viewport
    #ppgplot.pgpanl(2,1)
    bbJmax=magcut[3]
    g.cutonlbj(bbJmax)
    c.measurecontam(nv,nr,g)
    sub1plotmagcuts(c.mass,c.contam,c.compl)
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    label="M\db\dJ\u\u < "+str(bbJmax)
    ppgplot.pgtext(xlabel,ylabel,label) 

    ppgplot.pgend()

def plotmagcutsgif():
    nv=3.
    nr=1.
    #bunch of pgplot initialization...
    ppgplot.pgbeg("compcontmhalo-magcut-gif.ps/vcps",1,1)  #color port.
    ppgplot.pgpap(8.,1.25)
    ppgplot.pgpage
    ppgplot.pgsch(1.2) #font size
    ppgplot.pgslw(3)   #line width

    #magcut=[-17.,-18.,-19.,-20.]
    magcut=[-18.,-19.,-20.]

    npx=3. #number of panels in x-direction
    npy=6. #number of panels in y-direction
    xpmin=0.1
    xpmax=.9
    ypmin=.1
    ypmax=.9
    dxp=(xpmax-xpmin)/npx
    dyp=(ypmax-ypmin)/npy
    xp=N.arange(xpmin,(xpmax+dxp),dxp)
    yp=N.arange(ypmin,(ypmax+dyp),dyp)
    x1=.1
    x2=.45
    x3=.6
    x4=.95
    y1=.15
    y2=.425
    y3=.575
    y4=.85
    xlabel=14.1-14.
    ylabel=1.15
    xl=[x1,x3,x1,x3]#left x value of viewport
    xr=[x2,x4,x2,x4]
    yb=[y3,y3,y1,y1]
    yt=[y4,y4,y2,y2]

    z=[0.00,0.2,0.4,0.6,0.8,1.05]
    zindex=[0,1,2,3,4,5]
    #zindex=[1]
    #cprefix=[cgif0,cgif02,cgif04,cgif06,cgif08,cgif1]
    #gprefix=[ggif0,ggif02,ggif04,ggif06,ggif08,ggif1]
    for j in zindex:
        #cc=cgif0#current cluster catalog
        #cg=ggif0#current galaxy catalog
        if j == 0:
            cc = cgif0
            cg = ggif0
        if j == 1:
            cc = cgif02
            cg = ggif02
        if j == 2:
            cc = cgif04
            cg = ggif04
        if j == 3:
            cc = cgif06
            cg = ggif06
        if j == 4:
            cc = cgif08
            cg = ggif08
        if j == 5:
            cc = cgif1
            cg = ggif1
        for i in range(len(magcut)):
            bbJmax=magcut[i] 
            ppgplot.pgsvp(xp[i],xp[i+1],yp[j],yp[j+1])  #sets viewport
            print "before magcut", len(cg.z),len(cg.x1),len(cgif02.x1)
            cg.cutonlbj(bbJmax)
            print i,j,"magcut = ",bbJmax, len(cg.z),len(cg.x1)
            cc.measurecontam(nv,nr,cg)
            flag=cc.contam + cc.compl #wacky clusters have either set to -999.
            #if i < 1:
            #    for j in range(20):
            #        print "dude ",cc.mass[j],cc.contam[j],cc.compl[j],flag[j]

            mass=N.compress(flag > 0,cc.mass)
            contam=N.compress(flag > 0,cc.contam)
            compl=N.compress(flag > 0,cc.compl)
            gsub1plotmagcuts(mass,contam,compl,i,j)
            ppgplot.pgsch(.8)
            ppgplot.pgslw(3)
            if j == 5:
                label="M\db\dJ\u\u<"+str(bbJmax)
                ppgplot.pgtext(xlabel,ylabel,label)
            if i == 0:
                label="z="+str(z[j])
                ppgplot.pgtext(-.3,.5,label)

    ppgplot.pgsch(1.2)
    ppgplot.pgslw(4)

    ppgplot.pgsvp(xp[0],xp[(len(xp)-1)],yp[0],yp[(len(yp)-1)])  #sets viewport
    ppgplot.pgmtxt('b',2.5,0.5,0.5,"log\d10\u(10\u14\d M\dhalo\u/h\u-1\d M\d\(2281)\u)")    #xlabel
    ppgplot.pgmtxt('l',2.5,0.5,0.5,"F\dt\u(F\di\u)")

    ppgplot.pgend()

def gsub1plotmagcuts(x,y1,y2,i,j):
    schdef=1.2
    slwdef=4
    ppgplot.pgsch(schdef)
    xmin=-.5
    xmax=1.1
    ymin=-.15
    ymax=1.3
    #nbin=5
    ppgplot.pgslw(slwdef)   #line width
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    if (i == 0) & (j == 0):
        ppgplot.pgbox('bcnlst',1.0,0,'bcnvst',1.,5)  #tickmarks and labeling
    if (i > 0) & (j == 0):
        ppgplot.pgbox('bclnst',1.0,0,'bcvst',1.,5)  #tickmarks and labeling
    if (i == 0) & (j > 0):
        ppgplot.pgbox('bclst',1.0,0,'bcnvst',1.,5)  #tickmarks and labeling
    if (i > 0) & (j > 0):
        ppgplot.pgbox('bclst',1.0,0,'bcvst',1.,5)  #tickmarks and labeling

    avecont=N.average(y1)
    avecomp=N.average(y2)
    (xbin,ybin,ybinerr)=my.biniterr(x,y1,nbin)
    xbin=N.log10(xbin)+12.-14.#add back 10^12 Msun, then divide by 10^14
    ppgplot.pgsch(1.5)
    ppgplot.pgpt(xbin,ybin,7)
    ppgplot.pgslw(3)
    ppgplot.pgerrb(6,xbin,ybin,ybinerr,2.)
    #my.errory(xbin,ybin,ybinerr)
    (xbin,ybin,ybinerr)=my.biniterr(x,y2,nbin)
    xbin=N.log10(xbin)+12.-14.#add back 10^12 Msun, then divide by 10^14
    ppgplot.pgsch(1.75)
    ppgplot.pgpt(xbin,ybin,17)
    ppgplot.pgsch(schdef)
    ppgplot.pgerrb(6,xbin,ybin,ybinerr,2.)
    ppgplot.pgslw(slwdef)
    x=N.arange((xmin-1),(xmax+1),.1)
    y=avecont*N.ones(len(x),'f')
    ppgplot.pgslw(2)
    ppgplot.pgsls(4)
    ppgplot.pgline(x,y)

    y=avecomp*N.ones(len(x),'f')
    ppgplot.pgsls(3)
    ppgplot.pgline(x,y)
    ppgplot.pgsls(1)
    ppgplot.pgslw(slwdef)
#my.errory(xbin,ybin,ybinerr)
    #ppgplot.pgsci(2)

def plotradcutsgif():
    nv=3.
    nr=1.
    #bunch of pgplot initialization...
    ppgplot.pgbeg("compcont-rad-z-gif.ps/vcps",1,1)  #color port.
    ppgplot.pgpap(8.,1.25)
    ppgplot.pgpage
    ppgplot.pgsch(1.2) #font size
    ppgplot.pgslw(3)   #line width

    #magcut=[-17.,-18.,-19.,-20.]
    magcut=[-18.]

    npx=2. #number of panels in x-direction
    npy=1. #number of panels in y-direction
    xpmin=0.1
    xpmax=.9
    ypmin=.1
    ypmax=.9
    dxp=(xpmax-xpmin)/npx
    dyp=(ypmax-ypmin)/npy
    xp=N.arange(xpmin,(xpmax+dxp),dxp)
    yp=N.arange(ypmin,(ypmax+dyp),dyp)

    xlabel=14.1-14.
    ylabel=1.15
    xv1=.1
    xv2=.45
    xv3=.6
    xv4=.95
    yv1=.35
    yv2=.6
    xlabel=0.05
    ylabel=.75
    dx=.18 #length of line
    ystep=.1 #step between labels
    dyl=.02 #offset b/w line and label in x-direction
    dxl=.015 #offset b/w line and label in y-direction

    #ppgplot.pgsvp(x1,x2,y1,y2)  #sets viewport

    z=N.array([0.00,0.2,0.42,0.62,0.82,1.05],'f')
    zindex=[0,1,2,3,4,5]
    complz=N.zeros(len(z),'f')
    errcomplz=N.zeros(len(z),'f')
    contamz=N.zeros(len(z),'f')
    errcontamz=N.zeros(len(z),'f')
    lcomplz=N.zeros(len(z),'f')
    lerrcomplz=N.zeros(len(z),'f')
    lcontamz=N.zeros(len(z),'f')
    lerrcontamz=N.zeros(len(z),'f')
    #zindex=[1]
    #cprefix=[cgif0,cgif02,cgif04,cgif06,cgif08,cgif1]
    #gprefix=[ggif0,ggif02,ggif04,ggif06,ggif08,ggif1]
    for j in zindex:
        #cc=cgif0#current cluster catalog
        #cg=ggif0#current galaxy catalog
        if j == 0:
            cc = cgif0
            cg = ggif0
        if j == 1:
            cc = cgif02
            cg = ggif02
        if j == 2:
            cc = cgif04
            cg = ggif04
        if j == 3:
            cc = cgif06
            cg = ggif06
        if j == 4:
            cc = cgif08
            cg = ggif08
        if j == 5:
            cc = cgif1
            cg = ggif1
        bbJmax=magcut[0] 

        #print "before magcut", len(cg.z),len(cg.x1),len(cgif02.x1)
        cg.cutonlbj(bbJmax)
        #print i,j,"magcut = ",bbJmax, len(cg.z),len(cg.x1)
        cc.measurecontam(nv,nr,cg)
        flag=cc.contam + cc.compl #wacky clusters have either set to -999.
        mass=N.compress(flag > 0,cc.mass)
        contam=N.compress(flag > 0,cc.contam)
        compl=N.compress(flag > 0,cc.compl)
        lcontam=N.compress(flag > 0,cc.lcontam)
        lcompl=N.compress((flag > 0) & (cc.lcompl > 0.),cc.lcompl)
        contamtype=N.compress(flag > 0,cc.contamtype)
        
        complz[j]=N.average(compl)
        #errcomplz[j]=scipy.stats.std(compl)
        errcomplz[j]=pylab.std(compl)
        contamz[j]=N.average(contam)
        #errcontamz[j]=scipy.stats.std(contam)
        errcontamz[j]=pylab.std(contam)
        try:
            lcomplz[j]=N.average(lcompl)
            lerrcomplz[j]=pylab.std(lcompl)
        except:
            lcomplz[j]=0.
            lerrcomplz[j]=0.
        lcontamz[j]=N.average(lcontam)
        lerrcontamz[j]=pylab.std(lcontam)

    for k in range(2):
        ppgplot.pgswin(-.1,1.25,-.1,1.3) #axes limits
        if k == 0:
            y1=contamz
            erry1=errcontamz
            y2=complz
            erry2=errcomplz
            ppgplot.pgsvp(xv1,xv2,yv1,yv2)  #sets viewport
            ppgplot.pgtext(0,1.1,"All")  #sets viewport
        if k == 1:
            y1=lcontamz
            erry1=lerrcontamz
            y2=lcomplz
            erry2=lerrcomplz
            ppgplot.pgsvp(xv3,xv4,yv1,yv2)  #sets viewport
            ppgplot.pgtext(0,1.1,"Late-type")  #sets viewport
        #ppgplot.pgsvp(xp[k],xp[k+1],yp[0],yp[len(yp)-1])  #sets viewport
        ppgplot.pgbox('bcnst',1.0,5,'bcnvst',1.,5)  #tickmarks and labeling
        ppgplot.pgsch(1.5)
        ppgplot.pgpt(z,y1,7)
        ppgplot.pgslw(3)
        ppgplot.pgerrb(6,z,y1,erry1,2.)
        ppgplot.pgsch(1.75)
        ppgplot.pgpt(z,y2,17)
        ppgplot.pgsch(schdef)
        ppgplot.pgerrb(6,z,y2,erry2,2.)
        ppgplot.pgsch(1.2)
        ppgplot.pgslw(4)
        ppgplot.pgmtxt('b',2.5,0.5,0.5,"z")    #xlabel
        ppgplot.pgmtxt('l',2.5,0.5,0.5,"F\dt\u(F\di\u)")
    for i in range(len(z)):
        print "z = %5.2f %5.2f %5.2f" %(z[i],contamz[i],lcontamz[i])
    ppgplot.pgend()

def plottotalsfrgif():
    nv=3.
    nr=1.
    sfroutput=open('totalsfr.dat','w')
    #bunch of pgplot initialization...
    ppgplot.pgbeg("totalsfr-mcl-z-gif.ps/vcps",1,1)  #color port.
    ppgplot.pgpap(8.,1.25)
    ppgplot.pgpage
    ppgplot.pgsch(1.2) #font size
    ppgplot.pgslw(3)   #line width

    #magcut=[-17.,-18.,-19.,-20.]
    magcut=[-19.]

    npx=2. #number of panels in x-direction
    npy=1. #number of panels in y-direction
    xpmin=0.1
    xpmax=.9
    ypmin=.1
    ypmax=.9
    dxp=(xpmax-xpmin)/npx
    dyp=(ypmax-ypmin)/npy
    xp=N.arange(xpmin,(xpmax+dxp),dxp)
    yp=N.arange(ypmin,(ypmax+dyp),dyp)

    xlabel=14.1-14.
    ylabel=1.15
    xv1=.1
    xv2=.45
    xv3=.6
    xv4=.95
    yv1=.35
    yv2=.6
    #yv1=.1
    #yv2=.45
    #yv3=.6
    #yv4=.95
    xlabel=0.05
    ylabel=.75
    dx=.18 #length of line
    ystep=.1 #step between labels
    dyl=.02 #offset b/w line and label in x-direction
    dxl=.015 #offset b/w line and label in y-direction

    ppgplot.pgsvp(xv1,xv4,yv1,yv2)  #sets viewport
    ppgplot.pgmtxt('b',2.5,0.5,0.5,"log\d10\u(10\u14\d M\dhalo\u/h\u-1\d M\d\(2281)\u)")    #xlabel

    ppgplot.pgmtxt('l',2.5,0.5,0.5,"Total SFR")
    ppgplot.pgsvp(xv1,xv2,yv1,yv2)  #sets viewport
    #ppgplot.pgbox('bcnst',1.0,5,'bcnvst',1.,5)  #tickmarks and labeling
    ppgplot.pgsch(1.)
    
    xmin=-.5
    xmax=1.
    ymin=0.
    ymax=50.
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('bcnlst',1.0,0,'bcvnst',10.,2)  #tickmarks and labeling
    xlabel=-.4
    ylabel=44.
    ppgplot.pgtext(xlabel,ylabel,'R\d200\u')
    ppgplot.pgsvp(xv3,xv4,yv1,yv2)  #sets viewport
    #ppgplot.pgbox('bcnst',1.0,5,'bcnvst',1.,5)  #tickmarks and labeling
    #ppgplot.pgsch(1.5)
    #ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('bcnlst',1.0,0,'bcvnst',10.,2)  #tickmarks and labeling
    #ppgplot.pgmtxt('b',2.5,0.5,0.5,"log\d10\u(10\u14\d M\dhalo\u/h\u-1\d M\d\(2281)\u)")    #xlabel
    ppgplot.pgtext(xlabel,ylabel,'0.5R\d200\u')

    #ppgplot.pgmtxt('l',2.5,0.5,0.5,"Total SFR(R<0.5R\d200\u)")
    z=N.array([0.00,0.2,0.42,0.62,0.82,1.05],'f')
    zindex=[0,1,2,3,4,5]
    
    for j in zindex:
        #cc=cgif0#current cluster catalog
        #cg=ggif0#current galaxy catalog
        if j == 0:
            cc = cgif0
            cg = ggif0
            ppgplot.pgsci(1)
        if j == 1:
            cc = cgif02
            cg = ggif02
            ppgplot.pgsci(2)
        if j == 2:
            cc = cgif04
            cg = ggif04
            ppgplot.pgsci(3)
        if j == 3:
            cc = cgif06
            cg = ggif06
            ppgplot.pgsci(4)
        if j == 4:
            cc = cgif08
            cg = ggif08
            ppgplot.pgsci(5)
        if j == 5:
            cc = cgif1
            cg = ggif1
            ppgplot.pgsci(6)
        #bbJmax=magcut[0] 
        #cg.cutonlbj(bbJmax)
        #print i,j,"magcut = ",bbJmax, len(cg.z),len(cg.x1)
        nv=3.
        nr=1.
        print "getting observed member ids for z = ",z[j]
        #cc.getobsmemberids(nv,nr,cg)
        cc.getmemberids(nv,nr,cg)
        #cc.measurecontam(nv,nr,cg)
        #flag=cc.contam + cc.compl #wacky clusters have either set to -999.
        #mass=N.compress(flag > 0,cc.mass)
        print "calculating total sfr for z = ",z[j]
        cc.calctotalsfr(cg)
        sfr1=cc.totalsfr
        nv=3.
        nr=.5
        #cc.getobsmemberids(nv,nr,cg)
        cc.getmemberids(nv,nr,cg)
        #cc.measurecontam(nv,nr,cg)
        #flag=cc.contam + cc.compl #wacky clusters have either set to -999.
        #mass=N.compress(flag > 0,cc.mass)
        cc.calctotalsfr(cg)
        sfr05=cc.totalsfr
        (xbin,ybin1,ybin1err)=my.biniterr(cc.mass,sfr1,nbin)
        (xbin,ybin05,ybin05err)=my.biniterr(cc.mass,sfr05,nbin)
        mass=N.log10(xbin)+12.-14.#add back 10^12 Msun, then divide by 10^14
        ppgplot.pgsvp(xv1,xv2,yv1,yv2)  #sets viewport
        ppgplot.pgpt(mass,ybin1,17)
        ppgplot.pgline(mass,ybin1)
        ppgplot.pgsvp(xv3,xv4,yv1,yv2)  #sets viewport
        ppgplot.pgpt(mass,ybin05,17)
        ppgplot.pgline(mass,ybin05)
        for l in range(len(mass)):
            sfroutput.write("%5.3f %5.3f %10.2f %10.2f \n" % (z[j],mass[l],sfr1[l],sfr05[l]))
            
    sfroutput.close()
    ppgplot.pgend()


def plottypemagcuts():
    nv=3.
    nr=1.

    ppgplot.pgbeg("typemhalo-magcut.ps/vcps",1,1)  #color port.
    ppgplot.pgpap(8.,1.25)
    ppgplot.pgpage
    ppgplot.pgsch(1.2) #font size
    ppgplot.pgslw(3)   #line width

    # 1st panel with symbols w/ stddev errorbars
    x1=.1
    x2=.45
    x3=.6
    x4=.95
    y1=.15
    y2=.425
    y3=.575
    y4=.85
    xlabel=14.1-14.
    ylabel=1.15
    ppgplot.pgsvp(x1,x2,y3,y4)  #sets viewport
    bbJmax=magcut[0]
    g.cutonlbj(bbJmax)
    c.measurecontam(nv,nr,g)
    sub1plottype(c.mass,c.contamtype)
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    label="M\db\dJ\u\u<"+str(bJmax)
    ppgplot.pgtext(xlabel,ylabel,label) 

    ppgplot.pgsvp(x3,x4,y3,y4)  #sets viewport
    #ppgplot.pgpanl(2,2)
    bbJmax=magcut[1]
    g.cutonlbj(bbJmax)
    c.measurecontam(nv,nr,g)
    sub1plottype(c.mass,c.contamtype)
    label="M\db\dJ\u\u < "+str(bbJmax)
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    ppgplot.pgtext(xlabel,ylabel,label) 

    ppgplot.pgsvp(x1,x2,y1,y2)  #sets viewport
    #ppgplot.pgpanl(1,1)
    bbJmax=magcut[2]
    g.cutonlbj(bbJmax)
    c.measurecontam(nv,nr,g)
    sub1plottype(c.mass,c.contamtype)
    label="M\db\dJ\u\u < "+str(bbJmax)
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    ppgplot.pgtext(xlabel,ylabel,label) 

    ppgplot.pgsvp(x3,x4,y1,y2)  #sets viewport
    #ppgplot.pgpanl(2,1)
    bbJmax=magcut[3]
    g.cutonlbj(bbJmax)
    c.measurecontam(nv,nr,g)
    sub1plottype(c.mass,c.contamtype)
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    label="M\db\dJ\u\u < "+str(bbJmax)
    ppgplot.pgtext(xlabel,ylabel,label) 

    ppgplot.pgend()

def plotsigcuts():
    nr=1.
    bJmax=-18.
    ppgplot.pgbeg("compcontmhalo-sigcut.ps/vcps",1,1)  #color port.
    ppgplot.pgpap(8.,1.25)
    ppgplot.pgpage
    ppgplot.pgsch(1.2) #font size
    ppgplot.pgslw(3)   #line width

    # 1st panel with symbols w/ stddev errorbars


    x1=.1
    x2=.45
    x3=.6
    x4=.95
    y1=.15
    y2=.425
    y3=.575
    y4=.85
    xlabel=14.3-14.
    ylabel=1.15
    ppgplot.pgsvp(x1,x2,y3,y4)  #sets viewport
    g.cutonlbj(bJmax)

    nv=6.
    c.measurecontam(nv,nr,g)
    print "dv < 6 sigma"
    print "mass    contam  compl"
    for i in range(len(c.mass)):
	print c.mass[i], c.contam[i], c.compl[i]
    sub1plotmagcuts(c.mass,c.contam,c.compl)
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    label="\gDv<"+str(nv)+"\gs"
    ppgplot.pgtext(xlabel,ylabel,label) 

    nv=3.
    ppgplot.pgsvp(x3,x4,y3,y4)  #sets viewport
    #ppgplot.pgpanl(2,2)
    c.measurecontam(nv,nr,g)
    print "dv < 3 sigma"
    print "mass    contam  compl"
    for i in range(len(c.mass)):
	print c.mass[i], c.contam[i], c.compl[i]

    sub1plotmagcuts(c.mass,c.contam,c.compl)
    label="\gDv<"+str(nv)+"\gs"
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    ppgplot.pgtext(xlabel,ylabel,label) 

    nv=2.
    ppgplot.pgsvp(x1,x2,y1,y2)  #sets viewport
    #ppgplot.pgpanl(1,1)
    c.measurecontam(nv,nr,g)
    print "dv < 3 sigma"
    print "mass    contam  compl"
    for i in range(len(c.mass)):
	print c.mass[i], c.contam[i], c.compl[i]

    sub1plotmagcuts(c.mass,c.contam,c.compl)
    label="\gDv<"+str(nv)+"\gs"
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    ppgplot.pgtext(xlabel,ylabel,label) 

    nv=1.
    ppgplot.pgsvp(x3,x4,y1,y2)  #sets viewport
    #ppgplot.pgpanl(2,1)
    c.measurecontam(nv,nr,g)
    sub1plotmagcuts(c.mass,c.contam,c.compl)
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    label="\gDv<"+str(nv)+"\gs"
    ppgplot.pgtext(xlabel,ylabel,label) 

    ppgplot.pgend()

def plotsigcutssdss():
    nr=1.
    bJmax=-18.
    ppgplot.pgbeg("compcontmhalo-sigcut-sdss.ps/vcps",1,1)  #color port.
    ppgplot.pgpap(8.,1.25)
    ppgplot.pgpage
    ppgplot.pgsch(1.2) #font size
    ppgplot.pgslw(3)   #line width

    # 1st panel with symbols w/ stddev errorbars


    x1=.1
    x2=.45
    x3=.6
    x4=.95
    y1=.15
    y2=.425
    y3=.575
    y4=.85
    xlabel=14.3-14.
    ylabel=1.15
    ppgplot.pgsvp(x1,x2,y3,y4)  #sets viewport
    g.cutonlbj(bJmax)

    nv=6.
    c.measurecontam(nv,nr,g)
    print "dv < 6 sigma"
    print "mass    contam  compl"
    for i in range(len(c.mass)):
	print c.mass[i], c.contam[i], c.compl[i]
    print "Average contamination within 6 sigma = ",N.average(c.contam)
    sub1plotmagcuts(c.mass,c.contam,c.compl)
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    label="\gDv<"+str(nv)+"\gs"
    ppgplot.pgtext(xlabel,ylabel,label) 

    nv=3.
    ppgplot.pgsvp(x3,x4,y3,y4)  #sets viewport
    #ppgplot.pgpanl(2,2)
    c.measurecontam(nv,nr,g)
    print "dv < 3 sigma"
    print "mass    contam  compl"
    for i in range(len(c.mass)):
	print c.mass[i], c.contam[i], c.compl[i]
    print "Average contamination within 3 sigma = ",N.average(c.contam)

    sub1plotmagcuts(c.mass,c.contam,c.compl)
    label="\gDv<"+str(nv)+"\gs"
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    ppgplot.pgtext(xlabel,ylabel,label) 

    nv=2.
    ppgplot.pgsvp(x1,x2,y1,y2)  #sets viewport
    #ppgplot.pgpanl(1,1)
    c.measurecontam(nv,nr,g)
    print "dv < 3 sigma"
    print "mass    contam  compl"
    for i in range(len(c.mass)):
	print c.mass[i], c.contam[i], c.compl[i]

    sub1plotmagcuts(c.mass,c.contam,c.compl)
    label="\gDv<"+str(nv)+"\gs"
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    ppgplot.pgtext(xlabel,ylabel,label) 

    nv=1.
    ppgplot.pgsvp(x3,x4,y1,y2)  #sets viewport
    #ppgplot.pgpanl(2,1)
    c.measurecontam(nv,nr,g)
    sub1plotmagcuts(c.mass,c.contam,c.compl)
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    label="\gDv<"+str(nv)+"\gs"
    ppgplot.pgtext(xlabel,ylabel,label) 

    ppgplot.pgend()

def plotngalmclradcuts():
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


def plotradcuts():
    nr=1.
    nv=3.
    bbJmax=-18.
    ppgplot.pgbeg("compcontmhalo-radcut.ps/vcps",1,1)  #color port.
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
    c.measurecontam(nv,nr,g)
    print "nr = ",nr, " ave contam = ",N.average(c.contam)
    sub1plotmagcuts(c.mass,c.contam,c.compl)
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    #label="R\dp\u < "+str(nr)+"R\dv\u"
    label=str1+str(nr)+str2
    ppgplot.pgtext(xlabel,ylabel,label) 

    nr=.75
    ppgplot.pgsvp(x3,x4,y3,y4)  #sets viewport
    #ppgplot.pgpanl(2,2)
    c.measurecontam(nv,nr,g)
    print "nr = ",nr, " ave contam = ",N.average(c.contam)
    sub1plotmagcuts(c.mass,c.contam,c.compl)
    label=str1+str(nr)+str2
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    ppgplot.pgtext(xlabel,ylabel,label) 

    nr=.5
    ppgplot.pgsvp(x1,x2,y1,y2)  #sets viewport
    #ppgplot.pgpanl(1,1)
    c.measurecontam(nv,nr,g)
    print "nr = ",nr, " ave contam = ",N.average(c.contam)
    sub1plotmagcuts(c.mass,c.contam,c.compl)
    label=str1+str(nr)+str2
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    ppgplot.pgtext(xlabel,ylabel,label) 

    nr=.25
    ppgplot.pgsvp(x3,x4,y1,y2)  #sets viewport
    #ppgplot.pgpanl(2,1)
    c.measurecontam(nv,nr,g)
    print "nr = ",nr, " ave contam = ",N.average(c.contam)
    sub1plotmagcuts(c.mass,c.contam,c.compl)
    ppgplot.pgsch(.8)
    ppgplot.pgslw(3)
    label=str1+str(nr)+str2
    ppgplot.pgtext(xlabel,ylabel,label) 

    ppgplot.pgend()

def sub1plotmagcuts(x,y1,y2):
    schdef=1.2
    slwdef=4
    ppgplot.pgsch(schdef)
    xmin=-.5
    xmax=1.
    ymin=-.15
    ymax=1.3
    #nbin=5
    ppgplot.pgslw(slwdef)   #line width
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('bcnlst',1.0,0,'bcvnst',1.,5)  #tickmarks and labeling
    ppgplot.pgmtxt('b',2.5,0.5,0.5,"log\d10\u(10\u14\d M\dhalo\u/h\u-1\d M\d\(2281)\u)")    #xlabel
    ppgplot.pgmtxt('l',2.1,0.5,0.5,"F\dt\u(F\di\u)")

    avecont=N.average(y1)
    avecomp=N.average(y2)
    (xbin,ybin,ybinerr)=my.biniterr(x,y1,nbin)
    xbin=N.log10(xbin)+12.-14.#add back 10^12 Msun, then divide by 10^14
    ppgplot.pgsch(1.5)
    ppgplot.pgpt(xbin,ybin,7)
    ppgplot.pgslw(3)
    ppgplot.pgerrb(6,xbin,ybin,ybinerr,2.)
    #my.errory(xbin,ybin,ybinerr)
    (xbin,ybin,ybinerr)=my.biniterr(x,y2,nbin)
    xbin=N.log10(xbin)+12.-14.#add back 10^12 Msun, then divide by 10^14
    ppgplot.pgsch(1.75)
    ppgplot.pgpt(xbin,ybin,17)
    ppgplot.pgsch(schdef)
    ppgplot.pgerrb(6,xbin,ybin,ybinerr,2.)
    ppgplot.pgslw(slwdef)
    x=N.arange((xmin-1),(xmax+1),.1)
    y=avecont*N.ones(len(x),'f')
    ppgplot.pgslw(2)
    ppgplot.pgsls(4)
    ppgplot.pgline(x,y)

    y=avecomp*N.ones(len(x),'f')
    ppgplot.pgsls(3)
    ppgplot.pgline(x,y)
    ppgplot.pgsls(1)
    ppgplot.pgslw(slwdef)
#my.errory(xbin,ybin,ybinerr)
    #ppgplot.pgsci(2)


def sub1plottype(x,y1):
    schdef=1.2
    slwdef=4
    ppgplot.pgsch(schdef)
    xmin=-.5
    xmax=.7
    ymin=-.15
    ymax=1.3
    #nbin=5
    x=N.compress(c.contam > 0,x)
    y1=N.compress(c.contam > 0,y1)
    y2=1-y1
    ppgplot.pgslw(slwdef)   #line width
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('bcnlst',1.0,0,'bcvnst',1.,5)  #tickmarks and labeling
    ppgplot.pgmtxt('b',2.5,0.5,0.5,"log\d10\u(10\u14\d M\dhalo\u/h\u-1\d M\d\(2281)\u)")    #xlabel
    ppgplot.pgmtxt('l',2.1,0.5,0.5,"Type Fraction")

    (xbin,ybin,ybinerr)=my.biniterr(x,y1,nbin)
    #print xbin
    #print ybin
    #print ybinerr
    xbin=N.log10(xbin)+12.-14.#add back 10^12 Msun, then divide by 10^14
    ppgplot.pgsch(1.5)
    ppgplot.pgpt(xbin,ybin,7)
    ppgplot.pgslw(3)
    ppgplot.pgerrb(6,xbin,ybin,ybinerr,2.)
    #my.errory(xbin,ybin,ybinerr)

    (xbin,ybin,ybinerr)=my.biniterr(x,y2,nbin)
    #print xbin
    #print ybin
    #print ybinerr

    xbin=N.log10(xbin)+12.-14.#add back 10^12 Msun, then divide by 10^14
    ppgplot.pgsch(1.75)
    ppgplot.pgpt(xbin,ybin,17)
    ppgplot.pgsch(schdef)
    ppgplot.pgerrb(6,xbin,ybin,ybinerr,2.)
    ppgplot.pgslw(slwdef)

    avecont=N.average(y1)
    avecomp=N.average(y2)

    x=N.arange((xmin-1),(xmax+1),.1)
    y=avecont*N.ones(len(x),'f')
    ppgplot.pgslw(2)
    ppgplot.pgsls(4)
    ppgplot.pgline(x,y)

    y=avecomp*N.ones(len(x),'f')
    ppgplot.pgsls(3)
    ppgplot.pgline(x,y)
    ppgplot.pgsls(1)
    ppgplot.pgslw(slwdef)
#my.errory(xbin,ybin,ybinerr)
    #ppgplot.pgsci(2)

def plotcummag():
    nv=3.
    nr=1.
    bjmin=-16.
    ppgplot.pgbeg("cumulativeMbj.ps/vcps",1,1)  #color port.
    ppgplot.pgpap(8.,1.25)
    ppgplot.pgpage
    ppgplot.pgsch(1.2) #font size
    ppgplot.pgslw(4)   #line width

    x1=.1
    x2=.45
    x3=.6
    x4=.95
    y1=.35
    y2=.6
    xlabel=-19
    ylabel=.5
    ppgplot.pgsvp(x1,x2,y1,y2)  #sets viewport

    g.cutonlbj(bjmin)
    c.measurecontam(nv,nr,g)

    #store all members in list
    obsmembid=[]
    membflag=[]
    gtype=[]
    for i in range(len(c.mass)):
        l1=list(c.obsmembids[i])
        l2=list(c.membflag[i])
        for j in range(len(l1)):
            obsmembid.append(l1[j])
            membflag.append(l2[j])
            gtype.append(g.type[int(l1[j])])
    obsmembid=N.array(obsmembid,'f')
    membflag=N.array(membflag,'f')
    gtype=N.array(gtype,'f')
    membid=N.compress((membflag > 0.),obsmembid)
    contid=N.compress((membflag < 1.),obsmembid)
    latemembid=N.compress((gtype > 0) & (membflag > 0.),obsmembid)
    latecontid=N.compress((gtype > 0) & (membflag < 1.),obsmembid)
    print "CUMULATIVE COMPARING MAG OF ALL GALAXIES"
    cumsub1(membid,contid)
    cumsub2(xlabel,ylabel,"All")
    ppgplot.pgsvp(x3,x4,y1,y2)  #sets viewport
    print "CUMULATIVE COMPARING MAG OF LATE-TYPE memb/contam GALAXIES"
    cumsub1(latemembid,latecontid)
    cumsub2(xlabel,ylabel,"Late-type")

    ppgplot.pgend()

def cumsub2(xlabel,ylabel,label):#draw key
    schdef=ppgplot.pgqch()
    ppgplot.pgsch(1.1)
    ppgplot.pgtext(xlabel,ylabel,label)
    ppgplot.pgsch(1.1)
    xs=N.arange(-20.,-18.,.1)
    ys=.2*N.ones(len(xs),'f')
    ppgplot.pgsls(1)
    ppgplot.pgline(xs,ys)
    ppgplot.pgtext(-20.1,.18,"Memb")
    ys=.1*N.ones(len(xs),'f')
    ppgplot.pgsls(4)
    ppgplot.pgline(xs,ys)
    ppgplot.pgtext(-20.1,.08,"Contam")
    ppgplot.pgsch(schdef)
    ppgplot.pgsls(1)
def cumsub1(membid,contid):

    xmax=-22.5
    xmin=-15.5
    ymax=1.1
    ymin=-.1
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('bcnst',2.0,2,'bcvnst',1.,5)  #tickmarks and labeling
    ppgplot.pgmtxt('b',2.5,0.5,0.5,"M\db\dJ\u\u")    #xlabel
    ppgplot.pgmtxt('l',2.1,0.5,0.5,"Cumulative Fraction")



    mag=N.zeros(len(membid),'f')
    #print "number of members = ",len(membid)
    for i in range(len(membid)):
        mag[i]=Mbsol - 2.5*g.lbj[int(membid[i])]
        #if i < 25:
        #    print i,membid[i],g.lbj[int(membid[i])],mag[i]
    mag=-1*mag
    (x,y)=my.cumulative(mag)
    x=-1*x
    k1=mag
    #y=1-y
    #print mag[i]
    #print mag[i]
    #xmax=max(x)
    #xmin=min(x)
    #ymax=max(y)
    #ymin=min(y)
    #ppgplot.pglab("M\db\dJ\u\u","Cumulative Fraction","")
    ppgplot.pgline(x,y)
    ppgplot.pgsls(4)
    mag=N.zeros(len(contid),'f')
    for i in range(len(contid)):
        mag[i]=Mbsol - 2.5*g.lbj[int(contid[i])]
    mag=-1.*mag
    (x,y)=my.cumulative(mag)
    x=-1.*x
    k2=mag
    ppgplot.pgline(x,y)
    ppgplot.pgsci(1)
    ppgplot.pgsls(1)
    output=open('data1.dat','w')
    for i in k1:
        i=str(i)
        #print i
        i=i+"\n"
        output.write(i)
    output.close()
    output=open('data2.dat','w')
    for i in k2:
        i=str(i)+"\n"
        output.write(str(i))
    output.close()
    os.system("ks \n")

def plotcumdv():#plot cumulative distribution of memb/contam versus dv/sigma
    nv=3.
    nr=2.
    bjmin=-17.
    ppgplot.pgbeg("cumulativedv.ps/vcps",1,1)  #color port.
    ppgplot.pgpap(8.,1.25)
    ppgplot.pgpage
    ppgplot.pgsch(1.2) #font size
    ppgplot.pgslw(4)   #line width

    x1=.1
    x2=.45
    x3=.6
    x4=.95
    y1=.35
    y2=.6
    xlabel=1.2
    ylabel=.4
    dx=.6 #length of line
    ystep=.1 #step between labels
    dyl=.02 #offset b/w line and label in x-direction
    dxl=.1 #offset b/w line and label in y-direction
    
    ppgplot.pgsvp(x1,x2,y1,y2)  #sets viewport

    g.cutonlbj(bjmin)
    c.measurecontam(nv,nr,g)

    #store all members in list
    obsmembid=[]
    obsmembdv=[]
    membflag=[]
    gtype=[]
    for i in range(len(c.mass)):
        l1=list(c.obsmembids[i])
        l2=list(c.membflag[i])
        l3=list(c.obsmembdv[i])
        for j in range(len(l1)):
            gindex=int(l1[j])
            obsmembid.append(l1[j])
            membflag.append(l2[j])
            gtype.append(g.type[int(l1[j])])
            obsmembdv.append(l3[j])
    obsmembid=N.array(obsmembid,'f')
    obsmembdv=N.array(obsmembdv,'f')
    membflag=N.array(membflag,'f')
    gtype=N.array(gtype,'f')
    membid=N.compress((membflag > 0.),obsmembid)
    membdv=N.compress((membflag > 0.),obsmembdv)
    contid=N.compress((membflag < 1.),obsmembid)
    contdv=N.compress((membflag < 1.),obsmembdv)
    latemembid=N.compress((gtype > 0) & (membflag > 0.),obsmembid)
    latemembdv=N.compress((gtype > 0) & (membflag > 0.),obsmembdv)
    latecontid=N.compress((gtype > 0) & (membflag < 1.),obsmembid)
    latecontdv=N.compress((gtype > 0) & (membflag < 1.),obsmembdv)
    print "CUMULATIVE COMPARING dv OF ALL GALAXIES"
    cumdvsub1(membid,membdv,contid,contdv)
    linelabel(xlabel,ylabel,dx,ystep,dxl,dyl,"All")
    ppgplot.pgsvp(x3,x4,y1,y2)  #sets viewport
    print "CUMULATIVE COMPARING dv OF LATE-TYPE memb/contam GALAXIES"
    cumdvsub1(latemembid,latemembdv,latecontid,latecontdv)
    #cumdvsub2(xlabel,ylabel,"Late-type")
    linelabel(xlabel,ylabel,dx,ystep,dxl,dyl,"Late")

    ppgplot.pgend()

def cumdvsub2(xlabel,ylabel,label):#draw key
    schdef=ppgplot.pgqch()
    ppgplot.pgsch(1.1)
    ppgplot.pgtext(xlabel,ylabel,label)
    ppgplot.pgsch(1.1)
    xs=N.arange(1.4,2.,.1)
    ys=.2*N.ones(len(xs),'f')
    ppgplot.pgsls(1)
    ppgplot.pgline(xs,ys)
    ppgplot.pgtext(2.1,.18,"Memb")
    ys=.1*N.ones(len(xs),'f')
    ppgplot.pgsls(4)
    ppgplot.pgline(xs,ys)
    ppgplot.pgtext(2.1,.08,"Contam")
    ppgplot.pgsch(schdef)
    ppgplot.pgsls(1)

def cumdvsub1(membid,membdv,contid,contdv):

    xmax=3.2
    xmin=-0.2
    ymax=1.1
    ymin=-.1
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('bcnst',1.0,2,'bcvnst',1.,5)  #tickmarks and labeling
    ppgplot.pgmtxt('b',2.5,0.5,0.5,"\gDv/\gs")    #xlabel
    ppgplot.pgmtxt('l',2.1,0.5,0.5,"Cumulative Fraction")


    dv=membdv
    dv=abs(dv)
    print "average dv = ",N.average(dv),min(dv),max(dv)
    (x,y)=my.cumulative(dv)
    k1=dv
    #y=1-y
    #print mag[i]
    #print mag[i]
    #xmax=max(x)
    #xmin=min(x)
    #ymax=max(y)
    #ymin=min(y)
    #ppgplot.pglab("M\db\dJ\u\u","Cumulative Fraction","")
    ppgplot.pgline(x,y)
    ppgplot.pgsls(4)
    dv=contdv
    dv=abs(dv)
    print "average dv = ",N.average(dv),min(dv),max(dv)
    (x,y)=my.cumulative(dv)
    k2=dv
    ppgplot.pgline(x,y)
    ppgplot.pgsci(1)
    ppgplot.pgsls(1)
    output=open('data1.dat','w')
    for i in k1:
        i=str(i)
        #print i
        i=i+"\n"
        output.write(i)
    output.close()
    output=open('data2.dat','w')
    for i in k2:
        i=str(i)+"\n"
        output.write(str(i))
    output.close()
    os.system("ks \n")

def plothistdv():#plot histogram of memb/contam versus dv/sigma
    nv=6.
    nr=2.
    bjmin=bjdefault
    ppgplot.pgbeg("histdv.ps/vcps",1,1)  #color port.
    ppgplot.pgpap(8.,1.25)
    ppgplot.pgpage
    ppgplot.pgsch(1.2) #font size
    ppgplot.pgslw(4)   #line width

    x1=.1
    x2=.45
    x3=.6
    x4=.95
    y1=.35
    y2=.6
    xlabel=0.15
    ylabel=.2
    dx=.5 #length of line
    ystep=.1 #step between labels
    dyl=.02 #offset b/w line and label in x-direction
    dxl=.1 #offset b/w line and label in y-direction

    ppgplot.pgsvp(x1,x2,y1,y2)  #sets viewport


    g.cutonlbj(bjmin)
    c.getallmembers(g,nv,nr)
    c.measurecontam(nv,nr,g)

    #store all members in list
    obsmembid=[]
    obsmembdv=[]
    membflag=[]
    gtype=[]
    for i in range(len(c.mass)):
        l1=list(c.obsmembids[i])
        l2=list(c.membflag[i])
        l3=list(c.obsmembdv[i])
        for j in range(len(l1)):
            gindex=int(l1[j])
            obsmembid.append(l1[j])
            membflag.append(l2[j])
            gtype.append(g.type[int(l1[j])])
            obsmembdv.append(l3[j])
    obsmembid=N.array(obsmembid,'f')
    obsmembdv=N.array(obsmembdv,'f')
    membflag=N.array(membflag,'f')
    gtype=N.array(gtype,'f')
    membid=N.compress((membflag > 0.),obsmembid)
    membdv=N.compress((membflag > 0.),obsmembdv)
    contid=N.compress((membflag < 1.),obsmembid)
    contdv=N.compress((membflag < 1.),obsmembdv)
    latemembid=N.compress((gtype > 0) & (membflag > 0.),obsmembid)
    latemembdv=N.compress((gtype > 0) & (membflag > 0.),obsmembdv)
    latecontid=N.compress((gtype > 0) & (membflag < 1.),obsmembid)
    latecontdv=N.compress((gtype > 0) & (membflag < 1.),obsmembdv)
    #print membdv
    histdvsub1(membid,membdv,contid,contdv)
    ppgplot.pgtext(xlabel,.95,"All")
    linelabel(xlabel,ylabel,dx,ystep,dxl,dyl,"")
    ppgplot.pgsvp(x3,x4,y1,y2)  #sets viewport
    histdvsub1(latemembid,latemembdv,latecontid,latecontdv)
    ppgplot.pgtext(xlabel,.95,"Late")
    linelabel(xlabel,ylabel,dx,ystep,dxl,dyl,"")

    ppgplot.pgend()

def linelabel(xlabel,ylabel,dx,ystep,dxl,dyl,label):#draw key
    schdef=ppgplot.pgqch()
    ppgplot.pgsch(1.1)
    ppgplot.pgtext(xlabel,ylabel,label)
    ppgplot.pgsch(1.1)
    ylabel=ylabel-2.*ystep
    xs=N.arange(xlabel,(xlabel+dx),.01)
    ys=ylabel*N.ones(len(xs),'f')
    ppgplot.pgsls(1)
    ppgplot.pgline(xs,ys)
    ppgplot.pgtext((xlabel+dxl+dx),(ylabel-dyl),"Memb")
    ylabel=ylabel-ystep
    ys=ylabel*N.ones(len(xs),'f')
    ppgplot.pgsls(4)
    ppgplot.pgline(xs,ys)
    ppgplot.pgtext((xlabel+dxl+dx),(ylabel-dyl),"Contam")
    ppgplot.pgsch(schdef)
    ppgplot.pgsls(1)

def histdvsub1(membid,membdv,contid,contdv):
    dvall=list(membdv)+list(contdv)
    dvall=abs(N.array(dvall,'f'))
    print "checking dvall ",len(dvall),len(membdv),len(contdv)
    xmax=3.2
    xmin=-.2
    ymax=1.1
    ymin=-.2
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('bcnst',1.0,2,'bcvnst',1.,5)  #tickmarks and labeling
    ppgplot.pgmtxt('b',2.5,0.5,0.5,"\gDv/\gs")    #xlabel
    ppgplot.pgmtxt('l',2.1,0.5,0.5,"Fraction")
    
    xmax=3.
    xmin=0.
    dbin=0.25
    bins=N.arange(xmin,xmax,dbin)#make array of bins in step of 0.5mag
    nbin=len(bins)
    dv=abs(membdv)
    #print dv
    print "average dv = ",N.average(dv),min(dv),max(dv)
    y=scipy.stats.histogram2(dv,bins)
    yall=scipy.stats.histogram2(dvall,bins)
    ybin=N.zeros(len(bins),'f')
    ybinerr=N.zeros(len(bins),'f')
    #y=y/float(len(dv))
    for i in range(len(bins)):
        try:
            #ybin[i]=float(y[i])/float(yall[i])
	    (ybin[i],ybinerr[i])=my.ratioerror(float(y[i]),float(yall[i]))
        except ZeroDivisionError:
            ybin[i]=0.
    #erry=N.sqrt(y)/float(len(dv))
    #ppgplot.pgerrb(6,bins,y,erry,2.)
    #ppgplot.pgpt(bins,y,3)
    my.drawhist(bins,ybin)
    my.errory(bins+(bins[1]-bins[0])/2.,ybin,ybinerr)
    ppgplot.pgsls(4)
    ppgplot.pgslw(4)
    #print "average dv = ",N.average(dv),min(dv),max(dv)
    dv=abs(contdv)
    y=scipy.stats.histogram2(dv,bins)
    ybin=N.zeros(len(bins),'f')
    ybinerr=N.zeros(len(bins),'f')
    #y=y/float(len(dv))
    for i in range(len(bins)):
        try:
            #ybin[i]=float(y[i])/float(yall[i])
	    (ybin[i],ybinerr[i])=my.ratioerror(float(y[i]),float(yall[i]))
        except ZeroDivisionError:
            ybin[i]=0.

    #y=y/yall
    #erry=N.sqrt(y)/float(len(dv))
    #ppgplot.pgpt(bins,y,3)
    #ppgplot.pgsci(2)
    #ppgplot.pgerrb(6,bins,y,erry,2.)
    my.drawhist(bins,ybin)
    my.errory(bins+(bins[1]-bins[0])/2.,ybin,ybinerr)
    ppgplot.pgsci(1)
    ppgplot.pgsls(1)

def plothistdr():#plot histogram of memb/contam versus dv/sigma
    nv=3.
    nr=1.
    bjmin=bjdefault
    c.getallmembers(g,nv,nr)
    ppgplot.pgbeg("histdr.ps/vcps",1,1)  #color port.
    ppgplot.pgpap(8.,1.25)
    ppgplot.pgpage
    ppgplot.pgsch(1.2) #font size
    ppgplot.pgslw(4)   #line width

    x1=.1
    x2=.45
    x3=.6
    x4=.95
    y1=.35
    y2=.6
    xlabel=0.05
    ylabel=.75
    dx=.18 #length of line
    ystep=.1 #step between labels
    dyl=.02 #offset b/w line and label in x-direction
    dxl=.015 #offset b/w line and label in y-direction

    ppgplot.pgsvp(x1,x2,y1,y2)  #sets viewport

    g.cutonlbj(bjmin)
    c.getallmembers(g,nv,nr)
    c.measurecontam(nv,nr,g)

    #store all members in list
    obsmembid=[]
    obsmembdr=[]
    membflag=[]
    gtype=[]
    for i in range(len(c.mass)):
        l1=list(c.obsmembids[i])
        l2=list(c.membflag[i])
        l3=list(c.obsmembdr[i])
        for j in range(len(l1)):
            gindex=int(l1[j])
            obsmembid.append(l1[j])
            membflag.append(l2[j])
            gtype.append(g.type[int(l1[j])])
            obsmembdr.append(l3[j])
    obsmembid=N.array(obsmembid,'f')
    obsmembdr=N.array(obsmembdr,'f')
    membflag=N.array(membflag,'f')
    gtype=N.array(gtype,'f')
    membid=N.compress((membflag > 0.),obsmembid)
    membdr=N.compress((membflag > 0.),obsmembdr)
    contid=N.compress((membflag < 1.),obsmembid)
    contdr=N.compress((membflag < 1.),obsmembdr)
    latemembid=N.compress((gtype > 0) & (membflag > 0.),obsmembid)
    latemembdr=N.compress((gtype > 0) & (membflag > 0.),obsmembdr)
    latecontid=N.compress((gtype > 0) & (membflag < 1.),obsmembid)
    latecontdr=N.compress((gtype > 0) & (membflag < 1.),obsmembdr)
    #print membdv
    histdrsub1(membid,membdr,contid,contdr)
    linelabel(xlabel,ylabel,dx,ystep,dxl,dyl,"All")
    ppgplot.pgsvp(x3,x4,y1,y2)  #sets viewport
    histdrsub1(latemembid,latemembdr,latecontid,latecontdr)
    linelabel(xlabel,ylabel,dx,ystep,dxl,dyl,"Late")

    ppgplot.pgend()

def histdrsub1(membid,membdv,contid,contdv):
    dvall=list(membdv)+list(contdv)
    xmax=1.1
    xmin=-.1
    ymax=1.1
    ymin=-.1
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('bcnst',1.0,5,'bcvnst',1.,5)  #tickmarks and labeling
    ppgplot.pgmtxt('b',2.5,0.5,0.5,"R\dp\u/R\dv\u")    #xlabel
    ppgplot.pgmtxt('l',2.1,0.5,0.5,"Fraction")
    
    xmax=1.
    xmin=0.
    bins=N.arange(xmin,xmax,0.1)#make array of bins in step of 0.5mag
    nbin=len(bins)
    dv=membdv
    print "average dv = ",N.average(dv),min(dv),max(dv)
    y=scipy.stats.histogram2(dv,bins)
    yall=scipy.stats.histogram2(dvall,bins)
    ybin=N.zeros(len(bins),'f')
    ybinerr=N.zeros(len(bins),'f')
    for i in range(len(bins)):
        try:
	    (ybin[i],ybinerr[i])=my.ratioerror(float(y[i]),float(yall[i]))
            #ybin[i]=float(y[i])/float(yall[i])
	    
        except ZeroDivisionError:
            ybin[i]=0.
    my.drawhist(bins,ybin)
    my.errory(bins+(bins[1]-bins[0])/2.,ybin,ybinerr)
    ppgplot.pgsls(4)
    dv=contdv
    y=scipy.stats.histogram2(dv,bins)
    ybin=N.zeros(len(bins),'f')
    ybinerr=N.zeros(len(bins),'f')
    for i in range(len(bins)):
        try:
	    (ybin[i],ybinerr[i])=my.ratioerror(float(y[i]),float(yall[i]))
            #ybin[i]=float(y[i])/float(yall[i])
        except ZeroDivisionError:
            ybin[i]=0.
    ppgplot.pgslw(4)
    my.drawhist(bins,ybin)
    my.errory(bins+(bins[1]-bins[0])/2.,ybin,ybinerr)
    ppgplot.pgsci(1)
    ppgplot.pgsls(1)


def plotsigcutsdr():#plot contamination vs dr for various sigma cuts
    nr=1.
    bJmax=-16.
    ppgplot.pgbeg("contdr-sigcut.ps/vcps",1,1)  #color port.
    ppgplot.pgpap(8.,1.25)
    ppgplot.pgpage
    ppgplot.pgsch(1.2) #font size
    ppgplot.pgslw(3)   #line width

    # 1st panel with symbols w/ stddev errorbars


    x1=.1
    x2=.45
    x3=.6
    x4=.95
    y1=.15
    y2=.425
    y3=.575
    y4=.85
    xp1=[x1,x3,x1,x3]
    xp2=[x2,x4,x2,x4]
    yp1=[y3,y3,y1,y1]
    yp2=[y4,y4,y2,y2]
    xlabel=2.
    ylabel=1.15

    g.cutonlbj(bJmax)

    dv=[6.,3.,2.,1.]
    dr=N.arange(0,3.5,.5)
    for j in range(len(dv)):
        nv=dv[j]
        ppgplot.pgsvp(xp1[j],xp2[j],yp1[j],yp2[j])  #sets viewport    
        lcontam=N.zeros(len(dr)-1,'f')
        lcompl=N.zeros(len(dr)-1,'f')
        econtam=N.zeros(len(dr)-1,'f')
        ecompl=N.zeros(len(dr)-1,'f')
        errlcontam=N.zeros(len(dr)-1,'f')
        errlcompl=N.zeros(len(dr)-1,'f')
        errecontam=N.zeros(len(dr)-1,'f')
        errecompl=N.zeros(len(dr)-1,'f')
        contam=N.zeros(len(dr)-1,'f')
        compl=N.zeros(len(dr)-1,'f')
        errcontam=N.zeros(len(dr)-1,'f')
        errcompl=N.zeros(len(dr)-1,'f')
        rave=N.zeros((len(dr)-1),'f')
        for i in range(len(dr)-1):
            nr1=dr[i]
            nr2=dr[i+1]
            rave[i]=0.5*(nr1+nr2)
            c.getobsmemberidsdr(nv,nr1,nr2)
            c.contam=N.array(c.contam,'f')
            c.lcontam=N.array(c.lcontam,'f')
            c.econtam=N.array(c.econtam,'f')
            c.compl =N.array(c.compl,'f')
            c.lcompl=N.array(c.lcompl,'f')
            c.ecompl=N.array(c.ecompl,'f')
            c.contam=N.compress(c.contam > - 1,c.contam)
            c.lcontam=N.compress(c.lcontam > -1,c.lcontam)
            c.econtam=N.compress(c.econtam > -1,c.econtam)
            c.compl=N.compress(c.compl > -1,c.compl)
            c.lcompl=N.compress(c.lcompl > -1,c.lcompl)
            c.ecompl=N.compress(c.ecompl > -1,c.ecompl)
            
            lcontam[i]=N.average(c.lcontam)
            errlcontam[i]=pylab.std(c.lcontam)
            lcompl[i]=N.average(c.lcompl)
            errlcompl[i]=pylab.std(c.lcompl)
            econtam[i]=N.average(c.econtam)
            errecontam[i]=pylab.std(c.econtam)
            ecompl[i]=N.average(c.ecompl)
            errecompl[i]=pylab.std(c.ecompl)
            contam[i]=N.average(c.contam)
            errcontam[i]=pylab.std(c.contam)
            compl[i]=N.average(c.compl)
            errcompl[i]=pylab.std(c.compl)

        #sub1plotdr(dr,lcontam,errlcontam,lcompl,errlcompl,econtam,errecontam,ecompl,errecompl)
        sub1plotdr(rave,contam,errcontam,compl,errcompl,econtam,errecontam,ecompl,errecompl)
        ppgplot.pgsch(.8)
        ppgplot.pgslw(3)
        label="\gDv < "+str(nv)+"\gs"
        ppgplot.pgtext(xlabel,ylabel,label) 

        
    ppgplot.pgend()

def sub1plotdr(x,y1,erry1,y2,erry2,y3,erry3,y4,erry4):
    schdef=1.2
    slwdef=4
    ppgplot.pgsch(schdef)
    xmin=-.1
    xmax=3.1
    ymin=-.15
    ymax=1.3
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('bcnst',1.0,0,'bcvnst',1.,5)  #tickmarks and labeling
    ppgplot.pgmtxt('b',2.5,0.5,0.5,"R\dp\u (Mpc)")    #xlabel
    ppgplot.pgmtxt('l',2.1,0.5,0.5,"Ft(Fi)")


    ppgplot.pgpt(x,y1,12)
    ppgplot.pgerrb(6,x,y1,erry1,2.)
    ppgplot.pgpt(x,y2,18)
    ppgplot.pgerrb(6,x,y2,erry2,2.)
    #ppgplot.pgpt(x,y3,17)
    #ppgplot.pgerrb(6,x,y3,erry3,2.)
    #ppgplot.pgpt(x,y1,21)
    #ppgplot.pgerrb(6,x,y4,erry4,2.)

def plotsffdr():
    sfrmin=1.
    my.psplotinit("sffdr200.ps")
    xmin=0.
    xmax=5.5
    ymin=-.1
    ymax=1.1
    nbin=5
    ppgplot.pgbox("",0.0,0,"",0.0,1)
    ppgplot.pgenv(0,xmax,ymin,ymax,0)
    ppgplot.pglab("dr/R\d200\u","SF Fraction","")
    x=N.arange((xmin-1.),(xmax+1.),1.,'f')
    y=0.352*N.ones(len(x),'f')
    ppgplot.pgsls(1)
    ppgplot.pgline(x,y)
    ppgplot.pgsls(4)
    y2=y-.005
    ppgplot.pgline(x,y2)
    y2=y+.005
    ppgplot.pgline(x,y2)
    ppgplot.pgsls(1)

    z=[0.,.2,.4,.6,.8,1.]
    zcolor=[1,2,3,4,5,6]
    for j in range(len(z)):
        if j == 0:
            c=cgif0
            g=ggif0
        if j == 1:
            c=cgif02
            g=ggif02
        if j == 2:
            c=cgif04
            g=ggif04
        if j == 3:
            c=cgif06
            g=ggif06
        if j == 4:
            c=cgif08
            g=ggif08
        if j == 5:
            c=cgif1
            g=ggif1
        ppgplot.pgsci(zcolor[j])
        sffdrsub(c,g)

    ppgplot.pgend()

def sffdrsub(c,g):
    sfrmin=0.
    dr=[]
    sfr=[]
    ew=[]
    final=[]
    finalsf=[]
    for j in range(len(c.mass)):
        obsmembids=list(c.obsmembids[j])
        #print c.obsmembids[j],len(c.obsmembids[j])
        z=len(c.obsmembids[j])
        #print z
        if z <= 1:
            continue
        #print obsmembids, obsmembids[0]
        for i in obsmembids:
            i=int(i)
            print 
            dr.append(g.drRvir[i])
            sfr.append(g.sfr[i])
            ew.append(g.ew[i])
            final.append(1.)
            if (g.sfr[i] > sfrmin):
                finalsf.append(1.)
            else:
                finalsf.append(0.)
    ytot=N.array(final,'f')
    ysf=N.array(finalsf,'f')
    dr=N.array(dr,'f')
    #sig=N.array(sig,'f')
    #print "ytot, ysf, dr"
    print ytot
    print ysf
    print dr
    (xbin,y1sum)=my.binitsumequal(dr,ysf,nbin)#get total sfr in each mag bin
    (xbin,y2sum)=my.binitsumequal(dr,ytot,nbin)#get total sfr in each mag bin
    (sffrac,err) = my.ratioerror(y1sum,y2sum)
    a1=N.sum(y1sum)
    b1=N.sum(y2sum)
    (a,b)=my.ratioerror(a1,b1)
    print "average fraction of sf galaxies (poisson error)= ",a,b
    a=N.average(sffrac)
    b=pylab.std(sffrac)
    print "average fraction of sf galaxies (std)= ",a,b
    sffrac = N.array(sffrac,'f')
    ppgplot.pgpt(xbin,sffrac,17)
    my.errory(xbin,sffrac,err)


class Cluster:
    def __init__(self):
        self.id = []
        self.x1 = []
        self.x2  = []
        self.x3  = []
        self.rvir = []
        self.mass = []
    def readfiles(self,halocat):
        self.v1 = []
        self.v2  = []
        self.v3  = []
        self.sigma = []

        infile=open(halocat,'r')
        i=0
        for line in infile:
            if line.find('#') > -1:
                continue
            fields=line.split()
            if len(fields) < 4:#skip line with 3 entries
                continue
            if (float(fields[0]) < Mclmin): #
                print "HEYYYYYYYY, number of clusters = ",len(self.mass)
                infile.close()
                break
            if float(fields[0]) > Mclmin:
                #print "heyah",fields[0],Mclmin
                self.id.append(float(i+1))
                self.mass.append(float(fields[0]))#10^12 Msun
                self.x1.append(float(fields[1]))
                self.x2.append(float(fields[2]))
                self.x3.append(float(fields[3]))
                self.v1.append(float(fields[4]))
                self.v2.append(float(fields[5]))
                self.v3.append(float(fields[6]))
                self.sigma.append(float(fields[7]))
            #print i,fields[0]
            i=i+1



    def convarray(self):
        self.id = N.array(self.id,'i')
        self.mass = N.array(self.mass,'f')
        self.x1 = simL*N.array(self.x1,'f')
        self.x2 = simL*N.array(self.x2,'f')
        self.x3 = simL*N.array(self.x3,'f')
        self.v1 = N.array(self.v1,'f')
        self.v2 = N.array(self.v2,'f')
        self.v3 = N.array(self.v3,'f')
        self.sigma = N.array(self.sigma,'f')
        #self.rvir=3*G*self.mass*1.e10/(self.sigma**2)
        self.rvir=1.73/N.sqrt(2.)*self.sigma/1000.
        #self.rvir=1.5*N.ones(len(self.mass),'f')
        self.x=self.x1
        self.y=self.x2
        #self.z=(self.x3*H0+self.v3)/3.e5

    def readGIFfiles(self,halocat):

        infile=open(halocat,'r')
        i=0
        for line in infile:
            if line.find('#') > -1:
                continue
            fields=line.split()
            if len(fields) < 4:#get halo index and move to next line
                continue
            mcl=10.**(float(fields[1])-12.)#halo mass in units of 10^12
            #if mcl > Mclmin:
            #print "heyah",fields[0],Mclmin
            self.id.append(float(i+1))
            self.mass.append(10.**(float(fields[1])-12.))#10^12 Msun
            self.rvir.append(float(fields[2]))#R200 in h^{-1} Mpc
            self.x1.append(float(fields[4]))
            self.x2.append(float(fields[5]))
            self.x3.append(float(fields[6]))
            #print i,fields[0]
            i=i+1
    def convGIFarray(self):
        self.id = N.array(self.id,'i')
        self.mass = N.array(self.mass,'f')
        self.x1 = N.array(self.x1,'f')
        self.x2 = N.array(self.x2,'f')
        self.x3 = N.array(self.x3,'f')
        self.v1 = N.zeros(len(self.id),'f')
        self.v2 = N.zeros(len(self.id),'f')
        self.v3 = N.zeros(len(self.id),'f')
        self.rvir=N.array(self.rvir,'f')
        #self.v3 = self.x3*H0
        #if zbox > .05:
        #    self.v3=(self.v3+1.)*3.e5*((1+zbox)**2-1)/((1+zbox)**2 + 1)
        self.x=self.x1
        self.y=self.x2
    def dogif(self):
        self.readGIFfiles(halocat)
        self.convGIFarray()
        self.getvobs()
        self.getgifsigma()

    def getvobs(self):
        vh=self.x3*H0
        self.vobs=(vh + self.v3)/(1+vh*self.v3/(3.e5)**2)#Hubble flow
        self.z=N.sqrt((1+self.vobs/3.e5)/(1-self.vobs/3.e5)) - 1.#observed redshift
        #print "vobs, z",self.vobs[0:5],self.z[0:5]
        #self.z=self.vobs/3.e5#observed redshift
    def getgifsigma(self):
        #self.sigma = self.rvir/1.73*N.sqrt(omegaL+omega0*(1+self.z)**3)*1000.
        self.sigma = N.sqrt(G*self.mass*1.e12/self.rvir)
        #for i in range(len(self.sigma)):
        #    if self.mass[i] > Mclmin:
        #        print i,self.mass[i],self.rvir[i],self.sigma[i]

    def cutmass(self):
        #if mcl > Mclmin:
        self.id = N.compress(self.mass > Mclmin,self.id)
        self.x1 = N.compress(self.mass > Mclmin,self.x1)
        self.x2 = N.compress(self.mass > Mclmin,self.x2)
        self.x3 = N.compress(self.mass > Mclmin,self.x3)
        self.v1 = N.compress(self.mass > Mclmin,self.v1)
        self.v2 = N.compress(self.mass > Mclmin,self.v2)
        self.v3 = N.compress(self.mass > Mclmin,self.v3)
        self.rvir=N.compress(self.mass > Mclmin,self.rvir)
        self.sigma = N.compress(self.mass > Mclmin,self.sigma)
        self.vobs = N.compress(self.mass > Mclmin,self.vobs)
        self.z = N.compress(self.mass > Mclmin,self.z)

        self.mass = N.compress(self.mass > Mclmin,self.mass)#LAST LINE!!!
        #print "dude ",len(self.x3),len(self.vobs)
        #for i in range(10):
        #    print "vobs %i %5i %6.2f %6.2f %6.2f" % (i,self.id[i],self.x3[i],(self.x3[i]*H0),self.vobs[i])
            #print "vobs",i,self.x3[i],(self.x3[i]*H0),self.vobs[i])

    def getmemberidsall(self,nv,nr,cg):#get galaxy members
        print "within getmemberids(), nv,nr,len(x1),len(c.mass) = ",nv,nr,len(cg.x1),len(self.mass)
        self.membids=N.zeros(len(self.mass),'f')
        self.membincut=N.zeros(len(self.mass),'f')
        self.membids=list(self.membids)
        for i in range(len(self.mass)):
            deltaid=0.1#tolerance for matching
            self.membids[i]=[]

            (self.membids[i],matchflag)=my.findmatch(self.id[i],cg.ih,deltaid)#g.ih is already sorted
            #if i < 5:
            #    print i,self.id[i],self.membids[i],matchflag
            cdtheta=self.rvir[i]/my.DA(self.z[i],h)
	    #print 'cdtheta = ',cdtheta
            try:
                a=len(self.membids[i])
            except TypeError:
                #if (self.membids[i]+999) < 1.:
                #print "no members in cluster ",i,self.membids[i]
                self.membincut[i]=0
                continue
            if len(self.membids[i]) == 1:
                #print "one member in cluster ",i,self.membids[i]
                membids=self.membids[i]
                j=int(membids[0])
                dr=(N.sqrt((cg.x1[j]-self.x1[i])**2+(cg.x2[j]-self.x2[i])**2))/self.rvir[i]
                if dr < nr:
                    #dv = abs(self.vobs[i]-cg.vobs[j])*(1.+(self.vobs[i]*cg.vobs[j]/((3.e5)**2)))
                    dv = abs(self.z[i]-cg.z[j])*3.e5/(1+self.z[i])#postman 2001
                    if dv <= nv*N.sqrt(3.)*cg.hsigmax[j]:
                        gdtheta=dr/my.DA(cg.z[j],h)

                        if gdtheta < cdtheta:#make sure galaxy is in cone field-of-view
                            self.membincut[i]=1

                continue
            else:
                try:
		    #print "looping over memberids"
                    for j in self.membids[i]:
                        j=int(j)
                        dr=(N.sqrt((cg.x1[j]-self.x1[i])**2+(cg.x2[j]-self.x2[i])**2))/self.rvir[i]
                        if dr < nr:
			    #print "made radial cut"
                            #dv = abs(self.vobs[i]-cg.vobs[j])/self.sigma[i]#/(1.+(self.vobs[i]*g.vobs[j]/((3.e5)**2)))
                            dv = abs(self.z[i]-cg.z[j])*3.e5/(1+self.z[i])#postman 2001
                            #dv=abs(self.x3[i]*H0-(g.v3[j]+g.x3[j]*H0)) 
			    #print "dv = ",dv, nv
                            if abs(dv) <= nv*N.sqrt(3.)*cg.hsigmax[j]:#self.sigma[i]:
				#print "made vel cut"
                                gdtheta=dr/my.DA(cg.z[j],h)
                                if gdtheta <= cdtheta:#make sure galaxy is in cone field-of-view

                                    self.membincut[i]=self.membincut[i]+1
                except TypeError:
                    continue
	    #print "cluster ",i," self.membincut[i] = ",self.membincut[i]
    def getobsmemberidsall(self,nv,nr,cg):#get galaxies within nr*dr and nv*dv
        #print "within getobsmemberids(), nv,nr,len(x1) = ",nv,nr,len(g.x1)
        print "within measurecontam(), nv,nr,len(x1) = ",nv,nr,len(cg.x1),len(self.mass)
        #print self.mass
        self.obsmembids=N.zeros(len(self.mass),'f')
        self.lcontam=N.zeros(len(self.mass),'f')
        self.lcompl=N.zeros(len(self.mass),'f')
        self.obsmembids=list(self.obsmembids)
        self.obsmembincut=N.zeros(len(self.mass),'f')
        self.membflag=N.zeros(len(self.mass),'f')
        self.membflag=list(self.membflag)
        self.typeflag=N.zeros(len(self.mass),'f')
        self.typeflag=list(self.typeflag)
        self.obsmembdv=N.zeros(len(self.mass),'f')
        self.obsmembdv=list(self.typeflag)
        self.obsmembdr=N.zeros(len(self.mass),'f')
        self.obsmembdr=list(self.typeflag)
        self.lbjcontam=N.zeros(len(self.mass),'f')
        self.lbjcontam=list(self.lbjcontam)
        vbox=my.vofz(zbox)#redshift of box center given epoch of observation
        #i.e. recession velocity corresponding to z=0.2, 0.4, 0.6, etc of GIF observations
        for i in range(len(self.mass)):
            #self.obsmembids[i]=[]
            vcl=H0*self.x3[i]
            deltar=nr*self.rvir[i]#tolerance for matching
            deltav=nv*self.sigma[i]#tolerance for matching
            #print i,"cluster id, mass, rvir, sigma = ",self.id[i],self.mass[i],self.rvir[i],self.sigma[i],Mclmin
            temp=[]
            #put cluster in middle of box
            x1=cg.x1-self.x1[i]+0.5*simL
            #flip coordinates so galaxy x1 ranges from 0 to simL 
            for j in range(len(x1)):
                if x1[j] > simL:
                    x1[j]=x1[j]-simL
                if x1[j] < 0:
                    x1[j]=x1[j]+simL
            x2=cg.x2-self.x2[i]+0.5*simL
            for j in range(len(x2)):
                if x2[j] > simL:
                    x2[j]=x2[j]-simL
                if x2[j] < 0:
                    x2[j]=x2[j]+simL
            x3=cg.x3-self.x3[i]+0.5*simL
            for j in range(len(x3)):
                if x3[j] > simL:
                    x3[j]=x3[j]-simL
                if x3[j] < 0:
                    x3[j]=x3[j]+simL
            gvobs=x3*H0 + cg.v3
            #print "min, max x3 = ",min(x3), max(x3),simL
            gvobs=my.relsumv(vbox,gvobs)
            gz=my.zofv(gvobs)
            (x1sort,x1index)=my.sortwindex(x1)
            #print "cluster ",i,self.id[i],self.x1[i],self.x2[i],self.x3[i],self.v1[i]
            #print "length of x1 ",len(x1),len(x1sort),len(x1index)
            #print "sorting galaxy catalog on x2"
            (x2sort,x2index)=my.sortwindex(x2)

            #temp=my.findmatch(self.x1[i],x1sort,deltar)
            #print "max x1 = ",max(x1sort),min(x1sort),0.5*simL,self.x1[i],self.id[i]
            (temp,matchflag)=my.findmatch(0.5*simL,x1sort,deltar)
            templist1=temp
            for j in range(len(temp)):
		#print j,temp[j],len(templist1),len(x1index)
                templist1[j]=x1index[temp[j]]
            #print i,"len(x2sort) = ",len(x2sort),deltar,self.x2[i]
            temp=[]
            #temp=my.findmatch(self.x2[i],x2sort,deltar)
            #print "length of x2 ",len(x2),len(x2sort),len(x2index),self.x2[i],self.id[i]
            #print "max min x2 = ",max(x2sort),min(x2sort),0.5*simL
                        
            (temp,matchflag)=my.findmatch(0.5*simL,x2sort,deltar)
            templist2=temp
            for j in range(len(temp)):
                templist2[j]=x2index[temp[j]]
            a=sets.Set(templist1)
            b=sets.Set(templist2)
            members=a&b
            members=list(members)
            #print i, "potential members = ",len(members)
            templist1=N.argsort(templist1)
            templist2=N.argsort(templist2)
            subset=[]
            mflag=[]
            tflag=[]
            ltemp=[]
            dvtemp=[]
            drtemp=[]
            #cdtheta=abs(self.rvir[i]/self.x3[i])
            #cdtheta=self.rvir[i]/my.DA(self.z[i],h)
            ccvobs=my.relsumv(vbox,0.5*simL*H0)
            ccz=my.zofv(ccvobs)
            #cdtheta=self.rvir[i]/my.DA(0.5*simL*H0,h)
            cdtheta=self.rvir[i]/my.DA(ccz,h)
            ccx1=0.5*simL
            ccx2=0.5*simL
            for k in members:
                #dr=(N.sqrt((x1[k]-self.x1[i])**2+(x2[k]-self.x2[i])**2))/self.rvir[i]
                dr=(N.sqrt((x1[k]-ccx1)**2+(x2[k]-ccx2)**2))/self.rvir[i]
                if dr < nr:
                    #dv = (self.vobs[i]-g.vobs[k])/(N.sqrt(3.)*g.hsigmax[k])
                    #dv = (self.vobs[i]-gvobs[k])/(self.sigma[i])
                    #dv = (ccvobs-gvobs[k])/(self.sigma[i])
                    dv = (gvobs[k]-ccvobs)/(1-gvobs[k]*ccvobs/(9.e10))/self.sigma[i]
                    #dv = 3.e5*(gz[k]-ccz)/(1+ccz)/(self.sigma[i])
                    #print "dv = ",dv,nv,self.sigma[i],ccvobs,gvobs[k],cg.x3[k],cg.v3[k],self.x3[i],x3[k],vbox,0.5*simL*H0
                    #dv=abs(vcl-(g.v3[k]+g.x3[k]*H0))
                    if abs(dv) <= nv:#deltav:
                        #print "made velocity cut with dv = ", dv
                        #gdtheta=abs(dr/g.x3[k])
                        gdtheta=dr*self.rvir[i]/my.DA(gz[k],h)
                        #print "gdtheta = ",gdtheta, cdtheta
                        if gdtheta <= cdtheta:#make sure galaxy is in cone field-of-view
                            #print "made angular cut with gdtheta = ",gdtheta
                            subset.append(k)
                            if abs(cg.ih[k]-self.id[i]) < .5:
                                mflag.append(1)
                                dvtemp.append(dv)
                                drtemp.append(dr)
                            else:
                                mflag.append(0)
                                tflag.append(cg.type[k])
                                ltemp.append(cg.lbj[k])
                                dvtemp.append(dv)
                                drtemp.append(dr)
            #print i,"number of obs members = ",len(subset)
            if len(subset) > 0:
                t=N.array(tflag,'f')
                s=N.array(subset,'f')
                lcontam=1.*N.sum(t)/(1.*len(s))
            else:
                lcontam=-999.

            self.obsmembids[i]=subset
            self.obsmembincut[i]=len(subset)
            #print i,"number of obs members = ",len(subset),self.obsmembincut[i]
            self.membflag[i]=mflag
            self.typeflag[i]=tflag
            self.lbjcontam[i]=ltemp
            self.obsmembdv[i]=dvtemp
            self.obsmembdr[i]=drtemp
            self.lcontam[i]=lcontam
	    #print "within getobsmemberidsall, max(obsmembids)=",max(subset)," len(cg.x3) = ",len(cg.x3)
    def getmemberids(self,nv,nr,cg):#get galaxy members
        print "within getmemberids(), nv,nr,len(x1),len(c.mass) = ",nv,nr,len(cg.x1),len(self.mass)
        for i in range(len(self.mass)):
            deltaid=0.1#tolerance for matching
            cdtheta=self.rvir[i]/my.DA(self.z[i],h)
	    allmembids=self.membidsall[i]
            try:
                a=len(allmembids)
            except TypeError:
                #if (self.membids[i]+999) < 1.:
                #print "no members in cluster ",i,self.membids[i]
                self.membincut[i]=0
                continue
            if len(allmembids) == 1:
                #print "one member in cluster ",i,self.membids[i]
                membids=allmembids
                j=int(membids[0])
                dr=(N.sqrt((cg.x1[j]-self.x1[i])**2+(cg.x2[j]-self.x2[i])**2))/self.rvir[i]
                if dr < nr:
                    #dv = abs(self.vobs[i]-cg.vobs[j])*(1.+(self.vobs[i]*cg.vobs[j]/((3.e5)**2)))
                    dv = abs(self.z[i]-cg.z[j])*3.e5/(1+self.z[i])#postman 2001
                    if dv <= nv*N.sqrt(3.)*cg.hsigmax[j]:
                        gdtheta=dr/my.DA(cg.z[j],h)

                        if gdtheta < cdtheta:#make sure galaxy is in cone field-of-view
                            self.membincut[i]=1

                continue
            else:
                try:
		    #print "looping over memberids"
                    for j in self.membids[i]:
                        j=int(j)
                        dr=(N.sqrt((cg.x1[j]-self.x1[i])**2+(cg.x2[j]-self.x2[i])**2))/self.rvir[i]
                        if dr < nr:
			    #print "made radial cut"
                            #dv = abs(self.vobs[i]-cg.vobs[j])/self.sigma[i]#/(1.+(self.vobs[i]*g.vobs[j]/((3.e5)**2)))
                            dv = abs(self.z[i]-cg.z[j])*3.e5/(1+self.z[i])#postman 2001
                            #dv=abs(self.x3[i]*H0-(g.v3[j]+g.x3[j]*H0)) 
			    #print "dv = ",dv, nv
                            if abs(dv) <= nv*N.sqrt(3.)*cg.hsigmax[j]:#self.sigma[i]:
				#print "made vel cut"
                                gdtheta=dr/my.DA(cg.z[j],h)
                                if gdtheta <= cdtheta:#make sure galaxy is in cone field-of-view

                                    self.membincut[i]=self.membincut[i]+1
                except TypeError:
                    continue
	    #print "cluster ",i," self.membincut[i] = ",self.membincut[i]
    def getobsmemberids(self,nv,nr,cg):#get galaxies within nr*dr and nv*dv
        #print "within getobsmemberids(), nv,nr,len(x1) = ",nv,nr,len(g.x1)
        print "within measurecontam(), nv,nr,len(x1) = ",nv,nr,len(cg.x1),len(self.mass)
        #print self.mass

        vbox=my.vofz(zbox)#redshift of box center given epoch of observation
        for i in range(len(self.mass)):
            #self.obsmembids[i]=[]
            vcl=H0*self.x3[i]
            deltar=nr*self.rvir[i]#tolerance for matching
            deltav=nv*self.sigma[i]#tolerance for matching
            #print i,"cluster id, mass, rvir, sigma = ",self.id[i],self.mass[i],self.rvir[i],self.sigma[i],Mclmin
            temp=[]
            #put cluster in middle of box
            x1=cg.x1-self.x1[i]+0.5*simL
            #flip coordinates so galaxy x1 ranges from 0 to simL 
            for j in range(len(x1)):
                if x1[j] > simL:
                    x1[j]=x1[j]-simL
                if x1[j] < 0:
                    x1[j]=x1[j]+simL
            x2=cg.x2-self.x2[i]+0.5*simL
            for j in range(len(x2)):
                if x2[j] > simL:
                    x2[j]=x2[j]-simL
                if x2[j] < 0:
                    x2[j]=x2[j]+simL
            x3=cg.x3-self.x3[i]+0.5*simL
            for j in range(len(x3)):
                if x3[j] > simL:
                    x3[j]=x3[j]-simL
                if x3[j] < 0:
                    x3[j]=x3[j]+simL
            gvobs=x3*H0 + cg.v3
            #print "min, max x3 = ",min(x3), max(x3),simL
            gvobs=my.relsumv(vbox,gvobs)
            gz=my.zofv(gvobs)

            members=self.obsmembidsall[i]
            #print i, "potential members = ",len(members)#,members
            #templist1=N.argsort(templist1)
            #templist2=N.argsort(templist2)
            subset=[]
            mflag=[]
            tflag=[]
            ltemp=[]
            dvtemp=[]
            drtemp=[]
            #cdtheta=abs(self.rvir[i]/self.x3[i])
            #cdtheta=self.rvir[i]/my.DA(self.z[i],h)
            ccvobs=my.relsumv(vbox,0.5*simL*H0)
            ccz=my.zofv(ccvobs)
            #cdtheta=self.rvir[i]/my.DA(0.5*simL*H0,h)
            cdtheta=self.rvir[i]/my.DA(ccz,h)
            ccx1=0.5*simL
            ccx2=0.5*simL
            for k in members:
                #dr=(N.sqrt((x1[k]-self.x1[i])**2+(x2[k]-self.x2[i])**2))/self.rvir[i]
		try:
		    dr=(N.sqrt((x1[k]-ccx1)**2+(x2[k]-ccx2)**2))/self.rvir[i]
		except IndexError:
		    print k,len(x1),len(x2),i,len(self.rvir)
                if dr < nr:
                    #dv = (self.vobs[i]-g.vobs[k])/(N.sqrt(3.)*g.hsigmax[k])
                    #dv = (self.vobs[i]-gvobs[k])/(self.sigma[i])
                    #dv = (ccvobs-gvobs[k])/(self.sigma[i])
		    try:
			dv = (gvobs[k]-ccvobs)/(1-gvobs[k]*ccvobs/(9.e10))/self.sigma[i]
		    except IndexError:
			print "IndexError 2",len(gvobs),k,len(self.sigma),i
                    #dv = 3.e5*(gz[k]-ccz)/(1+ccz)/(self.sigma[i])
                    #print "dv = ",dv,nv,self.sigma[i],ccvobs,gvobs[k],cg.x3[k],cg.v3[k],self.x3[i],x3[k],vbox,0.5*simL*H0
                    #dv=abs(vcl-(g.v3[k]+g.x3[k]*H0))
                    if abs(dv) <= nv:#deltav:
                        #print "made velocity cut with dv = ", dv
                        #gdtheta=abs(dr/g.x3[k])
                        gdtheta=dr*self.rvir[i]/my.DA(gz[k],h)
                        #print "gdtheta = ",gdtheta, cdtheta
                        if gdtheta <= cdtheta:#make sure galaxy is in cone field-of-view
                            #print "made angular cut with gdtheta = ",gdtheta
                            subset.append(k)
                            if abs(cg.ih[k]-self.id[i]) < .5:
                                mflag.append(1)
                                dvtemp.append(dv)
                                drtemp.append(dr)
                            else:
                                mflag.append(0)
                                tflag.append(cg.type[k])
                                ltemp.append(cg.lbj[k])
                                dvtemp.append(dv)
                                drtemp.append(dr)
            #print i,"number of obs members = ",len(subset)
            if len(subset) > 0:
                t=N.array(tflag,'f')
                s=N.array(subset,'f')
                lcontam=1.*N.sum(t)/(1.*len(s))
            else:
                lcontam=-999.

            self.obsmembids[i]=subset
            self.obsmembincut[i]=len(subset)
            #print i,"number of obs members = ",len(subset),self.obsmembincut[i]
            self.membflag[i]=mflag
            self.typeflag[i]=tflag
            self.lbjcontam[i]=ltemp
            self.obsmembdv[i]=dvtemp
            self.obsmembdr[i]=drtemp
            self.lcontam[i]=lcontam
    def calccontam(self,cg):
            self.contam=N.zeros(len(self.mass),'f')
            #self.lcontam=N.zeros(len(self.mass),'f')
            self.compl=N.zeros(len(self.mass),'f')
            self.contamtype=N.zeros(len(self.mass),'f')
            for i in range(len(self.mass)):
                #if i < 5:
                #    print i,len(self.membids[i]),self.membincut[i],len(self.obsmembids[i]),N.sum(self.membflag[i]),self.sigma[i],self.rvir[i]

                try:
                    self.compl[i]=float(N.sum(self.membflag[i]))/float(len(self.membids[i]))
                except:
                    self.compl[i]=-999.
                try:
                    self.contam[i]=float(len(self.obsmembids[i]) - N.sum(self.membflag[i]))/float(len(self.obsmembids[i]))
                except ZeroDivisionError:
                    self.contam[i]=-999.
                    

                try:
                    self.contamtype[i]=N.average(self.typeflag[i])
                except ZeroDivisionError:
                    self.contamtype[i]=-999
                #print "%i %5.2f %5.2f %5.2f" % (i,self.compl[i],self.contam[i],self.contamtype[i])
                memb=self.membids[i]
                try:
                    t=len(memb)
                    suml=0.
                    for j in range(len(memb)):
                        if cg.type[int(j)] > 0.1:
                            suml+=1 #number of late-type members
                    #print "Number of member late types = ",suml
                    obsmemb=self.obsmembids[i]
                    a=sets.Set(memb)
                    b=sets.Set(obsmemb)
                    templist1=a&b#members that are observed
                    sumol=0.#sum of observed late types
                    for j in templist1:
                        if cg.type[int(j)] > 0.1:
                            sumol+=1. #number of late-type members
                    #print "Number of observed, member late-types",sumol
                    self.lcompl[i]=1.*sumol/(1.*suml)
                    #print "Compl for late-types obs,memb,frac",sumol,suml,self.lcompl[i]
                except TypeError:
                    self.lcompl[i]=-999.
                except ZeroDivisionError:
                    self.lcompl[i]=-999.



            compl=N.compress(self.contam > 0.,self.compl)
            print "average completeness = %5.2f" % (N.average(compl))
            contam=N.compress(self.contam > 0.,self.contam)
            print "average contamination = %5.2f" % (N.average(contam))
            contamtype=N.compress(self.contam > 0.,self.contamtype)
            (a,b)=my.ratioerror(N.sum(contamtype),float(len(contamtype)))
            print "average type of contaminating galaxies = %5.2f +/- %5.2f" % (a,b)
	    try:
		print "average type (0/early,1/late) of contaminating galaxies = %5.2f +/- %5.2f" % (N.average(contamtype),pylab.std(contamtype))
	    except ZeroDivisionError:
		print "average type (0/early,1/late) of contaminating galaxies = %5.2f" % (N.average(contamtype))


    def calctotalsfr(self,cg):
            self.totalsfr=N.zeros(len(self.mass),'f')
            for i in range(len(self.mass)):
                try:
                    nmemb=len(self.membids[i])
                    obsmemb=self.membids[i]
                    for j in obsmemb:
                        j = int(j)
                        self.totalsfr[i] = self.totalsfr[i] + cg.sfr[j]
                except:
                    self.totalsfr[i]=self.totalsfr[i]

    def getobsmemberidsdr(self,nv,nr1,nr2):#get galaxies within nr*dr and nv*dv
        #nr2 > nr1
        self.obsmembids=N.zeros(len(self.mass),'f')
        self.obsmembids=list(self.obsmembids)
        self.membids=N.zeros(len(self.mass),'f')
        self.membids=list(self.obsmembids)
        self.membflag=N.zeros(len(self.mass),'f')
        self.membflag=list(self.membflag)
        self.typeflag=N.zeros(len(self.mass),'f')
        self.typeflag=list(self.typeflag)
        self.lbjcontam=N.zeros(len(self.mass),'f')
        self.lbjcontam=list(self.lbjcontam)
        self.contam=N.zeros(len(self.mass),'f')
        self.contam=list(self.contam)
        self.compl=N.zeros(len(self.mass),'f')
        self.compl=list(self.compl)
        self.lcontam=N.zeros(len(self.mass),'f')#completeness for late type
        self.lcontam=list(self.contam)
        self.lcompl=N.zeros(len(self.mass),'f')
        self.lcompl=list(self.compl)
        self.econtam=N.zeros(len(self.mass),'f')#completeness for early type
        self.econtam=list(self.contam)
        self.ecompl=N.zeros(len(self.mass),'f')
        self.ecompl=list(self.compl)
        nol=N.zeros(len(self.mass),'f')#numb observed late-type
        noe=N.zeros(len(self.mass),'f')
        nml=N.zeros(len(self.mass),'f')#number of actual members
        nme=N.zeros(len(self.mass),'f')
        noml=N.zeros(len(self.mass),'f')#numb of observed members
        nome=N.zeros(len(self.mass),'f')
        (x1sort,x1index)=my.sortwindex(g.x1)
        (x2sort,x2index)=my.sortwindex(g.x2)
        for i in range(len(self.mass)):
            #self.obsmembids[i]=[]
            vcl=H0*self.x3[i]
            deltar1=nr1*self.rvir[i]#tolerance for matching
            deltar2=nr2*self.rvir[i]#tolerance for matching
            deltav=nv*self.sigma[i]#tolerance for matching

            if (nr1 > 0.):
                print "in here!"
                temp=[]
                (temp,matchflag)=my.findmatch(self.x1[i],x1sort,deltar1)
                templist1a=temp
                for j in range(len(temp)):
                    templist1a[j]=x1index[temp[j]]
                    
                
                temp=[]
                (temp,matchflag)=my.findmatch(self.x1[i],x1sort,deltar2)
                templist1b=temp
                for j in range(len(temp)):
                    templist1b[j]=x1index[temp[j]]

                a=sets.Set(templist1a)
                b=sets.Set(templist1b)
                templist1=a&b
                templist1=list(templist1)#number in x-direction slices
                print "len templist1 = ",len(templist1)
                
                temp=[]
                (temp,matchflag)=my.findmatch(self.x2[i],x2sort,deltar1)
                templist2a=temp
                for j in range(len(temp)):
                    templist2a[j]=x2index[temp[j]]

                temp=[]
                (temp,matchflag)=my.findmatch(self.x2[i],x2sort,deltar2)
                templist2b=temp
                for j in range(len(temp)):
                    templist2b[j]=x2index[temp[j]]

                a=sets.Set(templist2a)
                b=sets.Set(templist2b)
                templist2=a&b
                templist2=list(templist2)
                print "len templist2 = ",len(templist2)

                a=sets.Set(templist1)
                b=sets.Set(templist2)
                members=a&b
                print "len members = ",len(members),self.x1[i],self.x2[i]
            if (nr1 < .01):
                deltar=nr2*self.rvir[i]
                temp=[]
                (temp,matchflag)=my.findmatch(self.x1[i],x1sort,deltar)
                templist1=temp
                for j in range(len(temp)):
                    templist1[j]=x1index[temp[j]]
                temp=[]
                (temp,matchflag)=my.findmatch(self.x2[i],x2sort,deltar)
                templist2=temp
                for j in range(len(temp)):
                    templist2[j]=x2index[temp[j]]
                a=sets.Set(templist1)
                b=sets.Set(templist2)
                members=a&b
            members=list(members)
            templist1=N.argsort(templist1)
            templist2=N.argsort(templist2)
            subset=[]
            mflag=[]
            tflag=[]
            ltemp=[]
            memb=[]
            mtype=[]#type for members
            #print nr1,nr2,"len of potential members = ",len(members),members,len(g.x1)
            for k in members:
                dr=N.sqrt((g.x1[k]-self.x1[i])**2+(g.x2[k]-self.x2[i])**2)
                if dr < deltar2:
                    if dr > deltar1:
                        dv = abs(self.vobs[i]-g.vobs[k])#/(1+(self.vobs[i]*g.vobs[k]/((3.e5)**2)))
                        drphys=N.sqrt((g.x1[k]-self.x1[i])**2+(g.x2[k]-self.x2[i])**2 + (g.x3[k]-self.x3[i])**2)
                        if (drphys > deltar1) & (drphys < deltar2):
                            memb.append(k)
                            mtype.append(g.type[k])
                            if g.type[k] > 0:
                                nml[i]=nml[i]+1
                            if g.type[k] < 1:
                                nme[i]=nme[i]+1
                        if dv < deltav:
                            subset.append(k)
                            if g.type[k] > 0:
                                nol[i]=nol[i]+1
                            if g.type[k] < 1:
                                noe[i]=noe[i]+1

                            if (drphys > deltar1) & (drphys < deltar2):
                                mflag.append(1)
                                tflag.append(g.type[k])
                                ltemp.append(g.lbj[k])
                                if g.type[k] > 0:
                                    noml[i]=noml[i]+1
                                if g.type[k] < 1:
                                    nome[i]=nome[i]+1

                            else:
                                mflag.append(0)
                                tflag.append(g.type[k])
                                ltemp.append(g.lbj[k])
            self.obsmembids[i]=subset
            self.membflag[i]=mflag
            self.typeflag[i]=tflag
            self.lbjcontam[i]=ltemp
            self.membids[i]=memb
            tflag=N.array(tflag,'f')
            mflag=N.array(mflag,'f')
            mtype=N.array(mtype,'f')
            if len(subset) == 0:
                print "Warning - no observed members!",i,nr1,nr2,deltar1,deltar2,len(members),len(subset),len(memb)
                self.contam[i]=-999
                self.compl[i]=-999
                self.lcontam[i]=-999
                self.lcompl[i]=-999
                self.econtam[i]=-999
                self.ecompl[i]=-999
                continue
            if len(memb) == 0:
                print "Warning - no actual members!",i,nr1,nr2,deltar1,deltar2,len(members),len(subset),len(memb)
                
                self.contam[i]=1.
                self.compl[i]=0.
                self.lcompl[i]=0.
                self.ecompl[i]=0.
                self.lcontam[i]=-999
                self.econtam[i]=-999
                if (nol[i] > 0):
                    self.lcontam[i]=1.
                if (noe[i] > 0):
                    self.econtam[i]=1.
                continue

            #self.contam[i]=float(len(self.obsmembids[i]) - N.sum(self.membflag[i]))/float(len(self.obsmembids[i]))
            try:
                #self.contam[i]= float(len(self.obsmembids[i]) - N.sum(self.membflag[i]))/float(len(self.obsmembids[i]))
                self.contam[i]=(nol[i]+noe[i]-noml[i]-nome[i])/(nol[i]+noe[i])
            except ZeroDivisionError:
                self.contam[i]=-999
            try: 
                #self.compl[i]=float(N.sum(self.membflag[i]))/float(len(self.membids[i]))
                self.compl[i]=(noml[i]+nome[i])/(nml[i]+nme[i])
            except ZeroDivisionError:
                self.compl[i]=0.
            #print "len tflag, mflag",len(tflag),len(mflag),nr1,nr2
            #print tflag
            #print mflag
            #print memb
            #memb=N.array(memb,'f')
            #a=float(N.sum(N.compress(tflag > 0,mflag)))
            #b=float(len(N.compress(mtype > 0,subset)))
            try:
                #self.lcontam[i]=float(len(N.compress(tflag > 0,self.obsmembids[i]))-N.sum(N.compress(tflag > 0,mflag)))/float(len(N.compress(mtype > 0,self.obsmembids[i])))
                self.lcontam[i]=(nol[i]-noml[i])/nml[i]

            except ZeroDivisionError:
                self.lcontam[i]=1.
            try:
                #self.lcompl[i]=float(N.sum(N.compress(tflag > 0,mflag)))/float(len(N.compress(mtype > 0,memb)))
                self.lcompl[i]=noml[i]/nml[i]
            except ZeroDivisionError:
                self.lcompl[i]=0.

            try:
                #self.econtam[i]=float(len(N.compress(tflag > 0,self.obsmembids[i]))- N.sum(N.compress(tflag < 1.,mflag)))/float(len(N.compress(mtype < 1.,memb)))
                self.econtam[i]=(noe[i]-nome[i])/nme[i]
            except ZeroDivisionError:
                self.econtam[i]=1.

            try:
                #self.ecompl[i]=float(N.sum(N.compress(tflag < 1.,mflag)))/float(len(N.compress(mtype < 1.,memb)))
                self.ecompl[i]=nome[i]/nme[i]
            except ZeroDivisionError:
                self.ecompl[i]=0.
            


    #def measurecontam(self,nv,nr,gx1,gx2,gih,gvobs,gz,ghsigmax,gtype,glbj):
    def getallmembers(self,cg,nv,nr):#get member and observed ids within 6 sig and 3 Rv
	print "in get all members, len(g.v3) = ",len(cg.v3)

	#nv=6
	#nr=2
        self.getmemberidsall(nv,nr,cg)
	print "after getmemberidsall , len(g.v3) = ",len(cg.v3)
        self.getobsmemberidsall(nv,nr,cg)
	print "in get all members, after getobsmemberidsall len(g.v3) = ",len(cg.v3)
	self.membidsall=self.membids
	print "after membidsall , len(g.v3) = ",len(cg.v3)
	self.obsmembidsall=self.obsmembids
	print "after obsmembidsall , len(g.v3) = ",len(cg.v3)
    def measurecontam(self,nv,nr,cg):
        myg=cg
        print "within measurecontam(), nv,nr,len(x1) = ",nv,nr,len(cg.x1),len(self.mass)
        #print "Getting members"
        #self.getmemberids(nv,nr,gx1,gx2,gih,gvobs,gz,ghsigmax)
        #print "before getmemberids ",len(cg.ih),len(cg.z),len(cg.x3)
        self.getmemberids(nv,nr,cg)
        #print "within measurecontam(), nv,nr,len(x1) = ",nv,nr,len(gx1)
        #print "Getting apparent members"
        #self.getobsmemberids(nv,nr,gx1,gx2,gz,gvobs,gtype,glbj)
        #self.getobsmemberids(nv,nr,cg.x1,cg.x2,cg.z,cg.vobs,cg.type,cg.lbj)
        #print self.x1[0:10]
        #print self.x2[0:10]
        self.getobsmemberids(nv,nr,cg)
        #for i in range(10):
        #    membids=self.membids[i]
        #    try:
        #        a=len(membids)
        #    except:
        #        continue
        #    for j in range(len(membids)):                
        #        k=int(membids[j])
        #        dr=N.sqrt((self.x1[i]-cg.x1[k])**2+(self.x2[i]-cg.x2[k])**2)
        #        dr=dr/self.rvir[i]
        #        dv=abs(self.vobs[i]-cg.vobs[k])
        #        dv=dv/self.sigma[i]
        #        print "member stats %5i %5.2f %5.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f" % (self.id[i],dr,dv,self.sigma[i],self.x3[i],self.x3[i]*H0,self.v3[i],self.vobs[i],cg.x3[k],cg.v3[k],cg.vobs[k])

        #print "Calculating contamination"
        self.calccontam(cg)

    def measurengalcontam(self,nv,nr,cg):
        myg=cg
        print "within measurecontam(), nv,nr,len(x1) = ",nv,nr,len(cg.x1),len(self.mass)
        #print "Getting members"
        #self.getmemberids(nv,nr,gx1,gx2,gih,gvobs,gz,ghsigmax)
        #print "before getmemberids ",len(cg.ih),len(cg.z),len(cg.x3)
        self.getmemberids(nv,nr,cg)
        #print "within measurecontam(), nv,nr,len(x1) = ",nv,nr,len(gx1)
        #print "Getting apparent members"
        #self.getobsmemberids(nv,nr,gx1,gx2,gz,gvobs,gtype,glbj)
        #self.getobsmemberids(nv,nr,cg.x1,cg.x2,cg.z,cg.vobs,cg.type,cg.lbj)
        #print self.x1[0:10]
        #print self.x2[0:10]
        self.getobsmemberids(nv,nr,cg)
        #for i in range(10):
        #    membids=self.membids[i]
        #    try:
        #        a=len(membids)
        #    except:
        #        continue
        #    for j in range(len(membids)):                
        #        k=int(membids[j])
        #        dr=N.sqrt((self.x1[i]-cg.x1[k])**2+(self.x2[i]-cg.x2[k])**2)
        #        dr=dr/self.rvir[i]
        #        dv=abs(self.vobs[i]-cg.vobs[k])
        #        dv=dv/self.sigma[i]
        #        print "member stats %5i %5.2f %5.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f" % (self.id[i],dr,dv,self.sigma[i],self.x3[i],self.x3[i]*H0,self.v3[i],self.vobs[i],cg.x3[k],cg.v3[k],cg.vobs[k])

        #print "Calculating contamination"
        #self.calccontam(cg)

    def readobsmembids(self,file):
        self.obsmembids=N.zeros(len(self.mass),'f')
        self.obsmembids=list(self.obsmembids)
        input = open(file,'r')
        i=0
        for line in input:
            self.obsmembids[i]=[]
            temp=line[1:(len(line)-2)]
            z=temp.split(',')
            if len(z) == 1:
                i=i+1
                continue
            id=[]
            #print "len(x) = ",len(z)
            for mid in z:
                #print mid#,self.obsmembids[i]
                id.append(int(mid))
                #self.obsmembids[i].append(int(mid))
            #print id
            self.obsmembids[i]=id
            i=i+1
        input.close()

class Galaxy:
    def __init__(self):
        self.igax = []#galaxy id
        self.igh = [] #galaxy id in that halo
        self.ih = [] #halo id
        self.type=[] #type 1/late 2/early
        self.x1 = []
        self.x2  = []
        self.x3  = []
        self.v1 = []
        self.v2  = []
        self.v3  = []
        self.lbj = [] #luminosity of galaxy (log10 (L/h^{-2}L_{\odot}) in b_J band
        self.MR = [] #absolute R-band luminosity
        self.hmass=[] #halo mass (log10 (M/h^{-1}M_{\odot})
        self.hsigmax=[] # line-of-sigt velocity dispersion (log10 (v / km/s))
        self.hrvir = []
        self.hsigma = []
    def readfiles(self,galcat):

        infile=open(galcat,'r')
        i=0
        for line in infile:
            if line.find('#') > -1:
                continue
            fields=line.split()
            if len(fields) < 4:#skip line with 3 entries
                continue
            if float(fields[10]) > lbjmin:
                self.igax.append(float(fields[0]))
                self.igh.append(float(fields[1]))
                self.ih.append(float(fields[2]))
                self.type.append(float(fields[3]))
                self.x1.append(float(fields[4]))
                self.x2.append(float(fields[5]))
                self.x3.append(float(fields[6]))
                self.v1.append(float(fields[7]))
                self.v2.append(float(fields[8]))
                self.v3.append(float(fields[9]))
                self.lbj.append(float(fields[10]))
                self.hmass.append(float(fields[11]))
                self.hsigmax.append(10.**float(fields[12]))
                #self.hsigma.append(float(c.sigma[int(fields[2])]))
                #self.hrvir.append(float(c.rvir[int(fields[2])]))
            i=i+1


    def readGIFfiles(self,galcat,cmass,csigma,crvir):
        self.bj=[]
        self.vbulge=[]
        self.mstar=[]
        self.mgas=[]
        self.sfr=[]

        infile=open(galcat,'r')
        i=0
        j=0
        for line in infile:
            if line.find('#') > -1:
                continue
            fields=line.split()
            if len(fields) < 4:#skip line with 3 entries
                haloindex=int(fields[0])#halo index
                continue
            if len(fields) > 4:
                if j == 0:
                    go=0
                    if float(fields[5]) < bJmax:
                        self.bj.append(float(fields[5]))
                        self.igh.append(float(fields[1]))
                        self.ih.append(haloindex)
                        t=(Mbsol-float(fields[5]))/2.5
                        #llbjmin=N.log10(10.**((Mbsol-bbJmax)/2.5))
                        self.lbj.append(t)#log10 L(B_J)
                        self.MR.append(float(fields[7]))#log10 L(B_J)
                        self.hmass.append(cmass[haloindex])
                        self.hsigma.append(csigma[haloindex])
                        self.hsigmax.append(csigma[haloindex]/N.sqrt(3.))
                        self.hrvir.append(crvir[haloindex])
                        #print haloindex,float(fields[0])
                        go=1
                    j+=1
                    continue
                if j == 1:
                    j=0
                    if go > 0:
                        self.vbulge.append(float(fields[0]))
                        self.mstar.append(float(fields[1]))
                        self.mgas.append(float(fields[2]))
                        self.sfr.append(float(fields[3]))
                        self.igax.append(float(fields[4]))
                        self.x1.append(float(fields[5]))
                        self.x2.append(float(fields[6]))
                        self.x3.append(float(fields[7]))
                        self.v1.append(float(fields[8]))
                        self.v2.append(float(fields[9]))
                        self.v3.append(float(fields[10]))
                        if float(fields[3]) > 0:#star-forming
                            self.type.append(1.)#late type, to match Mo catalogs
                        else:
                            self.type.append(2.)#early type
                        go = 0
        self.sfr=N.array(self.sfr,'f')
        self.MR=N.array(self.MR,'f')
        self.ew=self.sfr*(3.3e7)*(10.**((self.MR-4.76)/2.5))
        print "length of sfr and hrvir = ",len(self.sfr),len(self.hrvir)
    def calcdr(self,c):
        self.x1=N.array(self.x1,'f')
        self.x2=N.array(self.x2,'f')
        self.dr=N.zeros(len(self.x1),'f')
        self.drRvir=N.zeros(len(self.x1),'f')
        for i in range(len(self.dr)):
            ci=(self.ih[i]-1)
            self.dr[i]=N.sqrt((self.x1[i]-c.x1[ci])**2 + (self.x2[i]-c.x2[ci])**2)
            self.drRvir[i]=N.sqrt((self.x1[i]-c.x1[ci])**2 + (self.x2[i]-c.x2[ci])**2)/self.hrvir[i]
    def doggif(self):
        self.convarray()
        #print "Number of z=0 clusters before mass cut = ",len(cgif0.mass)
        #if zbox > .05:
        #    v=3.e5*((1+zbox)**2-1)/((1+zbox)**2 + 1)
        #    self.v3=(self.v3+v)/(1+self.v3*v/(3.e5)**2)#relativistic addition of velocities

        self.getvobs()

    def convarray(self):
        self.igax = N.array(self.igax,'f')
        self.igh = N.array(self.igh,'f')
        self.ih = N.array(self.ih,'f')
        self.type = N.array(self.type,'f')
        self.type = abs(self.type - 2)#change type so that 0/early,1/late
        self.x1 = N.array(self.x1,'f')
        self.x2 = N.array(self.x2,'f')
        self.x3 = N.array(self.x3,'f')
        self.v1 = N.array(self.v1,'f')
        self.v2 = N.array(self.v2,'f')
        self.v3 = N.array(self.v3,'f')
        self.lbj = N.array(self.lbj,'f')
        self.hmass = N.array(self.hmass,'f')
        self.hsigmax = N.array(self.hsigmax,'f')
        self.hrvir = N.array(self.hrvir,'f')
        self.oigax = copy.copy(self.igax)#save tuple of original arrays
        self.oigh = copy.copy(self.igh)
        self.oih = copy.copy(self.ih)
        self.otype = copy.copy(self.type)
        self.ox1 = copy.copy(self.x1)
        self.ox2 = copy.copy(self.x2)
        self.ox3 = copy.copy(self.x3)
        self.ov1 = copy.copy(self.v1)
        self.ov2 = copy.copy(self.v2)
        self.ov3 = copy.copy(self.v3)
        self.olbj = copy.copy(self.lbj)
        self.ohmass = copy.copy(self.hmass)
        self.ohsigmax = copy.copy(self.hsigmax)
        self.ohrvir = copy.copy(self.hrvir)
                
	if (gif > 0):
	    self.odr = copy.copy(self.dr)
	    self.odrRvir = copy.copy(self.drRvir)

	    self.osfr = copy.copy(self.sfr)
	    self.oew = copy.copy(self.ew)
	    self.odr = copy.copy(self.dr)
	    self.odrRvir = copy.copy(self.drRvir)
	    self.hsigma = N.array(self.hsigma,'f')
	    self.hrvir = N.array(self.hrvir,'f')
                
        #self.hsigma = N.array(self.hsigma,'f')
        #self.hrvir = N.array(self.hrvir,'f')
        
    def getvobs(self):
        self.vobs=self.x3*H0+self.v3#Hubble flow
        #should evolve H0
        Hz=H0
        v=self.x3*Hz
        self.vobs=(v+self.v3)/(1+v*self.v3/(3.e5)**2)#Hubble flow
        self.ovobs=copy.copy(self.vobs)
        #self.z=N.sqrt((1+self.vobs/3.e5)/(1+self.vobs/3.e5)) - 1.#observed redshift
        self.z=N.sqrt((1+self.vobs/3.e5)/(1-self.vobs/3.e5))-1.#observed redshift
        self.oz=copy.copy(self.z)
    def cutonlbj(self,bbJmax):
        llbjmin=((Mbsol-bbJmax)/2.5)
        #print "within cutonlbj, Mbsol,bbJmax,llbjmin",Mbsol,bbJmax,llbjmin,len(self.oih)
        #for i in range(20):
        #    m=Mbsol-self.olbj[i]*2.5
        #    print "cut on lbj ",i,llbjmin,self.olbj[i],bbJmax,m
        self.igax = N.compress(self.olbj>llbjmin,self.oigax)
        self.igh = N.compress(self.olbj>llbjmin,self.oigh)
        self.ih = N.compress(self.olbj>llbjmin,self.oih)
        self.type =N.compress(self.olbj>llbjmin,self.otype)
        self.x1 = N.compress(self.olbj>llbjmin,self.ox1)
        self.x2 = N.compress(self.olbj>llbjmin,self.ox2)
        self.x3 = N.compress(self.olbj>llbjmin,self.ox3)
        self.v1 = N.compress(self.olbj>llbjmin,self.ov1)
        self.v2 = N.compress(self.olbj>llbjmin,self.ov2)
        self.v3 = N.compress(self.olbj>llbjmin,self.ov3)
        self.vobs = N.compress(self.olbj>llbjmin,self.ovobs)
        self.z = N.compress(self.olbj>llbjmin,self.oz)
        self.hmass = N.compress(self.olbj>llbjmin,self.ohmass)
        self.lbj = N.compress(self.olbj>llbjmin,self.olbj)
	if (gif > 0):
	    self.hsigmax = N.compress(self.olbj>llbjmin,self.ohsigmax)
	    self.hrvir = N.compress(self.olbj>llbjmin,self.ohrvir)

	    self.sfr = N.compress(self.olbj>llbjmin,self.osfr)
	    self.ew = N.compress(self.olbj>llbjmin,self.oew)
	    self.dr = N.compress(self.olbj>llbjmin,self.odr)
	    self.drRvir = N.compress(self.olbj>llbjmin,self.odrRvir)
#read in galaxy catalog, keeping those with bJ < bJmax

def writefiles():
#for each cluster, find indices of all member galaxies in galaxy array
    output=open("myclusters.cat",'w')
    output2=open("pca.dat",'w')
    header="# "+str(len(c.z))+ " clusters \n"
    output.write(header)
    header="#id z rvir r200 sigma totrlum lx stellmass05 stellmass1 stellmass2 stellmass3 sfr05 sfr1 sfr2 sfr3 sffrac05 sffrac1 sffrac2 sub r200bin sigmabin sfr05bin sfr1bin sfr2bin sfr3bin stellmass05bin stellmass1bin stellmass2bin stellmass3bin chisq(diff from gauss) sfr05hiz kauffmann ra dec ngal05 ngal1 ngal2 compl05 compl1 compl2 voffset"
    output.write("%s \n" % (header))
    for i in range(len(c.z)):
        #print i,c.id[i],c.z[i],c.rvir[i],c.r200[i],c.sigma[i],c.Lx[i],c.xrayflag[i],c.stellmass05[i],c.stellmass1[i],c.stellmass2[i],c.stellmass3[i],c.sfr05[i],c.sfr1[i],c.sfr2[i],c.sfr3[i],c.sffrac05[i],c.sffrac1[i],c.sffrac2[i]
        #             id    z     rvir  r200  sigma totrl lx    xflg  ste05 stel1 stel2 stel3 sfr05 sfr1  sfr2  sfr3  sff05 sff1  sff2  sub   totl
        output.write('%5.0f %5.4f %5.3f %5.3f %6.1f %6.1f %6.3f %2.0f %3.2e %3.2e %3.2e %3.2e %5.1f %5.1f %5.1f %5.1f %4.2f %4.2f %4.2f %2.0f %4.3e %4.3e %4.3e %5.3f %5.3f %5.1f %5.1f %5.1f %5.1f %3.2e %3.2e %3.2e %3.2e %5.2f %3.2e %5.3f %12.8f %12.8f %6.1f  %6.1f  %6.1f  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n' %(c.id[i],c.z[i],c.rvir[i],c.r200[i],c.sigma[i],c.totrlum[i],c.Lx[i],c.xrayflag[i],c.stellmass05[i],c.stellmass1[i],c.stellmass2[i],c.stellmass3[i],c.sfr05[i],c.sfr1[i],c.sfr2[i],c.sfr3[i],c.sffrac05[i],c.sffrac1[i],c.sffrac2[i],c.sub[i],c.totlum05[i],c.totlum1[i],c.totlum2[i],c.r200bin[i],c.sigmabin[i],c.sfr05bin[i],c.sfr1bin[i],c.sfr2bin[i],c.sfr3bin[i],c.stellmass05bin[i],c.stellmass1bin[i],c.stellmass2bin[i],c.stellmass3bin[i],c.chisq[i],c.sfr05hiz[i],c.kauffmann[i],c.ra[i],c.dec[i],c.ngal05[i],c.ngal1[i],c.ngal2[i],c.compl05[i],c.compl1[i],c.compl2[i],c.voffset[i],c.nsf05hiz[i],c.n05hiz[i]))
	output2.write('%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n' %(c.rvir[i],c.r200[i],c.sigma[i],c.totrlum[i],c.stellmass1[i],c.sfr1[i],c.sffrac1[i]))
    output.close()
    output2.close()

def readfiles():
    output=open("myclusters.cat",'w')
    output2=open("pca.dat",'w')
    header="# "+str(len(c.z))+ " clusters \n"
    output.write(header)
    header="#id z rvir r200 sigma totrlum lx stellmass05 stellmass1 stellmass2 stellmass3 sfr05 sfr1 sfr2 sfr3 sffrac05 sffrac1 sffrac2 sub r200bin sigmabin sfr05bin sfr1bin sfr2bin sfr3bin stellmass05bin stellmass1bin stellmass2bin stellmass3bin chisq(diff from gauss) sfr05hiz kauffmann ra dec ngal05 ngal1 ngal2 compl05 compl1 compl2 voffset"
    output.write("%s \n" % (header))
    for i in range(len(c.z)):
        #print i,c.id[i],c.z[i],c.rvir[i],c.r200[i],c.sigma[i],c.Lx[i],c.xrayflag[i],c.stellmass05[i],c.stellmass1[i],c.stellmass2[i],c.stellmass3[i],c.sfr05[i],c.sfr1[i],c.sfr2[i],c.sfr3[i],c.sffrac05[i],c.sffrac1[i],c.sffrac2[i]
        #             id    z     rvir  r200  sigma totrl lx    xflg  ste05 stel1 stel2 stel3 sfr05 sfr1  sfr2  sfr3  sff05 sff1  sff2  sub   totl
        output.write('%5.0f %5.4f %5.3f %5.3f %6.1f %6.1f %6.3f %2.0f %3.2e %3.2e %3.2e %3.2e %5.1f %5.1f %5.1f %5.1f %4.2f %4.2f %4.2f %2.0f %4.3e %4.3e %4.3e %5.3f %5.3f %5.1f %5.1f %5.1f %5.1f %3.2e %3.2e %3.2e %3.2e %5.2f %3.2e %5.3f %12.8f %12.8f %6.1f  %6.1f  %6.1f  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n' %(c.id[i],c.z[i],c.rvir[i],c.r200[i],c.sigma[i],c.totrlum[i],c.Lx[i],c.xrayflag[i],c.stellmass05[i],c.stellmass1[i],c.stellmass2[i],c.stellmass3[i],c.sfr05[i],c.sfr1[i],c.sfr2[i],c.sfr3[i],c.sffrac05[i],c.sffrac1[i],c.sffrac2[i],c.sub[i],c.totlum05[i],c.totlum1[i],c.totlum2[i],c.r200bin[i],c.sigmabin[i],c.sfr05bin[i],c.sfr1bin[i],c.sfr2bin[i],c.sfr3bin[i],c.stellmass05bin[i],c.stellmass1bin[i],c.stellmass2bin[i],c.stellmass3bin[i],c.chisq[i],c.sfr05hiz[i],c.kauffmann[i],c.ra[i],c.dec[i],c.ngal05[i],c.ngal1[i],c.ngal2[i],c.compl05[i],c.compl1[i],c.compl2[i],c.voffset[i],c.nsf05hiz[i],c.n05hiz[i]))
	output2.write('%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n' %(c.rvir[i],c.r200[i],c.sigma[i],c.totrlum[i],c.stellmass1[i],c.sfr1[i],c.sffrac1[i]))
    output.close()
    output2.close()

#Main


sim=3
if sim == 1:
    halocat='/Users/rfinn/mock-clusters/halo/gcatb02.8211'
    galcat='/Users/rfinn/mock-clusters/mock/m8211'
    simL=100.#100 Mpc cube
    nbin=5
    magcut=[-16.,-17.,-18.,-19.]
if sim == 3:
    halocat='/Users/rfinn/mock-clusters/halo/gcatb02.8110'
    galcat='/Users/rfinn/mock-clusters/mock/m8110'
    simL=300.#100 Mpc cube
    nbin=7
    magcut=[-16.,-17.,-18.,-19.]
#sim=3
#initialize
df=1
gif=0
if (df == 1):
    zbox=0.
    c=Cluster()
    print "Reading halo file"
    c.readfiles(halocat)
    print "converting to arrays"
    c.convarray()
    c.getvobs()
    g=Galaxy()
    print "Reading galaxy file"
    g.readfiles(galcat)
    print "converting to arrays"
    g.convarray()
    g.getvobs()

    #nv=3.
    #nr=2.
    #bJmax=-16.

    nv=3.
    nr=1.
    bjmin=-17.

    #g.cutonlbj(bJmax)
    #c.measurecontam(nv,nr,g)
    print "before get all members, len(g) = ",len(g.x1)
    c.getallmembers(g,nv,nr)


    #plotmagcuts()
    #print "PLOTTING SIGMA CUTS"
    #plotsigcuts()
    #plotsigcutssdss()


    #print "PLOTTING RADIAL CUTS"
    #plotradcuts()

    #print "PLOTTING NGAL VS Mcl"
    #plotngalmclradcuts()

    #print "PLOTTING TYPE VS MAG CUTS"
    #plottypemagcuts()

    #print "PLOTTING CUMULATIVE VS MAG"
    #plotcummag()

    #print "PLOTTING CUMULATIVE VS dv"
    #plotcumdv()

    print "PLOTTING HISTOGRAM of dv"
    plothistdv()

    print "PLOTTING HISTOGRAM of dr"
    plothistdr()
    #print "PLOTTING dV/dv/dOmega vs z"
    #plotdVdz()
    #print "PLOTTING dL vs z"
    #plotdLdz()
    #plotsigcutsdr()
    
endtime=time.clock()
print "end time = ",endtime
print "elapsed time = ",endtime-starttime

if (gif > 0):
    print "Starting GIF plots..."
    simL=141.
    nr=2.
    nv=3.
    nbin=4
    bJmax=-19.
    Mclmin=100.#min cluster mass in units of 10^12 Msun
    getmemb=0 #set to 1 to calculate observed members, 0 to read from files
    zbox=0.
    print "zbox = ",zbox
    halocat="/Users/rfinn/mock-clusters/GIF/halos/halos.propl_1178"
    galcat="/Users/rfinn/mock-clusters/GIF/galaxies/lcdm_galaxy_cat.z0.00"
    cgif0=Cluster()
    cgif0.dogif()


    ggif0=Galaxy()
    ggif0.readGIFfiles(galcat,cgif0.mass,cgif0.sigma,cgif0.rvir)
    ggif0.calcdr(cgif0)
    ggif0.doggif()
    print "len before cut mass = ",len(cgif0.mass)
    cgif0.cutmass()
    print "len after cut mass = ",len(cgif0.mass)
    #print "for z=0, self.x2[0:10] = ",cgif0.x2[0:10]
    #cgif0.measurecontam(nv,nr,ggif0)
    print "number of clusters with M > ",Mclmin," = ",len(cgif0.mass)
    zbox=0.2
    print "zbox = ",zbox
    halocat="/Users/rfinn/mock-clusters/GIF/halos/halos.propl_0938"
    galcat="/Users/rfinn/mock-clusters/GIF/galaxies/lcdm_galaxy_cat.z0.20"
    cgif02=Cluster()
    cgif02.dogif()


    ggif02=Galaxy()
    ggif02.readGIFfiles(galcat,cgif02.mass,cgif02.sigma,cgif02.rvir)
    ggif02.calcdr(cgif02)
    ggif02.doggif()
    cgif02.cutmass()
    print "number of clusters with M > ",Mclmin," = ",len(cgif02.mass)
    #cgif02.measurecontam(nv,nr,ggif02)

    zbox=0.42
    print "zbox = ",zbox
    halocat="/Users/rfinn/mock-clusters/GIF/halos/halos.propl_0731"
    galcat="/Users/rfinn/mock-clusters/GIF/galaxies/lcdm_galaxy_cat.z0.42"
    cgif04=Cluster()
    cgif04.dogif()

    ggif04=Galaxy()
    ggif04.readGIFfiles(galcat,cgif04.mass,cgif04.sigma,cgif04.rvir)
    ggif04.calcdr(cgif04)
    ggif04.doggif()
    cgif04.cutmass()
    print "number of clusters with M > ",Mclmin," = ",len(cgif04.mass)

    zbox=0.62
    print "zbox = ",zbox
    halocat="/Users/rfinn/mock-clusters/GIF/halos/halos.propl_0612"
    galcat="/Users/rfinn/mock-clusters/GIF/galaxies/lcdm_galaxy_cat.z0.62"
    cgif06=Cluster()
    cgif06.dogif()

    ggif06=Galaxy()
    ggif06.readGIFfiles(galcat,cgif06.mass,cgif06.sigma,cgif06.rvir)
    ggif06.calcdr(cgif06)
    ggif06.doggif()
    cgif06.cutmass()

    print "number of clusters with M > ",Mclmin," = ",len(cgif06.mass)
    zbox=0.82
    print "zbox = ",zbox
    halocat="/Users/rfinn/mock-clusters/GIF/halos/halos.propl_0519"
    galcat="/Users/rfinn/mock-clusters/GIF/galaxies/lcdm_galaxy_cat.z0.82"
    cgif08=Cluster()
    cgif08.dogif()

    ggif08=Galaxy()
    ggif08.readGIFfiles(galcat,cgif08.mass,cgif08.sigma,cgif08.rvir)
    ggif08.calcdr(cgif08)
    ggif08.doggif()
    cgif08.cutmass()

    print "number of clusters with M > ",Mclmin," = ",len(cgif08.mass)
    zbox=1.05
    print "zbox = ",zbox
    halocat="/Users/rfinn/mock-clusters/GIF/halos/halos.propl_0439"
    galcat="/Users/rfinn/mock-clusters/GIF/galaxies/lcdm_galaxy_cat.z1.05"
    cgif1=Cluster()
    cgif1.dogif()

    ggif1=Galaxy()
    ggif1.readGIFfiles(galcat,cgif1.mass,cgif1.sigma,cgif1.rvir)
    ggif1.calcdr(cgif1)
    ggif1.doggif()
    cgif1.cutmass()

    if getmemb > 0:
        z=['0','02','04','06','08','1']
        print "getting observed members, once and for all"

        zi=N.array([0,1,2,3,4,5],'i')
        #zi=N.array([],'i')
        #for i in range(len(z)):            
        #for i in range(1):
        for i in zi:
            file='obsmemb-cgif'+str(z[i])
            output=open(file,'w')
            if i == 0:
                ggif0.cutonlbj(bJmax)
                cgif0.getobsmemberids(nv,nr,ggif0)#get galaxies within nr*dr and nv*dv
                c=cgif0
                g=ggif0
            if i == 1:
                ggif02.cutonlbj(bJmax)
                cgif02.getobsmemberids(nv,nr,ggif02)#get galaxies within nr*dr and nv*dv
                c=cgif02
                g=ggif02
            if i == 2:
                ggif04.cutonlbj(bJmax)
                cgif04.getobsmemberids(nv,nr,ggif04)#get galaxies within nr*dr and nv*dv
                c=cgif04
                g=ggif04
            if i == 3:
                ggif06.cutonlbj(bJmax)
                cgif06.getobsmemberids(nv,nr,ggif06)#get galaxies within nr*dr and nv*dv
                c=cgif06
                g=ggif06
            if i == 4:
                ggif08.cutonlbj(bJmax)
                cgif08.getobsmemberids(nv,nr,ggif08)#get galaxies within nr*dr and nv*dv
                c=cgif08
                g=ggif08
            if i == 5:
                ggif1.cutonlbj(bJmax)
                cgif1.getobsmemberids(nv,nr,ggif1)#get galaxies within nr*dr and nv*dv
                c=cgif1
                g=ggif1
            for j in range(len(c.obsmembids)):
                out=str(c.obsmembids[j])+'\n'
                output.write(out)
            output.close()
    if getmemb < 1:
        print "reading observed members z=0 from file"
        cgif0.readobsmembids('obsmemb-cgif0')
        print "reading observed members z=0.2 from file"
        cgif02.readobsmembids('obsmemb-cgif02')
        print "reading observed members z=0.4 from file"
        cgif04.readobsmembids('obsmemb-cgif04')
        print "reading observed members z=0.6 from file"
        cgif06.readobsmembids('obsmemb-cgif06')
        print "reading observed members z=0.8 from file"
        cgif08.readobsmembids('obsmemb-cgif08')
        print "reading observed members z=1. from file"
        cgif1.readobsmembids('obsmemb-cgif1')
    #print "Number of z=0 clusters = ",len(cgif0.mass)
    endtime=time.clock()
    print "elapsed time = ",endtime-starttime
    #plotmagcutsgif()
    plotradcutsgif()
    #plottotalsfrgif()
    #plotsffdr()
