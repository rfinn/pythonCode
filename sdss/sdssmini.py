#!/usr/bin/env python
"""
mini routine to fit extinction vs Mr and stellar mass vs Mr

"""
import sys, glob
#import Numeric as N
import numarray as N
import scipy
import pylab #use matplotlib instead of scipy
from math import *
print "loading ppgplot"
import ppgplot
print "got here"
import time, os
import mystuff as my
import matplotlib

def fitting_func(x,a,b):
    return a*x**b

def objective_func(parameters,x,y):
    a,b = parameters
    sum_of_squares = 0.
    data = zip(x,y)
    for x,y in data:
	f_of_x = fitting_func(x,a,b)
	dy=f_of_x - y
	sum_of_squares += dy*dy
    return sum_of_squares

def leastsq_func(parameters,x,y):
    a,b = parameters
    f_of_x = fitting_func(x,a,b)
    dy=f_of_x - y
	
    return dy

def fitpowerlaw(x,y):
    p_init = scipy.array([0.,1.])
    p_opt = scipy.optimize.fmin(objective_func,p_init,(x,y),full_output=1,maxiter=1000)
    print p_opt
    t=p_opt[0]
    aopt=t[0]
    bopt=t[1]
    return aopt,bopt


class galaxy:
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
                self.ar=N.zeros(ncl,'f')
                self.mpaflag=N.zeros(ncl,'f')
                i += 1
                continue
            if line.find('#') > -1:#skip header
                continue
            f=line.split()
            for j in range(len(f)):
                f[j]=float(f[j])
            #(self.clusterid[i],self.clusterz[i],self.clusterrvir[i],self.clusterr200[i],self.clustersigma[i],self.dr[i],self.dv[i],self.ew[i],self.lHa[i],self.sfr[i],self.stellarmass[i],self.agn[i],self.memb[i],self.Mabs[i],self.ra[i],self.dec[i],self.sigma10[i])=(f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16])
            #print i, len(f)
            (self.clusterid[i],self.clusterz[i],self.clusterrvir[i],self.clusterr200[i],self.clustersigma[i],self.dr[i],self.dv[i],self.ew[i],self.lHa[i],self.sfr[i],self.stellarmass[i],self.agn[i],self.memb[i],self.Mabs[i],self.ra[i],self.dec[i],self.ar[i],self.mpaflag[i])=(f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17])
            
            i += 1
        output.close()
    def azmr(self):
	x=N.compress((self.mpaflag > 0.1) & (self.ew > 4.) & (self.Mabs < -18.),self.Mabs)
	y=N.compress((self.mpaflag > 0.1) & (self.ew > 4.) & (self.Mabs < -18.),self.ar)
	x1=N.compress((self.mpaflag > 0.1) & (self.ew > 4.) & (self.Mabs < -20.38),self.Mabs)
	y1=N.compress((self.mpaflag > 0.1) & (self.ew > 4.) & (self.Mabs < -20.38),self.ar)
	y=2.5*N.log10(y)
	#pylab.plot(x,y,'k.',markersize=0.1,zorder=1)

	print "average Ar for Mr < -20.38 = %5.2f +/- %5.2f"%(N.average(y1),pylab.std(y1))
	(xbin,ybin)=my.binit(x1,y1,20)
	#(xbin,ybin,ybinerr)=my.biniterr(x,y,20)
	for i in range(len(xbin)):
	    print i,xbin[i],ybin[i]
	print "Average of binned values = ",N.average(ybin)
	print "average Ar for Mr < -20.38 = %5.2f +/- %5.2f"%(N.average(N.log10(y1)),pylab.std(N.log10(y1)))
	#pylab.axis([-26.,-12.,0.1,30.])
	pylab.xlabel(r'$\rm{M_r}$',fontsize=28.)
	pylab.ylabel(r'$\rm{A_r}$',fontsize=28.)
	(xbin,ybin)=my.binit(x,y,20)
	#(xbin,ybin,ybinerr)=my.biniterr(x,y,20)
	for i in range(len(xbin)):
	    print i,xbin[i],ybin[i]

	pylab.plot(xbin,ybin,'r-',lw=5)
	ax=pylab.gca()
	xmin=-24.
	xmax=-18.
	ymin=-1.
	ymax=3.
	my.contourf(x,y,xmin,xmax,ymin,ymax)
	pylab.axvline(x=-20.6,linewidth=3,ls='--',c='g')
	xl=N.arange(-23.,-20.5,.2)
	yl=0.76*N.ones(len(xl),'f')
	pylab.plot(xl,yl,'b-',lw=3)


	pylab.axis([-24.,-18,-1.,2.4])
	#ax.set_yscale('log')
	#pylab.show()
	pylab.savefig('armr.eps')
	print "fraction w/MPA stellar mass and Az = ",N.sum(self.mpaflag)/(1.*len(self.mpaflag))

    def azhist(self):
	x=N.compress((self.mpaflag > 0.1) & (self.ew > 4.),self.Mabs)
	y=N.compress((self.mpaflag > 0.1) & (self.ew > 4.),self.ar)
	y=2.5*N.log10(y)
	m1=-18.
	self.azhistsub(x,y,m1,'b')
	#m1=-19.
	#self.azhistsub(x,y,m1,'r')
	m1=-20.
	self.azhistsub(x,y,m1,'g')
	#m1=-21.
	#self.azhistsub(x,y,m1,'k')
	m1=-22.
	self.azhistsub(x,y,m1,'y')
	#m1=-23.
	#self.azhistsub(x,y,m1,'c')
	#pylab.show()
	#ax=pylab.gca()
	#ax.set_xscale('log')
	pylab.legend(loc='upper right')
	pylab.savefig('azhist.eps')
    def azhistsub(self,x,y,m1,color):
	m2=m1-2.
	y1=N.compress((x< m1) & (x > m2),y)
	(n,bin,patches)=pylab.hist(y1,bins=10,normed=1,fill=False,ec=color,fc='w',alpha=1.,visible=False)
	s=color+'-'
	s1=str(m1)+'> M > '+str(m2)
	pylab.plot(bin,n,s,label=s1)

    def stellarmr(self):
	pylab.cla()
	x=N.compress((self.mpaflag > 0.1),self.Mabs)
	y=N.compress((self.mpaflag > 0.1),self.stellarmass)
	print "len x,y,Mabs = ",len(x),len(y),len(self.Mabs)
	#x=[]
	#y=[]
	#for i in range(len(self.Mabs)):
	#    if self.mpaflag[i] > 0.1:
	#	x.append(self.Mabs[i])
	#	y.append(self.stellarmass[i])
	#print "len x,y,Mabs =",len(x),len(y),len(self.Mabs)
	#x=N.array(x,'d')
	#y=N.array(y,'d')
	y=N.log10(y/1.e11)
	xmin=-24.
	xmax=-18.
	ymin=-3.
	ymax=2.

	#my.contour(x,y,xmin,xmax,ymin,ymax)
	#pylab.plot(x,y,'k.',markersize=.01,zorder=1)

	pylab.xlabel(r'$\rm{M_r}$',fontsize=28.)
	pylab.ylabel(r'$\rm{log_{10}(M_* / 10^{11} \ M_\odot)}$',fontsize=28.)
	#(a,b)=fitpowerlaw(x,y)

	xtrans=-21.5
	b=-13.8                                          
	m=-.63
	xl=N.arange(xtrans,-18.,.05)
	xl=N.array(xl,'f')
	yl=N.zeros(len(xl),'f')
	yl=m*(xl)+b
	#xl=-1.*xl
	#print xl
	#print yl
	pylab.plot(xl,yl,'b-',lw=4,label='_nolegend_')


	m2=-.47
	b2=(m-m2)*xtrans+b
	print "b2 = ",b2
	xl=N.arange(-24.,(xtrans+.05),.05)
	xl=N.array(xl,'f')
	yl=N.zeros(len(xl),'f')
	yl=m2*(xl)+b2
	#xl=-1.*xl
	#print xl
	#print yl
	pylab.plot(xl,yl,'c-',lw=4,label='_nolegend_')
	pylab.axvline(x=-20.6,linewidth=3,ls='--',c='g')
	#(xbin,ybin,ybinerr)=my.binitbinsfix(xl,(xl+0.5),x,y)
	#pylab.plot(xbin,ybin,'b-',linewidth=4.)
	pylab.hold()
	my.contourf(x,y,xmin,xmax,ymin,ymax)
	ax=pylab.gca()
	#ax.set_yscale('log')

	pylab.axis([-24.5,-17.5,-3.,1.2])
	pylab.savefig('stellarmr.eps')
	#pylab.show()

g=galaxy()
g.readdatafile("/Users/rfinn/SDSS/baloghDR5/mygalaxies.cat")


g.azmr()
#g.stellarmr()
#g.azhist()
#f.readdatafile("/Users/rfinn/SDSS/fieldDR5/mygalaxies.cat")


