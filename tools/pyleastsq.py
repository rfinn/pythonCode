#!/usr/bin/env python
"""
useage

pyleastsq.py dr clust

where dr = 4 for DR4
           5 for DR5

and clust = 1 for clusters
            0 for field
"""

import scipy
import scipy.optimize
import sys

dr = int(sys.argv[1])

class cluster:
    def readtable(self,file):
	infile = open(file,'r')	
	#if abs(dr - 4.) < .1:
	#infile = open('/Users/rfinn/SDSS/fieldDR5/myclusters.cat','r')
	#print "reading fieldDR5"
	#if abs(dr - 5.) < .1:
	#infile = open('/Users/rfinn/SDSS/baloghDR5/myclusters.cat','r')
	#print "reading baloghDR5"
	sigma=[]
	ngal=[]
	mstar=[]
	sfr=[]
	ngal05=[]
	mstar05=[]
	sfr05=[]
	ngal2=[]
	mstar2=[]
	sfr2=[]
	sffrac05=[]
	sffrac=[]
	sffrac2=[]
	for line in infile:
	    if line.find('#') > -1:
		continue
	    t=line.split()
	    sigma.append(float(t[4]))
	    mstar.append(float(t[9]))
	    mstar05.append(float(t[8]))
	    mstar2.append(float(t[10]))
	    ngal.append(float(t[39]))
	    ngal05.append(float(t[38]))
	    ngal2.append(float(t[40]))
	    sfr.append(float(t[13]))
	    sfr05.append(float(t[12]))
	    sfr2.append(float(t[14]))
	    sffrac05.append(float(t[16]))
	    sffrac.append(float(t[17]))
	    sffrac2.append(float(t[18]))
	self.sigma=scipy.array(sigma,'f')
	self.mstar=scipy.array(mstar,'f')/1.e11
	self.ngal=scipy.array(ngal,'f')
	self.sfr=scipy.array(sfr,'f')
	self.mstar05=scipy.array(mstar05,'f')/1.e11
	self.ngal05=scipy.array(ngal05,'f')
	self.sfr05=scipy.array(sfr05,'f')
	self.mstar2=scipy.array(mstar2,'f')/1.e11
	self.ngal2=scipy.array(ngal2,'f')
	self.sfr2=scipy.array(sfr2,'f')
	self.sffrac05=scipy.array(sffrac05,'f')
	self.sffrac=scipy.array(sffrac,'f')
	self.sffrac2=scipy.array(sffrac2,'f')
	#return sigma,ngal05,ngal,ngal2,mstar05,mstar,mstar2,sfr05,sfr,sfr2
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


def bootstrap(x,y):
    avea=[]
    aveb=[]
    for i in range(100):
	(xrand,yrand)=getrand(x,y)
	p_opt = scipy.optimize.fmin(objective_func,p_init,(xrand,yrand),full_output=0,disp=0,maxiter=1000)
	#print p_opt
	t=p_opt
	avea.append(t[0])
	aveb.append(t[1])

    print "bootstrap intercept ave/std = %1.4e %1.4e"%(scipy.mean(avea), scipy.std(avea))
    print "bootstrap power ave/std = %1.4e %1.4e"%(scipy.mean(aveb), scipy.std(aveb))
    bsort=scipy.sort(aveb)
    #print bsort
    print "2.5 percent value = ",bsort[int(len(bsort)*.025)]
    print "97.5 percent value = ",bsort[int(len(bsort)*.975)]
    return scipy.mean(avea),scipy.std(avea),scipy.mean(aveb),scipy.std(aveb)
def getrand(x,y):
    xrand=scipy.zeros(len(x),'f')
    yrand=scipy.zeros(len(x),'f')
    for i in range(len(x)):
	j=int(len(x)*scipy.rand())
	xrand[i]=x[j]
	yrand[i]=y[j]
    return xrand,yrand

def fitting_func3(x,a):
    #return a*x**2.35

    return a*x**2.97

def objective_func3(parameters,x,y):
    a = parameters
    sum_of_squares = 0.
    data = zip(x,y)
    for x,y in data:
	f_of_x = fitting_func3(x,a)
	dy=f_of_x - y
	sum_of_squares += dy*dy
    return sum_of_squares



def fitpowerlaw3(x,y):
    p_init = scipy.array([0.])
    p_opt = scipy.optimize.fmin(objective_func3,p_init,(x,y),full_output=1,maxiter=1000)
    print p_opt
    t=p_opt[0]
    aopt=t[0]
    #bopt=t[1]
    return aopt#,bopt

def bootstrap3(x,y):
    avea=[]
    #aveb=[]
    p_init = scipy.array([0.])
    for i in range(100):
	(xrand,yrand)=getrand(x,y)
	p_opt = scipy.optimize.fmin(objective_func3,p_init,(xrand,yrand),full_output=0,disp=0,maxiter=1000)
	#print p_opt
	t=p_opt
	avea.append(t[0])
	#aveb.append(t[1])

    print "bootstrap intercept ave/std = %1.4e %1.4e"%(scipy.mean(avea), scipy.std(avea))
    #print "bootstrap power ave/std = %1.4e %1.4e"%(scipy.mean(aveb), scipy.std(aveb))
    #bsort=scipy.sort(aveb)
    #print bsort
    #print "2.5 percent value = ",bsort[int(len(bsort)*.025)]
    #print "97.5 percent value = ",bsort[int(len(bsort)*.975)]
    return scipy.mean(avea),scipy.std(avea)#,scipy.mean(aveb),scipy.std(aveb)


p_init = scipy.array([0.,1.])

runit=0.
if runit > 0.:
    c=cluster()
    #c.readtable('/Users/rfinn/SDSS/baloghDR5/myclusters.cat')
    c.readtable('/Users/rfinn/SDSS/fieldDR5/myclusters.cat')
    
#(sigma,ngal05,ngal,ngal2,mstar05,mstar,mstar2,sfr05,sfr,sfr2)=(c.sigma,
#minimize by simplex method
    
#print "least-squares"
    p_init = scipy.array([0.,1.5])
#minimize via conjugate gradient
#p_opt = scipy.optimize.fmin_cg(objective_func, p_init, fprime)
#p_opt = scipy.optimize.leastsq(leastsq_func,p_init,args=(x,y),full_output=0,col_deriv=1)
#print p_opt
    
    
#list=[ngal,mstar,sfr]
#list=[c.ngal05,c.ngal,c.ngal2,c.mstar05,c.mstar,c.mstar2,c.sfr05,c.sfr,c.sfr2]
#name=['ngal05','ngal1','ngal2','mstar05','mstar1','mstar2','sfr05','sfr1','sfr2']
    list=[c.sffrac05*c.ngal05,c.sffrac*c.ngal,c.sffrac2*c.ngal2]
    name=['Nsf05','Nsf1','Nsf2']
    x=c.sigma
    for i in range(len(list)):
	print name[i],": fitting powerlaw to "
	y=list[i]
	p_opt = scipy.optimize.fmin(objective_func,p_init,(x,y),full_output=1,maxiter=1000)
	print p_opt
	t=p_opt[0]
	aopt=t[0]
	bopt=t[1]

	bootstrap(x,y)
	print "optimum a,b = %1.4e %1.4e"%(aopt,bopt)
