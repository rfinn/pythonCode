#!/usr/bin/env python

from pylab import *
import atpy
import numpy as np
import os

from lf_double_powerlaw import lf_double_powerlaw
from lf_double_powerlaw import intlf_double_powerlaw

def schechter(x,alpha,Lstar,phistar):
    y=phistar*(x/Lstar)**alpha*exp(-1.*x/Lstar)
    return y
def schechterlogL(x,alpha,Lstar,phistar):
    y=phistar*(x/Lstar)**(alpha+1)*exp(-1.*x/Lstar)
    return y



# set number of points
npoints=5000000

simschechter=0
simdoublepowerlaw=1
if simschechter:
    # LF parameters
    alpha=-1.5
    logLstar=0
    phistar=1
    eta=0.01 # min L/Lstar cut so that LF w/negative alpha does not diverge

    # draw from gamma function
    lum=np.random.gamma(alpha+2,10.**logLstar,npoints)
    
    # draw another set of random numbers
    test=np.random.uniform(size=npoints)
    keepflag=(lum/10.**logLstar > eta) & (test < eta*10.**logLstar/lum)
    ngood=sum(keepflag)
    logL=log10(lum[keepflag])
    #yfit=schechterlogL(10.**bincenters,alpha,10.**logLstar,phistar)
    outfile='testLF3.fits'
if simdoublepowerlaw:
    alpha=.36
    beta=2.17
    logLstar=9.6
    phistar=1.2e-1
    eta=0.01 # min L/Lstar cut so that LF w/negative alpha does not diverge
    logLmin=eta

    # draw from luminosity
    loglum=np.random.uniform(low=logLstar-2,high=logLstar+2,size=npoints)
    # calculate the prob corresponding to that luminosity
    prob=lf_double_powerlaw(loglum,phistar,logLstar,alpha,beta)/intlf_double_powerlaw(logLmin,phistar,logLstar,alpha,beta)
    # draw random number
    test=np.random.uniform(size=npoints)
    # if prob > than test, keep lum
    keepflag=prob>test
    logL=loglum[keepflag]
    ngood=sum(keepflag)

    outfile='test_double_powerlaw.fits'
figure()
dl=.2
lbins=arange(logLstar-2,logLstar+1+dl,dl)
t=hist(logL)
clf()
binedges=t[1]
bincenters=zeros(len(t[0]),'f')
for i in range(len(bincenters)):
    bincenters[i]=0.5*(binedges[i]+binedges[i+1])
    
plot(bincenters,log10(t[0]/(dl*ngood)),'ko')
errup=log10(t[0]+sqrt(t[0])) - log10(t[0])
errdown=log10(t[0]) - log10(t[0]-sqrt(t[0]))
yerror=array(zip(errdown,errup),'f').T
errorbar(bincenters,log10(t[0]/(dl*ngood)),yerr=yerror,fmt=None)
if simdoublepowerlaw:
    yfit=lf_double_powerlaw(bincenters,phistar,logLstar,alpha,beta)
if simschechter:
    yfit=schechterlogL(10.**bincenters,alpha,10.**logLstar,phistar)
plot(bincenters,log10(yfit),'r-')
show()
lf=atpy.Table()
lf.add_column('logL',logL,unit='log10(L/Lsun)')
lf.write(outfile,overwrite='yes')
