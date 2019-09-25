#!/usr/bin/env python
from pylab import *


lookbackt=array([6.404,5.047,0.779,.312],'f')
flirg=array([19.,8.,.00347,0.],'f')
erryup=array([4.,2.,.008,.0067],'f')
errydown=array([4.,2.,.003,.00],'f')
plot(lookbackt,flirg,'ko',markersize=12,label=r'$\rm Cluster$')
yerrall=zeros([2,len(erryup)],'f')
yerrall[0]=erryup#+pylab.log10(scale)
yerrall[1]=errydown#+pylab.log10(scale)
errorbar(lookbackt,flirg,yerr=yerrall,fmt=None,ecolor='k',label='_nolegend_')

lookbackt=array([6.404,5.047],'f')
flirg=array([31.,14.],'f')
erryup=array([1.,1.],'f')
errydown=array([1.,1.],'f')

plot(lookbackt,flirg,'k^',markersize=12,label=r'$\rm Field$')
yerrall=zeros([2,len(erryup)],'f')
yerrall[0]=erryup#+pylab.log10(scale)
yerrall[1]=errydown#+pylab.log10(scale)
errorbar(lookbackt,flirg,yerr=yerrall,fmt=None,ecolor='k',label='_nolegend_')


xlabel(r'$\rm Lookback \ Time \ (Gyr)$')
ylabel(r'$\rm LIRG \ Fraction (\%) $')
x=arange(0.,8,.1)
y=exp(x/2.2)-1.
plot(x,y,'k--',color='0.5',label=r'$\rm y = e^{t/2.2} -1$')

plot(x,1.63*y,'k:',color='0.5',label=r'$\rm y = 1.6e^{t/2.2} -1$')

y=exp((x+1)/2.2)-1.
plot(x,y,'k-.',color='0.5',label=r'$\rm y = e^{(t+1)/2.2} -1$')


legend(loc='upper left',numpoints=2)

#firfield=array([
axis([0.,7.5,-1.,35])
savefig('sfrdecline.eps')
