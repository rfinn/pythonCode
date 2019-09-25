#!/usr/bin/env python

import pylab 
import Numeric as N

f24IRz=[]
f24IRconv=[]
errf24IRconv=[]
#infile2=open('/Users/rfinn/clusters/spitzer/Mastertables/fluxconv24toIR.v2.dat','r')
infile2=open('/Users/rfinn/clusters/spitzer/Mastertables/fluxconv24toIR.dat','r')
for line in infile2:
	if line.find('#') > -1:
		continue
	t=line.split()
	for j in range(len(t)):
		t[j]=float(t[j])
	f24IRz.append(t[0])
	f24IRconv.append(t[6])
	errf24IRconv.append(t[7])
	#f24IRconv.append(t[1])
	#errf24IRconv.append(t[2])
f24IRz=N.array(f24IRz,'f')
f24IRconv=N.array(f24IRconv,'f')
errf24IRconv=N.array(errf24IRconv,'f')
infile2.close()
pylab.cla()
pylab.clf()
pylab.plot(f24IRz,f24IRconv,'k-')
y=f24IRconv+errf24IRconv
pylab.plot(f24IRz,y,'k--')
y=f24IRconv-errf24IRconv
pylab.plot(f24IRz,y,'k--')
pylab.axis([0.2,1.,3.,17.])
pylab.xlabel(r'$\rm z$',fontsize=32)
pylab.ylabel(r'$\rm F(3-1100 \mu m)/\nu F_\nu (24\mu m)$',fontsize=32)
pylab.savefig('Ave24toIRconvvsz.eps')
