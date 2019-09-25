#!/usr/bin/env python
import numarray as N
import mystuff as my

h=0.7#H100
id=[]
z=[]
fHa=[]

infile=open('targets','r')
for line in infile:
    if line.find('#') > -1:
	continue
    t=line.split()
    id.append(t[0])
    z.append(float(t[6]))
    fHa.append(float(t[5]))
z=N.array(z,'f')
fHa=N.array(fHa,'d')*1.e-16#f(HA) 1e-16 erg/s/cm2
dL=N.zeros(len(z),'d')
for i in range(len(dL)):
    dL[i]=my.dLcm(z[i],h)
LHa=fHa*4.*N.pi*dL**2
sfr=LHa*7.9e-42
print "id     z      SFR  SFR/dL^2 f/(1Msun/yr@z=0.35) expt   Ncycles (70um)"
for i in range(len(z)):
    r=1.e6*sfr[i]/(my.dL(z[i],h)**2)
    r2=r/(3.12/11.5*10.)
    exp=1/N.sqrt(r2)
    print "%s %5.4f %5.1f   %5.2f         %5.2f        %5.2f   %2i"%(id[i], z[i], sfr[i],r,r2,exp,round(10.*exp))
zmin=min(z)
zmax=max(z)
r=(my.dL(zmax,h)/my.dL(zmin,h))**2
print "(max dL/min dL)^2 = ",r
print "Average SFR = ",N.average(sfr)
