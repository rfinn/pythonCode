#!/usr/bin/env python

import numarray as N
#check translation of ediscs ra and dec 
era=[]
edec=[]
gra=[]
gdec=[]
infile1=open('cl1018radec.dat','r')
for line in infile1:
    t=line.split()
    era.append(float(t[0]))
    edec.append(float(t[1]))
infile1.close()
era=N.array(era,'d')
edec=N.array(edec,'d')

infile1=open('cl1018-GMOSradec.dat','r')
for line in infile1:
    t=line.split()
    gra.append(float(t[0]))
    gdec.append(float(t[1]))
infile1.close()
gra=N.array(gra,'d')
gdec=N.array(gdec,'d')

dr=N.sqrt((gra-era)**2 + (gdec-edec)**2)
dr=dr*3600.
for i in range(len(dr)):
    print dr[i]
print "Average difference = ",N.average(dr),max(dr),min(dr)
