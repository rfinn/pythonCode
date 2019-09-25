#!/usr/bin/env python
import sys, os
import Numeric as N
import scipy
from math import *
from mystuff import *
import ppgplot
import random
numcl=N.array([0,1,2],'i')
dA=N.array([5.02,5.14,5.25],'f')
dA=dA/1000.
nneigh=5
for ncl in numcl:
    if ncl == 0:
        infile=open('cl1040-1155.RF_cluster_Mr','r')
        outfile=open('cl1040-1155.sigma10','w')
    if ncl == 1:
        infile=open('cl1054-1245.RF_cluster_Mr','r')
        outfile=open('cl1054-1245.sigma10','w')
    if ncl == 2:
        infile=open('cl1216-1201.RF_cluster_Mr','r')
        outfile=open('cl1216-1201.sigma10','w')
    
    name=[]
    oldname=[]
    ra=[]
    dec=[]
    Mr=[]
    
    for line in infile:
        if line.find('#') > -1: #skip lines with '#' in them
            continue
        t=line.split()
        name.append(t[0])
        oldname.append(t[1])
        ra.append(float(t[2]))
        dec.append(float(t[3]))
        Mr.append(float(t[3]))
    ra=N.array(ra,'f')
    dec=N.array(dec,'f')
    Mr=N.array(Mr,'f')
    infile.close()
    dNN = N.zeros(len(ra),'f')
    sigma10 = N.zeros(len(ra),'f')
    for i in range(len(dNN)):
        temp=N.sqrt((ra[i]-ra)**2+(dec[i]-dec)**2)
        temp=N.sort(temp)
        dNN[i]=temp[nneigh]#distance to 5th NN, no edge effect accounted for
        sigma10[i]=(1.*nneigh)/(dNN[i]*3600.*dA[ncl])**2
    for i in range(len(sigma10)):
        outfile.write("%s %8.2f \n" % (name[i],sigma10[i]))
        
    outfile.close()
