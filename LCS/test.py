#!/usr/bin/env python

from pylab import *

def readmorphcat(prefix):
    infile='/home/rfinn/research/LocalClusters/Morphologies2010/Final.'+prefix+'.txt'
    infile1=open(infile,'r')
    name=[]
    morph=[]
    disturb=[]
    for line in infile1:
        t=line.split()
        n=t[0]
        myn=n.split('-')
        print len(myn),myn
        if len(myn) > 4:
            s='-'+myn[2]
            name.append(s)
        else:
            name.append(myn[1])
            morph.append(t[1])
            try:
                disturb.append(t[2])
            except:
                print 'trouble in paradise'
                disturb.append('-99')
    return name,morph,disturb


def localdensity(x,y,xref,yref,n1,n2):#find local density, using nearest neighbors n1 through n2
    DA=my.DA(self.cz,h)
    sigma=zeros(len(x),'d')
    for i in range(len(x)):
        d=sqrt((x[i]-xref)**2+(y[i]-yref)**2)*3600.*DA#d in kpc
        d=d/1000.#d  in Mpc
        d.sort()#sort in ascending order, zeroth element is distance from galaxy to itself
        sig=0
        for j in range(n1,n2+1):
            sig += (1.*j)/(d[j])**2
        sigma[i]=1./(4.*pi)*sig
    return sigma
