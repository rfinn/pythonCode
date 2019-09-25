#!/usr/bin/env python
import sys, os
import Numeric as N
import scipy
from math import *
import mystuff as my
import ppgplot
import random

f0=1600.e-23#Jband ZP in ergs/s/cm^2/Hz
h=1.
pi=3.14159
kenn=7.9e-42
dust=2.5
NII=.77
ZP = N.array([24.94, 24.93, 24.93, 22.04],'f')

fzp=f0*10.**(-1*ZP/2.5)

z = N.array([.704,.748,.794,.832],'f')#cluster redshift

r = N.array([.0604,.06126,.07382,.065],'f')#NB/J throughput

l = .6563*(1.+z)*1.e-6#central wavelength in m

dl=0.02*l#bandwidth

dnu=3.e8*dl/l**2#bandwidth in Hz

dLcm=N.zeros(len(z),'f')
for i in range(len(z)):
    d=my.dL(z[i],h)
    dLcm=d*3.08e24
    
fzpNB=fzp/r*dnu

sfrconv=fzpNB*4*pi*dLcm**2*NII*dust*kenn

for i in range(len(z)):
    print i,z[i],sfrconv[i],fzpNB[i]
