#!/usr/bin/env python

import glob
import sys, glob, os
import numarray as n
import Numeric as N
#import scipy
from math import *
from pyraf import iraf
import pyfits
import mystuff as my
#import ppgplot
#pisces noise stats
READNOISE=4.
GAIN=4.35


prefix=['D17s','EGSdJ','EGSdn','EGSsn','standA']

iraf.images()
iraf.images.imutil()
iraf.immatch()


for name in prefix:
    i=0
    s=name+"*.fits"
    files=glob.glob(s)
    for file in files:
	(pre,junk)=file.split('.')
	index=pre[(len(pre)-4):len(pre)]

	iraf.imgets(image=file,param='FILTER')#get RA of image
	filter=iraf.imgets.value
	if my.beginsWith('D17s',file) > .1:
	    outname='D17s'+str(file[7])+"-"+str(filter)+"-"+str(index)+".fits"
	    s="cp "+file+" "+outname+" \n"
	    os.system(s)
    

	if my.beginsWith('EGSd',file) > .1:
	    outname='EGS'+"-"+str(filter)+"-n1-"+str(index)+".fits"
	    s="cp "+file+" "+outname+" \n"
	    os.system(s)
    
	if my.beginsWith('EGSs',file) > .1:
	    outname='EGSs'+str(int(file[7])+1)+"-"+str(filter)+"-"+str(index)+".fits"
	    s="cp "+file+" "+outname+" \n"
	    os.system(s)
    
	if my.beginsWith('stand',file) > .1:
	    outname='s852c'+"-"+str(filter)+"-"+str(index)+".fits"
	    s="cp "+file+" "+outname+" \n"
	    os.system(s)

files=glob.glob("*-*.fits")
allprefix=[]
for file in files:
    t=file.split('-')
    s=str(t[0])+'-'+str(t[1])
    allprefix.append(s)
    
prefix=[allprefix[0]]#keep only unique prefixes
for i in range(1,len(allprefix)):
    t=allprefix[i-1]
    t1=allprefix[i]
    if t1.find(t) > -1:
	continue
    prefix.append(t1)
print "unique prefixes = ",prefix


for pre in prefix:
    s=pre+"*"
    files=glob.glob(s)
    revfiles=[]
    for i in range(len(files)):
	revfiles.append(files[len(files)-(i+1)])
    if pre.find('EGS-J') > -1:
	revfiles=files
    for i in range(len(files)):
	n=len(revfiles)-i
	if pre.find('EGS-J') > -1:
	    n=i+1
	if n < 10:
	    index='000'+str(n)
	if n > 9:
	    index='00'+str(n)
	if n > 99:
	    index='0'+str(n)

	file=revfiles[i]#start w/max index first so I don't write over original files when adding 1 to index
	
	(t1,t2)=file.split('.')
	firstpart=str(t1[0:(len(t1)-4)])
	s=firstpart+str(index)+'.fits'

	s1='mv '+str(file)+' '+s
	print s1,t1,t2,firstpart
	os.system(s1)

os.system('mv *-* renamed/.')
os.system('ls *.fits > inlist')
os.system('ls renamed/*.fits > outlist')
os.system('wc -l *list')
