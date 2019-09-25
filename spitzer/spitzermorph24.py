#!/usr/bin/env python
import sys, os
import Numeric as N
import numarray as n
#import scipy
from math import *
import mystuff as my
#import ppgplot
#import random
import sets
#import steve
import pylab
from pyraf import iraf
#from matplotlib import rc
#rc('font',family='serif', style='normal', variant='normal',weight='bold', stretch='normal', size='large')

#rcfigure.figsize   : 8, 7    # figure size in inches
class cluster:
    def readdata(self,input1,input2):
	infile1=open(input1,'r')
	self.morphall=[]
	for line in infile1:
	    if line.find('#') > -1:
		continue
	    f=line.split()
	    self.morphall.append(float(f[1]))
	self.morphall=N.array(self.morphall,'f')
	infile1.close()

	infile2=open(input2,'r')
	self.morph24=[]
	self.f24=[]
	self.errf24=[]
	for line in infile2:
	    if line.find('#') > -1:
		continue
	    f=line.split()
	    self.morph24.append(float(f[1]))
	    self.f24.append(float(f[2]))
	    self.errf24.append(float(f[3]))
	self.morph24=N.array(self.morph24,'f')
	self.f24=N.array(self.f24,'f')
	self.errf24=N.array(self.errf24,'f')
	infile2.close()


def plotmorphall():

    pylab.subplots_adjust(left=0.1, right=.9,bottom=.1,top=0.9,wspace=.3,hspace=.2)
    p1=pylab.subplot(321)
    plotmorph(cl1232.morphall,cl1232.morph24,'CL1232, z=0.54')
    p1=pylab.subplot(322)
    plotmorph(cl1037.morphall,cl1037.morph24,'CL1037, z=0.58')
    p1=pylab.subplot(323)
    plotmorph(cl1227.morphall,cl1227.morph24,'CL1227, z=0.63')
    pylab.ylabel('Ngal (%)',fontsize=20)
    p1=pylab.subplot(324)
    plotmorph(cl1040.morphall,cl1040.morph24,'CL1040, z=0.70')
    p1=pylab.subplot(325)
    plotmorph(cl1354.morphall,cl1354.morph24,'CL1040, z=0.76')
    p1=pylab.subplot(326)
    plotmorph(cl1216.morphall,cl1216.morph24,'CL1216, z=0.79')

    pylab.text(-4.,-10.,r'Morphology',horizontalalignment='center',verticalalignment='center',fontsize=20)

    pylab.savefig('morph24all.eps')
    #pylab.show()


def plotmorph(mall,m24,prefix):
    xbins=N.array([-7,-6,-5,-2,1,2,3,4,5,6,7,8,9,10,11,66,111],'f')
    nmall=N.zeros(len(xbins),'f')
    nm24=N.zeros(len(xbins),'f')
    ntot=len(mall)
    ntot24=len(m24)
    i=0
    for x in xbins:
	for m in mall:
	    if abs(m-x) < .1:
		nmall[i]=nmall[i]+1.
	for m in m24:
	    if abs(m-x) < .1:
		nm24[i]=nm24[i]+1.
	i=i+1
    Nbin = len(xbins)
    nmallerr=N.sqrt(nmall)/N.sum(nmall)*100.
    nmall=nmall/N.sum(nmall)*100.
    nm24err=N.sqrt(nm24)/N.sum(nm24)*100.
    nm24=nm24/N.sum(nm24)*100.
    #ind = N.arange(0,2*Nbin,2)  # the x locations for the groups
    ind = N.arange(Nbin)  # the x locations for the groups
    width = 0.35       # the width of the bars
    p1 = pylab.bar(ind, nmall, width, color='r', yerr=nmallerr)

    p2 = pylab.bar(ind+width, nm24, width, color='b', yerr=nm24err)

    pylab.text(.5,34.,prefix,fontsize=12.)
    pylab.xlim(-width,len(ind))
    pylab.ylim(0.,40.)
    #pylab.xticks(ind+width, ('St', 'C', 'E', 'S0', 'Sa', 'Sab', 'Sb','Sbc','Sc','Scd','Sd','Sdm','Sm','Im','Irr','?','N/A') ,fontsize=10)
    ytick=N.arange(0.,40.,10.)
    #pylab.yticks(ytick)
    pylab.xticks(ind+width, ('St', 'C', 'E', 'S0', 'Sa', '', 'Sb','','Sc','','Sd','','Sm','Im','Irr','?','N/A') ,fontsize=10)

    pylab.legend( (p1[0], p2[0]), ('All ('+str(ntot)+' gal)', '24um ('+str(ntot24)+' gal)'), shadow=True)
    #pylab.show()
    #pylab.savefig(prefix+'morph.eps')



cl1040=cluster()
cl1040.readdata('cl1040AllMembers.dat','cl1040Members24.dat')
#plotmorph(cl1040.morphall,cl1040.morph24,'cl1040')

cl105412=cluster()
cl105412.readdata('cl105412AllMembers.dat','cl105412Members24.dat')
#plotmorph(cl1040.morphall,cl1040.morph24,'cl1040')

cl1216=cluster()
cl1216.readdata('cl1216AllMembers.dat','cl1216Members24.dat')
#plotmorph(cl1216.morphall,cl1216.morph24,'cl1216')

cl1354=cluster()
cl1354.readdata('cl1354AllMembers.dat','cl1354Members24.dat')
#plotmorph(cl1354.morphall,cl1354.morph24,'cl1354')

cl1227=cluster()
cl1227.readdata('cl1227AllMembers.dat','cl1227Members24.dat')
#plotmorph(cl1227.morphall,cl1227.morph24,'cl1227')

cl1232=cluster()
cl1232.readdata('cl1232AllMembers.dat','cl1232Members24.dat')
#plotmorph(cl1232.morphall,cl1232.morph24,'cl1232')

cl1037=cluster()
cl1037.readdata('cl1037AllMembers.dat','cl1037Members24.dat')
#plotmorph(cl1037.morphall,cl1037.morph24,'cl1037')

plotmorphall()
