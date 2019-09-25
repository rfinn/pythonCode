#!/usr/bin/env python
"""
useage 
spitzer1.py mode

mode = 0 : calculate quantities, match ediscs, write mastertable
mode = 1 : read master table

"""
import sys, os
import Numeric as N
import numarray as n
import glob
#import scipy
from math import *
import mystuff as my
#import ppgplot
#import random
import sets
#import steve
import matplotlib
import pylab
from pyraf import iraf
iraf.imcoords()
#frame=int(sys.argv[1])
frame=1
catalogpath='/Users/rfinn/clusters/ediscs/catalogs/'
delta=2. #max allowed offset (in arcseconds) between matched sources
delta=delta/3600.
h=.7#H100
Lsol=3.9e33

#mode=float(sys.argv[1])

mode=1

matplotlib.rc('legend',numpoints=1,fontsize=12,markerscale=1)
def findnearest(x1,y1,x2,y2,delta):
	dmin = 100000000000000000000000000000000000.
	matchflag=1
	nmatch=0
	for i in range(len(x2)):
		d = N.sqrt((x1-x2[i])**2 + (y1-y2[i])**2)
		if d < delta:
			nmatch=nmatch+1
		if d < dmin:
			dmin = d
			imatch = i

	
	if dmin > delta:
		imatch = 0
		matchflag = 0
	return imatch, matchflag,nmatch
def ratioerror(a,b):
    ratio=a/b
    try:
        err = a/b*N.sqrt(1./a + 1./b)
    except:
        err = 0.
    return ratio,err


class ediscs24:
    def __init__(self):#individual galaxy properties
        print "dude - 24 micron data!"
	self.prefix=prefix
	self.idra=idra
	self.iddec=iddec
	self.mipsimage='/Users/rfinn/clusters/spitzer/'+prefix+'/mips/24/'+rvalue+'/ch1/bcd/pbcd/Combine/mosaic.fits'
	self.z=z
	self.rac=float(rac)
	self.decc=float(decc)
	self.DL = my.dL(z,h)


    def readmaster(self):	
	mastertable=self.prefix+'mastertable24.dat'
	self.readmastertable(mastertable)#comment this line and uncomment line above

    def readmastertable(self,file):
	
	nmemb=0.
	nspecmemb=0.
        input=open(file,'r')
        #get number of galaxies
        ngal=0
	for line in input:
		if line.find('#') > -1: #skip lines with '#' in them
			continue
		if line.find('\\') > -1: #skip lines with '#' in them
			continue
		if line.find('|') > -1: #skip lines with '#' in them
			continue
		t=line.split()
		if (float(t[32]) > 23.): #skip of R > 23
			continue

		if (float(t[32]) < 5.): #estimate R from 0.5(V + I) for galaxies w/no R data
			#print "Rmag < 5: ",t[0]
			testR = 0.5*(float(t[30])+float(t[34]))-0.2#0.2 is fudge factor to make sure all galaxies with R > 23 are included
			if testR > 23.:
				continue
		s= self.prefix#do custom photoz cut for fields w/multiple structures
		if (s.find('1103') > -1):
			if (float(t[48]) < .32) or (float(t[48]) > 1.25):
				continue
		elif (s.find('1301') > -1):
			if (float(t[48]) < .1) or (float(t[48]) > .78):
				continue
		elif (s.find('1037') > -1):
			if (float(t[48]) < .12) or (float(t[48]) > .88):
				continue
		elif (s.find('1232') > -1):
			if (float(t[48]) < .28) or (float(t[48]) > .93):
				continue	
		elif abs(float(t[48])-self.z) > 0.3:#48 - self.bestz
			#print float(t[48]),self.z
			continue
		if float(t[41]) > -.9:#if newmatchflag > -1, then it already has spectroscopy
			continue
		ngal=ngal+1
        input.close()


	self.specmembflag = N.zeros(ngal,'f')
	self.matchflagspec = N.zeros(ngal,'f')
	self.fflat24 = N.zeros(ngal,'f')
	self.errfflat24 = N.zeros(ngal,'f')
	self.matchflagediscsflat24 = N.zeros(ngal,'f')
	

        self.ediscsID = []
	self.ra = N.zeros(ngal,'f')
	self.dec = N.zeros(ngal,'f')
	self.xcorr = N.zeros(ngal,'f')
	self.ycorr = N.zeros(ngal,'f')
	self.starflag = N.zeros(ngal,'f')
	self.EW = N.zeros(ngal,'f')
	self.EWerr = N.zeros(ngal,'f')
	self.SFR = N.zeros(ngal,'f')
	self.SFRerr = N.zeros(ngal,'f')
	self.matchflagha = N.zeros(ngal,'f')
	self.SFflag = N.zeros(ngal,'f')
	self.gimtype = N.zeros(ngal,'f')
	self.matchflagmorphgimtype = N.zeros(ngal,'f')	
	self.vistype = N.zeros(ngal,'f')
	self.matchflagvistype = N.zeros(ngal,'f')
	self.matchflag24 = N.zeros(ngal,'f')
	self.flux24 = N.zeros(ngal,'f')
	self.flux24err = N.zeros(ngal,'f')
	self.nmatchediscs24 = N.zeros(ngal,'f')
	self.misoV = N.zeros(ngal,'f')
	self.misoeVapsim = N.zeros(ngal,'f')
	self.misoR = N.zeros(ngal,'f')
	self.misoeRapsim = N.zeros(ngal,'f')
	self.misoI = N.zeros(ngal,'f')
	self.misoeIapsim = N.zeros(ngal,'f')
	self.misoJ = N.zeros(ngal,'f')
	self.misoeJapsim = N.zeros(ngal,'f')	
	self.misoK = N.zeros(ngal,'f')
	self.misoeKapsim = N.zeros(ngal,'f')
	self.magV = N.zeros(ngal,'f')
	self.mageVapsim = N.zeros(ngal,'f')
	self.magR = N.zeros(ngal,'f')
	self.mageRapsim = N.zeros(ngal,'f')
	self.magI = N.zeros(ngal,'f')
	self.mageIapsim = N.zeros(ngal,'f')
	self.magJ = N.zeros(ngal,'f')
	self.mageJapsim = N.zeros(ngal,'f')
	self.magK = N.zeros(ngal,'f')
	self.mageKapsim = N.zeros(ngal,'f')
	self.membflag = N.zeros(ngal,'f')
	self.newmatchflag = N.zeros(ngal,'f')	
	self.defmembflag = N.zeros(ngal,'f')
	self.specz = N.zeros(ngal,'f')
	self.spectype = N.zeros(ngal,'f')
	self.specEWOII = N.zeros(ngal,'f')
	self.matchflagspecediscs = N.zeros(ngal,'f')
	self.specEWOIIflag = N.zeros(ngal,'f')
	self.bestz = N.zeros(ngal,'f')
	self.lowz = N.zeros(ngal,'f')
	self.highz = N.zeros(ngal,'f')
	self.wmin = N.zeros(ngal,'f')
	self.Pclust = N.zeros(ngal,'f')
	self.LUlowzclust = N.zeros(ngal,'f')
	self.LUbestzclust = N.zeros(ngal,'f')
	self.LUhighzclust = N.zeros(ngal,'f')	
	self.LBlowzclust = N.zeros(ngal,'f')
	self.LBbestzclust = N.zeros(ngal,'f')
	self.LBhighzclust = N.zeros(ngal,'f')	
	self.LVlowzclust = N.zeros(ngal,'f')
	self.LVbestzclust = N.zeros(ngal,'f')
	self.LVhighzclust = N.zeros(ngal,'f')	
	self.LRlowzclust = N.zeros(ngal,'f')
	self.LRbestzclust = N.zeros(ngal,'f')
	self.LRhighzclust = N.zeros(ngal,'f')	
	self.LIlowzclust = N.zeros(ngal,'f')
	self.LIbestzclust = N.zeros(ngal,'f')
	self.LIhighzclust = N.zeros(ngal,'f')	
	self.LJlowzclust = N.zeros(ngal,'f')
	self.LJbestzclust = N.zeros(ngal,'f')
	self.LJhighzclust = N.zeros(ngal,'f')	
	self.LKlowzclust = N.zeros(ngal,'f')
	self.LKbestzclust = N.zeros(ngal,'f')
	self.LKhighzclust = N.zeros(ngal,'f')
	self.fluxK = N.zeros(ngal,'f')
	self.UBlowzclust = N.zeros(ngal,'f')
	self.UBbestzclust = N.zeros(ngal,'f')
	self.UBhighzclust = N.zeros(ngal,'f')
	self.BVlowzclust = N.zeros(ngal,'f')
	self.BVbestzclust = N.zeros(ngal,'f')
	self.BVhighzclust = N.zeros(ngal,'f')
	self.UVlowzclust = N.zeros(ngal,'f')
	self.UVbestzclust = N.zeros(ngal,'f')
	self.UVhighzclust = N.zeros(ngal,'f')
	self.matchflagediscsirac = N.zeros(ngal,'f')
	self.iracf1 = N.zeros(ngal,'f')
	self.iracf2 = N.zeros(ngal,'f')
	self.iracf3 = N.zeros(ngal,'f')
	self.iracf4 = N.zeros(ngal,'f')
	self.erriracf1 = N.zeros(ngal,'f')
	self.erriracf2 = N.zeros(ngal,'f')
	self.erriracf3 = N.zeros(ngal,'f')
	self.erriracf4 = N.zeros(ngal,'f')
	self.iracsexflag0 = N.zeros(ngal,'f')
	self.iracsexflag1 = N.zeros(ngal,'f')
	self.iracwch1 = N.zeros(ngal,'f')
	self.iracwch2 = N.zeros(ngal,'f')
	self.iracwch3 = N.zeros(ngal,'f')
	self.iracwch4 = N.zeros(ngal,'f')
	self.iracwmin = N.zeros(ngal,'f')
	self.nmatchediscsirac = N.zeros(ngal,'f')
	self.L24 = N.zeros(ngal,'d')
	self.errL24 = N.zeros(ngal,'d')
	self.LHa = N.zeros(ngal,'d')
	self.errLHa = N.zeros(ngal,'d')
	self.snr24 = N.zeros(ngal,'d')
	self.imagex24 = N.zeros(ngal,'d')
	self.imagey24 = N.zeros(ngal,'d')
	

        input=open(file,'r')
        i=0
	k=0
        for line in input:
		if line.find('#') > -1: #skip lines with '#' in them
			continue
		if line.find('\\') > -1: #skip lines with '#' in them
			continue
		if line.find('|') > -1: #skip lines with '#' in them
			continue


		t=line.split()
		if (float(t[32]) > 23.): #skip of R > 23
			continue

		if (float(t[32]) < 5.): #estimate R from 0.5(V + I) for galaxies w/no R data
			#print "Rmag < 5: ",t[0]
			testR = 0.5*(float(t[30])+float(t[34]))-0.2#0.2 is fudge factor to make sure all galaxies with R > 23 are included
			if testR > 23.:
				continue
		s= self.prefix#do custom photoz cut for fields w/multiple structures
		if (s.find('1103') > -1):
			if (float(t[48]) < .32) or (float(t[48]) > 1.25):
				continue
		elif (s.find('1301') > -1):
			if (float(t[48]) < .1) or (float(t[48]) > .78):
				continue
		elif (s.find('1037') > -1):
			if (float(t[48]) < .12) or (float(t[48]) > .88):
				continue
		elif (s.find('1232') > -1):
			if (float(t[48]) < .28) or (float(t[48]) > .93):
				continue	
		elif abs(float(t[48])-self.z) > 0.3:#48 - self.bestz
			#print float(t[48]),self.z
			continue
		#if abs(float(t[48])-self.z) > 0.3:#48 - self.bestz
		#	print "abs(photoz - zcl) > 0.3: ",float(t[48]),self.z
		#	continue
		if float(t[41]) > -.9:#if newmatchflag > -1, then it already has spectroscopy
			print self.prefix,": gal already has spec z ",self.z,float(t[43]),float(t[41])
			continue

		self.ediscsID.append(t[0])

		for j in range(1,len(t)):
			t[j]=float(t[j])

			

		(self.ra[i],self.dec[i],self.xcorr[i],self.ycorr[i],self.starflag[i],self.EW[i],self.EWerr[i],self.SFR[i],self.SFRerr[i],self.matchflagha[i],self.SFflag[i],self.gimtype[i],self.matchflagmorphgimtype[i],self.vistype[i],self.matchflagvistype[i],self.matchflag24[i],self.flux24[i],self.flux24err[i],self.nmatchediscs24[i],self.misoV[i],self.misoeVapsim[i],self.misoR[i],self.misoeRapsim[i],self.misoI[i],self.misoeIapsim[i],self.misoJ[i],self.misoeJapsim[i],self.misoK[i],self.misoeKapsim[i],self.magV[i],self.mageVapsim[i],self.magR[i],self.mageRapsim[i],self.magI[i],self.mageIapsim[i],self.magJ[i],self.mageJapsim[i],self.magK[i],self.mageKapsim[i],self.membflag[i],self.newmatchflag[i],self.defmembflag[i],self.specz[i],self.spectype[i],self.specEWOII[i],self.matchflagspecediscs[i],self.specEWOIIflag[i],self.bestz[i],self.lowz[i],self.highz[i],self.wmin[i],self.Pclust[i],self.LUlowzclust[i],self.LUbestzclust[i],self.LUhighzclust[i],self.LBlowzclust[i],self.LBbestzclust[i],self.LBhighzclust[i],self.LVlowzclust[i],self.LVbestzclust[i],self.LVhighzclust[i],self.LRlowzclust[i],self.LRbestzclust[i],self.LRhighzclust[i],self.LIlowzclust[i],self.LIbestzclust[i],self.LIhighzclust[i],self.LJlowzclust[i],self.LJbestzclust[i],self.LJhighzclust[i],self.LKlowzclust[i],self.LKbestzclust[i],self.LKhighzclust[i],self.fluxK[i],self.UBlowzclust[i],self.UBbestzclust[i],self.UBhighzclust[i],self.BVlowzclust[i],self.BVbestzclust[i],self.BVhighzclust[i],self.UVlowzclust[i],self.UVbestzclust[i],self.UVhighzclust[i],self.matchflagediscsirac[i],self.iracf1[i],self.iracf2[i],self.iracf3[i],self.iracf4[i],self.erriracf1[i],self.erriracf2[i],self.erriracf3[i],self.erriracf4[i],self.iracsexflag0[i],self.iracsexflag1[i],self.iracwch1[i],self.iracwch2[i],self.iracwch3[i],self.iracwch4[i],self.iracwmin[i],self.nmatchediscsirac[i],self.L24[i],self.errL24[i],self.LHa[i],self.errLHa[i],self.snr24[i],self.imagex24[i],self.imagey24[i])=t[1:(len(t)-20)]
		If self.magR[i] < 0.5:
			self.magR[i]= 0.5*(self.magI[i] + self.magV[i])-0.2
	    #print i,self.defmembflag[i],self.matchflagvistype[i],self.vistype[i],self.matchflag24[i]
		if (self.defmembflag[i] > 0.) & (self.matchflagvistype[i] > 0.):
			k=k+1
		#print k,self.vistype[i],self.matchflag24[i]
		i=i+1


	input.close()

	self.photmembflag=N.zeros(len(self.membflag),'f')
	for i in range(len(self.membflag)):
	    if (self.membflag[i] > 0.) & (self.wmin[i] > 0.4):
		self.photmembflag[i]=1.

	pre=self.prefix
	#if pre.find('1216') > -1:#scale GTO data
	#    self.flux24=self.flux24*146.
	#    self.flux24err=self.flux24err*146.
	#    self.L24=self.L24*146.
	#    self.errL24=self.errL24*146.
	
	self.L24=self.L24/Lsol
	self.errL24=self.errL24/Lsol





    def readoptIRphot(self,file):#reading in ediscs photometry file
	    input=open(file,'r')
            #get number of galaxies
	    ngal=0
	    for line in input:
		    if line.find('#') > -1:
			    continue
		    ngal=ngal+1
	    input.close()
	    print "Number of galaxies in optIR file = ",ngal
		    
	    self.ediscsID = []
	    self.ra = N.zeros(ngal,'f')
	    self.dec = N.zeros(ngal,'f')
	    self.misoB = N.zeros(ngal,'f')
	    self.misoeBapsim = N.zeros(ngal,'f')
	    self.misoV = N.zeros(ngal,'f')
	    self.misoeVapsim = N.zeros(ngal,'f')
	    self.misoR = N.zeros(ngal,'f')
	    self.misoeRapsim = N.zeros(ngal,'f')
	    self.misoI = N.zeros(ngal,'f')
	    self.misoeIapsim = N.zeros(ngal,'f')
	    self.misoJ = N.zeros(ngal,'f')
	    self.misoeJapsim = N.zeros(ngal,'f')
	    self.misoK = N.zeros(ngal,'f')
	    self.misoeKapsim = N.zeros(ngal,'f')
	    self.EWha =-99.*N.ones(ngal,'f')
	    self.EWhaerr =-99.*N.ones(ngal,'f')
	    self.SFR =-99.*N.ones(ngal,'f')
	    self.SFRerr =-99.*N.ones(ngal,'f')
	    self.lha =-99.*N.ones(ngal,'f')
	    self.lhaerr =-99.*N.ones(ngal,'f')
	    self.matchhaediscs = N.zeros(ngal,'i')
	    self.matchflaghaediscs = N.zeros(ngal,'f')
	    self.mautoB= N.zeros(ngal,'f')
	    self.mautoeBapsim= N.zeros(ngal,'f')
	    self.mautoV= N.zeros(ngal,'f')
	    self.mautoeVapsim= N.zeros(ngal,'f')
	    self.mautoR= N.zeros(ngal,'f')
	    self.mautoeRapsim= N.zeros(ngal,'f')
	    self.mautooI= N.zeros(ngal,'f')
	    self.mautoeIapsim= N.zeros(ngal,'f')
	    self.mautoJ= N.zeros(ngal,'f')
	    self.mautoeJapsim= N.zeros(ngal,'f')
	    self.mautoK= N.zeros(ngal,'f')
	    self.mautoeKapsim= N.zeros(ngal,'f')

	    self.magB = N.zeros(ngal,'f')
            self.mageBapsim = N.zeros(ngal,'f')

	    self.magV = N.zeros(ngal,'f')
            self.mageVapsim = N.zeros(ngal,'f')
            self.magR = N.zeros(ngal,'f')
            self.mageRapsim = N.zeros(ngal,'f')
            self.magI = N.zeros(ngal,'f')
            self.mageIapsim = N.zeros(ngal,'f')
            self.magJ = N.zeros(ngal,'f')
            self.mageJapsim = N.zeros(ngal,'f')
            self.magK = N.zeros(ngal,'f')
            self.mageKapsim = N.zeros(ngal,'f')
	    self.SF_flag = N.zeros(ngal,'f')
	    #initialize other flags
	    self.hstvisnumtype = N.zeros(ngal,'f')
	    self.matchmorphvisediscs = N.zeros(ngal, 'i') 
	    self.matchflagmorphvisediscs = N.zeros(ngal, 'i')
	    self.hstgimtype = N.zeros(ngal,'f')
	    self.matchmorphgimediscs = N.zeros(ngal,'i')
	    self.matchflagmorphgimediscs = N.zeros(ngal,'i')
	    self.fflat24 = N.zeros(ngal,'f')
            self.errfflat24 = N.zeros(ngal,'f')
            self.matchflagediscsflat24 = N.zeros(ngal, 'i')
	    self.matchediscsflat24 = N.zeros(ngal,'i')
	    self.fluxK = N.zeros(ngal,'d')

	    if file.find('vrijk') > -1:

		    input=open(file,'r')
		    i=0
		    for line in input:
			    if line.find('#') > -1: #skip lines with '#' in them
				    continue
			    if line.find('\\') > -1: #skip lines with '#' in them
				    continue
			    if line.find('|') > -1: #skip lines with '#' in them
				    continue
			    t=line.split()
			    self.ediscsID.append(t[0])
			    self.ra[i]=float(t[6])*15.
			    self.dec[i]=float(t[7])

			    self.misoV[i]=float(t[128])
			    self.misoeVapsim[i]=float(t[134])
			    self.misoR[i]=float(t[136])
			    self.misoeRapsim[i]=float(t[142])
			    self.misoI[i]=float(t[144])
			    self.misoeIapsim[i]=float(t[150])
			    self.misoJ[i]=float(t[152])
			    self.misoeJapsim[i]=float(t[158])
			    self.misoK[i]=float(t[160])
			    self.misoeKapsim[i]=float(t[166])
			    self.mautoV[i]=float(t[168])
			    self.mautoeVapsim[i]=float(t[174])
			    self.mautoR[i]=float(t[176])
			    self.mautoeRapsim[i]=float(t[182])
			    self.mautooI[i]=float(t[184])
			    self.mautoeIapsim[i]=float(t[190])
			    self.mautoJ[i]=float(t[192])
			    self.mautoeJapsim[i]=float(t[198])
			    self.mautoK[i]=float(t[200])
			    self.mautoeKapsim[i]=float(t[206])
		    
			    self.magV[i] = float(t[8])
			    self.mageVapsim[i] = float(t[14])
			    self.magR[i] = float(t[16])
			    self.mageRapsim[i] = float(t[22])
			    self.magI[i] = float(t[24])
			    self.mageIapsim[i] = float(t[30])
			    self.magJ[i] = float(t[32])
			    self.mageJapsim[i] = float(t[38])
			    self.magK[i] = float(t[40])
			    self.mageKapsim[i] = float(t[46])
		    #print self.magK[i]
			    i = i+1
		    input.close()
		    self.nediscs=ngal

	    if file.find('bvik') > -1:

		    input=open(file,'r')
		    i=0
		    for line in input:
			    if line.find('#') > -1: #skip lines with '#' in them
				    continue
			    if line.find('\\') > -1: #skip lines with '#' in them
				    continue
			    if line.find('|') > -1: #skip lines with '#' in them
				    continue
			    t=line.split()
			    self.ediscsID.append(t[0])
			    self.ra[i]=float(t[6])*15.
			    self.dec[i]=float(t[7])

			    self.misoB[i]=float(t[104])
			    self.misoeBapsim[i]=float(t[110])

			    self.misoV[i]=float(t[112])
			    self.misoeVapsim[i]=float(t[118])
			    self.misoI[i]=float(t[120])
			    self.misoeIapsim[i]=float(t[126])
			    self.misoK[i]=float(t[128])
			    self.misoeKapsim[i]=float(t[134])

			    self.mautoB[i]=float(t[136])
			    self.mautoeBapsim[i]=float(t[142])
			    self.mautoV[i]=float(t[144])
			    self.mautoeVapsim[i]=float(t[150])
			    self.mautooI[i]=float(t[152])
			    self.mautoeIapsim[i]=float(t[158])
			    self.mautoK[i]=float(t[160])
			    self.mautoeKapsim[i]=float(t[166])

			    self.magB[i] = float(t[8])
			    self.mageBapsim[i] = float(t[14])
			    self.magV[i] = float(t[16])
			    self.mageVapsim[i] = float(t[22])
			    self.magI[i] = float(t[24])
			    self.mageIapsim[i] = float(t[30])
			    self.magK[i] = float(t[32])
			    self.mageKapsim[i] = float(t[38])
		    #print self.magK[i]
			    i = i+1
		    input.close()
		    self.nediscs=ngal

	    if file.find('bvijk') > -1:

		    input=open(file,'r')
		    i=0
		    for line in input:
			    if line.find('#') > -1: #skip lines with '#' in them
				    continue
			    if line.find('\\') > -1: #skip lines with '#' in them
				    continue
			    if line.find('|') > -1: #skip lines with '#' in them
				    continue
			    t=line.split()
			    self.ediscsID.append(t[0])
			    self.ra[i]=float(t[6])*15.
			    self.dec[i]=float(t[7])

			    self.misoB[i]=float(t[128])
			    self.misoeBapsim[i]=float(t[134])
			    self.misoV[i]=float(t[136])
			    self.misoeVapsim[i]=float(t[142])
			    self.misoI[i]=float(t[144])
			    self.misoeIapsim[i]=float(t[150])
			    self.misoJ[i]=float(t[152])
			    self.misoeJapsim[i]=float(t[158])
			    self.misoK[i]=float(t[160])
			    self.misoeKapsim[i]=float(t[166])

			    self.mautoB[i]=float(t[168])
			    self.mautoeBapsim[i]=float(t[174])
			    self.mautoV[i]=float(t[176])
			    self.mautoeVapsim[i]=float(t[182])
			    self.mautooI[i]=float(t[184])
			    self.mautoeIapsim[i]=float(t[190])
			    self.mautoJ[i]=float(t[192])
			    self.mautoeJapsim[i]=float(t[198])
			    self.mautoK[i]=float(t[200])
			    self.mautoeKapsim[i]=float(t[206])
		    
			    self.magB[i] = float(t[8])
			    self.mageBapsim[i] = float(t[14])
			    self.magV[i] = float(t[16])
			    self.mageVapsim[i] = float(t[22])
			    self.magI[i] = float(t[24])
			    self.mageIapsim[i] = float(t[30])
			    self.magJ[i] = float(t[32])
			    self.mageJapsim[i] = float(t[38])
			    self.magK[i] = float(t[40])
			    self.mageKapsim[i] = float(t[46])
		    #print self.magK[i]
			    i = i+1
		    input.close()
		    self.nediscs=ngal
 
	    	
    def read24(self,file,mipsimage):
        input=open(file,'r')
        #get number of galaxies
        ngal=0
	for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue
            ngal=ngal+1
        input.close()


        self.id24 = N.zeros(ngal,'f')
        self.imagex24 = N.zeros(ngal,'f')
        self.imagey24  = N.zeros(ngal,'f')
        self.ra24 = N.zeros(ngal,'f')
        self.dec24 = N.zeros(ngal,'f')
        self.f24 = N.zeros(ngal,'d')#flux
        self.errf24 = N.zeros(ngal,'d')
        self.fap1 = N.zeros(ngal,'d')#flux in aperture 1 (1,1.5,2,2.6,3,3.5,4,4.5,5.,5.5) pixels
        self.fap2 = N.zeros(ngal,'d')#flux
        self.fap3 = N.zeros(ngal,'d')#flux
        self.fap4 = N.zeros(ngal,'d')#flux in ap 4 - this is one w/ap cor of 1.67 (Calzetti et al 2007)
        self.fap5 = N.zeros(ngal,'d')#flux
        self.fap6 = N.zeros(ngal,'d')#flux
        self.fap7 = N.zeros(ngal,'d')#flux
        self.fap8 = N.zeros(ngal,'d')#flux
        self.fap9 = N.zeros(ngal,'d')#flux
        self.fap10 = N.zeros(ngal,'d')#flux
        self.errfap1 = N.zeros(ngal,'d')#flux in aperture 1 (1,1.5,2,2.6,3,3.5,4,4.5,5.,5.5) pixels
        self.errfap2 = N.zeros(ngal,'d')#flux
        self.errfap3 = N.zeros(ngal,'d')#flux
        self.errfap4 = N.zeros(ngal,'d')#flux in ap 4 - this is one w/ap cor of 1.67 (Calzetti et al 2007)
        self.errfap5 = N.zeros(ngal,'d')#flux
        self.errfap6 = N.zeros(ngal,'d')#flux
        self.errfap7 = N.zeros(ngal,'d')#flux
        self.errfap8 = N.zeros(ngal,'d')#flux
        self.errfap9 = N.zeros(ngal,'d')#flux
        self.errfap10 = N.zeros(ngal,'d')#flux
        self.snr24 = N.zeros(ngal,'d')#SNR calculated by mopex
        self.deblend = N.zeros(ngal,'f')#SNR calculated by mopex


        input=open(file,'r')
        i=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue

            t=line.split()
            #(self.id24[i],self.imagex24[i],self.imagey24[i],self.ra24[i],self.dec24[i],self.f24[i],self.errf24[i])=(float(t[0]),float(t[1]),float(t[2]),float(t[3]),float(t[4]),float(t[5]),float(t[6]))
	    #print 'length of t = ',len(t)
	    #print t
	    (self.id24[i],self.ra24[i],self.dec24[i],self.imagex24[i],self.imagey24[i],self.f24[i],self.errf24[i],self.snr24[i],self.deblend[i],self.fap1[i],self.fap2[i],self.fap3[i],self.fap4[i],self.fap5[i],self.fap6[i],self.fap7[i],self.fap8[i],self.fap9[i],self.fap10[i],self.errfap1[i],self.errfap2[i],self.errfap3[i],self.errfap4[i],self.errfap5[i],self.errfap6[i],self.errfap7[i],self.errfap8[i],self.errfap9[i],self.errfap10[i])=(float(t[0]),float(t[3]),float(t[5]),float(t[8]),float(t[10]),float(t[13]),float(t[14]),float(t[18]),float(t[2]),float(t[23]),float(t[24]),float(t[25]),float(t[26]),float(t[27]),float(t[28]),float(t[29]),float(t[30]),float(t[31]),float(t[32]),float(t[33]),float(t[34]),float(t[35]),float(t[36]),float(t[37]),float(t[38]),float(t[39]),float(t[40]),float(t[41]),float(t[42]))
            i=i+1
        input.close()
	print 'number of galaxies in 24 micron catalog = ',i,len(self.f24),len(self.fap4)
	#pylab.errorbar(self.f24,self.fap4*1.67,yerr=self.errf24,xerr=self.errfap4*1.67)
	pylab.plot(self.f24,self.fap4*1.67,'bo')
	pylab.plot(self.f24,self.fap4,'ro')
	pylab.plot(self.f24,self.fap3,'go')
	pylab.plot(self.f24,self.fap2,'yo')
	pylab.plot(self.f24,self.fap1,'co')
	x=N.arange(0.,300.,10.)
	y=x
	pylab.plot(x,y,'k-')
	y=2.*x
	pylab.plot(x,y,'k-')
	y=3.*x
	pylab.plot(x,y,'k-')
	y=4.*x
	pylab.plot(x,y,'k-')
	y=5.*x
	pylab.plot(x,y,'k-')
	pylab.xlabel('Mopex F(24)')
	pylab.ylabel('Ap Flux (r=2.6 pix) w/apcor=1.67')
	pylab.axis([0.,50.,0.,50.])
	s=str(self.prefix)+'fluxcomp.eps'
	pylab.savefig(s)
        outfile=open('xy24.dat','w')
        for i in range(len(self.imagex24)):
            #print self.imagex24[i],self.imagey24[i],self.f24[i],self.errf24[i]
            string="%6.2f %6.2f %8.1f %6.2f \n" % (self.imagex24[i],self.imagey24[i],self.f24[i],self.errf24[i])
            outfile.write(string)
        outfile.close()

	dL=my.dLcm(z,.7)
	self.l24=self.f24*1.e-6#convert to Jy
	self.l24=self.l24*1.e-23#convert to erg/s/cm^2/Hz
	self.l24=self.l24*4*N.pi*dL**2#convert to erg/s/Hz
	self.l24=self.l24*2.448e12#multiply by bandwidth in Hz to get erg/s
	self.errl24=self.l24*(self.errf24/self.f24)
	#delta=2.5 #max allowed offset (in arcseconds) between matched sources
	#delta=delta/3600.
	x1=self.ra
	y1=self.dec
	x2=self.ra24
	y2=self.dec24


	
	
	self.matchediscs24 = N.zeros(len(x1), 'i')
	self.matchflagediscs24 = N.zeros(len(x1), 'i')
	self.nmatchediscs24 = N.zeros(len(x1), 'i')
	outfile=open('wcsin','w')
	for i in range(len(x1)):
		s="%8.2f %8.2f \n"%(x1[i],y1[i])
		outfile.write(s)
	outfile.close()
	try:
		os.system('rm wcsout')
	except:
		print "running wcsctran for first time"
	iraf.wcsctran('wcsin','wcsout',image=mipsimage,inwcs='world',outwcs='physical',verbose='no')
	xpix=N.zeros(len(x1),'f')
	ypix=N.zeros(len(x1),'f')
	outfile=open('wcsout','r')
	i=0
	for line in outfile:
		f=line.split()
		#print i, line, f,len(x1),len(xpix)
		xpix[i]=float(f[0])
		ypix[i]=float(f[1])
		i=i+1
	outfile.close()
	buffer=0.
	xmin=1.
	iraf.imgets(image=mipsimage,param='naxis1')#get RA of image
	xmax=float(iraf.imgets.value)
	ymin=1.
	iraf.imgets(image=mipsimage,param='naxis2')#get RA of image
	ymax=float(iraf.imgets.value)

	for i in range(len(x1)):
		if (xpix[i] < (1.+buffer)):
			self.matchflagediscs24[i]=-1.
			continue
		if (xpix[i] > (xmax-buffer)):
			self.matchflagediscs24[i]=-1.
			continue
		if (ypix[i] < (1.+buffer)):
			self.matchflagediscs24[i]=-1.
			continue
		if (ypix[i] > (ymax-buffer)):
			self.matchflagediscs24[i]=-1.
			continue
		(self.matchediscs24[i],self.matchflagediscs24[i],self.nmatchediscs24[i]) = findnearest(x1[i],y1[i],x2,y2,delta)

    def readirac(self,prefix):
	self.matchediscsirac = N.zeros(len(self.ra), 'i')
	self.matchflagediscsirac = N.zeros(len(self.ra), 'i')
	self.nmatchediscsirac = N.zeros(len(self.ra), 'i')

        print 'Reading IRAC data'
        file='/Users/rfinn/clusters/spitzer/irac-catalogs/'+str(prefix)+'.irac.cat'
	try:
		input=open(file,'r')
		ngal=0        #get number of galaxies

		for line in input:
			if line.find('#') > -1: #skip lines with '#' in them
				continue
			if line.find('\\') > -1: #skip lines with '#' in them
				continue
			if line.find('|') > -1: #skip lines with '#' in them
				continue
			ngal=ngal+1
		input.close()


		self.iracra = N.zeros(ngal,'f')#irac ra
		self.iracdec = N.zeros(ngal,'f')#irac dec
		self.iracf1 = N.zeros(ngal,'f')#irac flux channel 1
		self.erriracf1 = N.zeros(ngal,'f')
		self.iracf2  = N.zeros(ngal,'f')#irac flux channel 2
		self.erriracf2 = N.zeros(ngal,'f')
		self.iracf3 = N.zeros(ngal,'f')#irac flux channel 3
		self.erriracf3 = N.zeros(ngal,'f')
		self.iracf4 = N.zeros(ngal,'f')#irac flux channel 4
		self.erriracf4 = N.zeros(ngal,'f')
		self.iracsexflag0 = N.zeros(ngal,'f')
		self.iracsexflag1 = N.zeros(ngal,'f')
		self.iracwch1 = N.zeros(ngal,'f')
		self.iracwch2 = N.zeros(ngal,'f')
		self.iracwch3 = N.zeros(ngal,'f')
		self.iracwch4 = N.zeros(ngal,'f')	
		self.iracwmin = N.zeros(ngal,'f')#relative exposure on irac image

		print "IRAC file = ",file
		input=open(file,'r')
		i=0
		for line in input:
			if line.find('#') > -1: #skip lines with '#' in them
				continue
			if line.find('\\') > -1: #skip lines with '#' in them
				continue
			if line.find('|') > -1: #skip lines with '#' in them
				continue

			t=line.split()
			(self.iracra[i],self.iracdec[i],self.iracf1[i],self.iracf2[i],self.iracf3[i],self.iracf4[i],self.erriracf1[i],self.erriracf2[i],self.erriracf3[i],self.erriracf4[i],self.iracsexflag0[i],self.iracsexflag1[i],self.iracwch1[i],self.iracwch2[i],self.iracwch3[i],self.iracwch4[i],self.iracwmin[i])=(float(t[1]),float(t[2]),float(t[5]),float(t[6]),float(t[7]),float(t[8]),float(t[9]),float(t[10]),float(t[11]),float(t[12]),float(t[13]),float(t[14]),float(t[15]),float(t[16]),float(t[17]),float(t[18]),float(t[19]))
			i=i+1
		input.close()

		x1=self.ra
		y1=self.dec
		x2=self.iracra
		y2=self.iracdec
	

		for i in range(len(x1)):
			(self.matchediscsirac[i],self.matchflagediscsirac[i],self.nmatchediscsirac[i]) = findnearest(x1[i],y1[i],x2,y2,delta)

	except IOError:
		print "No IRAC data yet!"
		ngal=len(self.ra)
		self.iracra = N.zeros(ngal,'f')#irac ra
		self.iracdec = N.zeros(ngal,'f')#irac dec
		self.iracf1 = N.zeros(ngal,'f')#irac flux channel 1
		self.erriracf1 = N.zeros(ngal,'f')
		self.iracf2  = N.zeros(ngal,'f')#irac flux channel 2
		self.erriracf2 = N.zeros(ngal,'f')
		self.iracf3 = N.zeros(ngal,'f')#irac flux channel 3
		self.erriracf3 = N.zeros(ngal,'f')
		self.iracf4 = N.zeros(ngal,'f')#irac flux channel 4
		self.erriracf4 = N.zeros(ngal,'f')
		self.iracsexflag0 = N.zeros(ngal,'f')
		self.iracsexflag1 = N.zeros(ngal,'f')
		self.iracwch1 = N.zeros(ngal,'f')
		self.iracwch2 = N.zeros(ngal,'f')
		self.iracwch3 = N.zeros(ngal,'f')
		self.iracwch4 = N.zeros(ngal,'f')	
		self.iracwmin = N.zeros(ngal,'f')#relative exposure on irac image

	#delta=2. #max allowed offset (in arcseconds) between matched sources


    def readwide24(self,file,file2):
        input=open(file,'r')
        #get number of galaxies
        ngal=0
	for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue
            ngal=ngal+1
        input.close()


        self.idwide24 = N.zeros(ngal,'f')
        self.imagexwide24 = N.zeros(ngal,'f')
        self.imageywide24  = N.zeros(ngal,'f')
        self.rawide24 = N.zeros(ngal,'f')
        self.decwide24 = N.zeros(ngal,'f')
        self.fwide24 = N.zeros(ngal,'f')#flux
        self.errfwide24 = N.zeros(ngal,'f')

        input=open(file,'r')
        i=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue

            t=line.split()

	    (self.idwide24[i],self.rawide24[i],self.decwide24[i],self.imagexwide24[i],self.imageywide24[i],self.fwide24[i],snr)=(float(t[0]),float(t[1]),float(t[2]),float(t[3]),float(t[4]),float(t[5]),float(t[6]))
	    #print "Number of fields in first line = ",len(t)," expecting 5"
	    self.errfwide24[i]=self.fwide24[i]/snr
	    #if i < 10:
	    #	    print line
	#	    for ta in t:
	#		    print ta
	#	    print i, self.fwide24[i]
		
            i=i+1
        input.close()

	delta=1. #max allowed offset (in arcseconds) between matched sources
	delta=delta/3600.
	x1=self.ra24
	y1=self.dec24
	x2=self.rawide24
	y2=self.decwide24
	
	self.matchwide24 = N.zeros(len(x1), 'i')
	self.nmatchwide24 = N.zeros(len(x1), 'i')
	self.matchflagwide24 = N.zeros(len(x1), 'i')
	j=0
	for i in range(len(x1)):
		(self.matchwide24[i],self.matchflagwide24[i],self.nmatchwide24[i]) = findnearest(x1[i],y1[i],x2,y2,delta)
		if self.matchflagwide24[i] > 0:
			if j < 10:
				print j, self.ra24[i],self.dec24[i],self.rawide24[self.matchwide24[i]],self.decwide24[self.matchwide24[i]]
				print j, self.f24[i], self.fwide24[self.matchwide24[i]]
				j=j+1
	i=0
	fw=[]
	fwerr=[]
	f=[]
	ferr=[]
	for i in range(len(self.ra24)):
	    if (self.matchflagwide24[i] > 0):
		f.append(self.f24[i])
		ferr.append(self.errf24[i])
		fw.append(self.fwide24[self.matchwide24[i]])
		fwerr.append(self.fwide24[self.matchwide24[i]])

	f=N.array(f,'d')
	ferr=N.array(ferr,'d')
	fw=N.array(fw,'d')
	fwerr=N.array(fwerr,'d')
	print "average ratio of core/wide flux = ",N.average(f/fw)," +/-",pylab.std(f/fw)
	pylab.plot(f,fw,'y^',label="",markersize=5)
	pylab.xlabel('F(24)-core')
	pylab.ylabel('F(24)-wide')
	#pylab.show()
	pylab.clf()
	pylab.cla()
	pylab.plot(f,fw,'b^',label="Scan vs Phot",markersize=5)

	pylab.xlabel('F(24)-core (uJy)')
	pylab.ylabel('F(24)-wide (uJy)')
	x=N.arange(0.,(max(f)+100.),100.)
	y=x
	pylab.plot(x,y,'k',label='y=x')
	y=0.7*x
	pylab.plot(x,y,'k--',label="y=0.7 x")
	t=pylab.legend(loc='upper left')				
	pylab.axis([100.,300.,100.,300.])
	d=N.compress(f < 200.,f)
	w=N.compress(f < 200.,fw)
	r=(w-d)/d
	print "STD for f < 250 uJy (wide-deep)/deep = ",pylab.std(r)
	#pylab.errorbar(f,fw,yerr=fwerr,xerr=ferr,fmt=None,ecolor='k')
	pylab.savefig(file2)

    def magnitudetofluxconversion(self,prefix,plothaflag):
	    i=0
	    pylab.cla()
	    pylab.clf()
	    log=open('Klog','w')
	    for i in range(len(self.magK)):
		    self.fluxK[i] = 640/10**(self.magK[i]/2.5)
		    print >>log, i, self.magK[i],self.fluxK[i]
		    i=i+1
	    log.close()
	    ra=[]
	    dec=[]
	    xbad = N.compress(self.magK < -90.,self.ra)
	    ybad = N.compress(self.magK < -90.,self.dec)
	    pylab.plot(self.ra,self.dec,'ko',label="full phot",markersize=3)
	    pylab.plot(xbad,ybad,'ro',label="No NIR",markersize=3)
	    pylab.plot(self.ra24,self.dec24,'bs',label="24um",markersize=7)
	    if (plothaflag > 0.1):
		    pylab.plot(self.raha,self.decha,'c.',label="Halpha",markersize=3)

	    i=0
	    for i in range(len(self.ra)):
		    if (self.matchflaghaediscs[i] > 0):
			    if (self.matchflagediscs24[i] > 0):
				    ra.append(self.ra[i])
				    dec.append(self.dec[i])
				    #print i,"got a match"
		    i=i+1

	    ra=N.array(ra,'d')
	    dec=N.array(dec,'d')
	    pylab.plot(ra,dec,'r^',label="24-Ha-Edi",markersize=9)
	    t=pylab.legend(loc='upper left')					       
	    outfile=str(prefix)+'.eps'
	    pylab.savefig(outfile)
		    
    def readsextractor(self,file):#read sextractor output
        input=open(file,'r')
        #get number of galaxies
        ngal=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue
            ngal=ngal+1
        input.close()

        self.x2 = N.zeros(ngal,'f')
        self.y2  = N.zeros(ngal,'f')
        self.f2 = N.zeros(ngal,'f')#flux
        self.errf2 = N.zeros(ngal,'f')

        input=open(file,'r')
        i=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue

            t=line.split()
	    self.f2[i]=float(t[8])
	    self.errf2[i]=float(t[9])
	    #self.f2[i]=float(t[35])
	    #self.errf2[i]=float(t[42])
	    self.x2[i]=float(t[10])
	    self.y2[i]=float(t[11])	

            i=i+1
        input.close()


    def readhafile(self,file):#read halpha output and match to ediscs
	input=open(file,'r')
	#get number of galaxies
	ngal=0
	for line in input:
	    if line.find('#') > -1:
		continue
	    ngal=ngal+1
	input.close()
		
	self.haID = []
	self.raha = N.zeros(ngal,'f')
	self.decha = N.zeros(ngal,'f')
	#self.Fluxn = N.zeros(ngal,'f')
	#self.Fluxnerr = N.zeros(ngal,'f')
	#self.FluxJ = N.zeros(ngal,'f')
	#self.FluxJerr = N.zeros(ngal,'f')
	#self.magJ = N.zeros(ngal,'f')
	#self.magJerr = N.zeros(ngal,'f')
	#self.ContSub = N.zeros(ngal,'f')
	#self.ContSuberr = N.zeros(ngal,'f')
	#self.SNR = N.zeros(ngal,'f')
	self.EWha = N.zeros(ngal,'f')
	self.EWhaerr = N.zeros(ngal,'f')
	self.SFR = N.zeros(ngal,'d')
	self.SFRerr = N.zeros(ngal,'d')
	self.lha = N.zeros(ngal,'d')
	self.lhaerr = N.zeros(ngal,'d')
	self.SF_flag = N.zeros(ngal,'f')

	input=open(file,'r')
	#get number of galaxies
	i=0
	for line in input:
		if line.find('#') > -1:
			continue
		t=line.split()
		shortname=t[0]
		self.haID.append(shortname)
		self.raha[i]=((float(shortname[0:2])+float(shortname[2:4])/60.+(float(shortname[4:7])/10.)/3660.)*15.)#convert from hours to degrees
		
		if float(shortname[7:10]) < 0:
			self.decha[i]=float(shortname[7:10])-float(shortname[10:12])/60.-(float(shortname[12:])/10.-0.)/3660#to align w/mips
		else:
			self.decha=float(shortname[7:10])+float(shortname[10:12])/60.+float(shortname[12:])/10./3660
		self.SFR[i]=float(t[14])
		self.SFRerr[i]=float(t[15])
		self.EWha[i]=float(t[12])
		self.EWhaerr[i] = float(t[13])
		self.SF_flag[i] = float(t[16])
		i=i+1
	input.close()
       
	self.matchhaediscs = N.zeros(len(self.ediscsID), 'i')
	self.matchflaghaediscs = -1.*N.ones(len(self.ediscsID), 'i')
	i=0
	for i in range(len(self.ediscsID)):
		j=0
		ediscsname=str(self.ediscsID[i])
		    #print "checking ha id \n"
		for j in range(len(self.haID)):
			haname=str(self.haID[j])
			if ediscsname.find(haname) > -1:
				#print "matched Halpha name to ediscs ID \n",(self.raha[j]-self.ra[i])*3600./15.,(self.decha[j]-self.dec[i])*3600.
				self.matchhaediscs[i]=j
				self.matchflaghaediscs[i]=1
				self.raha[j]=self.ra[i]
				self.decha[j]=self.dec[i]
				if (self.EWha[j] < 0.):
					self.matchflaghaediscs[i]=0.
				break
	self.lha=self.SFR/(7.9e-42)
	self.lhaerr=self.SFRerr/7.9e-42


    def readha2file(self,file):#read h alpha output of cl0023
	    input=open(file,'r')
	#get number of galaxies
	    ngal=0
	    for line in input:
		    if line.find('#') > -1:
			    continue
		    ngal=ngal+1
		    input.close()
		    racenter=(0.+23.0/60.+51.8/60./60.)*15.
		    deccenter=4.+22./60.+41./60./60.
		    self.raha = N.zeros(ngal,'f')
		    self.decha = N.zeros(ngal,'f')
	#self.x = N.zeros(ngal,'f')
	#self.y = N.zeros(ngal,'f')
	#self.Fluxn = N.zeros(ngal,'f')
	#self.Fluxnerr = N.zeros(ngal,'f')
	#self.FluxJ = N.zeros(ngal,'f')
	#self.FluxJerr = N.zeros(ngal,'f')
	#self.Contsub = N.zeros(ngal,'f')
	#self.Contsuberr = N.zeros(ngal,'f')
	#self.EW(Ha) = N.zeros(ngal,'f')
	#self.EW(Ha)err = N.zeros(ngal,'f')
	#self.Flux(Ha) =  N.zeros(ngal,'f')
	#self.Flux(Ha)err = N.zeros(ngal,'f')
	#self.L(Ha) = N.zeros(ngal,'f')
	#self.L(Ha)err = N.zeros(ngal,'f')
		    self.SFR = N.zeros(ngal,'f')
		    self.SFRerr = N.zeros(ngal,'f')
		    
	    input=open(file,'r')
	    i=0
	    for line in input:
		    if line.find('#') > -1: #skip lines with '#' in them
			    continue
		    if line.find('\\') > -1: #skip lines with '#' in them
			    continue
		    if line.find('|') > -1: #skip lines with '#' in them
			    continue
			    
		    t=line.split()
		    dra=float(t[1])/60./60.
		    dec=float(t[2])/60./60.
		    self.raha[i]=racenter+dra
		    self.decha[i]=deccenter+dec
		    self.SFR[i]=float(t[17])
		    self.SFRerr[i]=float(t[18])
			    
		    i=i+1
		    input.close()


    def readphotoz(self,file):
        input=open(file,'r')
	#get number of galaxies
	ngal=0
	for line in input:
	    if line.find('#') > -1:
		continue
	    ngal=ngal+1
	input.close()

        self.idnew =[]
	self.idold =[]
	self.bestz = N.zeros(ngal,'f')
 	self.membflag = N.zeros(ngal,'f')
	self.lowz = N.zeros(ngal,'f')
        self.highz = N.zeros(ngal,'f')
        self.wmin = N.zeros(ngal,'f')
        self.Pclust = N.zeros(ngal,'f')
        self.LUlowzclust = N.zeros(ngal,'f')
        self.LUbestzclust = N.zeros(ngal,'f')
        self.LUhighzclust = N.zeros(ngal,'f')
        self.LBlowzclust = N.zeros(ngal,'f')
        self.LBbestzclust = N.zeros(ngal,'f')
        self.LBhighzclust = N.zeros(ngal,'f')
        self.LVlowzclust = N.zeros(ngal,'f')
        self.LVbestzclust = N.zeros(ngal,'f')
        self.LVhighzclust = N.zeros(ngal,'f')
        self.LRlowzclust = N.zeros(ngal,'f')
        self.LRbestzclust = N.zeros(ngal,'f')
        self.LRhighzclust = N.zeros(ngal,'f')
        self.LIlowzclust = N.zeros(ngal,'f')
        self.LIbestzclust = N.zeros(ngal,'f')
        self.LIhighzclust = N.zeros(ngal,'f')
        self.UBlowzclust = N.zeros(ngal,'f')
        self.UBbestzclust = N.zeros(ngal,'f')
        self.UBhighzclust = N.zeros(ngal,'f')
        self.BVlowzclust = N.zeros(ngal,'f')
        self.BVbestzclust = N.zeros(ngal,'f')
        self.BVhighzclust = N.zeros(ngal,'f')
        self.UVlowzclust = N.zeros(ngal,'f')
        self.UVbestzclust = N.zeros(ngal,'f')
        self.UVhighzclust = N.zeros(ngal,'f')
        self.Membflag = N.zeros(ngal,'f')
        self.LJlowzclust = N.zeros(ngal,'f')
        self.LJbestzclust = N.zeros(ngal,'f')
        self.LJhighzclust = N.zeros(ngal,'f')
        self.LKlowzclust = N.zeros(ngal,'f')
        self.LKbestzclust = N.zeros(ngal,'f')
        self.LKhighzclust = N.zeros(ngal,'f')
        self.xcorr = N.zeros(ngal,'f')
        self.ycorr = N.zeros(ngal,'f')
        self.starflag = N.zeros(ngal,'f')

	input=open(file,'r')

        i=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue

	    t=line.split()
	    self.idnew.append(t[0])
	    self.idold.append(t[1])
	    self.bestz[i] = float(t[3])
	    self.membflag[i]= float(t[56])
	 
	    self.lowz[i] = float(t[2])
            self.highz[i] = float(t[4])
            self.wmin[i] = float(t[30])
            self.Pclust[i] = float(t[31])
            self.LUlowzclust[i] = float(t[32])
            self.LUbestzclust[i] = float(t[33])
            self.LUhighzclust[i] = float(t[34])
            self.LBlowzclust[i] = float(t[35])
            self.LBbestzclust[i] = float(t[36])
            self.LBhighzclust[i] = float(t[37])
            self.LVlowzclust[i] = float(t[38])
            self.LVbestzclust[i] = float(t[39])
            self.LVhighzclust[i] = float(t[40])
            self.LRlowzclust[i] = float(t[41])
            self.LRbestzclust[i] = float(t[42])
            self.LRhighzclust[i] = float(t[43])
            self.LIlowzclust[i] = float(t[44])
            self.LIbestzclust[i] = float(t[45])
            self.LIhighzclust[i] = float(t[46])
            self.UBlowzclust[i] = float(t[47])
            self.UBbestzclust[i] = float(t[48])
            self.UBhighzclust[i] = float(t[49])
            self.BVlowzclust[i] = float(t[50])
            self.BVbestzclust[i] = float(t[51])
            self.BVhighzclust[i] = float(t[52])
            self.UVlowzclust[i] = float(t[53])
            self.UVbestzclust[i] = float(t[54])
            self.UVhighzclust[i] = float(t[55])
            self.LJlowzclust[i] = float(t[57])
            self.LJbestzclust[i] = float(t[58])
            self.LJhighzclust[i] = float(t[59])
            self.LKlowzclust[i] = float(t[60])
            self.LKbestzclust[i] = float(t[61])
            self.LKhighzclust[i] = float(t[62])
            self.xcorr[i] = float(t[102])
            self.ycorr[i] = float(t[103])
            self.starflag[i] = float(t[104])
	    i = i+1
	    
	input.close()


    def match24sextractor(self,file):
        input=open(file,'r')
	#get number of galaxies
	ngal=0
	for line in input:
	    if line.find('#') > -1:
		continue
	    ngal=ngal+1
	input.close()

    def readgimmorphology(self,file1,file2): #read gim morphology and match to ediscs
        input=open(file1,'r')
	#get number of galaxies
	ngal=0
	for line in input:
	    if line.find('#') > -1:
		continue
	    ngal=ngal+1
	input.close()

	self.hstgimID = [] #id for gim
	self.hstgimstargal = N.zeros(ngal,'f') #star(0) or galaxy (1)
	self.hstgimtype = N.zeros(ngal,'f') #galaxy type

	input=open(file1,'r') #opens hstgim id file
        i=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue   

	    t=line.split()
	    self.hstgimID.append(t[0])

	    i=i+1
	input.close()

	input=open(file2,'r') #opens hstgim morphology data
        i=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue   

	    t=line.split()
	    self.hstgimstargal[i]=float(t[5])
	    self.hstgimtype[i]=float(t[12])
	    
	    i=i+1
	input.close()

        self.matchmorphgimediscs = N.zeros(len(self.ediscsID), 'i') #match morphology to ediscs
        self.matchflagmorphgimediscs = N.zeros(len(self.ediscsID), 'i')
	i=0
	for i in range(len(self.ediscsID)):
		j=0
		ediscsname=str(self.ediscsID[i])
		    #print "checking ha id \n"
		for j in range(len(self.hstgimID)):
			hstgimname=str(self.hstgimID[j])
			if ediscsname.find(hstgimname) > -1:
				self.matchmorphgimediscs[i]=j
				self.matchflagmorphgimediscs[i]=1
				break



    def readvismorphology(self,file): #read vis morphology and match to ediscs
        input=open(file,'r')
        #get number of galaxies
        ngal=0
        for line in input:
            if line.find('#') > -1:
                continue
            ngal=ngal+1
        input.close()

        self.hstvisID = []
        self.hstvisra = N.zeros(ngal,'f')
        self.hstvisdec = N.zeros(ngal,'f')
        self.hstvisIauto =  N.zeros(ngal,'f')
        self.hstvisnumtype = N.zeros(ngal,'f')
        self.hstvisflagfeatureless = N.zeros(ngal,'f')
        self.hstvisflagedgeon = N.zeros(ngal,'f')
        self.hstvisflagtoosmall = N.zeros(ngal,'f')
        self.hstvisflagLSB = N.zeros(ngal,'f')
        self.hstvisflagcosray = N.zeros(ngal,'f')
        self.hstvisflagdust = N.zeros(ngal,'f')
        self.hstvisflagdistorted = N.zeros(ngal,'f')


        input=open(file,'r')
        i=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue  

            t=line.split()
            self.hstvisID.append(t[0])
            self.hstvisra[i] = float(t[1])
            self.hstvisdec[i] = float(t[2])
            self.hstvisIauto[i] = float(t[3])
            self.hstvisnumtype[i] = float(t[4])
            self.hstvisflagfeatureless[i] = float(t[5])
            self.hstvisflagedgeon[i] = float(t[6])
            self.hstvisflagtoosmall[i] = float(t[7])
            self.hstvisflagLSB[i] = float(t[8])
            self.hstvisflagcosray[i] = float(t[9])
            self.hstvisflagdust[i] = float(t[10])
            #self.hstvisdistorted[i] = float(t[11])

            i = i+1
        input.close()
	
        self.matchmorphvisediscs = N.zeros(len(self.ediscsID), 'i') #match morphology to ediscs
        self.matchflagmorphvisediscs = N.zeros(len(self.ediscsID), 'i')
	i=0
	for i in range(len(self.ediscsID)):
		j=0
		ediscsname=str(self.ediscsID[i])
		    #print "checking ha id \n"
		for j in range(len(self.hstvisID)):
			hstvisname=str(self.hstvisID[j])
			if ediscsname.find(hstvisname) > -1:
				self.matchmorphvisediscs[i]=j
				self.matchflagmorphvisediscs[i]=1
				break


    def readspectroscopy(self,file):
        input=open(file,'r')
        #get number of galaxies
	ngal=0
	for line in input:
		if line.find('#') > -1:
			continue
		if line.find('\\') > -1: #skip lines with '#' in them
			continue
		if line.find('|') > -1: #skip lines with '#' in them
			continue
		if line[0:4].find('c') < 0: #skip lines with out 'cl1040' in them
			continue
		ngal=ngal+1
       	input.close()
	
	self.specediscsID = []
	self.specmembership = N.zeros(ngal,'i')
	self.specz = N.zeros(ngal,'f')
	self.spectype = N.zeros(ngal,'f')
	self.specEWOII = N.zeros(ngal,'f')
	self.specEWOIIflag = N.zeros(ngal,'f')

		
		
	input=open(file,'r')
	i=0
	for line in input:
		if line.find('#') > -1: #skip lines with '#' in them
			continue
		if line.find('\\') > -1: #skip lines with '#' in them
			continue
		if line.find('|') > -1: #skip lines with '#' in them
			continue
		if line[0:4].find('c') < 0: #skip lines with out 'cl1040' in them
			continue
		t = line.split()
		name=t[4]
		if name.find(':') > -1:
			self.specediscsID.append(name[:(len(name)-1)])
		else:
			self.specediscsID.append(name)
		try:
			self.specz[i]=float(t[6])
		except:
			self.specz[i]=-99.
		#memb=t[3]
		try:
			self.specmembership[i] = float(t[3])
		except:
			self.specmembership[i] = 0.
		try:
			self.spectype[i] = float(t[9])
		except:
			self.spectype[i] = -99.
		try:			
			self.specEWOII[i]=float(line[80:85])
			self.specEWOIIflag[i]=1.
		except:
			self.specEWOII[i]=-99.	
		i = i+1
       	input.close()

   

	self.matchspecediscs = N.zeros(len(self.ediscsID),'i')
	self.matchflagspecediscs = N.zeros(len(self.ediscsID),'i')

	i=0
	for i in range(len(self.idold)):
		j=0
		ediscsoldname=str(self.idold[i])
		    #print "checking ha id \n"
		for j in range(len(self.specediscsID)):
			specname=str(self.specediscsID[j])
			if ediscsoldname.find(specname) > -1:
				try:
					self.matchspecediscs[i]=j
					self.matchflagspecediscs[i]=1
					break
				except IndexError:
					print "Index error in match spec"
					print "len(self.matchspecediscs)",len(self.matchspecediscs),len(self.ediscsID)
					print "i, j",i,j
    def readspectroscopyv2(self,file):
	    #update this for reading v24

	    input=open(file,'r')
         #get number of galaxies
	    ngal=0
	    for line in input:
		    if line.find('#') > -1:
			    continue
		    if line.find('\\') > -1: #skip lines with '#' in them
			    continue
		    if line.find('|') > -1: #skip lines with '#' in them
			    continue
		    if line[0:4].find('cl') < 0: #skip lines with out 'cl1040' in them
			    continue
		    ngal=ngal+1
	    input.close()
       
	    self.specediscsID = []
	    self.specmembership = N.zeros(ngal,'i')
	    self.specz = N.zeros(ngal,'f')
	    self.spectype = N.zeros(ngal,'f')
	    self.specEWOII = N.zeros(ngal,'f')
	    self.specEWOIIflag = N.zeros(ngal,'f')

               
               
	    input=open(file,'r')
	    i=0
	    for line in input:
		    if line.find('#') > -1: #skip lines with '#' in them
			    continue
		    if line.find('\\') > -1: #skip lines with '#' in them
			    continue
		    if line.find('|') > -1: #skip lines with '#' in them
			    continue
		    if line[0:4].find('c') < 0: #skip lines with out 'cl1040' in them
			    continue
		    t = line.split()
		    name=t[3]
		    if name.find(':') > -1:
			    self.specediscsID.append(name[:(len(name)-1)])
		    else:
			    self.specediscsID.append(name)
		    try:
			    self.specz[i]=float(t[9])
		    except:
			    self.specz[i]=-99.
		    memb=t[8]
		    try:
			    memb=int(memb)
			    self.specmembership[i]=memb
		    except:
			    self.specmembership[i]=0.
		    try:
			    self.spectype[i] = float(t[10])
		    except:
			    self.spectype[i] = -99.
		    self.specEWOII[i]=-99. 
		    i = i+1
	    input.close()


	    self.matchspecediscs = N.zeros(len(self.ediscsID),'i')
	    self.matchflagspecediscs = N.zeros(len(self.ediscsID),'i')

	    i=0
	    for i in range(len(self.idold)):
                j=0
                ediscsoldname=str(self.idold[i])
                    #print "checking ha id \n"
                for j in range(len(self.specediscsID)):
                        specname=str(self.specediscsID[j])
                        if ediscsoldname.find(specname) > -1:
                                self.matchspecediscs[i]=j
                                self.matchflagspecediscs[i]=1
                                break
	    print "Number of matches with spectroscopy = ",N.sum(self.matchflagspecediscs)


    def memb(self):
	    nspec=0.
	    nspecphotagree=0
	    nspecphotdisagree=0

	    self.defmembership=N.zeros(self.nediscs,'i')
	    self.newspecmatchflag=N.zeros(self.nediscs,'i')
	    for i in range(len(self.ediscsID)):
		    if (self.membflag[i] > 0) & (self.matchflagspecediscs[i] < 1.):
			    self.defmembership[i] = 1
			    continue
		    if (self.matchflagspecediscs[i] > 0.):
			    #print "got here",self.matchflagediscs24[i]
			    if (self.matchflagediscs24[i] > 0.):
				    nspec=nspec + 1
			    j=self.matchspecediscs[i]
			    diff=abs(self.specmembership[j]-self.membflag[i])
			    if (self.matchflagediscs24[i] > 0.):
				    if (diff < .1):
					    nspecphotagree=nspecphotagree + 1.
				    else:
					    nspecphotdisagree=nspecphotdisagree + 1.
			    if self.specmembership[j] > 0:
				    self.defmembership[i]=1.
			            continue
		    #print i, "defmembership = ",self.defmembership[i]
		    

	    for i in range(len(self.newspecmatchflag)):
		    if (self.matchflagspecediscs[i] > 0.):
			    if (self.specmembership[self.matchspecediscs[i]] > 0.):
				    self.newspecmatchflag[i] = 1.
				    continue
		    if (self.matchflagspecediscs[i] > 0.):
			    if (self.specmembership[self.matchspecediscs[i]] < 1.):
				    self.newspecmatchflag[i] = 0.
				    continue
		    if (self.matchflagspecediscs[i] < 1.):
			    self.newspecmatchflag[i] = -1.

	    print "24 micron sources only:  Nspec, nphotspecagree, %, nphotspecdisagree, %"
	    try:
		    print nspec, nspecphotagree, nspecphotagree/(1.*nspec), nspecphotdisagree,nspecphotdisagree/(1.*nspec)
	    except ZeroDivisionError:
		    print "Warning:  nspec = 0!!!"
	    nspecmemb=len(N.compress((self.matchflagediscs24>0.) & (self.newspecmatchflag > 0.),self.newspecmatchflag))
	    nspecmembNphotmemb=len(N.compress((self.matchflagediscs24>0.) & (self.newspecmatchflag > 0.) & (self.membflag > 0.),self.newspecmatchflag))
	    (a,b)=ratioerror(1.*nspecmembNphotmemb,1.*nspecmemb)
	    print "Fraction of spec members & 24micron that are photoz members = %5.3f +/- %5.3f"%(a,b)
	    nspecmemb=len(N.compress((self.matchflagediscs24>0.) & (abs(self.newspecmatchflag) < .1) ,self.newspecmatchflag))
	    nspecmembNphotmemb=len(N.compress((self.matchflagediscs24>0.) & (abs(self.newspecmatchflag) < .1) & (self.membflag > 0.),self.newspecmatchflag))
	    (a,b)=ratioerror(1.*nspecmembNphotmemb,1.*nspecmemb)
	    print "Fraction of spec non-members & 24 micron that are photoz members = %5.3f +/- %5.3f"%(a,b)," Nspecmemb = ",nspecmemb

	    nspecmemb=len(N.compress( (abs(self.newspecmatchflag) < .1) ,self.newspecmatchflag))
	    nspecmembNphotmemb=len(N.compress( (abs(self.newspecmatchflag) < .1) & (self.membflag > 0.),self.newspecmatchflag))
	    (a,b)=ratioerror(1.*nspecmembNphotmemb,1.*nspecmemb)
	    print "Fraction of ALL spec non-members that are photoz members = %5.3f +/- %5.3f"%(a,b)," Nspec nonmemb = ",nspecmemb
    def readflat24(self,file):
            input=open(file,'r')
        #get number of galaxies
            ngal=0
            for line in input:
                    if line.find('#') > -1: #skip lines with '#' in them
                            continue
                    if line.find('\\') > -1: #skip lines with '#' in them
                            continue
                    if line.find('|') > -1: #skip lines with '#' in them
                            continue
                    ngal=ngal+1
            input.close()
                   
                   
            self.idflat24 = N.zeros(ngal,'f')
            self.imagexflat24 = N.zeros(ngal,'f')
            self.imageyflat24  = N.zeros(ngal,'f')
            self.raflat24 = N.zeros(ngal,'f')
            self.decflat24 = N.zeros(ngal,'f')
            self.fflat24 = N.zeros(ngal,'f') #flux
            self.errfflat24 = N.zeros(ngal,'f')
                   
            input=open(file,'r')
            i=0
            for line in input:
                    if line.find('#') > -1: #skip lines with '#' in them
                            continue
                    if line.find('\\') > -1: #skip lines with '#' in them
                            continue
                    if line.find('|') > -1: #skip lines with '#' in them
                            continue
                           
                    t=line.split()
                    (self.idflat24[i],self.imagexflat24[i],self.imageyflat24[i],self.raflat24[i],self.decflat24[i],self.fflat24[i],self.errfflat24[i])=(float(t[0]),float(t[1]),float(t[2]),float(t[3]),float(t[4]),float(t[5]),float(t[6]))
                    i=i+1
            input.close()
            outfile=open('xyflat24.dat','w')
            for i in range(len(self.imagex24)):
            #print self.imagex24[i],self.imagey24[i],self.f24[i],self.errf24[i]
                    string="%6.2f %6.2f %8.1f %6.2f \n" % (self.imagex24[i],self.imagey24[i],self.f24[i],self.errf24[i])
                    outfile.write(string)
            outfile.close()
                                   
            #delta=3.#max allowed offset (in arcseconds) between matched sources
            #delta=delta/3600.
            x1=self.ra
            y1=self.dec
            x2=self.raflat24
            y2=self.decflat24
                                   
            self.matchediscsflat24 = N.zeros(len(x1), 'i')
            self.nmatchediscsflat24 = N.zeros(len(x1), 'i')
            self.matchflagediscsflat24 = N.ones(len(x1), 'i')
            for i in range(len(x1)):
                    (self.matchediscsflat24[i],self.matchflagediscsflat24[i],self.nmatchediscsflat24[i]) = findnearest(x1[i],y1[i],x2,y2,delta)

    def mastertable(self,name):
	    output=open('mastertable.fmt','w')
	    output.write('#1 - ediscsID \n')
	    output.write('#2 - ediscs RA \n')
	    output.write('#3 - ediscs Dec \n')
	    output.write('#4 - Halpha EW (A) \n')
	    output.write('#5 - err Halpha EW (A) \n')
	    output.write('#6 - Halpha SFR (Msol/yr, 1 mag extinct) \n')
	    output.write('#7 - Halpha err SFR (Msol/yr, 1 mag extinct) \n')
	    output.write('#8 - match flag Halpha: -1, 0, 1  \n')
	    output.write('#9 - GIM2D Type of HST image \n')
	    output.write('#10 - matchflag GIM2D Type of HST image \n')
	    output.write('#11 - HST visual morph \n')
	    output.write('#12 - HST visual morph matchflag\n')
	    output.write('#13 - 24 um total flux (uJy) \n')
	    output.write('#14 - error in 24 um total flux (uJy) \n')
	    output.write('#15 - matchflag 24 um \n')
	    output.write('#16-25 - isophotal mag and err: V,R,I,J,K \n')
	    output.write('#26-35 - 1" aper mag and err: V,R,I,J,K \n')
	    output.write('#36 - photoz memb flag \n')
	    output.write('#37 - spec membership flag: -1, 0, 1\n')
	    output.write('#38 - photoz + spec membership flag: -0, 1\n')
	    output.write('#39 - spec redshift\n')
	    output.write('#40 - spec type:1=abs-line; 2=abs-line+weak emission; 3=emission w/EW(OII)<25; 4=strong emission w/EW(OII)\n')
	    output.write('#41 - EW OII (A)\n')
	    output.write('#42 - match flag spectroscopy\n')
	    output.write('#43 - best photoz \n')
	    output.write('#44 - low photoz \n')
	    output.write('#45 - high photoz \n')
	    output.write('#46 - wmin: fraction of total exposure in near-IR\n')
	    output.write('#47 - pclust: probability of member based on photozs \n')
	    output.write('#48-62 - low, best, high rest-frame luminosity in U,B,V,R,I \n')
	    output.write('#63-71 - low, best, high rest-frame colors  U-B,B-V,U-V \n')
	    output.write('#72-77 - low, best, high rest-frame luminosity in J,K \n')
	    output.write('#78 - x position on geo corrected vlt image \n')
	    output.write('#79 - y position on geo corrected vlt image \n')
	    output.write('#80 - starflag \n')
	    output.write('#81 - K-band flux \n')
	    output.write('#82 - Star-forming flag based on Halpha criteria \n')
	    output.write('#83 - match flag for ediscs spectroscopy \n')
	    output.write('#84 - EW OII flag for ediscs spectroscopy \n')
	    output.write('#85 - number of optical ediscs galaxies within 2arcsec of 24um source \n')
	    output.write('#86 - matchflagediscsirac\n') 
	    output.write('#87-90 - iracf1 iracf2 iracf3 iracf4 (uJy) \n')
	    output.write('#91-94 - erriracf1 erriracf2 erriracf3 erriracf4 (uJy) \n') 
	    output.write('#95 - iracsexflag0 \n') 
	    output.write('#96 - iracsexflag1 \n') 
	    output.write('#97-100 - iracwch1 iracwch2 iracwch3 iracwch4 \n') 
	    output.write('#101 - iracwmin \n') 
	    output.write('#102 - nmatchediscsirac: number of optical sources w/in 2arcsec radius of irac source\n') 
	    output.write('#103 - 24 micron luminosity (H0=70 km/s/Mpc)\n') 
	    output.write('#104 - error in 24 micron luminosity (H0=70 km/s/Mpc)\n') 
	    output.write('#105 - Halpha luminosity (H0=70 km/s/Mpc)\n') 
	    output.write('#106 - error in Halpha luminosity (H0=70 km/s/Mpc)\n') 
	    #output.write('#106 - error in Halpha luminosity (H0=70 km/s/Mpc)\n') 
	    output.close()
	    outfile=str(name)+'mastertable.dat'
	    output= open(outfile,'w')
	    outfile2=str(name)+'mastertable24.dat'
	    output2= open(outfile2,'w')
	    outfile3=str(name)+'AllMembers.dat'
	    output3=open(outfile3,'w')
	    outfile4=str(name)+'Members24.dat'
	    output4=open(outfile4,'w')
	    outfile5=str(name)+'ediscs24.dat'#tables for ediscs collaboration
	    output5=open(outfile5,'w')
            output5.write('#ediscsID     matchflag24  f24(uJy) errf24(uJy) Vmag(1arc) errVmag photmembflag specmembflag photo+specmemb wmin spectype OIIEW HaFlag HaEW HaEWerr HaSFR HaSFRerr\n')

            output5.write('#matchflag24: -1=not on 24um image; 0=on 24um image but no 24um detectioin; 1=24um source \n')
            output5.write('#NOTE: f24 is meaningful only if matchflag24 > 0 \n')
            output5.write('#photmembflag: 0 = not a member; 1 = a member according to photoz \n')
	    output5.write('#specmembflag: -1 = no spec; 0 = spec but not a member; 1 = spec and a member \n')
            output5.write('#phot+specmemb: 0 = if photoz or spec says not a member; 1 = if spec or photoz says is a member \n')
            output5.write('#wmin = fraction of full near-IR exposure.  photozs are unreliable for wmin < 0.3\n')
	    output.write('#ediscsID  ra  dec  xcorr  ycorr  starflag  EWha  EWhaerr  SFR  SFRerr  matchflaghaediscs SF_flag hstgimtype  matchflagmorphgimediscs  hstvisnumtype  matchflagmorphvisediscs  matchflagediscs24  f24  errf24  nmatchediscs24  misoV  misoeVapsim  misoR  misoeRapsim  misoI  misoeIapsim  misoJ  misoeJapsim  misoK  misoeKapsim   magV_1  mageVapsim  magR_1  mageRapsim  magI_1  mageIapsim  magJ_1  mageJapsim  magK_1  mageKapsim  membflag  newspecmatchflag defmembership  specz  spectype  specEWOII matchflagspecediscs specEWOIIflag  bestz   lowz   highz   wmin   Pclust  LUlowzclust  LUbestzclust  LUhighzclust  LBlowzclust  LBbestzclust  LBhighzclust  LVlowzclust  LVbestzclust  LVhighzclust  LRlowzclust  LRbestzclust  LRhighzclust  LIlowzclust  LIbestzclust  LIhighzclust  LJlowzclust  LJbestzclust  LJhighzclust  LKlowzclust  LKbestzclust  LKhighzclust  fluxK  UBlowzclust  UBbestzclust  UBhighzclust  BVlowzclust  BVbestzclust  BVhighzclust  UVlowzclust  UVbestzclust  UVhighzclust  matchflagediscsirac iracf1 iracf2  iracf3 iracf4 erriracf1 erriracf2 erriracf3 erriracf4 iracsexflag0 iracsexflag1 iracwch1 iracwch2 iracwch3 iracwch4 iracwmin nmatchediscsirac L24 errL24 LHa errLHa \n')
	    for i in range(len(self.ediscsID)):
		    try:
			    ii=self.matchediscsirac[i]
	                    ih=self.matchhaediscs[i]
			    i24=self.matchediscs24[i]
			    ispec=self.matchspecediscs[i]
			    string="%s %12.6f %12.6f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %8.4e %8.4e %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.4f %7.2f %7.2f %7.2f %7.2f %7.2f %7.4f %7.3f %7.3f %7.2f %7.2f %7.2f %7.2f %7.2f %7.3f %7.3f %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %7.2f %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %7.2f %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e \n" % (self.ediscsID[i],self.ra[i],self.dec[i],self.xcorr[i],self.ycorr[i],self.starflag[i],self.EWha[ih],self.EWhaerr[ih],self.SFR[ih],self.SFRerr[ih],self.matchflaghaediscs[i],self.SF_flag[ih], self.hstgimtype[self.matchmorphgimediscs[i]],self.matchflagmorphgimediscs[i],self.hstvisnumtype[self.matchmorphvisediscs[i]],self.matchflagmorphvisediscs[i],self.matchflagediscs24[i],self.f24[i24],self.errf24[i24],self.nmatchediscs24[i],self.misoV[i],self.misoeVapsim[i],self.misoR[i],self.misoeRapsim[i],self.misoI[i],self.misoeIapsim[i],self.misoJ[i],self.misoeJapsim[i],self.misoK[i],self.misoeKapsim[i], self.magV[i],self.mageVapsim[i],self.magR[i],self.mageRapsim[i],self.magI[i],self.mageIapsim[i],self.magJ[i],self.mageJapsim[i],self.magK[i],self.mageKapsim[i], self.membflag[i],self.newspecmatchflag[i],self.defmembership[i],self.specz[ispec],self.spectype[ispec],self.specEWOII[ispec],self.matchflagspecediscs[i],self.specEWOIIflag[ispec],self.bestz[i], self.lowz[i], self.highz[i], self.wmin[i], self.Pclust[i],self.LUlowzclust[i],self.LUbestzclust[i],self.LUhighzclust[i],self.LBlowzclust[i],self.LBbestzclust[i],self.LBhighzclust[i],self.LVlowzclust[i],self.LVbestzclust[i],self.LVhighzclust[i],self.LRlowzclust[i],self.LRbestzclust[i],self.LRhighzclust[i],self.LIlowzclust[i],self.LIbestzclust[i],self.LIhighzclust[i],self.LJlowzclust[i],self.LJbestzclust[i],self.LJhighzclust[i],self.LKlowzclust[i],self.LKbestzclust[i],self.LKhighzclust[i],self.fluxK[i],self.UBlowzclust[i],self.UBbestzclust[i],self.UBhighzclust[i],self.BVlowzclust[i],self.BVbestzclust[i],self.BVhighzclust[i],self.UVlowzclust[i],self.UVbestzclust[i],self.UVhighzclust[i],self.matchflagediscsirac[i],self.iracf1[ii],self.iracf2[ii],self.iracf3[ii],self.iracf4[ii],self.erriracf1[ii],self.erriracf2[ii],self.erriracf3[ii],self.erriracf4[ii],self.iracsexflag0[ii],self.iracsexflag1[ii],self.iracwch1[ii],self.iracwch2[ii],self.iracwch3[ii],self.iracwch4[ii],self.iracwmin[ii],self.nmatchediscsirac[i],self.l24[i24],self.errl24[i24],self.lha[ih],self.lhaerr[ih],self.snr24[i24],self.imagex24[i24],self.imagey24[i24],self.fap1[i24],self.fap2[i24],self.fap3[i24],self.fap4[i24],self.fap5[i24],self.fap6[i24],self.fap7[i24],self.fap8[i24],self.fap9[i24],self.fap10[i24],self.errfap1[i24],self.errfap2[i24],self.errfap3[i24],self.errfap4[i24],self.errfap5[i24],self.errfap6[i24],self.errfap7[i24],self.errfap8[i24],self.errfap9[i24],self.errfap10[i24])
			    
			    output.write(string)
			    #print i,self.matchflaghaediscs[i],self.SFR[ih],self.SFRerr[ih],self.lha[ih],self.lhaerr[ih]
		    except IndexError:
			    print  "Index Error in mastertable",i,len(self.ediscsID)

		    #print "%s \n" %(self.ediscsID[i])
		    #print string

		    if self.matchflagediscs24[i] > 0:
			    string5="%s %7.2f %8.2f %8.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f \n"%(self.ediscsID[i],self.matchflagediscs24[i],self.f24[self.matchediscs24[i]],self.errf24[self.matchediscs24[i]],self.magV[i],self.mageVapsim[i],self.membflag[i],self.newspecmatchflag[i],self.defmembership[i],self.wmin[i],self.spectype[self.matchspecediscs[i]],self.specEWOII[self.matchspecediscs[i]],self.matchflaghaediscs[i],self.EWha[self.matchhaediscs[i]],self.EWhaerr[self.matchhaediscs[i]],self.SFR[self.matchhaediscs[i]],self.SFRerr[self.matchhaediscs[i]])
		    else:
			    string5="%s %7.2f     0.00     0.00 %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f    0.00    0.00    0.00    0.00    0.00    0.00    0.00\n"%(self.ediscsID[i],self.matchflagediscs24[i],self.magV[i],self.mageVapsim[i],self.membflag[i],self.newspecmatchflag[i],self.defmembership[i],self.wmin[i])
		    output5.write(string5)
		    if self.matchflagediscs24[i] > 0:
			    output2.write(string)

		    if self.defmembership[i] > 0:
			    if self.matchflagmorphvisediscs[i] > 0.:
				    string2="%s %7.2f \n" % (self.ediscsID[i],self.hstvisnumtype[self.matchmorphvisediscs[i]])
				    output3.write(string2)
				    if self.matchflagediscs24[i] > 0:
					    string2="%s %7.2f %7.2f %7.2f \n" % (self.ediscsID[i],self.hstvisnumtype[self.matchmorphvisediscs[i]],self.f24[self.matchediscs24[i]],self.errf24[self.matchediscs24[i]])
					    output4.write(string2)

		    i=i+1
	    output.close()
	    output2.close()
	    output3.close()
	    output4.close()
	    output5.close()

    def count(self):
	total = 0.
	photo = 0.
	total24 = 0.
	photo24 = 0.
	den = 0.
	num = 0.
	den24 = 0.
	num24 = 0.
	N = 0.
	N24 = 0.
	Nha = 0.
	totalha = 0.
	b = 0.
	bspec = 0.
	spec24 = 0.
	photoha = 0.
	photob = 0.
	
	self.complete = 0.
	self.contam = 0.
	N = 0.
	N24 = 0.
	for i in range(len(self.membflag)):
	    if self.wmin[i] > .3:
		if self.newmatchflag[i] > 0:
		    total = total + 1.
		    if self.membflag[i] > 0:
			photo = photo + 1.
		if self.newmatchflag[i] > -1:
		    if self.membflag[i] > 0:
			den = den + 1.
			if self.newmatchflag[i] == 0:
			    num = num + 1.
		 
	complete = my.ratioerror(photo,total)
	contam = my.ratioerror(num,den)
	print complete, contam
	for i in range(len(self.membflag)):
	    if self.wmin[i] > .3:
		if self.matchflag24[i] > 0:
		    if self.newmatchflag[i] > 0:
			total24 = total24 + 1.
			if self.membflag[i] > 0:
			    photo24 = photo24 + 1.
		    if self.newmatchflag[i] > -1:
			if self.membflag[i] > 0:
			    den24 = den24 + 1.
			    if self.newmatchflag[i] == 0:
				num24 = num24 + 1.
    	complete24 = my.ratioerror(photo24,total24)
	self.complete = complete24[0]
	contam24 = my.ratioerror(num24,den24)
	self.contam = contam[0]
	print  complete24, contam24
	for i in range(len(self.membflag)):
	    if self.wmin[i] > .3:
		if self.membflag[i] > 0:
		    N = N + 1.		    
		if self.matchflag24[i] > 0:
		    if self.membflag[i] > 0:
			N24 = N24 + 1.
		    if self.newmatchflag[i] > 0:
			spec24 = spec24 + 1
		if self.matchflagha[i] > 0:
		    if self.membflag[i] > 0:
			Nha = Nha + 1.
		    if self.newmatchflag[i] > 0:
			totalha = totalha + 1.
			if self.membflag[i] > 0:
			    photoha = photoha + 1
		    if self.matchflag24[i] > 0:
			if self.membflag[i] > 0:
			    b = b + 1.
			if self.newmatchflag[i] > 0:
			    bspec = bspec + 1.
			    if self.membflag[i] > 0:
				photob = photob + 1.
			    
			
	print total, N, photo
	print spec24, N24, photo24
	print totalha, Nha, photoha
	print bspec, b, photob

    def L24LIsub(self,mark):
	l24=[]
	li=[]
	LIbestzclust=self.LIbestzclust*1.e10/h**2#convert to solar luminosities & cosmol

	for i in range(len(self.L24)):
	    if (self.matchflag24[i] >0.) & (self.photmembflag[i] > 0.):
		l24.append(self.L24[i])
		li.append(LIbestzclust[i])
	l24=N.array(l24,'d')
	li=N.array(li,'d')
	pylab.plot(li,l24,mark,label=self.prefix,markersize=6.)

    def LIRLIsub(self,mark):
	l24=[]
	li=[]
	LIbestzclust=self.LIbestzclust*1.e10/h**2#convert to solar luminosities & cosmol
	Lir=self.L24/.09#approx conversion from Dale models

	for i in range(len(self.L24)):
	    if (self.matchflag24[i] >0.) & (self.photmembflag[i] > 0.):
		    #print "got one"
		    l24.append(Lir[i])
		    li.append(LIbestzclust[i])
	l24=N.array(l24,'d')
	li=N.array(li,'d')
	print self.prefix,"len(l24) = ",len(l24)
	s='CL'+str(self.idra)+'-'+str(self.iddec)
	pylab.plot(li,l24,mark,label=s,markersize=6.)
    def LIRhistsub(self,mark):
	lir=[]
	morph=[]
	Lir=self.L24/.09#approx conversion from Dale models
	for i in range(len(self.L24)):
	    if (int(self.matchflag24[i]) & int(self.photmembflag[i]) & int(self.matchflagvistype[i]) > 0):
		lir.append(Lir[i])
		morph.append(self.vistype[i])
	return lir,morph
    def LIRhistsubspec(self,mark):
	lir=[]
	morph=[]
	Lir=self.L24/.09#approx conversion from Dale models
	for i in range(len(self.L24)):
	    if (int(self.matchflag24[i]) & int(self.matchflagvistype[i]) > 0) & (self.newmatchflag[i] > 0.) :
		lir.append(Lir[i])
		morph.append(self.vistype[i])
	return lir,morph

    def plotmorph(self):
	prefix=self.prefix
	mall=[]
	m24=[]
	for i in range(len(self.L24)):
	    if (int(self.photmembflag[i]) & int(self.matchflagvistype[i]) > 0.):
		mall.append(self.vistype[i])
		if self.matchflag24[i] > 0.:
		    m24.append(self.vistype[i])
	mall=N.array(mall,'f')
	m24=N.array(m24,'f')

	#xbins=N.array([-7,-6,-5,-2,1,2,3,4,5,6,7,8,9,10,11,66,111],'f')
	xbins=N.array([-6,-5,-2,1,2,3,4,5,6,7,8,9,10,11],'f')
	nmall=N.zeros(len(xbins),'f')
	nm24=N.zeros(len(xbins),'f')
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
	ntot=int(N.sum(nmall))
	ntot24=int(N.sum(nm24))

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
	ytick=N.arange(0.,40.,10.)
	#pylab.xticks(ind+width, ('St', 'C', 'E', 'S0', 'Sa', '', 'Sb','','Sc','','Sd','','Sm','Im','Irr','?','N/A') ,fontsize=10)
	pylab.xticks(ind+width, ('C', 'E', 'S0', 'Sa', '', 'Sb','','Sc','','Sd','','Sm','Im','Irr') ,fontsize=12)
	
	pylab.legend( (p1[0], p2[0]), ('All ('+str(ntot)+' gal)', '24um ('+str(ntot24)+' gal)'), shadow=True)


    def plotpositions(self,mark):
	prefix=self.prefix
	ra=[]
	dec=[]
	ra24=[]
	dec24=[]
	DA=my.DA(self.z,h)#kpc/arcsec
	for i in range(len(self.L24)):
	    if (int(self.newmatchflag[i]) > 0):
		dra=(self.ra[i]-self.rac)*3600.*DA
		ddec=(self.dec[i]-self.decc)*3600.*DA
		ra.append(dra)
		dec.append(ddec)

		if (self.matchflag24[i] > 0.):
		    ra24.append(dra)
		    dec24.append(ddec)
	pylab.plot(ra,dec,'ko',markersize=4)
	pylab.plot(ra24,dec24,'ro',markersize=6)
	#return ra,dec,ra24,dec24

    def colormag(self):
	    pylab.cla()
	    pylab.clf()
	    name='CL'+str(self.ra)+'-'+str(self.dec)

	    self.colormagallsub()
	    pylab.legend(loc='upper right')

	    pylab.ylabel('V-I',fontsize=26)
	    pylab.xlabel('I',fontsize=26)
	    title='V-I vs. I'+str(name)

	    #pylab.axis([19.,24.,0.,5,])
	    file=self.prefix+'colormag.eps'
	    pylab.savefig(file)

    def colormagallsub(self):
	    name='CL'+str(self.idra)+'-'+str(self.iddec)
	    print name
	    dm=2.5*N.log10(self.DL/cl1216.DL)

	    flux1=[]
	    errflux1=[]
	    flux2=[]
	    errflux2=[]	
	    flux1h=[]
	    errflux1h=[]
	    flux2h=[]
	    errflux2h=[]	
	    flux124=[]
	    errflux124=[]
	    flux224=[]
	    errflux224=[]	
	    
	    pflux1=[]
	    perrflux1=[]
	    pflux2=[]
	    perrflux2=[]	
	    pflux1h=[]
	    perrflux1h=[]
	    pflux2h=[]
	    perrflux2h=[]	
	    pflux124=[]
	    perrflux124=[]
	    pflux224=[]
	    perrflux224=[]	
	    n03=0#number w/photoz from zc-0.3 to zc+.3
	    n031=0#number w/photoz from zc-0.3 to z=1
	    n02=0#number w/photoz from zc+/-0.2
	    n03a=0#number w/photoz from zc-0.3  to a max zbest=1
	    for i in range(len(self.ra)):
	    #print name,' phot,spec,wmin= ',self.membflag[i],self.newmatchflag[i],self.wmin[i]
		    if self.matchflag24[i] > 0.:#count how many 24 micron sources are within broad photoz cuts from cluster
			    if abs(self.z-self.bestz[i]) < 0.3:
				    n03=n03+1
				    if self.bestz[i] < 1:
					    n03a=n03a+1
			    if abs(self.z-self.bestz[i]) < 0.2:
				    n02=n02+1

			    if (self.bestz[i] > 0.3) & (self.bestz[i] < 1.):
				    n031=n031+1
			    
		    if self.newmatchflag[i] > 0:
			    flux2.append(self.magV[i]-self.magI[i])
			    flux1.append(self.magI[i])
			    if self.matchflagha[i] > 0:
				    if self.SFflag[i] > 0:
					    flux2h.append(self.magV[i]-self.magI[i])
					    flux1h.append(self.magI[i])

			    if self.matchflag24[i] > 0:
				    flux224.append(self.magV[i]-self.magI[i])
				    flux124.append(self.magI[i])
			    #print "spec+24 ",self.ediscsID[i],N.log10(self.L24[i]),self.vistype[i]
			    s='circle(%12.8f,%12.8f,4.0" \n'%(self.ra[i],self.dec[i])
			    #dsfile.write(s)
		    elif (self.photmembflag[i] > 0.):
			    pflux2.append(self.magV[i]-self.magI[i])
			    pflux1.append(self.magI[i])
			    if self.matchflagha[i] > 0:
				    if self.SFflag[i] > 0:
					    pflux2h.append(self.magV[i]-self.magI[i])
					    pflux1h.append(self.magI[i])
			    
			    if self.matchflag24[i] > 0:
				    pflux224.append(self.magV[i]-self.magI[i])
				    pflux124.append(self.magI[i])

	
	    flux1=N.array(flux1,'d')
	    flux2=N.array(flux2,'d')
	    flux1h=N.array(flux1h,'d')
	    flux2h=N.array(flux2h,'d')
	    flux124=N.array(flux124,'d')
	    flux224=N.array(flux224,'d')

	    pflux1=N.array(pflux1,'d')
	    pflux2=N.array(pflux2,'d')
	    pflux1h=N.array(pflux1h,'d')
	    pflux2h=N.array(pflux2h,'d')
	    pflux124=N.array(pflux124,'d')
	    pflux224=N.array(pflux224,'d')
		
	    #pylab.plot(flux1,flux2,'bo')
	    title='V-I vs. I'+str(name)
	    #pylab.title(name)

	    pylab.plot(pflux1,pflux2, 'wo',markeredgecolor='k',markeredgewidth=1.,markersize=3,label="Photo-z Members")
	    
	    pylab.plot(pflux124,pflux224,'wo',markeredgecolor='r',mew=1.,markersize=5,label="Photo-z + 24um")
	    #pylab.plot(pflux1h,pflux2h,'wo',markeredgecolor='b',mew=1.,markersize=3,label="Photo-z + Halpha")
	    
	    pylab.plot(flux1,flux2, 'ko',markersize=6,label="Spec Members")
	    
	    pylab.plot(flux124,flux224,'ro',markersize=10,label="Spec + 24um")
	    #pylab.plot(flux1h,flux2h,'bo',markersize=6,label="Spec +Halpha")
	

	    print self.prefix, "z = ",self.z
	    print "Nspec+24, Nphot+24 = ",len(flux124),len(pflux124)
	    print "N 24 with abs(bestz-z) < 0.2 = ",n02
	    print "N 24 with abs(bestz-z) < 0.3 = ",n03
	    print "N 24 with abs(bestz-z) < 0.3, but capping at zbest < 1 = ",n03a
	    print "N 24 with bestz> 0.3 and bestz < 1 = ",n031

	#pylab.legend(loc='lower left')


	    dm=2.5*N.log10(self.DL/cl1216.DL)
	    x=N.arange(19.8,24.4,1.)
	    y=-.09*(x-20)+2.8+dm
	    pylab.plot(x,y,'k-',label="_nolegend_")
	    
	    pylab.xticks(N.arange(20,24))
	    pylab.yticks(N.arange(1,4))
	    pylab.axis([19.5,23.5,0.25,3.9])

	    s=name+', z=%5.3f, Nsp=%2d, Nph=%2d'%(self.z,len(flux124),len(pflux124))
	    pylab.text(19.8,3.2,s,fontsize=14)
#	pylab.axis([0,3.5,14,22])
	#file=str(name)+'colormag.jpg'
	#pylab.savefig(file)

    def colormagallsubspec(self,name,dm):
	flux1=[]
	errflux1=[]
	flux2=[]
	errflux2=[]	
	flux1h=[]
	errflux1h=[]
	flux2h=[]
	errflux2h=[]	
	flux124=[]
	errflux124=[]
	flux224=[]
	errflux224=[]	

	pflux1=[]
	perrflux1=[]
	pflux2=[]
	perrflux2=[]	
	pflux1h=[]
	perrflux1h=[]
	pflux2h=[]
	perrflux2h=[]	
	pflux124=[]
	perrflux124=[]
	pflux224=[]
	perrflux224=[]	

	for i in range(len(self.ra)):
	    #print name,' phot,spec,wmin= ',self.membflag[i],self.newmatchflag[i],self.wmin[i]
	    if self.newmatchflag[i] > 0:
		    flux2.append(self.magV[i]-self.magI[i])
		    flux1.append(self.magI[i])
		    if self.matchflagha[i] > 0:
			if self.SFflag[i] > 0:
			    flux2h.append(self.magV[i]-self.magI[i])
			    flux1h.append(self.magI[i])

		    if self.matchflag24[i] > 0:
			flux224.append(self.magV[i]-self.magI[i])
			flux124.append(self.magI[i])
			#print "spec+24 ",self.ediscsID[i],N.log10(self.L24[i]),self.vistype[i]
			#s='circle(%12.8f,%12.8f,4.0" \n'%(self.ra[i],self.dec[i])
			#dsfile.write(s)
	    elif (self.photmembflag[i] > 0.):
		    pflux2.append(self.magV[i]-self.magI[i])
		    pflux1.append(self.magI[i])
		    if self.matchflagha[i] > 0:
			if self.SFflag[i] > 0:
			    pflux2h.append(self.magV[i]-self.magI[i])
			    pflux1h.append(self.magI[i])
			    
		    if self.matchflag24[i] > 0:
			pflux224.append(self.magV[i]-self.magI[i])
			pflux124.append(self.magI[i])

	
	flux1=N.array(flux1,'d')
	flux2=N.array(flux2,'d')
	flux1h=N.array(flux1h,'d')
	flux2h=N.array(flux2h,'d')
	flux124=N.array(flux124,'d')
	flux224=N.array(flux224,'d')

	pflux1=N.array(pflux1,'d')
	pflux2=N.array(pflux2,'d')
	pflux1h=N.array(pflux1h,'d')
	pflux2h=N.array(pflux2h,'d')
	pflux124=N.array(pflux124,'d')
	pflux224=N.array(pflux224,'d')
		
	    #pylab.plot(flux1,flux2,'bo')
	title='V-I vs. I'+str(name)
	#pylab.title(title)

	#pylab.plot(pflux1,pflux2, 'wo',markeredgecolor='k',markeredgewidth=1.,markersize=3,label="Photo-z Members")

	#pylab.plot(pflux124,pflux224,'wo',markeredgecolor='r',mew=1.,markersize=5,label="Photo-z + 24um")
	#pylab.plot(pflux1h,pflux2h,'wo',markeredgecolor='b',mew=1.,markersize=3,label="Photo-z + Halpha")

	#pylab.plot(flux1,flux2, 'ko',markersize=16,label="Spec Members")


	pylab.plot(flux1,flux2, 'ko',markeredgecolor='k',markeredgewidth=2,markersize=8,label="Spec Members")
	pylab.plot(flux1h,flux2h,'ws',markersize=22,markeredgecolor='k',markeredgewidth=2,label="Spec +Halpha")
	pylab.plot(flux124,flux224,'k^',markersize=22,label="Spec + 24um")	



	print self.prefix
	print "Nspec+24, Nphot+24 = ",len(flux124),len(pflux124)
	#pylab.legend(loc='lower left')



	x=N.arange(20.2,23.4,1.)
	dm=5.*N.log10(self.DL/cl1216.DL)
	y=-.09*(x-20)+2.8+dm
	pylab.plot(x,y,'k-')
	
	pylab.xticks(N.arange(20,24),fontsize=20.)
	pylab.yticks(N.arange(1,4),fontsize=20.)
	pylab.axis([20.,23.5,0.5,3.8])

#	pylab.axis([0,3.5,14,22])
	#file=str(name)+'colormag.jpg'
	#pylab.savefig(file)


    def gotoitHa(self,catalogpath,prefix,ra,dec,speccat,morphidcat,mipsimage):
	    if mode < .1:
		    self.readoptIRphot(catalogpath+'photometry/cl'+str(ra)+'-'+str(dec)+'vrijk_optIRphot.v23.cat')
		    self.read24('/Users/rfinn/clusters/spitzer/MasterTables/'+str(prefix)+'mosaic_extract_final.tbl',mipsimage)
		    self.readirac(prefix)
		    self.readphotoz(catalogpath+'photoz/v241/zcat.final.cl'+str(ra)+'-'+str(dec)+'vrijk.V2.4.1.dat')
		    self.readhafile(catalogpath+'hafiles/'+str(prefix)+'sfrtableapj.dat')
		    self.readgimmorphology(catalogpath+'morphology/HST/GIM2D/'+str(morphidcat),catalogpath+'morphology/HST/GIM2D/cl'+str(ra)+'-'+str(dec)+'_drz_sci.struct2.cat')
		    self.readvismorphology(catalogpath+'morphology/HST/Visual/cl'+str(ra)+'_'+str(dec)+'.f.cat')
		    self.readspectroscopy(catalogpath+'spectroscopy/'+str(speccat))
		    self.memb()
		    self.magnitudetofluxconversion(prefix,1.)
		    self.mastertable(prefix)
	    if mode > .1:
		    self.mode1()

    def gotoit(self):
	    if mode < .1:
		    prefix=self.prefix
		    ra=self.idra
		    dec=self.iddec
		    cluster=str(ra)+'-'+str(dec)
		    s=catalogpath+'photometry/v23/cl'+cluster+'*optIRphot*'
		    files=glob.glob(s)
		    print s
		    print files
		    phot=files[0]
		    s=catalogpath+'photoz/v241/zcat.final.cl'+cluster+'*'
		    files=glob.glob(s)
		    photz=files[0]
		    self.readoptIRphot(phot)
		    #print "Length of ediscsID = ",len(self.ediscsID)

		    s='/Users/rfinn/clusters/spitzer/MasterTables/'+str(prefix)+'mosaic_extract_final.tbl'
		    self.read24(s,self.mipsimage)
		    self.readirac(prefix)
		    
		    self.readphotoz(photz)
		    s=catalogpath+'morphology/HST/GIM2D/cl'+cluster+'*'
		    files=glob.glob(s)
		    if len(files) > 0:
			    morphcat=files[0]
			    s=catalogpath+'morphology/HST/GIM2D/'+prefix+'*.id'
			    files=glob.glob(s)
			    #print s
			    #print files
			    morphid=files[0]

			    self.readgimmorphology(morphid,morphcat)
			    self.readvismorphology(catalogpath+'morphology/HST/Visual/cl'+str(ra)+'_'+str(dec)+'.f.cat')
			    
		    s=catalogpath+'spectroscopy/final/'+prefix+'*'
		    files=glob.glob(s)
		    speccat=files[0]
		    if speccat.find('v1') > -1:
			    self.readspectroscopy(speccat)
		    if speccat.find('v24') > -1:
			    self.readspectroscopyv2(speccat)
					    
		    self.memb()
		    self.magnitudetofluxconversion(prefix,0.)
		    self.mastertable(prefix)

	    if mode > 0.1:
		    self.mode1()
    def mode1(self):
	    self.readmaster()
	    #self.makeplots()

    def makeplots(self):
	#self.plothavs24(prefix)
	#self.plothavs24morph(prefix)
	self.count()
	self.colormag()
	#self.plotston(prefix)
	#self.hist(prefix)
	#self.morphhist(prefix)
	#self.plot24vsEWOII(prefix)


def plotcolormagall():
	matplotlib.rc('xtick',labelsize=18)
	matplotlib.rc('ytick',labelsize=18)

	xlabx=23.5
	xlaby=-1.
	ylabx=18.75
	ylaby=5.5
	pylab.cla()
	pylab.clf()
	pylab.subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.9,wspace=0.001,hspace=0.001)
	pylab.subplot(321)
	cl1216.colormagallsub()

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))

	pylab.subplot(322)
	cl1354.colormagallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(323)
	cl105412.colormagallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(324)
	cl1040.colormagallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(325)
	cl105411.colormagallsub()
	pylab.text(xlabx,xlaby,'I-band Magnitude',fontsize=24,horizontalalignment='center')
	pylab.text(ylabx,ylaby,'V-I',fontsize=24,verticalalignment='center',rotation='vertical')

	pylab.subplot(326)
	cl1227.colormagallsub()

	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.savefig('colormagalla.eps')

	pylab.cla()
	pylab.clf()
	pylab.subplot(321)
	cl1103.colormagallsub()

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(322)
	cl1353.colormagallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(323)
	cl1037.colormagallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))

	pylab.subplot(324)
	cl1232.colormagallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(325)
	cl1411.colormagallsub()
	pylab.text(xlabx,xlaby,'I-band Magnitude',fontsize=24,horizontalalignment='center')
	pylab.text(ylabx,ylaby,'V-I',fontsize=24,verticalalignment='center',rotation='vertical')


	pylab.subplot(326)
	cl1420.colormagallsub()
	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.savefig('colormagallb.eps')

	pylab.cla()
	pylab.clf()
	pylab.subplot(321)
	cl1301.colormagallsub()

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))

	pylab.subplot(322)
	cl1138.colormagallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	pylab.subplot(323)
	cl1018.colormagallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))

	pylab.subplot(324)
	cl1059.colormagallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	pylab.subplot(325)
	cl1202.colormagallsub()
	pylab.text(xlabx,xlaby,'I-band Magnitude',fontsize=24,horizontalalignment='center')
	pylab.text(ylabx,ylaby,'V-I',fontsize=24,verticalalignment='center',rotation='vertical')



	#pylab.subplot(326)
	#pylab.legend(loc='upper right')
	pylab.savefig('colormagallc.eps')

def plotLIRLI():
    pylab.cla()
    pylab.clf()
    cl1216.LIRLIsub('bo')
    cl1354.LIRLIsub('ro')
    cl105412.LIRLIsub('go')
    cl1040.LIRLIsub('ko')
    cl105411.LIRLIsub('yo')
    cl1227.LIRLIsub('co')
    #cl1216.LIRLIsub('bo')
    pylab.xlabel(r'$\rm{L(I)/L_\odot}$',fontsize=24)
    pylab.ylabel(r'$\rm{L(IR)/L_\odot}$',fontsize=24)
    pylab.legend(loc='upper left')
    x=N.arange(1.,2.e12,1.e11)
    y=1.e11*N.ones(len(x),'d')
    pylab.plot(x,y,'k--')
    pylab.text(2.e8,1.1e11,'LIRGS',fontsize=18)
    pylab.axis([1.e8,1.e12,9.e9,1.e12])
    ax=pylab.gca()
    ax.set_yscale("log")
    ax=pylab.gca()
    ax.set_xscale("log")

    pylab.savefig('LIRLI.eps')


def plotpositions():
	pylab.cla()
	pylab.clf()
	c=cl1037
	pylab.plot(c.imagex24,c.imagey24,'b.')
	c=cl1040
	pylab.plot(c.imagex24,c.imagey24,'g.')
	c=cl105411
	pylab.plot(c.imagex24,c.imagey24,'r.')
	c=cl105412
	pylab.plot(c.imagex24,c.imagey24,'c.')
	c=cl1216
	pylab.plot(c.imagex24,c.imagey24,'k.')
	c=cl1227
	pylab.plot(c.imagex24,c.imagey24,'bs')
	c=cl1232
	pylab.plot(c.imagex24,c.imagey24,'gs')	
	c=cl1354
	pylab.plot(c.imagex24,c.imagey24,'rs')
	c=cl1018
	pylab.plot(c.imagex24,c.imagey24,'cs')
	c=cl1059
	pylab.plot(c.imagex24,c.imagey24,'ks')
	c=cl1103
	pylab.plot(c.imagex24,c.imagey24,'b^')
	c=cl1138
	pylab.plot(c.imagex24,c.imagey24,'g^')
	c=cl1202
	pylab.plot(c.imagex24,c.imagey24,'r^')
	c=cl1301
	pylab.plot(c.imagex24,c.imagey24,'c^')
	c=cl1353
	pylab.plot(c.imagex24,c.imagey24,'k^')
	c=cl1411
	pylab.plot(c.imagex24,c.imagey24,'bo')
	c=cl1420
	pylab.plot(c.imagex24,c.imagey24,'go')
	pylab.savefig('position24um.eps')

def plotsnrR():
	pylab.cla()
	pylab.clf()
	c=cl1037
	snrsubplot(c,'b.')
	c=cl1040
	snrsubplot(c,'g.')
	c=cl105411
	snrsubplot(c,'r.')
	c=cl105412
	snrsubplot(c,'c.')
	c=cl1216
	snrsubplot(c,'k.')
	c=cl1227
	snrsubplot(c,'bs')
	c=cl1232
	snrsubplot(c,'gs')	
	c=cl1354
	snrsubplot(c,'rs')
	c=cl1018
	snrsubplot(c,'cs')
	c=cl1059
	snrsubplot(c,'ks')
	c=cl1103
	snrsubplot(c,'b^')
	c=cl1138
	snrsubplot(c,'g^')
	c=cl1202
	snrsubplot(c,'r^')
	c=cl1301
	snrsubplot(c,'c^')
	c=cl1353
	snrsubplot(c,'k^')
	c=cl1411
	snrsubplot(c,'bo')
	c=cl1420
	snrsubplot(c,'go')
	pylab.savefig('snr24R.eps')

def snrsubplot(c,symbol):
	pylab.plot(c.magR,c.snr24,symbol)
	

def plotRIvsI():
	pylab.cla()
	pylab.clf()
	c=cl1037
	RIvsIsubplot(c,'b.')
	c=cl1040
	RIvsIsubplot(c,'g.')
	c=cl105411
	RIvsIsubplot(c,'r.')
	c=cl105412
	RIvsIsubplot(c,'c.')
	c=cl1216
	RIvsIsubplot(c,'k.')
	c=cl1227
	RIvsIsubplot(c,'bs')
	c=cl1232
	RIvsIsubplot(c,'gs')	
	c=cl1354
	RIvsIsubplot(c,'rs')
	c=cl1018
	RIvsIsubplot(c,'cs')
	c=cl1059
	RIvsIsubplot(c,'ks')
	c=cl1103
	RIvsIsubplot(c,'b^')
	c=cl1138
	RIvsIsubplot(c,'g^')
	c=cl1202
	RIvsIsubplot(c,'r^')
	c=cl1301
	RIvsIsubplot(c,'c^')
	c=cl1353
	RIvsIsubplot(c,'k^')
	c=cl1411
	RIvsIsubplot(c,'bo')
	c=cl1420
	RIvsIsubplot(c,'go')
	pylab.axis([15.,23.5,-.5,2.0])
	pylab.xlabel('R mag')
	pylab.ylabel('R - I')
	pylab.savefig('RIvsI.eps')

def RIvsIsubplot(c,symbol):
	pylab.plot(c.magR,(c.magR-c.magI),symbol)
	

def plotRest():#estimate R from V and I and plot vs R mag
	pylab.cla()
	pylab.clf()
	c=cl1037
	Restsubplot(c,'b.')
	c=cl1040
	Restsubplot(c,'g.')
	c=cl105411
	Restsubplot(c,'r.')
	c=cl105412
	Restsubplot(c,'c.')
	c=cl1216
	Restsubplot(c,'k.')
	c=cl1227
	Restsubplot(c,'bs')
	c=cl1232
	Restsubplot(c,'gs')	
	c=cl1354
	Restsubplot(c,'rs')
	c=cl1018
	Restsubplot(c,'cs')
	c=cl1059
	Restsubplot(c,'ks')
	c=cl1103
	Restsubplot(c,'b^')
	c=cl1138
	Restsubplot(c,'g^')
	c=cl1202
	Restsubplot(c,'r^')
	c=cl1301
	Restsubplot(c,'c^')
	c=cl1353
	Restsubplot(c,'k^')
	c=cl1411
	Restsubplot(c,'bo')
	c=cl1420
	Restsubplot(c,'go')
	xl=N.arange(10.,30.,1.)
	yl=0.*xl
	pylab.plot(xl,yl,'k--')
	pylab.axis([15.,23.5,-.6,.6])
	pylab.xlabel('R mag')
	pylab.ylabel('0.5*(magI + magV) - R')
	pylab.savefig('Rest.eps')

def Restsubplot(c,symbol):

	#y=2.5*N.log10(.5*10.**(c.magI/2.5)+0.5*10.**(c.magV/2.5))-.25
	y=0.5*(c.magI+c.magV)
	pylab.plot(c.magR,(y-c.magR),symbol)
	print c.prefix, " number of sources = ",len(c.imagex24)
	

gocl105412=0.
gocl105411=0.
	
	
gocl1037=0.
gocl1040=0.
gocl1216=0.    
gocl1227=1.
gocl1232=0.
gocl1354=0.

gocl1018=0.
gocl1059=0.
gocl1103=0.
gocl1138=0.
gocl1202=0.
gocl1301=0.
gocl1353=0.
gocl1411=0.
gocl1420=0.


allclusters=1.
if allclusters > .1:
	gocl105412=1.
	gocl105411=1.
	
	
	gocl1037=1.
	gocl1040=1.
	gocl1216=1.    
	gocl1227=1.
	gocl1232=1.
	gocl1354=1.
	
	gocl1018=1.
	gocl1059=1.
	gocl1103=1.
	gocl1138=1.
	gocl1202=1.
	gocl1301=1.
	gocl1353=1.
	gocl1411=1.
	gocl1420=1.

if (gocl1216 > 0.1):
	z=0.794
	prefix='cl1216'
	idra='1216'
	iddec='1201'
	rvalue='r13810688'
	print "CL1216"
	mipsimage='/Users/rfinn/clusters/spitzer/cl1216/mips/24/r13810688/ch1/bcd/pbcd/Combine/mosaic.fits'
	#mipsimage='/Users/rfinn/clusters/spitzer/cl1216/hinz/mosaic.fits'#use Joannah's image
	rac='184.1879'
	decc='-12.0214'
	cl1216 = ediscs24()
	cl1216.gotoitHa(catalogpath,'cl1216','1216','1201','cl1216.v1.spe','cl1216.id',mipsimage)
	#cl1216.readoptIRphot(catalogpath+'photometry/cl1216-1201vrijk_optIRphot.v23.cat')
	##cl1216.read24('/Users/rfinn/clusters/spitzer/MasterTables/cl1216mosaic_extract.tbl',mipsimage)
	#cl1216.read24('/Users/rfinn/clusters/spitzer/cl1216/hinz/output_apex_step2/mosaic_extract.tbl',mipsimage)#use Joannah's image
	#print "reading irac"
	#cl1216.readirac('cl1216')
	##cl1216.read24('/Users/rfinn/clusters/spitzer/cl1216/mips/24/r13810688/ch1/pbcd/output/mosaic_extract.tbl')#SSC pbcd version
	#print "reading photoz catalog"
	#cl1216.readphotoz(catalogpath+'photoz/zcat.final.cl1216-1201vrijk.V2.4.1.dat')
	#print "reading H-alpha file"
	#cl1216.readhafile(catalogpath+'hafiles/cl1216sfrtableapj.dat')
        ## ##cl1216.x2=cl1216.x2+40.
        ## ##cl1216.y2=cl1216.y2+38.
        ## ##cl1216.idra=cl1216.ra+0.8/3600.*15.
        ## ##cl1216.raha=cl1216.raha+0.8/3600.*15.
	#print "reading gim2d morphology"
	#cl1216.readgimmorphology(catalogpath+'morphology/HST/GIM2D/cl1216.id',catalogpath+'morphology/HST/GIM2D/cl1216-1201_drz_sci.struct2.cat')
	#print "reading visible morphology"
	#cl1216.readvismorphology(catalogpath+'morphology/HST/Visual/cl1216_1201.f.cat')
	#print "reading spectroscopy"
	#cl1216.readspectroscopy(catalogpath+'spectroscopy/cl1216.v1.spe')
	#print "running memb()"
	#cl1216.memb()
	#print "magnitude to flux conversion"
	#cl1216.magnitudetofluxconversion('cl1216',1.)
	#print "writing master table"
	#cl1216.mastertable('cl1216')
	##cl1216.readwide24('/Users/rfinn/clusters/spitzer/desai-mosaic_extract.tbl','cl1216wide-SSC.eps')
	##cl1216.readwide24('/Users/rfinn/clusters/spitzer/desai-mosaic_extract.tbl','cl1216wide.eps')
if (gocl1040 > 0.1):
	print "CL1040"
	z=0.704
	prefix='cl1040'
	idra='1040'
	iddec='1155'
	mipsimage='/Users/rfinn/clusters/spitzer/cl1040/mips/24/r13810432/ch1/bcd/pbcd/Combine/mosaic.fits'
	rac='160.1733'
	decc='-11.9308'
	cl1040 = ediscs24()
	cl1040.gotoitHa(catalogpath,'cl1040','1040','1155','cl1040.v1.spe','cl1040.id',mipsimage)
if (gocl105412 > 0.1):
	print "CL105412"
	z=.748
	prefix='cl105412'
	idra='1054'
	iddec='1245'
	mipsimage='/Users/rfinn/clusters/spitzer/cl105412/mips/24/r13812480/ch1/bcd/pbcd/Combine/mosaic.fits'
	rac='163.6813'
	decc='-12.7639'
	cl105412 = ediscs24()
	cl105412.gotoitHa(catalogpath,'cl105412','1054','1245','cl105412.v1.spe','cl1054-12.id',mipsimage)


if (gocl105411 > 0.1):
	print "CL105411"
	z=.697
	prefix='cl105411'
	idra='1054'
	iddec='1146'
	rvalue='r13810176'
	rac='163.6008'
	decc='-11.7717'

	cl105411 = ediscs24()
	cl105411.gotoit()

if (gocl1354 > 0.1):
	print "CL1354"
	z=.762
	prefix='cl1354'
	idra='1354'
	iddec='1230'
	rvalue='r13810944'
	rac='208.5396'
	decc='-12.5164'

	cl1354 = ediscs24()
	#(self,catalogpath,prefix,ra,dec,speccat,morphidcat,mipsimage):
	cl1354.gotoit()

if (gocl1037 > 0.1):
	print "CL1037"
	z=.578
	prefix='cl1037'
	idra='1037'
	iddec='1243'
	rvalue='r17053184'
	rac='0.'
	decc='0.0'
	cl1037 = ediscs24()
	cl1037.gotoit()
if (gocl1227 > 0.1):
	print "CL1227"
	z=.636
	prefix='cl1227'
	idra='1227'
	iddec='1138'
	rvalue='r17053952'
	rac='186.9746'
	decc='-11.6389'
	cl1227 = ediscs24()
	cl1227.gotoit()

if (gocl1232 > 0.1):
	print "CL1232"
	z=.541
	prefix='cl1232'
	idra='1232'
	iddec='1250'
	rvalue='r17053696'
	rac='0.'
	decc='0.'

	cl1232 = ediscs24()
	cl1232.gotoit()



if (gocl1018 > .1):
	print "CL1018"
	z=.473
	prefix='cl1018'
	idra='1018'
	iddec='1211'
	rvalue='r17783552'
	rac='0'
	decc='0.'

	cl1018 = ediscs24()
	cl1018.gotoit()


if (gocl1059 > .1):
	print "CL1059"
	z=.456
	prefix='cl1059'
	idra='1059'
	iddec='1253'
	rvalue='r17783808'
	rac='0'
	decc='0.'
	cl1059 = ediscs24()
	cl1059.gotoit()
if (gocl1103 > .1):
	print "CL1103"
	z=.626
	prefix='cl1103'
	idra='1103'
	iddec='1245'
	rvalue='r17785344'
	rac='0'
	decc='0.'
	cl1103 = ediscs24()
	cl1103.gotoit()
if (gocl1138 > .1):
	print "CL1138"
	z=.480
	prefix='cl1138'
	idra='1138'
	iddec='1133'
	rvalue='r17785600'
	rac='0'
	decc='0.'
	cl1138 = ediscs24()
	cl1138.gotoit()
if (gocl1202 > .1):
	print "CL1202"
	z=.424
	prefix='cl1202'
	idra='1202'
	iddec='1224'
	rvalue='r17784064'
	rac='0'
	decc='0.'
	cl1202 = ediscs24()
	cl1202.gotoit()
if (gocl1301 > .1):
	print "CL1301"
	z=.483
	prefix='cl1301'
	idra='1301'
	iddec='1139'
	rvalue='r17784320'
	rac='0'
	decc='0.'
	cl1301 = ediscs24()
	cl1301.gotoit()
if (gocl1353 > .1):
	print "CL1353"
	z=.588
	prefix='cl1353'
	idra='1353'
	iddec='1137'
	rvalue='r17784576'
	rac='0'
	decc='0.'
	cl1353 = ediscs24()
	cl1353.gotoit()
if (gocl1411 > .1):
	print "CL1411"
	z=.520
	prefix='cl1411'
	idra='1411'
	iddec='1148'
	rvalue='r17784832'
	rac='0'
	decc='0.'
	cl1411 = ediscs24()
	cl1411.gotoit()
if (gocl1420 > .1):
	print "CL1420"
	z=.496
	prefix='cl1420'
	idra='1420'
	iddec='1236'
	rvalue='r17785088'
	rac='0'
	decc='0.'
	cl1420 = ediscs24()
	cl1420.gotoit()



#cl0023 = ediscs24()
#cl0023.read24('/Users/Corey/spitzer/cl0023/mips/24/r13809664/ch1/pbcd/output/mosaic_extract.tbl')
#cl0023.readha2file('/Users/Corey/spitzer/cl0023/mips/24/r13809664/ch1/pbcd/output/cl0023sfrtable2.dat')


#plotcolormagall()


#GMOS plots

plotpositions()
plotsnrR()
#plotRIvsI()
plotRest()

#plotLIRLI()


