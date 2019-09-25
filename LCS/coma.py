#!/usr/bin/env python
import sqlcl
from pylab import *


def findnearest(x1,y1,x2,y2,delta):
	dmin = 100000000000000000000000000000000000.
	matchflag=1
	nmatch=0
	for i in range(len(x2)):
		d = sqrt((x1-x2[i])**2 + (y1-y2[i])**2)
		if d < delta:
			nmatch=nmatch+1
		if d < dmin:
			dmin = d
			imatch = i

	
	if dmin > delta:
		imatch = 0
		matchflag = 0
	return imatch, matchflag,nmatch

class galaxy:
    def __init__(self):
	print "initiating galaxy class"
	self.ra=194.92940   
	self.dec=27.93860 
	self.z=0.02310

    def readBai(self):
	ngal=0
	infile=open('/Users/rfinn/clusters/spitzer/coma/mips/BaiCatalog/Bei_IRsource_shape.dat','r')
	for line in infile:
	    if line.find('#') > -1:
		continue
	    ngal=ngal+1
	infile.close()


	self.Bra=zeros(ngal,'f')#quantities from Lei Bai's catalog
	self.Bdec=zeros(ngal,'f')
	self.Bf24=zeros(ngal,'f')
	self.BsexA=zeros(ngal,'f')
	self.BsexB=zeros(ngal,'f')
	self.Br=zeros(ngal,'f')#radius from center of cluster
	infile=open('/Users/rfinn/clusters/spitzer/coma/mips/BaiCatalog/Bei_IRsource_shape.dat','r')
	i=0
	for line in infile:
	    if line.find('#') > -1:
		continue
	    f=line.split()
	    for t in range(len(f)):
		f[t]=float(f[t])
	    (self.Bra[i],self.Bdec[i],self.Bf24[i],self.BsexA[i],self.BsexB[i],self.Br[i])=f
	    i += 1
	infile.close()
	self.Brmax=max(self.Br)
	print "Max radius of spec memb (arcmin) = ",self.Brmax

    def getsdss(self):
	#dA=DA(self.z[i],h100)
	#r200arcmin=self.r200[i]*1000./dA/60.
	drsearch=self.Brmax+1.#2xR200 in arcmin for sdss query
	query="select g.ra, g.dec, g.isoA_r, g.isoB_r, n.distance from galaxy g, dbo.fGetNearbyObjEq(%12.8f,%12.8f,%8.3f) n where g.objID = n.objID and  (g.PrimTarget & 0x00000040) > 0 order by distance" % (self.ra,self.dec,drsearch)#added flags to get rid of saturated objects, stars, etc
	try:
	    lines=sqlcl.query(query).readlines()
	except IOError:
	    print "IOError for cluster",self.id[i],i," trying phot query again"
	lines=sqlcl.query(query).readlines()
	print "got number+1 phot objects = ",len(lines)
	out1=open('coma-sdss-phot.dat','w')
	for line in lines:
	    out1.write(line)
	out1.close()
    
	query="select g.ra,g.dec, g.isoA_r, g.isoB_r, n.distance,  s.z from galaxy g, specobj s, dbo.fGetNearbyObjEq(%12.8f,%12.8f,%8.3f) n where g.objid = s.bestobjid and g.objID = n.objID and (g.PrimTarget & 0x00000040) > 0 order by distance" % (self.ra,self.dec,drsearch)
    
	try:
	    lines=sqlcl.query(query).readlines()
	except IOError:
	    print "IOError for cluster",self.id[i],i," trying spec query again"
	lines=sqlcl.query(query).readlines()
	out2=open('coma-sdss-spec.dat','w')
	for line in lines:
	    out2.write(line)
	out2.close()
    

	    #dr=[]
	    #print "got number + 1 of spec objects = ",len(lines)
	    #if (len(lines) > 1.):
	#	for line in lines[1:]:
	#	    f=line.split(',')
	#	    dr.append(float(f[0]))
	#    dr=array(dr,'f')

    def readsdssphot(self):
	print "Reading in sdss-phot file"
	infile=open('coma-sdss-phot.dat','r')
	ngal=0
	for line in infile:
	    ngal += 1
	infile.close()
	ngal = ngal-1
	sra=zeros(ngal,'f')
	sdec=zeros(ngal,'f')
	sisoAr=zeros(ngal,'f')
	sisoBr=zeros(ngal,'f')
	sdr=zeros(ngal,'f')

	infile=open('coma-sdss-phot.dat','r')
	i=0
	for line in infile:
	    try:
		f=line.split(',')
		for j in range(len(f)):
		    f[j]=float(f[j])
	    except ValueError:
		continue
	    (sra[i],sdec[i],sisoAr[i],sisoBr[i],sdr[i])=f
	    i += 1
	infile.close()

	delta=1.#search radius in arcsec
	delta=delta/3600.#convert to degrees
	self.sra=zeros(len(self.Bra),'f')
	self.sdec=zeros(len(self.Bra),'f')
	self.sisoAr=zeros(len(self.Bra),'f')
	self.sisoBr=zeros(len(self.Bra),'f')
	self.sdr=zeros(len(self.Bra),'f')
	self.sphotmatchflag=zeros(len(self.Bra),'f')
	self.sphotnmatch=zeros(len(self.Bra),'f')
	for i in range(len(self.Bra)):
	    (imatch, matchflag,nmatch)=findnearest(self.Bra[i],self.Bdec[i],sra,sdec,delta)
	    self.sphotmatchflag[i]=matchflag
	    self.sphotnmatch[i]=nmatch
	    #print "found a match?", matchflag
	    if matchflag > 0:
		#print "found a match with ",nmatch," matches w/in 2 arcsec"
		self.sra[i]=sra[imatch]
		self.sdec[i]=sdec[imatch]
		self.sisoAr[i]=sisoAr[imatch]
		self.sisoBr[i]=sisoBr[imatch]
		self.sdr[i]=sdr[imatch]
	outfile=open('mips-sdss-phot.dat','w')
	for i in range(len(self.sra)):
	    s='%12.8f %12.8f %8.4f %8.4f %8.4f %5.2f %5.2f \n'%(self.Bra[i],self.Bdec[i],self.sisoAr[i],self.sisoBr[i],self.sdr[i],self.sphotmatchflag[i],self.sphotnmatch[i])
	    outfile.write(s)
	outfile.close()
	rat1=[]
	rat2=[]
	dr=[]
	f24=[]
	file2='coma.reg'
	output1=open(file2,'w')
	output1.write("global color=red font=\"helvetica 10 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\nfk5\n")

	for i in range(len(self.BsexA)):
	    if self.sphotmatchflag[i] > 0:
	    #if (self.sphotmatchflag[i] > 0) & (self.Bf24[i] > 3.):
		#rat1.append(self.BsexA[i]/self.sisoAr[i])
		rat1.append((self.BsexA[i]**2+self.BsexB[i]**2)/(self.sisoAr[i]**2+self.sisoBr[i]**2))
		rat2.append(self.BsexB[i]/self.sisoBr[i])
		dr.append(self.Br[i])
		f24.append(self.Bf24[i])

		string1 = "circle(%12.8f, %12.8f, 20\") \n"%(self.Bra[i],self.Bdec[i])
		output1.write(string1)
	output1.close()
	cla()
	clf()
	plot(dr,rat1,'bo')
	#plot(dr,rat2,'go')
	savefig('ratios.eps')

	cla()
	clf()
	plot(f24,rat1,'bo')
	#plot(dr,rat2,'go')
	savefig('f24ratios.eps')
class sextractor:
    def __init__(self):
	print "reading sextractor"

    def readcat(self):
	infile=open('test.cat','r')
	ngal=0
	for line in infile:
	    if line.find('#') > -1: 
		continue
	
	    ngal=ngal+1
	infile.close()
    
	number=zeros(ngal,'f')
	fauto=zeros(ngal,'f')
	fautoerr=zeros(ngal,'f')	
	fbest=zeros(ngal,'f')
	fbesterr=zeros(ngal,'f')
	isoarea=zeros(ngal,'f')
	xim=zeros(ngal,'f')
	yim=zeros(ngal,'f')
	ra=zeros(ngal,'f')
	dec=zeros(ngal,'f')
	elong=zeros(ngal,'f')
	ellip=zeros(ngal,'f')
	fwhmim=zeros(ngal,'f')
	fwhmworld=zeros(ngal,'f')
	flags=zeros(ngal,'f')
	classstar=zeros(ngal,'f')
	infile=open('test.cat','r')
	i=0
	for line in infile:
	    if line.find('#') > -1: 
		continue
	    f=line.split()
	    for j in range(len(f)):
		f[j]=float(f[j])
	    (number[i],fauto[i],fautoerr[i],fbest[i],fbesterr[i],isoarea[i],xim[i],yim[i],ra[i],dec[i],elong[i],ellip[i],fwhmim[i],fwhmworld[i],flags[i],classstar[i])=f
	    i += 1
	infile.close()

	plot(fwhmworld,classstar,'bo')
	xlabel('FWHM (pixels)')
	ylabel('SExtractor Classifier')
	axis([0.,10.,0.,1.])
	savefig('fwhm.eps')

#mips=sextractor()
#mips.readcat()

coma=galaxy()
coma.readBai()
#coma.getsdss()
coma.readsdssphot()
