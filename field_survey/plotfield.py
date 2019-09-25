#!/usr/bin/env python
import pylab
import Numeric as N
#write class filter to read catalog of each filter
class field:
    def __init__(self):
        self.ra = []
        self.dec  = []
    def readallcats(self,fieldid):
	self.readjcat(fieldid)
	self.read13cat(fieldid)
	self.read84cat(fieldid)

    def readcat(self,cat):
	    infile=open(cat,'r')
	    ngal=0
	    for line in infile:
		if line.find('#') > -1.:
		    continue
		ngal=ngal+1
	    self.m1=N.zeros(ngal,'f')
	    self.m1err=N.zeros(ngal,'f')
	    self.m2=N.zeros(ngal,'f')
	    self.m2err=N.zeros(ngal,'f')
	    self.flag=N.zeros(ngal,'f')
	    self.x=N.zeros(ngal,'f')
	    self.y=N.zeros(ngal,'f')
	    self.starclass=N.zeros(ngal,'f')
	    self.fwhm=N.zeros(ngal,'f')
	    self.mbest=N.zeros(ngal,'f')
	    self.mbesterr=N.zeros(ngal,'f')
	    self.ra=N.zeros(ngal,'f')
	    self.dec=N.zeros(ngal,'f')
	    infile.close()
	    infile=open(cat,'r')
	    i=0
	    for line in infile:
		if line.find('#') > -1.:
		    continue
		t=line.split()
		for k in range(len(t)):
		    t[k]=float(t[k])
		#print i,len(t),len(self.x),t[48],t[10],self.x[i]
		(self.x[i],self.y[i],self.mbest[i],self.mbesterr[i],self.fwhm[i],self.flag[i],self.starclass[i],self.m1[i],self.m1err[i],self.m2[i],self.m2err[i],self.ra[i],self.dec[i])=(t[10],t[11],t[12],t[13],t[16],t[17],t[18],t[33],t[40],t[34],t[41],t[47],t[48])
		#(self.x[i],self.y[i],self.mbest[i],self.mbesterr[i],self.fwhm[i],self.flag[i],self.starclass[i],self.m1[i],self.m1err[i],self.m2[i],self.m2err[i],self.ra[i])=(t[10],t[11],t[12],t[13],t[16],t[17],t[18],t[33],t[40],t[34],t[41],t[47])#,self.mbest[i],self.mbesterr[i],self.fwhm[i],self.flag[i],self.starclass[i],self.m1[i],self.m1err[i],self.m2[i],self.m2err[i],self.ra[i],self.dec[i])=(t[10],t[11],t[12],t[13],t[16],t[17],t[18],t[33],t[40],t[34],t[41],t[47],t[48])
		#print self.mbest[i]
		i=i+1
	    print "got ",len(self.mbest),' objects in ',cat

    def readjcat(self,fieldid):
	    cat=str(fieldid)+'-J-final.cat'
	    self.readcat(cat)
	    self.m1j=self.m1
	    self.m1jerr=self.m1err
	    self.m2j=self.m2
	    self.m2jerr=self.m2err
	    self.flagj=self.flag
	    self.xj=self.x
	    self.yj=self.y
	    self.classj=self.starclass
	    self.fwhmj=self.fwhm
	    self.mbestj=self.mbest
	    self.mbestjerr=self.mbesterr

    def read13cat(self,fieldid):
	    cat=str(fieldid)+'-1113-final.cat'
	    self.readcat(cat)
	    self.m113=self.m1
	    self.m113err=self.m1err
	    self.m213=self.m2
	    self.m213err=self.m2err
	    self.flag13=self.flag
	    self.x13=self.x
	    self.y13=self.y
	    self.class13=self.starclass
	    self.fwhm13=self.fwhm
	    self.mbest13=self.mbest
	    self.mbest13err=self.mbesterr
    def read84cat(self,fieldid):
	    cat=str(fieldid)+'-1184-final.cat'
	    self.readcat(cat)
	    self.m184=self.m1
	    self.m184err=self.m1err
	    self.m284=self.m2
	    self.m284err=self.m2err
	    self.flag84=self.flag
	    self.x84=self.x
	    self.y84=self.y
	    self.class84=self.starclass
	    self.fwhm84=self.fwhm
	    self.mbest84=self.mbest
	    self.mbest84err=self.mbesterr
    def plotcolormag(self,fieldid):
	#x=self.mbestj
	pylab.cla()
	pylab.clf()

	ymin=-5.
	ymax=5.
	xmin=8.
	xmax=22.
	pylab.subplot(221)
	x=self.mbest84
	y=self.mbestj-self.mbest84

	pylab.plot(x,y,'bo')
	pylab.xlabel('1184')
	pylab.ylabel('J - 1184')
	xl=N.arange((xmin-1),(xmax+1),.5)
	yl=0*xl
	pylab.plot(xl,yl,'k--')
	pylab.axis([8.,22.,ymin,ymax])
	pylab.title(fieldid)
	pylab.subplot(222)
	x=self.mbest13
	y=self.mbestj-self.mbest13

	pylab.plot(x,y,'bo')
	pylab.xlabel('1113')
	pylab.ylabel('J - 1113')
	xl=N.arange((xmin-1),(xmax+1),.5)
	yl=0*xl
	pylab.plot(xl,yl,'k--')
	pylab.axis([8.,22.,ymin,ymax])
	pylab.subplot(223)
	x=self.mbest84
	y=self.mbest13-self.mbest84

	pylab.plot(x,y,'bo')
	pylab.xlabel('1184')
	pylab.ylabel('1113 - 1184')
	xl=N.arange((xmin-1),(xmax+1),.5)
	yl=0*xl
	pylab.plot(xl,yl,'k--')
	pylab.axis([8.,22.,ymin,ymax])
	#pylab.show()
	file=str(fieldid)+'colormag.eps'
	pylab.savefig(file)
	#for i in range(len(self.mbestj)):
	 #   print i,"%8.3f %8.3f %8.3f %8.3f %8.3f"%(self.mbestj[i],self.mbest13[i],self.mbest84[i],self.mbestj[i]-self.mbest84[i],self.mbestj[i]-self.mbest13[i])

    def plotstarclass(self,fieldid):
	#x=self.mbestj
	pylab.cla()
	pylab.clf()


	xmin=8.
	xmax=24.
	ymin=-0.05
	ymax=1.05
	pylab.subplot(221)
	x=self.mbest84
	y=self.class84

	pylab.plot(x,y,'bo')
	pylab.xlabel('1184 (mag-best)')
	pylab.ylabel('SE Classifier')
	pylab.axis([xmin,xmax,ymin,ymax])
	pylab.title(fieldid)
	pylab.subplot(222)
	x=self.mbest13
	y=self.class13

	pylab.plot(x,y,'bo')
	pylab.xlabel('1113 (mag-best)')
	#pylab.ylabel('SE C')
	pylab.axis([xmin,xmax,ymin,ymax])
	pylab.subplot(223)
	x=self.mbestj
	y=self.classj
	pylab.plot(x,y,'bo')
	pylab.xlabel('J (mag-best)')
	pylab.ylabel('SE Classifier')
	pylab.axis([xmin,xmax,ymin,ymax])
	pylab.show()

    def plotfwhmmag(self,fieldid):
	#x=self.mbestj
	pylab.cla()
	pylab.clf()


	xmin=8.
	xmax=24.
	ymin=0.
	ymax=15.
	pylab.subplot(221)
	x=self.mbest84
	y=self.fwhm84

	pylab.plot(x,y,'bo')
	pylab.xlabel('1184 (mag-best)')
	pylab.ylabel('FWHM')
	pylab.axis([xmin,xmax,ymin,ymax])
	pylab.title(fieldid)
	pylab.subplot(222)
	x=self.mbest13
	y=self.fwhm13

	pylab.plot(x,y,'bo')
	pylab.xlabel('1113 (mag-best)')
	#pylab.ylabel('SE C')
	pylab.axis([xmin,xmax,ymin,ymax])
	pylab.subplot(223)
	x=self.mbestj
	y=self.fwhmj
	pylab.plot(x,y,'bo')
	pylab.xlabel('J (mag-best)')
	pylab.ylabel('FWHM (pixels)')
	pylab.axis([xmin,xmax,ymin,ymax])
	pylab.show()
    def plotclassfwhm(self,fieldid):
	#x=self.mbestj
	pylab.cla()
	pylab.clf()


	xmin=-.05
	xmax=15.
	ymin=-0.05
	ymax=1.05
	pylab.subplot(221)
	x=self.fwhm84
	y=self.class84

	pylab.plot(x,y,'bo')
	pylab.xlabel('1184 FWHM (pixels)')
	pylab.ylabel('SE Classifier')
	pylab.axis([xmin,xmax,ymin,ymax])
	pylab.title(fieldid)
	pylab.subplot(222)
	x=self.fwhm13
	y=self.class13

	pylab.plot(x,y,'bo')
	pylab.xlabel('1113 FWHM (pixels)')
	#pylab.ylabel('SE C')
	pylab.axis([xmin,xmax,ymin,ymax])
	pylab.subplot(223)
	x=self.fwhmj
	y=self.classj
	pylab.plot(x,y,'bo')
	pylab.xlabel('J FWHM (pixels)')
	pylab.ylabel('SE Classifier')
	pylab.axis([xmin,xmax,ymin,ymax])
	pylab.show()
    def makeplots(self,fieldid):
	self.plotcolormag(fieldid)
	self.plotstarclass(fieldid)
	self.plotfwhmmag(fieldid)
	self.plotclassfwhm(fieldid)

fieldid='EGSs23'
e23=field()
e23.readallcats(fieldid)
e23.makeplots(fieldid)
fieldid='EGSs24'
e24=field()
e24.readallcats(fieldid)
e24.makeplots(fieldid)

