#!/usr/bin/env python

import Numeric as N
import pylab
import os
from pyraf import iraf
#import getfieldnoise

def findnearest(x1,y1,x2,y2,delta):
    dmin = 100000000000000000000000000000000000.
    matchflag=1
    for i in range(len(x2)):
	d = N.sqrt((x1-x2[i])**2 + (y1-y2[i])**2)
	if d < dmin:
	    dmin = d
	    imatch = i
	
    if dmin > delta:
	imatch = 0
	matchflag = 0
    return imatch, matchflag



class Galaxy:
    def __init__(self):#individual galaxy properties
        print "MIPS data!  MmmmmmmIPS"

    def read24(self,file):
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
	    if (float(t[13]) < 0.):
		continue
            ngal=ngal+1
        input.close()


        self.id24 = N.zeros(ngal,'d')
        self.imagex24 = N.zeros(ngal,'d')
        self.imagey24  = N.zeros(ngal,'d')
        self.ra24 = N.zeros(ngal,'d')
        self.dec24 = N.zeros(ngal,'d')
        self.f24 = N.zeros(ngal,'d')#flux
        self.errf24 = N.zeros(ngal,'d')
        self.f1 = N.zeros(ngal,'d')
        self.f2 = N.zeros(ngal,'d')
        self.f3  = N.zeros(ngal,'d')
        self.f4 = N.zeros(ngal,'d')
        self.f5 = N.zeros(ngal,'d')

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
	    if (float(t[13]) < 0.):
		continue

	    for k in range(len(t)):
		try:
		    t[k]=float(t[k])
		except:
		    t[k]=t[k]
	    #print t
            #(self.id24[i],self.imagex24[i],self.imagey24[i],self.ra24[i],self.dec24[i],self.f24[i],self.errf24[i],self.f1[i],self.f2[i],self.f3[i],self.f4[i],self.f5[i])=(float(t[0]),float(t[1]),float(t[2]),float(t[3]),float(t[4]),float(t[5]),float(t[6]),t[23],t[24],t[25],t[26],t[27])
	    (self.id24[i],self.imagex24[i],self.imagey24[i],self.ra24[i],self.dec24[i],self.f24[i],self.errf24[i],self.f1[i],self.f2[i],self.f3[i],self.f4[i],self.f5[i])=(float(t[0]),float(t[8]),float(t[10]),float(t[3]),float(t[5]),float(t[13]),float(t[14]),float(t[23]),float(t[24]),float(t[25]),float(t[26]),float(t[27]))
	    print i,self.id24[i],self.f24[i],self.f2[i],self.f3[i]
            i=i+1
        input.close()
	#for i in range(20):
	#    print i, self.id24[i],self.f24[i],self.f2[i]

    def scaledata(self,scale):
	self.f24=scale*self.f24
	self.errf24=scale*self.errf24
	self.f1=scale*self.f1
	self.f2=scale*self.f2
	self.f3=scale*self.f3
	self.f4=scale*self.f4
	self.f5=scale*self.f5
	self.errf1=scale*self.errf1
	self.errf2=scale*self.errf2
	self.errf3=scale*self.errf3
	self.errf4=scale*self.errf4
	self.errf5=scale*self.errf5

    def matchwide(self):

	delta=1.5 #max allowed offset (in arcseconds) between matched sources
	delta=delta/3600.
	x1=self.ra24
	y1=self.dec24
	x2=wide.ra24
	y2=wide.dec24

	self.matchwide24 = N.zeros(len(x1), 'i')
	self.matchflagwide24 = N.ones(len(x1), 'i')
	for i in range(len(x1)):
	    (self.matchwide24[i],self.matchflagwide24[i]) = findnearest(x1[i],y1[i],x2,y2,delta)
	print "total number of matches w/wide field data = ",N.sum(self.matchflagwide24)
	    

    def matchGTO(self):

	delta=1.5 #max allowed offset (in arcseconds) between matched sources
	delta=delta/3600.
	x1=self.ra
	y1=self.dec
	x2=gto.ra
	y2=gto.dec
	f=[]
	fsn=[]
	fg=[]
	f1=[]
	f1g=[]
	f2=[]
	f2g=[]
	f3=[]
	f3g=[]
	f4=[]
	f4g=[]
	self.matchGTO = N.zeros(len(x1), 'i')
	self.matchflagGTO = N.zeros(len(x1), 'i')
	for i in range(len(x1)):
	    (self.matchGTO[i],self.matchflagGTO[i]) = findnearest(x1[i],y1[i],x2,y2,delta)
	    if (self.matchflagGTO[i] > 0.):
		if (self.y[i] > 10.):
		    if (self.x[i]>10.):
			if (self.x[i]<123.):
			    if (self.y[i]<135.):

				f.append(self.f24[i])
				fsn.append(abs(self.f24[i]/self.errf24[i]))
				fg.append(gto.f24[self.matchGTO[i]])
				f1.append(self.f1[i])
				f1g.append(gto.f1[self.matchGTO[i]])
				
				f2.append(self.f2[i])
				f2g.append(gto.f2[self.matchGTO[i]])
				
				f3.append(self.f3[i])
				f3g.append(gto.f3[self.matchGTO[i]])
				f4.append(self.f4[i])
				f4g.append(gto.f4[self.matchGTO[i]])
	#if matchflagwide24[i] > 0:
	    #    print "found a match!\n"
	f=N.array(f,'f')
	fg=N.array(fg,'f')
	f1=N.array(f1,'f')
	f1g=N.array(f1g,'f')
	f2=N.array(f2,'f')
	f2g=N.array(f2g,'f')
	f3=N.array(f3,'f')
	f3g=N.array(f3g,'f')
	f4=N.array(f4,'f')
	f4g=N.array(f4g,'f')
	fsn=N.array(fsn,'f')

	pylab.plot(f1,f1g,'ro',label="Ap 1")	
	pylab.plot(f2,f2g,'go',label="Ap 2")
	pylab.plot(f3,f3g,'ko',label="Ap 3")
	pylab.plot(f4,f4g,'bo',label="Ap 4")
	pylab.legend(loc="upper left")
	#pylab.plot(f1,f1g,'r^')
	x=N.arange(0,(max(f4)+200.),100.)
	y=x
	pylab.plot(x,y,'k')
	pylab.xlabel('Mopex Flux')
	pylab.ylabel('DAT flux')
	r=[]
	r2=[]
	r3=[]
	r4=[]
	for i in range(len(f1)):
	    try:
		r.append(float(f1g[i])/float(f1[i]))
		r2.append(float(f2g[i])/float(f2[i]))
		r3.append(float(f3g[i])/float(f3[i]))
		r4.append(float(f4g[i])/float(f4[i]))
	    except ZeroDivisionError:
		continue
	r=N.array(r,'f')
	r2=N.array(r2,'f')
	r3=N.array(r3,'f')
	r4=N.array(r4,'f')
	print r4
	print "average value of GTO/mopex ap1 = ",N.average(r),pylab.std(r)
	print "average value of GTO/mopex ap2 = ",N.average(r2),pylab.std(r2)
	print "average value of GTO/mopex ap3= ",N.average(r3),pylab.std(r3)
	print "average value of GTO/mopex ap4= ",N.average(r4),pylab.std(r4)
	r=N.compress(fsn > 3,r)
	print "average value of GTO/mopex (SNR > 3) = ",N.average(r),pylab.std(r)

	pylab.savefig('matched.eps')
	r=fg/f
	pylab.cla()
	pylab.clf()
	pylab.plot(f,r,'bo')
	pylab.savefig('matchedratio.eps')
	print "total number of matches w/GTO data = ",N.sum(self.matchflagGTO)

	pylab.cla()
	pylab.clf()
	#f1g=f1g*146.
	pylab.plot(f1,f1g,'bo')
	#pylab.plot(f1,f1g,'r^')
	x=N.arange(0,max(f1),100.)
	y=x
	pylab.plot(x,y,'k')
	y=1.05*x
	pylab.plot(x,y,'k--')
	y=.95*x
	pylab.plot(x,y,'k--')
	pylab.xlabel('Mopex Aperture Flux')
	pylab.ylabel('DAT Aperture flux')
	print "average value of GTO/mopex apflux = ",N.average(f1g)/N.average(f1)
	pylab.savefig('matchedap.eps')

    def plotsnr(self,file,title):
	pylab.cla()
	pylab.clf()
	x=self.f24
	y=abs(self.f24/self.errf24)

	pylab.plot(x,y,'bo')
	pylab.xlabel('F(24) (uJy)')
	pylab.ylabel('SNR')
	pylab.title(title)
	ax=pylab.gca()
	ax.set_yscale("log")
	ax.set_xscale("log")
	pylab.savefig(file)
	    
    def readsexcat(self,file):
	ngal=0
	infile = open(file,'r')
	for line in infile:
	    if line.find('#') > -1:
		continue
	    ngal=ngal + 1
	infile.close()
	self.x=N.zeros(ngal,'f')
	self.y=N.zeros(ngal,'f')
	self.f1=N.zeros(ngal,'f')
	self.f2=N.zeros(ngal,'f')
	self.f3=N.zeros(ngal,'f')
	self.f4=N.zeros(ngal,'f')
	self.f5=N.zeros(ngal,'f')
	self.f6=N.zeros(ngal,'f')
	self.f7=N.zeros(ngal,'f')
	self.f24=N.zeros(ngal,'f')
	self.errf1=N.zeros(ngal,'f')
	self.errf2=N.zeros(ngal,'f')
	self.errf3=N.zeros(ngal,'f')
	self.errf4=N.zeros(ngal,'f')
	self.errf5=N.zeros(ngal,'f')
	self.errf6=N.zeros(ngal,'f')
	self.errf7=N.zeros(ngal,'f')
	self.errf24=N.zeros(ngal,'f')
	self.flag=N.zeros(ngal,'f')
	self.ra=N.zeros(ngal,'f')
	self.dec=N.zeros(ngal,'f')
	infile = open(file,'r')
	i=0
	for line in infile:
	    if line.find('#') > -1:
		continue
	    t = line.split()
	    for j in range(len(t)):
		t[j]=float(t[j])
	    if t[17] > 0.1:
		continue
	    self.x[i]=t[10]
	    self.y[i]=t[11]
	    self.f1[i]=t[19]
	    self.f2[i]=t[20]
	    self.f3[i]=t[21]
	    self.f4[i]=t[22]
	    self.f5[i]=t[23]
	    self.f6[i]=t[24]
	    self.f7[i]=t[25]
	    self.f24[i]=t[49]
	    self.errf1[i]=t[26]
	    self.errf2[i]=t[27]
	    self.errf3[i]=t[28]
	    self.errf4[i]=t[29]
	    self.errf5[i]=t[30]
	    self.errf6[i]=t[31]
	    self.errf7[i]=t[32]
	    self.errf24[i]=t[50]
	    self.flag[i]=t[17]
	    self.ra[i]=t[47]
	    self.dec[i]=t[48]
	    #print self.flag[i]
	    i=i+1
	self.f24=self.f1
def plotsnr(file):
    pylab.cla()
    pylab.clf()

    f=wide.f24
    errf=wide.errf24
    x=f
    y=abs(f/errf)
    pylab.plot(x,y,'k.',label="Wide",markersize=1.)

    f=mopex.f24
    errf=mopex.errf24
    x=f
    y=abs(f/errf)
    pylab.plot(x,y,'bo',label="Mopex Deep")

    f=gto.f24
    errf=gto.errf24
    x=f
    y=abs(f/errf)
    pylab.plot(x,y,'ro',label="GTO Deep")


    pylab.axhline(y=1,linewidth=2,color='k',label='SNR=1')
    pylab.axvline(x=17.8,linewidth=2,color='b',label='17.8 uJy')

    pylab.xlabel('F(24) (uJy)')
    pylab.ylabel('SNR')
    pylab.title('CL1216-1201 (LCDCS 0504)', fontsize=14)
    pylab.legend(loc="upper left")
    pylab.axis([1.,1.e5,.01,1000.])
    ax=pylab.gca()
    ax.set_yscale("log")
    ax.set_xscale("log")
    pylab.savefig(file)
	    
	    

def matchcat(file,scalex,scaley,xlab,ylab):
    file2=file+'match.eps'
    i=0
    fw=[]
    fwerr=[]
    f=[]
    ferr=[]
    f1=[]
    fw1=[]
    for i in range(len(mopex.f24)):
	if (mopex.matchflagGTO[i] > 0):
	    f.append(mopex.f24[i])
	    ferr.append(mopex.errf24[i])
	    fw.append(gto.f24[mopex.matchGTO[i]])
	    fwerr.append(gto.f24[mopex.matchGTO[i]])
	    f1.append(mopex.f3[i])
	    fw1.append(gto.f3[mopex.matchGTO[i]])

    f=scalex*N.array(f,'f')
    ferr=scalex*N.array(ferr,'f')
    fw=scaley**N.array(fw,'f')
    fwerr=scaley*N.array(fwerr,'f')
    f1=scalex*N.array(f1,'f')
    fw1=scaley*N.array(fw1,'f')
    print "number of matches points = ",len(f)
    #for i in range(10):
#	print i,f[i],fw[i]
    print "average ratio of mopex/gto flux = ",N.average(f/fw)

    print "average ratio of mopex/gto flux = ",N.average(f/fw)," +/-",pylab.std(f/fw)
    print "average ratio of mopex/gto aperture flux = ",N.average(f1/fw1)," +/-",pylab.std(f1/fw1)
    pylab.clf()
    pylab.cla()

    #pylab.subplot(221)
    #pylab.plot(f,fw,'y^',label="",markersize=5)
    lab=ylab+" vs "+xlab
    pylab.plot(f,fw,'b^',label=lab,markersize=5)
    pylab.plot(f1,fw1,'r^',label=lab,markersize=5)

    

    #ylabel=xlab+"-"+ylab+" (uJy)"
    #pylab.ylabel(ylabel)
    #x=N.arange(0.,(max(f)+100.),100.)
    #y=x
    #pylab.plot(x,y,'k',label='y=x')
    #y=1/0.868*x
    #pylab.plot(x,y,'k--',label="y=0.7 x")
    #t=pylab.legend(loc='upper left')				

    #d=N.compress(f < 200.,f)
    #w=N.compress(f < 200.,fw)
    #r=(w-d)/d
    #ax=pylab.gca()
    #ax.set_yscale("log")
    #ax.set_xscale("log")

    #print "STD for f < 250 uJy = ",pylab.std(r)
    pylab.savefig(file2)

def matchcatold():
    pylab.subplot(222)
    d=N.compress(f > 500.,f)
    w=N.compress(f > 500.,fw)
    #scale=N.average(d/w)
    #scale=5.58
    y=(f-fw)/f
    
    pylab.plot(f,y,'b^',label=lab,markersize=5)
    xlabel='F(24) '+xlab+' (uJy)'
    pylab.xlabel(xlabel,fontsize=12)
    ylabel=xlab+'-'+ylab+'/'+xlab
    pylab.ylabel(ylabel,fontsize=12)

    ax=pylab.gca()
    #ax.set_yscale("log")
    ax.set_xscale("log")

    pylab.subplot(223)
    pylab.plot(f1,fw1,'r^',label="Aperture Flux",markersize=5)
    ax=pylab.gca()
    ax.set_yscale("log")
    ax.set_xscale("log")
    #ax.set_xtick(labelsize=12)
    xlabel='F(24) '+xlab+' (uJy)'
    pylab.xlabel(xlabel,fontsize=12)
    ylabel='F(24) '+ylab+' (uJy)'
    pylab.ylabel(ylabel,fontsize=12)
    t=pylab.legend(loc='upper left')				

    pylab.savefig(file2)
def runsextractor(image,auto):
    s="sex "+str(image)
    for i in range(1000):
	os.system(s)
	os.system('getxycat.pl')
	iraf.display(image,1,fill='yes')
	iraf.display(image,2,fill='yes')
	iraf.display('check.fits',3,fill='yes')
	#iraf.display('background.fits',4)
	iraf.tvmark(2,'testxy.cat',color=204,radii=2)
	#iraf.tvmark(4,'testxy.cat',color=204,radii=2)
	if auto > 0.:
		break
	try:
	    flag= raw_input("how does it look?  1=keep, 0=adjust default.sex \n")
	    flag=float(flag)
	except ValueError:
	    print "Sorry, didn't get that.  Let's try one more time."
	    flag= raw_input("how does it look?  1=keep, 0=adjust default.sex \n")
	    flag=float(flag)
	if flag > .1:
	    break
	flag= raw_input("Edit default.sex.  Hit any key when ready \n")	

def runsextractorold(im1,im2):
    #im1='cl1216mopex-sky'
    #im2='scl1216dat2'
    #im2='cl1216mopex-sky'
    #im2='mosaic-new'
    os.system('cd /Users/rfinn/clusters/spitzer/DAT-mopex/') 
    s='sex '+im1+'.fits'
    os.system(s)
    file1=im1+'-test.cat'
    #(a1,b1)=getfieldnoise.getnoise(im1,'test.cat')
    s='cp test.cat '+file1
    os.system(s)
    s='sex '+im1+'.fits,'+im2+'.fits'
    os.system(s)
    file2=im2+'-test.cat'
    s='cp test.cat '+file2
    os.system(s)
    #(a2,b2)=getfieldnoise.getnoise(im2,'test.cat')
    a1=0
    b1=0
    a2=0
    b2=0
    return a1,b1,a2,b2
	      
def plotcompflux():
    pylab.subplot(221)
    x=mopex.f2
    y=gto.f2
    labels=['4 pixel diameter', '5 pixel diameter', '6 pixel diameter', '7 pixel diameter']
    plotcompfluxsub(x,y,'mopex flux aper 2','DAT flux aper 2',labels[0])
    pyxlabel='mopex flux auto'
    ylabel='DAT flux auto'

    pylab.subplot(222)
    x=mopex.f3
    y=gto.f3
    plotcompfluxsub(x,y,'mopex flux aper 3','DAT flux aper 3',labels[1])

    pylab.subplot(223)
    x=mopex.f4
    y=gto.f4
    plotcompfluxsub(x,y,'mopex flux aper 4','DAT flux aper 4',labels[2])

    pylab.subplot(224)
    x=mopex.f5
    y=gto.f5
    plotcompfluxsub(x,y,'mopex flux aper 5','DAT flux aper 5',labels[3])
    
    pylab.savefig('compflux.jpg')
def plotcompfluxsub(x,y,xlabel,ylabel,mylabels):
    scale=N.average(y/x)
    #y=(y-x)/y
    print scale
    pylab.plot(x,y,'bo',markersize=5,label=mylabels)    
    ystat=N.compress( (x > .1),y)
    xstat=N.compress((x > .1),x)
    diff=ystat-xstat
    print "average difference = %7.3f +/- %6.3f"%(N.average(diff),pylab.std(diff))
    s="ave(DAT-SSC) = %6.3f+/-%4.3f"%(N.average(diff),pylab.std(diff))
    pylab.text(.22,7,s)
    pylab.xlabel(xlabel)
    #pylab.ylabel('(DAT - mopex)/DAT')
    pylab.ylabel(ylabel)
    pylab.legend(loc='upper left')
    x=N.arange(.01,100.,10.)
    y=N.zeros(len(x),'f')
    y=x
    pylab.plot(x,y,'k-')
    y=x+.1*x
    pylab.plot(x,y,'r--')
    y=x-.1*x
    pylab.plot(x,y,'r--')
    #pylab.axis([-5.,50.,-.5,.5])
    pylab.axis([.2,20.,.3,12.])
    ax=pylab.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
#mopex.read24('/Users/rfinn/clusters/spitzer/MasterTables/cl1216mosaic_extract.tbl','Users/rfinn/clusters/spitzer/MasterTables/cl1216mosaic_aperture.tbl')

#im1='mosaic-sky-final'
im1='cl1216_final24-scaled'
im2='mosaic-DAT-scaled'

runsex=0
if runsex > 0:
    auto=0
    file=im1+'.fits'
    runsextractor(file,auto)
    s='cp test.cat '+im1+'-test.cat'
    os.system(s)
    
    file=im2+'.fits'
    runsextractor(file,auto)
    s='cp test.cat '+im2+'-test.cat'
    os.system(s)



mopex=Galaxy()

cat=im1+'-test.cat'
mopex.readsexcat(cat)

gto=Galaxy()
cat=im2+'-test.cat'
gto.readsexcat(cat)
mopex.matchGTO()


#plotcompflux()
#mopex.matchwide()


#mopex.plotsnr('mopexSNR.eps','Mopex Reduction of Deep Data')
#wide.plotsnr('wideSNR.eps','Wide-field Data')
#gto.plotsnr('gtoSNR.eps','GTO Reduction of Deep Data')

#plotsnr('allsnr.eps')
#matchcat()
#matchcat('apflux',1,1,'ap3 flux','ap 3 flux')
#mopex.matchwide()
#mopex.matchgto()



