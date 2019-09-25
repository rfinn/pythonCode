#!/usr/bin/env python

#writing some scripts to aid in the geometric distortion correction 

#import Numeric as N
import numarray as N
from pyraf import iraf
import glob
import os
import pylab
import getfieldnoise
iraf.images()
iraf.images.imutil()
iraf.immatch()


RAcenter=132.82483
DECcenter=11.8000
#dimensions of pisces image
imxmax=1024
imymax=1024

def extractphotflux():
    input=open("phot.dat",'r')
    aperture=[]
    counts=[]
    area=[]
    #j=0
    ndat=0
    for line in input:
        if line.find('#') > -1: #skip lines with '#' in them
            continue
        if line.find('.fits') > -1: #skip lines with '#' in them
            ndat=ndat+1
            continue
    input.close()
    phot=N.zeros([ndat,8],'f')
    area=N.zeros([ndat,8],'f')
    input=open("phot.dat",'r')
    i=-1
    for line in input:
        if line.find('#') > -1: #skip lines with '#' in them
            continue
        if line.find('.fits') > -1: #skip lines with '#' in them
            j=0
	    i=i+1
	    nap=0
            continue
        j=j+1
        if (j > 3):
            t = line.split()
            phot[i,nap]=float(t[1])
            area[i,nap]=float(t[2])
	    nap=nap+1
    input.close()

    return phot,area


def autoalign(images):
	ra=[]
	dec=[]
	xcenter=[]
	ycenter=[]
	dx=[]
	dy=[]
	xpscale=[]
	ypscale=[]
	for im in  images:
		iraf.imgets(image=im,param='CRVAL1')#get RA of image in deg
		t=iraf.imgets.value
		ra.append(float(t))
		iraf.imgets(image=im,param='CRPIX1')#get pix of RA 
		t=iraf.imgets.value
		xcenter.append(float(t))
		iraf.imgets(image=im,param='CD1_1')#get pix of RA 
		t=iraf.imgets.value
		xpscale.append(float(t))
		iraf.imgets(image=im,param='CRVAL2')#get RA of image
		t=iraf.imgets.value
		dec.append(float(t))
		iraf.imgets(image=im,param='CRPIX2')#get RA of image
		t=iraf.imgets.value
		ycenter.append(float(t))
		iraf.imgets(image=im,param='CD2_2')#get RA of image
		t=iraf.imgets.value
		ypscale.append(float(t))

	dra=(ra-ra[0])*N.cos(N.pi/180.*dec)#correct delta ra for cos declination
	ddec=(dec-dec[0])
	dx=dra/xpscale+(xcenter-xcenter[0])
	dy=dra/ypscale+(ycenter-ycenter[0])
	xshift=dx
	yshift=-1.*dy
	outfile=open('shifts','w')
	for i in range(len(xshift)):
		s='%8.3f %8.3f \n'%(xshift[i],yshift[i])
		outfile.write(s)
	outfile.close()

def addzp(images,zp):#add zp to header	
	for im in images:
		iraf.hedit(image=im,fields='ZEROPOINT',value=zp,add='yes',verify='no',update='yes')#add zp to image header

def addfwhm(images):#add zp to header	
	for im in images:
		print "Measure fwhm on image"
		iraf.imexam(im,1)
		fwhm=raw_input('Enter FWHM in arcsec\n')
		try:
			fwhm=float(fwhm)
		except:
			print "Sorry, didn't understand ",fwhm
			fwhm=raw_input('Enter FWHM in arcsec\n')
		iraf.hedit(image=im,fields='FWHM',value=fwhm,add='yes',verify='no',update='yes')#add fwhm to image header
def normexp(images):#normalize by exposure time
	for im in images:
		outim='n'+str(im)
		iraf.imgets(image=im,param='EXPTIME')#get RA of image
		t=iraf.imgets.value
		expt=float(t)
		iraf.imarith(im,'/',expt,outim)
def wcsalign(images):
	print images
	outfile=open('imshiftin','w')
	outfile2=open('imshiftout','w')
	for im in images:
		s=str(im)+'\n'
		outfile.write(s)
		s='s'+str(im)+'\n'
		outfile2.write(s)
	outfile.close()
	outfile2.close()
	ref=str(images[1])
	#iraf.imshift('imshiftin','imshiftout',shifts_file='shifts')
	iraf.sregister('@imshiftin',ref,'@imshiftout')
def trimimages(images):
	for i in range(len(images)):
		frame=i+1
		iraf.display(images[i],frame)
	try:
		xmin=raw_input("enter min x value (X first!) \n")
		xmin=int(xmin)
	except:
		print "Sorry, I didn't get that."
		xmin=raw_input("enter min x value \n")
		xmin=int(xmin)

	try:
		xmax=raw_input("enter max x value \n")
		xmax=int(xmax)
	except:
		print "Sorry, I didn't get that."
		xmax=raw_input("enter max x value \n")
		xmax=int(xmax)

	try:
		ymin=raw_input("enter min y value \n")
		ymin=int(ymin)
	except:
		print "Sorry, I didn't get that."
		ymin=raw_input("enter min y value \n")
		ymin=int(ymin)

	try:
		ymax=raw_input("enter max y value \n")
		ymax=int(ymax)
	except:
		print "Sorry, I didn't get that."
		ymax=raw_input("enter max y value \n")
		ymax=int(ymax)

	
	for i in range(len(images)):
		dimen='['+str(int(xmin))+':'+str(int(xmax))+','+str(int(ymin))+':'+str(int(ymax))+']'
		tim='t'+str(images[i])
		im=str(images[i])+dimen
		iraf.imcopy(im,tim)
	for i in range(len(images)):
		frame=i+1
		im='t'+str(images[i])
		iraf.display(im,frame)
	print "hope you like your trimmed images!"

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


def plotsdsspos(im):
    infile=open('opencl.astrom.star21','r')
    dra=[]
    ddec=[]
    iraf.imgets(image=im,param='RA')#get RA of image
    t=iraf.imgets.value
    t1=t.split(':')
    for i in range(len(t1)):#convert to floats
	t1[i]=float(t1[i])
    RAcenter=15.*(t1[0]+t1[1]/60.+t1[2]/3600.)#convert to degrees
    iraf.imgets(image=im,param='DEC')#get RA of image
    t=iraf.imgets.value
    t1=t.split(':')
    for i in range(len(t1)):#convert to floats
	t1[i]=float(t1[i])
    DECcenter=(t1[0]+t1[1]/60.+t1[2]/3600.)#convert to degrees
    for line in infile:
	t=line.split()
	ra=(float(t[0])-RAcenter)*N.cos(N.pi/180.*float(t[1]))#correct delta ra for cos declination
	dec=(DECcenter-float(t[1]))#makes east to the left
	ra=ra*3600.#convert to arcsec
	dec=dec*3600.
	dra.append(ra)#save offsets from center position
	ddec.append(dec)#save offsets from center
    dra=N.array(dra,'f')
    ddec=N.array(ddec,'f')
    infile.close()

    xcenter=512
    ycenter=512
    plate=0.5#arcsec/pixel
    #convert offsets to pixels on pisces at 90"
    x=dra/plate
    y=ddec/plate
    #center at (512,512) in pixel coordinates
    x=x+xcenter
    y=y+ycenter
#reflect about y axis and then swap x & y to match the orientation of pisces
    x=imxmax-x#reflect about x
    yt=x#swap x & t
    xt=y
    xt=imxmax-xt
    
    dy=-26.
    dx=20.

    for i in range(1000):
	yt=yt+dy
	xt=xt+dx
	outfile=open('sdsscoords.dat','w')
	for i in range(len(x)):
	    s="%8.2f %8.2f \n"%(xt[i],yt[i])
	    outfile.write(s)
	outfile.close()

	iraf.display(im,2)
	iraf.display(im,1)
	iraf.tvmark(2,'testxy.cat',color=207,radii=7)
	iraf.tvmark(2,'sdsscoords.dat',color=204,radii=17)
	try:
	    flag= raw_input("how does it look?  1=keep, 0=adjust positions \n")
	    flag=float(flag)
	except ValueError:
	    print "Sorry, didn't get that.  Let's try one more time."
	    flag= raw_input("how does it look?  1=keep, 0=adjust positions \n")
	    flag=float(flag)
	if flag > .1:
	    break
	dx=raw_input("enter new x offset for sdss coords\n")
	dx=float(dx)
	dy=raw_input("enter new y offset for sdss coords\n")
	dy=float(dy)

def plot2massposgeo(im):#for geometric distortion correction
    #infile=open('fp_2mass.fp_psc1374.tbl','r')
    #infile=open('/Users/rfinn/field/2MASSPSCats/2massEGS.tbl','r')
    #infile=open('/Users/rfinn/field/2MASSPSCats/2massEGS1.tbl','r')
    infile=open('/Users/rfinn/field/2MASSPSCats/EGSs01.tbl','r')
    geocor=0.
    dra=[]
    ddec=[]
    j2=[]
    j2err=[]
    id2=[]
    iraf.imgets(image=im,param='RA')#get RA of image
    t=iraf.imgets.value
    t1=t.split(':')
    for i in range(len(t1)):#convert to floats
	t1[i]=float(t1[i])
    RAcenter=15.*(t1[0]+t1[1]/60.+t1[2]/3600.)#convert to degrees
    iraf.imgets(image=im,param='DEC')#get RA of image
    t=iraf.imgets.value
    t1=t.split(':')
    for i in range(len(t1)):#convert to floats
	t1[i]=float(t1[i])
    DECcenter=(t1[0]+t1[1]/60.+t1[2]/3600.)#convert to degrees
    for line in infile:
	if line.find('\\') > -1:
	    continue
	if line.find('=') > -1:
	    continue
	if line.find('|') > -1:
	    continue
	t=line.split()
	try:
		if (float(t[6]) >15.):
			continue
	except IndexError:
		continue
	ra=(RAcenter-float(t[0]))*N.cos(N.pi/180.*float(t[1]))#correct delta ra for cos declination
	dec=(float(t[1])-DECcenter)#makes east to the left
	ra=ra*3600.#convert to arcsec
	dec=dec*3600.
	try:
		j2.append(float(t[6]))
		j2err.append(float(t[7]))
	except ValueError:
		continue

	dra.append(ra)#save offsets from center position
	ddec.append(dec)#save offsets from center
	id2.append(t[5])
    dra=N.array(dra,'f')
    ddec=N.array(ddec,'f')
    infile.close()

    iraf.imgets(image=im,param='naxis1')#get RA of image
    t=float(iraf.imgets.value)
    xcenter=t/2.
    iraf.imgets(image=im,param='naxis2')#get RA of image
    t=float(iraf.imgets.value)
    ycenter=t/2.
    plate=0.5#arcsec/pixel
    #convert offsets to pixels on pisces at 90"
    x=dra/plate
    y=ddec/plate
    dx=0.
    dy=0.
    if (geocor > 0.1):
	    x=x+xcenter#center at (512,512) in pixel coordinates
	    y=y+ycenter
#reflect about y axis and then swap x & y to match the orientation of pisces
	    x=imxmax-x#reflect about x
	    yt=x#swap x & t
	    xt=y
	    xt=imxmax-xt
    
	    dy=-26.
	    dx=20.
    else:
	    yt=y+ycenter
	    xt=x+xcenter

    for i in range(1000):
	yt=yt+dy
	xt=xt+dx
	outfile=open('2masscoords.dat','w')
	for i in range(len(x)):
	    s="%8.2f %8.2f %8.2f %8.2f %s\n"%(xt[i],yt[i],j2[i],j2err[i],id2[i])
	    outfile.write(s)
	outfile.close()

	iraf.display(im,2)
	iraf.tvmark(2,'testxy.cat',color=207,radii=5)
	iraf.tvmark(2,'2masscoords.dat',color=206,radii=7.)
	try:
	    flag= raw_input("how does it look?  1=keep, 0=adjust positions \n")
	    flag=float(flag)
	except ValueError:
	    print "Sorry, didn't get that.  Let's try one more time."
	    flag= raw_input("how does it look?  1=keep, 0=adjust positions \n")
	    flag=float(flag)
	if flag > .1:
	    break
	dx=raw_input("enter new x offset for sdss coords\n")
	dx=float(dx)
	dy=raw_input("enter new y offset for sdss coords\n")
	dy=float(dy)

def plot2masspos(im,auto):
    #infile=open('fp_2mass.fp_psc1374.tbl','r')
    #infile=open('/Users/rfinn/field/2MASSPSCats/2massEGS.tbl','r')
    #infile=open('/Users/rfinn/field/2MASSPSCats/EGSs01.tbl','r')
    imname=str(im)
    if imname.find('EGS') > -1:
	    infile=open('/Users/rfinn/field/2MASSPSCats/EGSsall.tbl','r')
    if imname.find('D17') > -1:
	    infile=open('/Users/rfinn/field/2MASSPSCats/D17all.tbl','r')
    geocor=0.
    dra=[]
    ddec=[]
    j2=[]
    j2err=[]
    id2=[]
    iraf.imgets(image=im,param='CRVAL1')#get RA of image in deg
    t=iraf.imgets.value
    RAcenter=float(t)
    iraf.imgets(image=im,param='CRVAL2')#get RA of image
    t=iraf.imgets.value
    DECcenter=float(t)
    i=0
    for line in infile:
	if line.find('\\') > -1:
	    continue
	if line.find('=') > -1:
	    continue
	if line.find('|') > -1:
	    continue
	t=line.split()
	try:
		if (float(t[6]) >15.):
			continue
	except IndexError:
		continue
	ra=(RAcenter-float(t[0]))*N.cos(N.pi/180.*float(t[1]))#correct delta ra for cos declination
	dec=(float(t[1])-DECcenter)#makes east to the left
	ra=ra*3600.#convert to arcsec
	dec=dec*3600.
	try:
		j2err.append(float(t[7]))

	except ValueError:
		continue

	j2.append(float(t[6]))
	dra.append(ra)#save offsets from center position
	ddec.append(dec)#save offsets from center
	id2.append(t[5])
	id=str(t[5])
	#if id.find('14150542') > -1:
		#print t[6],len(id2),len(j2),len(j2err),len(dra)

    dra=N.array(dra,'f')
    ddec=N.array(ddec,'f')
    infile.close()

    iraf.imgets(image=im,param='CRPIX1')#get RA of image
    t=float(iraf.imgets.value)
    xcenter=t
    iraf.imgets(image=im,param='CRPIX2')#get RA of image
    t=float(iraf.imgets.value)
    ycenter=t
    plate=0.5#arcsec/pixel
    #convert offsets to pixels on pisces at 90"
    x=dra/plate
    y=ddec/plate
    dx=0.
    dy=0.
    if (geocor > 0.1):
	    x=x+xcenter#center at (512,512) in pixel coordinates
	    y=y+ycenter
#reflect about y axis and then swap x & y to match the orientation of pisces
	    x=imxmax-x#reflect about x
	    yt=x#swap x & t
	    xt=y
	    xt=imxmax-xt
    
	    dy=-26.
	    dx=20.
    else:
	    yt=y+ycenter
	    xt=x+xcenter
	    #yt=ycenter-y
	    #xt=xcenter-x

    for i in range(1000):
	yt=yt+dy
	xt=xt+dx
	outfile=open('2masscoords.dat','w')
	for j in range(len(x)):
	    s="%8.2f %8.2f %8.3f %8.3f %s\n"%(xt[j],yt[j],j2[j],j2err[j],id2[j])
	    outfile.write(s)
	    id=str(id2[j])
	    if id.find('14150542') > -1:
		    print id2[j],j2[j],j2err[j]
	outfile.close()

	if auto > 0.:
		break
	iraf.display(im,2)
	iraf.tvmark(2,'testxy.cat',color=207,radii=5)
	iraf.tvmark(2,'2masscoords.dat',color=206,radii=7.)
	try:
	    flag= raw_input("how does it look?  1=keep, 0=adjust positions \n")
	    flag=float(flag)
	except ValueError:
	    print "Sorry, didn't get that.  Let's try one more time."
	    flag= raw_input("how does it look?  1=keep, 0=adjust positions \n")
	    flag=float(flag)
	if flag > .1:
	    break
	dx=raw_input("enter new x offset for sdss coords\n")
	dx=float(dx)
	dy=raw_input("enter new y offset for sdss coords\n")
	dy=float(dy)
def getxyfromradec(im,ra,dec):
    iraf.imgets(image=im,param='CRVAL1')#get RA of image in deg
    t=iraf.imgets.value
    RAcenter=float(t)
    iraf.imgets(image=im,param='CRVAL2')#get RA of image
    t=iraf.imgets.value
    DECcenter=float(t)
    dra=(RAcenter-ra)*N.cos(N.pi/180.*dec)#correct delta ra for cos declination
    ddec=(dec-DECcenter)
    #ra=ra*3600.#convert to arcsec
    #dec=dec*3600.

    #dra=N.array(dra,'f')
    #ddec=N.array(ddec,'f')
    iraf.imgets(image=im,param='CRPIX1')#get x value corresponding to RA 
    t=float(iraf.imgets.value)
    xcenter=t
    iraf.imgets(image=im,param='CRPIX2')#get y value corresponding to dec
    t=float(iraf.imgets.value)
    ycenter=t
    iraf.imgets(image=im,param='CD1_1')#get x value corresponding to RA 
    xplate=abs(float(iraf.imgets.value))#deg/pixel
    iraf.imgets(image=im,param='CD2_2')#get x value corresponding to RA 
    yplate=abs(float(iraf.imgets.value))#deg/pixel
    #convert offsets to pixels on pisces at 90"
    x=dra/xplate+xcenter
    y=ddec/yplate+ycenter
    return x,y

def getorigxyfromradec(im,ra,dec):#get xy on original image
    iraf.imgets(image=im,param='CRVAL1')#get RA of image in deg
    t=iraf.imgets.value
    RAcenter=float(t)
    iraf.imgets(image=im,param='CRVAL2')#get RA of image
    t=iraf.imgets.value
    DECcenter=float(t)
    dra=(RAcenter-ra)*N.cos(N.pi/180.*dec)#correct delta ra for cos declination
    ddec=(dec-DECcenter)
    #ra=ra*3600.#convert to arcsec
    #dec=dec*3600.

    xcenter=600.#reference point for ra
    ycenter=600.#reference point for dec
    iraf.imgets(image=im,param='CD1_1')#get x value corresponding to RA 
    xplate=abs(float(iraf.imgets.value))#deg/pixel
    iraf.imgets(image=im,param='CD2_2')#get x value corresponding to RA 
    yplate=abs(float(iraf.imgets.value))#deg/pixel
    #convert offsets to pixels on pisces at 90"
    x=dra/xplate+xcenter
    y=ddec/yplate+ycenter
    return x,y

def checksexcoords(im):
    #infile=open('fp_2mass.fp_psc1374.tbl','r')
    #infile=open('/Users/rfinn/field/2MASSPSCats/2massEGS.tbl','r')
    #infile=open('/Users/rfinn/field/2MASSPSCats/EGSs01.tbl','r')
    iraf.imgets(image=im,param='naxis1')#get x value corresponding to RA 
    t=float(iraf.imgets.value)
    xcenter=t/2.
    iraf.imgets(image=im,param='naxis2')#get y value corresponding to dec
    t=float(iraf.imgets.value)
    ycenter=t/2.
    print "xcenter = ",xcenter," ycenter = ",ycenter
    infile=open('test.cat','r')
    geocor=0.
    ra=[]
    dec=[]
    sflag=[]
    r=[]
    i=0
    for line in infile:
	if line.find('\\') > -1:
	    continue
	if line.find('=') > -1:
	    continue
	if line.find('|') > -1:
	    continue
	if line.find('#') > -1:
	    continue
	t=line.split()
	ra.append(float(t[47]))
	dec.append(float(t[48]))
	sflag.append(float(t[17]))
	r.append(N.sqrt((float(t[10])-xcenter)**2+(float(t[11])-ycenter)**2))
    ra=N.array(ra,'d')
    dec=N.array(dec,'d')
    infile.close()

    (x,y)=getxyfromradec(im,ra,dec)
    dx=0.
    dy=0.
    rmax=510.
    for i in range(1000):
	yt=y+dy
	xt=x+dx
	outfile=open('checksexcoords.dat','w')
	for j in range(len(x)):
		if sflag[j] < 1.:
			if r[j] < rmax:
				s="%8.2f %8.2f \n"%(xt[j],yt[j])
				outfile.write(s)
	outfile.close()
	
	iraf.display(im,2)
	iraf.tvmark(2,'testxy.cat',color=207,radii=5)
	iraf.tvmark(2,'checksexcoords.dat',color=206,radii=7.)
	try:
	    flag= raw_input("how does it look?  1=keep, 0=adjust positions \n")
	    flag=float(flag)
	except ValueError:
	    print "Sorry, didn't get that.  Let's try one more time."
	    flag= raw_input("how does it look?  1=keep, 0=adjust positions \n")
	    flag=float(flag)
	if flag > .1:
	    break
	dx=raw_input("enter new x offset for sex coords\n")
	dx=float(dx)
	dy=raw_input("enter new y offset for sex coords\n")
	dy=float(dy)
def match2mass(im,auto,filter):

    xy=1.
    ra2=[]
    dec2=[]
    j2=[]
    j2err=[]
    iraf.imgets(image=im,param='CRVAL1')#get RA of image in deg
    t=iraf.imgets.value
    RAcenter=float(t)
    iraf.imgets(image=im,param='CRVAL2')#get RA of image
    t=iraf.imgets.value
    DECcenter=float(t)

    #infile=open('/Users/rfinn/field/2MASSPSCats/EGSs01.tbl','r')
    infile=open('/Users/rfinn/field/2MASSPSCats/EGSsall.tbl','r')
    for line in infile:
	if line.find('\\') > -1:
	    continue
	if line.find('=') > -1:
	    continue
	if line.find('|') > -1:
	    continue
	t=line.split()
	try:
		if (float(t[6]) >15.):
			continue
	except IndexError:
		continue
	ra=(RAcenter-float(t[0]))*N.cos(N.pi/180.*float(t[1]))#correct delta ra for cos declination
	dec=(float(t[1])-DECcenter)#makes east to the left
	ra=ra*3600.#convert to arcsec
	dec=dec*3600.
	try:
		j2.append(float(t[6]))
		j2err.append(float(t[7]))
	except ValueError:
		continue

	ra2.append(ra)#save offsets from center position
	dec2.append(dec)#save offsets from center


    if (xy > 0.):#use 2mass coords in image xy coords until I get right RA and Dec output from sextractor
	    x2=[]
	    y2=[]
	    j2=[]
	    j2err=[]
	    #x2=N.zeros(len(ra2),'f')
	    #y2=N.zeros(len(ra2),'f')
	    id2=[]
	    infile=open('2masscoords.dat','r')
	    i=0
	    for line in infile:
		    t=line.split()
		    for j in range(len(t)-1):
			    t[j]=float(t[j])
		    #(x2[i],y2[i],j2[i],j2err[i])=t[0:(len(t)-1)]
		    x2.append(t[0])
		    y2.append(t[1])
		    j2.append(t[2])
		    j2err.append(t[3])
		    id2.append(t[len(t)-1])
		    i=i+1
	    infile.close()

    ra2=N.array(ra2,'f')
    dec2=N.array(dec2,'f')
    infile.close()
    outfile=open('temp2masspos','w')
    for i in range(len(ra2)):
	s="%15.12f %15.12f"%(ra2[i],dec2[i])
	outfile.write(s)
    outfile.close()
    os.system('rm temp2massxy')
    iraf.wcsctran(input='temp2masspos',output='temp2massxy',image=im,inwcs='world',outwcs='logical')
    #iraf.display(im,2)
    #iraf.tvmark(2,'testxy.cat',color=207,radii=2)

    #iraf.tvmark(2,'temp2massxy',color=206,radii=7.)
    ra=[]
    dec=[]
    j=[]
    jerr=[]
    x=[]
    y=[]
    infile=open('test.cat','r')
    for line in infile:
	if line.find('\\') > -1:
	    continue
	if line.find('=') > -1:
	    continue
	if line.find('|') > -1:
	    continue
	if line.find('#') > -1:
	    continue
	t=line.split()
	j.append(float(t[12]))
	jerr.append(float(t[13]))
	ra.append(float(t[47]))
	dec.append(float(t[48]))#makes east to the left
	x.append(float(t[10]))
	y.append(float(t[11]))#makes east to the left


    mp=[]#pisces mag
    mperr=[]
    m2=[]#2mass mag
    m2err=[]
    delta=10./3600.#min matching radius in degrees
    rp=[]
    xim=[]
    yim=[]
    if (xy > 0.):
	    ra=N.array(x,'f')
	    dec=N.array(y,'f')
	    ra2=N.array(x2,'f')
	    dec2=N.array(y2,'f')
	    delta=5.
	    iraf.imgets(image=im,param='naxis1')#get RA of image
	    t=float(iraf.imgets.value)
	    xcenter=t/2.
	    iraf.imgets(image=im,param='naxis2')#get RA of image
	    t=float(iraf.imgets.value)
	    ycenter=t/2.
	    r=N.sqrt((ra-xcenter)**2 + (dec-ycenter)**2)
    for i in range(len(ra)):
	(match,matchflag,nmatch) = findnearest(ra[i],dec[i],ra2,dec2,delta)
	#print i, match,matchflag,nmatch,ra[i],dec[i],delta
	if matchflag > 0.:
	    print "found a match"#,i,match,ra[i],ra2[int(match)],j[i],j2[int(match)],id2[int(match)]
	    if nmatch > 1.:
		print "Warning:  Multiple matches to 2MASS Point Source Catalog, Found ",nmatch
	    if j2[int(match)] > 12.:#cut out brightest stars that are saturated in pisces
		    
		    #if r[i] < 500.: #avoid edges
		    mp.append(j[i])
		    mperr.append(jerr[i])
		    m2.append(j2[int(match)])
		    m2err.append(j2err[int(match)])
		    rp.append(r[i])
		    xim.append(x[i])
		    yim.append(y[i])
			    
    mp=N.array(mp,'f')
    mperr=N.array(mperr,'f')
    m2=N.array(m2,'f')
    m2err=N.array(m2err,'f')
    #zp=N.average(m2-mp)
    #zp=1.90-1.4

    #mp=mp+zp-m2
    #mperr=N.sqrt(mperr**2 + m2err**2)
    rperr=N.zeros(len(rp),'f')
    xim=N.array(xim,'f')
    yim=N.array(yim,'f')
    if auto < 1.:
	    pylab.cla()
	    pylab.clf()
	    pylab.subplot(221)
	    pylab.errorbar(m2,mp,mperr,m2err,'bo')
	    pylab.xlabel('2MASS J mag')
	    #s='(Pisces J mag + %6.2f'%(zp)+') - 2MASSJ'
	    s='Pisces J mag - 2MASSJ'
	    pylab.ylabel(s)
	    
	    x=N.arange(12.,16.,1.)
	    y=0*x
	    
	    pylab.plot(x,y,'k--')
            #pylab.axis([12.0,15.,-,.2])
	    pylab.subplot(222)
	    pylab.errorbar(rp,mp,mperr,rperr,'bo')
	    pylab.xlabel('Distance from image center (pixels)')
	    #s='(Pisces '+str(filter)+' mag + %6.2f'%(zp)+') - 2MASSJ'
	    s='Pisces '+str(filter)+'mag - 2MASSJ'
	    #pylab.ylabel(s)
	    x=N.arange(0.,520.,1.)
	    y=0*x
	    pylab.plot(x,y,'k--')
	    pylab.subplot(223)
	    pylab.errorbar(xim,mp,mperr,rperr,'bo')
	    pylab.xlabel('x (pixels)')
	    #s='Pisces '+str(filter)+' mag + %6.2f'%(zp)+') - 2MASSJ'
	    s='Pisces '+str(filter)+' mag - 2MASSJ'
	    pylab.ylabel(s)
	    x=N.arange(0.,520.,1.)
	    y=0*x
	    pylab.plot(x,y,'k--')
	    pylab.subplot(224)
	    pylab.errorbar(yim,mp,mperr,rperr,'bo')
	    pylab.xlabel('y (pixels)')
	    #s='(Pisces '+str(filter)+' mag + %6.2f'%(zp)+') - 2MASSJ'
	    #pylab.ylabel(s)
	    x=N.arange(0.,520.,1.)
	    y=0*x
	    pylab.plot(x,y,'k--')
            #pylab.axis([0.,512.,-1,1])
	    pylab.show()
    return mp,mperr,m2,m2err,rp,rperr,xim,yim
def runsextractor(image,auto):
    s="sex "+str(image)
    for i in range(1000):
	os.system(s)
	os.system('getxycat.pl')
	iraf.display(image,1)
	iraf.display(image,2)
	iraf.display('check.fits',3)
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
def runsextractorJ():
	images=glob.glob('n*-J.fits')
	for i in range(len(images)):
		im=images[i]
		iraf.imgets(image=im,param='ZEROPOINT')#get RA of image in deg
		t=iraf.imgets.value
		zp=float(t)
		fwhm='1.6'
		s='sex '+str(images[i])+' -MAG_ZEROPOINT '+str(zp)+' -SEEING_FWHM '+str(fwhm)
		os.system(s)
		t=im.split('.')
		im=t[0]
		t=im.split('-')
		filter=t[1]
		pre=t[0]
		fieldid=pre[3:]
		s='cp test.cat '+str(fieldid)+'-J.cat'
		#print images[i],im,pre,fieldid
		#print s
		os.system(s)

def runsextractor2im(detectim,images,detectfilter,fieldid):
    for i in range(len(images)):
	    im=images[i]
	    print "image 1 = ",detectim, "image 2 = ",im
	    iraf.imgets(image=im,param='ZEROPOINT')#get RA of image in deg
	    t=iraf.imgets.value
	    zp=float(t)
	    #iraf.imgets(image=im,param='FWHM')#get RA of image in deg
	    #t=iraf.imgets.value
	    #fwhm=float(t)
	    fwhm=1.6
	    s='sex '+str(detectim)+','+str(images[i])+' -MAG_ZEROPOINT '+str(zp)+' -SEEING_FWHM '+str(fwhm)
	    os.system(s)
	    im=str(images[i])

	    t=im.split('.')
	    im=t[0]
	    t=im.split('-')
	    filter=t[1]
	    #fieldid=t[1]
	    iraf.imgets(image=im,param='naxis1')#get RA of image
	    t=float(iraf.imgets.value)
	    xc=t/2.
	    iraf.imgets(image=im,param='naxis2')#get RA of image
	    t=float(iraf.imgets.value)
	    yc=t/2.
	    #define useable field
	    #1113

	    apertures=N.array([3.,4.,5.,6.,7.,8.,9.,10.],'f')
	    area=N.pi*apertures**2
	    linarea=N.sqrt(area)
	    iraf.imgets(image=im,param='ANOISE')#get RA of image in deg
	    t=iraf.imgets.value
	    a=float(t)
	    iraf.imgets(image=im,param='BNOISE')#get pix of RA 
	    t=iraf.imgets.value
	    b=float(t)
	    apnoise=linarea*(a+b*a*linarea)
	    apnoise=N.array(apnoise,'d')

		    

	    #(a,b)=getfieldnoise.getnoise('check.fits','test.cat')
	    infile=open('test.cat','r')
	    ofile=str(fieldid)+'-'+str(detectfilter)+'-'+str(filter)+'-test.cat'
	    s='cp test.cat '+str(ofile)
	    os.system(s)#keep entire catalog just in case
	    ofile=str(fieldid)+'-'+str(detectfilter)+'-'+str(filter)+'-final.cat'#final.cat is cut for r< 500 pix
	    outfile=open(ofile,'w')
	    s='# a = %8.4f \n'%(a)
	    outfile.write(s)
	    s='# b = %8.4f \n'%(b)
	    outfile.write(s)
	    s='# noise = x*(a+b*a*x), where x=linear size of aperture = sqrt(area) \n '
	    outfile.write(s)
	    s='#apertures = 3,4,5,6,7,8,9,10 pixels radius\n'
	    outfile.write(s)

	    iraf.imgets(image=im,param='CRVAL1')#get RA of image in deg
	    t=iraf.imgets.value
	    RAcenter=float(t)
	    iraf.imgets(image=im,param='CRVAL2')#get RA of image
	    t=iraf.imgets.value
	    DECcenter=float(t)
	    xcenter=600.#reference point for ra
	    ycenter=600.#reference point for dec
	    iraf.imgets(image=im,param='CD1_1')#get x value corresponding to RA 
	    xplate=abs(float(iraf.imgets.value))#deg/pixel
	    iraf.imgets(image=im,param='CD2_2')#get x value corresponding to RA 
	    yplate=abs(float(iraf.imgets.value))#deg/pixel

	    iraf.imgets(image=im,param='CRVAL1')#get RA of image in deg
	    t=iraf.imgets.value
	    dRAcenter=float(t)
	    iraf.imgets(image=im,param='CRVAL2')#get RA of image
	    t=iraf.imgets.value
	    dDECcenter=float(t)
	    xcenter=600.#reference point for ra
	    ycenter=600.#reference point for dec
	    iraf.imgets(image=im,param='CD1_1')#get x value corresponding to RA 
	    dxplate=abs(float(iraf.imgets.value))#deg/pixel
	    iraf.imgets(image=im,param='CD2_2')#get x value corresponding to RA 
	    dyplate=abs(float(iraf.imgets.value))#deg/pixel
	    print "rewriting catalog for image = ",im
	    for line in infile:#keep only targets w/in 500 pixels of center
		    if line.find('#') > -1:
			    outfile.write(line)
			    continue
		    t=line.split()
		    ra=float(t[51])
		    dec=float(t[52])
		    dra=(dRAcenter-ra)*N.cos(N.pi/180.*dec)#correct delta ra for cos declination
		    ddec=(dec-dDECcenter)
		    x=dra/dxplate+xcenter
		    y=ddec/dyplate+ycenter
		    detectimflag=checkposition(x,y,detectfilter)
		    if detectimflag < 1.:#only keep objects w/in useable field of detection image
			continue
		    #check to see if object w/in useable field on image 2
		    dra=(RAcenter-ra)*N.cos(N.pi/180.*dec)#correct delta ra for cos declination
		    ddec=(dec-DECcenter)
		    x=dra/xplate+xcenter
		    y=ddec/yplate+ycenter
		    im2flag=checkposition(x,y,filter)
		    ### end check image 2
		    for j in range(len(t)):
			t[j]=float(t[j])
		    (number,magiso,magerriso,magisocor,magerrisocor,magauto,magerrauto,isoarea,fluxiso,fluxerriso,xim,yim,magbest,magerrbest,elong,ellip,fwhm,flags,seclass,fluxap1,fluxap2,fluxap3,fluxap4,fluxap5,fluxap6,fluxap7,fluxap8,errap1,errap2,errap3,errap4,errap5,errap6,errap7,errap8,map1,map2,map3,map4,map5,map6,map7,map8,merrap1,merrap2,merrap3,merrap4,merrap5,merrap6,merrap7,merrap8,raim,decim,fluxbest,fluxerrbest,fluxauto,fluxerrauto)=t[:]
		    errap1=apnoise[0]
		    errap2=apnoise[1]
		    errap3=apnoise[2]
		    errap4=apnoise[3]
		    errap5=apnoise[4]
		    errap6=apnoise[5]
		    errap7=apnoise[6]
		    errap8=apnoise[7]

		    merrap1=2.5*N.log10(N.abs(fluxap1/(fluxap1+apnoise[0])))
		    merrap2=2.5*N.log10(N.abs(fluxap2/(fluxap2+apnoise[1])))
		    merrap3=2.5*N.log10(N.abs(fluxap3/(fluxap3+apnoise[2])))
		    merrap4=2.5*N.log10(N.abs(fluxap4/(fluxap4+apnoise[3])))
		    merrap5=2.5*N.log10(N.abs(fluxap5/(fluxap5+apnoise[4])))
		    merrap6=2.5*N.log10(N.abs(fluxap6/(fluxap6+apnoise[5])))
		    merrap7=2.5*N.log10(N.abs(fluxap7/(fluxap7+apnoise[6])))
		    merrap8=2.5*N.log10(N.abs(fluxap8/(fluxap8+apnoise[7])))
		    fluxerrauto=isoarea*(a+b*a*isoarea)
		    magerrauto=2.5*N.log10(N.abs(fluxauto/(fluxauto+fluxerrauto)))

		    newfields="%12.0f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %4.1f\n"%(number,magiso,magerriso,magisocor,magerrisocor,magauto,magerrauto,isoarea,fluxiso,fluxerriso,xim,yim,magbest,magerrbest,elong,ellip,fwhm,flags,seclass,fluxap1,fluxap2,fluxap3,fluxap4,fluxap5,fluxap6,fluxap7,fluxap8,errap1,errap2,errap3,errap4,errap5,errap6,errap7,errap8,map1,map2,map3,map4,map5,map6,map7,map8,merrap1,merrap2,merrap3,merrap4,merrap5,merrap6,merrap7,merrap8,raim,decim,fluxbest,fluxerrbest,fluxauto,fluxerrauto,im2flag)
		    outfile.write(newfields)
	    outfile.close()

def checkposition(x,y,filter):
	flag=1.
	if filter.find('1113') > -1:#define useable field
		xbox=485.#box(503.36558,523.38613,972.20034,968.69692,0)
		ybox=523.
		dxbox=966.
		dybox=968.
		xc=505.##circle(505,520,510)
		yc=520.
		rc=500.
	if filter.find('1184') > -1:
		xbox=497.#box(497.23459,512,963.44178,966.94521,0)
		ybox=513.
		dxbox=963.
		dybox=966.
		xc=506.##circle(506.24144,520,510)#1184
		yc=520.
		rc=500.
	if filter.find('J') > -1:
		xbox=489.#box(489.37413,514.36259,978.11085,987.56119,0)
		ybox=514.
		dxbox=978.
		dybox=987.
		xc=497.##circle(497.52182,513.75172,510)#J
		yc=514.
		rc=500.
	xmin=xbox-dxbox/2.
	xmax=xbox+dxbox/2.
	ymin=ybox-dybox/2.
	ymax=ybox+dybox/2.
	if x < xmin:
		flag=0.
		return flag
	if x > xmax:
		flag=0.
		return flag
	if y < ymin:
		flag=0.
		return flag
	if y > ymax:
		flag=0.
		return flag
	d=N.sqrt((x-xc)**2+(y-yc)**2)
	if d > rc:
		flag=0.
		return flag
	return flag

	
def runsextractorassoc(detectim,images,detectfilter,fieldid):
	im=detectim
	iraf.imgets(image=im,param='ZEROPOINT')#get RA of image in deg
	t=iraf.imgets.value
	zp=float(t)
	fwhm='1.6'#use median seeing for now
	print "running sextractor on ",detectim
	s='sex '+str(detectim)+' -MAG_ZEROPOINT '+str(zp)+' -SEEING_FWHM '+str(fwhm)
	os.system(s)
	ra=[]
	dec=[]
	infile=open('test.cat','r')#keep only sources w/in useable field
	for line in infile:#keep only targets w/in 500 pixels of center
		if line.find('#') > -1:
			continue
		t=line.split()
		flag=checkposition(float(t[10]),float(t[11]),detectfilter)#check to see if (x,y) is w/in useable field
		if flag > .1:
			ra.append(float(t[51]))
			dec.append(float(t[52]))
	infile.close()

	ra=N.array(ra,'f')
	dec=N.array(dec,'f')

	#for i in range(len(images)):
	for i in range(1):
		im=images[i]
		iraf.imgets(image=im,param='ZEROPOINT')#get RA of image in deg
		t=iraf.imgets.value
		zp=float(t)
	        #iraf.imgets(image=im,param='FWHM')#get RA of image in deg
	        #t=iraf.imgets.value
	        #fwhm=float(t)
		fwhm='1.6'#use median seeing for now
	        #convert ra, dec list to x,y positions using image header
			      
		(x,y)=getxyfromradec(im,ra,dec)
		outfile=open('assoc.list','w')
		for j in range(len(x)):
			s="%8.2f %8.2f %i\n"%(x[j],y[j],j)
			outfile.write(s)
		outfile.close()
		iraf.display(im,1)
		iraf.display(im,2)
		iraf.tvmark(2,'assoc.list',color=204,radii=3)
		print "running sextractor on ",im
		s='sex '+str(images[i])+' -c default.sex.assoc -PARAMETERS_NAME default.param.assoc -MAG_ZEROPOINT '+str(zp)+' -SEEING_FWHM '+str(fwhm)+' -ASSOC_NAME assoc.list -ASSOC_RADIUS 3.0 -ASSOC_TYPE NEAREST -ASSOC_PARAMS 1,2 -ASSOCSELEC_TYPE MATCHED -ASSOC_DATA 1,2,3' 
		os.system(s)
		os.system('getxycat.pl')
		iraf.tvmark(2,'testxy.cat',color=207,radii=10)
		im=str(images[i])

		t=im.split('.')
		im=t[0]
		t=im.split('-')
		filter=t[1]

		iraf.imgets(image=im,param='ANOISE')#get RA of image in deg
		t=iraf.imgets.value
		a=float(t)
		iraf.imgets(image=im,param='BNOISE')#get pix of RA 
		t=iraf.imgets.value
		b=float(t)

		#(a,b)=getfieldnoise.getnoise('check.fits','test.cat')
		infile=open('test.cat','r')
		ofile=str(fieldid)+'-'+str(detectfilter)+'-'+str(filter)+'-test.cat'
		s='cp test.cat '+str(ofile)
		os.system(s)#keep entire catalog just in case
		ofile=str(fieldid)+'-'+str(detectfilter)+'-'+str(filter)+'-final.cat'#final.cat is cut for r< 500 pix
		outfile=open(ofile,'w')
		s='# a = %8.4f \n'%(a)
		outfile.write(s)
		s='# b = %8.4f \n'%(b)
		outfile.write(s)
		s='# noise = x*(a+b*a*x), where x=linear size of aperture = sqrt(area) \n '
		apertures=N.array([3.,4.,5.,6.,7.,8.,9.,10.],'f')
		linarea=N.sqrt(N.pi)*apertures

		apnoise=linarea*(a+b*a*linarea)
		apnoise=N.array(apnoise,'d')
		outfile.write(s)
		for line in infile:#keep only targets w/in 500 pixels of center
			if line.find('#') > -1:
				outfile.write(line)
				continue
			t=line.split()
			flag=checkposition(float(t[10]),float(t[11]),filter)
			oline=line[:(len(line)-1)]
			newfields=" %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %4.1f\n"%(apnoise[0],apnoise[1],apnoise[2],apnoise[3],apnoise[4],apnoise[5],apnoise[6],apnoise[7],flag)
			#print apnoise
			#print "%8.2f \n"%(apnoise[0])
			#newfields=" %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e \n"%(apnoise[:])#,flag)
			oline=oline+newfields
			outfile.write(oline)
		outfile.close()

def runsextractorphot(detectim,images,detectfilter,fieldid):
	im=detectim
	iraf.imgets(image=im,param='ZEROPOINT')#get RA of image in deg
	t=iraf.imgets.value
	zp=float(t)
	fwhm='1.6'#use median seeing for now
	print "running sextractor on ",detectim
	s='sex '+str(detectim)+' -MAG_ZEROPOINT '+str(zp)+' -SEEING_FWHM '+str(fwhm)
	os.system(s)
	ra=[]
	dec=[]
	infile=open('test.cat','r')#keep only sources w/in useable field
	f=fieldid+'-'+detectfilter+'-final.cat'
	outfile=open(f,'w')
	for line in infile:#keep only targets w/in 500 pixels of center
		if line.find('#') > -1:
			continue
		t=line.split()
		flag=checkposition(float(t[10]),float(t[11]),detectfilter)#check to see if (x,y) is w/in useable field
		if flag > .1:
			outfile.write(line)
			ra.append(float(t[51]))
			dec.append(float(t[52]))
	infile.close()

	ra=N.array(ra,'f')
	dec=N.array(dec,'f')
	npoints=len(ra)


	sky = open("sky",'w')
	for i in range(npoints):
		sky.write("0.0 \n")
	sky.close()
	aps = open("apertures",'w')
	#aps.write("3,4,5,6,7,8,9,10")
	aps.write("1.5,2.,2.5,3.,3.5,4.,4.5,10.")
	aps.close()

	iraf.digiphot()
	iraf.daophot()



	for i in range(len(images)):
	#for i in range(1):
		im=images[i]
		iraf.imgets(image=im,param='ZEROPOINT')#get RA of image in deg
		t=iraf.imgets.value
		zp=float(t)
	        #iraf.imgets(image=im,param='FWHM')#get RA of image in deg
	        #t=iraf.imgets.value
	        #fwhm=float(t)
		fwhm='1.6'#use median seeing for now
	        #convert ra, dec list to x,y positions using image header
			      
		(x,y)=getxyfromradec(im,ra,dec)
		outfile=open('photpos.in','w')
		for j in range(len(x)):
			s="%8.2f %8.2f %i\n"%(x[j],y[j],j)
			outfile.write(s)
		outfile.close()
		iraf.display(im,1)
		iraf.display(im,2)
		iraf.tvmark(2,'photpos.in',color=204,radii=3)
		print "running sextractor on ",im,' to create background-subtracted image'
		s='sex '+str(images[i])+' -MAG_ZEROPOINT '+str(zp)+' -SEEING_FWHM '+str(fwhm)
		os.system(s)#create background-subtracted image, check.fits
		#photfile=str(fieldid)+'-'+str(detectfilter)+'-'+str(filter)+'-phot.dat'
		os.system('rm phot.dat')
		print 'running phot on ',im
		iraf.phot('check.fits',coords='photpos.in',output='phot.dat',skyfile='sky',salgori='file',aperture='apertures',interactive='no',calgorith='none',verify='no',verbose='no')
		print 'running extraphotflux'
		(phot,area)=extractphotflux()

		iraf.imgets(image=im,param='ANOISE')#get RA of image in deg
		t=iraf.imgets.value
		a=float(t)
		iraf.imgets(image=im,param='BNOISE')#get pix of RA 
		t=iraf.imgets.value
		b=float(t)

		linarea=N.sqrt(area)

		apnoise=linarea*(a+b*a*linarea)
		apnoise=N.array(apnoise,'d')

		t=im.split('.')
		im=t[0]
		t=im.split('-')
		filter=t[1]


		ofile=str(fieldid)+'-'+str(detectfilter)+'-'+str(filter)+'-phot.dat'
		print 'writing output file ',ofile
		outfile=open(ofile,'w')
		s='#flux-ap1 errflux-ap1 flux-ap2 err flux-ap3 err flux-ap4 err flux-ap5 err flux-ap6 err flux-ap7 err flux-ap8 err flag\n'
		outfile.write(s)
		s='#apertures = 3,4,5,6,7,8,9,20 pixels DIAMETER\n'
		outfile.write(s)
		s='#catalog is already cut to include only objects w/in useable area of detection image\n'
		outfile.write(s)
		s='#flag = 0 if object is outside useable area of 2nd/phot image\n'
		outfile.write(s)
		s='# noise = x*(a+b*a*x), where x=linear size of aperture = sqrt(area) \n '
		outfile.write(s)
		s='# a = %8.4f \n'%(a)
		outfile.write(s)
		s='# b = %8.4f \n'%(b)
		outfile.write(s)
		for j in range(len(phot)):
			flag=checkposition(x[j],y[j],str(filter))
			newfields="%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %4.1f\n"%(phot[j,0],apnoise[j,0],phot[j,1],apnoise[j,1],phot[j,2],apnoise[j,2],phot[j,3],apnoise[j,3],phot[j,4],apnoise[j,4],phot[j,5],apnoise[j,5],phot[j,6],apnoise[j,6],phot[j,7],apnoise[j,7],flag)
			outfile.write(newfields)
		outfile.close()

		    
def runall(prefixes,cat):
    for prefix in prefixes:
	images=glob.glob(prefix)
	iraf.imgets(images[0],'FILTER')
	filter=str(iraf.imgets.value)

	for im in images:
	    runsextractor(im)

	    t=im.split('.')
	    if cat < 0.1:
		out='sdssmatch-'+str(t[0])
		refcat='sdsscoords.dat'
		s='cat '+out+' >> sdssmatch'+filter
		plotsdsspos(im)
	    if cat > .1:
		out='2massmatch-'+str(t[0])
		refcat='2masscoords.dat'
		s='cat '+out+' >> 2massmatch'+filter
		plot2masspos(im)
	    #iraf.xyxymatch(input='testxy.cat',reference='sdsscoords.dat',output=out,tolerance=20.,refpoints='refpoints',interactive='yes')
	    iraf.xyxymatch(input='testxy.cat',reference=refcat,output=out,tolerance=10.,refpoints='refpoints',interactive='no')
	    os.system(s)

	infile='sdssmatch'+filter
	infile='2massmatch'+filter
	iraf.geomap(input=infile,database='PiscesBok',transform=filter,xmin=1.,xmax=1024,ymin=1.,ymax=1024.)

def testone():
    im='gmqopen0004.fits'
    runsextractor(im)
    plot2masspos(im)

def dophot(matchstring,filename,auto,filter):
	#auto=0.
	mpall=[]
	mpallerr=[]
	m2all=[]
	m2allerr=[]
	rpall=[]
        xim=[]
        yim=[]
	rpallerr=[]
	files=glob.glob(matchstring)
	for im in files:
		runsextractor(im,auto)		
		plot2masspos(im,auto)
		(mp,mperr,m2,m2err,rp,rperr,x,y)=match2mass(im,auto,filter)
		for k in range(len(mp)):
			mpall.append(mp[k])
			mpallerr.append(mp[k])
			m2all.append(m2[k])
			m2allerr.append(m2err[k])
			rpall.append(rp[k])
			rpallerr.append(rperr[k])
			xim.append(x[k])
			yim.append(y[k])


	mpall=N.array(mpall,'f')
	mpallerr=N.array(mpallerr,'f')
	m2all=N.array(m2all,'f')
	m2allerr=N.array(m2allerr,'f')
	rpall=N.array(rpall,'f')
	rpallerr=N.array(rpallerr,'f')

	file=str(filename)+'.dat'
	outfile=open(file,'w')
	for i in range(len(mpall)):
		s='%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n'%(mpall[i],mpallerr[i],m2all[i],m2allerr[i],rpall[i],rpallerr[i],xim[i],yim[i])
		outfile.write(s)
	outfile.close()


def plotallphot(filename,filter):
	pylab.cla()
	pylab.clf()
	mpall=[]
	mpallerr=[]
	m2all=[]
	m2allerr=[]
	rpall=[]
	rpallerr=[]

	file=str(filename)+'.dat'
	infile=open(file,'r')
	ngal=0
	for line in infile:
		ngal=ngal+1
	infile.close()
	mpall=N.zeros(ngal,'f')
	mpallerr=N.zeros(ngal,'f')
	m2all=N.zeros(ngal,'f')
	m2allerr=N.zeros(ngal,'f')
	rpall=N.zeros(ngal,'f')
	rpallerr=N.zeros(ngal,'f')
	xim=N.zeros(ngal,'f')
	yim=N.zeros(ngal,'f')

	infile=open(file,'r')
	i=0
	for line in infile:
		t=line.split()
		for j in range(len(t)):
			t[j]=float(t[j])
		(mpall[i],mpallerr[i],m2all[i],m2allerr[i],rpall[i],rpallerr[i],xim[i],yim[i])=t
		i=i+1
		#mpall.append(t[0])
		#mpallerr.append(t[1])
		#m2all.append(t[2])
		#m2allerr.append(t[3])
		#rpall.append(t[4])
		#rpallerr.append(t[5])
	infile.close()
	zp=N.average(m2all-mpall)
	zperr=pylab.std(m2all-mpall)
	print "estimate of zp = %8.3f +/- %8.3f"%(zp,zperr)
	#zp=0.

	for l in range(1000):
		plotallphotsub(mpall,mpallerr,zp,m2all,m2allerr,rpall,rpallerr,xim,yim,filter)
		pylab.show()
		zpflag=raw_input('Adjust ZP (0=no,1=yes)\n')
		zpflag=float(zpflag)
		if zpflag > 0.:
			print "current zp = ",zp
			zp=raw_input('Enter new zp \n')
			zp=float(zp)
		if zpflag < 1.:
			plotallphotsub(mpall,mpallerr,zp,m2all,m2allerr,rpall,rpallerr,xim,yim,filter)
			file=str(filename)+'.eps'

			pylab.savefig(file)
			break
	return zp,zperr
def plotallphotsub(mpall,mpallerr,zp,m2all,m2allerr,rpall,rpallerr,xim,yim,filter):
	pylab.cla()
	pylab.clf()
	pylab.subplot(221)
	mp=mpall+zp
	dmag=mp-m2all
	dmagerr=N.sqrt(mpallerr**2 + m2allerr**2)
	dy=.3
	m2all=N.array(m2all,'d')
	m2allerr=N.array(m2allerr,'d')
	try:
	    dmag=N.array(dmag,'d')
	except:
	    print "Error casting dmag as array on L1324"
	    print dmag
	try:
	    dmagerr=N.array(dmagerr,'d')
	except:
	    print "Error casting dmag as array on L1324"
	    print dmagerr
	try:
	    pylab.errorbar(m2all,dmag,dmagerr,m2allerr,'bo')
	except:
	    print "error printing at L1326"
	    print m2all
	    print dmag
	    print dmagerr
	    print m2allerr
	pylab.xlabel('2MASS J mag')
	s='(Pisces '+str(filter)+' mag + %6.2f'%(zp)+') - 2MASSJ'
	pylab.ylabel(s)
	x=N.arange(12.,16.,1.)
	y=N.zeros(len(x),'f')
	print x
	print y
	pylab.plot(x,y,'k--')
	pylab.axis([12.0,15.,-1*dy,dy])
	pylab.subplot(222)
	pylab.errorbar(rpall,dmag,dmagerr,rpallerr,'bo')
	#pylab.plot(rpall,dmag,'bo')
	pylab.xlabel('Distance from image center (pixels)')
	x=N.arange(0.,520.,1.)
	y=0*x
	pylab.plot(x,y,'k--')
	pylab.axis([0.,512.,-1.*dy,dy])
	
	pylab.subplot(223)
	pylab.errorbar(xim,dmag,dmagerr,rpallerr,'bo')
	pylab.xlabel('x (pixels)')
	s='(Pisces '+str(filter)+' mag + %6.2f'%(zp)+') - 2MASSJ'
	pylab.ylabel(s)
	x=N.arange(0.,520.,1.)
	y=0*x
	pylab.plot(x,y,'k--')
	pylab.axis([0.,1024.,-1.*dy,dy])
	pylab.subplot(224)
	pylab.errorbar(yim,dmag,dmagerr,rpallerr,'bo')
	pylab.xlabel('y (pixels)')
	s='(Pisces '+str(filter)+' mag + %6.2f'%(zp)+') - 2MASSJ'
	x=N.arange(0.,520.,1.)
	y=0*x 
	pylab.plot(x,y,'k--')
	pylab.axis([0.,1024.,-1.*dy,dy])


#cat=1#0=sdss, 1=2mass
#prefixes=['mqopen3*']#,'mqopen1*','mqopen2*','mqopen3*','mqopen4*','mqopen6*']
#runall(prefixes,cat)

#testone()
#runsextractor('mopexmosaic.fits')
im='cfdq-EGSs01-J.fits'

def calibratephot():
	images=glob.glob('c*.fits')
	normexp(images)
	auto=1.
	os.system('cp /Users/rfinn/field/data/sexfiles/default.sex.2mass default.sex')
	os.system('cp /Users/rfinn/field/data/sexfiles/default.param .')
	os.system('cp /Users/rfinn/field/data/sexfiles/gauss_3.0_5x5.conv .')
	os.system('cp /Users/rfinn/field/data/sexfiles/default.nnw .')
	#namematch=['n*EGS*J.fits','n*EGS*1184.fits','n*EGS*1113.fits']
	#filename=['EGS-J-photall','EGS-1184-photall','EGS-1113-photall']
	#namematch=['n*J*.fits','n*1184*.fits','n*1113*.fits']
	#filename=['photall-J','photall-1184','photall-1113']
	#filter=['J','1184','1113']
	namematch=['n*1184*.fits','n*1113*.fits']
	filename=['photall-1184','photall-1113']
	filter=['1184','1113']
	#namematch=['n*1184*.fits','n*1113*.fits','n*J*.fits']
	#filename=['photall-1184','photall-1113','photall-J']
	#filter=['1184','1113','J']
	#namematch=['n*J*.fits']
	#filename=['photall-J']
	#filter=['J']
	#i=2
	#calibratephotsub(namematch[i],filename[i],auto,filter[i])
	for i in range(len(namematch)):
		calibratephotsub(namematch[i],filename[i],auto,filter[i])
def calibratephotsub(namematch,filename,auto,filter):
        dophot(namematch,filename,auto,filter)
	(zp,zperr)=plotallphot(filename,filter)
	images=glob.glob(namematch)
	addzp(images,zp)

def finalsteps():
	#fieldids=['EGSs23','EGSs2']
	fieldids=['EGSs24']
	cdir=str(os.getcwd())
	if cdir.find('EGSs') > -1:
	    ids=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]
	    #ids=[11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]
	    #ids=[11]
	    for i in range(len(ids)):
		ids[i]=str(ids[i])
	    field=['EGSs']
	if cdir.find('D17') > -1:
	    ids=[1,2,3,4,5,6,13,14]#D17
	    field=['D17s']

	fieldids=[]
	for i in ids:
		s=field[0]+str(i)
		fieldids.append(s)

	for fieldid in fieldids:
		fieldid=str(fieldid)
		print "running getcatalogs for ",fieldid
		s='n*'+str(fieldid)+'-*.fits'
		images=glob.glob(s)
		print images
		os.system('cp /Users/rfinn/field/data/sexfiles/* .')
		#noise2header(images)
		#addfwhm(images)
		alignimages(fieldid)#align using wcs
		if fieldid.find('EGSs28') > -1:
		    alignimagesv28(fieldid)
		    print "running 1113 only for EGSs28"
		    getcatalogsassoc(fieldid,'1113')#use 1184 as detection image
		else:
		    alignimagesv2(fieldid)#align w/geotran
		    getcatalogsassoc(fieldid,'1184')#use 1184 as detection image
		    getcatalogsassoc(fieldid,'1113')#use 1184 as detection image


def finalsteps2im():
	#fieldids=['EGSs23','EGSs2']
	fieldids=['EGSs24']
	cdir=str(os.getcwd())
	if cdir.find('EGSs') > -1:
	    #ids=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]
	    ids=[7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]
	    #ids=[11]
	    for i in range(len(ids)):
		ids[i]=str(ids[i])
	    field=['EGSs']
	if cdir.find('D17') > -1:
	    ids=[1,2,3,4,5,6,13,14]#D17
	    field=['D17s']

	fieldids=[]
	for i in ids:
		s=field[0]+str(i)
		fieldids.append(s)

	for fieldid in fieldids:
		fieldid=str(fieldid)
		print "running getcatalogs for ",fieldid
		s='n*'+str(fieldid)+'-*.fits'
		images=glob.glob(s)
		print images
		os.system('cp /Users/rfinn/field/data/sexfiles/* .')
		alignimages(fieldid)
		if fieldid.find('EGSs28') > -1:
		    print "running 1113 only for EGSs28"
		    getcatalogs(fieldid,'1113')#use 1184 as detection image
		else:
		    getcatalogs(fieldid,'1184')#use 1184 as detection image
		    getcatalogs(fieldid,'1113')#use 1184 as detection image


checksex=0.
if checksex > 0.:
	im='cfdq-EGSs23-J.fits'
	auto=0.
	runsextractor(im,auto)
	checksexcoords(im)
def noise2header(images):
	for i in range(len(images)):
		im=images[i]
		s='sex '+str(im)+'\n'
		os.system(s)
		(a,b)=getfieldnoise.getnoise('check.fits','test.cat')
		iraf.hedit(image=im,fields='ANOISE',value=a,add='yes',verify='no',update='yes')#add zp to image header
		iraf.hedit(image=im,fields='BNOISE',value=b,add='yes',verify='no',update='yes')#add zp to image header
def alignimagesv2(fieldid):
    print "aligning images for ",fieldid
    s='n*'+str(fieldid)+'-*.fits'
    images=glob.glob(s)
    #wcsalign(images)

    s='sncq'+fieldid+'-1184*.fits'
    im=glob.glob(s)
    im1=im[0]

    s='sncq'+fieldid+'-1113*.fits'
    im=glob.glob(s)
    im2=im[0]


    s='sncq'+fieldid+'-J*.fits'
    im=glob.glob(s)
    try:
	im3=im[0]
	flagj=1
    except IndexError:
	flagj=0

    #im2='sncq'+fieldid+'-1113.fits'
    iraf.imgets(image=im1,param='ZEROPOINT')#get RA of image in deg
    t=iraf.imgets.value
    zp=float(t)
    fwhm='1.6'#use median seeing for now
    print "running sextractor on ",im1
    filter='1184'
    s='sex '+str(im1)+' -MAG_ZEROPOINT '+str(zp)+' -SEEING_FWHM '+str(fwhm)
    os.system(s)
    infile=open('test.cat','r')#keep only sources w/in useable field
    f=fieldid+'-'+filter+'-xy.cat'
    outfile=open(f,'w')
    for line in infile:#keep only targets w/in 500 pixels of center
	if line.find('#') > -1:
	    continue
	t=line.split()
	flag=checkposition(float(t[10]),float(t[11]),filter)#check to see if (x,y) is w/in useable field
	if flag > .1:
	    out='%8.4f %8.4f \n'%(float(t[10]),float(t[11]))
	    outfile.write(out)
    infile.close()
    outfile.close()

    print "running sextractor on ",im2
    filter='1113'
    s='sex '+str(im2)+' -MAG_ZEROPOINT '+str(zp)+' -SEEING_FWHM '+str(fwhm)
    os.system(s)
    infile=open('test.cat','r')#keep only sources w/in useable field
    f2=fieldid+'-'+filter+'-xy.cat'
    outfile=open(f2,'w')
    for line in infile:#keep only targets w/in 500 pixels of center
	if line.find('#') > -1:
	    continue
	t=line.split()
	flag=checkposition(float(t[10]),float(t[11]),filter)#check to see if (x,y) is w/in useable field
	if flag > .1:
	    out='%8.4f %8.4f \n'%(float(t[10]),float(t[11]))
	    outfile.write(out)
    infile.close()
    outfile.close()

    match=fieldid+'-match'
    dbase=fieldid+'-match-geomap'
    iraf.xyxymatch(f2,f,match,tolerance=10,refpoints='refpoints-1113-1184',interactive='no')
    #iraf.xyxymatch(f2,f,match,tolerance=15,interactive='yes')
    iraf.geomap(input=match,database=dbase,xmin=1.,xmax=1024.,ymin=1.,ymax=1024.,transform='1113',interactive='no')
    outim2='sg'+im2
    iraf.geotran(im2,output=outim2,database=dbase,transform='1113')
    iraf.display(im2,1)
    iraf.display(outim2,2)
    iraf.display(im1,3)
    s='mv '+im2+' oldshifts/'
    os.system(s)
    iraf.imrename(outim2,im2)


    if flagj > 0:
	print "running sextractor on ",im3
	filter='J'
	s='sex '+str(im3)+' -MAG_ZEROPOINT '+str(zp)+' -SEEING_FWHM '+str(fwhm)
	os.system(s)
	infile=open('test.cat','r')#keep only sources w/in useable field
	f3=fieldid+'-'+filter+'-xy.cat'
	outfile=open(f3,'w')
	for line in infile:#keep only targets w/in 500 pixels of center
	    if line.find('#') > -1:
		continue
	    t=line.split()
	    flag=checkposition(float(t[10]),float(t[11]),filter)#check to see if (x,y) is w/in useable field
	    if flag > .1:
		out='%8.4f %8.4f \n'%(float(t[10]),float(t[11]))
		outfile.write(out)
	infile.close()
	outfile.close()

	match=fieldid+'-match-J'
	dbase=fieldid+'-match-geomap'
	iraf.xyxymatch(f3,f,match,tolerance=5,refpoints='refpoints-1113-1184',interactive='no')
	iraf.geomap(input=match,database=dbase,xmin=1.,xmax=1024.,ymin=1.,ymax=1024.,transform='J',interactive='no')
	outim3='sg'+im3
	iraf.geotran(im3,output=outim3,database=dbase,transform='J')
	s='mv '+im3+' oldshifts/'
	os.system(s)
	iraf.imrename(outim3,im3)

def alignimagesv28(fieldid):#for EGSs28
    s='sncq'+fieldid+'-1113*.fits'
    im=glob.glob(s)
    im2=im[0]

    s='sncq'+fieldid+'-J*.fits'
    im=glob.glob(s)
    im3=im[0]

    #im2='sncq'+fieldid+'-1113.fits'
    iraf.imgets(image=im2,param='ZEROPOINT')#get RA of image in deg
    t=iraf.imgets.value
    zp=float(t)
    fwhm='1.6'#use median seeing for now
    print "running sextractor on ",im2
    filter='1113'
    s='sex '+str(im2)+' -MAG_ZEROPOINT '+str(zp)+' -SEEING_FWHM '+str(fwhm)
    os.system(s)
    infile=open('test.cat','r')#keep only sources w/in useable field
    f2=fieldid+'-'+filter+'-xy.cat'
    outfile=open(f2,'w')
    for line in infile:#keep only targets w/in 500 pixels of center
	if line.find('#') > -1:
	    continue
	t=line.split()
	flag=checkposition(float(t[10]),float(t[11]),filter)#check to see if (x,y) is w/in useable field
	if flag > .1:
	    out='%8.4f %8.4f \n'%(float(t[10]),float(t[11]))
	    outfile.write(out)
    infile.close()
    outfile.close()

    print "running sextractor on ",im3
    filter='J'
    s='sex '+str(im3)+' -MAG_ZEROPOINT '+str(zp)+' -SEEING_FWHM '+str(fwhm)
    os.system(s)
    infile=open('test.cat','r')#keep only sources w/in useable field
    f3=fieldid+'-'+filter+'-xy.cat'
    outfile=open(f3,'w')
    for line in infile:#keep only targets w/in 500 pixels of center
	if line.find('#') > -1:
	    continue
	t=line.split()
	flag=checkposition(float(t[10]),float(t[11]),filter)#check to see if (x,y) is w/in useable field
	if flag > .1:
	    out='%8.4f %8.4f \n'%(float(t[10]),float(t[11]))
	    outfile.write(out)
    infile.close()
    outfile.close()
    
    match=fieldid+'-match-J'
    dbase=fieldid+'-match-geomap'
    iraf.xyxymatch(f3,f2,match,tolerance=5,refpoints='refpoints-1113-1184',interactive='no')
    #iraf.xyxymatch(f2,f,match,tolerance=15,interactive='yes')
    iraf.geomap(input=match,database=dbase,xmin=1.,xmax=1024.,ymin=1.,ymax=1024.,transform='J',interactive='no')
    outim3='sg'+im3
    iraf.geotran(im3,output=outim3,database=dbase,transform='J')
    s='mv '+im3+' oldshifts/'
    os.system(s)
    iraf.imrename(outim3,im3)

def alignimages(fieldid):
	s='n*'+str(fieldid)+'-*.fits'
	images=glob.glob(s)
        wcsalign(images)
	s1='s'+s
	images=glob.glob(s1)
	#trimimages(images)


def getcatalogs(fieldid,filter):
	auto=1.
	s='s*'+str(fieldid)+'-'+str(filter)+'*.fits'#use 1184 as detection image
	detectim=glob.glob(s)
	#print s,detectim,fieldid
	detectim=detectim[0]
	os.system('cp /Users/rfinn/field/data/sexfiles/* .')
	#s='cp /Users/rfinn/field/data/sexfiles/default.sex default.sex.'+str(fieldid)+'-'+str(filter)
	s='s*'+str(fieldid)+'-*.fits'#use 1184 as detection image
	images=glob.glob(s)
	runsextractor2im(detectim,images,filter,fieldid)

def getcatalogsassoc(fieldid,filter):
	auto=1.
	s='sn*'+str(fieldid)+'-'+str(filter)+'*.fits'#detection image
	detectim=glob.glob(s)
	detectim=detectim[0]
	os.system('cp /Users/rfinn/field/data/sexfiles/* .')
	s='sn*'+str(fieldid)+'-*.fits'#use 1184 as detection image
	images=glob.glob(s)
	#runsextractorassoc(detectim,images,filter,fieldid)
	runsextractorphot(detectim,images,filter,fieldid)





#calibratephot()
finalsteps()
#runsextractorJ()
#runsextractor('ncqEGSs10-1113-n3.fits',0)
#finalsteps2im()


