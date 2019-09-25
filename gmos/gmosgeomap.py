#!/usr/bin/env python
"""
writing some scripts to aid in the geometric distortion correction 
to register I-band ediscs images to gmos pre-imaging

run sextractor on gmos preimage and vlt i-band image

create refpoints file to feed into xyxymatch

match coordinates from sextractor and run geotran

run geomap to make transformed image

get ra and dec of pixel coords using gmos image header

writes file clXXXXmastertable24GMOS.dat

give this to Vandana and she will generate fits table 

use fits table to design Nod & Shuffle mask

"""

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
iraf.imcoords()


RAcenter=132.82483
DECcenter=11.8000
#dimensions of pisces image
imxmax=1024
imymax=1024


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

def alignimages(im1,im2):
    runsex=1#run sextractor to get object catalogs
    auto=0.
    t=im1.split('.')
    prefixim1=t[0]
    t=prefixim1.split('pre')
    fieldid=t[0]#saves cl1018 or whatever as fieldid
    t=im2.split('.')
    prefixim2=t[0]

    filexy1=prefixim1+'-xy.cat'
    filexy2=prefixim2+'-xy.cat'
    print "aligning images for ",fieldid
    if runsex > 0:
	
	runsextractor(im1,auto)
	infile=open('test.cat','r')#keep only sources w/in useable field

	outfile=open(filexy1,'w')
	for line in infile:#keep only targets w/in 500 pixels of center
	    if line.find('#') > -1:
		continue
	    t=line.split()
	    out='%8.4f %8.4f \n'%(float(t[10]),float(t[11]))
	    outfile.write(out)
	infile.close()
	outfile.close()

	print "running sextractor on ",im2
	runsextractor(im2,auto)

	infile=open('test.cat','r')#keep only sources w/in useable field

	outfile=open(filexy2,'w')
	for line in infile:#keep only targets w/in 500 pixels of center
	    if line.find('#') > -1:
		continue
	    t=line.split()
	    out='%8.4f %8.4f \n'%(float(t[10]),float(t[11]))
	    outfile.write(out)
	infile.close()
	outfile.close()

    match=fieldid+'-match'
    dbase=fieldid+'-match-geomap'
    iraf.display(im1,1,fill='yes')
    iraf.display(im2,2,fill='yes')
    reffile=fieldid+'refpoints'#make these individually.  First line w/3points in gmos image (x1 y1 x2 y2 x3 y3). 2nd line w/same points in ediscs I-band
    flag= raw_input("Should we create a refpoints file (0=no, 1=yes)?\n")
    flag=int(flag)
    if flag > 0.1:
	s='cp cl1018refpoints '+reffile
	os.system(s)
	print "open ",reffile," in emacs, delete 2 lines"
	print "Enter x1 y1 x2 y2  x3 y3 for ref image (gmos)"
	print "Enter x1 y1 x2 y2  x3 y3 for 2nd image (vlt)"
	iraf.imexam()

    iraf.xyxymatch(filexy2,filexy1,match,tolerance=15,refpoints=reffile,matching='tolerance',interactive='no')
    iraf.imgets(image=im1,param='i_naxis1')
    t=iraf.imgets.value
    im1xmax=(float(t))
    iraf.imgets(image=im1,param='i_naxis2')
    t=iraf.imgets.value
    im1ymax=(float(t))

    iraf.geomap(input=match,database=dbase,function='polynomial',xmin=1,xmax=im1xmax,ymin=1,ymax=im1ymax,xxorder=4,xyorder=4,yyorder=4,yxorder=4,transform='gmos',interactive='yes')
    transformim=1
    if transformim > 0:
	outim2='g'+im2
	iraf.geotran(im2,output=outim2,database=dbase,transform='gmos')
	iraf.display(im2,1,fill='yes')
	iraf.display(outim2,2,fill='yes')
	iraf.display(im1,3,fill='yes')

def transformcoords(im1,im2):
    checkxy=1
    infile='/Users/rfinn/clusters/spitzer/Mastertables/'+prefix+'mastertable24.dat'
    input=open(infile,'r')
    i=0
    k=0
    era=[]#ediscs ra
    edec=[]#ediscs dec
    excorr=[]#ediscs x pixel position
    eycorr=[]#ediscs y pixel position
    for line in input:
	if line.find('#') > -1: #skip lines with '#' in them
	    continue
	if line.find('\\') > -1: #skip lines with '#' in them
	    continue
	if line.find('|') > -1: #skip lines with '#' in them
	    continue
	t=line.split()#break up columns & save ra and dec
	era.append(float(t[1]))
	edec.append(float(t[2]))
	excorr.append(float(t[3]))
	eycorr.append(float(t[4]))
    input.close()
    era=N.array(era,'d')
    edec=N.array(edec,'d')
    excorr=N.array(excorr,'d')
    eycorr=N.array(eycorr,'d')
    output=prefix+'radec.dat'
    out1=open(output,'w')
    for i in range(len(era)):
	s="%14.10f %14.10f \n"%(era[i],edec[i])
	out1.write(s)
    out1.close()
    output=prefix+'-ediscsxy.dat'#now contains xcorr, ycorr
    out1=open(output,'w')
    for i in range(len(era)):
	s="%14.10f %14.10f \n"%(excorr[i],eycorr[i])
	out1.write(s)
    out1.close()
    

    #transform ra and dec to x and y coords on ediscs image using ediscs i band image
    #write era and edec to a file
    outfile=output
    #outfile=prefix+'-ediscsxy.dat'#file containing x and y pixels values on ediscs iband image
    #iraf.wcsctran(output,outfile,image=im2,inwcs='world',outwcs='physical',columns='1 2', verbose='no')

    if checkxy > 0:
	iraf.display(im2,1,fill='yes')
	iraf.tvmark(1,outfile,color=204,radii=15)
	#iraf.tvmark(1,output,color=204,radii=15)#plot xcorr and ycorr instead of transformed coordinates
	flag= raw_input("Enter 0 to terminate now or 1 to continue \n")	
	if flag < 0.1:
	    return

    dbase=prefix+'-match-geomap'
    outfile2=prefix+'-GMOSxy.dat'#file containing x and y pixels values on gmos image
    iraf.geoxytran(outfile,output=outfile2,database=dbase,transform='gmos', direction='backward')

    if checkxy > 0:
	iraf.display(im1,2,fill='yes')
	iraf.tvmark(2,outfile2,color=204,radii=15)
	im3='g'+im2
	iraf.display(im3,3,fill='yes')
	iraf.tvmark(3,outfile2,color=204,radii=15)
	flag= raw_input("Enter anything to continue \n")	

def writegmoscat(gmosim):
    iraf.imgets(image=gmosim,param='i_naxis1')
    t=iraf.imgets.value
    gxmax=(float(t))
    iraf.imgets(image=gmosim,param='i_naxis2')
    t=iraf.imgets.value
    gymax=(float(t))


    gxmin=1
    gymin=1

    outfile2=prefix+'-GMOSxy.dat'#file containing x and y pixels values on gmos image
    gmosx=[]
    gmosy=[]
    input=open(outfile2,'r')
    for line in input:
	t=line.split()
	gmosx.append(float(t[0]))
	gmosy.append(float(t[1]))
    input.close()
    gmosx=N.array(gmosx,'f')
    gmosy=N.array(gmosy,'f')

    #get gmos ra and dec positions from x,y pixels values
    radecfile=prefix+'-GMOSradec.dat'
    iraf.wcsctran(outfile2,radecfile,image=gmosim,inwcs='physical',outwcs='world',columns='1 2', verbose='no')
    infile=open(radecfile,'r')
    gmosra=[]
    gmosdec=[]
    for line in infile:
	    if line.find('#') > -1:
		    continue
	    t=line.split()
	    gmosra.append(float(t[0]))
	    gmosdec.append(float(t[1]))
    infile.close()
    gmosra=N.array(gmosra,'f')
    gmosdec=N.array(gmosdec,'f')

    s=prefix+'.reg'
    dsfile=open(s,'w')
    s="global color=green font='helvetica 10 normal' select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n"
    dsfile.write(s)
    s="fk5 \n"
    dsfile.write(s)

    infile='/Users/rfinn/clusters/spitzer/Mastertables/'+prefix+'mastertable24.dat'
    gmoscat=prefix+'mastertable24GMOS.dat'
    input=open(infile,'r')
    output=open(gmoscat,'w')
    i=0
    for line in input:
	if line.find('#') > -1: #skip lines with '#' in them
	    output.write(line)

	    continue
	if line.find('\\') > -1: #skip lines with '#' in them
	    output.write(line)
	    continue
	if line.find('|') > -1: #skip lines with '#' in them
	    output.write(line)
	    continue
	if gmosx[i] < 1.:
	    print "bad x value: ",gmosx[i],1,gxmax
	    i=i+1
	    continue
	if gmosx[i] > gxmax:
	    print "bad x value: ",gmosx[i],1,gxmax
	    i=i+1
	    continue
	if gmosy[i] < 1.:
	    print "bad y value (too small): ",i,gmosy[i],1,gymax
	    i=i+1
	    continue
	if gmosy[i] > gymax:
	    print "bad y value: ",i,gmosy[i],1,gymax
	    i=i+1
	    continue
	line=line[0:(len(line)-1)]#get rid of \n at end of line 
	s=" %10.3f %10.3f %13.8f %13.8f \n"%(gmosx[i],gmosy[i],gmosra[i],gmosdec[i])
	outline=line+s#append gmos x and y pixels to the end of the line
	output.write(outline)
	string1 = "circle(%12.8f, %12.8f, 3\") \n"%(gmosra[i],gmosdec[i])
	dsfile.write(string1)
	i=i+1
    input.close()

    dsfile.close()


def prelimstuff():#pull out ext 1 as science image
    images=glob.glob('*preimage.fits')
    for im in images:
	input = im+'[1]'
	prefix=im.split('.')
	pre=prefix[0]
	newname=pre+'sci.fits'
	iraf.imcopy(input,newname)
	

auto=0.

#prelimstuff()
#os.system('cp /Users/rfinn/field/data/sexfiles/default.sex.2mass default.sex')
#os.system('cp /Users/rfinn/field/data/sexfiles/default.param .')
#os.system('cp /Users/rfinn/field/data/sexfiles/gauss_3.0_5x5.conv .')
#os.system('cp /Users/rfinn/field/data/sexfiles/default.nnw .')


#gmosimages=glob.glob('*preimagesci.fits')
#edisimages=glob.glob('*i.fss_avg.fits')

gmosimages=[]
edisimages=[]

inlist=open('vltlist','r')
for line in inlist:
    if line.find('#') > -1:
	continue
    t=line.split()
    gmosimages.append(t[0])
    edisimages.append(t[1])

#for i in range(len(gmosimages)):
#for i in range(1,len(gmosimages)):
#for i in range(2,3):
#for i in range(1,2):#cl1037
for i in range(len(gmosimages)):
    t=gmosimages[i]#get prefix from gmos filename
    f=t.split('.')
    f=f[0]
    p=f.split('pre')
    prefix=p[0]
    print "Working on cluster ",prefix
    #alignimages(gmosimages[i],edisimages[i])
    #transformcoords(gmosimages[i],edisimages[i])
    writegmoscat(gmosimages[i])

