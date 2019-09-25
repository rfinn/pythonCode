#!/usr/bin/env python
"""
run from /Users/rfinn/clusters/spitzer/GroupLirgs

"""

from pylab import *
from pyraf import iraf
import pyfits
import glob


delta=10.#width of cutouts in arcsec

def getxyfromradec(im,ra,dec):
    iraf.imgets(image=im,param='CRVAL1')#get RA of image in deg
    t=iraf.imgets.value
    RAcenter=float(t)
    iraf.imgets(image=im,param='CRVAL2')#get RA of image
    t=iraf.imgets.value
    DECcenter=float(t)
    dra=(RAcenter-ra)*cos(pi/180.*dec)#correct delta ra for cos declination
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


infile=open('/Users/rfinn/clusters/spitzer/MasterTables/GroupLirgs.dat','r')
ids=[]
for line in infile:
    t=line.split()
    ids.append(t[0])
infile.close()


infile=open('/Users/rfinn/clusters/spitzer/MasterTables/GroupLirgs.reg','r')
ra=[]
dec=[]
o2flag=[]
for line in infile:
    if line.find('circle')>-1:
	ra.append(float(line[7:18]))
	dec.append(float(line[21:32]))
	if line.find('blue')>-1:
	    o2flag.append(1)
	else:
	    o2flag.append(0)
infile.close()
ra=array(ra,'f')
dec=array(dec,'f')

images=[]
images24=[]
rimages24=[]
flag=ones(len(ra),'f')
i=0
out1=open('cl1040lirgs.radec','w')
out2=open('cl1103lirgs.radec','w')
out3=open('cl1054-12lirgs.radec','w')
out4=open('cl1054-11lirgs.radec','w')
images.append('/Users/rfinn/clusters/ediscs/HST-images/cl1040_drz_gsci.fits')
images.append('/Users/rfinn/clusters/ediscs/HST-images/cl1103_drz_gsci.fits')
images.append('/Users/rfinn/clusters/ediscs/HST-images/cl1054-12_drz_gsci.fits')
images.append('/Users/rfinn/clusters/ediscs/HST-images/cl1054-11_drz_gsci.fits')
path24='/Users/rfinn/clusters/spitzer/final-images/'
prefix=['cl1040','cl1103','cl105412','cl105411']
for pre in prefix:
    s=path24+pre+'_final24.fits'
    images24.append(s)
    s=path24+'r'+pre+'_final24.fits'
    rimages24.append(s)
incoords=['cl1040lirgs.radec','cl1103lirgs.radec','cl1054-12lirgs.radec','cl1054-11lirgs.radec']
prefix=['cl1040','cl1103','cl1054-12','cl1054-11']
for i in range(len(ids)):
    name=ids[i]
    s='%f %f %s\n'%(ra[i],dec[i],ids[i])
    if name.find('J1040') > -1:
	out1.write(s)
    if name.find('J110') > -1:
	out2.write(s)
    if name.find('J1301') > -1:
	flag[i]=0
    if name.find('J1054') > -1:
	    if name.find('-12') > -1:
		out3.write(s)
	    else:
		out4.write(s)
    i=i+1
out1.close()
out2.close()
out3.close()
out4.close()

def transcoords():
    for i in range(len(images)):
	print i,prefix,prefix[i]
	outcoords=str(prefix[i])+'.xy'
	print images[i],incoords[i],outcoords
	iraf.imcoords.wcsctran(image=images[i],input=incoords[i],output=outcoords,inwcs='world',outwcs='logical')



def makecutouts():
    cutouts=[]
    for i in range(len(images)):
    #for i in range(1):
	#print i,prefix,prefix[i]
	outcoords=str(prefix[i])+'.xy'

	iraf.imgets(image=images[i],param='CD2_1')#get x value corresponding to RA 
	xplate=abs(float(iraf.imgets.value))#deg/pixel
	dpix=delta/3600./xplate/2.
	print dpix
	#delta=100.

	infile=open(outcoords,'r')
	for line in infile:
	    if line.find('#') > -1:
		continue
	    if len(line)<2:
		continue
	    x,y,id=line.split()
	    x=float(x)
	    y=float(y)

	    xmin=int(round(x-dpix))
	    xmax=int(round(x+dpix))
	    ymin=int(round(y-dpix))
	    ymax=int(round(y+dpix))
	    s=images[i]+'[%i:%i,%i:%i]'%(xmin,xmax,ymin,ymax)
	    print s
	    #print ra[i],dec[i]
	    outim=id+'cutout.fits'
	    print outim
	    iraf.imcopy(s,outim)
	    cutouts.append(outim)
	infile.close()
    return cutouts

def transcoords24():
    images24=rimages24
    for i in range(len(images24)):
    #for i in range(1):
	print i,prefix,prefix[i]
	outcoords=str(prefix[i])+'.xy24'
	print images24[i],incoords[i],outcoords
	iraf.imcoords.wcsctran(image=images24[i],input=incoords[i],output=outcoords,inwcs='world',outwcs='logical')


def makecutouts24(delta):
    cutouts24=[]
    images24=rimages24
    for i in range(len(images24)):
    #for i in range(1):
	#print i,prefix,prefix[i]
	outcoords=str(prefix[i])+'.xy24'
	#iraf.imgets(image=images24[i],param='CDELT1')#get x plate scale
	iraf.imgets(image=images24[i],param='CD2_1')#get x plate scale on rotated image
	print iraf.imgets.value
	xplate=abs(float(iraf.imgets.value))#deg/pixel

	dpix=delta/3600./xplate/2.

	iraf.imgets(image=images24[i],param='naxis1')#get x value corresponding to RA 
	xpixmax=(int(iraf.imgets.value))#deg/pixel
	iraf.imgets(image=images24[i],param='naxis2')#get x value corresponding to RA 
	ypixmax=(int(iraf.imgets.value))#deg/pixel


	infile=open(outcoords,'r')
	#print outcoords
	for line in infile:
	    #print images24[i],line,outcoords
	    if line.find('#') > -1:
		continue
	    if len(line)<2:
		continue
	    x,y,id=line.split()
	    x=float(x)
	    y=float(y)
	    #delta=100.
	    xmin=int(round(x-dpix))
	    xmax=int(round(x+dpix))
	    ymin=int(round(y-dpix))
	    ymax=int(round(y+dpix))
	    if xmin < 1:
		xmin=1
	    if ymin < 1:
		ymin=1
	    if xmax > xpixmax:
		xmax=xpixmax
	    if ymax > ypixmax:
		ymax=ypixmax
	    s=images24[i]+'[%i:%i,%i:%i]'%(xmin,xmax,ymin,ymax)
	    print s
	    #print ra[i],dec[i]
	    outim=id+'cutout24.fits'
	    print outim
	    iraf.imcopy(s,outim)
	    cutouts24.append(outim)
	infile.close()
    return cutouts24



def plotcutouts():
    o2lirgs=['EDCSNJ1040533-1153541', 'EDCSNJ1040532-1153540','EDCSNJ1040532-1154044', 'EDCSNJ1104001-1242319','EDCSNJ1104001-1248465', 'EDCSNJ1103599-1242595','EDCSNJ1054418-1144511', 'EDCSNJ1054418-1145508', 'EDCSNJ1054416-1148043', 'EDCSNJ1040532-1155482']
    clf()
    cla()
    figure(figsize=(12,16))
    subplots_adjust(left=0.1, right=.95,bottom=.1,top=0.95,wspace=0.001,hspace=0.001)
    hstfiles=glob.glob('EDCS*cutout.fits')
    mipsfiles=glob.glob('EDCS*cutout24.fits')
    nx=6
    ny=8
    nplot=0
    for i in range(len(hstfiles)):
	name=hstfiles[i]
	if name.find('EDCSNJ1104001-1243540')>-1:#not on HST image
	    continue
	if name.find('EDCSNJ1104000-1247179')>-1:#no real 24um source
	    continue
	if name.find('EDCSNJ1104000-1243555')>-1:#no real 24um source
	    continue
	for j in range(2):
	    #k=j+2*i+1
	    nplot=nplot+1
	    k=nplot
	    print k
	    subplot(ny,nx,nplot)
	    if (j == 0):
		fits=pyfits.open(hstfiles[i])
		name=hstfiles[i]
		t=name.split('cutout')
		im=fits[0].data.copy()

		#ave=average(average(im))
		ave=median(median(im))
		stddev=std(std(im))
		nave=5
		#print max(max(im))
		if name.find('EDCSNJ1104002-1248394')>-1:
		    print 'found the rogue!'
		    nave=.8
	
		#im[where(im>(ave+nave*stddev))]=(ave+nave*stddev)
		lim=.07
		im[where(im>(lim))]=lim
		nave=3.5
		#im[where(im<(ave-nave*stddev))]=(ave-nave*stddev)
		im[where(im<0)]=0
		myvmin=None
		myvmax=None


	    else:
		fits=pyfits.open(mipsfiles[i])
		name=mipsfiles[i]
		im=fits[0].data.copy()
		myvmin=0.01
		myvmax=.03
		im[where(im<0)]=0
		im[where(im>(.08))]=.08
		axis([1.,5.,1.,5.])
	
	    #imshow(-1*log(im+1))
	    fits.close()
	    axis('equal')
	    imshow(-1.*(im),interpolation='nearest',origin='upper',cmap='gray')#,vmin=myvmin,vmax=myvmax)
	    #(ymin,ymax)=ylim()
	    #ylim(ymax,ymin)
	    ax=gca()
	    ax.set_xticklabels(([]))
	    ax=gca()
	    ax.set_yticklabels(([]))
	    id=t[0]
	    if id in o2lirgs:
		col='b'
		#text(t[0],fontsize=4,color='b')
	    else:
		#title(t[0],fontsize=4,color='k')
		col='k'
	    #print col
	    junk,id=id.split('J')
	    text(.5,.82,id,fontsize=10,transform=ax.transAxes,color=col,horizontalalignment='center')
    savefig('AllCutouts.eps')

def rotate24():

    for i in range(len(images24)):
    #for i in range(1):
	iraf.imgets(image=images24[i],param='CROTA2')#[deg] Orientation of axis 2 (W of N, +=CW)
	print iraf.imgets.value
	rot24=float(iraf.imgets.value)
	iraf.imgets(image=images[i],param='ORIENTAT')#position angle of image y axis (deg. e of n)
	print iraf.imgets.value
	rotacs=float(iraf.imgets.value)
	angle=-1*rot24-rotacs
	outimage=rimages24[i]
	iraf.rotate(input=images24[i],output=outimage,rotation=angle)

#rotate24()

#transcoords()
#cutouts=makecutouts()

#transcoords24()
#cutouts24=makecutouts24(delta)


plotcutouts()
