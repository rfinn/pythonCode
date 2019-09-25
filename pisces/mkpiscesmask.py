#!/usr/bin/env python

#writing some scripts to aid in the geometric distortion correction 

import Numeric as N
from pyraf import iraf

iraf.images()
iraf.images.imutil()

#search y=0 - 300, then flip regions to 3 other quadrants based on circle
#centered at (512,512) with a radius=555
#set pixels to 0 if N.sqrt((x-512)**2 + (y-512)**2) > 555

imxmax=1024.
imymax=1024.
#circle parameters
xc=512.
yc=512.
r=555.
#template image
template='qopen0004.fits'
try:
    iraf.imarith(template,'/',template,'mask')
except:
    iraf.imdelete('mask')
    iraf.imarith(template,'/',template,'mask')


for y in range(1,yc):
    xmax=int(512. - N.sqrt(r**2-(1.*y-yc)**2))
    if xmax < 1.:
	break
    s='mask[1:'+str(xmax)+","+str(y)+":"+str(y)+"]"
    print y,s
    iraf.imreplace(s,0.)
    x1=int(imxmax-xmax)
    x2=int(imxmax)
    s='mask['+str(x1)+':'+str(x2)+','+str(y)+':'+str(y)+']'
    print y,s
    iraf.imreplace(s,0.)
    y1=int(imymax-y+1)
    y2=int(imymax)
    s='mask[1:'+str(xmax)+","+str(y1)+":"+str(y1)+"]"
    print y,s
    iraf.imreplace(s,0.)
    s='mask['+str(x1)+':'+str(x2)+','+str(y1)+':'+str(y1)+']'
    print y,s
    iraf.imreplace(s,0.)

iraf.display('mask',2)
