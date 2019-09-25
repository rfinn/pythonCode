#!/usr/bin/env python


from pylab import *

def getagc(file):
    agcnumber=[]
    awhich=[]
    ngcic=[]
    rah=[]
    ram=[]
    ras10=[]
    decd=[]
    decm=[]
    decs=[]
    a100=[]
    b100=[]
    mag10=[]
    bsteintype=[]
    vopt=[]
    vsource=[]
    flux100=[]
    rms100=[]
    v21=[]
    width=[]
    widthcode=[]
    telcode=[]
    detcode=[]
    hisource=[]
    ibandqual=[]
    ibandsrc=[]
    posang=[]
    iposition=[]
    ipalomar=[]

    #file='gals.txt'
    infile=open(file,'r')
    for line in infile:
    #idl code that writes gals.txt
#
#                    printf, gals_txt_lun, agccat.AGCNUMBER[i],agccat.WHICH[i],agccat.NGCIC[i],agccat.RAH[i],agccat.RAM[i],agccat.RAS10[i],agccat.SIGN[i],$
#                      agccat.DECD[i],agccat.DECM[i],$
#                      agccat.DECS[i],agccat.A100[i],agccat.B100[i],agccat.MAG10[i],agccat.BSTEINTYPE[i],$
#                      agccat.VOPT[i],agccat.VSOURCE[i],agccat.FLUX100[i],agccat.RMS100[i],$
#                      agccat.V21[i],agccat.WIDTH[i],agccat.WIDTHCODE[i],agccat.TELCODE[i],agccat.DETCODE[i],agccat.HISOURCE[i],$
#                      agccat.IBANDQUAL[i],agccat.IBANDSRC[i],agccat.POSANG[i],agccat.IPOSITION[i],agccat.IPALOMAR[i],$
#                      format = '(I6,A1,A8,1x,2I2.2,I3.3,A1,3I2.2,I5,2I4,i6,2x,I5,I3,2x,I5,I7,1x,I7,I4,1x,A4,A4,i1,1x,I3,2i3,i5,2i3)'

	agcnumber.append(line[0:7])
	awhich.append(line[6:7])
	ngcic.append(line[7:15])
	rah.append(float(line[16:18]))
	ram.append(float(line[18:20]))
	ras10.append(float(line[20:23]))
	decd.append(float(line[24:26]))
	decm.append(float(line[26:28]))
	decs.append(float(line[28:30]))
	a100.append(float(line[30:35]))
	b100.append(float(line[35:39]))
	mag10.append(float(line[39:43]))
	bsteintype.append(line[43:49])
	vopt.append(float(line[51:56]))
	vsource.append(float(line[56:59]))
	flux100.append(float(line[61:66]))
	rms100.append(float(line[66:73]))
	v21.append(float(line[74:81]))
	width.append(float(line[81:85]))
	widthcode.append(line[86:90])
	telcode.append(line[90:94])
	detcode.append(line[94:95])
	hisource.append(line[96:99])
	ibandqual.append(line[99:102])
	ibandsrc.append(line[102:105])
	posang.append(line[105:110])
	iposition.append(line[110:113])
	ipalomar.append(line[113:])


#print rah
    rah=array(rah,'f')
    ram=array(ram,'f')
    ras10=array(ras10,'f')
    decd=array(decd,'f')
    decm=array(decm,'f')
    decs=array(decs,'f')
    radeg=15*(rah+ram/60.+ras10/3600./10.)
    decdeg=decd+decm/60.+decs/3600.
    vopt=array(vopt,'f')
    v21=array(v21,'f')


    return(radeg,decdeg,vopt,v21)

#cla()
#clf()
#plot(radeg,decdeg,'r+')
#savefig('test.eps')

#clf()
#cla()
#a=vopt[where(vopt > 5000.)]
#print a, vopt[where(vopt > 5000.)]
#hist(a,20,facecolor='r')
#b=v21[where(v21 > 5000.)]
#hist(b,20)
#savefig('testhist.eps')


def transcoords(prefix,imagein,incoords):#get x,y from ra,dec
    outcoords=str(prefix)+'.xy'
    iraf.imcoords.wcsctran(image=imagein,input=incoords,output=outcoords,inwcs='world',outwcs='logical')
    return outcoords

def makecutouts(imagein,x,y,delta):#image, array containing x coords, array containing y coords
    cutouts=[]
    iraf.imgets(image=imagein,param='CD2_1')#get x value corresponding to RA 
    xplate=abs(float(iraf.imgets.value))#deg/pixel
    dpix=delta/3600./xplate/2.

    iraf.imgets(image=imagein,param='naxis1')#get x value corresponding to RA 
    xpixmax=(int(iraf.imgets.value))#deg/pixel
    iraf.imgets(image=imagein,param='naxis2')#get x value corresponding to RA 
    ypixmax=(int(iraf.imgets.value))#deg/pixel

    for i in range(len(x)):
	xmin=int(round(x[i]-dpix))
	xmax=int(round(x[i]+dpix))
	ymin=int(round(y[i]-dpix))
	ymax=int(round(y[i]+dpix))
	if xmin < 1:
	    xmin=1
	if ymin < 1:
	    ymin=1
	if xmax > xpixmax:
	    xmax=xpixmax
	if ymax > ypixmax:
	    ymax=ypixmax

	s=images[i]+'[%i:%i,%i:%i]'%(xmin,xmax,ymin,ymax)
	outim=id+'cutout.fits'
	print outim
	iraf.imcopy(s,outim)
	cutouts.append(outim)
    infile.close()
    return cutouts,OnImageFlag

def checklocation(imagein,x,y,delta):#image, array containing x coords, array containing y coordswidth and heigh of cutout in arcsec
#check to see if x,y is on the image
    iraf.imgets(image=imagein,param='CD2_1')#get x value corresponding to RA 
    xplate=abs(float(iraf.imgets.value))#deg/pixel
    dpix=delta/3600./xplate/2.
    iraf.imgets(image=imagein,param='naxis1')#get x value corresponding to RA 
    xpixmax=(int(iraf.imgets.value))#deg/pixel
    iraf.imgets(image=imagein,param='naxis2')#get x value corresponding to RA 
    ypixmax=(int(iraf.imgets.value))#deg/pixel

    for i in range(len(x)):
	xmin=int(round(x[i]-dpix))
	xmax=int(round(x[i]+dpix))
	ymin=int(round(y[i]-dpix))
	ymax=int(round(y[i]+dpix))
	if xmax < 1.:
	    OnImageFlag[i]=0
	    continue
	if xmin > xpixmax:
	    OnImageFlag[i]=0
	    continue
	if ymax < 1.:
	    OnImageFlag[i]=0
	    continue
	if ymin > ypixmax:
	    OnImageFlag[i]=0
	    continue
    return OnImageFlag



