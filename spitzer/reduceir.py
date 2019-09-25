#!/usr/bin/env python
"""
rewriting reduce_ir scripts in python to try to improve flatfielding technique

assumes corquad is already done

assumes darks are combined

"""

import sys, glob, os
import numarray as n
import Numeric as N
import scipy
from math import *
from pyraf import iraf
import pyfits
#import ppgplot
#pisces noise stats
READNOISE=4.
GAIN=4.35

npass = 0

iraf.images()
iraf.images.imutil()
iraf.immatch()

def combinedarks():
    darktimes= [10,120,20,30,300,360,4,45,5,60,90]
    for t in darktimes:
        prefix="qdark_"+str(t)+"????.fits"
        out="dk"+str(t)
        iraf.imcombine(prefix,output=out,combine="ave", reject="ccdclip",rdnoise=READNOISE, gain=GAIN, lsigma=3,hsigma=3)

def corquad():
    allfiles=glob.glob('*fits')
    for file in allfiles:
	s="corquad "+str(file)+" \n"
	os.system(s)

def darksubtract(files,dark):
    for image in files:
        out="d"+image
        print image,out
        iraf.imarith(image,"-",dark,out)
def domask(files,invmask):
    current=files
    #for current in files:
    mfile="m"+current
    #iraf.mask(current,mask=invmask,stat="[400:600,400:600]")
    mode=iraf.imstatistics(current,fields="mode",lower=1,format=0,Stdout=1)
    mode=float(mode[0])
    print "masking ",current, " with mode = ",mode
    iraf.imarith(1,"-",invmask,"tempmask")
    iraf.imarith(current,"*","tempmask",mfile)
    iraf.imarith(invmask,"*",mode,"tempmask2")
    iraf.imarith(mfile,"+","tempmask2",mfile)
    iraf.imdel("tempmask")
    iraf.imdel("tempmask2")

def domask999():
    files=glob.glob('fm*.fits')
    mask='mask2007.pl'#bad values = 1, good values=0
    iraf.imarith(1,"-",mask,"tempmask")#good values=1, badvalues=0
    mode=-9999.
    for current in files:
	mfile="m"+current
	iraf.imarith(current,"*","tempmask",mfile)
	iraf.imarith(mask,"*",mode,"tempmask2")
	print "masking ",current, " with mode = ",mode
	iraf.imarith(mfile,"+","tempmask2",mfile)
	iraf.imdel("tempmask2")
    iraf.imdel("tempmask")

def FlattenData():
    filters=['J','1113','1184']
    for filt in filters:
	s='m*'+filt+'*.fits'
	files=glob.glob(s)
	flat='flat'+filt
	for file in files:
	    out='f'+file
	    print file," -> ",out
	    iraf.imarith(file,'/',flat,out)
	    iraf.display(out,1)

def GeoTran():
    filters=['J','1113','1184']
    for filt in filters:
	s='mfm*'+filt+'*.fits'
	files=glob.glob(s)
	flist='flist'+str(filt)
	glist='glist'+str(filt)
	fout=open(flist,'w')
	gout=open(glist,'w')
	for file in files:
	    ffile=str(file)#+'\n'
	    gfile='g'+str(file)#+'\n'
	    fout.write(ffile+'\n')
	    gout.write(gfile+'\n')
	    iraf.geotran(input=ffile, output=gfile, database='PiscesBok',transforms=filt,fluxconserve='no')
	    iraf.rotate(input=gfile, output=gfile,rotation=90)
	    gfilein='g'+str(file)+'[-*,*]'
	    gfileout='g'+str(file)
	    iraf.imcopy(gfilein, gfileout)

	infiles='@'+str(flist)
	outfiles='@'+str(glist)
	#infiles=str(flist)
	#outfiles=str(glist)
	fout.close()
	gout.close()

def makeflats():
    filters=['J','1113','1184']
    flatfiles=['@FlatfilesJ','@Flatfiles1113','@Flatfiles1184']
    for i in range(len(filters)):
	filt=filters[i]
	files=flatfiles[i]
	#files='*'+filt+'*.fits'
	flat='flat'+filt
	iraf.imcombine(files,output=flat,combine="median",reject="minmax",scale='median',weight='exposure',statsec="[400:600,400:600]",nlow=2,nhigh=4)
	input=flat+'.fits[300:800,300:800]'
	stats=iraf.imstatistics(input,fields="mean",lower=1,format=0,Stdout=1)
	print 'stats = ',stats,stats[0]
	ave=float(stats[0])
	iraf.imarith(flat,"/",ave,flat)#normalize flat

def MoveDarks():
    os.system('mkdir darks')
    os.system('mv *dark*fits darks/.')
    os.system('mv *dk*fits darks/.')
def MoveMdata():
    try:
	os.system('mkdir mdata')
    except:
	print "I think there is already an mdata directory"
    os.system('mv m*fits mdata/.')

def MoveQdata():
    try:
	os.system('mkdir qdata')
    except:
	print "I think there is already a qdata directory"
    os.system('mv q*fits qdata/.')

def MoveQdataBack():
    os.system('mv qdata/q*fits .')

def MoveRawdata():
    os.system('mkdir rawdata')
    os.system('mv *fits rawdata/.')

def MoveStandards():
    os.system('mkdir standards')
    try:
	os.system('mv mqp*.fits standards/.')
    except:
	print "apparently no p standards"
    try:
	os.system('mv mqs*.fits standards/.')
    except:
	print "apparently no s standards"

def MoveStandardsBack():
    try:
	os.system('cp standards/* . ')
    except:
	print "didn't find any standards to move back"


def writedimbatch(i,npass,prefix):
    #print 'writedimbatch ',i,npass
    if npass == 1:
        output = open('dimbatch1'+str(i)+'.cl','w')
        dimcommand='reduce @diminlist'+str(i)+' out_img=dim'+str(npass)+str(i)+' ref_img='+reference[i]+' prefix=ss_ chk+ fp_xslm+ fp_fixpix- fp_xzap+ fp_badpixupd- mk_shifts+ fp_registar+ masking+ maskdereg+ mp_xslm- mp_fixpix- mp_xzap- mp_badpixupd- mp_registar- bpimage='+mask+' \n'
        output.write('dimsum\n')
        output.write(dimcommand)
        output.write('!cp objmask* backup/. \n')
    if npass == 2:
        output = open('dimbatch2'+str(i)+'.cl','w')
        #now rerun all of dimsum on distortion-corrected images
        dimcommand='reduce @diminlist'+str(i)+' out_img=dim'+str(npass)+str(i)+' ref_img='+reference2[i]+' prefix=ss_ chk+ fp_xslm+ fp_fixpix- fp_xzap+ fp_badpixupd- mk_shifts+ fp_registar+ masking+ maskdereg+ mp_xslm+ mp_fixpix- mp_xzap+ mp_badpixupd- mp_registar+ bpimage=g'+mask+' \n'
        output.write('dimsum\n')
        output.write(dimcommand)
        output.write('!cp objmask* backup/. \n')
        output.close()
    if i == len(prefix):
        output.write('!cp backup/objmask* . \n')
        output.close()

def sumimages(files):
    imsum = 0
    #for filename in files:
    for i in range(len(files)):
        print i,files[i]
        fitsobj = pyfits.open(files[i])
        imsum = imsum + fitsobj[0].data
    outim = pyfits.HDUList()
    outhdu = pyfits.PrimaryHDU()
    outhdu.data = imsum
    outim.append(outhdu)
#    outim.writeto('total.fits')    
    try:
        outim.writeto('total.fits')
    except:
        os.system('rm total.fits')
        outim.writeto('total.fits')
def imageprocess():
    #corquad()
    #MoveQdata()
    #MoveRawdata()
    #MoveQdataBack()
    #domask()
    #MoveQdata()
    #MoveStandards()
    #makeflats()
    #MoveStandardsBack()
    #FlattenData()
    domask999()
    GeoTran()

imageprocess()
#narrow band MS1054
#prefix=['qdata*']
#darks=['dk300.fits']
#flats=['qflatn.fits']
#flats2=['qflatn2.fits']
#reference=['fmdqdata0085.fits']
#reference2=['gfmdqdata0085.fits']
#finalname=['ms1054']

#j-band MS1054
#prefix=['qdata*']
#darks=['dk120.fits']
#flats=['qflatj.fits']
#flats2=['qflatj2.fits']
#reference=['fmdqdata0108.fits']
#reference2=['gfmdqdata0108.fits']

#CL1040 J-band
prefix=['qc104j*']
darks=['dk60.fits']
flats=['flatj.fits']
flats2=['flatj2.fits']
reference=['fmdqc104j009.fits']
reference2=['gfmdc104j009.fits']
finalname=['cl1040j']

#CL1040 narrow-band
prefix=['qc104n*']
darks=['dk600.fits']
flats=['flatn.fits']
flats2=['flatn2.fits']
reference=['fmdqc104n003.fits']
reference2=['gfmdqc104n003.fits']
finalname=['cl1040n']


#CL1054 J-band
prefix=['qc105j*']
darks=['dk120.fits']
flats=['flatj.fits']
flats2=['flatj2.fits']
reference=['fmdqc105j106.fits']
reference2=['gfmdqc105j106.fits']
finalname=['cl1054j']

#CL1054 narrow-band
prefix=['qc105n*']
darks=['dk600.fits']
flats=['flatn.fits']
flats2=['flatn2.fits']
reference=['fmdqc105n000.fits']
reference2=['gfmdqc105n000.fits']
finalname=['cl1054n']

#CL1216 narrow-band
prefix=['qc12n*']
darks=['dk300.fits']
flats=['flatn.fits']
flats2=['flatn2.fits']
reference=['fmdqc12n006.fits']
reference2=['gmfmdqc12n006.fits']
finalname=['cl1216n']

#CL1216 narrow-band
prefix=['qc12j*']
darks=['dk120.fits']
flats=['flatj.fits']
flats2=['flatj2.fits']
reference=['fmdqc12j009.fits']
reference2=['gmfmdqc12j009.fits']
finalname=['cl1216j']

#CL0152 j-band
prefix=['qj015j*']
darks=['dk120.fits']
flats=['flatj.fits']
flats2=['flatj2.fits']
reference=['fmdqj015j000.fits']
reference2=['gmfmdqj015j000.fits']
finalname=['cl0152j']

#CL0152 narrow-band
prefix=['qj015n*']
darks=['dk300.fits']
flats=['flatn.fits']
flats2=['flatn2.fits']
reference=['fmdqj015n004.fits']
reference2=['gmfmdqj015n004.fits']
finalname=['cl0152n']

#HDFN1 narrow-band 03/25/05
prefix=['qdata*']

flats=['flatn.fits']
flats2=['flatn2.fits']
#reference=['fmdqdata0062.fits'] #3/25
#reference2=['gmfmdqdata0062.fits']
darks=['/home/rfinn/field/data/data0328/dk300.fits']
reference=['fmdqdata0035.fits'] #03/26/05
reference2=['gmfmdqdata0035.fits']
finalname=['HDFN1n']
#J-band data
darks=['/home/rfinn/field/data/data0328/dk120.fits']
reference=['fmdqdata0005.fits'] #03/26/05
reference2=['gmfmdqdata0005.fits']
finalname=['HDFN1j']


mask="/home/rfinn/field/data/jammask_2004.pl"
invmask="/home/rfinn/field/data/jammask_inv.pl"

geodb="geodb4"
geotrans="mmtj"

#90-inch
geodb="/home/rfinn/field/data/J_9904_geomapout"
geotrans="J_9904_geomapout"

#combinedarks()

i=0
for name in prefix:
    print i,name
    files = glob.glob(name)

#    dfiles = 'd'+files
#    mfiles = 'md'+files
#    ffiles = 'fmd'+files
#    gfiles = 'gfmd'+files
    dname='d'+name
    mname='md'+name
    fname='fmd'+name
    gname='gfmd'+name



    writedimbatch(i,npass,prefix)
    if npass == 1:
        flat=flats[i]
        darksubtract(files,darks[i])
        for file in files:
            dfile='d'+file
        #    iraf.mask(dfile,invmask)
            domask(dfile,invmask)
        #iraf.mask('@dfiles',invmask,stat='[400:600,400:600]')
        iraf.imcombine(mname,output=flat,combine="median",reject="none",zero="none",scale="median",statsec="[400:600,400:600]",nlow="1",nhigh="2")

    if npass == 2:
        flat=flats2[i]
        print dname
        for file in files:
            dfile = 'd'+file
            infile = 'md' + file
            mfile='objmask_ss_fmd'+ file[:-4] + 'pl'
            print "mask = ",mfile
            iraf.hedit(infile,fields="BPM",value=mfile,add="yes",delete="no",verify="no",update="yes")
        #iraf.imcombine(mname,output=flat,combine="median",reject="none",zero="none",scale="median",statsec="[400:600,400:600]",masktype="goodvalue",maskvalue="0",nlow="1",nhigh="2")
        iraf.imcombine(mname,output=flat,combine="median",scale="median",stats="[400:600,400:600]", reject="ccdclip",lsigma=5.5,hsig=4.5, rdnoise=17.5,gain=4.35, maskt="goodval", lthres=500,hthres=32000,nkeep=2,blank=-999,grow=3)

    stats=iraf.imstatistics(flat,fields="mean",lower=1,format=0,Stdout=1)
    ave=float(stats[0])
    iraf.imarith(flat,"/",ave,flat)#normalize flat

    for file in files:
        dfile = 'd' + file
        mfile = 'md' + file
        ffile = 'fmd' + file
        m999file = 'mfmd' + file
        gfile = 'gmfmd' + file
        if npass == 1:
            iraf.imarith(mfile,"/",flat,ffile)
            diminput=fname
        if npass == 2:
            iraf.imdelete(ffile)
            iraf.imarith(mfile,"/",flat,ffile)
            domask999(ffile,invmask)
#            iraf.geotran(ffile,gfile,database="geodb4",trans="mmtj",geometry="geometric")


            iraf.geotran(m999file,gfile,database=geodb,trans=geotrans,geometry="geometric",interpolant="nearest",boundary="nearest",fluxconserve="no")
            diminput=gname
    makedim = 'ls '+ diminput +' > diminlist'+str(i)
    os.system(makedim)    
    i = i + 1
gmask = 'g' + mask
try:
    iraf.geotran(mask,gmask,database=geodb,trans=geotran,geometry="geometric",)
except:
    print "does ",gmask," already exist?"
    print "why are you asking me to make it again?"

if npass == 1:
    print "Now run dimsum through first pass to generate object masks."
    print "Copy objmask* to directory backup/"
    print "Repeat for ",prefix
    print "Copy backup/objmask* to ."
    print "Come back here when you are finished."
    print "Good Luck!"
if npass == 2:
    print "You did it!"
    #print "Now rerun disum on second pass."
    print "Now let's run align to create final image!"


    #iraf.align(refer=reference2, images = "@diminlist", outimage=finalname)

    print "Have fun, and come back soon!"
"""
** after first flattening, run dimsum
- the point is to just generate object masks.
- repeat for all objects.

"""


