#! /usr/bin/env python
from pylab import *
import sys, glob, os
# to run
# [rb227rfmbp151:irac/ch1/bcd] rfinn% pwd 
# /Users/rfinn/research/KISS/data/K0225/irac/ch1/bcd
# [rb227rfmbp151:irac/ch1/bcd] rfinn% reduceirac.py K0225 1

kgal=sys.argv[1]
channel=sys.argv[2]
print kgal,channel

def writefileswithpath(files,outf):
    outfile=open(outf,'w')
    path=os.getcwd()
    for file in files:
        newfile=path+'/'+file+'\n'
        outfile.write(newfile)
    outfile.close()

def initialize():
    os.system('mkdir cdf')
    os.system('mkdir cal')
    os.system('mkdir output')
    #os.system('cp /Applications/mopex/cal/prfmap.tbl output/.')
    os.system('cp /Users/rfinn/research/KISS/data/cdf/FIF.tbl output/.')
    os.system('mkdir FirstFrames')
    os.system('mv *0000_0000*.fits FirstFrames/.')
    #os.system('cp /Users/rfinn/clusters/spitzer/flatfield_24_ediscs.nl cdf/flatfield_24_ediscs.nl')
#    os.system('cp /Users/rfinn/clusters/spitzer/mosaic_24_ediscs.nl cdf/mosaic_24_ediscs.nl')
    os.system('ln -s ~/mopex_030106/cal cal')
#    os.system('mkdir firstframes')
    files=glob.glob('SPITZER*cbcd.fits')
    outf='imageList.txt'
    writefileswithpath(files,outf)


    files=glob.glob('SPITZER*cbunc.fits')
    outf='sigmaList.txt'
    writefileswithpath(files,outf)


    files=glob.glob('SPITZER*bimsk.fits')
    outf='maskList.txt'
    writefileswithpath(files,outf)

    #os.system('ls SPITZER*cbcd.fits  > InputImageList.txt')
    #os.system('ls SPITZER*cbunc.fits > SigmaList.txt')
    #os.system('ls SPITZER*bimsk.fits > DmaskList.txt')


    s="sed -e 's/K0847/"+kgal+"/g' /Applications/mopex/cdf/FinnIRACch"+channel+".nl > /Applications/mopex/cdf/FinnIRACch"+channel+kgal+".nl"
    os.system(s)
    print "Date of Observation = "
    os.system("gethead -f DATE_OBS ../pbcd/SPIT*maic.fits")
    print "Find pmask file with the first date AFTER the date of observation"
    print " "
    print "Loading pbcd mosaic in ds9 for your viewing pleasure"
    os.system("ds9 ../pbcd/SPIT*maic.fits -zscale -regions /Users/rfinn/research/KISS/RegionsFiles/kiss.reg&")

initialize()
