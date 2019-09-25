#!/usr/bin/env python
import glob
import os

os.system('mkdir ../firstframes')
os.system('mv SPITZER*_000*_0000*.fits ../firstframes')
currentDir = os.getcwd()
def writefile(matchstring,output):
    #bcdfiles=glob.glob('SPITZER*bcd.fits')
    bcdfiles=glob.glob(matchstring)
    outfile=open(output,'w')
    for file in bcdfiles:
	fullpath=currentDir+'/'+file+'\n'
	outfile.write(fullpath)
    outfile.close()

writefile('SPITZER*bcd.fits','InputImageList.txt')
writefile('SPITZER*bunc.fits','SigmaList.txt')
writefile('SPITZER*bbmsk.fits','DmaskList.txt')
os.system('ln -s /home/rfinn/mopex/cal cal')
os.system('mkdir cdf')
os.system('cp ~/mopex/mycdf/flatfield_24_LC.nl cdf/')
#os.system('echo COPY flatfield_24_LC.nl to cdf directory!!!')
