#!/usr/bin/env python

# Written by Rose A. Finn, July 1, 2010
#
# program reads in the agcnorth.sav file using idlsave
# and then cuts the agc according to RA, Dec, and redshift required
# for the Local Cluster Survey
#
# usage:
# import ReadAGCsav
# agc=ReadAGCsav.agc()
#
# to then access the agc data, for example:
# plot(agc.radeg,agc.decdeg,'k.',markersize=1)
# 
# if you don't want to cut the AGC based on RA, Dec, and recession vel,
# remove or comment the line self.cull() at the end of the __init__ routine.
#
# a full description of the AGC can be found at
# http://caborojo.astro.cornell.edu/alfalfalog/idldocs/agcinfo.php


import idlsave
from pylab import *
import atpy
import os
mypath=os.getcwd()
if mypath.find('Users') > -1:
    print "Running on Rose's mac pro or laptop"
    homedir='/Users/rfinn/'
elif mypath.find('home') > -1:
    print "Running on coma"
    homedir='/home/rfinn/'

class agc:
    def __init__(self,infile):
        a=idlsave.read(infile)
        if infile.find('agctotal') > -1:
            a.agc=a.agctotal
        elif infile.find('agcnorthminus1') > -1:
            a.agc=a.agcnorthminus1
        elif infile.find('agcnorth') > -1:
            a.agc=a.agcnorth

        self.outfile=infile.split('.sav')[0]+'.fits'
        self.agcnumber=a.agc.agcnumber[0]
        self.which=a.agc.which[0]
        self.rah=array(a.agc.rah[0],'f')
        self.ram=array(a.agc.ram[0],'f')
        self.ras10=array(a.agc.ras10[0],'f')
        self.radeg=zeros(len(self.ras10),'d')
        self.radeg=15.*(self.rah+self.ram/60.+self.ras10/10./3600.)
        self.decd=array(a.agc.decd[0],'f')
        self.decm=array(a.agc.decm[0],'f')
        self.decs=array(a.agc.decs[0],'f')
        self.decdeg=zeros(len(self.ras10),'d')
        self.decdeg=self.decd+self.decm/60.+self.decs/3600.
        self.a100=a.agc.a100[0]
        self.b100=a.agc.b100[0]
        self.mag10=a.agc.mag10[0]
        self.inccode=a.agc.inccode[0]
        self.posang=a.agc.posang[0]
        self.description=a.agc.description[0]
        self.bsteintype=a.agc.bsteintype[0]
        self.vopt=a.agc.vopt[0]
        self.verr=a.agc.verr[0]
        self.extrc3=a.agc.extrc3[0]
        self.extdirbe=a.agc.extdirbe[0]
        self.vsource=a.agc.vsource[0]
        self.ngcic=a.agc.ngcic[0]
        self.flux100=a.agc.flux100[0]
        self.rms100=a.agc.rms100[0]
        self.v21=a.agc.v21[0]
        self.width=a.agc.width[0]
        self.widtherr=a.agc.widtherr[0]
        self.telcode=a.agc.telcode[0]
        self.detcode=a.agc.detcode[0]
        self.hisource=a.agc.hisource[0]
        self.statuscode=a.agc.statuscode[0]
        self.snratio=a.agc.snratio[0]
        self.ibandqual=a.agc.ibandqual[0]
        self.ibandsrc=a.agc.ibandsrc[0]
        self.irasflag=a.agc.irasflag[0]
        self.icluster=a.agc.icluster[0]
        self.hidata=a.agc.hidata[0]
        self.iposition=a.agc.iposition[0]
        self.ipalomar=a.agc.ipalomar[0]
        self.rc3flag=a.agc.rc3flag[0]
        #self.cull()

    def cull(self,ramin=170,ramax=250,decmin=0,decmax=38.,vmin=0.018*3.e5,vmax=0.05*3.e5):

        posflag=(self.radeg>ramin)&(self.radeg < ramax)&(self.decdeg>decmin)&(self.decdeg<decmax)
        voptflag=(self.vopt > vmin)&(self.vopt<vmax)
        v21flag=(self.v21>vmin)&(self.v21<vmax)
        vflag=voptflag|v21flag
        keepflag=posflag&vflag
        self.agcnumber=self.agcnumber[keepflag]
        self.which=self.which[keepflag]
        self.radeg=self.radeg[keepflag]
        self.decdeg=self.decdeg[keepflag]
        self.a100=self.a100[keepflag]
        self.b100=self.b100[keepflag]
        self.mag10=self.mag10[keepflag]
        self.inccode=self.inccode[keepflag]
        self.posang=self.posang[keepflag]
        self.description=self.description[keepflag]
        self.bsteintype=self.bsteintype[keepflag]
        self.vopt=self.vopt[keepflag]
        self.verr=self.verr[keepflag]
        self.extrc3=self.extrc3[keepflag]
        self.extdirbe=self.extdirbe[keepflag]
        self.vsource=self.vsource[keepflag]
        self.ngcic=self.ngcic[keepflag]
        self.flux100=self.flux100[keepflag]
        self.rms100=self.rms100[keepflag]
        self.v21=self.v21[keepflag]
        self.width=self.width[keepflag]
        self.widtherr=self.widtherr[keepflag]
        self.telcode=self.telcode[keepflag]
        self.detcode=self.detcode[keepflag]
        self.hisource=self.hisource[keepflag]
        self.statuscode=self.statuscode[keepflag]
        self.snratio=self.snratio[keepflag]
        self.ibandqual=self.ibandqual[keepflag]
        self.ibandsrc=self.ibandsrc[keepflag]
        self.irasflag=self.irasflag[keepflag]
        self.icluster=self.icluster[keepflag]
        self.hidata=self.hidata[keepflag]
        self.iposition=self.iposition[keepflag]
        self.ipalomar=self.ipalomar[keepflag]
        self.rc3flag=self.rc3flag[keepflag]
        self.agcdict=dict((a,b) for a,b in zip(self.agcnumber,arange(len(self.agcnumber))))
    
    def writefitsfile(self,outfile=None):
        lcstab=atpy.Table()
        lcstab.add_column('AGCNUMBER',self.agcnumber)
        lcstab.add_column('WHICH',self.which)
        lcstab.add_column('RA',self.radeg,unit='deg')
        lcstab.add_column('DEC',self.decdeg,unit='deg')
        lcstab.add_column('A100',self.a100,unit='arcsec')
        lcstab.add_column('B100',self.b100,unit='arcsec')
        lcstab.add_column('MAG10',self.mag10,unit='mag')
        lcstab.add_column('INCCODE',self.inccode,unit='')
        lcstab.add_column('POSANG',self.posang,unit='')
        lcstab.add_column('DESCRIPTION',self.description,unit='')
        lcstab.add_column('BSTEINTYPE',self.bsteintype,unit='')
        lcstab.add_column('VOPT',self.vopt,unit='km/s')
        lcstab.add_column('VERR',self.verr,unit='km/s')
        lcstab.add_column('EXTRC3',self.extrc3,unit='')
        lcstab.add_column('EXTDIRBE',self.extdirbe,unit='')
        lcstab.add_column('VSOURCE',self.vsource,unit='')
        lcstab.add_column('NGCIC',self.ngcic,unit='')
        lcstab.add_column('FLUX100',self.flux100,unit='')
        lcstab.add_column('RMS100',self.rms100,unit='')
        lcstab.add_column('V21',self.v21,unit='')
        lcstab.add_column('WIDTH',self.width,unit='')
        lcstab.add_column('WIDTHERR',self.widtherr,unit='')
        lcstab.add_column('TELCODE',self.telcode,unit='')
        lcstab.add_column('DETCODE',self.detcode,unit='')
        lcstab.add_column('HISOURCE',self.hisource,unit='')
        lcstab.add_column('STATUSCODE',self.statuscode,unit='')
        lcstab.add_column('SNRATIO',self.snratio,unit='')
        lcstab.add_column('IBANDQUAL',self.ibandqual,unit='')
        lcstab.add_column('IBANDSRC',self.ibandsrc,unit='')
        lcstab.add_column('IRASFLAG',self.irasflag,unit='')
        lcstab.add_column('ICLUSTER',self.icluster,unit='')
        lcstab.add_column('HIDATA',self.hidata,unit='')
        lcstab.add_column('IPOSITION',self.iposition,unit='')
        lcstab.add_column('IPALOMAR',self.ipalomar,unit='')
        lcstab.add_column('RC3FLAG',self.rc3flag,unit='')
        #lcstab.write(homedir+'idl/programs/idl_alfa/agcnorthLCS.fits')
        #lcstab.write(homedir+'idl/programs/idl_alfa/'+outfile+'.fits')
        if outfile == None:
            outfile = self.outfile
        if os.path.exists(outfile):
            os.remove(outfile)
        lcstab.write(outfile)
if __name__ == '__main__':
    #infile=homedir+'idl/programs/idl_alfa/agcnorthminus1.sav'
    #infile=homedir+'idl/programs/idl_alfa/agctotal.sav'
    #infiles=[homedir+'research/AGC/agctotal.sav',homedir+'research/AGC/agcnorth.sav',homedir+'research/AGC/agcnorthminus1.sav']
    infiles=[homedir+'research/AGC/agctotal.sav',homedir+'research/AGC/agcnorthminus1.sav']
    #infiles=[homedir+'research/AGC/agcnorthminus1.sav']
    for file in infiles:
        a=agc(file)
        a.writefitsfile()
    a.cull()
    a.writefitsfile(outfile=homedir+'research/AGC/agcnorthLCSregion.fits')
