#!/usr/bin/env python
'''
writing code to do a quick analysis of the depth of the 09/2009 data

RXJ1821 looks pretty good

RXJ1716 is compromised by two bright stars on the field

'''
from z08common import *
from pylab import *
import os
import asciidata
import atpy

class baseCluster:
    def __init__(self,clustername):
        self.radeg=clusterRAdeg[clustername]
        self.decdeg=clusterDecdeg[clustername]
        #self.readSextractor('testJ.cat')
        #self.readSextractorNB('testNB.cat')
        self.prefix=clustername
        #self.getNewfirmCenters()
    def readSextractor(self,filename):
        sexJ=atpy.Table(filename,type='ascii')
        self.x=sexJ.col2
        self.y=sexJ.col3
        self.fluxbest=sexJ.col22
        self.fluxbesterr=sexJ.col23
        self.magbest=sexJ.col28
        self.magbesterr=sexJ.col29
        self.fluxauto=sexJ.col22
        self.fluxautoerr=sexJ.col23
        self.magauto=sexJ.col24
        self.magautoerr=sexJ.col25
        self.fluxiso=sexJ.col10
        self.fluxisoerr=sexJ.col11
        self.magiso=sexJ.col12
        self.magisoerr=sexJ.col13
        self.classstar=sexJ.col95
        self.snr=abs(self.fluxbest/self.fluxbesterr)
        self.snrauto=abs(self.fluxauto/self.fluxautoerr)
    def readSextractorNB(self,filename):
        sexJ=atpy.Table(filename,type='ascii')
        self.NBx=sexJ.col2
        self.NBy=sexJ.col3
        self.NBfluxbest=sexJ.col22
        self.NBfluxbesterr=sexJ.col23
        self.NBmagbest=sexJ.col28
        self.NBmagbesterr=sexJ.col29
        self.NBfluxauto=sexJ.col22
        self.NBfluxautoerr=sexJ.col23
        self.NBmagauto=sexJ.col24
        self.NBmagautoerr=sexJ.col25
        self.NBfluxiso=sexJ.col10
        self.NBfluxisoerr=sexJ.col11
        self.NBmagiso=sexJ.col12
        self.NBmagisoerr=sexJ.col13
        self.NBclassstar=sexJ.col95
        self.NBsnr=abs(self.NBfluxbest/self.NBfluxbesterr)
        self.NBsnrauto=abs(self.NBfluxauto/self.NBfluxautoerr)

        
        #rxj18=cluster('RXJ1821')
    def getNewfirmCenters(self):
        newfirm_offset_ra_deg=newfirm_offset[self.prefix][0]/3600.
        newfirm_offset_dec_deg=newfirm_offset[self.prefix][1]/3600.
        self.newfirm_ra=self.radeg+newfirm_offset_ra_deg
        self.newfirm_dec=self.decdeg+newfirm_offset_dec_deg
        # print in hr:mm:ss
        rah=(self.newfirm_ra/15.)
        hr=floor(rah)
        min=floor((rah-floor(rah))*60.)
        sec=(rah-(hr+min/60.))*3600.
        
        ddeg=floor(self.newfirm_dec)
        dmin=floor((self.newfirm_dec-ddeg)*60.)
        dsec=(self.newfirm_dec-(ddeg+dmin/60.))*3600.

        print 'RA = ',hr,min,sec,' Dec=',ddeg,dmin,dsec

    def makeHoldenRegionFile(self):
        if self.prefix.find('RDCSJ1317') > -1:
            infile=homedir+'research/z08clusters/AncillaryData/Holden/r1317/c1317+29.skycat.v2'
        in1=asciidata.open(infile)
        self.hdata=in1
        ncol=len(in1)
        nrow=len(in1[0])
        redshift=array(in1[11],'f')
        zindex=where(redshift > -1)
        out1=homedir+'research/z08clusters/RegionsFiles/'+self.prefix+'_holdenRADec.reg'
        outfile=open(out1,'w')
        outfile.write('global color=green width=2\n')
        outfile.write('fk5 \n')
        zi=zindex[0]
        for i in zi:
            print i,zindex,redshift[i], float(in1[11][i])
            if (redshift[i] > zmin) & (redshift[i] < zmax):
                ccolor='blue'
            else:
                ccolor='cyan'
            s='circle(%12.8f,%12.8f,5\") # color= %s text={%s}\n'%(float(in1[1][i]),float(in1[2][i]),ccolor,in1[11][i])
            outfile.write(s)
        outfile.close()
        in1.close()
        
