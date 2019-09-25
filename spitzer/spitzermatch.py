#!/usr/bin/env python
import sys, os
import Numeric as N
#import scipy
from math import *
import mystuff as my
import ppgplot
import random
import sets
import pylab
from pyraf import iraf

deflw=3
ewlabel = "EW H\ga+[NII] (\(2078))"
ewo2label = "EW [OII] (\(2078))"
sfrlabel = "SFR (h\d100\u\u-2\d M\d\(2281)\u yr\u-1\d)"
sfrlabel70 = "SFR (M\d\(2281)\u yr\u-1\d)"#label for ms1054, h100=.7
sumsfrlabel = "\gSSFR (h\d100\u\u-2\d M\d\(2281)\u yr\u-1\d)"
sfrmlabel = "\gSSFR/M\dcl\u (h\d100\u\u-3\d M\d\(2281)\u yr\u-1\d/10\u14\d M\d\(2281)\u)"
mcllabel = "M\dcl\u (h\d100\u\u-1\d 10\u14\d M\d\(2281)\u)"
omega0=0.3
omegaL=0.7
h100=.7
dlambdaJ=2500 #bandwidth of J filter
nstar=0
pscale=0.18
sfrmin=.5
sfrcomp=.5#sfr cut for comparision w/literature
#parameters for comparison w/literature
f200a=.5 #fraction of R200 to include in integrated SFR and femiss
f200b=2. #fraction of R200 to include in integrated SFR and femiss
lhamin=sfrmin*1.27#lhamin in units of 1d40 ergs/s
snmin=3. #min s/n cut for continuum-subtracted flux
ewmin=10.
ewcomp=20#ew cut for comparison w/low-z studies
dm=0.
ewburst=40 #ew cut for defining starburst
sfrcorr=2.5*0.77#sfr correction factors for dust and NII


def psplotinit(output):
    file=output+"/vcps"
    ppgplot.pgbeg(file,1,1)
    ppgplot.pgpap(8.,1.)
    ppgplot.pgsch(1.7) #font size
    ppgplot.pgslw(6)  #line width

def plothistsfr():
    DATAMIN = -4.
    DATAMAX = 15.
    NBIN = int((DATAMAX-DATAMIN)*2.)    
    #print "ngal = ",len(g0.sfr)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ppgplot.pgenv(DATAMIN,DATAMAX,0,45,0)
    ppgplot.pglab("SFR (h\d100\u\u-2\d M\d\(2281)\u yr\u-1 \d)","Number of Galaxies","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    #ppgplot.pgsci(4)
    #x=N.compress((abs(g0.ew) > ewmin),g0.sfr)
    x=N.compress((g0.final > 0),g0.sfrc)
    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    xlabel = 6.5
    ylabel = 38.
    ystep = 3.
    dy=.4
    dxl=3
    dxr=.5
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1040")
    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgline(xlin,ylin)

    ppgplot.pgslw(5)
    ppgplot.pgsls(3)#dot-dash-dot-dash
    #ppgplot.pgsci(3)
    #x=N.compress((abs(g1.ew) > ewmin),g1.sfr)
    x=N.compress((g1.final > 0),g1.sfrc)
    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)

    ylabel = ylabel - ystep

    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgline(xlin,ylin)
    ppgplot.pgsls(1)
    ppgplot.pgslw(deflw)
    ppgplot.pgtext(xlabel,ylabel,"CL1054-12")
    
    ppgplot.pgsls(1)#dot-dash-dot-dash
    #ppgplot.pgsci(2)
    ppgplot.pgslw(2)  #line width
    #x=N.compress((abs(g2.ew) > ewmin),g2.sfr)
    x=N.compress((g2.final > 0),g2.sfrc)
    ppgplot.pghist(len(x),x,DATAMIN,DATAMAX,NBIN,5)
    ylabel = ylabel - ystep
    ppgplot.pgslw(deflw)  #line width
    ppgplot.pgtext(xlabel,ylabel,"CL1216")
    xlin = N.array([xlabel-dxl,xlabel-dxr],'f')
    ylin = N.array([ylabel+dy,ylabel+dy],'f')
    ppgplot.pgslw(2)  #line width
    ppgplot.pgline(xlin,ylin)


    #print "Number in g2.ratios = ",len(g2.ratio)
    #for ratio in g2.ratio:
    #    print ratio
    #drawbinned(x,y,5)
    ppgplot.pgsci(1)
    #ppgplot.pgend()


def plotgal(xg,yg,final,finalsf,ncl):
    ppgplot.pgslw(6)
    ppgplot.pgsls(1)

    xmin = (-1.*c.pscale[ncl]*(c.xc[ncl]-1))
    xmax = (c.pscale[ncl]*(c.xmax[ncl]-c.xc[ncl]))    
    ymin = (-1.*c.pscale[ncl]*(c.yc[ncl]-1))
    ymax = (c.pscale[ncl]*(c.ymax[ncl]-c.yc[ncl]))
    ppgplot.pgbox("",0.0,0,"L",0.0,0)

    dx=5.
    ppgplot.pgenv(xmin-dx,xmax+dx,ymin-dx,ymax+dx,0)
    ppgplot.pglab("\gD DEC (\")","\gD RA (\")","")
    ppgplot.pgtext(-4,-4,"X")
    r = (0.5*c.r200pix[ncl]*c.pscale[ncl])
    ppgplot.pgslw(1)
    ppgplot.pgsls(2)
    ppgplot.pgsfs(2)
    ppgplot.pgcirc(0,0,r)
    #print "cluster ",ncl," r200: ",c.r200pix[ncl],c.r200Mpc[ncl], " Mpc"

    ppgplot.pgslw(3)
    ppgplot.pgsls(1)

    x = (xg - c.xc[ncl])*c.pscale[ncl]
    y = (yg - c.yc[ncl])*c.pscale[ncl]
    x = (N.compress((final > 0) & (finalsf < 1), xg) - c.xc[ncl])*c.pscale[ncl]
    y = (N.compress((final > 0) & (finalsf < 1), yg) - c.yc[ncl])*c.pscale[ncl]
    ppgplot.pgpt(x,y,22)
    x = (N.compress((final > 0) & (finalsf > 0), xg) - c.xc[ncl])*c.pscale[ncl]
    y = (N.compress((final > 0) & (finalsf > 0), yg) - c.yc[ncl])*c.pscale[ncl]
    ppgplot.pgpt(x,y,18)


def getnearestgen(self,ra1,dec1,ra2,dec2,j):#measure distances from ra1, dec1 to members in catalog ra2, dec2
    sig5=N.zeros(len(ra1),'f')
    sig10=N.zeros(len(ra1),'f')
    for i in range(len(ra1)):
        angdist=my.DA(c.z[j],h100)#kpc/arcsec
        dspec=N.sqrt((ra1[i]-ra2)**2+(dec1[i]-dec2)**2)#sorted array of distances in degrees
        dspecsort=N.take(dspec,N.argsort(dspec))
        sig5[i]=5./(N.pi)/(dspecsort[5]*3600.*angdist/1000.)**2#convert from deg to arcsec, multiply by DA (kpc/arcsec), divide by 1000 to convert to Mpc, index 5 element b/c 0 is itself
        sig10[i]=10./(N.pi)/(dspecsort[10]*3600.*angdist/1000.)**2#convert from deg to arcsec, multiply by DA (kpc/arcsec), divide by 1000 to convert to Mpc, index 5 element b/c 0 is itself
    return sig5, sig10


def findnearest(self,ra1,dec1,ra2,dec2,delta):#measure distances from ra1, dec1 to members in catalog ra2, dec2
    ramatch=my.findmatch(ra1,ra2,delta)
    decmatch=my.findmatch(dec1,dec2,delta)
    match=ramatch&decmatch
    for i in range(len(ra1)):
        angdist=my.DA(c.z[j],h100)#kpc/arcsec
        dspec=N.sqrt((ra1[i]-ra2)**2+(dec1[i]-dec2)**2)#sorted array of distances in degrees
        dspecsort=N.take(dspec,N.argsort(dspec))
        sig5[i]=5./(N.pi)/(dspecsort[5]*3600.*angdist/1000.)**2#convert from deg to arcsec, multiply by DA (kpc/arcsec), divide by 1000 to convert to Mpc, index 5 element b/c 0 is itself
        sig10[i]=10./(N.pi)/(dspecsort[10]*3600.*angdist/1000.)**2#convert from deg to arcsec, multiply by DA (kpc/arcsec), divide by 1000 to convert to Mpc, index 5 element b/c 0 is itself
    return sig5, sig10



class Galaxy:
    def __init__(self):#individual galaxy properties
        print "dude - a galaxy!"
        
    def readhafile(self,file):
        input=open(file,'r')
        #get number of galaxies
        ngal=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            ngal=ngal+1
        input.close()

        self.clusterid = N.zeros(ngal,'f')
        self.clusterr200 = N.zeros(ngal,'f')
        self.clustersig  = N.zeros(ngal,'f')
        self.clusterz = N.zeros(ngal,'f')
        self.clusterx = N.zeros(ngal,'f')
        self.clustery = N.zeros(ngal,'f')

        self.shortname = []
        self.dra = N.zeros(ngal,'f')
        self.ddec  = N.zeros(ngal,'f')
        self.fn = N.zeros(ngal,'f')
        self.errfn = N.zeros(ngal,'f')
        self.fj = N.zeros(ngal,'f')
        self.errfj = N.zeros(ngal,'f')
        self.magj = N.zeros(ngal,'f')
        self.errmagj = N.zeros(ngal,'f')
        self.contsub = N.zeros(ngal,'f')
        self.errcontsub = N.zeros(ngal,'f')
        self.snc = N.zeros(ngal,'f')

        self.ew = N.zeros(ngal,'f')
        self.errew = N.zeros(ngal,'f')
        
        self.sfr = N.zeros(ngal,'f')
        self.errsfr = N.zeros(ngal,'f')
        self.sf = N.zeros(ngal,'f')

        self.ra = N.zeros(ngal,'f')#ra in deg
        self.dec  = N.zeros(ngal,'f')#dec in deg

        input=open(file,'r')
        for i in range(ngal):
            if line.find('#') > -1: #skip lines with '#' in them
                continue

            t=line.split()
            (self.dra[i],self.ddec[i],self.fn[i],self.errfn[i],self.fj[i],self.errfj[i],self.magj[i],self.errmagj[i],self.contsub[i],self.errcontsub[i],self.snc[i],self.ew[i],self.errew[i],self.sfr[i],self.errsfr[i],self.sf[i])=(float(t[1]),float(t[2]),float(t[3]),float(t[4]),float(t[5]),float(t[6]),float(t[7]),float(t[8]),float(t[9]),float(t[10]),float(t[11]),float(t[12]),float(t[13]),float(t[14]),float(t[15]),float(t[16]))
            self.shortname.append(t[0])
            shortname=t[0]
            self.ra[i]=(float(shortname[0:2])+float(shortname[2:4])/60.+float(shortname[4:7])/10./3660.)*15.#convert from hours to degrees
            self.dec[i]=float(shortname[7:10])+float(shortname[10:12])/60.+float(shortname[12:])/10./3660.
    

class Spitzer24:
    def __init__(self):#individual galaxy properties
        print "dude - 24 micron data!"

    def readfile(self,file):
        input=open(file,'r')
        #get number of galaxies
        ngal=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue
            ngal=ngal+1
        input.close()


        self.id = N.zeros(ngal,'f')
        self.imagex = N.zeros(ngal,'f')
        self.imagey  = N.zeros(ngal,'f')
        self.ra = N.zeros(ngal,'f')
        self.dec = N.zeros(ngal,'f')
        self.f = N.zeros(ngal,'f')
        self.errf = N.zeros(ngal,'f')

        input=open(file,'r')
        i=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue

            t=line.split()
            (self.id[i],self.imagex[i],self.imagey[i],self.ra[i],self.dec[i],self.f[i],self.errf[i])=(float(t[0]),float(t[1]),float(t[2]),float(t[3]),float(t[4]),float(t[5]),float(t[6]))
            i=i+1
        input.close()
        outfile=open('xy24.dat','w')
        for i in range(len(self.imagex)):
            #print self.imagex[i],self.imagey[i],self.f[i],self.errf[i]
            string="%8.2f %8.2f %8.1f %8.2f \n" % (self.imagex[i],self.imagey[i],self.f[i],self.errf[i])
            outfile.write(string)
        outfile.close()


    def matchha(self):#match coordinates with ha data
        self.haflag=N.zeros(len(self.ra),'f')#flag for halpha
        self.hamatch=N.zeros(len(self.ra),'i')
        (x1sort,x1index)=my.sortwindex(gha.ra)
        (x2sort,x2index)=my.sortwindex(gha.dec)
        deltar=100.#matching tolerance in arcsec
        deltar=deltar/3600.#convert to degrees
        for i in range(len(self.ra)):
            temp=[]
            (temp,flag)=my.findmatch(self.ra[i],x1sort,deltar)#returns indices of sorted array
            if flag < 1:
                print i, "No match to Halpha data"
                continue
            templist1=temp
            for j in range(len(temp)):
                templist1[j]=x1index[temp[j]]
            
            temp=[]
            (temp,flag)=my.findmatch(self.dec[i],x2sort,deltar)
            if flag < 1:
                print i, "No match to Halpha data"
                continue

            templist2=temp
            for j in range(len(temp)):
                templist2[j]=x2index[temp[j]]
            a=sets.Set(templist1)
            b=sets.Set(templist2)
            members=a&b
            members=list(members)#contains all galaxies w/in a 4"x4" square
            min=100.
            for j in members:
                d=N.sqrt((self.ra[i]-gha.ra[j])**2+(self.dec[i]-gha.dec[j])**2)
                if d < min:
                    min=d
                    minindex=j
            if min < 2.:
                self.haflag[i]=1.
                self.hamatch[i]=minindex
                
            


gha = Galaxy()
gha.readhafile('cl1216sfrtableapj.dat')

g24 = Spitzer24()
g24.readfile('mosaic_extract.tbl')
g24.matchha()

iraf.images.tv.display('../mosaic.fits',1,zscale='no',zrange='no',z1=33.,z2=40.)
iraf.images.tv.tvmark(1,'xy24.dat',mark='circle', radii=3,color=204)
#pylab.plot(g24.f, g24.errf,'bo')
#pylab.show()
#pylab.savefig("test.ps")

