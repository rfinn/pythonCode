#!/usr/bin/env python
"""useage
sdssinter.py mode completenesscorr, where mode
= 0 for updating all
= 1 for updating plots only

and completenesscorr
=0 for no completeness correction for fiber sampling
=1 for yes
"""
import sys, glob
import numarray as N
#import scipy
import pylab #use matplotlib instead of scipy
from math import *
import ppgplot
import time, os
import mystuff as my
import matplotlib
matplotlib.use('Agg')

from matplotlib import rc
rc('font',family='serif', style='normal', variant='normal',weight='bold', stretch='normal', size='small')
import sqlcl

#import matplotlib.matlab as mat
#import ediscssfr
#from sdssplots import *
starttime=time.clock()
print "start time = ",starttime

#mode=int(sys.argv[1])#
#print "mode = ",mode
omega0=0.3
omegaL=0.7
zminp=.04#min z for plots
zmaxp=.1#max z for plots
h100=.7
H0=100.*h100
c=3.e5
mabscut=-19.8#M_V cut
#mabscut=-16.

def psplotinit(output):
    file=output+"/vcps"
    ppgplot.pgbeg(file,1,1)
    ppgplot.pgpap(8.,1.)
    ppgplot.pgsch(1.7) #font size
    ppgplot.pgslw(7)  #line width


class Cluster:
    def __init__(self):
        #columns of Bianca's file: clusters.spec.dat
        self.id = []#Abell Number
        self.zSR = []#redshift from Struble & Rood 1991
        self.sigmav = []#velocity dispersion from Struble & Rood 1991
        self.Nz = [] #number of galaxies used to determine z
        self.raSR = [] #RA from Struble & Rood
        self.decSR  = [] #Dec from Struble & Rood
        self.counts = []#number of cluster members b/w m3 and m3+2
        self.richness = []#richness class
        self.distance = []#distance class
        self.m10 = []#m10 class?
        self.A1 = []#2.0 Abell radii
        self.ra = []#ra manually from sky server
        self.dec = []#dec manually from sky server
        self.z = []#from SDSS galaxies
        self.sigmaz=[]#from SDSS galaxies
        self.sigmazerr=[]
        self.sigma = []#sigma_z/(1+z)*300000 km/s
        self.sigmaerr = []#error from bootstrapping
        self.r200 = []#in Mpc, Finn 2005
        self.r200err=[]
        self.r200deg = []#R200 in degrees
        self.r200degerr=[]
        self.Ngal = []#number of galaxies
    def creadfiles(self,clusterfile):
        for line in open(clusterfile):
            if line.find('#') > -1:
                continue
            fields=line.split()
            sum=0
            self.id.append(int(fields[0]))
            self.zSR.append(float(fields[1]))
            self.sigmav.append(float(fields[2]))
            self.Nz.append(float(fields[3]))
            self.raSR.append(float(fields[4]))
            self.decSR.append(float(fields[5]))
            self.counts.append(float(fields[6]))
            self.richness.append(float(fields[7]))
            self.distance.append(float(fields[8])) #distance class
            self.m10.append(float(fields[9]))#m10 class?
            self.A1.append(float(fields[10])) #2.0 Abell radii
            self.ra.append(float(fields[11])) #ra manually from sky server
            self.dec.append(float(fields[12])) #dec manually from sky server
            self.z.append(float(fields[13])) #from SDSS galaxies
            self.sigmaz.append(float(fields[14])) #from SDSS galaxies
            self.sigma.append(float(fields[16])) #sigma_z/(1+z)*300000 km/s
            self.r200.append(float(fields[18])) #in Mpc, Finn 2005
            self.r200deg.append(float(fields[20])) #R200 in degrees
            self.Ngal.append(float(fields[22]))#number of galaxies
        print "got ",len(self.z)," clusters!"
    def convarray(self):
        self.id = N.array(self.id,'i')#Abell Number
        self.zSR = N.array(self.zSR,'f')#redshift from Struble & Rood 1991
        self.sigmav = N.array(self.sigmav,'f')#redshift from Struble & Rood 1991
        self.Nz = N.array(self.Nz,'f')#number of galaxies used to determine z
        self.raSR = N.array(self.raSR,'f')#RA from Struble & Rood
        self.decSR  = N.array(self.decSR,'f')#Dec from Struble & Rood
        self.counts = N.array(self.counts,'f')#number of cluster members b/w m3 and m3+2
        self.richness =N.array(self.richness,'f')#richness class
        self.distance = N.array(self.distance,'f')#distance class
        self.m10 = N.array(self.m10,'f')#m10 class?
        self.A1 = N.array(self.A1,'f')#2.0 Abell radii
        self.ra = N.array(self.ra,'f')#ra manually from sky server
        self.dec = N.array(self.dec,'f')#dec manually from sky server
        self.z = N.array(self.z,'f')#from SDSS galaxies
        self.sigmaz=N.array(self.sigmaz,'f')#from SDSS galaxies
        self.sigma = N.array(self.sigma,'f')#sigma_z/(1+z)*300000 km/s
        self.r200 = N.array(self.r200,'f')#in Mpc, Finn 2005
        self.r200deg = N.array(self.r200deg,'f')#R200 in degrees
        self.Ngal = N.array(self.Ngal,'f')#number of galaxies
    def Kcorr(self):#calculate K_u(z) (Blanton et al 2003) for each cluster
        self.kcorr=N.zeros(len(self.z),'f')
        self.dL=N.zeros(len(self.z),'f')
        z=N.arange(0.,1.2,.1)
        #kc=N.array([0.,0.,.05,.05,.1,.15,.2,.225,.25,.3,.35,.35],'f')#Ku(z)
        kc=N.array([0.,0.,.025,.05,.07,.1,.14,.16,.2,.25,.25,.3],'f')#Kg(z)
        r=self.z/.01
        for i in range(len(self.z)):
            self.kcorr[i]=kc[int(r[i])]+(r[i]-int(r[i]))*(kc[int(r[i]+1)]-kc[int(r[i])])
            self.dL[i] = my.dL(self.z[i],h100)
            
    def getsdssphotcats(self):  #get photometric sources within 2R200 
        print "elapsed time = ",time.clock()-starttime
        self.mcut=N.zeros(len(self.z),'f')
        cl=N.arange(17,len(self.z),1)
        for i in range(len(self.z)):
        #for i in cl:
            dL = self.dL[i]
            print "getting phot cat for cluster abell",self.id[i]
            r200arcmin=self.r200deg[i]*60.
            #drsearch=2.*r200arcmin#2xR200 in arcmin for sdss query
            drsearch=3.*r200arcmin#2xR200 in arcmin for sdss query
            #Vg=0.3556-0.7614*((self.avegr)-0.6148)#(V-g) from Blanton et al 2003
            mr=mabscut - 0.1331 + 5.*N.log10(dL)+25.+self.kcorr[i]
            print i, self.z[i], dL, mr
            self.mcut[i]=mr
            print "ra, dec, dr, mr = %12.8f %12.8f %8.3f %5.2f" % (self.ra[i],self.dec[i],drsearch,mr)
            #query="select g.ra, g.dec, g.u, g.g, g.r, g.i, g.z, g.plate_ID, g.MJD,  from galaxy g, dbo.fGetNearbyObjEq(%12.8f,%12.8f,%8.3f) n where g.objID = n.objID and (g.g < %5.2f) and ((0.384*g.g + 0.716*g.r)< %5.2f)" % (self.ra[i],self.dec[i],drsearch,(mr+1.5),mr)
            query="select g.ra, g.dec, g.u, g.g, g.r, g.i, g.z, g.objid, g.specObjID,g.extinction_u, g.extinction_g, g.extinction_r, g.extinction_i, g.extinction_z from galaxy g, dbo.fGetNearbyObjEq(%12.8f,%12.8f,%8.3f) n where g.objID = n.objID and (g.g < %5.2f) and  (g.PrimTarget & 0x00000040) > 0 " % (self.ra[i],self.dec[i],drsearch,(mr))
	    #line from sdssinter.py code
	    #query="select n.distance from galaxy g, dbo.fGetNearbyObjEq(%12.8f,%12.8f,%8.3f) n where g.objID = n.objID and g.r< %5.2f and  (g.PrimTarget & 0x00000040) > 0 order by distance" % (self.ra[i],self.dec[i],drsearch,mr)#added flags to get rid of saturated objects, stars, etc


            #query="select g.ra, g.dec, g.u, g.g, g.r, g.i, g.z from galaxy g, dbo.fGetNearbyObjEq(%12.8f,%12.8f,%8.3f) n where g.objID = n.objID" % (self.ra[i],self.dec[i],drsearch)#no mag cut

            #query="select g.ra, g.dec, g.u, g.g, g.r, g.i, g.z from galaxy g, dbo.fGetNearbyObjEq(%12.8f,%12.8f,%8.3f) n where g.objID = n.objID and g.r<17.7 and g.g<18.0" % (self.ra[i],self.dec[i],drsearch)
            lines=sqlcl.query(query).readlines()
            #print query
            print "got number+1 phot objects = ",len(lines)
            #print lines
            output="abell"+str(self.id[i])+".phot.dat"
            outfile=open(output,'w')
            outfile.write("#%s " % (lines[0]))
            for line in lines[1:]:
                outfile.write("%s " % (line))
            outfile.close()
            #calculate Mv
    def getsdssspeccats(self):  #get photometric sources within 2R200 
        print "elapsed time = ",time.clock()-starttime
        self.mcut=N.zeros(len(self.z),'f')
        for i in range(len(self.z)):
            dL = self.dL[i]
            print "getting spec cat for cluster abell",self.id[i]
            r200arcmin=self.r200deg[i]*60.
            drsearch=3.*r200arcmin#2xR200 in arcmin for sdss query
            #Vg=0.3556-0.7614*((self.avegr)-0.6148)#(V-g) from Blanton et al 2003
            mr=mabscut - 0.1331 + 5.*N.log10(dL)+25.+self.kcorr[i]
            print i, self.z[i], dL, mr
            self.mcut[i]=mr
            dz=3*self.sigma[i]/(3.e5)*(1+self.z[i])
            zmax=self.z[i]+.5*dz
            zmin=self.z[i]-.5*dz
            print "ra, dec, dr, mr = %12.8f %12.8f %8.3f %5.2f" % (self.ra[i],self.dec[i],drsearch,mr)
            query="select g.ra, g.dec, g.u, g.g, g.r, g.i, g.z, g.objid, g.specObjID,g.extinction_u, g.extinction_g, g.extinction_r, g.extinction_i, g.extinction_z, l.ew, l.ewErr, l2.ew, l2.ewErr from galaxy g, specobj s, SpecLine l, SpecLine l2, dbo.fGetNearbyObjEq(%12.8f,%12.8f,%8.3f) n where g.objID = n.objID and g.objID = s.bestobjid and s.specobjID=l.specobjID and s.specobjID=l2.specobjID and (g.g < %5.2f) and  (g.PrimTarget & 0x00000040) > 0 and (s.z > %6.4f) and (s.z < %6.4f) and l.LineID = 3727 and l2.LineID = 6565" % (self.ra[i],self.dec[i],drsearch,(mr+3),zmin,zmax)
#	    query="select g.ra, g.dec, g.u, g.g, g.r, g.i, g.z, g.objid, g.specObjID,g.extinction_u, g.extinction_g, g.extinction_r, g.extinction_i, g.extinction_z, l.ew, l.ewErr, l2.ew, l2.ewErr from galaxy g, specobj s, SpecLine l, SpecLine l2, dbo.fGetNearbyObjEq(%12.8f,%12.8f,%8.3f) n where g.objID = n.objID and g.objID = s.bestobjid and s.specobjID=l.specobjID and s.specobjID=l2.specobjID and (g.g < %5.2f) and  (g.PrimTarget & 0x00000040) > 0 and l.LineID = 3727 and l2.LineID = 6565" % (self.ra[i],self.dec[i],drsearch,(mr))


            

            lines=sqlcl.query(query).readlines()
            #print query
            print "got number+1 spec objects w/in 2R200= ",len(lines)
            #print lines
            output="abell"+str(self.id[i])+".spec2r200.dat"
            outfile=open(output,'w')
            outfile.write("#%s " % (lines[0]))
            for line in lines[1:]:
                outfile.write("%s " % (line))
            outfile.close()
            #calculate Mv

                
    def readsdsscompleteness(self):
        self.nspec05=N.zeros(len(self.z),'f')
        self.nphot05=N.zeros(len(self.z),'f')
        self.compl05=N.zeros(len(self.z),'f')
        self.nspec1= N.zeros(len(self.z),'f')
        self.nphot1= N.zeros(len(self.z),'f')
        self.compl1= N.zeros(len(self.z),'f')
        self.nspec2= N.zeros(len(self.z),'f')
        self.nphot2= N.zeros(len(self.z),'f')
        self.compl2= N.zeros(len(self.z),'f')
        complout=open('sdsscompleteness.dat','r')
        i=0
        for line in complout:
            f=line.split()
            #print line
            #print f
            j=0
            for j in range(len(f)):
                f[j]=float(f[j])
            print f
            (self.nspec05[i],self.nphot05[i],self.compl05[i],self.nspec1[i],self.nphot1[i],self.compl1[i],self.nspec2[i],self.nphot2[i],self.compl2[i])=(f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8])
            #print f[0],f[1],f[2],self.nspec05[i],self.nphot05[i],self.compl05[i]
            #print line
            i=i+1
        complout.close()
        print "got ",i," clusters from sdsscompleteness.dat"
        print "min of c.compl05 = ",min(self.compl05)
        #print self.compl05
    def assigncolor(self,gr):#assign average g-r color to each cluster based on observed colors of spec members
        self.avegr=gr
        print "make sure these are equal:"
        print "length of clusters = ",len(self.z)
        print "length of gr = ",len(gr)

class Galaxy:
    def __init__(self):
        self.Nr = []
        self.plateID = []
        self.MJD = []
        self.fiberID = []
        self.z = []
        self.ra = []
        self.dec= []
        self.o2=[]#EW [OII]
        self.erro2=[]#err in EW [OII]
        self.u=[]#u mag, k-corrected to z=0.1
        self.g=[]#""
        self.r=[]#""
        self.i=[]#""
        self.zm=[]#""
        self.gpetrext=[]#g_petr corrected for extinction
        self.rpetrext=[]
        self.ipetrext=[]
        self.zpetrext=[]
        self.errgpetrext=[]#error in g_petr corrected for extinction
        self.errrpetrext=[]
        self.erripetrext=[]
        self.errzpetrext=[]
        self.V=[]#Vega mag, z=0.0, Blanton
        self.distBCG = []#distance to BCG in degrees
        self.distBCGR200=[]#distance to BCG in units of R200
        self.dz=[]#redshift - clusterz
    def greadfiles(self,clusters,i):
        infile=open(clusters,'r')
        distmod=(5.*N.log10(c.dL[i])+25.)
        for line in infile:
            if line.find('#') > -1:
                continue
            fields=line.split()
            #print fields[25]
            dLgal=my.dL(float(fields[4]),h100)
            dm=5.*N.log10(dLgal/c.dL[i])#correct abs mag to cluster redshift
            m=float(fields[22])+dm
            m=m-0.77#convert to h100=0.7
             

            if (m > mabscut):#keep only galaxies that meet mag cut
                continue


            #apply color selection to mimic ediscs selection
            distmod=5.*N.log10(dLgal)+25.
            g=float(fields[10])-0.77+ distmod
            r=float(fields[11])-0.77+ distmod
            if (g > 18.0):
                #print "didn't make g cut, g = ",g,fields[10],m,distmod
                continue
            if (r > 17.7):
                #print "didn't make r cut, g = ",r,fields[11],m
                continue
            
            if (abs(float(fields[25])/(c.z[i]+1)*3.e5) > 3*c.sigma[i]):#keep only galaxies w/in 3 sigma
                continue
            self.z.append(float(fields[4]))
            self.ra.append(float(fields[5]))
            self.dec.append(float(fields[6]))
            self.o2.append(float(fields[7]))
            self.erro2.append(float(fields[8]))
            self.u.append(float(fields[9]))
            self.g.append(float(fields[10]))
            self.r.append(float(fields[11]))
            self.i.append(float(fields[12]))
            self.zm.append(float(fields[13]))
            self.gpetrext.append(float(fields[14]))
            self.rpetrext.append(float(fields[15]))
            self.ipetrext.append(float(fields[16]))
            self.zpetrext.append(float(fields[17]))
            self.V.append(m)
            self.distBCG.append(float(fields[23]))
            self.distBCGR200.append(float(fields[23]))
            self.dz.append(float(fields[25]))

        infile.close()
    def convarray(self):
        self.z = N.array(self.z,'f')
        self.ra = N.array(self.ra,'f')
        self.dec = N.array(self.dec,'f')
        self.distBCG = N.array(self.distBCG,'f')
        self.distBCGR200 = N.array(self.distBCGR200,'f')
        self.dz = N.array(self.dz,'f')
        self.o2 = N.array(self.o2,'f')
        self.erro2 = N.array(self.erro2,'f')
        self.u = N.array(self.u,'f')
        self.g = N.array(self.g,'f')
        self.r = N.array(self.r,'f')
        self.i = N.array(self.i,'f')
        self.zm = N.array(self.zm,'f')
        self.V = N.array(self.V,'f')
        #self.V=self.V-0.77#convert to h=0.7
        self.memb = N.zeros(len(self.distBCGR200),'f')
        self.sf = N.zeros(len(self.distBCGR200),'f')
        self.tot = N.ones(len(self.distBCGR200),'f')
        for i in range(len(self.memb)):
            if (self.distBCGR200[i] < 1):
                self.memb[i]=1.
            if (self.o2[i] < -4.):
                self.sf[i]=1.
    def cullmembers(self):
        self.z = N.compress(self.memb>0,self.z)
        self.ra = N.compress(self.memb>0,self.ra)
        self.dec =N.compress(self.memb>0,self.dec)
        self.distBCG = N.compress(self.memb>0,self.distBCG)
        self.distBCGR200 = N.compress(self.memb>0,self.distBCGR200)
        self.dz = N.compress(self.memb>0,self.dz)
        self.o2 = N.compress(self.memb>0,self.o2)
        self.erro2 = N.compress(self.memb>0,self.erro2)
        self.u = N.compress(self.memb>0,self.u)
        self.g = N.compress(self.memb>0,self.g)
        self.r = N.compress(self.memb>0,self.r)
        self.i = N.compress(self.memb>0,self.i)
        self.zm = N.compress(self.memb>0,self.zm)

        self.sf = N.compress(self.memb>0,self.sf)
        self.tot =N.compress(self.memb>0,self.tot)
        self.memb = N.compress(self.memb>0,self.memb)
    def greadphotfiles(self,gphotfile,dL,kcorr):
        self.photra=[]
        self.photdec=[]
        self.phu=[]
        self.phMv=[]
        for line in open(gphotfile):
            if line.find('#') > -1:
                continue
            if len(line) <10:#phot files have blank line at end - skip it
                continue
            #print line
            fields=line.split(',')
            #print fields[0], fields[1], line
            #need to implement a mag cut
            g=float(fields[3])
            r=float(fields[4])
            if (g > 18.0):
                continue
            if (r > 17.7):
                continue
            g01=g+0.3134+0.4608*((g-r)* -0.6148)
            g01=g01+.01#convert to vega
            r01=g-0.4118-.8597*((g-r)* -0.6148)
            r01=r01-.04#convert to vega
            V=g01-0.6689-0.9264*((g01-r01)-0.7252)

            g=g-(-0.08)#convert to vega
            r=r-(.16)#convert to vega
            #g=g-float(fields[10])#correct for extinction
            #r=r-float(fields[11])#correct for extinction

            V=g-0.3556-0.7614*((g-r)-0.6148)

            #mv=V - (5.*N.log10(dL)+25.+kcorr)
            mv=V- (5.*N.log10(dL)+25.)
            if ((mv) < mabscut):
                self.photra.append(float(fields[0]))
                self.photdec.append(float(fields[1]))
                self.phu.append(float(fields[2]))
                self.phMv.append(float(mv))
            #self.photra.append(float(fields[0]))
            #self.photdec.append(float(fields[1]))
            #self.phu.append(float(fields[2]))
            #self.phMv.append(float(mv))

        self.photra=N.array(self.photra,'f')
        self.photdec=N.array(self.photdec,'f')
        self.phu=N.array(self.phu,'f')
        self.phMv=N.array(self.phMv,'f')
    def greadspecfiles(self,gphotfile,dL,kcorr,j):
        self.sra=[]
        self.sdec=[]
        self.su=[]
        self.sg=[]
        self.sr=[]
        self.sV=[]
        self.so2=[]
        self.sha=[]
        for line in open(gphotfile):
            if line.find('#') > -1:
                continue
            if len(line) <10:#phot files have blank line at end - skip it
                continue
            #print line
            fields=line.split(',')
            #print fields[0], fields[1]
            #need to implement a mag cut
            g=float(fields[3])
            r=float(fields[4])
            if (g > 18.0):
                continue
            if (r > 17.7):
                continue

            g01=g+0.3134+0.4608*((g-r)* -0.6148)
            g01=g01+.01#convert to vega
            r01=g-0.4118-.8597*((g-r)* -0.6148)
            r01=r01-.04#convert to vega
            V=g01-0.6689-0.9264*((g01-r01)-0.7252)

            g=g-(-0.08)#convert to vega
            r=r-(.16)#convert to vega
            #g=g-float(fields[10])#correct for extinction
            #r=r-float(fields[11])#correct for extinction

            V=g-0.3556-0.7614*((g-r)-0.6148)

            #mv=V - (5.*N.log10(dL)+25.+kcorr)
            mv=V- (5.*N.log10(dL)+25.)
            if ((mv) < mabscut):
                self.sra.append(float(fields[0]))
                self.sdec.append(float(fields[1]))
                self.su.append(float(fields[2]))
                self.sg.append(float(fields[3]))
                self.sr.append(float(fields[4]))
                self.sV.append(float(mv))
                self.so2.append(float(fields[14]))#14-o2, 16 - Halpha
                self.sha.append(float(fields[16]))#14-o2, 16 - Halpha

        self.sra=N.array(self.sra,'f')
        self.sdec=N.array(self.sdec,'f')
        self.su=N.array(self.su,'f')
        self.sg=N.array(self.sg,'f')
        self.sr=N.array(self.sr,'f')
        self.sV=N.array(self.sV,'f')
        self.so2=N.array(self.so2,'f')
        self.sha=N.array(self.sha,'f')
        self.sdistBCGR200 = N.sqrt((c.ra[j]-self.sra)**2 + (c.dec[j]-self.sdec)**2)/c.r200deg[j]
        self.smemb = N.zeros(len(self.sdistBCGR200),'f')
        self.ssf = N.zeros(len(self.sdistBCGR200),'f')
        self.stot = N.ones(len(self.sdistBCGR200),'f')
        for i in range(len(self.smemb)):
            if (self.sdistBCGR200[i] < 4):
                self.smemb[i]=1.
            if (self.so2[i]> 4.):
                self.ssf[i]=1.
            if (self.so2[i]< -50.):
                self.smemb[i]=0.

    def getnearest(self,j):
        self.sig5=N.zeros(len(self.ra),'f')
        self.sig10=N.zeros(len(self.ra),'f')
        self.sig5phot=N.zeros(len(self.ra),'f')
        self.sig10phot=N.zeros(len(self.ra),'f')
        self.nearest=N.zeros(len(self.ra),'f')#distance to nearest neighbor
        self.dmagnearest=N.zeros(len(self.ra),'f')#mag diff b/w spec obj and nearest
        for i in range(len(self.ra)):
            angdist=my.DA(c.z[j],h100)#kpc/arcsec
            if self.memb[i] > 0:
                dspec=N.sqrt((self.ra[i]-self.ra)**2+(self.dec[i]-self.dec)**2)#sorted array of distances in degrees
                dphot=N.sqrt((self.ra[i]-self.photra)**2+(self.dec[i]-self.photdec)**2)#sorted array of distances in degrees
                Mvsort=N.take(self.V,N.argsort(dspec))
                phMvsort=N.take(self.phMv,N.argsort(dphot))
                dspecsort=N.take(dspec,N.argsort(dspec))
                self.sig5[i]=5./(N.pi)/(dspecsort[5]*3600.*angdist/1000.)**2#convert from deg to arcsec, multiply by DA (kpc/arcsec), divide by 1000 to convert to Mpc, index 5 element b/c 0 is itself
                if (len(dspecsort) > 10):
                    self.sig10[i]=10./(N.pi)/(dspecsort[10]*3600.*angdist/1000.)**2
                else:
                    self.sig10[i]=0.
                dphotsort=N.take(dphot,N.argsort(dphot))
                self.nearest[i]=dphotsort[0]
                self.dmagnearest[i]=phMvsort[0]-Mvsort[0]
                self.sig5phot[i]=5./(N.pi)/(dphotsort[5]*3600.*angdist/1000.)**2
                self.sig10phot[i]=10./(N.pi)/(dphotsort[10]*3600.*angdist/1000.)**2
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
    def gwritefiles(self,clusters,j):
        outfile=clusters[:(len(clusters)-4)]+'2.dat'
        out=open(outfile,'w')
        infile=open(clusters,'r')
        i=0
        ds5=[]
        ds10=[]
        k=0
        for line in infile:
            if line.find('#') > -1:
                out.write(line)
                continue
            if k < 1:
                out.write('#  27 Sigma_5 calculated from spec sample \n')
                out.write('#  28 Sigma_10 calculated from spec sample \n')
                out.write('#  29 Sigma_5 calculated from phot sample \n')
                out.write('#  30 Sigma_10 calculated from phot sample \n')
                out.write('#  31 Sigma_5 calculated from spec sample extended to 2xR200 \n')
                out.write('#  32 Sigma_10 calculated from spec sample extended to 2xR200\n')
                out.write('#  33 M_V (h100 = 0.7) corrected to cluster redshift \n')
                k=1
            fields=line.split()
            dLgal=my.dL(float(fields[4]),h100)
            dm=5.*N.log10(dLgal/c.dL[j])#correct abs mag to cluster redshift
            m=float(fields[22])+dm
            m=m-0.77#convert to h100=0.7

            if (m > mabscut):#keep only galaxies that meet mag cut
                out.write(line)
                #print "didn't make mag cut",m,mabscut,
                continue
            distmod=5.*N.log10(dLgal)+25.
            g=float(fields[10])-0.77+ distmod
            r=float(fields[11])-0.77+ distmod
            if (g > 18.0):
                #print "didn't make g cut, g = ",g,fields[10],m,distmod
                continue
            if (r > 17.7):
                #print "didn't make r cut, g = ",r,fields[11],m
                continue

            if (abs(float(fields[25])/(c.z[j]+1)*3.e5) > 3*c.sigma[j]):#keep only galaxies w/in 3 sigma
                out.write(line)
                #print "didn't make vel cut"
                continue

            #print i, 'writing output ',len(self.sig5),len(self.sig5),len(self.sig10),len(self.sig5phot),len(self.sig10phot),len(line)
            nline=line[:(len(line)-2)]+" %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n"%(self.sig5[i],self.sig10[i],self.sig5phot[i],self.sig10phot[i], self.sig52r200[i],self.sig102r200[i], self.V[i])
            d5=(self.sig5phot[i]-self.sig5[i])
            d10=(self.sig10phot[i]-self.sig10[i])
            #print " %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f \n"%((self.sig5phot[i]-self.sig5[i]),(self.sig10phot[i]-self.sig10[i]),self.sig5[i],self.sig10[i],self.sig5phot[i],self.sig10phot[i])
            out.write(nline)
            ds5.append(d5)
            ds10.append(d10)
            i=i+1
        
        out.close()
        infile.close()
        return ds5, ds10
            
"""
plot fraction of OII EW > 4 versus sigma 5, where sigma 5 is measured from spec members only
"""
def plotsigsff(sig,sf,file,nbin):

    psplot=file+".ps"
    psplotinit(psplot)
    tot=N.ones(len(sf),'f')
    (sigbin,sfbin)=my.binitsumequal(sig,sf,nbin)
    (sigbin,totbin)=my.binitsumequal(sig,tot,nbin)
    print sfbin
    print totbin
    (sff,sfferr)=my.ratioerror(sfbin,totbin)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ymin=-.05
    ymax=1.05
    xmin=min(sig)-10.
    #xmax=max(sig)-200.
    xmax=350.
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0)
    ppgplot.pglab("\gS\d5\u (gal/Mpc\u2\d)","Fraction EW([OII])>4 \(2078)","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    sig=N.array(sig,'f')
    sff=N.array(sff,'f')
    ppgplot.pgsci(2)
    ppgplot.pgline(sigbin,sff)
    ppgplot.pgsci(1)

    ppgplot.pgpt(sigbin,sff,17)
    my.errory(sigbin,sff,sfferr)
    ppgplot.pgend()
def plotsigsffall(sigspec,sigphot,sf,file,nbin):

    psplot=file+".ps"
    psplotinit(psplot)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ymin=-.01
    ymax=1.01
    #xmin=min(sigspec)-10.
    #xmax=max(sig)-200.
    #xmax=400.
    xmin=0.2
    xmax=2.7
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,10)
    ppgplot.pglab("\gS\d5\u (gal/Mpc\u2\d)","Fraction EW([OII])>4 \(2078)","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    tot=N.ones(len(sf),'f')
    (sigbin,sfbin)=my.binitsumequal(sigspec,sf,nbin)
    (sigbin,totbin)=my.binitsumequal(sigspec,tot,nbin)
    (sff,sfferr)=my.ratioerror(sfbin,totbin)
    #sig=N.array(sig,'f')
    #sff=N.array(sff,'f')
    ppgplot.pgsci(2)
    sigbin=N.log10(sigbin)
    ppgplot.pgline(sigbin,sff)
    ppgplot.pgsci(1)

    ppgplot.pgpt(sigbin,sff,17)
    my.errory(sigbin,sff,sfferr)
    
    (sigbin,sfbin)=my.binitsumequal(sigphot,sf,nbin)
    (sigbin,totbin)=my.binitsumequal(sigphot,tot,nbin)
    (sff,sfferr)=my.ratioerror(sfbin,totbin)
    #sig=N.array(sig,'f')
    #sff=N.array(sff,'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgsci(4)
    sigbin=N.log10(sigbin)
    ppgplot.pgline(sigbin,sff)
    ppgplot.pgsci(1)

    ppgplot.pgpt(sigbin,sff,21)
    #my.errory(sigbin,sff,sfferr)
    ppgplot.pgend()
def plotsig10sffall(sigspec,sigphot,sf,file,nbin):

    psplot=file+".ps"
    psplotinit(psplot)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ymin=-.01
    ymax=1.01
    #xmin=min(sigspec)-10.
    #xmax=max(sig)-200.
    #xmax=400.
    xmin=-1.
    xmax=2.7
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,10)
    ppgplot.pglab("\gS\d10\u (gal/Mpc\u2\d)","Fraction EW([OII])>4 \(2078)","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    tot=N.ones(len(sf),'f')
    (sigbin,sfbin)=my.binitsumequal(sigspec,sf,nbin)
    (sigbin,totbin)=my.binitsumequal(sigspec,tot,nbin)
    (sff,sfferr)=my.ratioerror(sfbin,totbin)
    #sig=N.array(sig,'f')
    #sff=N.array(sff,'f')
    ppgplot.pgsci(2)
    sigbin=N.log10(sigbin)
    ppgplot.pgline(sigbin,sff)
    ppgplot.pgsci(1)

    ppgplot.pgpt(sigbin,sff,17)
    my.errory(sigbin,sff,sfferr)
    
    (sigbin,sfbin)=my.binitsumequal(sigphot,sf,nbin)
    (sigbin,totbin)=my.binitsumequal(sigphot,tot,nbin)
    (sff,sfferr)=my.ratioerror(sfbin,totbin)
    #sig=N.array(sig,'f')
    #sff=N.array(sff,'f')
    ppgplot.pgslw(4)  #line width
    ppgplot.pgsci(4)
    sigbin=N.log10(sigbin)
    ppgplot.pgline(sigbin,sff)
    ppgplot.pgsci(1)

    ppgplot.pgpt(sigbin,sff,21)
    #my.errory(sigbin,sff,sfferr)
    ppgplot.pgend()

def plotsigo2(sig,o2,file,nbin):

    psplot=file+".ps"
    psplotinit(psplot)
    (sigbin,o2bin)=my.binit(sig,o2,nbin)
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
    ymin=-10.
    ymax=2.
    xmin=min(sig)-10.
    #xmax=max(sig)-200.
    xmax=350.
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0)
    ppgplot.pglab("\gS\d5\u (gal/Mpc\u2\d)","EW([OII]) (\(2078))","")
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    sig=N.array(sig,'f')
    o2=N.array(o2,'f')
    ppgplot.pgpt(sig,o2,1)
    ppgplot.pgsci(2)
    ppgplot.pgline(sigbin,o2bin)
    ppgplot.pgsci(1)
    ppgplot.pgend()
def plotsigo2all(sig,psig,o2b,file,nbin):
    #o2=N.zeros(len(o2b),'f')
    #for i in range(len(o2b)):
        #print i, sig[i], psig[i], o2b[i]
    #    if o2b[i] < 0:

    #        o2[i]=-1*o2b[i]
            #print "hey", o2[i] 
    o2=o2b
    psplot=file+".ps"
    psplotinit(psplot)
    ppgplot.pgsch(0.7)
    (sigbin,o2bin)=my.binit(sig,o2,nbin)
    #print 'dude', sigbin, o2bin
    sigbin=N.log10(sigbin)
    ppgplot.pgswin(-1.,3.,-.5,10.)
    ppgplot.pgbox('bcnst',0.0,0.0,'bcvnst',0.0,0.0)  #tickmarks and labeling
    ppgplot.pgsch(1.0)
    ppgplot.pgmtxt('b',2.5,0.5,0.5,"\gS\d10\u (gal/Mpc\u2\d)" )    #xlabel
    ppgplot.pgsch(1.2)
    ppgplot.pgmtxt('l',2.6,0.5,0.5,'EW([OII]) (\(2078))')

    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    ppgplot.pgsci(2)
    ppgplot.pgpt(sigbin,o2bin,17)
    #print 'dude2', sigbin, o2bin
    ppgplot.pgline(sigbin,o2bin)
    (sigbin,o2bin)=my.binit(psig,o2,nbin)
    #print 'dude', sigbin, o2bin
    sigbin=N.log10(sigbin)
    ppgplot.pgsci(4)
    ppgplot.pgpt(sigbin,o2bin,21)
    ppgplot.pgline(sigbin,o2bin)
    ppgplot.pgsci(1)
    ppgplot.pgend()
def plotsighaall(sig,psig,o2b,file,nbin):
    o2b=N.array(o2b,'f')
    sig=N.array(sig,'f')
    psig=N.array(psig,'f')
    #o2b=o2b+4.
    o2=N.compress(o2b > -500., o2b)
    sig=N.compress(o2b > -500., sig)
    psig=N.compress(o2b > -500., psig)
    
    psplot=file+".ps"
    psplotinit(psplot)
    #ppgplot.pgsch(0.7)
    ppgplot.pgslw(7)
    (sigbin,o2bin)=my.binit(sig,o2,nbin)
    #print 'dude', sigbin, o2bin
    sigbin=N.log10(sigbin)
    ppgplot.pgswin(-2.,2.,-5.,20.)
    ppgplot.pgbox('blcnst',0.0,0.0,'bcvnst',0.0,0.0)  #tickmarks and labeling
    ppgplot.pgsch(1.0)
    ppgplot.pgmtxt('b',2.5,0.5,0.5,"\gS\d10\u (gal/Mpc\u2\d)" )    #xlabel
    ppgplot.pgsch(1.2)
    ppgplot.pgmtxt('l',2.6,0.5,0.5,'EW(H\ga) (\(2078))')

    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width

    ppgplot.pgpt(sigbin,o2bin,17)
    ppgplot.pgpt(N.log10(sig),o2,1)
    #my.errory(sigbin,o2bin,yerr)
    #print 'dude2', sigbin, o2bin
    ppgplot.pgsci(2)
    ppgplot.pgline(sigbin,o2bin)
    (sigbin,o2bin)=my.binit(psig,o2,nbin)

    #print 'dude', sigbin, o2bin
    sigbin=N.log10(sigbin)
    ppgplot.pgsci(1)
    ppgplot.pgpt(sigbin,o2bin,21)
    #my.errory(sigbin,o2bin,yerr)
    ppgplot.pgsci(4)
    ppgplot.pgline(sigbin,o2bin)
    ppgplot.pgsci(1)
    ppgplot.pgend()


def gotoit():
    nbin=10
    #c=Cluster()
    #g=Galaxy()
    clusterfile="clusters.spec.dat"
    print "reading in cluster file to get cluster parameters"
    c.creadfiles(clusterfile)
    print "got ",len(c.z)," clusters"
    c.convarray()
    c.Kcorr()

    go2=[]#combined arrays containing all galaxies
    gsf=[]#combined arrays containing all galaxies
    gsig5=[]
    gsig10=[]
    gsig52r200=[]#spec catalogs extended out to 2xR200
    gsig102r200=[]#spec catalogs extended out to 2xR200
    gsig5phot=[]
    gsig10phot=[]
    sgo2=[]#combined arrays containing all galaxies
    sgha=[]#combined arrays containing all galaxies
    sgsf=[]#combined arrays containing all galaxies
    sgsig5=[]
    sgsig10=[]
    sgsig52r200=[]#spec catalogs extended out to 2xR200
    sgsig102r200=[]#spec catalogs extended out to 2xR200
    sgsig5phot=[]
    sgsig10phot=[]

    if (mode < 1):
        c.getsdssphotcats()
        c.getsdssspeccats()

    gr=[]#list of median g-r colors
    psplotinit('summary.ps')
    x1=.1
    x2=.45
    x3=.6
    x4=.95
    y1=.15
    y2=.45
    y3=.55
    y4=.85
    ppgplot.pgsch(1.2) #font size
    ppgplot.pgslw(2)
    #for i in range(len(c.z)):
    cl=[10]
    (xl,xu,yl,yu)=ppgplot.pgqvp(0)
    print "viewport = ",xl,xu,yl,yu
    complall=[]
    for i in range(len(c.z)):
    #for i in cl:
        gname="g"+str(i)
        gname=Galaxy()
        gspecfile="abell"+str(c.id[i])+".spec.dat"
        gname.greadfiles(gspecfile,i)
        print "number of members = ", len(gname.z)
        if len(gname.z) < 10:
            print "less than 10 members", len(gname.z)
            continue
        gname.convarray()
        #gname.cullmembers()
        #gname.getmemb()#get members w/in R200
        #gr.append(N.average(gname.g-gname.r))

        gspec2r200file="abell"+str(c.id[i])+".spec2r200.dat"
        gname.greadspecfiles(gspec2r200file,c.dL[i],c.kcorr[i],i)
        print i, c.id[i]," getnearest, first call", len(gname.ra),len(gname.sra),sum(gname.smemb)
        #gname.getnearest(i)
        (gname.sig52r200,gname.sig102r200)=gname.getnearestgen(gname.ra,gname.dec,gname.sra,gname.sdec,i)#measure distances from ra1, dec1 to members in catalog ra2, dec2
        sig52r200=N.compress(gname.memb > 0,gname.sig52r200)
        gsig52r200[len(gsig5phot):]=sig52r200
        sig102r200=N.compress(gname.memb > 0,gname.sig102r200)
        gsig102r200[len(gsig10phot):]=sig102r200


        gphotfile="abell"+str(c.id[i])+".phot.dat"
        gname.greadphotfiles(gphotfile,c.dL[i],c.kcorr[i])
        gname.getnearest(i)
        #print "len of local density arrays = ",len(gname.sig5),len(gname.sig5phot)
        #print gspecfile, c.z[i],c.kcorr[i]
        (ds5,ds10)=gname.gwritefiles(gspecfile,i)
        o2=N.compress(gname.memb > 0,gname.o2)
        go2[len(go2):]=o2
        sf=N.compress(gname.memb > 0,gname.sf)
        gsf[len(gsf):]=sf
        sig5=N.compress(gname.memb > 0,gname.sig5)
        gsig5[len(gsig5):]=sig5
        sig10=N.compress(gname.memb > 0,gname.sig10)
        gsig10[len(gsig10):]=sig10
        sig5phot=N.compress(gname.memb > 0,gname.sig5phot)
        gsig5phot[len(gsig5phot):]=sig5phot
        sig10phot=N.compress(gname.memb > 0,gname.sig10phot)
        gsig10phot[len(gsig10phot):]=sig10phot


        ds5=N.array(ds5,'f')
        ds10=N.array(ds10,'f')
        #print len(ds5),len(ds10)
        #ppgplot.pgsvp(xl,xu,yl,yu)
        ppgplot.pgsvp(0.1,.9,.08,.92)
        ppgplot.pgslw(7)
        label='Abell '+str(c.id[i])+' (z=%5.2f, \gs=%3.0f km/s)'%(c.z[i],c.sigma[i])
        ppgplot.pgtext(0.,1.,label)
        ppgplot.pgslw(2)
        ppgplot.pgsvp(x1,x2,y1,y2)  #sets viewport
        #ppgplot.pgbox("",0.0,0,"",0.0)
        ppgplot.pgswin(-1.,3.,-1.,3.) #axes limits
        ppgplot.pgbox('bcnst',1,2,'bcvnst',1,2)  #tickmarks and labeling
        ppgplot.pgmtxt('b',2.5,0.5,0.5,"\gS\d10\u(phot) (gal/Mpc\u2\d)")    #xlabel
        ppgplot.pgmtxt('l',2.6,0.5,0.5,"\gS\d10\u(spec) (gal/Mpc\u2\d)")

        x=N.arange(-5.,10.,.1)
        y=x
        ppgplot.pgsls(1)#dotted
        ppgplot.pgslw(4)  #line width
        ppgplot.pgline(x,y)
        x=N.log10(sig10phot)
        y=N.log10(sig10)
        ppgplot.pgsch(.7)
        ppgplot.pgpt(x,y,17)
        xp=N.array([-0.5],'f')
        yp=N.array([2.5],'f')
        ppgplot.pgpt(xp,yp,17)
        ppgplot.pgtext((xp+.1),yp,'spec(1.2xR200) vs phot')
        ppgplot.pgsci(4)
        xp=N.array([-0.5],'f')
        yp=N.array([2.2],'f')
        ppgplot.pgpt(xp,yp,21)
        ppgplot.pgtext((xp+.1),yp,'spec(2xR200) vs phot')

        y=N.log10(sig102r200)

        ppgplot.pgsch(.9)
        ppgplot.pgpt(x,y,21)
        ppgplot.pgsch(1.2)
        ppgplot.pgslw(2)  #line width
        ppgplot.pgsci(1)

        #ppgplot.pgenv(-200.,200.,-1.,20.,0,0)
        #ppgplot.pgsci(2)
        #ppgplot.pghist(len(ds5),ds5,-200.,200.,30,1)
        #ppgplot.pgsci(4)
        #ppgplot.pghist(len(ds10),ds10,-200.,200.,30,1)
        #ppgplot.pgsci(1)
        #ppgplot.pglab("\gD\gS","Ngal",gspecfile)
        #ppgplot.pgpanl(1,2)
        g=N.compress(gname.memb>0,gname.g)
        r=N.compress(gname.memb>0,gname.r)
        V=N.compress(gname.memb>0,gname.V)
        dmag=N.compress(gname.memb>0,gname.dmagnearest)
        dnearest=N.compress(gname.memb>0,gname.nearest)
        dz=N.compress(gname.memb>0,gname.dz)
        #ppgplot.pgsvp(x3,x4,y1,y2)  #sets viewport
        #ppgplot.pgenv(-.5,3.,-1.,5.,0,0)
        #ppgplot.pgpt((g-V),(g-r),17)
        #ppgplot.pgsci(1)
        #ppgplot.pglab("g - M\dV\u",'g-r',gspecfile)
        ppgplot.pgsvp(x1,x2,y3,y4)  #sets viewport
        #ppgplot.pgbox("",0.0,0,"",0.0)
        ppgplot.pgswin((c.ra[i]+2.*c.r200deg[i]/N.cos(c.dec[i]*N.pi/180.)),(c.ra[i]-2*c.r200deg[i]/N.cos(c.dec[i]*N.pi/180.)),(c.dec[i]-2.*c.r200deg[i]),(c.dec[i]+2.*c.r200deg[i]))
        ppgplot.pgbox('bcnst',0.0,0.0,'bcvnst',0.0,0.0)  #tickmarks and labeling
        ppgplot.pgmtxt('b',2.5,0.5,0.5,"RA")    #xlabel
        ppgplot.pgmtxt('l',2.6,0.5,0.5,"Dec")

        #ppgplot.pglab("RA",'Dec',gspecfile)
        ppgplot.pgsfs(2)
        ppgplot.pgcirc(c.ra[i],c.dec[i],c.r200deg[i])
        ppgplot.pgsls(4)
        ppgplot.pgcirc(c.ra[i],c.dec[i],1.2*c.r200deg[i])
        ppgplot.pgsls(1)
        #ppgplot.pgcirc(c.ra[i],c.dec[i],c.r200deg[i]/N.cos(c.dec[i]*N.pi/180.))
        ppgplot.pgsci(2)
        ppgplot.pgpt(gname.ra,gname.dec,17)
        ppgplot.pgsci(4)
        ppgplot.pgpt(gname.photra,gname.photdec,21)
        ppgplot.pgsci(1)

	#calculate completeness w/in R200

	dspec=N.sqrt((gname.ra-c.ra[i])**2 +(gname.dec-c.dec[i])**2)
	dphot=N.sqrt((gname.photra-c.ra[i])**2 +(gname.photdec-c.dec[i])**2)
	nphot=1.*len(N.compress(dphot < c.r200deg[i],dphot))
	nspec=1.*len(N.compress(dspec < c.r200deg[i],dspec))
	s="Completeness for cluster Abell %s = %6.2f (nspec=%6.1f,nphot= %6.1f)" % (str(c.id[i]), float(nspec/nphot),nspec,nphot)
	print s
	complall.append(float(nspec/nphot))
        ppgplot.pgsvp(x3,x4,y3,y4)  #sets viewport
        #ppgplot.pgsvp(x1,x2,y3,y4)  #sets viewport
        #ppgplot.pgbox("",0.0,0,"",0.0)
        ppgplot.pgswin(-0.005,.05,-1.,1.)
        ppgplot.pgbox('bcnst',.02,2,'bcvnst',1,4)  #tickmarks and labeling
        ppgplot.pgsch(1.0)
        ppgplot.pgmtxt('b',2.5,0.5,0.5,"Dist to nearest phot neighbor (deg)" )    #xlabel
        ppgplot.pgsch(1.2)
        ppgplot.pgmtxt('l',2.6,0.5,0.5,'M\dV\u(phot) - M\dV\u(spec)')
        ppgplot.pgsci(2)
        ppgplot.pgpt(dnearest,dmag,17)
        ppgplot.pgsci(1)
        x=N.arange(-30.,30.,1.)
        y=0*x
        ppgplot.pgsci(1)
        ppgplot.pgsls(2)
        ppgplot.pgline(x,y)
        ppgplot.pgsls(1)
        ppgplot.pgsci(1)
        dm=N.compress(dnearest < 0.01, dmag)
        std='%5.3f (%5.3f)'%(pylab.mean(dm),pylab.std(dm))
        #ppgplot.pgslw(7)
        #label='Abell '+str(c.id[i])
        #ppgplot.pgtext(0.,1.,label)
        ppgplot.pgslw(2)
        label = '\gDM\dV\u(err) = '+std
        ppgplot.pgsch(.9)
        ppgplot.pgtext(0.,.8,label)
        #label = "z = %5.2f"%(c.z[i])
        #ppgplot.pgtext(0.,.8,label)
        ppgplot.pgsch(1.2)
        #ppgplot.pgsvp(x3,x4,y3,y4)  #sets viewport
        #ppgplot.pgenv(-.15,.15,-3.,3.,0,0)
        #ppgplot.pgsci(2)
        #ppgplot.pgpt(dz,dmag,17)
        #ppgplot.pgsci(1)
        #ppgplot.pglab("z-z\dcl\u",'\gD Mag',gspecfile)
        ppgplot.pgsvp(x3,x4,y1,y2)  #sets viewport
        ppgplot.pgswin(-3.,3.,-1.,1.)
        ppgplot.pgbox('bcnst',1,2,'bcvnst',1,4)  #tickmarks and labeling
        ppgplot.pgmtxt('b',2.5,0.5,0.5,"\gDv/\gs")    #xlabel
        ppgplot.pgmtxt('l',2.6,0.5,0.5,'M\dV\u(phot) - M\dV\u(spec)')

        ppgplot.pgsci(2)
        dv=dz/(1+c.z[i])*3.e5/c.sigma[i]
        ppgplot.pgpt(dv,dmag,17)
        ppgplot.pgsci(1)
        x=N.arange(-30.,30.,1.)
        y=0*x
        ppgplot.pgsci(1)
        ppgplot.pgsls(2)
        ppgplot.pgline(x,y)
        ppgplot.pgsls(1)
        ppgplot.pgsci(1)
        #ppgplot.pgsvp(x1,x2,y1,y2)  #sets viewport
        #ppgplot.pgenv(0.,3.5,-3.,3.,0,0)
        #ppgplot.pgsci(4)
        #ppgplot.pgpt((g-r),dmag,17)
        #ppgplot.pgsci(1)
        #ppgplot.pglab("g-r",'\gD Mag',gspecfile)

        #ppgplot.pgsvp(x1,x2,y1,y2)  #sets viewport
        #ppgplot.pgenv(-25.,-18.,-1.,1.,0,0)
        #ppgplot.pgsci(4)
        #ppgplot.pgpt((V),dmag,17)
        #x=N.arange(-30.,30.,1.)
        #y=0*x
        #ppgplot.pgsci(1)
        #ppgplot.pgsls(2)
        #ppgplot.pgline(x,y)
        #ppgplot.pgsls(1)
        #ppgplot.pgsci(1)
        #ppgplot.pglab("M\dV\u(spec)",'M\dV\u(phot) - M\dV\u(spec)',gspecfile)
        #ppgplot.pgpage()
        #ppgplot.pgpage()
        #combine galaxy data
        ppgplot.pgpage()

        (sssig5,sssig10)=gname.getnearestgen(gname.sra,gname.sdec,gname.sra,gname.sdec,i)#get spec-spec local density
        (spsig5,spsig10)=gname.getnearestgen(gname.sra,gname.sdec,gname.photra,gname.photdec,i)#get spec-phot local density


        o2=N.compress(gname.smemb > 0,gname.so2)
        sgo2[len(sgo2):]=o2
        ha=N.compress(gname.smemb > 0,gname.sha)
        sgha[len(sgha):]=ha
        sf=N.compress(gname.smemb > 0,gname.ssf)
        sgsf[len(sgsf):]=sf
        sig5=N.compress(gname.smemb > 0,sssig5)
        sgsig5[len(sgsig5):]=sig5
        sig10=N.compress(gname.smemb > 0,sssig10)
        sgsig10[len(sgsig10):]=sig10
        sig5phot=N.compress(gname.smemb > 0,spsig5)
        sgsig5phot[len(sgsig5phot):]=sig5phot
        sig10phot=N.compress(gname.smemb > 0,spsig10)
        sgsig10phot[len(sgsig10phot):]=sig10phot

    #gr=N.array(gr,'f')
    #c.assigncolor(gr)

    #for i in range(len(c.z)):
    #    print c.id[i],c.z[i],c.r200[i],c.r200deg[i]

    print "Average Completeness w/in R200 = ",N.average(N.array(complall,'f'))
    print "sig o2",len(gsig10),len(gsig10phot),len(go2)
    print "sig o2 large",len(sgsig10),len(sgsig10phot),len(sgo2)
    plotsigo2all(gsig10,gsig10phot,go2,'o2vsig10spec',nbin)
    #plotsigo2(gsig5phot,-1*go2,'o2vsig5phot',nbin)
    plotsigsff(gsig5,gsf,'sffvsig5spec',nbin)#sf frac versus sigma
    plotsigsff(gsig5phot,gsf,'sffvsig5phot',nbin)#sf frac versus sigma
    plotsigsffall(gsig5,gsig5phot,gsf,'sffvsig5all',nbin)#sf frac versus sigma
    plotsig10sffall(gsig10,gsig10phot,gsf,'sffvsig10all',nbin)#sf frac versus sigma
    #plotsighaall(gsig10,gsig10phot,gha,'havsig10spec',20)
    #plotsigo2all(sgsig10,sgsig10phot,sgo2,'o2vsig10spec.large',30)
    plotsighaall(sgsig10,sgsig10phot,sgha,'havsig10spec.large',10)
    #plotsigsffall(sgsig5,sgsig5phot,sgsf,'sffvsig5.large',nbin)#sf frac versus sigma
    #plotsig10sffall(sgsig10,sgsig10phot,sgsf,'sffvsig10.large',nbin)#sf frac versus sigma
    psplotinit('one2one.ps')
    ppgplot.pgenv(-1.5,2.5,-1.5,2.5,0)
    ppgplot.pglab("\gS\d10\u(phot) (gal/Mpc\u2\d)","\gS\d10\u(spec) (gal/Mpc\u2\d)","")
    x=N.arange(-5.,10.,.1)
    y=x
    ppgplot.pgsls(1)#dotted
    ppgplot.pgslw(4)  #line width
    ppgplot.pgline(x,y)
    x=N.log10(gsig10phot)
    y=N.log10(gsig10)
    ppgplot.pgsch(.7)
    ppgplot.pgpt(x,y,17)
    ppgplot.pgsch(1.)
    ppgplot.pgsci(1)
    ppgplot.pgend()

    #do large area


c=Cluster()
g=Galaxy()
mode=0
gotoit()
endtime=time.clock()
print "end time = ",endtime
print "elapsed time = ", endtime-starttime
    


