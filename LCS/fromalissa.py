#!/usr/bin/env python

from pylab import *
import pyfits

clusterra={'MKW11':202.3800, 'MKW8':220.1796, 'AWM4':241.2375, 'A2063':230.7578, 'A2052':229.1896, 'NGC6107':244.333750, 'Coma':194.9531, 'A1367':176.1231, 'Hercules':241.3125}
clusterdec={'MKW11':11.78861, 'MKW8':3.4530, 'AWM4':23.9206, 'A2063':8.6394, 'A2052':7.0003, 'NGC6107':34.901389, 'Coma':27.9807, 'A1367':19.8391, 'Hercules':17.7485}
#clusterV={'MKW11':.022849*3.0e5, 'MKW8':.027000*3.0e5, 'AWM4':.031755*3.0e5,'A2063':.034937*3.0e5, 'A2052':.0354918*3.0e5,'NGC6107':.030658*3.0e5, 'Coma':.023*3.0e5, 'A1367':.028*3.0e5, 'Hercules':.037*3.0e5}#original
#clusterV={'MKW11':6850, 'MKW8':8094, 'AWM4':9520, 'A2063':10474, 'A2052':10640, 'NGC6107':9191, 'Coma':6925, 'A1367':6595, 'Hercules':10972}#NED
clusterV={'MKW11':7569.7,'MKW8':8296.1,'AWM4':10117.45,'A2063':10300.5,'A2052':10719.5,'NGC6107':9607.76,'Coma':7014.03,'A1367':6893.12,'Hercules':10934.95}#Cbi from hist program

######Numbers taken from dispersion.py
#clustercentralvelocity={'MKW11':6854, 'MKW8':8100, 'AWM4':9526, 'A2063':10481, 'A2052':10647, 'NGC6107':9197, 'Coma':6900, 'A1367':8400, 'Hercules':11100}#original
clustercentralvelocity={'MKW11':6850, 'MKW8':8094, 'AWM4':9520, 'A2063':10474, 'A2052':10640, 'NGC6107':9191, 'Coma':6925, 'A1367':6595, 'Hercules':10972}#NED

#clusterlow8={'MKW11':5816.85, 'MKW8':5758.61,'AWM4':6344.15,'A2063':8551,'A2052':8514.45,'NGC6107':7269.02, 'Coma':4534.67,'A1367':4008.8,'Hercules':8648.51}#updated dispersion July 6th

#clusterhigh8={'MKW11':7891.15, 'MKW8':10441.39,'AWM4':12707.85,'A2063':12411,'A2052':12779.55,'NGC6107':11124.989, 'Coma':9265.33,'A1367':9181.2,'Hercules':13551.49}#updated dispersion July 6th

#clusterlow8={'MKW11':5990,'MKW8':7540,'AWM4':8540,'A2063':8610,'A2052':9230,'NGC6107':8280,'Coma':4410,'A1367':4820,'Hercules':9830}#eyeball v cut
#clusterhigh8={'MKW11':8070,'MKW8':8880,'AWM4':10600,'A2063':12400,'A2052':12300,'NGC6107':10700,'Coma':8910,'A1367':8730,'Hercules':12800}#eyeball v cut

#clusterlow8={'MKW11':5816.85, 'MKW8':5758.61,'AWM4':6344.15,'A2063':8551,'A2052':8514.45,'NGC6107':7269.02, 'Coma':4534.67,'A1367':4008.8,'Hercules':8648.51}#NED

#clusterhigh8={'MKW11':7891.15, 'MKW8':10441.39,'AWM4':12707.85,'A2063':12411,'A2052':12779.55,'NGC6107':11124.989, 'Coma':9265.33,'A1367':9181.2,'Hercules':13551.49}#NED

#clusterlow8={'MKW11':5816.85, 'MKW8':5758.61,'AWM4':6344.15,'A2063':8551,'A2052':8514.45,'NGC6107':7269.02, 'Coma':4534.67,'A1367':3231.03,'Hercules':8648.51}#original
#clusterhigh8={'MKW11':7891.15, 'MKW8':10441.39,'AWM4':12707.85,'A2063':12411,'A2052':12779.55,'NGC6107':11124.989, 'Coma':9265.33,'A1367':13568.97,'Hercules':13551.49}#original

clustersbi={'MKW11':392.37,'MKW8':491.32,'AWM4':476.67,'A2063':727.06,'A2052':626.32,'NGC6107':616.86,'Coma':937.03,'A1367':794.61,'Hercules':772.74}


###### Cut for figure 8 & 9
clusterDFCCut={'MKW11':1.2, 'MKW8':.9, 'AWM4':1.25, 'A2063':.8, 'A2052':.75, 'NGC6107':1.1, 'Coma':1, 'A1367':.9, 'Hercules':1}


class cluster:
   def __init__(self,clustername):
       infile='/home/rfinn/research/LocalClusters/MasterTables/'+clustername+'mastertable.fits'
       #infile='/home/alissa/LocalClusters/'+clustername+'mastertable.fits'
       tb=pyfits.open(infile)
       tbdata=tb[1].data
       tb.close()
       self.agcflag=tbdata.field('AGCflag')
       self.sdssflag=tbdata.field('SDSSflag')
       self.sexsdssflag=tbdata.field('SEXSDSSflag')
       self.sex24flag=tbdata.field('SEX24FLAG')
       self.agcnumber=tbdata.field('AGCNUMBER')
       self.ra=tbdata.field('AGC-RA')
       self.dec=tbdata.field('AGC-DEC')
       self.a100=tbdata.field('A100')
       self.b100=tbdata.field('B100')
       self.mag10=tbdata.field('MAG10')
       self.posang=tbdata.field('POSANG')
       self.bsteintype=tbdata.field('BSTEINTYPE')
       self.vopt=tbdata.field('VOPT')#optical velocity
       self.verr=tbdata.field('VERR')
       self.vsource=tbdata.field('VSOURCE')#HI velocity
       self.flux100=tbdata.field('FLUX100')
       self.rms100=tbdata.field('RMS100')
       self.v21=tbdata.field('V21')
       self.width=tbdata.field('WIDTH')
       self.widtherr=tbdata.field('WIDTHERR')
       self.sdssra=tbdata.field('SDSSRA')
       self.sdssdec=tbdata.field('SDSSDEC')
       self.sdssu=tbdata.field('SDSSU')
       self.sdssg=tbdata.field('SDSSG')
       self.sdssr=tbdata.field('SDSSR')
       self.sdssi=tbdata.field('SDSSI')
       self.sdssz=tbdata.field('SDSSZ')
       self.sdssspecz=tbdata.field('SDSSSPECZ')
       self.sdsshaew=tbdata.field('SDSSHAEW')
       self.sdsshaewerr=tbdata.field('SDSSHAEWERR')
       self.numberser=tbdata.field('NUMBERSER')
       self.ximageser=tbdata.field('XIMAGESER')
       self.yimageser=tbdata.field('YIMAGESER')
       self.xminimageser=tbdata.field('XMINIMAGESER')
       self.xmaximageser=tbdata.field('XMAXIMAGESER')
       self.yminimageser=tbdata.field('YMINIMAGESER')
       self.raser=tbdata.field('RASER')
       self.decser=tbdata.field('DECSER')
       self.fluxisoser=tbdata.field('FLUXISOSER')
       self.fluxerrisoser=tbdata.field('FLUXERRISOSER')
       self.magisoser=tbdata.field('MAGISOSER')
       self.magerrisoser=tbdata.field('MAGERRISOSER')
       self.fluxautoser=tbdata.field('FLUXAUTOSER')
       self.fluxerrautoser=tbdata.field('FLUXERRAUTOSER')
       self.magautoser=tbdata.field('MAGAUTOSER')
       self.magerrautoser=tbdata.field('MAGERRAUTOSER')
       self.fluxpetroser=tbdata.field('FLUXPETROSER')
       self.fluxerrpetroser=tbdata.field('FLUXERRPETROSER')
       self.magpetroser=tbdata.field('MAGPETROSER')
       self.magerrpetroser=tbdata.field('MAGERRPETROSER')
       self.kronradser=tbdata.field('KRONRADSER')#kron radius
       self.petroradser=tbdata.field('PETRORADSER')#petrosian radius
       self.fluxradser=tbdata.field('FLUXRADSER')#1/2 light radius
       self.isoareaser=tbdata.field('ISOAREASER')
       self.aworldser=tbdata.field('AWORLDSER')
       self.bworldser=tbdata.field('BWORLDSER')
       self.thetaser=tbdata.field('THETASER')
       self.errthetaser=tbdata.field('ERRTHETASER')
       self.thetaj2000ser=tbdata.field('THETAJ2000SER')
       self.errthetaj2000ser=tbdata.field('ERRTHETAJ2000SER')
       self.elongser=tbdata.field('ELONGATIONSER')
       self.elliptser=tbdata.field('ELLIPTICITYSER')
       self.fwhmser=tbdata.field('FWHMSER')
       self.flagsser=tbdata.field('FLAGSSER')
       self.classstarser=tbdata.field('CLASSSTARSER')
       self.numberse24=tbdata.field('NUMBERSE24')
       self.ximagese24=tbdata.field('XIMAGESE24')
       self.yimagese24=tbdata.field('YIMAGESE24')
       self.xminimagese24=tbdata.field('XMINIMAGESE24')
       self.xmaximagese24=tbdata.field('XMAXIMAGESE24')
       self.xminimagese24=tbdata.field('YMINIMAGESE24')
       self.rase24=tbdata.field('RASE24')
       self.decse24=tbdata.field('DECSE24')
       self.fluxisose24=tbdata.field('FLUXISOSE24')
       self.fluxerrisose24=tbdata.field('FLUXERRISOSE24')
       self.magisose24=tbdata.field('MAGISOSE24')
       self.magerrisose24=tbdata.field('MAGERRISOSE24')
       self.fluxautose24=tbdata.field('FLUXAUTOSE24')
       self.fluxerrautose24=tbdata.field('FLUXERRAUTOSE24')
       self.magautose24=tbdata.field('MAGAUTOSE24')
       self.magerrautose24=tbdata.field('MAGERRAUTOSE24')
       self.fluxpetrose24=tbdata.field('FLUXPETROSE24')
       self.fluxerrpetrose24=tbdata.field('FLUXERRPETROSE24')
       self.magpetrose24=tbdata.field('MAGPETROSE24')
       self.magerrpetrose24=tbdata.field('MAGERRPETROSE24')
       self.kronradse24=tbdata.field('KRONRADSE24')
       self.petroradse24=tbdata.field('PETRORADSE24')
       self.fluxradse24=tbdata.field('FLUXRADSE24')
       self.isoarease24=tbdata.field('ISOAREASE24')
       self.aworldse24=tbdata.field('AWORLDSE24')
       self.bworldse24=tbdata.field('BWORLDSE24')
       self.thetase24=tbdata.field('THETASE24')
       self.errthetase24=tbdata.field('ERRTHETASE24')
       self.thetaj2000se24=tbdata.field('THETAJ2000SE24')
       self.errthetaj2000se24=tbdata.field('ERRTHETAJ2000SE24')
       self.elongse24=tbdata.field('ELONGATIONSE24')
       self.elliptse24=tbdata.field('ELLIPTICITYSE24')
       self.fwhmse24=tbdata.field('FWHMSE24')
       self.flagsse24=tbdata.field('FLAGSSE24')
       self.classstarse24=tbdata.field('CLASSSTARSE24')
       self.f24dist=self.fluxautose24[self.sex24flag]
       #used to plot the velocity vs. distance from center
       self.clusterra=clusterra[clustername]
       self.clusterdec=clusterdec[clustername]
       self.alldfc=sqrt(((self.dec-self.clusterdec)**2)+((self.ra-self.clusterra)**2))
       self.allvelocity=(3e5)*self.sdssspecz
       self.clusterV=clusterV[clustername]
       self.clustercentralvelocity=clustercentralvelocity[clustername]
      # self.clusterMinV=clusterlow8[clustername]
      # self.clusterMaxV=clusterhigh8[clustername]
       self.clusterMinV=self.clusterV-(3*clustersbi[clustername])
       self.clusterMaxV=self.clusterV+(3*clustersbi[clustername])


       self.centerV=arange(len(self.alldfc),dtype=float)
       self.centerMinV=arange(len(self.alldfc),dtype=float)
       self.centerMaxV=arange(len(self.alldfc),dtype=float)
       self.dfcCut=clusterDFCCut[clustername]
       p=0
       while p<len(self.alldfc):
           self.centerV[p]=self.clusterV
           self.centerMinV[p]=self.clusterMinV
           self.centerMaxV[p]=self.clusterMaxV
           p=p+1
       self.r=size(self.vopt)
       a=0
       while a<self.r:
           if self.sdssflag[a]<=0 and self.v21[a]>0:
               self.allvelocity[a]=self.v21[a]
               a=a+1
           elif self.sdssflag[a]<=0 and self.v21[a]<=0:
               self.allvelocity[a]=self.vopt[a]
              # print(self.agcnumber[a])
               a=a+1
           else:
               a=a+1
#########used to make distance cut color coding for figure 5
       self.raCut=[]
       self.raCut=array(self.raCut,'f')
       self.decCut=[]
       self.decCut=array(self.decCut,'f')
       self.raCut1=[]
       self.raCut1=array(self.raCut1,'f')
       self.decCut1=[]
       self.decCut1=array(self.decCut1,'f')
       self.raCut2=[]
       self.raCut2=array(self.raCut2,'f')
       self.decCut2=[]
       self.decCut2=array(self.decCut2,'f')
       cut=0
       while cut<len(self.alldfc):
           if self.alldfc[cut]<=1: #change this value to change cut
               self.raCut=append(self.raCut,self.ra[cut])
               self.decCut=append(self.decCut,self.dec[cut])
               cut=cut+1
           elif self.alldfc[cut]<=2:
               self.raCut1=append(self.raCut1,self.ra[cut])
               self.decCut1=append(self.decCut1,self.dec[cut])
               cut=cut+1
           elif self.alldfc[cut]<=3:
               self.raCut2=append(self.raCut2,self.ra[cut])
               self.decCut2=append(self.decCut2,self.dec[cut])
               cut=cut+1
           else:
               cut=cut+1
#########used to make velocity cut color coding for figure 7
       self.raV=[]
       self.raV=array(self.raV,'f')
       self.decV=[]
       self.decV=array(self.decV,'f')
       self.raV1=[]
       self.raV1=array(self.raV1,'f')
       self.decV1=[]
       self.decV1=array(self.decV1,'f')
       self.raV2=[]
       self.raV2=array(self.raV2,'f')
       self.decV2=[]
       self.decV2=array(self.decV2,'f')
       self.dfcV=[]
       self.dfcV=array(self.dfcV,'f')
       cut=0
       while cut<len(self.alldfc):
           #diff=self.allvelocity[cut]-self.clustercentralvelocity
           diff=self.allvelocity[cut]
           if diff<=self.clusterMaxV and diff>=self.clusterMinV: #change this value to change cut
               self.raV=append(self.raV,self.ra[cut])
               self.decV=append(self.decV,self.dec[cut])
               self.dfcV=append(self.dfcV,self.alldfc[cut])
               cut=cut+1
           elif diff>self.clusterMaxV:
               self.raV1=append(self.raV1,self.ra[cut])
               self.decV1=append(self.decV1,self.dec[cut])
               cut=cut+1
           elif diff<self.clusterMinV:
               self.raV2=append(self.raV2,self.ra[cut])
               self.decV2=append(self.decV2,self.dec[cut])
               cut=cut+1
           else:
               cut=cut+1

##########used for figure 8- velocity cut and dfc cut
       self.raCutV=[]
       self.raCutV=array(self.raCutV,'f')
       self.decCutV=[]
       self.decCutV=array(self.decCutV,'f')
       cut=0
       while cut<len(self.dfcV):
           if self.dfcV[cut]<=self.dfcCut: #change this value to change cut
               self.raCutV=append(self.raCutV,self.raV[cut])
               self.decCutV=append(self.decCutV,self.decV[cut])
               cut=cut+1
           else:
               cut=cut+1


       #outfile=clustername+'f24.dat'
       #out1=open(outfile,'w')
       #for i in range(len(self.f24dist)):
       #    s='%3.2e \n'%(self.f24dist[i])
       #    out1.write(s)
       #out1.close()

       #make plot of RA and Dec
  # def plotpositions(self):
   #    figure()
   #    clf()
   #    plot(self.ra,self.dec,'k.')
   #    xlabel('RA (deg)')
   #    ylabel('Dec (deg)')
       #make plot of
  # def plotlf(self):
      # figure(1)

       #y=hist(self.f24dist,histtype='step')
       #ngal=y[0]
       #x=y[1]
       #xbin=zeros(len(ngal))
       #for i in range(len(xbin)):
        #   xbin[i]=0.5*(x[i]+x[i+1])
       #clf()
       #self.xbin=xbin
       #self.ngal=ngal
       #figure(2)
       #plot(xbin,ngal,'ro')
       #errorbar(xbin,ngal,sqrt(ngal))
       #ax=gca()
       #ax.set_yscale('log')
       #ax.set_xscale('log')
       #xlabel('24um Flux')


mkw11=cluster('MKW11')
coma=cluster('Coma')
herc=cluster('Hercules')
awm4=cluster('AWM4')
a1367=cluster('A1367')
a2052=cluster('A2052')
a2063=cluster('A2063')
ngc=cluster('NGC6107')
mkw8=cluster('MKW8')
clusters=[mkw11,mkw8,awm4,a2063,a2052,ngc,coma,a1367,herc]
clusternames=['MKW11','MKW8','AWM4','Abell2063','Abell2052','NGC6107','Coma','Abell1367','Hercules']

#RA vs. Dec plots using subplot
figure(4,figsize=(16,14))
clf()
subplots_adjust(left=0.1, right=.9, bottom=.1, wspace=.27, hspace=.22)
#clf()
f=1
while f<10:
    subplot(3,3,f)
    xlabel('RA (deg)')
    ylabel('Dec (deg)')
    plot(clusters[f-1].ra,clusters[f-1].dec,'k.')
    xc=clusters[f-1].clusterra
    yc=clusters[f-1].clusterdec
    plot([xc],[yc],'go')
    gca().add_patch(Circle((xc,yc),radius=.5,alpha=1,ec='m', lw=1.5, fill=False))
    gca().add_patch(Circle((xc,yc),radius=1,alpha=1,ec='b', lw=1.5, fill=False))
    gca().add_patch(Circle((xc,yc),radius=1.5,alpha=1,ec='c', lw=1.5, fill=False))
    gca().add_patch(Circle((xc,yc),radius=2,alpha=1,ec='g', lw=1.5, fill=False))
    gca().add_patch(Circle((xc,yc),radius=2.5,alpha=1,ec='y', lw=1.5, fill=False))
    gca().add_patch(Circle((xc,yc),radius=3,alpha=1,ec='r', lw=1.5, fill=False))
    title(clusternames[f-1])
    f=f+1
subplot(3,3,4)
plot([a2052.clusterra],[a2052.clusterdec],'mo')
subplot(3,3,5)
plot([a2063.clusterra],[a2063.clusterdec],'mo')
show()

#RA vs. Dec plots with cutoff
figure(5,figsize=(16,14))
clf()
subplots_adjust(left=0.1, right=.9, bottom=.1, wspace=.27, hspace=.22)
#clf()
f=1
while f<10:
    subplot(3,3,f)
    xlabel('RA (deg)')
    ylabel('Dec (deg)')
    plot(clusters[f-1].ra,clusters[f-1].dec,'k.')
    plot(clusters[f-1].raCut,clusters[f-1].decCut,'b.')
    plot(clusters[f-1].raCut1,clusters[f-1].decCut1,'m.')
    plot(clusters[f-1].raCut2,clusters[f-1].decCut2,'g.')
    xc=clusters[f-1].clusterra
    yc=clusters[f-1].clusterdec
    plot([xc],[yc],'go')
    #gca().add_patch(Circle((xc,yc),radius=1,alpha=1,ec='m', lw=1.5, fill=False))
    title(clusternames[f-1])
    f=f+1
subplot(3,3,4)
plot([a2052.clusterra],[a2052.clusterdec],'mo')
subplot(3,3,5)
plot([a2063.clusterra],[a2063.clusterdec],'mo')
show()




#velocity vs. Distance from center of all the clusters using subplot(vopt/sdssspecz)
figure(6,figsize=(16,14))
clf()
subplots_adjust(left=0.1, right=.9, bottom=.1, wspace=.27, hspace=.22)
#clf()
d=1
while d<10:
    subplot(3,3,d)
    xlabel('Distance from Center')
    ylabel('Velocity')
    plot(clusters[d-1].alldfc,clusters[d-1].allvelocity,'g.')
    plot(clusters[d-1].alldfc,clusters[d-1].centerV,'k-')
    plot(clusters[d-1].alldfc,clusters[d-1].centerMinV,'b-')
    plot(clusters[d-1].alldfc,clusters[d-1].centerMaxV,'b-')
    title(clusternames[d-1])
    d=d+1
show()



#RA vs. Dec colorcoding by velocity
figure(7,figsize=(16,14))
clf()
subplots_adjust(left=0.1, right=.9, bottom=.1, wspace=.27, hspace=.22)
f=1
while f<10:
    subplot(3,3,f)
    xlabel('RA (deg)')
    ylabel('Dec (deg)')
    plot(clusters[f-1].ra,clusters[f-1].dec,'k.')
    plot(clusters[f-1].raV,clusters[f-1].decV,'b.')
    plot(clusters[f-1].raV1,clusters[f-1].decV1,'m.')
    plot(clusters[f-1].raV2,clusters[f-1].decV2,'g.')
    xc=clusters[f-1].clusterra
    yc=clusters[f-1].clusterdec
    plot([xc],[yc],'go')
    #gca().add_patch(Circle((xc,yc),radius=1,alpha=1,ec='m', lw=1.5, fill=False))
    title(clusternames[f-1])
    f=f+1
subplot(3,3,4)
plot([a2052.clusterra],[a2052.clusterdec],'mo')
subplot(3,3,5)
plot([a2063.clusterra],[a2063.clusterdec],'mo')
show()

#RA vs. Dec removing galaxies with high/low velocities based on calculations from dispersion.py
figure(8,figsize=(16,14))
clf()
subplots_adjust(left=0.1, right=.9, bottom=.1, wspace=.27, hspace=.22)
f=1
while f<10:
    subplot(3,3,f)
    xlabel('RA (deg)')
    ylabel('Dec (deg)')
    plot(clusters[f-1].raV,clusters[f-1].decV,'k.')
    #plot(clusters[f-1].raCutV,clusters[f-1].decCutV,'m.')
    xc=clusters[f-1].clusterra
    yc=clusters[f-1].clusterdec
    plot([xc],[yc],'go')
    #gca().add_patch(Circle((xc,yc),radius=.5,alpha=1,ec='g', lw=1.5, fill=False))
    #gca().add_patch(Circle((xc,yc),radius=1,alpha=1,ec='b', lw=1.5, fill=False))
    #gca().add_patch(Circle((xc,yc),radius=1.5,alpha=1,ec='c', lw=1.5, fill=False))
    #gca().add_patch(Circle((xc,yc),radius=2,alpha=1,ec='g', lw=1.5, fill=False))
    #gca().add_patch(Circle((xc,yc),radius=2.5,alpha=1,ec='y', lw=1.5, fill=False))
    #gca().add_patch(Circle((xc,yc),radius=3,alpha=1,ec='r', lw=1.5, fill=False))
    title(clusternames[f-1])
    f=f+1
subplot(3,3,4)
plot([a2052.clusterra],[a2052.clusterdec],'bo')
subplot(3,3,5)
plot([a2063.clusterra],[a2063.clusterdec],'bo')
show()

#Same as Figure 8 but w/out black dots for cut off galaxies
figure(9,figsize=(16,14))
clf()
subplots_adjust(left=0.1, right=.9, bottom=.1, wspace=.27, hspace=.22)
f=1
while f<10:
    subplot(3,3,f)
    xlabel('RA (deg)')
    ylabel('Dec (deg)')
    plot(clusters[f-1].raCutV,clusters[f-1].decCutV,'m.')
    xc=clusters[f-1].clusterra
    yc=clusters[f-1].clusterdec
    plot([xc],[yc],'go')
    title(clusternames[f-1])
    f=f+1
show()
