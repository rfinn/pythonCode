#!/usr/bin/env python
import pyfits
from pylab import *

lowerLimit={'MKW11':-10, 'MKW8':-6, 'AWM4':-5, 'A2063':-10, 'A2052':-6, 'NGC6107':-10, 'Coma':-6, 'A1367':-5, 'Hercules':-8}
#read in data
class cluster:
   def __init__(self,clustername):
       infile='/home/rfinn/LocalClusters/MasterTables/'+clustername+'mastertable.fits'
       #infile='/home/alissa/LocalClusters/'+clustername+'mastertable.fits'
       tb=pyfits.open(infile)
       tbdata=tb[1].data
       tb.close()
       self.agcflag=tbdata.field('AGCflag')
       self.sdssflag=tbdata.field('SDSSflag')
       self.sexsdssflag=tbdata.field('SEXSDSSflag')
       self.sex24flag=tbdata.field('SEX24FLAG')
       self.mpaflag=tbdata.field('MPAFLAG')
       self.agcnumber=tbdata.field('AGCNUMBER')
       self.ra=tbdata.field('AGCRA')
       self.dec=tbdata.field('AGCDEC')
       self.hA=tbdata.field('MPAHALPHA')
       self.hB=tbdata.field('MPAHBETA')
       self.oIII=tbdata.field('MPAOIII')
       self.nII=tbdata.field('MPANII')
       count=0
       self.xaxis=array([],'f')
       self.yaxis=array([],'f')
       self.agcnumbers=array([],'f')
       while count<len(self.agcnumber):
           if self.sdssflag[count]==True and self.mpaflag[count]==True:
               self.xaxis=append(self.xaxis,log(self.nII[count]/self.hA[count]))
               self.yaxis=append(self.yaxis,log(self.oIII[count]/self.hB[count]))
               self.agcnumbers=append(self.agcnumbers,self.agcnumber[count])
               count=count+1
           else:
               count=count+1
       #create array of x values to plot line
       self.lowerLimit=lowerLimit[clustername]
       self.xsort=arange(self.lowerLimit,12,.1)
       self.xsortK=arange(self.lowerLimit,10,.1)
       self.xsortS=arange(self.lowerLimit,12,.1)

       #equations for lines
       self.line=(.61/(self.xsort-.05)+1.3)#Kauffman 2003?
       self.lineK=(.61/(self.xsort-.47)+1.19)#Kewely
       self.lineS=((-30.787+(1.1358*self.xsort)+((.27297)*(self.xsort)**2))*tanh(5.7409*self.xsort))-31.093 #Stasinska 2006

       #cut used to remove outliers and top half of hyperbola to make plots easier to view
       self.xaxisCut=array([],'f')
       self.lineCut=array([],'f')
       cc=0
       while cc<len(self.xsort):
           if self.line[cc]<1.25 and self.line[cc]>-4:
               self.xaxisCut=append(self.xaxisCut,self.xsort[cc])
               self.lineCut=append(self.lineCut,self.line[cc])
               cc=cc+1
           else:
               cc=cc+1
       self.lineCutK=array([],'f')
       self.xaxisCutK=array([],'f')
       cc=0
       while cc<len(self.xsortK):
           if self.lineK[cc]<1.25 and self.lineK[cc]>-4:
               self.xaxisCutK=append(self.xaxisCutK,self.xsortK[cc])
               self.lineCutK=append(self.lineCutK,self.lineK[cc])
               cc=cc+1
           else:
               cc=cc+1
       self.lineCutS=array([],'f')
       self.xaxisCutS=array([],'f')
       cc=0
       while cc<len(self.xsortS):
           if self.lineS[cc]<1.25 and self.lineS[cc]>-4:
               self.xaxisCutS=append(self.xaxisCutS,self.xsortS[cc])
               self.lineCutS=append(self.lineCutS,self.lineS[cc])
               cc=cc+1
           else:
               cc=cc+1

       #makes arrays with points above Stasinska line so that they can be plotted seperatly
       self.xagn=array([],'f')
       self.yagn=array([],'f')
       self.agcAgnList=array([],'f')
       countAgn=0
       while countAgn<len(self.xaxis):
           if self.yaxis[countAgn]>((-30.787+(1.1358*self.xaxis[countAgn])+((.27297)*(self.xaxis[countAgn])**2))*tanh(5.7409*self.xaxis[countAgn]))-31.093:#this equation can be edited to match a different reference
               self.xagn=append(self.xagn, self.xaxis[countAgn])
               self.yagn=append(self.yagn, self.yaxis[countAgn])
               self.agcAgnList=append(self.agcAgnList,self.agcnumbers[countAgn])
               countAgn=countAgn+1
           else:
               countAgn=countAgn+1



#set cluster names
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



#Plots to determine agn  using subplot
figure(figsize=(16,14))
#clf()
subplots_adjust(left=0.1, right=.9, bottom=.1, wspace=.27, hspace=.22)
f=1
while f<10:
    subplot(3,3,f)
    v=[-3,1,-2,2]
    xlabel('log([NII]/[Halpha])')
    ylabel('log([OIII]/[Hbeta])')
    print(clusternames[f-1])
    i=0
    while i<len(clusters[f-1].xaxis):
#         try:
        plot(array([clusters[f-1].xaxis[i]],'f'),array([clusters[f-1].yaxis[i]],'f'),'b.')
        i=i+1
#         except OverflowError:
#             print i,clusters[f-1].xaxis[i],clusters[f-1].yaxis[i]
#             i=i+1
    plot(clusters[f-1].xagn,clusters[f-1].yagn,'c.')
    plot(clusters[f-1].xaxisCut,clusters[f-1].lineCut,'g-')
    plot(clusters[f-1].xaxisCutK,clusters[f-1].lineCutK,'k-')
    plot(clusters[f-1].xaxisCutS,clusters[f-1].lineCutS,'r-')
    axis(v)
    title(clusternames[f-1])
    f=f+1
show()



#Plots to determine agn  using subplot
#figure(figsize=(16,14))
#clf()
#subplots_adjust(left=0.1, right=.9, bottom=.1, wspace=.27, hspace=.22)
#f=1
#while f<10:
#     subplot(3,3,f)
#     v=[-3,1,-2,2]
#     xlabel('log([NII]/[Halpha])')
#     ylabel('log([OIII]/[Hbeta])')
#     plot(clusters[f-1].xaxis,clusters[f-1].yaxis,'b.')
#     plot(clusters[f-1].xagn,clusters[f-1].yagn,'c.')
#     plot(clusters[f-1].xaxisCut,clusters[f-1].lineCut,'g-')
#     plot(clusters[f-1].xaxisCutK,clusters[f-1].lineCutK,'k-')
#     plot(clusters[f-1].xaxisCutS,clusters[f-1].lineCutS,'r-')
#     axis(v)
#     title(clusternames[f-1])
#     f=f+1
#show()
