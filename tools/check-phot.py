#!/usr/bin/env python
"""
checking photometry from flatfielding with object masking with flatfielding
w/out object masking
"""

import sys, os
import Numeric as N
import scipy
from math import *
import ppgplot

def xysimple(x,y,xlabel,ylabel):
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0)
    ppgplot.pglab(xlabel,ylabel,"")
    ppgplot.pgpt(x,y,3)
def compdiff(x,y,xlabel,ylabel):
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ave=N.average(N.compress((x < 21) & (im1.sn > 3),y))
    std=scipy.stats.std(N.compress((x < 21) & (im1.sn > 3),y))
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,-.8,.8,0)
    ppgplot.pglab(xlabel,ylabel,"")
    ppgplot.pgpt(x,y,3)
    x=N.arange((int(xmin)-1),(int(xmax)+2),1)
    y=ave*x/x
    ppgplot.pgsci(2)
    ppgplot.pgline(x,y)
    ppgplot.pgsci(1)
    y1=y-std
    ppgplot.pgline(x,y1)
    y1=y+std
    ppgplot.pgline(x,y1)

def plotoutlierpos():
    dm=0.01
    xlabel="x position of outliers"
    ylabel="y position of outliers"
    delta=abs(im1.magap5 - im2.magap5)
    x=N.compress(delta > dm,im1.x)
    y=N.compress(delta > dm,im1.y)
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    #ppgplot.pgenv(xmin,xmax,ymin,ymax,0)
    ppgplot.pgenv(0,920,0,875,0)
    ppgplot.pglab(xlabel,ylabel,"")
    ppgplot.pgpt(x,y,3)
    x=N.compress(delta > dm,im2.x)
    y=N.compress(delta > dm,im2.y)
    diff=N.compress(delta > dm,delta)
    ppgplot.pgpt(x,y,4)
    for i in range(len(x)):
        print x[i], y[i],diff[i]
def checkposition(x,y):
    pos=1
    if x < 10:
        pos = 0
    if (y > 500):
	d=sqrt((x-430)**2 + (y-423)**2)
	if(d > 495):
            pos=0
    if (y < 300):
	d=sqrt((x-465)**2 + (y-518)**2)
	if(d > 500):
            pos=0
      #check for bottom right corner
    if ((x > 514) & (y < 384)):
	d=sqrt((x-514)**2 + (y-384)**2)
	if(d > 400):
            pos=0
      #bottom edge
    if (y < 10):
	pos=0
      #reject stars
    d=sqrt((x-498.5)**2 + (y-43.3)**2)
    if (d < 4.):
	pos=0
    d=sqrt((x-100.7)**2 + (y-91.5)**2)
    if (d < 6.):
        pos=0
    d=sqrt((x-716)**2 + (y-125)**2)
    if (d < 6.):
        pos=0
    d=sqrt((x-186.9)**2 + (y-134.5)**2)
    if (d < 4.):
	pos=0
    d=sqrt((x-299)**2 + (y-141)**2)
    if (d < 6.):
        pos=0
    d=sqrt((x-69)**2 + (y-147)**2)
    if (d < 4.):
	pos=0
    d=sqrt((x-389)**2 + (y-219.6)**2)
    if (d < 4.):
	pos=0
    d=sqrt((x-520.9)**2 + (y-235)**2)
    if (d < 4.):
	pos=0
    d=sqrt((x-745.7)**2 + (y-311.6)**2)
    if (d < 4.):
        pos=0
    d=sqrt((x-704.)**2 + (y-342)**2)
    if (d < 6.):
	pos=0
    d=sqrt((x-248.2)**2 + (y-422.5)**2)
    if (d < 6.):
	pos=0
    d=sqrt((x-441.8)**2 + (y-484)**2)
    if (d < 4.):
	pos=0
    d=sqrt((x-292.)**2 + (y-574)**2)
    if (d < 6.):
        pos=0
    d=sqrt((x-879)**2 + (y-598)**2)
    if (d < 4.):
	pos=0
    d=sqrt((x-203.6)**2 + (y-618.7)**2)
    if (d < 6.):
	pos=0
    d=sqrt((x-349)**2 + (y-620)**2)
    if (d < 6.):
	pos=0
    d=sqrt((x-469.6)**2 + (y-676)**2)
    if (d < 6.):
	pos=0
    d=sqrt((x-676.4)**2 + (y-659.7)**2)
    if (d < 4.):
	pos=0
    d=sqrt((x-766.5)**2 + (y-704)**2)
    if (d < 4.):
	pos=0      
    d=sqrt((x-510.5)**2 + (y-828.7)**2)
    if (d < 4.):
	pos=0
    d=sqrt((x-323.19)**2 + (y-634.62)**2)
    if (d < 4.):
	pos=0
    d=sqrt((x-297.)**2 + (y-84.)**2)
    if (d < 4.):
	pos=0
    d=sqrt((x-827.)**2 + (y-653.)**2)
    if (d < 4.):
	pos=0
    #ghost image
    d=sqrt((x-704.)**2 + (y-286.)**2)
    if (d < 4.):
	pos=0
    d=sqrt((x-715.)**2 + (y-68.)**2)
    if (d < 4.):
	pos=0
    return pos
      

class Catalog:
    def __init__(self):
        self.number = []
        self.magiso = []
        self.errmagiso = []
        self.magisocor = []
        self.errmagisocor = []
        self.magauto = []
        self.errmagauto = []
        self.isoarea= []
        self.fluxiso = []
        self.errfluxiso = []
        self.x = []
        self.y = []
        self.magbest = []
        self.errmagbest = []
        self.elong = []
        self.ellip = []
        self.fwhm = []
        self.flags = []
        self.classstar = []
        self.fluxap1 = []
        self.errfluxap1 = []
        self.fluxap2 = []
        self.errfluxap2 = []
        self.fluxap3 = []
        self.errfluxap3 = []
        self.fluxap4 = []
        self.errfluxap4 = []
        self.fluxap5 = []
        self.errfluxap5 = []
        self.fluxap6 = []
        self.errfluxap6 = []
        self.fluxap7 = []
        self.errfluxap7 = []
        self.magap1 = []
        self.errmagap1 = []
        self.magap2 = []
        self.errmagap2 = []
        self.magap3 = []
        self.errmagap3 = []
        self.magap4 = []
        self.errmagap4 = []
        self.magap5 = []
        self.errmagap5 = []
        self.magap6 = []
        self.errmagap6 = []
        self.magap7 = []
        self.errmagap7 = []

    def readcat(self):
        i = 0
        for line in open("test.cat",'r'):
            if line.find('#') > -1:
                continue
            fields=line.split()
            sum=0
            x=float(fields[10])
            y=float(fields[11])
            #if i < 10:
            #    print i,x,y,fields[10],fields[11],fields[0]
            #    i = i + 1
            pos=checkposition(x,y)#reject vignetted area and stars
            if pos < 1:
                continue
            self.number.append(int(fields[0]))
            self.magiso.append(float(fields[1]))
            self.errmagiso.append(float(fields[2]))
            self.magisocor.append(float(fields[3]))
            self.errmagisocor.append(float(fields[4]))
            self.magauto.append(float(fields[5]))
            self.errmagauto.append(float(fields[6]))
            self.isoarea.append(float(fields[7]))
            self.fluxiso.append(float(fields[8]))
            self.errfluxiso.append(float(fields[9]))
            self.x.append(float(fields[10]))
            self.y.append(float(fields[11]))
            self.magbest.append(float(fields[12]))
            self.errmagbest.append(float(fields[13]))
            self.elong.append(float(fields[14]))
            self.ellip.append(float(fields[15]))
            self.fwhm.append(float(fields[16]))
            self.flags.append(float(fields[17]))
            self.classstar.append(float(fields[18]))
            self.fluxap1.append(float(fields[19]))
            self.fluxap2.append(float(fields[20]))
            self.fluxap3.append(float(fields[21]))
            self.fluxap4.append(float(fields[22]))
            self.fluxap5.append(float(fields[23]))
            self.fluxap6.append(float(fields[24]))
            self.fluxap6.append(float(fields[25]))
            self.errfluxap1.append(float(fields[26]))
            self.errfluxap2.append(float(fields[27]))
            self.errfluxap3.append(float(fields[28]))
            self.errfluxap4.append(float(fields[29]))
            self.errfluxap5.append(float(fields[30]))
            self.errfluxap6.append(float(fields[31]))
            self.errfluxap7.append(float(fields[32]))
            self.magap1.append(float(fields[33]))
            self.magap2.append(float(fields[34]))
            self.magap3.append(float(fields[35]))
            self.magap4.append(float(fields[36])) 
            self.magap5.append(float(fields[37]))
            self.magap6.append(float(fields[38]))
            self.magap7.append(float(fields[39]))
            self.errmagap1.append(float(fields[40]))
            self.errmagap2.append(float(fields[41]))
            self.errmagap3.append(float(fields[42]))
            self.errmagap4.append(float(fields[43]))
            self.errmagap5.append(float(fields[44]))
            self.errmagap6.append(float(fields[45]))
            self.errmagap7.append(float(fields[46]))
            #def convarray(self): #convert to Numeric arrays
        self.number = N.array(self.number,'i')
        self.magiso = N.array(self.magiso,'f')
        self.errmagiso = N.array(self.errmagiso,'f')
        self.magisocor = N.array(self.magisocor,'f')
        self.errmagisocor = N.array(self.errmagisocor,'f')
        self.magauto = N.array(self.magauto,'f')
        self.errmagauto = N.array(self.errmagauto,'f')
        self.isoarea = N.array(self.isoarea,'f')
        self.fluxiso = N.array(self.fluxiso,'f')
        self.errfluxiso = N.array(self.errfluxiso,'f')
        self.x = N.array(self.x,'f')
        self.y = N.array(self.y,'f')
        self.magbest = N.array(self.magbest,'f')
        self.errmagbest = N.array(self.errmagbest,'f')
        self.elong = N.array(self.elong,'f')
        self.ellip = N.array(self.ellip,'f')
        self.fwhm = N.array(self.fwhm,'f')
        self.flags = N.array(self.flags,'f')
        self.classstar = N.array(self.classstar,'f')
        self.fluxap1 = N.array(self.fluxap1,'f')
        self.errfluxap1 = N.array(self.errfluxap1,'f')
        self.fluxap2 = N.array(self.fluxap2,'f')
        self.errfluxap2 = N.array(self.errfluxap2,'f')
        self.fluxap3 = N.array(self.fluxap3,'f')
        self.errfluxap3 = N.array(self.errfluxap3,'f')
        self.fluxap4 = N.array(self.fluxap4,'f')
        self.errfluxap4 = N.array(self.errfluxap4,'f')
        self.fluxap5 = N.array(self.fluxap5,'f')
        self.errfluxap5 = N.array(self.errfluxap5,'f')
        self.fluxap6 = N.array(self.fluxap6,'f')
        self.errfluxap6 = N.array(self.errfluxap6,'f')
        self.fluxap7 = N.array(self.fluxap7,'f')
        self.errfluxap7 = N.array(self.errfluxap7,'f')
        self.magap1 = N.array(self.magap1,'f')
        self.errmagap1 = N.array(self.errmagap1,'f')
        self.magap2 = N.array(self.magap2,'f')
        self.errmagap2 = N.array(self.magap2,'f')
        self.magap3 = N.array(self.magap3,'f')
        self.errmagap3 = N.array(self.errmagap3,'f')
        self.magap4 = N.array(self.magap4,'f')
        self.errmagap4 = N.array(self.errmagap4,'f')
        self.magap5 = N.array(self.magap5,'f')
        self.errmagap5 = N.array(self.errmagap5,'f')
        self.magap6 = N.array(self.magap6,'f')
        self.errmagap6 = N.array(self.errmagap6,'f')
        self.magap7 = N.array(self.magap7,'f')
        self.errmagap7 = N.array(self.errmagap7,'f')
        self.sn = abs(self.fluxiso/self.errfluxiso)
#os.system("sex cnc105jalign2.fits \n")
#os.system("sex j+n.fits,cnc105jalign2.fits \n")
#os.system("sex j+n.fits,ncc105jalign3.fits \n")
os.system("sex ncc105jalign3.fits \n")
#os.system("sex c1054j2cna.fits \n")
im1=Catalog()
im1.readcat()
#os.system("sex cnc105jalign2.fits,c1054j2cna.fits \n")
#os.system("sex j+n.fits,c1054j2cna.fits \n")
#os.system("sex c1054j2cna.fits,gediscsj.fits -MAG_ZEROPOINT 23.88\n")
#os.system("sex ncc105jalign3.fits,gediscsj.fits -MAG_ZEROPOINT 22.84\n")
os.system("sex ncc105jalign3.fits,ngdim20_mp.fits\n")
im2=Catalog()
im2.readcat()


ppgplot.pgbeg("all.ps/vcps",2,2)
ppgplot.pgsch(2.) #font size
ppgplot.pgslw(4)  #line width

xysimple(im1.magauto,im2.magauto,"magauto mask","magauto")
ppgplot.pgpage
y=im1.magauto - im2.magauto
ave=N.average(y)
std=scipy.stats.std(y)
print "Ave diff in mag auto = ",ave,"+/-",std 
xysimple(im1.magauto,y,"magauto mask","magauto mask - magauto")
ppgplot.pgpage
y=im1.magiso - im2.magiso
ave=N.average(y)
std=scipy.stats.std(y)
print "Ave diff in iso mags = ",ave,"+/-",std 
xysimple(im1.magiso,y,"magiso mask","magiso mask - magiso")
ppgplot.pgpage
y=im1.sn-im2.sn
x=N.compress(im1.magiso > 15,im1.magiso)
y=N.compress(im1.magiso > 15,y)
ave=N.average(y)
std=scipy.stats.std(y)
print "Ave diff in sn = ",ave,"+/-",std 

xysimple(x,y,"magiso mask","sn mask - sn")
ppgplot.pgpage
x=im1.magap3
y=im1.magap3 - im2.magap3
compdiff(x,y,"magap3 w/mask","magap3 w/m - magap5")
ppgplot.pgpage
x=im1.magap4
y=im1.magap4 - im2.magap4
compdiff(x,y,"magap4 w/mask","magap4 w/m - magap5")
ppgplot.pgpage
x=im1.magap5
y=im1.magap5 - im2.magap5
compdiff(x,y,"magap5 w/mask","magap5 w/m - magap5")
ppgplot.pgpage
x=im1.magap6
y=im1.magap6 - im2.magap6
compdiff(x,y,"magap6 w/mask","magap6 w/m - magap5")
ppgplot.pgpage
x=im1.magap7
y=im1.magap7 - im2.magap7
compdiff(x,y,"magap7 w/mask","magap7 w/m - magap5")
ppgplot.pgpage
plotoutlierpos()
#xysimple(im1.fluxiso,im2.fluxiso,"fluxiso mask","fluxiso")
y=im1.magap6 - im2.magap6
ave=N.average(y)
std=scipy.stats.std(y)
print "Ave diff in mag ap 6 = ",ave,"+/-",std 
    
ppgplot.pgend()
