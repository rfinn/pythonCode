#!/usr/bin/env python
import matplotlib 
#matplotlib.rc('font',family='times', style='normal', variant='normal',weight='ultrabold', stretch='normal', size=24)
#matplotlib.rc('text', usetex=True)
import pylab

import numarray as N
#import pylab #use matplotlib instead of scipy
from math import *
import time, os
import ppgplot

#matplotlib.use('Agg')

h=0.7

def psplotinit(output):
    file=output+"/vcps"
    ppgplot.pgbeg(file,1,1)
    ppgplot.pgpap(8.,1.)
    ppgplot.pgsch(1.7) #font size
    ppgplot.pgslw(7)  #line width

class ps:
    def __init__(self):
	print "getting some data!"
    def readpsoutput(self,file,z):
	input=open(file,'r')
	mass=[]
	frac=[]

	for line in input:
	    t=line.split()
	    mass.append(float(t[2]))
	    frac.append(float(t[3]))

	    #try:
		#frac.append(float(t[3]))
	    #except:
	#	frac.append(1.e-10)
	input.close()

	mass=N.array(mass,'d')
	sigma=((mass/(1.2e15/h)*N.sqrt(.7+.3*(1+z)**3))**(1./3.))*1000.
	frac=N.array(frac,'d')
	sigma=N.log10(sigma)


	#maccret=(frac*mass)
	maccret=N.zeros(len(mass),'d')
	for i in range(len(mass)):
	    maccret[i]=mass[i]*frac[i]
	    print i,mass[i],frac[i],maccret[i]
	frac=N.log10(frac)
	self.sigma=sigma
	self.mass=mass
	self.frac=frac
	self.maccret=maccret

print 'psmass_m_lowz_3Gyr_1e13.dat'
lz3lm=ps()
#lz3.readpsoutput('psmass_m_z0_3Gyr.dat',0.)
lz3lm.readpsoutput('psmass_m_lowz_3Gyr_1e13.dat',0.07)

#print lz3lm.maccret
print 'psmass_m_lowz_3Gyr_1e14.dat'
lz3hm=ps()
lz3hm.readpsoutput('psmass_m_lowz_3Gyr_1e14.dat',0.07)

print 'psmass_m_lowz_1Gyr_1e13.dat'
lz1lm=ps()
#lz1.readpsoutput('psmass_m_z0_1Gyr.dat',0.)
lz1lm.readpsoutput('psmass_m_lowz_1Gyr_1e13.dat',0.07)

print 'psmass_m_lowz_1Gyr_1e14.dat'
lz1hm=ps()
lz1hm.readpsoutput('psmass_m_lowz_1Gyr_1e14.dat',0.07)

print 'psmass_m_hiz_3Gyr_1e13.dat'
hz3lm=ps()
#hz3.readpsoutput('psmass_m_3Gyr.dat',0.75)
hz3lm.readpsoutput('psmass_m_hiz_3Gyr_1e13.dat',0.75)

print 'psmass_m_hiz_3Gyr_1e14.dat'
hz3hm=ps()
hz3hm.readpsoutput('psmass_m_hiz_3Gyr_1e14.dat',0.75)

print 'psmass_m_hiz_1Gyr_1e13.dat'
hz1lm=ps()
#hz1.readpsoutput('psmass_m_1Gyr.dat',0.75)
hz1lm.readpsoutput('psmass_m_hiz_1Gyr_1e13.dat',0.75)

print 'psmass_m_hiz_1Gyr_1e14.dat'
hz1hm=ps()
hz1hm.readpsoutput('psmass_m_hiz_1Gyr_1e14.dat',0.75)

#lz05=ps()
#lz05.readpsoutput('psmass_m_z0_0.5Gyr.dat',0.)
#hz05=ps()
#hz05.readpsoutput('psmass_m_0.5Gyr.dat',0.75)
#lz01=ps()
#lz01.readpsoutput('psmass_m_z0_0.1Gyr.dat',0.)
#hz01=ps()
#hz01.readpsoutput('psmass_m_0.1Gyr.dat',0.75)
#pylab.plot(mass,frac,'b')
#pylab.xlabel(r'Cluster Mass (1e14 $M_\odot$)')
#pylab.ylabel('fS(1e11:1e13)')
#pylab.show()

#pylab.plot(sigma,frac,'b')
#pylab.xlabel(r'Velocity Dispersion')
#pylab.ylabel('fS(1e11:1e13)')
#pylab.show()

def mratio():
    r3=hz3.maccret/lz3.maccret
    for i in range(len(r3)):
	print i,lz3.sigma[i],hz3.sigma[i],lz3.mass[i],hz3.mass[i]
	print i,lz01.sigma[i],hz01.sigma[i],lz01.mass[i],hz01.mass[i]
    r1=hz1.maccret/lz1.maccret
    r05=hz05.maccret/lz05.maccret
    ra=N.array(hz01.maccret,'d')
    rb=N.array(lz01.maccret,'d')
    r01=ra/rb
    for i in range(len(r01)):
	print "ratio ",hz01.maccret[i],lz01.maccret[i],ra[i],rb[i],r01[i]
    pylab.plot(lz3.mass,r3,'b-',label="3 Gyr",linewidth=4)
    pylab.plot(lz1.mass,r1,'c--',label="1 Gyr",linewidth=4)
    pylab.plot(lz05.mass,r05,'r-.',label="0.5 Gyr",linewidth=4)
    pylab.plot(lz01.mass,r01,'g:',label="0.3 Gyr",linewidth=4)
    pylab.axis([1.e14,3.e15,2.8,7.])
    ax=pylab.gca()
    ax.set_xscale('log')
    pylab.legend(loc='lower right')
    pylab.xlabel(r"$M_{cl} (M_\cdot)$",fontsize=20)
    pylab.ylabel(r"$M_{acc}(z=0.75)\ / \M_{acc}(z=0.07)$",fontsize=20)
    pylab.savefig('psratio.eps')

def mrationew():
    #for i in range(len(lz1lm.mass)):
#	m=lz1lm.mass[i]
#	l=lz1lm.maccret[i]
#	h=hz1lm.maccret[i]
#	r=h/l
#	print i,m,l,h,r
    #print lz1lm.maccret
    #print hz1lm.maccret
    #print hz3lm.maccret
    r3lm=(hz3lm.maccret)/(lz3lm.maccret)
    r3hm=(hz3hm.maccret)/(lz3hm.maccret)
    #for i in range(len(r3)):
#	print i,lz3.sigma[i],hz3.sigma[i],lz3.mass[i],hz3.mass[i]
#	print i,lz01.sigma[i],hz01.sigma[i],lz01.mass[i],hz01.mass[i]
    r1lm=hz1lm.maccret/lz1lm.maccret
    r1hm=hz1hm.maccret/lz1hm.maccret
    #ra=N.array(hz01.maccret,'d')
    #rb=N.array(lz01.maccret,'d')
    #r01=ra/rb
    #for i in range(len(r01)):
	#print "ratio ",hz01.maccret[i],lz01.maccret[i],ra[i],rb[i],r01[i]
    pylab.plot(lz3lm.mass,r3lm,'b-',label="3 Gyr,lowm",linewidth=2)
    pylab.plot(lz3hm.mass,r3hm,'b--',label="3 Gyr, highm",linewidth=2)
    pylab.plot(lz1lm.mass,r1lm,'k-',label="1 Gyr,lowm",linewidth=2)
    pylab.plot(lz1hm.mass,r1hm,'k--',label="1 Gyr, highm",linewidth=2)

    pylab.axis([8.e13,3.e15,.7,6.])

    ax=pylab.gca()
    ax.set_xscale('log')
    #pylab.xticks((1.e14,1.e15),(r'$10^{14}$',r'10$^{15}$'),color='k',size=24,weight='ultrabold')
    ax.set_yscale('log')
    pylab.legend(loc='lower right')
    pylab.xlabel(r'{M$\rm{_{Cl}}$ (M$_\odot$)}',fontsize=30,fontweight=7)
    pylab.ylabel(r'{M$\rm{_{acc}}$(z=0.75) / M$\rm{_{acc}}$(z=0.07)}',fontsize=30,fontweight='bold')
    pylab.savefig('psratio2.eps')

def mass():
    pylab.cla()
    pylab.clf()
    for i in range(len(lz1lm.mass)):
	m=hz1hm.mass[i]
	l=lz1hm.maccret[i]
	h=hz1hm.maccret[i]
	r=hz1hm.frac[i]
	print i,m,h,r,m*r
    #print lz1lm.maccret
    #print hz1lm.maccret
    #print hz3lm.maccret
    #r3lm=(hz3lm.maccret)/(lz3lm.maccret)
    #r3hm=(hz3hm.maccret)/(lz3hm.maccret)
    #for i in range(len(r3)):
#	print i,lz3.sigma[i],hz3.sigma[i],lz3.mass[i],hz3.mass[i]
#	print i,lz01.sigma[i],hz01.sigma[i],lz01.mass[i],hz01.mass[i]
    #r1lm=hz1lm.maccret/lz1lm.maccret
    #r1hm=hz1hm.maccret/lz1hm.maccret
    #ra=N.array(hz01.maccret,'d')
    #rb=N.array(lz01.maccret,'d')
    #r01=ra/rb
    #for i in range(len(r01)):
	#print "ratio ",hz01.maccret[i],lz01.maccret[i],ra[i],rb[i],r01[i]
    pylab.subplot(221)
    masssub(lz3lm.mass,lz3lm.maccret,hz3lm.mass,hz3lm.maccret,'3Gyr, lowm')
    pylab.subplot(222)
    masssub(lz3hm.mass,lz3hm.maccret,hz3hm.mass,hz3hm.maccret,'3Gyr, highm')
    pylab.subplot(223)
    masssub(lz1lm.mass,lz1lm.maccret,hz1lm.mass,hz1lm.maccret,'1Gyr, lowm')
    pylab.subplot(224)
    masssub(lz1hm.mass,lz1hm.maccret,hz1hm.mass,hz1hm.maccret,'1Gyr, highm')
    #pylab.plot(lz3lm.mass,lz3lm.maccret,'b-',label="lowz, 3 Gyr,lowm",linewidth=2)
    #pylab.plot(lz3hm.mass,lz3hm.maccret,'b--',label="lowz, 3 Gyr, highm",linewidth=2)
    #pylab.plot(lz1lm.mass,lz1lm.maccret,'k-',label="lowz, 1 Gyr,lowm",linewidth=2)
    #pylab.plot(lz1hm.mass,lz1hm.maccret,'k--',label="lowz, 1 Gyr, highm",linewidth=2)

    #pylab.plot(hz3lm.mass,hz3lm.maccret,'r-',label="highz, 3 Gyr,lowm",linewidth=2)
    #pylab.plot(hz3hm.mass,hz3hm.maccret,'r--',label="highz, 3 Gyr, highm",linewidth=2)
    #pylab.plot(hz1lm.mass,hz1lm.maccret,'g-',label="highz, 1 Gyr,lowm",linewidth=2)
    #pylab.plot(hz1hm.mass,hz1hm.maccret,'g--',label="highz, 1 Gyr, highm",linewidth=2)

    pylab.text(7.e13,1.e11,r'{M$\rm{_{Cl}}$ (M$_\odot$)}',fontsize=30,fontweight=7,horizontalalignment='center')
    pylab.text(5.e11,2.e15,r'{M$\rm{_{acc}}$}',fontsize=30,rotation='vertical',fontweight='bold',verticalalignment='center')
    pylab.savefig('psmass.eps')
def masssub(x1,y1,x2,y2,label):
    pylab.plot(x1,y1,'k-')
    pylab.plot(x2,y2,'k--')
    pylab.text(2.e14,3.e14,label,fontsize=20.)
    pylab.axis([1.e14,3.e15,1.e12,1.e15])

    ax=pylab.gca()
    ax.set_xscale('log')
    pylab.xticks((1.e14,1.e15),(r'$10^{14}$',r'$10^{15}$'),color='k',size=24,weight='ultrabold')
    ax.set_yscale('log')
    pylab.legend(('Low-z','High-z'),loc='lower right')



def plotold():
    xmin=2.2
    xmax=3.2
    ymin=-2.5
    ymax=-.5
    psplotinit('fSsigma3Gyr.ps')
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,30)
    ppgplot.pglab("\gs (km/s)",'fS(10\u11\d:10\u13\d)',"")
    ppgplot.pgsci(1)
    ppgplot.pgline(sigma,frac)
    ppgplot.pgsls(2)
    ppgplot.pgsci(2)
    ppgplot.pgline(sigma08,frac08)
    ppgplot.pgsls(1)
    ppgplot.pgsci(1)
    
    ppgplot.pgend()


    xmin=2.2
    xmax=3.2
    ymin=11.
    ymax=14.2
    psplotinit('maccretsigma3Gyr.ps')
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0,30)
    ppgplot.pglab("\gs (km/s)",'M\dacc\u (M\d\(2281)\u)',"")
    ppgplot.pgsci(1)
    ppgplot.pgline(sigma,maccret)
    ppgplot.pgsls(2)
    ppgplot.pgsci(2)
    ppgplot.pgline(sigma08,maccret08)
    ppgplot.pgsls(1)
    ppgplot.pgsci(1)
    
    mylines=N.arange(-20.,20.,.4)
    mylineswidth=3
    ppgplot.pgsls(4)
    ppgplot.pgslw(mylineswidth)
    x=N.arange(0.,5.,1.)
    lines=mylines
    for y0 in lines:  
	y=3*x +y0 
	ppgplot.pgline(x,y)
	
	ppgplot.pgsls(1)
	ppgplot.pgend()
    os.system('cp maccretsigma.ps /Users/rfinn/SDSS/paper/.')
    os.system('cp fSsigma.ps /Users/rfinn/SDSS/paper/.')


def mratiopg():
    ppgplot.pgbeg("maccratio.ps/vcps",1,1)  #color port.
    ppgplot.pgpap(8.,1.)
    ppgplot.pgpage
    ppgplot.pgsch(1.3) #font size
    ppgplot.pgslw(7)   #line width

    # 1st panel with symbols w/ stddev errorbars
    #ylabel="SFR (M\d\(2281) \u yr\u-1\d)"
    ylabel="L(H\ga) (10\u41\d  erg s\u-1\d)"
    xlabel="M\dr\u "
    x1=.15
    x2=.5
    x3=.5
    x4=.85
    y1=x1
    y2=x2
    y3=x3
    y4=x4
    emarker=18
    smarker=23
    xmin=N.log10(1.e14)
    xmax=N.log10(2.5e15)
    #ymin=-1.
    #ymax=3.
    ymin=0.
    ymax=25.
    ppgplot.pgsvp(x1,x4,y1,y4)  #sets viewport
    ppgplot.pgswin(xmin,xmax,ymin,ymax) #axes limits
    ppgplot.pgbox('blncst',1.,2,'bcvnst',2.,2)  #tickmarks and labeling


    for i in range(len(lz1lm.mass)):
	m=lz1lm.mass[i]
	l=lz1lm.maccret[i]
	h=hz1lm.maccret[i]
	r=h/l
	print i,m,l,h,r
    #print lz1lm.maccret
    #print hz1lm.maccret
    #print hz3lm.maccret
    r3lm=(hz3lm.maccret)/(lz3lm.maccret)
    r3hm=(hz3hm.maccret)/(lz3hm.maccret)
    #for i in range(len(r3)):
#	print i,lz3.sigma[i],hz3.sigma[i],lz3.mass[i],hz3.mass[i]
#	print i,lz01.sigma[i],hz01.sigma[i],lz01.mass[i],hz01.mass[i]
    r1lm=hz1lm.maccret/lz1lm.maccret
    r1hm=hz1hm.maccret/lz1hm.maccret
    #ra=N.array(hz01.maccret,'d')
    #rb=N.array(lz01.maccret,'d')
    #r01=ra/rb
    #for i in range(len(r01)):
	#print "ratio ",hz01.maccret[i],lz01.maccret[i],ra[i],rb[i],r01[i]
    ppgplot.pgsci(14)
    ppgplot.pgsls(1)
    ppgplot.pgline(N.log10(lz3lm.mass),r3lm)
    ppgplot.pgsls(2)
    ppgplot.pgline(N.log10(lz3hm.mass),r3hm)

    ppgplot.pgsci(1)
    ppgplot.pgsls(1)
    ppgplot.pgline(N.log10(lz1lm.mass),r1lm)
    ppgplot.pgsls(2)
    ppgplot.pgline(N.log10(lz1hm.mass),r1hm)

    xlabel='M\dcl\u (M\d\(2281)\u)'
    ylabel='M\dacc\u(z=0.75) / M\dacc\u(z=0.07)'

    ppgplot.pgsch(1.8)
    ppgplot.pgslw(7)
    ppgplot.pgmtxt('b',2.2,0.5,0.5,ylabel)    #xlabel
    ppgplot.pgmtxt('l',2.5,0.5,0.5,xlabel)

    ppgplot.pgend()



    #for i in range(len(ybin)):
#	try:
#	    ybin[i]=N.log10(ybin[i])
#	except:
#	    ybin[i]=-99.


mrationew()
mass()
#mratiopg()
print "the end"


