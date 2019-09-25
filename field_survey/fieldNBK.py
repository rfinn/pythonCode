#!/usr/bin/env python
import pylab as py
import matplotlib
import Numeric as N
import mystuff as my

z1=2.2
z1=2.19
z1=2.193
#z1=2.194
z2=2.258
#z1=z2

def intfunc(wave,flux,wmin,wmax):
    sum=0.
    for i in range(len(wave)):
	if wave[i] > wmin:
	    #print wave[i],wmin,wmax,flux[i]
	    if wave[i] < wmax:
		#print 'dude!',wave[i],wmin,wmax,flux[i]
		sum=sum+(wave[i+1]-wave[i])*flux[i]
	if wave[i] > wmax:
	    break
    return sum
	
def avefunc(wave,flux,wmin,wmax):
    sum=0.
    npt=0
    for i in range(len(wave)):
	if wave[i] > wmin:
	    #print wave[i],wmin,wmax,flux[i]
	    if wave[i] < wmax:
		#print 'dude!',wave[i],wmin,wmax,flux[i]
		sum=sum+flux[i]
		npt+=1
	if wave[i] > wmax:
	    break
    sum=sum/(1.*npt)
    return sum
	

matplotlib.rc('legend',fontsize=14)
lmin=1.9
lmax=2.4
readR=0
if readR > 0:
    waveo=[]
    fluxo=[]
    infile=open('Rousselot-OH-spectrum_v2.0.dat','r')
    for line in infile:
	if line.find('#') > -1:
	    continue
	t=line.split()
	if float(t[0])/10000. < lmin:
	    continue
	if float(t[0])/10000. > lmax:
	    break

	waveo.append(float(t[0])/10000.)#convert from A to um
	fluxo.append(float(t[1]))
    infile.close()
    waveo=py.array(waveo,'f')
    fluxo=py.array(fluxo,'f')
    fluxo=fluxo/max(fluxo)
    py.plot(waveo,fluxo,'r-')

wavea=[]
fluxa=[]
infile=open('atmosph_trans-0.90-6.0um.dat','r')
for line in infile:
    if line.find('#') > -1:
	continue
    t=line.split()
    if float(t[0]) < lmin:
	continue
    if float(t[0]) > lmax:
	break

    wavea.append(float(t[0]))#wavelength in um
    fluxa.append(float(t[1]))
infile.close()
wavea=py.array(wavea,'f')
fluxa=py.array(fluxa,'f')

fluxa=N.take(fluxa,N.argsort(wavea))#sort y according to x rankings
wavea=N.take(wavea,N.argsort(wavea))


waveoh=[]
fluxoh=[]
infile=open('nearIR_skybg_2.dat','r')
for line in infile:
    if line.find('#') > -1:
	continue
    t=line.split()
    if float(t[0])/1000. < lmin:
	continue
    if float(t[0])/1000. > lmax:
	break

    waveoh.append(float(t[0])/1000.)#convert from nanometers to um
    fluxoh.append(float(t[1]))
infile.close()
waveoh=py.array(waveoh,'f')
fluxoh=py.array(fluxoh,'f')
fluxoh=fluxoh/max(fluxoh)

fluxoh=N.take(fluxoh,N.argsort(waveoh))#sort y according to x rankings
waveoh=N.take(waveoh,N.argsort(waveoh))


waveNBJ=[]#narrowband J at 1187 (Halpha at z=0.804)
fluxNBJ=[]
infile=open('1187-filter.dat','r')
for line in infile:
    if line.find('#') > -1:
	continue
    t=line.split()

    waveNBJ.append(float(t[0])/1000.)
    fluxNBJ.append(float(t[2]))#column 1 is warm data, column 2 is cold
infile.close()
waveNBJ=py.array(waveNBJ,'f')
fluxNBJ=py.array(fluxNBJ,'f')

fluxNBJ=fluxNBJ/max(fluxNBJ)
#assumed NB filter detects OII - what is observed wavelength of Halpha?
waveNBJ=0.6563/0.3727*waveNBJ


waveK=[]#narrowband J at 1187 (Halpha at z=0.804)
fluxK=[]
infile=open('Newfirm-Kband-filter.dat','r')
for line in infile:
    if line.find('#') > -1:
	continue
    t=line.split()

    waveK.append(float(t[0])/1000.)
    fluxK.append(float(t[2]))#column 1 is warm data, column 2 is cold
infile.close()
waveK=py.array(waveK,'f')
fluxK=py.array(fluxK,'f')

fluxK=fluxK/max(fluxK)
#sort arrays
fluxK=N.take(fluxK,N.argsort(waveK))#sort y according to x rankings
waveK=N.take(waveK,N.argsort(waveK))


def plottrans():

    py.plot(waveoh,fluxoh,'k-',label='OH')
    py.plot(wavea,fluxa,'k--',label="Atmos Trans")
    
    py.plot(waveNBJ,fluxNBJ,'b-',label="NB J OII-Ha")
    
    py.plot(waveK,fluxK,'r-',label="K-band")

    width=.01#percent filter width

    lobs=.6563*(1.+z1)
#py.axvline(x=lobs,color='r',label='z=2.17')
    l2=(1.-width/2.)*lobs
    s='z='+str(z1)
    py.axvline(x=l2,color='c',label=s,ls='-')
    l2=(1.+width/2.)*lobs
    py.axvline(x=l2,color='c',label='_nolegend_',ls='-')
    l2=0.992*(1.-width/2.)*lobs#blue shift corresponding to incidence angle of 13.5 degrees
    py.axvline(x=l2,color='c',label=r'$\theta_i$=13.5',ls='--')
    
    

    s='z='+str(z2)
    lobs=.6563*(1.+z2)
#py.axvline(x=lobs,color='g',label=s)
    l2=(1.-width/2.)*lobs
    py.axvline(x=l2,color='g',label=s,ls='-')
    l2=(1.+width/2.)*lobs
    py.axvline(x=l2,color='g',label='_nolegend_',ls='-')
    
    l2=0.992*(1.-width/2.)*lobs#blue shift corresponding to incidence angle of 13.5 degrees
    py.axvline(x=l2,color='g',label=r'$\theta_i$=13.5',ls='--')
    
#py.legend(loc='center')
    py.legend(loc='best')
    py.xlabel(r"Wavelength ($\mu$m)",fontsize=24)
    py.ylabel("Arbitrary",fontsize=24)
    py.axis([2.05, 2.18, -0.05, 1.1])
    py.savefig('Kplot.eps')

def plotsenstheta():
    py.cla()
    py.clf()
    lcenter=.6563*(1.+z1)
    width=0.01
    dtheta=0.1
    tolerance=0.
    theta=N.arange(0.,(13.+dtheta),dtheta)
    (sumsky1,trans1,sens1,sumatm1,sumK1)=plotsenssub(lcenter,width,theta,)

    lcenter=.6563*(1.+z2)
    (sumsky2,trans2,sens2,sumatm2,sumK2)=plotsenssub(lcenter,width,theta)

    sens1=sens1/max(sens1)
    sens2=sens2/max(sens2)
    s="z="+str(z1)
    py.plot(theta,sens1,'c-',label=s)
    s="z="+str(z2)
    py.plot(theta,sens2,'g-',label=s)
    sens2dL=sens2*(my.dL(z1,.7)/my.dL(z2,.7))**2
    s="z="+str(z2)+' dL weight'
    py.plot(theta,sens2dL,'g--',label=s)
    py.axvline(x=9.545,color='r',label=r'Edge of Field',ls='--')
    py.legend(loc='best')
    py.axis([0.,14.,0.3,1.05])
    py.xlabel(r"$\theta$ (degrees)",fontsize=24)
    py.ylabel(r"Relative Sensitivity",fontsize=24)

    py.savefig('senstheta.eps')

    #for i in range(len(theta)):
	#print i,theta[i],sumsky1[i],sumatm1[i],sumK1[i],trans1[i],sens1[i],sumsky2[i],sumatm2[i],sumK2[i],trans2[i],sens2[i]


def plotsenssub(lcenter,width,theta):
    n=1.9#index of refraction
    #tolerance=4.2#in nm
    tolerance=0.0#in nm
    tolerance=tolerance/1000.#convert to um
    lmin=(1-width/2.)*lcenter-tolerance
    lmax=(1+width/2.)*lcenter+tolerance
    sumsky=N.zeros(len(theta),'f')
    sumatm=N.zeros(len(theta),'f')
    sumK=N.zeros(len(theta),'f')
    trans=N.zeros(len(theta),'f')
    sens=N.zeros(len(theta),'f')
    for i in range(len(theta)):
	thet=theta[i]*N.pi/180.
	l1=lmin*N.cos(thet/n)
	l2=lmax*N.cos(thet/n)
	#sumsky[i]=intfunc(waveoh,fluxoh,l1,l2)
	#a=intfunc(wavea,fluxa,l1,l2)
	#b=intfunc(waveK,fluxK,l1,l2)
	sumsky[i]=intfunc(waveoh,fluxoh,l1,l2)
	sumatm[i]=avefunc(wavea,fluxa,l1,l2)
	sumK[i]=avefunc(waveK,fluxK,l1,l2)
	
	trans[i]=sumatm[i]*sumK[i]
	#print 'theta = ',theta[i],'sumsky=',sumsky[i],trans[i],a,b
    sens=trans/N.sqrt(10000.*sumsky)

    return sumsky,trans,sens,sumatm,sumK
	
def plotsensthetatol():
    py.cla()
    py.clf()
    lcenter=.6563*(1.+z1)
    width=0.01
    dtheta=0.1
    tolerance=0.
    theta=N.arange(0.,(13.+dtheta),dtheta)
    (sumsky1,trans1,sens1,sumatm1,sumK1)=plotsenssub2(lcenter,width,theta,tolerance)

    #lcenter=.6563*(1.+z2)
    tolerance=4.2
    tolerance=3.
    #tolerance=0.
    (sumsky2,trans2,sens2,sumatm2,sumK2)=plotsenssub2(lcenter,width,theta,tolerance)

    scale=max(sens1)
    sens1=sens1/scale
    sens2=sens2/scale
    s="z="+str(z1)
    py.plot(theta,sens1,color='0.4',ls='-',label=s)
    s=s+' w/tolerance'
    py.plot(theta,sens2,'k-',label=s)
    #sens2dL=sens2*(my.dL(z1,.7)/my.dL(z2,.7))**2
    s="z="+str(z2)+' dL weight'
    #py.plot(theta,sens2dL,'g--',label=s)
    py.axvline(x=9.545,color='r',label=r'Edge of Field',ls='--')
    py.legend(loc='best')
    s='$\lambda_c$ = %6.4f $\mu$m ,~ width = %5.3f,~ tolerance = %4.1f nm'%(lcenter,width,tolerance)
    print s
    py.title(s,fontsize=24)
    py.axis([0.,14.,0.1,1.05])
    py.xlabel(r"$\theta$ (degrees)",fontsize=24)
    py.ylabel(r"Relative Sensitivity",fontsize=24)

    py.savefig('sensthetatolerance.eps')

    #for i in range(len(theta)):
	#print i,theta[i],sumsky1[i],sumatm1[i],sumK1[i],trans1[i],sens1[i],sumsky2[i],sumatm2[i],sumK2[i],trans2[i],sens2[i]

def plotsenssub2(lcenter,width,theta,tolerance):
    n=1.9#index of refraction
    #tolerance=4.2#in nm
    #tolerance=0.0#in nm
    tolerance=tolerance/1000.#convert from nm to um
    lmin=(1-width/2.)*lcenter-tolerance
    lmax=(1+width/2.)*lcenter+tolerance
    sumsky=N.zeros(len(theta),'f')
    sumatm=N.zeros(len(theta),'f')
    sumK=N.zeros(len(theta),'f')
    trans=N.zeros(len(theta),'f')
    sens=N.zeros(len(theta),'f')
    for i in range(len(theta)):
	thet=theta[i]*N.pi/180.
	l1=lmin*N.cos(thet/n)
	l2=lmax*N.cos(thet/n)
	#sumsky[i]=intfunc(waveoh,fluxoh,l1,l2)
	#a=intfunc(wavea,fluxa,l1,l2)
	#b=intfunc(waveK,fluxK,l1,l2)
	sumsky[i]=intfunc(waveoh,fluxoh,l1,l2)
	sumatm[i]=avefunc(wavea,fluxa,l1,l2)
	sumK[i]=avefunc(waveK,fluxK,l1,l2)
	
	trans[i]=sumatm[i]*sumK[i]
	#print 'theta = ',theta[i],'sumsky=',sumsky[i],trans[i],a,b
    sens=trans/N.sqrt(10000.*sumsky)

    return sumsky,trans,sens,sumatm,sumK

#plottrans()
#plotsenstheta()
plotsensthetatol()
	
