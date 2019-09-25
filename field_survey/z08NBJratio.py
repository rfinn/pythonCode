#! /usr/bin/env python
from pylab import *
from scipy.interpolate import splrep,splev,splint

nbx=[]
nby=[]
infile=open('/Users/rfinn/FieldSurvey/filter/transmission.dat','r')
for line in infile:
    t=line.split()
    nbx.append(float(t[0]))#wavelength in nm
    nby.append(float(t[1]))#transmission at 77K
#reverse order of arrays so that wavelength is increasing
nbx.reverse()
nby.reverse()
nbx=array(nbx,'f')
nby=array(nby,'f')

jx=[]
jy=[]
infile=open('/Users/rfinn/proposals/observing/NOAO2009B/newfirmJ.dat','r')
for line in infile:
    t=line.split()
    jx.append(float(t[0]))#wavelength in nm
    jy.append(float(t[1]))#transmission at 77K
jx=array(jx,'f')
jy=array(jy,'f')*100

ix=[]
iy=[]
infile=open('/Users/rfinn/FieldSurvey/filter/sdssi.dat','r')
for line in infile:
    if line.find('#') > -1:
	continue
    t=line.split()
    ix.append(float(t[0]))#wavelength in nm
    iy.append(float(t[1]))#transmission at 77K
ix=array(ix,'f')/10.
iy=array(iy,'f')

sbx=[]
sby=[]
infile=open('/Users/rfinn/snr/narrowband/nearIR_skybg_2.dat','r')
for line in infile:
    t=line.split()
    if float(t[0]) > 1400.:
	break
    sbx.append(float(t[0]))#wavelength in nm
    sby.append(float(t[1]))#transmission at 77K
sbx=array(sbx,'f')
sby=array(sby,'f')

ax=[]
ay=[]
infile=open('/Users/rfinn/proposals/observing/NOAO2009B/atmosph_short.lines.dat','r')
for line in infile:
    if line.find('#') > -1:
	continue
    t=line.split()
    ax.append(float(t[0]))#wavelength in nm
    ay.append(float(t[1]))#transmission at 77K
ax=array(ax,'f')*1000
ay=array(ay,'f')*100

#read in galaxy SEDs
#lambda E Sbc Scd Im SB2 SB1 10Myr

galx=[]
gal1=[]
gal2=[]
gal3=[]
gal4=[]
gal5=[]
gal6=[]
gal7=[]

infile=open('/Users/rfinn/SEDs/cww_extended.dat','r')
for line in infile:
    if line.find('#') > -1:
	continue
    t=line.split()
    galx.append(float(t[0]))#wavelength in A
    gal1.append(float(t[1]))#galaxy SED
    gal2.append(float(t[2]))#galaxy SED
    gal3.append(float(t[3]))#galaxy SED
    gal4.append(float(t[4]))#galaxy SED
    gal5.append(float(t[5]))#galaxy SED
    gal6.append(float(t[6]))#galaxy SED
    gal7.append(float(t[7]))#galaxy SED

galx=array(galx,'d')/10.
gal1=array(gal1,'d')
gal2=array(gal2,'d')
gal3=array(gal3,'d')
gal4=array(gal4,'d')
gal5=array(gal5,'d')
gal6=array(gal6,'d')
gal7=array(gal7,'d')
gals=zeros([7,len(galx)],'f')
#for i in range(7):
#    g='gal'+str(i+1)
#    gals[i]=g
scale=2.5/1.e26
gals[0]=gal1*scale
gals[1]=gal2*scale
gals[2]=gal3*scale
gals[3]=gal4*scale
gals[4]=gal5*scale
gals[5]=gal6*scale
gals[6]=gal7*scale
z=0.81
galxz=galx*(1+z)


def fitsplines():
    global tckJ, tckNB, tckSkyTran,tckSkyEm,tcki
    tckJ=splrep(jx,jy)
    tckNB=splrep(nbx,nby)
    tckSkyEm=splrep(sbx,sby)
    tckSkyTran=splrep(ax,ay)
    tcki=splrep(ix,iy)
    #return(tckJ,tckNB,tckSkyTran,tckSkyEm,tcki)


def plotseds():
    figure(7)
    clf()
    cla()
    for i in range(7):
	subplot(4,2,i+1)
	plot(galx,gals[i])



def getfluxes():
#fit splines

#redshift galaxy seds
    global fJ,fnb,fi
    fJ=zeros(7,'d')
    fnb=zeros(7,'d')
    fi=zeros(7,'d')

    for k in range(len(galx)):
	wave=galxz[k]
	if (wave > 1100.) & (wave < 1400.):
	    sky=splev(wave,tckSkyEm)
	    j=splev(wave,tckJ)/100.
	    atten=j*splev(wave,tckSkyTran)/100.
	    for i in range(len(fJ)):
		dwave=wave-galxz[k-1]
		#fJ[i]=fJ[i]+(gals[i,k]+sky)*atten*dwave
		fJ[i]=fJ[i]+(gals[i,k]+sky)*j
		if (wave > 1170.) & (wave < 1200):
		    nb=splev(wave,tckNB)/100.
		    #fnb[i]=fnb[i]+(gals[i,k]+sky)*atten*nb*dwave
		    fnb[i]=fnb[i]+(gals[i,k]+sky)*nb*j
		    print i,fJ[i],fnb[i],nb,j,dwave
	if (wave > 650.) & (wave < 880.):
	    dwave=wave-galxz[k-1]
	    sky=0#splev(wave,tckSkyEm)
	    atten=splev(wave,tcki)#*splev(wave,tckSkyTran)
	    for i in range(len(fJ)):
		fi[i]=fi[i]+(gals[i,k]+sky)*atten*dwave
    figure(6)
    cla()
    clf()
    print fJ,fnb,fi
    print fnb/fJ
    plot(2.5*log10(fi/fJ),fnb/fJ,'bo')
    xlabel('2.5 log(F(i)/F(J))')
    #ylabel('2.5 log(F(J)/F(NB))')
    ylabel('F(NB)/F(J)')


def getfluxes2():
    """
    z(l)=1.1/(1.0*zstep)*(l-1)
    c         write(*,*) l,z(l)
    do i=1,nfwav
    l1obs=fwav(i)/(1+z(l))
    l2obs=fwav(i+1)/(1+z(l))
    scale=(filter(i)+filter(i+1))/2
    do j=1,nwav
    if(j .eq. nwav)then
    dwav=wav(j)-wav(j-1)
    else
    dwav=wav(j+1)-wav(j)
    endif
    if(wav(j) .ge. l1obs .and. wav(j) .lt. l2obs)then
    c..flux is normalized by flux/A of vega
    fobsE=fobsE+E(j)*dwav*scale/f0
    fobsSbc=fobsSbc+Sbc(j)*dwav*scale/f0
    fobsScd=fobsScd+Scd(j)*dwav*scale/f0
    fobsIm=fobsIm+Im(j)*dwav*scale/f0
    c                  if (i .lt. 2) then
    c                     write(*,*) "hey",i,j,E(j),Sbc(j),Scd(j)
    c                 
    """
    global fJ,fnb,fi
    fJ=zeros(7,'d')
    fnb=zeros(7,'d')
    fi=zeros(7,'d')

    xmin=1050.
    xmax=1450.
    gx=galxz[(galxz>xmin)&(galxz<xmax)]
    gs=gals[:,(galxz>xmin)&(galxz<xmax)]
    for k in range(len(jx)-1):
	wave=jx[k]
	dwave=jx[k+1]-wave
	#j=(jy[k]+jy[k+1])/2.
	j=jy[k]/100.
	sky=splev(wave,tckSkyEm)
	atten=j*splev(wave,tckSkyTran)/100.
        mini=1000
	maxi=0
	for i in range(len(gx)):
	    if gx[i] > wave:
		mini=i-1
		maxi=i
		break
	for i in range(len(fJ)):
	    galflux=gals[i,mini]+(gs[i,maxi]-gs[i,mini])/(gx[maxi]-gx[mini])*(wave-gx[mini])
	    fJ[i]=fJ[i]+(galflux+sky)*atten
	    if (wave > 1170.) & (wave < 1200):
		nb=splev(wave,tckNB)/100.
		fnb[i]=fnb[i]+(galflux+sky)*nb*atten
		#print i,fJ[i],fnb[i],nb,j,dwave


    xmin=650.
    xmax=880.
    gx=galxz[(galxz>xmin)&(galxz<xmax)]
    gs=gals[:,(galxz>xmin)&(galxz<xmax)]
    for k in range(len(ix)-1):
	wave=ix[k]
	dwave=ix[k+1]-wave
	#j=(jy[k]+jy[k+1])/2.
	transi=iy[k]/100.
	#sky=splev(wave,tckSkyEm)
	#atten=transi*splev(wave,tckSkyTran)/100.
	atten=transi
	for i in range(len(gx)):
	    if gx[i] > wave:
		mini=i-1
		maxi=i
		break
	for i in range(len(fi)):
	    galflux=gs[i,mini]+(gs[i,maxi]-gs[i,mini])/(gx[maxi]-gx[mini])*(wave-gx[mini])
	    fi[i]=fi[i]+(galflux+sky)*atten

    figure(6)
    cla()
    clf()
    #print fJ,fnb,fi
    #print fnb/fJ
    plot(2.5*log10(fi/fJ),fnb/fJ,'bo')
    xlabel('2.5 log(F(i)/F(J))')
    #ylabel('2.5 log(F(J)/F(NB))')
    ylabel('F(NB)/F(J)')


def getnbj(z):
    xmin=1050.
    xmax=1450.
    galxz=galx*(1+z)
    gx=galxz[(galxz>xmin)&(galxz<xmax)]
    gs=gals[:,(galxz>xmin)&(galxz<xmax)]
    fJ=zeros(7,'d')
    fnb=zeros(7,'d')
    fi=zeros(7,'d')

    for k in range(len(jx)-1):
	wave=jx[k]
	dwave=jx[k+1]-wave
	#j=(jy[k]+jy[k+1])/2.
	j=jy[k]/100.
	sky=splev(wave,tckSkyEm)
	atten=j*splev(wave,tckSkyTran)/100.
        mini=1000
	maxi=0
	for i in range(len(gx)):
	    if gx[i] > wave:
		mini=i-1
		maxi=i
		break
	for i in range(len(fJ)):
	    galflux=gs[i,mini]+(gs[i,maxi]-gs[i,mini])/(gx[maxi]-gx[mini])*(wave-gx[mini])
	    fJ[i]=fJ[i]+(galflux+sky)*atten
	    if (wave > 1170.) & (wave < 1200):
		nb=splev(wave,tckNB)/100.
		fnb[i]=fnb[i]+(galflux+sky)*nb*atten
		#print i,fJ[i],fnb[i],nb,j,dwave

    return fJ,fnb

def plotnbjz():
    z=arange(0.1,2.5,.01)
    nbj=zeros([7,len(z)],'f')
    for i in range(len(z)):
	fJ,fnb=getnbj(z[i])
	for j in range(len(fJ)):
	    nbj[j,i]=fnb[j]/fJ[j]
    figure()
    colors=['r-','m-','y-','g-','c-','b-','k-']
    for i in range(7):
	plot(z,nbj[i],colors[i],lw=2)
    
    lab=('E','Sbc','Scd','Im','SB1','SB2','10Myr')
    legend(lab,loc='upper left')
    lobs=1187.
    lem=array([372.7,656.3,486.8,500.7])
    lab=('[OII]',r'$H\alpha$',r'$H\beta$','[OIII]')
    y=array([.2,.2,.17,.2])
    ax=gca()
    for i in range(len(lem)):
	zobs=lobs/lem[i]-1.
	print 'zobs = ',zobs
	axvline(x=zobs,color='k',ls='--')
	text(zobs,y[i],lab[i],fontsize=14)
    axis([0.,2.5,0.,.25])
    xticks(fontsize=14)
    yticks(fontsize=14)
    xlabel('Redshift',fontsize=18)
    ylabel('F(1187)/F(J)',fontsize=18)
    savefig('NB1187Jz.eps')

    
def plotfig1():
    figure(1)
    cla()
    clf()
    plot(nbx,nby,'b-')
    plot(jx,jy,'b-',lw=1)
    plot(sbx,sby,'k-',lw=1)
    plot(ax,ay,'r-',lw=1)
    xlabel('Wavelength (nm)',fontsize=26)
    ylabel('Transmission',fontsize=26)
    axis([1100,1400,0,105])
    savefig('plot1a.eps')

def plotfig2():
    figure(2)
    cla()
    clf()
    plot(nbx,nby,'b-',label='_nolegend_')
    plot(sbx,sby,'k-',lw=1,label='_nolegend_')
    plot(ax,ay,'r-',lw=1,label='_nolegend_')
    z=array([.813,.816,.805,.807,.794],'f')
    sigma=array([1000.,800.,1000.,1000.,1018.],'f')
    colors=['k','c','g','r','b']
    clnames=['RXJ1716.6+6708','RXJ1821', 'RDCS J1317+2911', 'RDCS J1350+6007', 'CL1216']
    x=arange(1000,1500,1)
    le=656.3
    c=3.e5
    a=-0.5*9e4
    ytext=.82
#ax=gca()
    for i in range(len(z)-1):
	if i == 1:
	    continue
	y=80.*exp(-0.5*(x-le*(1+z[i]))**2/(sigma[i]/c*le*(1+z[i]))**2)
	plot(x,y,color=colors[i],ls='--',label=clnames[i])
	
    leg=legend(loc='upper left',numpoints=4)
    for t in leg.get_texts():
	t.set_fontsize('14')
    axis([1170,1200,0,110])
    xlabel('Wavelength (nm)',fontsize=26)
    ylabel('Transmission',fontsize=26)

    savefig('plot1b.eps')

def plotfig3():
    figure(3)
    cla()
    clf()
    plot(jx,jy,'b-')
    wave=arange(1110.,1400.,1.)
    model=splev(wave,tckJ)
    plot(wave,model,'r-')


def plotfig4():
    figure(4)
    cla()
    clf()
    plot(galx,gal1,'r-')
    plot(galx,gal2,'b-')
    plot(galx,gal3,'k-')
    plot(galx,gal4,'g-')
    plot(galx,gal5,'c-')
    plot(galx,gal6,'k--')
    ax=gca()
    ax.set_yscale('log')
    xlabel('Wavelength (nm)')
    
fitsplines()
plotnbjz()
