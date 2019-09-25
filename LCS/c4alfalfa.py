#!/usr/bin/env python


from pylab import *
import sqlcl
import numarray as N
import mystuff as my
from matplotlib import rc
cra=[]
cdec=[]
cz=[]
csig=[]
cid=[]
cLx=[]
cTx=[]
ra1min=7.5*15.#8hr
ra1max=16.5*15#15hr
ra2min=0.*15.#0hr
ra2max=3.*15#3hr

ra3min=23.*15.#23hr
ra3max=0.*15#0hr

ncl=0


class Galaxy:
    def __init__(self):
        self.clusterid = []
        self.clusterrvir = []
        self.clusterr200 = []
        self.clustersigma = []
        self.clusterz = []
        self.ra = []
        self.dec  = []
        self.z = []
        self.Mabs=[]
        self.dr=[]
        self.dv=[]
        self.dLMpc=[]
        self.fHa=[]
        self.lHa=[]
        self.apcor=[]
        self.stellarmass=[]
        self.az=[]
        self.ar=[]
        self.O3Hb=[]
        self.N2Ha=[]
        self.sfr=[]
        self.ew=[]
        self.agn1=[]
        self.agn2=[]
        self.agn3=[]
        self.agn4=[]
        self.agn5=[]        

    def readdatafile(self,file):
        output=open(file,'r')
        i=-1
        for line in output:
            if (i == -1):
                fields=line.split()
                ncl=int(fields[1])#number of clusters
                print "number of galaxies = ",ncl
                self.clusterid=N.zeros(ncl,'f')
                self.clusterz=N.zeros(ncl,'f')
                self.clusterrvir=N.zeros(ncl,'f')
                self.clusterr200=N.zeros(ncl,'f')
                self.clustersigma=N.zeros(ncl,'f')
                self.dr=N.zeros(ncl,'f')
                self.dv=N.zeros(ncl,'f')
                self.ew=N.zeros(ncl,'f')
                self.lHa=N.zeros(ncl,'f')
                self.sfr=N.zeros(ncl,'f')
                self.stellarmass=N.zeros(ncl,'f')
                self.agn=N.zeros(ncl,'f')
                self.memb=N.zeros(ncl,'f')
                self.Mabs=N.zeros(ncl,'f')
                self.ra=N.zeros(ncl,'f')
                self.dec=N.zeros(ncl,'f')
                self.sigma10=N.zeros(ncl,'f')
                i += 1
                continue
            if line.find('#') > -1:#skip header
                continue
            f=line.split()
            for j in range(len(f)):
                f[j]=float(f[j])
            #(self.clusterid[i],self.clusterz[i],self.clusterrvir[i],self.clusterr200[i],self.clustersigma[i],self.dr[i],self.dv[i],self.ew[i],self.lHa[i],self.sfr[i],self.stellarmass[i],self.agn[i],self.memb[i],self.Mabs[i],self.ra[i],self.dec[i],self.sigma10[i])=(f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16])
            #print i, len(f)
            (self.clusterid[i],self.clusterz[i],self.clusterrvir[i],self.clusterr200[i],self.clustersigma[i],self.dr[i],self.dv[i],self.ew[i],self.lHa[i],self.sfr[i],self.stellarmass[i],self.agn[i],self.memb[i],self.Mabs[i],self.ra[i],self.dec[i])=(f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15])
            
            i += 1
        output.close()



def getclusters(catalog):
    ncl=0
    if abs(catalog-0) < .01:
	in1=open('/Users/rfinn/SDSS/SDSS-ALFALFA/myclusters.cat','r')


    if abs(catalog-1) < .01:
	in1=open('/Users/rfinn/proposals/Spitzer-Cycle5/cirstab1.dat','r')

    if abs(catalog-2) < .01:
	in1=open('/Users/rfinn/proposals/Spitzer-Cycle5/fastmfn.rose.dat','r')

    if abs(catalog-3) < .01:
	in1=open('/Users/rfinn/proposals/Spitzer-Cycle5/fastmfnsouth.rose.dat','r')
    if abs(catalog-4) < .01:
	in1=open('/Users/rfinn/proposals/Spitzer-Cycle5/groups.dat','r')
    for line in in1:
	
	if line.find('#') > -1:#skip header
	    continue
	f=line.split()
	if (catalog < .1):
	    for j in range(len(f)):
		f[j]=float(f[j])

	    cl=str(f[0])
	    if cl.find('1221') > -1:
		continue
	    if cl.find('1574') > -1:
		continue
	    if cl.find('3934') > -1:
		continue
	    id=f[0]
	    ra=f[36]
	    dec=f[37]
	    z=f[1]
	    sig=f[4]
	    Lx=0
	    Tx=0
	if (abs(catalog -1)< .01):
	    id=f[0]
	    ra=float(f[1])
	    dec=float(f[2])
	    z=float(f[3])
	    sig=float(line[63:66])
	    Lx=float(f[4])
	    Tx=float(line[58:62])
	if (catalog > 1.9) & (catalog < 4):
	    try:
		id=f[6]
	    except:
		id=f[5]

	    ra=float(f[0])
	    dec=float(f[1])
	    z=float(f[2])
	    Lx=float(f[4])
	    Tx=0.
	    sig=-1
	if (catalog > 3.1):
	    id=f[0]
	    rah=f[1]
	    rah=rah.split(':')

	    ra=15.*(float(rah[0])+float(rah[1])/60.+float(rah[2])/3600.)
	    rah=f[2]
	    rah=rah.split(':')
	    dec=(float(rah[0])+float(rah[1])/60.+float(rah[2])/3600.)
	    z=float(f[3])
	    #print id, ra,dec,z
	    Lx=0.
	    Tx=0.
	    sig=-1
	s=str(id)
#	if s.find('A2199') > -1:
#	    ncl += 1
#	    
#	    cra.append(ra)
#	    cdec.append(dec)
#	    cz.append(z)
#	    csig.append(sig)
#	    cid.append(id)
#	    cLx.append(Lx)
#	    cTx.append(Tx)
#	    continue
	
	if (dec > 0.) & (dec < 36.):
	    if (z < .039):
		if (z > .02):
		    if (ra > ra1min) & (ra < ra1max):
			ncl += 1
	    
			#print "%2i Cluster: id=%13s, ra=%6.2f (%6.3f hr), dec=%5.2f, z=%4.3f, sigma=%6.1f, Lx=%5.3f, Tx=%3.2f"%(ncl,str(id),ra,ra/15.,dec,z,sig,Lx,Tx)
			#print "%2i %13s %s %7.3f %s %7.3f %s %4.3f %s %6.1f %s %5.3f %s %3.2f"%(ncl,str(id),'%',ra,'%',dec,'%',z,'%',sig,'%',Lx,'%',Tx)
			cra.append(ra)
			cdec.append(dec)
			cz.append(z)
			csig.append(sig)
			cid.append(id)
			cLx.append(Lx)
			cTx.append(Tx)
			continue
		    elif (ra > ra2min) & (ra < ra2max):
			ncl +=1
			#print "%2i Cluster: id=%13s, ra=%6.2f (%6.3f hr), dec=%5.2f, z=%4.3f, sigma=%6.1f, Lx=%5.3f, Tx=%3.2f"%(ncl,str(id),ra,ra/15.,dec,z,sig,Lx,Tx)
			#print "%2i %13s %s %7.3f %s %7.3f %s %4.3f %s %6.1f %s %5.3f %s %3.2f"%(ncl,str(id),'%',ra,'%',dec,'%',z,'%',sig,'%',Lx,'%',Tx)
			cra.append(ra)
			cdec.append(dec)
			cz.append(z)
			csig.append(sig)
			cid.append(id)
			cLx.append(Lx)
			cTx.append(Tx)
			
			continue
		    elif (ra > ra3min) & (ra < ra3max):
			ncl += 1
			#print "%2i Cluster: id=%13s, ra=%6.2f (%6.3f hr), dec=%5.2f, z=%4.3f, sigma=%6.1f, Lx=%5.3f, Tx=%3.2f"%(ncl,str(id),ra,ra/15.,dec,z,sig,Lx,Tx)
			#print "%2i %13s %s %7.3f %s %7.3f %s %4.3f %s %6.1f %s %5.3f %s %3.2f"%(ncl,str(id),'%',ra,'%',dec,'%',z,'%',sig,'%',Lx,'%',Tx)
			cra.append(ra)
			cdec.append(dec)
			cz.append(z)
			csig.append(sig)
			cid.append(id)
			cLx.append(Lx)
			cTx.append(Tx)

			
			continue


    in1.close()



def getsdsscatalogs():
    for i in range(len(cz)):
    #for i in range(1):
	cname=str(cid[i])
	if cname.find('1367') < 0:
	    continue
	print i, cid[i]

	drsearch=5.*60.#search radius in arcmin for sdss query
	zmin=cz[i]-.005
	zmax=cz[i]+.005
	#print i,cid[i]," ra, dec, dr, mr = %12.8f %12.8f %8.3f %5.2f" % (cra[i],cdec[i],drsearch)
	query="select n.distance,g.ra,g.dec, g.u, g.g, g.r, g.i, g.z, s.z,l.ew,l.ewErr from galaxy g, specobj s, specline l, dbo.fGetNearbyObjEq(%12.8f,%12.8f,%8.3f) n where g.objid = s.bestobjid and g.objID = n.objID and l.specobjid = s.specobjid and s.z < %5.4f and s.z > %5.4f and (g.PrimTarget & 0x00000040) > 0 and l.LineId = 6565 order by distance" % (cra[i],cdec[i],drsearch,zmax,zmin)

	#print query
	try:
	    lines=sqlcl.query(query).readlines()
	except IOError:
	    print "IOError for cluster",self.id[i],i," trying spec query again"
	    lines=sqlcl.query(query).readlines()
	print "got number + 1 of spec objects = ",len(lines)
	n=str(cid[i])+'galaxy.dat'
	outfile=open(n,'w')
	j=0
	if (len(lines) > 1.):
	    for line in lines[1:]:
		if j < 0:
		    print line
		    j=j+1
		outfile.write(line)
	outfile.close()


def plotsdssgals():
    clf()
    cla()
    l=0
    npanel=1
    rmax=N.sqrt(2.)*60.
    subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.9,wspace=0.01,hspace=0.01)
    for j in range(len(cz)):
	n=str(cid[j])+'galaxy.dat'

	ngal=0
	infile=open(n,'r')
	#print cid[j],cra[j],cdec[j]
	for line in infile:
	    f=line.split(',')
	    if abs(float(f[1])-cra[j]) > .71:
		continue
	    if abs(float(f[2])-cdec[j]) > 1.:
		continue
	    ngal=ngal+1
	infile.close()


	gdr=   zeros(ngal,'f')
	gra=   zeros(ngal,'f')
	gdec=  zeros(ngal,'f')
	gz=    zeros(ngal,'f')
	gew=   zeros(ngal,'f')
	gewerr=zeros(ngal,'f')
	mu=    zeros(ngal,'f')
	mg=    zeros(ngal,'f')
	mr=    zeros(ngal,'f')
	mi=    zeros(ngal,'f')
	mz=    zeros(ngal,'f')
    

	infile=open(n,'r')
	i=0
	for line in infile:
	    f=line.split(',')
	    #if float(f[0]) > rmax:
	#	continue
	    if abs(float(f[1])-cra[j]) > .71:
		continue
	    if abs(float(f[2])-cdec[j]) > 1.:
		continue

	    for k in range(len(f)):
		f[k] = float(f[k])
	    #print i,len(gra),len(f),f
	    #(gdr[i],gra[i],gdec[i],mu[i],mg[i],mr[i],mi[i],mz[i],gz[i],gew[i],gewerr[i])=f    
	    gra[i]=f[1]
	    gdec[i]=f[2]
	    gew[i]=f[9]
	    i=i+1
	dra=cra[j]-gra
	ddec=gdec-cdec[j]
	hra=[]
	hdec=[]
	for i in range(len(gdr)):
	    if gew[i] > 1.:
		hra.append(dra[i])
		hdec.append(ddec[i])
	hra=array(hra,'f')
	hdec=array(hdec,'f')
#	print "%13s %s %7.3f %s %7.3f %s %4.3f %s %6.1f %s %4i %s %4i %s %3.2f %s %5.3f %s %3.2f \\\\ "%(str(cid[j]),'&',cra[j],'&',cdec[j],'&',cz[j],'&',csig[j],'&',len(dra),'&',len(hra),'&',1.*len(hra)/(1.*len(dra)),'&',cLx[j],'&',cTx[j])
	#print "%13s %s %4.3f %s %6.1f %s %4i %s %4i %s %3.2f %s %5.3f %s %3.2f \\\\ "%(str(cid[j]),'&',cz[j],'&',csig[j],'&',len(dra),'&',len(hra),'&',1.*len(hra)/(1.*len(dra)),'&',cLx[j],'&',cTx[j])
	try:
	    print "%13s %s %4.3f %s %6.1f %s %4i %s %4i %s %3.2f %s %5.3f \\\\ "%(str(cid[j]),'&',cz[j],'&',csig[j],'&',len(dra),'&',len(hra),'&',1.*len(hra)/(1.*len(dra)),'&',cLx[j])
	except ZeroDivisionError:
	    print "%13s %s %4.3f %s %6.1f %s %4i %s %4i %s %3.2f %s %5.3f \\\\ "%(str(cid[j]),'&',cz[j],'&',csig[j],'&',len(dra),'&',len(hra),'&',1.*len(hra),'&',cLx[j])
	if npanel == 1:
	    subplot(231)
	if npanel == 2:
	    subplot(232)
	if npanel == 3:
	    subplot(233)
	if npanel == 4:
	    subplot(234)
	if npanel == 5:
	    subplot(235)
	if npanel == 6:
	    subplot(236)
	if npanel == 7:
	    subplot(237)
	if npanel == 8:
	    subplot(238)
	if npanel == 9:
	    subplot(239)


	plot(dra,ddec,'k.',markersize=4)
	plot(hra,hdec,'bo',markersize=6)

	x=(-.71,-.71,.71,.71)
	y=(-1,1,1,-1)
	fill(x,y,fc=None,ec='k',lw=1)

	#x=array([-.5,-.5,.5,.5],'f')
	#y=array([-1.5,1.5,1.5,-1.5],'f')
	#fill(x,y,fc=None,ec='r',lw=1)
	#fill(y,x,fc=None,ec='g',lw=1)


	#if (npanel == 5):
	#    x=x+.5
	#    y=y+.5
	#    fill(x,y,fc=None,ec='r',lw=1)
	#    fill(y,x,fc=None,ec='g',lw=1)

	#elif (npanel == 7):
	#    x1=x-.5
	#    y1=y+.5
	#    fill(x1,y1,fc=None,ec='r',lw=1)
	#    x1=x+.5
	#    y1=y-.5
	#    fill(y1,x1,fc=None,ec='g',lw=1)

	#else:
	#    fill(x,y,fc=None,ec='r',lw=1)
	#    fill(y,x,fc=None,ec='g',lw=1)


	max=3.2
	limits=max*array([-1,1,-1,1],'f')
	axis(limits)

	ticks=("","","","","","","")
	#spot=(npanel/3.-int(npanel/3.))*3.
	if npanel == 1:
	    xticks(arange(-3,4,1),ticks)
	    yticks(arange(-3,4,1))#

	if npanel == 2:
	    xticks(arange(-3,4,1),ticks)
	    yticks(arange(-3,4,1),ticks)#

	if npanel == 3:
	    xticks(arange(-3,4,1),ticks)
	    yticks(arange(-3,4,1),ticks)
	if npanel == 4:
	    xticks(arange(-3,4,1),ticks)
	    yticks(arange(-3,4,1))

	if npanel == 5:
	    xticks(arange(-3,4,1),ticks)
	    yticks(arange(-3,4,1),ticks)

	if npanel == 6:
	    xticks(arange(-3,4,1),ticks)
	    yticks(arange(-3,4,1),ticks)
	if npanel == 7:
	    xticks(arange(-3,4,1))
	    yticks(arange(-3,4,1))

	if npanel == 8:
	    xticks(arange(-3,4,1))
	    yticks(arange(-3,4,1),ticks)
	
	if npanel == 9:
	    xticks(arange(-3,4,1))
	    yticks(arange(-3,4,1),ticks)
	try:
	    nsf=len(hra)
	except:
	    nsf=0

	try:
	    ntot=len(dra)
	except:
	    ntot=0
	s=r"\rm $\sigma$=%5.1f \ L$_x$=%5.3f \ T$_x$=%4.3f"% (csig[j],cLx[j],cTx[j])
	s2=r"\rm Nsp=%i \ Nsf=%i \ z=%5.3f" % (ntot,nsf,cz[j])
	text(-3,-2.8,s,fontsize=14)
	text(-3,1.8,s2,fontsize=14)
	text(-3,2.5,str(cid[j]),fontsize=14)
	
	if npanel == 6:
	    text(-8,-6,r'$\Delta$RA (deg)',fontsize=24,horizontalalignment='center')
	    text(-19.5,7,r'$\Delta$Dec (deg)',fontsize=24,verticalalignment='center',rotation='vertical')
	    l=l+1
	    p='area'+str(l)+'.eps'
	    savefig(p)
	    npanel=0
	    cla()
	    clf()
	npanel += 1
    if npanel > 1:
	l=l+1
	p='area'+str(l)+'.eps'
	savefig(p)
	npanel=0
	cla()
	clf()


def plotvelhist():
    clf()
    cla()
    l=0
    npanel=1
    xmin=-4000.
    xmax=4000.
    nbin=40
    subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.9,wspace=0.1,hspace=0.1)
    for j in range(len(cz)):
	n=str(cid[j])+'galaxy.dat'

	ngal=0
	infile=open(n,'r')
	for line in infile:
	    f=line.split(',')
	    if float(f[0]) > 100.:
		continue
	    ngal=ngal+1
	infile.close()


	gdr=   zeros(ngal,'f')
	gra=   zeros(ngal,'f')
	gdec=  zeros(ngal,'f')
	gz=    zeros(ngal,'f')
	gew=   zeros(ngal,'f')
	gewerr=zeros(ngal,'f')
	mu=    zeros(ngal,'f')
	mg=    zeros(ngal,'f')
	mr=    zeros(ngal,'f')
	mi=    zeros(ngal,'f')
	mz=    zeros(ngal,'f')
    

	infile=open(n,'r')
	i=0
	for line in infile:
	    f=line.split(',')
	    if float(f[0]) > 100.:
		continue

	    for k in range(len(f)):
		f[k] = float(f[k])
	    #print i,len(gra),len(f),f
	    #(gdr[i],gra[i],gdec[i],mu[i],mg[i],mr[i],mi[i],mz[i],gz[i],gew[i],gewerr[i])=f    
	    gra[i]=f[1]
	    gdec[i]=f[2]
	    gew[i]=f[9]
	    gz[i]=f[8]
	    i=i+1
	dvel=(gz-cz[j])*3.e5#convert dz to dv
	if npanel == 1:
	    subplot(231)
	if npanel == 2:
	    subplot(232)
	if npanel == 3:
	    subplot(233)
	if npanel == 4:
	    subplot(234)
	if npanel == 5:
	    subplot(235)
	if npanel == 6:
	    subplot(236)
	if npanel == 7:
	    subplot(237)
	if npanel == 8:
	    subplot(238)
	if npanel == 9:
	    subplot(239)
	#print dvel
	(xbin,ybin,ybinerr)=my.horizontalhist(dvel,xmin,xmax,nbin)
	my.drawhistpylab(xbin,ybin,'b-')



	xmin=-2000.
	xmax=-1.*xmin
	xgauss=N.arange(xmin,xmax,100.)
	sig=300.
	ygauss=max(ybin)*exp(-1.*(xgauss/(1.414*sig))**2)
	plot(xgauss,ygauss,'r-',lw=1)

	sig=500.
	ygauss=max(ybin)*exp(-1.*(xgauss/(1.414*sig))**2)
	plot(xgauss,ygauss,'r-',lw=1)

	sig=700.
	ygauss=max(ybin)*exp(-1.*(xgauss/(1.414*sig))**2)
	plot(xgauss,ygauss,'r-',lw=1)

	sig=1000.
	ygauss=max(ybin)*exp(-1.*(xgauss/(1.414*sig))**2)
	plot(xgauss,ygauss,'r-',lw=1)
	limits=array([xmin,xmax,0.,30],'f')
	axis(limits)

	xtickl=N.arange(xmin,xmax,500.)
	ticks=("","","","","","","")

	ticks=("","","","","","","")

	#spot=(npanel/3.-int(npanel/3.))*3.
	text(-3,2.5,str(cid[j]),fontsize=14)
	
	if npanel == 9:
	    text(-8,-6,r'$\Delta$RA (deg)',fontsize=24,horizontalalignment='center')
	    text(-19.5,7,r'$\Delta$Dec (deg)',fontsize=24,verticalalignment='center',rotation='vertical')
	    l=l+1
	    p='velhist'+str(l)+'.eps'
	    savefig(p)
	    npanel=0
	    cla()
	    clf()
	npanel += 1
    if npanel > 1:
	l=l+1
	p='velhist'+str(l)+'.eps'
	savefig(p)
	npanel=0
	cla()
	clf()


def plotsample():
    subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.9,wspace=0.3,hspace=0.3)
    subplot(221)
    plot(cz,cLx,'bo')
    xlabel('z',fontsize=20)
    ylabel('$L_x$',fontsize=20)
    ax=gca()
    ax.set_yscale('log')
    yticks(fontsize=16)
    xticks(arange(0.01,.035,.005),('0.01','','0.02','','0.03'),fontsize=16)
    #xticks([0.0,.005,.01,.015,.02,.025])
    subplot(222)
    plot(cz,csig,'bo')
    ax=gca()
    #ax.set_yscale('log')
    yticks(fontsize=16)
    xticks(arange(0.01,.035,.005),('0.01','','0.02','','0.03'),fontsize=16)
    xlabel('z',fontsize=20)
    ylabel(r'$\sigma$',fontsize=20)
    
    subplot(223)
    plot(csig,cLx,'bo')
    xlabel('$\sigma$',fontsize=20)
    ylabel('$L_x$',fontsize=20)
    ax=gca()
    ax.set_yscale('log')
    yticks(fontsize=16)
    xticks(fontsize=16)

    savefig('sample.eps')

def plotsiglx():
    cla()
    clf()
    rc('ytick',labelsize=2)
    rc('xtick',labelsize=20)
    rc('legend',numpoints=1)  # the number of points in the legend line
    lx=N.array([.033,.096,.083,.55,2.58,1.94],'f')
    sig=N.array([361.,325.,464.,500.,562.,660.],'f')


    plot(sig,lx,'bD',label="Spitzer Cycle 5 (PI Finn)",markersize=16)
    xlabel('$\sigma$ (km/s)',fontsize=26)
    ylabel('L$_x$ (10$^{43}$ erg/s)',fontsize=26)
    
    #coma(A1656), virgo, hercules(A2151), A1367, A1185
    lxgto=N.array([7.01,0.251,0.110, 1.51, .11],'f')
    #txgto=N.array([7.01,0.251,0.110, 5.26],'f')
    siggto=N.array([1000.,776.,689., 745.,489],'f')

    #coma(A1656), hercules(A2151), A1367, A1185
    lxgto=N.array([7.01,0.110, 1.51, .11],'f')
    #txgto=N.array([7.01,0.251,0.110, 5.26],'f')
    siggto=N.array([1000.,689., 745.,489],'f')

    plot(siggto,lxgto,'ko',mfc='w',label="Spitzer Archive",markersize=12,mew=2)

    #A0779 A1142 MKW4 NGC4325 
    lx2=N.array([.024,.071,.091,.051],'f')
    sig2=N.array([528.,579,577,484],'f')
    #plot(sig2,lx2,'r^',label="Other potential targets",markersize=12,mew=2)


    ax=gca()
    axis([250.,1050.,.01,15.])
    ax.set_yscale('log')
    yticks(fontsize=30)
    xticks(arange(300.,1100.,100.),("$300$","","$500$","","$700$","","$900$",""),fontsize=30)
    legend(loc='upper left')



    savefig('siglx.eps')



print "C4 + ALFALFA Clusters"
catalog=0
getclusters(catalog)

print "CIRS + ALFALFA Clusters"
catalog=1
getclusters(catalog)

print "RINES + ALFALFA Clusters"
catalog=2
getclusters(catalog)

#print "RINES-South + ALFALFA Clusters"
#catalog=3
#getclusters(catalog)

#print "Groups + ALFALFA Clusters"
#catalog=4
#getclusters(catalog)


cz=array(cz,'f')
cLx=array(cLx,'f')
cTx=array(cTx,'f')
cra=array(cra,'f')
cdec=array(cdec,'f')
csig=array(csig,'f')


#getsdsscatalogs()
print "cluster ids = ",cid
plotsdssgals()
print "Average redshift = ",N.average(cz)
#plotvelhist()

plotsiglx()


def sfrhist():
    g=Galaxy()
    
    print "reading in galaxy data file"
    g.readdatafile("/Users/rfinn/SDSS/baloghDR5/mygalaxies.cat")
    
    sfr=N.compress(g.ew > 10.,g.sfr)
    
    xmin=0.1 
    xmax=100.1
    nbin=1000
    (xbin,ybin,ybinerr)=my.horizontalhist(sfr,xmin,xmax,nbin)
    my.drawhistpylab(xbin,ybin,'b-')
    xlabel(r'SFR (h$^{-2}_{70}$ M$_\odot$/yr)',fontsize=24)
    ylabel(r'Ngal',fontsize=24)
    
    
#hist(sfr,bins=1000,fc=None,ec='k')
    ax=gca()
    axvline(1.,linewidth=2,color='r',ls='--')
    axis([.1,22.,1.,5000.])
    
    ax.set_xscale("log")
    ax.set_yscale("log")
    savefig('spitzersfrs.eps')
