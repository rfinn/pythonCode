#!/usr/bin/env python

'''
read in 'cluster'mastertable.fits 

(1) change mastertable string to the directory where you are storing the fits files.

(2) usage from w/in ipython -pylab
%run edsicsreadspitzercatalogs.py
cl = baseCluster('cl1216')

You must feed it the ediscs prefix as a string.  the possible names are
'cl1018'
'cl1037'
'cl1040'
'cl105412'
'cl105411'
'cl1059'
'cl1103'
'cl1138'
'cl1202'
'cl1216' 
'cl1227'
'cl1232'
'cl1301'
'cl1353'
'cl1354'
'cl1411'
'cl1420'

#then for example
figure()
plot(cl.ra,cl.dec,'k.')

(3) to use as an inherited class, include 

from ediscsreadspitzercatalogs import *

then when defining a cluster class:

class cluster(baseCluster):

when initiated, it will establish a new class and read in the mastertable

'''
from pylab import *
import pyfits, os
from ediscsCommon import *
import mystuff as my
mypath=os.getcwd()
if mypath.find('rfinn') > -1:

    if mypath.find('Users') > -1:
        print "Running on Rose's mac pro"
        homedir='/Users/rfinn/Dropbox/research-macbook/'
    elif mypath.find('home') > -1:
        print "Running on coma"
        homedir='/home/rfinn/research/'
    mastertablepath=homedir+'clusters/spitzer/MasterTables/'
else:
    mastertablepath=''

class baseCluster:
    def __init__(self,prefix):
        self.prefix=prefix
        self.fullname=fullname[self.prefix]
        self.cra=racenter[self.prefix]
        self.cdec=deccenter[self.prefix]
        self.cz=redshift[self.prefix]
        self.f80=F80[self.prefix]
        self.errf80=ErrorF80[self.prefix]
        self.csigma=sigma[self.prefix]
        self.csigmaerrplus=errsigmaplus[self.prefix]
        self.csigmaerrminus=errsigmaminus[self.prefix]
        self.r200=2.02*(self.csigma)/1000./sqrt(OmegaL+OmegaM*(1.+self.cz)**3)*H0/70. # in Mpc
        if self.r200 < .5:
            print 'WARNING:  R200 unrealistically small for ',self.prefix
            print 'resetting R200 to 0.5 Mpc'
            self.r200=.5
        self.r200deg=self.r200*1000./my.DA(self.cz,h)/3600.
        self.mcl=my.clusterMass(self.csigma,self.cz,h)

        mastertable=mastertablepath+self.fullname+'mastertable.fits'
        tb=pyfits.open(mastertable)
        tbdata=tb[1].data
        tb.close()

        self.ediscsID      =tbdata.field('EDISCS-ID')
        self.ediscsIDold  =tbdata.field('EDISCS-ID-OLD')
        self.ra             =tbdata.field('RA')
        self.dec            =tbdata.field('DEC')
        self.xcorr          =tbdata.field('xcorr')
        self.ycorr          =tbdata.field('ycorr')
        self.starflag       =tbdata.field('starflag')
        self.EW             =tbdata.field('EW')
        self.EWerr          =tbdata.field('EWerr')
        self.SFR            =tbdata.field('SFR') 
        self.SFRerr         =tbdata.field('SFRerr') 
        self.matchflaghalpha=tbdata.field('matchflaghalpha') 
        self.onHaimageflag  =tbdata.field('onHaimageflag') 
        self.sfflag         =tbdata.field('sfflag') 
        self.matchflag24    =tbdata.field('matchflag24') 
        self.on24imageflag  =tbdata.field('on24imageflag') 
        self.flux24         =tbdata.field('flux24') 
        self.flux24err      =tbdata.field('flux24err') 
        self.nmatchediscs24 =tbdata.field('nmatchediscs24') 
        self.snr24          =tbdata.field('snr24') 
        self.imagex24       =tbdata.field('imagex24') 
        self.imagey24       =tbdata.field('imagey24') 
        self.ra24           =tbdata.field('ra24') 
        self.dec24          =tbdata.field('dec24')
        self.flux80flag     =tbdata.field('flux80flag')
        self.L24            =tbdata.field('L24') 
        self.L24err         =tbdata.field('L24err') 
        self.Lir            =tbdata.field('Lir') 
        self.errLir         =tbdata.field('errLir') 
        self.SFRir          =tbdata.field('SFRir') 
        self.SFRirerr       =tbdata.field('SFRirerr') 
        self.matchflagediscsirac=tbdata.field('matchflagediscsirac')
        self.iracf1       =tbdata.field('iracf1')
        self.iracf2       =tbdata.field('iracf2')
        self.iracf3       =tbdata.field('iracf3')
        self.iracf4       =tbdata.field('iracf4')
        self.erriracf1    =tbdata.field('erriracf1')
        self.erriracf2    =tbdata.field('erriracf2')
        self.erriracf3    =tbdata.field('erriracf3')
        self.erriracf4    =tbdata.field('erriracf4')
        self.iracsexflag0 =tbdata.field('iracsexflag0')
        self.iracsexflag1 =tbdata.field('iracsexflag1')
        self.iracwch1     =tbdata.field('iracwch1')
        self.iracwch2     =tbdata.field('iracwch2')
        self.iracwch3     =tbdata.field('iracwch3')
        self.iracwch4     =tbdata.field('iracwch4')
        self.iracwmin     =tbdata.field('iracwmin')
        self.nmatchediscsirac     =tbdata.field('nmatchediscsirac') 
        self.matchflagmorphgimtype=tbdata.field('matchflagmorphgimtype') 
        self.gimtype              =tbdata.field('gimtype') 
        self.matchflagvistype     =tbdata.field('matchflagvistype') 
        self.vistype      =tbdata.field('vistype')
        self.misoV        =tbdata.field('misoV')
        self.misoeVapsim  =tbdata.field('misoeVapsim') 
        self.misoR        =tbdata.field('misoR')
        self.misoeRapsim  =tbdata.field('misoeRapsim') 
        self.misoI        =tbdata.field('misoI')
        self.misoeIapsim  =tbdata.field('misoeIapsim') 
        self.misoJ        =tbdata.field('misoJ')
        self.misoeJapsim  =tbdata.field('misoeJapsim') 
        self.misoK        =tbdata.field('misoK')
        self.misoeKapsim  =tbdata.field('misoeKapsim') 
        self.magV         =tbdata.field('magV') 
        self.mageVapsim   =tbdata.field('mageVapsim') 
        self.magR         =tbdata.field('magR') 
        self.mageRapsim   =tbdata.field('mageRapsim') 
        self.magI         =tbdata.field('magI') 
        self.mageIapsim   =tbdata.field('mageIapsim') 
        self.magJ         =tbdata.field('magJ') 
        self.mageJapsim   =tbdata.field('mageJapsim') 
        self.magK         =tbdata.field('magK') 
        self.mageKapsim   =tbdata.field('mageKapsim') 
        self.membflag     =tbdata.field('membflag') 
        self.newspecmatchflag=tbdata.field('newspecmatchflag')
        self.defmembflag     =tbdata.field('defmembflag') 
        self.photmembflag    =tbdata.field('photmembflag') 
        self.supermembflag   =tbdata.field('supermembflag')
        self.specz           =tbdata.field('specz') 
        self.spectype        =tbdata.field('spectype') 
        self.specEWOII       =tbdata.field('specEWOII')
        self.matchflagspecediscs  =tbdata.field('matchflagspecediscs') 
        self.specEWOIIflag =tbdata.field('specEWOIIflag') 
        self.bestz         =tbdata.field('bestz')
        self.lowz          =tbdata.field('lowz') 
        self.highz         =tbdata.field('highz')
        self.wmin          =tbdata.field('wmin')
        self.Pclust        =tbdata.field('Pclust')
        self.MR            =tbdata.field('MR') 
        self.MU            =tbdata.field('MU') 
        self.MV            =tbdata.field('MV') 
        self.MB            =tbdata.field('MB') 
        self.stellarmass     =tbdata.field('stellmass') 
        self.redflag       =tbdata.field('redflag') 
        self.LUlowzclust   =tbdata.field('LUlowzclust')
        self.LUbestzclust  =tbdata.field('LUbestzclust') 
        self.LUhighzclust  =tbdata.field('LUhighzclust')
        self.LBlowzclust   =tbdata.field('LBlowzclust') 
        self.LBbestzclust  =tbdata.field('LBbestzclust')
        self.LBhighzclust  =tbdata.field('LBhighzclust')
        self.LVlowzclust   =tbdata.field('LVlowzclust ')
        self.LVbestzclust  =tbdata.field('LVbestzclust')
        self.LVhighzclust  =tbdata.field('LVhighzclust')
        self.LRlowzclust   =tbdata.field('LRlowzclust') 
        self.LRbestzclust  =tbdata.field('LRbestzclust')
        self.LRhighzclust  =tbdata.field('LRhighzclust')
        self.LIlowzclust   =tbdata.field('LIlowzclust') 
        self.LIbestzclust  =tbdata.field('LIbestzclust')
        self.LIhighzclust  =tbdata.field('LIhighzclust')
        self.LJlowzclust   =tbdata.field('LJlowzclust') 
        self.LJbestzclust  =tbdata.field('LJbestzclust')
        self.LJhighzclust  =tbdata.field('LJhighzclust')
        self.LKlowzclust   =tbdata.field('LKlowzclust') 
        self.LKbestzclust  =tbdata.field('LKbestzclust') 
        self.LKhighzclust  =tbdata.field('LKhighzclust')


        # some extra quantities that are not included in the mastertable
        dr=sqrt((self.ra-self.cra)**2+(self.dec-self.cdec)**2)
        self.drflag = (dr < self.r200deg)


	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
