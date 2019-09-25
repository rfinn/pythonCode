#!/usr/bin/env python

import atpy
import os
from pylab import *
from ediscsCommon import *
mypath=os.getcwd()
if mypath.find('Users') > -1:
    print "Running on Rose's mac pro"
    homedir='/Users/rfinn/'
elif mypath.find('home') > -1:
    print "Running on coma"
    homedir='/home/rfinn/'

from LCSanalyzespirals import *
class cluster:
    def __init__(self,prefix):
        self.prefix=prefix
        gregfile=homedir+'research/WifiEdiscs/rudnick_catalog/'+self.prefix+'_megacat.v5.4.fits'
        self.gdat=atpy.Table(gregfile)
        # read in Finn mips catalog

        rfile=homedir+'research/ediscsClusters/MIPS_catalogs/'+fullname[self.prefix]+'mastertable.fits'
        self.rdat=atpy.Table(rfile)
        self.zcluster=redshift[self.prefix]
        self.cRA=racenter[self.prefix]
        self.cDec=deccenter[self.prefix]
        self.csigma=sigma[self.prefix]
        self.r200=2.02*(self.csigma)/1000./sqrt(OmegaL+OmegaM*(1.+self.zcluster)**3)*H0/70. # in Mpc
        self.r200deg=self.r200*1000./my.DA(self.zcluster,h)/3600.

        #self.plotLIRI()
    def plotLIRI(self,plotsingle=1):
        if plotsingle:
            figure()
            title(self.prefix)
            xlabel('$ R $',fontsize=20)
            ylabel('$ L_{IR} $',fontsize=20)
        distance=sqrt((self.gdat.RA-self.cRA)**2+(self.gdat.Dec-self.cDec)**2)
        keepflag=(abs(self.zcluster - self.gdat.zLDP) < .03) &( distance < 1.4*7./60)
        plot(self.gdat.Rauto[keepflag],self.gdat.LIR_ZCLUST[keepflag],'ko',mfc='w',mec='k',markersize=10,label='_nolegend_')
        n1=sum(keepflag & (self.gdat.LIR_ZCLUST>1.e11)& (self.gdat.Rauto < 23.))

        keepflag=self.rdat.matchflag24 & self.rdat.supermembflag & (self.rdat.newspecmatchflag < .1)
        plot(self.rdat.magR[keepflag],self.rdat.Lir[keepflag],'ko',mfc='w',markersize=10,label='Photoz')
        n2=sum(keepflag & (self.rdat.Lir>1.e11) & (self.rdat.magR < 23.))
        keepflag=self.rdat.matchflag24 & (self.rdat.newspecmatchflag > 0.1)
        plot(self.rdat.magR[keepflag],self.rdat.Lir[keepflag],'k^',markersize=8,label='Spec Memb')

        gca().set_yscale('log')
        axhline(y=10.**11,color='k',ls='--',lw=1)
        axvline(x=23.,color='k',ls='--',lw=1)
        axis([20.5,24.5,5.e10,2.e12])
        xticks(arange(21,25))
        if self.prefix.find('cl1216') > -1:
            legend(loc='upper right',numpoints=1,prop={'size':12})
        s=fullname[self.prefix]+'\nz=%4.3f \nN=%2i'%(self.zcluster,n1+n2)
        text(.1,.95,s,verticalalignment='top',horizontalalignment='left',transform=gca().transAxes,fontsize=12)
        #s='N Target=%i'%(n1+n2)
        #text(.1,.8,s,horizontalalignment='left',transform=gca().transAxes,fontsize=12)
        print self.prefix,': number of targets with log10(L_IR) > 11 = ',n1+n2


    def plotSFRstellmass(self,plotsingle=1):
        if plotsingle:
            figure()
            title(self.prefix)
            xlabel('$ R $',fontsize=20)
            ylabel('$ L_{IR} $',fontsize=20)
        distance=sqrt((self.gdat.RA-self.cRA)**2+(self.gdat.Dec-self.cDec)**2)
        keepflag=(abs(self.zcluster - self.gdat.zLDP) < .03) &( distance < 1.4*7./60)
        #plot(self.gdat.Rauto[keepflag],self.gdat.LIR_ZCLUST[keepflag],'ko',mfc='w',mec='k',markersize=10,label='_nolegend_')
        n1=sum(keepflag & (self.gdat.LIR_ZCLUST>1.e11)& (self.gdat.Rauto < 23.))

        keepflag=self.rdat.matchflag24 & self.rdat.supermembflag & (self.rdat.newspecmatchflag < .1)
        plot(self.rdat.stellmass[keepflag],self.rdat.SFRir[keepflag],'ko',mfc='w',markersize=10,label='EDisCS Photoz')
        n2=sum(keepflag & (self.rdat.Lir>1.e11) & (self.rdat.magR < 23.))
        keepflag=self.rdat.matchflag24 & (self.rdat.newspecmatchflag > 0.1)
        plot(self.rdat.stellmass[keepflag],self.rdat.SFRir[keepflag],'k^',markersize=8,label='EDisCS Spec-z')
        keepflag=(self.rdat.supermembflag) & (self.rdat.sfflag > .1) & (self.rdat.stellmass < 1.e11)
        #plot(self.rdat.stellmass[keepflag],self.rdat.SFR[keepflag],'gs',markersize=8,label='EDisCS Halpha')

        gca().set_yscale('log')
        gca().set_xscale('log')
        #axhline(y=10.**11,color='k',ls='--',lw=1)
        #axvline(x=23.,color='k',ls='--',lw=1)
        #axis([20.5,24.5,5.e10,2.e12])
        #xticks(arange(21,25))
        if self.prefix.find('cl1216') > -1:
            legend(loc='upper left',numpoints=1,prop={'size':12})
        #s=fullname[self.prefix]+'\nz=%4.3f \nN=%2i'%(self.zcluster,n1+n2)
        #text(.1,.95,s,verticalalignment='top',horizontalalignment='left',transform=gca().transAxes,fontsize=12)
        #s='N Target=%i'%(n1+n2)
        #text(.1,.8,s,horizontalalignment='left',transform=gca().transAxes,fontsize=12)
        print self.prefix,': number of targets with log10(L_IR) > 11 = ',n1+n2

    def plotsSFRstellmass(self,plotsingle=1):
        if plotsingle:
            figure()
            title(self.prefix)
            xlabel('$ R $',fontsize=20)
            ylabel('$ L_{IR} $',fontsize=20)

        distance=sqrt((self.gdat.RA-self.cRA)**2+(self.gdat.Dec-self.cDec)**2)
        keepflag=(abs(self.zcluster - self.gdat.zLDP) < .03) &( distance < 1.4*7./60)
        #plot(self.gdat.Rauto[keepflag],self.gdat.LIR_ZCLUST[keepflag],'ko',mfc='w',mec='k',markersize=10,label='_nolegend_')
        n1=sum(keepflag & (self.gdat.LIR_ZCLUST>1.e11)& (self.gdat.Rauto < 23.))

        keepflag=self.rdat.matchflag24 & self.rdat.supermembflag & (self.rdat.newspecmatchflag < .1)
        plot(self.rdat.stellmass[keepflag],self.rdat.SFRir[keepflag],'ko',mfc='w',markersize=10,label='EDisCS Photoz')
        n2=sum(keepflag & (self.rdat.Lir>1.e11) & (self.rdat.magR < 23.))
        keepflag=self.rdat.matchflag24 & (self.rdat.newspecmatchflag > 0.1)
        plot(self.rdat.stellmass[keepflag],self.rdat.SFRir[keepflag],'k^',markersize=8,label='EDisCS Spec-z')
        keepflag=(self.rdat.supermembflag) & (self.rdat.sfflag > .1) & (self.rdat.stellmass < 1.e11)
        #plot(self.rdat.stellmass[keepflag],self.rdat.SFR[keepflag],'gs',markersize=8,label='EDisCS Halpha')

        gca().set_yscale('log')
        gca().set_xscale('log')
        #axhline(y=10.**11,color='k',ls='--',lw=1)
        #axvline(x=23.,color='k',ls='--',lw=1)
        #axis([20.5,24.5,5.e10,2.e12])
        #xticks(arange(21,25))
        if self.prefix.find('cl1216') > -1:
            legend(loc='upper left',numpoints=1,prop={'size':12})
        #s=fullname[self.prefix]+'\nz=%4.3f \nN=%2i'%(self.zcluster,n1+n2)
        #text(.1,.95,s,verticalalignment='top',horizontalalignment='left',transform=gca().transAxes,fontsize=12)
        #s='N Target=%i'%(n1+n2)
        #text(.1,.8,s,horizontalalignment='left',transform=gca().transAxes,fontsize=12)
        print self.prefix,': number of targets with log10(L_IR) > 11 = ',n1+n2

    def plotpositions(self,plotsingle=1):
        if plotsingle:
            figure()
            title(self.prefix)
            xlabel('$ RA $',fontsize=20)
            ylabel('$ Dec $',fontsize=20)

        #cir=Circle((self.clusterra-racenter,self.clusterdec-deccenter),radius=1.3*self.r200deg,color='None',ec='0.2',alpha=0.2)
        cir=Circle((self.cRA,self.cDec),radius=1.*self.r200deg,color='None',ec='0.2',alpha=0.2)
        gca().add_patch(cir)

        #keepflag=(abs(self.zcluster - self.gdat.zLDP) < .03) &( distance < 1.4*7./60)
        distance=sqrt((self.rdat.RA-self.cRA)**2+(self.rdat.DEC-self.cDec)**2)
        keepflag=self.rdat.matchflag24 & self.rdat.supermembflag & (self.rdat.newspecmatchflag < .1)
        keepflag=self.rdat.matchflag24 &  (abs(self.zcluster - self.rdat.bestz) < .2) & ( distance < 2.*self.r200deg) & (self.rdat.magR < 23.) & (self.rdat.newspecmatchflag < 0.1)
        plot(self.rdat.RA[keepflag],self.rdat.DEC[keepflag],'r.',mec='r',markersize=10,label='EDisCS Photoz')

        keepflag=self.rdat.matchflag24 & (self.rdat.newspecmatchflag > 0.1)
        plot(self.rdat.RA[keepflag],self.rdat.DEC[keepflag],'k^',markersize=8,label='EDisCS Spec-z')

        keepflag=(self.rdat.newspecmatchflag > 0.1)
        plot(self.rdat.RA[keepflag],self.rdat.DEC[keepflag],'k.',markersize=8,label='_nolegend_')

        keepflag=self.rdat.matchflag24 & ( distance < 2.*self.r200deg) & (abs(self.zcluster - self.rdat.bestz) > .2) & (self.rdat.magR < 23.) & (self.rdat.newspecmatchflag < 0.1)
        plot(self.rdat.RA[keepflag],self.rdat.DEC[keepflag],'ro',mfc='w',mec='r',markersize=10,label='EDisCS Photoz')

        axis('equal')
        ax=gca()
        text(.5,.9,fullname[self.prefix],transform=ax.transAxes,horizontalalignment='center')
        #plot(self.rdat.RA[keepflag],self.rdat.DEC[keepflag],'rs',markersize=8,label='MIPS Spec-z')     



cl1216=cluster('cl1216')
cl1354=cluster('cl1354')
cl105412=cluster('cl105412')
cl1040=cluster('cl1040')
cl105411=cluster('cl105411')
cl1227=cluster('cl1227')
cl1353=cluster('cl1353')
cl1037=cluster('cl1037')

eclusters=[cl1216,cl1354,cl105412,cl1040,cl105411,cl1227,cl1353,cl1037]

def makeplot():
    figure(figsize=(10,8))
    subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.95,wspace=0.001,hspace=0.001)
    i=0
    for cl in eclusters:
        subplot(2,3,i+1)
        cl.plotLIRI(plotsingle=0)
        i+=1
        ax=gca()
        if i < 4:
            ax.set_xticklabels(([]))
            if i > 1:
                ax.set_yticklabels(([]))
        if i > 4:
            ax.set_yticklabels(([]))
    text(.5,-.2,'$ R $',horizontalalignment='center',transform=ax.transAxes,fontsize=20)
    text(-1.3,1,'$ L_{IR}/L_\odot $',verticalalignment='center',rotation=90,transform=ax.transAxes,fontsize=20)
    savefig(homedir+'proposals/observing/ESO2014Spring/LIRvsI.pdf')

def plotsfrstellmass():
    figure(figsize=(10,8))
    #subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.95,wspace=0.001,hspace=0.001)
    i=0
    for cl in eclusters:
        cl.plotSFRstellmass(plotsingle=0)
        i+=1
        ax=gca()
        if i < (len(eclusters)-1):
            ax.set_xticklabels(([]))
            ax.set_yticklabels(([]))
    axis([2.e9,1.e12,.01,300])
    text(.5,-.1,'$ Stellar \ Mass \ (M_\odot) $',horizontalalignment='center',transform=ax.transAxes,fontsize=20)
    text(-.1,.5,'$ SFR_{IR} \ (M_\odot/yr) $',verticalalignment='center',rotation=90,transform=ax.transAxes,fontsize=20)
    s.plotSFRvsStellarmassv2(plotsingle=0)
    #except:
    #    print 'LCS not loaded'
    #    print '%run ~/svnRepository/pythonCode/LCSanalyzespirals.py'
    xl=arange(2.e9,2.e12,5.e9)
    yl=xl/6.e10
    plot(xl,yl,'k',label='z=0 (Finn+ in prep)')
    plot(xl,7.*yl,'k--',label='7x (Ly+2011)')
    plot(xl,30.*yl,'k:',label='30x')
    legend()
    handles,labels=ax.get_legend_handles_labels()
    legend(handles[-6:len(labels)],labels[-6:len(labels)],numpoints=1,loc='lower right',prop={'size':12})
    plot(xl,20.*yl,'g:',label='30x')
    plot(xl,40.*yl,'g:',label='30x')
    savefig(homedir+'proposals/observing/ESO2014Spring/SFRstellmass.pdf')

def plotpositions():
    figure(figsize=(12,4))
    
    subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.95,wspace=0.2,hspace=0.001)
    subplot(1,3,1)
    cl1040.plotpositions(plotsingle=0)
    subplot(1,3,2)
    cl105411.plotpositions(plotsingle=0)
    subplot(1,3,3)
    cl105412.plotpositions(plotsingle=0)
    savefig(homedir+'proposals/observing/ESO2014Spring/positions.pdf')
