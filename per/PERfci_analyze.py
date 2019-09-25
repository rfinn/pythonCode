#!/usr/bin/env python

'''
written by Rose A. Finn


'''

import pylab as pl
import xlrd
#import mystuff as my
from scipy.stats.stats import spearmanr
from scipy.stats import ks_2samp

max_score=30

def spearman(x,y):
    rho,pvalue=spearmanr(x,y)
    print 'Spearman Rank Test:'
    print 'rho = %6.2f'%(rho)
    print 'p-vale = %6.5f (prob that samples are uncorrelated)'%(pvalue)
    return rho,pvalue

def ks(x,y):
    D,pvalue=ks_2samp(x,y)
    print 'KS Test:'
    print 'D = %6.2f'%(D)
    print 'p-vale = %6.5f (prob that samples are from same distribution)'%(pvalue)
    return D,pvalue


def calc_gini(x):
    x=list(x)
    N = len(x)
    x.sort()
    G = sum( x[i]*(N-i) for i in xrange(N))
    G = 2.*G/(N*sum(x))
    return (1 + 1./N) - G

def read_data(file):
    xdat=xlrd.open_workbook(file)
    sheet=xdat.sheet_by_index(0)
    section=sheet.col_values(0)[1:]
    id=sheet.col_values(1)[1:]
    pretest=pl.array(sheet.col_values(2)[1:],'f')
    posttest=pl.array(sheet.col_values(3)[1:],'f')
    lawson=pl.array(sheet.col_values(4)[1:],'f')
    medsker_gain=pl.array(sheet.col_values(5)[1:],'f')
    return section,id,pretest,posttest,lawson,medsker_gain

infile='/Users/rfinn/Copy/PER/SummaryFCIDataLM20Feb2014.xlsx'
section,id,pretest,posttest,lawson,medsker_gain = read_data(infile)


#dictionary to connect section number to time/day

class classdata():
    def __init__(self,classid):
        matchindex=[i for i,s in enumerate(section) if classid in s]
        self.pretest=pretest[matchindex]
        self.posttest=posttest[matchindex]
        self.lawson=posttest[matchindex]
        self.medsker_gain=medsker_gain[matchindex]
        self.classid=classid
        self.calc_gain()
        self.calc_G()
    def calc_gain(self):

        self.gain=(self.posttest-self.pretest)/(max_score - self.pretest)
        self.ave_gain = pl.mean(self.gain)
        self.std_gain=pl.std(self.gain)/pl.sqrt(1.*len(self.gain))
        # calculate gini coefficient
        self.gini=calc_gini(self.gain)
    def calc_G(self):
        self.ave_pretest=pl.mean(self.pretest)
        self.ave_posttest=pl.mean(self.posttest)
        self.G=(self.ave_posttest - self.ave_pretest)/(max_score - self.ave_pretest)

    def pre_hist(self):
        mybins=pl.arange(0,34,4)
        pl.hist(self.pretest,bins=mybins,color='0.5')
        pl.xlim(0,32)
        pl.ylim(0,12)
        ax=pl.gca()
        pl.text(0.9,0.8,self.classid,horizontalalignment='center',transform=ax.transAxes,fontsize=14)
        pl.axvline(x=pl.mean(self.pretest),color='k',ls='--')

    def post_hist(self):
        mybins=pl.arange(0,34,4)
        pl.hist(self.posttest,bins=mybins,color='0.5')
        pl.xlim(0,32)
        pl.ylim(0,7)
        ax=pl.gca()
        pl.text(0.9,0.8,self.classid,horizontalalignment='center',transform=ax.transAxes,fontsize=14)
        pl.axvline(x=pl.mean(self.posttest),color='k',ls='--')
    def gain_hist(self):
        mybins=pl.arange(-.4,1.3,.15)
        pl.hist(self.gain,bins=mybins,color='0.5')
        pl.xlim(-.4,1.1)
        pl.ylim(0,8)
        ax=pl.gca()
        pl.text(0.9,0.8,self.classid,horizontalalignment='center',transform=ax.transAxes,fontsize=14)
        pl.axvline(x=pl.mean(self.gain),color='k',ls='--')


course_sections=['CA1','CA2','CA3','UP1','UP2','UP3']
ca1=classdata(course_sections[0])
ca2=classdata(course_sections[1])
ca3=classdata(course_sections[2])
up1=classdata(course_sections[3])
up2=classdata(course_sections[4])
up3=classdata(course_sections[5])

courses=[ca1,ca2,ca3,up1,up2,up3]

def plot_whiskers():
    pl.figure()
    pl.clf()
    bp=pl.boxplot([ca1.gain,ca2.gain,ca3.gain,up1.gain,up2.gain,up3.gain],vert=False)#,bootstrap=1000)
    pl.setp(bp['boxes'],color='black')
    pl.setp(bp['whiskers'],color='black')
    pl.setp(bp['fliers'],color='black')
    pl.setp(bp['medians'],color='0.5')
    pl.xlabel('Normalized Gain',fontsize=14)
    pl.ylabel('Class Section',fontsize=14)
    pl.xlim(-.4,1.05)
    a,b=pl.yticks()
    pl.yticks(a,course_sections)
    pl.savefig('whisker_plot.png')
    pl.savefig('whisker_plot.eps')
def gain_hist():
    pl.figure(figsize=(12,7))
    i=1
    mybins=pl.arange(-.4,1.3,.15)
    print mybins
    pl.subplots_adjust(left=.1,right=0.95,wspace=.25,hspace=0.2)
    for c in courses:
        pl.subplot(2,3,i)
        pl.hist(c.gain,bins=mybins,color='k')
        pl.xlim(-.4,1.1)
        pl.ylim(0,8)
        ax=pl.gca()
        pl.text(0.9,0.9,c.classid,horizontalalignment='center',transform=ax.transAxes,fontsize=14)
        pl.axvline(x=pl.mean(c.gain),color='k',ls='--')
        i+= 1

    pl.text(-.7,-.2,'Normalized Gain',transform=ax.transAxes,horizontalalignment='center',fontsize=20)
    pl.savefig('gain_hist.png')
    pl.savefig('gain_hist.eps')

def post_hist():
    pl.figure(figsize=(12,7))
    i=1
    mybins=pl.arange(0,34,4)
    print mybins
    pl.subplots_adjust(left=.1,right=0.95,wspace=.25,hspace=0.2)
    for c in courses:
        pl.subplot(2,3,i)
        pl.hist(c.posttest,bins=mybins)
        pl.xlim(0,32)
        pl.ylim(0,7)
        ax=pl.gca()
        pl.text(0.9,0.9,c.classid,horizontalalignment='center',transform=ax.transAxes,fontsize=14)
        pl.axvline(x=pl.mean(c.posttest),color='k',ls='--')
        i+= 1

    pl.text(-.7,-.2,'Post-Test Score',transform=ax.transAxes,horizontalalignment='center',fontsize=20)
    pl.savefig('posttest_hist.png')
    pl.savefig('posttest_hist.eps')

def allhists():
    pl.figure(figsize=(10,12))
    pl.subplots_adjust(left=.05,right=0.95,wspace=.25,hspace=0.25,top=0.95)
    i=1
    for c in courses:

        pl.subplot(6,3,i)
        c.pre_hist()
        i+= 1
        pl.subplot(6,3,i)
        c.post_hist()
        i+=1
        pl.subplot(6,3,i)
        c.gain_hist()
        i+=1
    pl.subplot(6,3,16)
    pl.xlabel('Pre-test Score',fontsize=16)
    pl.subplot(6,3,17)
    pl.xlabel('Post-test Score',fontsize=16)
    pl.subplot(6,3,18)
    pl.xlabel('Normalized Gain',fontsize=16)
    pl.savefig('allhists.png')
    pl.savefig('allhists.eps')

def compare_sections():
    group1=ca1.gain.tolist()+ca2.gain.tolist()

    group2=up1.gain.tolist()+up2.gain.tolist()+up3.gain.tolist()

    print 'Compare CA1-CA2 to UP sections'
    ks(group1,group2)
    print '****************** \n'
    print 'Compare CA3 to CA1-CA2 '
    ks(ca3.gain,group1)
    print '****************** \n'
    print 'Compare CA3 to UP sections'
    ks(ca3.gain,group2)
    print '****************** \n'
    print 'Compare UP1 to UP2-UP3'
    ks(up1.gain,up2.gain.tolist()+up3.gain.tolist())
    print '****************** \n'
