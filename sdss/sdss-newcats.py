#!/usr/bin/env python
#read in cluster files from mike balogh.
#match RA, Dec, z of each galaxy with dr2 catalogs
#write one new file for all cluster galaxies containing:
#ra, dec, z, EW(Ha), errEW(Ha), EW(OII), errEW(OII), AGN flags from Miller, line# in dr files
import sys, glob
import Numeric as N
import scipy
from math import *
import ppgplot
import mx.DateTime,time
import mystuff as my

starttime=time.clock()
print "start time = ",starttime

clusters = glob.glob("*_DR*")#get list of cluster/field catalogs
ra=[]
dec=[]
z=[]
clusterid=[]
id=[]
ewha=[]
errewha=[]
agnflag1=[]
agnflag2=[]
agnflag3=[]
catra=[]
catdec=[]
catz=[]
catline=[]
catewha=[]
caterrewha=[]
nomatch=0
print "Reading DR files"
for i in range(len(clusters)):
    prefix=clusters[i].split('.')
    j=0

    file=open(clusters[i],'r')
    for line in open(clusters[i]):
        if line.find('#') > -1:
            continue
        if j == 0:
            fields=line.split()
            id.append(int(fields[0]))
            tid=float(fields[0])
            j=1
            #print "cluster id = ",fields[0]
            continue
        if line.find('*') > -1: #skip any galaxies w/ ****** in one or more fields
            continue
        fields=line.split()
        ra.append(float(fields[0]))
        dec.append(float(fields[1]))
        z.append(float(fields[2]))
        clusterid.append((tid))
print "Reading Ha EW file"
file=open("/home/rfinn/SDSS/cluster_catalogs/dr2/LINES/sdss_dr2_spec6564_61.dat",'r')
i=0
for line in file:
    if line.find('#') > -1:
        continue
    fields=line.split()
    catra.append(float(fields[0]))
    catdec.append(float(fields[1]))
    catz.append(float(fields[2]))
    catline.append(i)
    catewha.append(float(fields[5]))
    caterrewha.append(float(fields[6]))
    i += 1
file.close()
#convert to arrays
ra=N.array(ra,'f')
dec=N.array(dec,'f')
z=N.array(z,'f')
clusterid=N.array(clusterid,'f')
id=N.array(id,'f')
catra=N.array(catra,'f')
catdec=N.array(catdec,'f')
catz=N.array(catz,'f')
catline=N.array(catline,'i')
catewha=N.array(catewha,'f')
caterrewha=N.array(caterrewha,'f')
index = N.zeros(len(ra),'i')
print "number of cluster galaxies = ",len(ra)
ntime=0
print "sorting catalog ra array"
(catrasort,catraindex)=my.sortwindex(catra)
deltara=.01#matching tolerance in ra in arcsec
deltara=deltara/3600.*15.#convert to degrees
print "matching cluster galaxies against dr2 catalog"
endloop=len(ra)
for i in range(len(ra)):
    if 1.*i/500. > ntime:
        print i, "galaxies", mx.DateTime.localtime(), "nomatch = ",nomatch
        ntime += 1
    match=my.findmatch(ra[i],catrasort,deltara)#returns indices of catrasort that lie w/in deltara
    #print "matched",i
    try:#if more than one match, find closest in terms of ra,dec,z
        if len(match) > 1:
            diffmin=1000000.
            for j in range(len(match)):
                catindex=catraindex[match[j]]#catindex is for unsorted cat arrays
                diff=N.sqrt((catra[catindex]-ra[i])**2+(catdec[catindex]-dec[i])**2+(catz[catindex]-z[i])**2)
                if diff < diffmin:
                    diffmin=diff
                    catmatchindex=catindex
            index[i]=catmatchindex
            ewha.append(catewha[index[i]])
            errewha.append(caterrewha[index[i]])
    except:#match is number, not list
        if match > 0:
            catmatchindex=catraindex[match]
            index[i]=catmatchindex
            ewha.append(catewha[index[i]])
            errewha.append(caterrewha[index[i]])
        else:
            catmatchindex=-999.
            index[i]=-999.
            ewha.append(-999.)
            errewha.append(-999.)
            nomatch=nomatch+1


print "Reading AGN file"
catagn1=[]
catagn2=[]
catagn3=[]
file=open("/home/rfinn/SDSS/cluster_catalogs/dr2/sdss_dr2_deriv_agnflag.dat",'r')
i=0
for line in file:
    if line.find('#') > -1:
        continue
    fields=line.split()
    catagn1.append(float(fields[3]))
    catagn2.append(float(fields[4]))
    catagn3.append(float(fields[5]))
file.close()

print "writing output, mycat"
output=open('mycat','w')

nomatchvalue=-999.
for i in range(len(ra)):    
#for i in range(10):
    #print >> output,'%16.10f %16.10f %16.10f %16.10f %16.10f %i' % (ra[i],dec[i],z[i],ewha[i],errewha[i],agn[index[i]],index[i])
    if index[i] < 0:
        print >> output, '%16.10f %16.10f %16.10f %16.10f %16.10f %16.10f %16.10f %16.10f %i' % (ra[i],dec[i],z[i],ewha[i],errewha[i],nomatchvalue,nomatchvalue,nomatchvalue,index[i])
    else:
        print >> output, '%16.10f %16.10f %16.10f %16.10f %16.10f %16.10f %16.10f %16.10f %i' % (ra[i],dec[i],z[i],ewha[i],errewha[i],catagn1[index[i]],catagn2[index[i]],catagn3[index[i]],index[i])
    #print >> output,'5%16.10f %i' % (ra[i],dec[i],z[i],ewha[i],errewha[i],index[i])
output.close()    
    
endtime=time.clock()    
print "end time = ",endtime
print "elapsed time = ",endtime-starttime
