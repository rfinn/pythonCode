#!/usr/bin/env python
import sqlcl

def getsdsscatalogs():
    drsearch=5.*60.#search radius in arcmin for sdss query
    zmin=.0
    zmax=6.
    ra=180.
    dec=0.
	#print i,cid[i]," ra, dec, dr, mr = %12.8f %12.8f %8.3f %5.2f" % (cra[i],cdec[i],drsearch)
    query="select n.distance,g.ra,g.dec, g.u, g.g, g.r, g.i, g.z, s.z from galaxy g, specobj s, dbo.fGetNearbyObjEq(%12.8f,%12.8f,%8.3f) n where g.objid = s.bestobjid and g.objID = n.objID and s.z < %5.4f and s.z > %5.4f and (g.PrimTarget & 0x00000040) > 0 order by distance" % (ra,dec,drsearch,zmax,zmin)

    print query
    try:
        lines=sqlcl.query(query).readlines()
    except IOError:
        print "IOError for cluster trying spec query again"
        lines=sqlcl.query(query).readlines()
    print "got number + 1 of spec objects = ",len(lines)
    n='JonGalaxy.dat'
    outfile=open(n,'w')
    j=0
    if (len(lines) > 1.):
        for line in lines[1:]:
            if j < 0:
                print line
                j=j+1
            outfile.write(line)
    outfile.close()

getsdsscatalogs()
