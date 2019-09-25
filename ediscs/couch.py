#!/usr/bin/env python
omegaL=0.7
omega0=0.3
snmin=2
c=3.e5
sfrmin=0.#min sfr in solar mass/yr
#f200=.75
#f200=.79
f200=.5
snmin=3
ewmin=20#min ew in Angstrom
dm=0.5
balogh()
balogh02()
couch()
#postman01()


os.system("cp *totalsfr /home/rfinn/clusters/final/.")



def getR200():
  R200=1.73*sig/1000./sqrt(omegaL+omega0*(1+z)**3) #Mpc/h
  R200=R200*1000.  #convert to kpc
  print "R200 (kpc) = R200 \n"  
  R200=R200/dA#convert to arcsec
  print "R200 (arcsec) = R200 \n"
  R200=R200*f200

def couch():
    IN1=open("/home/rfinn/clusters/final/lit_data/couchdatafile1.txt",'r') 
    output="couchtotalsfr"
    open(OUT100,">output") || die "can't open outfile !"
    #  sig=1660 #mahdavi and geller
    sig=1390
    errsigp=200
    errsigm=200
    ra_c=15*(22+58./60.+52.3/3600.)#RA center from NED,convert to deg
    dec_c=-1*(34+46/60. +55./3600.)#DEC center from NED
    i=0
    sfrtot=0
    sfrtot2=0
    nel=0
    ntot=0
    fel=0
    errfel=0
    zcl=0.32
    dv=4*sig/c#apply +/- 3sigma cut in velocity
    dz=dv*(1+zcl)
    zmin=zcl-dz
    zmax=zcl+dz
    dL=1172. #Mpc/h
    dA=3.26 #kpc/arcsec
    getR200()
    mImax=3.21+5*log(dL)/log(10)#see notebook for calculation of limiting I mag
    mImax=mImax+dm
    while (line = <IN1>):
        chomp(line)
        if (i > 33):
            z=substr(line,52,7)
            if((z > zmin) & (z < zmax)):
                rah=substr(line,5,2)
                ram=substr(line,8,2)
                ras=substr(line,11,7)
                signdec=substr(line,19,1)
                decd=substr(line,19,3)
                decm=substr(line,23,2)
                decs=substr(line,26,6)
                ra=15*(rah+ram/60. + ras/3600.)
                dec=-1*(-1*decd+decm/60. + decs/3600.)
                d=sqrt((ra_c - ra)**2+(dec_c - dec)**2)*3600.
                sfr=substr(line,84,7)
                f=substr(line,60,7)
                errf=substr(line,68,7)
                sfr=f*((dL)**2)*9.48*10**(-9)
                errsfr=errf*((dL)**2)*9.48*10**(-9)
                sfr=sfr*2.5#correct for dust
                errsfr=errsfr*2.5
                sn = abs(f/errf)
                sfr2=substr(line,84,7)
                sfr3=4*sfr
                mI=substr(line,33,5)
                if (d < R200):
                    if(sfr > sfrmin & sn > snmin):
                        sfrtot=sfrtot+sfr
                        errsfrtot=errsfrtot+errsfr**2
                        if (mI < mImax):
                            ntot++
                            if(sfr > sfrmin && sn > snmin):
                                nel++
                                
    i+=1
    close(IN1)
    errsfrtot=sqrt(errsfrtot)
    print "Total SFR (Ho=100) = sfrtot \n"
    errfel=sqrt(nel+nel/ntot**2)/ntot
    fel=nel/ntot
    print "Fraction of Emission Line galaxies = fel \n"
    writeoutput()
    close(OUT100)

def balogh():
    open(IN1,"/home/rfinn/clusters/final/lit_data/balogh00table2.dat") || die "can't open balogh00table2.dat"
    output="baloghmorristotalsfr"
    open(OUT100,">output") || die "can't open outfile !"
    print "Balogh & Morris \n"
    i=0
    nel
    fel=0
    errfel=0
    ntot=0
    sfrtot=0
    errsfrtot=0
    
    sig=1023 #mahdavi and geller
    errsigp=200
    errsigm=200
    zcl=0.23
    const=0.152#see notebook re conversion, from flux to sfr
    
    dL=795#Mpc
    dA=2.56 #kpc/arcsec
    getR200()
    fmin=sfrmin/dL**2/0.001/const
    fmin=fmin*6400*0.23/ewmin
    mRmax=3.86+5*log(dL)/log(10)
    mRmax=mRmax+dm
    print "fmin = fmin\n"
    print "R max = mRmax \n"
    while (line = <IN1>):
        chomp(line)
        l=length(line)
        if (l > 110):
            member=substr(line,113,3)
            ew=substr(line,87,6)
            errew=substr(line,94,4)
            f=substr(line,72,7) #flux in mJy
            errf=substr(line,80,6)#errflux in mJy
            sn=abs(f/errf)
            sfr=f*.001*dL**2*const
            errsfr=errf*.001*dL**2*const      
            sfr=sfr*2.5#correct for dust
            errsfr=errsfr*2.5
            dra=substr(line,37,7)
            ddec=substr(line,45,6)
            d=sqrt(dra**2 + ddec**2)
            mR=substr(line,52,6)
            if (member =~ /y/):
                if (d < R200):
                    if(sfr > sfrmin && sn > snmin):
                        sfrtot=sfrtot+sfr
                        errsfrtot=errsfrtot+errsfr**2
                        if (mR < mRmax):
                            ntot++
                            if(sfr > sfrmin && sn > snmin):
                                nel++
    i++
    close(IN1)
    errsfrtot=sqrt(errsfrtot)
    errfel=sqrt(nel+nel/ntot**2)/ntot
    fel=nel/ntot
    print "Fraction of Emission Line galaxies = fel \n"
    writeoutput()
    close(OUT100)

sub balogh02{
  #H=100, omega_m=0.3, omega_L=0.7 - yeah!!
  open(IN1,"/home/rfinn/clusters/final/lit_data/balogh02table1.dat") || die "can't open balogh00table1.dat"
  output="balogh02totalsfr"
  open(OUT100,">output") || die "can't open outfile !"
  i=0
  fel=0
  nel
  ntot=0
  sfrtot=0
  errtot=0
#  sig=1819
  sig=1273
  errsigp=200
  errsigm=200
  zcl=0.183
  dv=3*sig/c
  dz=dv*(1+zcl)
  zmin=zcl-dz
  zmax=zcl+dz
  ra_c=15*(13+11./60.+34.2/3600.)#RA center from NED,convert to deg
  dec_c=-1*(1+21/60. +56./3600.)#DEC center from NED

  dL=622. #Mpc/h
  dA=2.15 #kpc/arcsec
  mImax=3.21+5*log(dL)/log(10)#see notebook for calculation of limiting I mag
  mImax=mImax+dm
  print "Max I mag = mImax \n"
  getR200()
  #print "z zmin=zmin zmax=zmax \n"
  i=0
  while (line = <IN1>){
    chomp(line)
    z=substr(line,39,5)
    sfr=substr(line,65,5)
    errsfr=substr(line,71,5)
    
    #print "i z sfr\n"
    #z = 0.159 - 0.206.

    if ((z >= zmin) && (z <= zmax)){
      #print "hey\n"
      rah=substr(line,5,2)
      ram=substr(line,8,2)
      ras=substr(line,11,6)
      #signdec=substr(line,18,1)
      decd=substr(line,18,2)
      decm=substr(line,21,2)
      decs=substr(line,24,5)
      decd=join "",(signdec,decd)
      ra=15*(rah+ram/60. + ras/3600.)
      dec=-1*(-1*decd+decm/60. + decs/3600.)
      #dec=join "",(signdec,dec)
      #print "ra dec decd ra_c dec_c \n"
      d=sqrt((ra_c - ra)**2+(dec_c - dec)**2)*3600.
      mI=substr(line,30,5)
      #print "mI mImax \n"
      if (errsfr > 0){
	sn=abs(sfr/errsfr)
      }else{
	sn=0
      }
      if (d < R200){
	if(sfr > sfrmin && sn > snmin){
	  sfrtot=sfrtot+sfr
	  errsfrtot=errsfrtot+errsfr**2
	}
	if (mI < mImax){
	  ntot++
	  if(sfr > sfrmin && sn > snmin){
	    nel++
	  }
	}
      }
      
    }
    i++
  }
  close(IN1)
  errsfrtot=sqrt(errtot)
  errfel=sqrt(nel+nel/ntot**2)/ntot
  fel=nel/ntot
  print "Fraction of Emission Line galaxies = fel \n"
  print "Integrated SFR = sfrtot +/- errtot \n"
  printf "nel ntot %8.3f %8.3f \n",fel, errfel
  writeoutput()
  close(OUT100)
}


sub writeoutput{
  #account for 30% uncertainty in halpha - SFR conversion
  errsfrtot=sqrt(errsfrtot**2 + (0.3*sfrtot)**2)
  #open(OUT100,">output") || die "can't open outfile !"
  printf OUT100 "%8.2f %8.2f %8.2f %8.2f %8.3f %8.3f sig errsigp errsigm zcl\n",sfrtot,errsfrtot,nel,ntot,fel,errfel
  printf "%8.2f %8.2f %8.2f %8.2f %8.3f %8.3f sig errsigp errsigm zcl\n",sfrtot,errsfrtot,nel,ntot,fel,errfel
  #close(OUT100)
}
