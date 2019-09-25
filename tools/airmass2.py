#!/usr/bin/env python
#started to transform to python, but can't get skycalendar to work
import ppgplot
import sys, glob
import Numeric as N
from math import *
import time, os

try:
    print "sys.argv"
except:
    printf "pgam [FILE] [YEAR] [MONTH] [DAY]\n"
    

# Constant
#################################################################
pi=3.1415

# Mauna kea (UKIRT)
#longitude='-155:28:23.6'  
#latitude='19:49:32.2'  

#MMT
#latitude = 31.6883
#longitude = 110.885
latitude = '31:41:19.6'
longitude ='110:53:04.4'

lambda=dec_to_deg(longitude) * pi/180      #!<----should be "dec_to_deg" 

phai=dec_to_deg(latitude) * pi/180

a_min=1
a_max=2.00001

# time

t_min=-8
t_max=8
t_delta=0.01

t_min_curve=-30
t_max_curve=30

# Plot A (Airmass-Time)
########################
# PS OUTPUT 
###############################################################
filename=sys.argv[1]
psfile = str(filename)+".ps" # print ("psfile\n")
ppgplot.pgbegin(0,"psfile/VCPS", 1, 1)
# pgbegin(0,"psfile/PS", 1, 1)

# Plot Setting
####################################################################
ppgplot.pgpaper(8,1.25) # window/paper size (width(inch), aspect)

ppgplot.pgscf(2)   # characte font (1: normal, 2: roman, 3: italic, 4: script)
ppgplot.pgslw(3) # line width
ppgplot.pgsvp(0.15, 0.9, 0.53, 0.89) # viewport in the window (relative)
ppgplot.pglab("", "", "Local Time [hour]")
ppgplot.pgsvp(0.12, 0.9, 0.53, 0.88) # viewport in the window (relative)
ppgplot.pglabel("", "Airmass", "") # label settingoto s

ppgplot.pgsch(1.0) # character height (size)
ppgplot.pgslw(3) # line width
ppgplot.pgsvp(0.15, 0.9, 0.53, 0.88) # viewport in the window (relative)
ppgplot.pgswin(t_min, t_max, a_max, a_min) # MIN,MAX of coordinate

ppgplot.pgbox('BCTS', 0.0, 0, 'BCTSNV1', 0.1, 0) # coordinate settings
ppgplot.pgbox('0', 0.0, 0, 'BCTSMV1', 0.1, 0) # coordinate settings

# Put Header/ Axes Label
#####################################################################
#####################################################################
# READ USER INPUT
#####################################################################
year  =  sys.argv[2]#ARGV[1]
month =  sys.argv[3]#ARGV[2]
day   =  sys.argv[4]#ARGV[3]

# print "TODAY   year month day\n "

# Read LMST 
########################################################

#calendar_buf =  `echo m 1 year | skycalendar`
#calendar_buf =~ m/month day\//g
#(calendar) = split (/\n/, ')

(buf, buf,buf, jdmid, hour,min,sec,ss_hour, ss_min, twie_hour,twie_min, twib_hour, twib_min, sr_hour, sr_min) = calendar.split()#/\s+/, calendar)
ppgplot.pgsch(0.8) # character height (size)
ppgplot.pgslw(2) # line width

ppgplot.pgtext(6.5,0.80,filename)

ppgplot.pgsch(0.6) # character height (size)
ppgplot.pgtext(-10.32,0.80,str(year))
ppgplot.pgtext(-9.5,0.80,str(month))
ppgplot.pgtext(-8.9,0.80,str(day))

ppgplot.pgtext(-8,0.80,"LMST Midnight")
ppgplot.pgtext(-4.8,0.80,"hour")
ppgplot.pgtext(-4.25,0.80,"min")
ppgplot.pgtext(-3.75,0.80,"sec")

ppgplot.pgtext(-8,0.84,"Sun Set")
ppgplot.pgtext(-4.8,0.84,"ss_hour")
ppgplot.pgtext(-4.25,0.84,"ss_min")

ppgplot.pgtext(-8,0.88,"Twilight End")
ppgplot.pgtext(-4.8,0.88,"twie_hour")
ppgplot.pgtext(-4.25,0.88,"twie_min")

ppgplot.pgtext(3.8,0.84,"Sun Rise")
ppgplot.pgtext(7.2,0.84,"sr_hour")
ppgplot.pgtext(7.55,0.84,"sr_min")

ppgplot.pgtext(3.8,0.88,"Twilight BEGIN")
ppgplot.pgtext(7.2,0.88,"twib_hour")
ppgplot.pgtext(7.55,0.88,"twib_min")


t_lmstmidn = hour + min/60 + sec/3600
twie_dn = twie_hour - 24 + twie_min/60 
twib_dn = twib_hour + twib_min/60 

ss_dn = ss_hour - 24 + ss_min/60 
sr_dn = sr_hour + sr_min/60 

# SHADE Twilight zone
#####################################################################
ppgplot.pgsfs(1)

ppgplot.pgsci(6) 
ppgplot.pgrect(t_min, twie_dn, a_min, a_max) 
ppgplot.pgrect(twib_dn, t_max,  a_min, a_max) 

ppgplot.pgsci(12)
ppgplot.pgrect(t_min, ss_dn, a_min, a_max) 
ppgplot.pgrect(sr_dn, t_max,  a_min, a_max) 

ppgplot.pgsci(1)

# Read (RA,DEC) FROM OBJECT LIST
############################################################
open(FILE, "filename") || die ("cannot open filename")
n=0

while(line=<FILE>) {
    if (line !~ /^#.*|^\n/) {
	line =~ tr/+/ /
	line =~ s/%3A/':'/eg
	line =~ s/%2B/'+'/eg
	(obj, id, ra, dec, note) = split(/\s+/, line)
	note =~ s/%[0-9A-F][0-9A-F]/'|'/eg

	ra_deg=ra_to_deg(ra)
        alpha=ra_deg*pi/180
	delta_deg=dec_to_deg(dec)
        delta=delta_deg*pi/180

	ra_deg=ra_to_deg("ra")
	alpha=ra_deg*pi/180
	delta_deg=dec_to_deg(dec)
	delta=delta_deg*pi/180

# Equinox conversion (e.g., 1950->2000)
# NOT YET

# Derive SideRealTime
#   then
# (RA,DEC)->(Az,h)

# Caliculate Airmass and Store them to T[i], secz[i]
#####################################################################
	i=0
	for (t=t_min_curve t<t_max_curve t+=t_delta) {

	    T[i]=t - t_lmstmidn 
	    Theta=t*15*pi/180 

	    sin_h= sin(phai)*sin(delta)+cos(phai)*cos(delta)*cos(Theta-alpha)
	    cos_h_cos_A = -cos(phai)*sin(delta)+sin(phai)*cos(delta)*cos(Theta-alpha)
	    cos_h_sin_A = cos(delta)*sin(Theta-alpha)

	    secz[i]=1/sin_h

	    t_abs = abs (T[i] )
            if (t_abs < 8 and (i > 1) )  {
#            print "<P> absolute t_abs: t_abs\n"
		if ((sin_h < 0.5) and (sin_h_buf > 0.5)){
		    t_lab = t - t_lmstmidn 
		    secz_lab = 1/sin_h 
		} elsif ((sin_h > 0.5) and (sin_h_buf < 0.5 )){
		    t_lab = t - t_lmstmidn 
		    secz_lab = 1/sin_h 
		}
	    }
#           print "<P> LABEL:  t:T[i] s:secz[i] i:i\n"
#            print "<P> LABEL:  t:t s:sin_h i:i\n"
            sin_h_buf = sin_h  
	    i++
	} # for loop T  

#####################################################################
# LEGEND 
#####################################################################

        ppgplot.pgsci(1) # color black
        ppgplot.pgslw(2) # line width
#       ppgplot.pgsch(1) # character height (size)
        ppgplot.pgsch(0.75) # character height (size)

        ppgplot.pgtext(-9.5,2.4+n*0.06,"obj")
        ppgplot.pgtext(-9,2.4+n*0.06,"id")
        ppgplot.pgtext(-6,2.4+n*0.06,"ra")
        ppgplot.pgtext(-3,2.4+n*0.06,"dec")
        ppgplot.pgtext(0,2.4+n*0.06,"note")


# Airmass Plot with Label
#####################################################################
	ppgplot.pgsch(1) # character height (size)

	if (obj =~ "o") {
	    ppgplot.pgsci (2)
            ppgplot.pgscf(2)   # (1: normal, 2: roman, 3: italic, 4: script)
            ppgplot.pgsch(0.8) # character height (size)
 	    ppgplot.pgslw(1) # line width
    	    ppgplot.pgptxt(t_lab, 2.02, -90.0, 0.0,"id") 
	    ppgplot.pgslw(6) # line width
	    ppgplot.pgpt(i, *T, *secz,-1)
            ppgplot.pgscf(2)   # (1: normal, 2: roman, 3: italic, 4: script)

	}else{
	    ppgplot.pgsci (4)
            ppgplot.pgscf(2)   # (1: normal, 2: roman, 3: italic, 4: script)
	    ppgplot.pgslw(1) # line width
            ppgplot.pgscf(1)   # (1: normal, 2: roman, 3: italic, 4: script)
            ppgplot.pgsch(0.6) # character height (size)
	    ppgplot.pgptxt(t_lab, 2.02, -90.0, 0.0,"id")
	    ppgplot.pgslw(2) # line width
	    ppgplot.pgpt(i, *T, *secz,-1)
            ppgplot.pgscf(2)   # (1: normal, 2: roman, 3: italic, 4: script)
	}
        n++

    }# if line

} # while

# Draw Axes and Label
#####################################################################
ppgplot.pgsci (1)

l_min=4
l_max=8

ppgplot.pgslw(3) # line width
ppgplot.pgslw(2) # line width
ppgplot.pgsch(1.0) # character height (size)

for (l=l_min l<=12 l+=1){
    ppgplot.pgtext(-12.012+l*0.974,0.96,"l")
}

for (l=1 l<=l_max l+=1){
    ppgplot.pgtext(-0.012+l*0.974,0.96,"l")
}

ppgplot.pgslw(1) # line width
ppgplot.pgsls(1) # line sytle 1: --- 2: --- 3: .-. : 4 ... 5: -...
ppgplot.pgbox('G', 1.0, 0, '0', 0.1, 0) # coordinate settings
ppgplot.pgsls(1) # line sytle 1: --- 2: --- 3: .-. : 4 ... 5: -...

ppgplot.pgslw(1) # line width
ppgplot.pgsls(4) # line sytle 1: --- 2: --- 3: .-. : 4 ... 5: -...
ppgplot.pgbox('0', 1.0, 0, 'G', 0.1, 0) # coordinate settings
ppgplot.pgsls(1) # line sytle 1: --- 2: --- 3: .-. : 4 ... 5: -...

ppgplot.pgslw(3) # line width
ppgplot.pgsls(1) # line sytle 1: --- 2: --- 3: .-. : 4 ... 5: -...
ppgplot.pgbox('G', 8.0, 0, 'G', 0.5, 0) # coordinate settings
ppgplot.pgsls(1) # line sytle 1: --- 2: --- 3: .-. : 4 ... 5: -...

ppgplot.pgend 

#####################################################################
# PLOT END
#####################################################################

system("gv -color -media letter psfile")
# system("rm  -f psfile")



###########################################################################
###########################################################################
# SUB ROUTINES
###########################################################################
###########################################################################

sub ra_to_deg {
    my (h,m,s) = split(':', @_[0])

    return h*15 + m*15/60 + s*15/3600
}

sub dec_to_deg {
    my (d,m,s) = split(':', @_[0])

    # NOTE: d should be string ex. "8","-8","+4","-0"

    if (d =~ /-.*/) {
        return -(-d + m/60 + s/3600)
    }
    else {
        return d + m/60 + s/3600
    }
}

###########################################################################
###########################################################################
# END
###########################################################################
###########################################################################

