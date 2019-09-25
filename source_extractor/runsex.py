#! /usr/bin/env python
'''


'''
from pyraf import iraf
import os, sys
#image=sys.argv[1]
#auto=sys.argv[2]

def runsextractor(image,auto):
    s="sex "+str(image)
    flag=0
    while ~flag:
        print s
	os.system(s)
	os.system('getxycat.pl')
	iraf.display(image,1,fill='yes')
	iraf.display(image,2,fill='yes')
	iraf.display('check.fits',3,fill='yes')
	iraf.display('background.fits',4)
	iraf.tvmark(2,'testxy.cat',color=204,radii=2)
	#iraf.tvmark(4,'testxy.cat',color=204,radii=2)
	if auto > 0.:
		break
	try:
	    flag= raw_input("how does it look?  1=keep, 0=adjust default.sex \n")
	    flag=int(flag)
	except ValueError:
	    print "Sorry, didn't get that.  Let's try one more time."
	    flag= raw_input("how does it look?  1=keep, 0=adjust default.sex \n")
	    flag=int(flag)
	if flag > .1:
	    break
	junk= raw_input("Edit default.sex.  Hit any key when ready \n")	

        #runsextractor(image,auto)

def writedefaultsex(ZP,gain,fwhm):
  outfile=open('default.sex','w')
  outfile.write("# Default configuration file for SExtractor V1.2b14 - > 2.0\n")
  outfile.write("# EB 26/10/97 \n")
  outfile.write("# (*) indicates parameters which can be omitted from this config file. \n")
  
  outfile.write("#-------------------------------- Catalog ------------------------------------ \n")
  outfile.write(" \n")
  outfile.write("CATALOG_NAME	testcat.fits	# name of the output catalog \n")
  outfile.write("CATALOG_TYPE	FITS_1.0	# \"ASCII_HEAD\",\"ASCII\",\"FITS_1.0\" or \"FITS_LDAC\" \n")
  outfile.write("PARAMETERS_NAME	default.param	# name of the file containing catalog contents \n")
  
  outfile.write("#------------------------------- Extraction ---------------------------------- \n")
  outfile.write("DETECT_TYPE	CCD             # \"CCD\" or \"PHOTO\" (*) \n")
  outfile.write("FLAG_IMAGE	flag.fits	# filename for an input FLAG-image \n")
  outfile.write("DETECT_MINAREA	20.		# minimum number of pixels above threshold \n")
  outfile.write("DETECT_THRESH	1.		# <sigmas> or <threshold>,<ZP> in mag.arcsec-2 \n")
  outfile.write("ANALYSIS_THRESH	1.		# <sigmas> or <threshold>,<ZP> in mag.arcsec-2 \n")
  outfile.write("FILTER		Y		# apply filter for detection (\"Y\" or \"N\")? \n")
  outfile.write("FILTER_NAME	tophat_3.0_3x3.conv	# name of the file containing the filter \n")
  outfile.write("DEBLEND_NTHRESH	32		# Number of deblending sub-thresholds \n")
  outfile.write("DEBLEND_MINCONT	0.005		# Minimum contrast parameter for deblending \n")
  outfile.write("CLEAN		Y		# Clean spurious detections? (Y or N)? \n")
  outfile.write("CLEAN_PARAM	1.0		# Cleaning efficiency \n")
  outfile.write("MASK_TYPE	CORRECT		# type of detection MASKing: can be one of \n")
  outfile.write("				# \"NONE\", \"BLANK\" or \"CORRECT\" \n")
  
  outfile.write("#------------------------------ Photometry ----------------------------------- \n")
  outfile.write("PHOT_APERTURES	2.2,3.1,4.4,6.3,8.9,12.6,17.8,25.1,35.5,50.3,71.1,100.5 \n")
  outfile.write("		# MAG_APER aperture diameter(s) in pixels \n")
  outfile.write("PHOT_AUTOPARAMS	2.5,3.5		# MAG_AUTO parameters: <Kron_fact>,<min_radius> \n")
  outfile.write("SATUR_LEVEL	55.0	# level (in ADUs) at which arises saturation \n")
  outfile.write("MAG_ZEROPOINT	%5.2f		# magnitude zero-point \n"%(ZP))
  outfile.write("MAG_GAMMA	3.0		# gamma of emulsion (for photographic scans) \n")
  outfile.write("GAIN		%5.2f		# detector gain in e-/ADU. \n"%(gain))
  outfile.write("PIXEL_SCALE	0.3		# size of pixel in arcsec (0=use FITS WCS info). \n")
  
  outfile.write("#------------------------- Star/Galaxy Separation ---------------------------- \n")
  outfile.write("SEEING_FWHM	%5.2f		# stellar FWHM in arcsec \n"%(fwhm))
  outfile.write("STARNNW_NAME	default.nnw	# Neural-Network_Weight table filename \n")
  
  outfile.write("#------------------------------ Background ----------------------------------- \n")
  outfile.write("BACK_SIZE	100		# Background mesh: <size> or <width>,<height> \n")
  outfile.write("BACK_FILTERSIZE	5		# Background filter: <size> or <width>,<height> \n")
  outfile.write("BACKPHOTO_TYPE	GLOBAL		# can be \"GLOBAL\" or \"LOCAL\" (*) \n")
  outfile.write("BACKPHOTO_THICK	24		# thickness of the background LOCAL annulus (*) \n")
  
  outfile.write("#------------------------------ Check Image ---------------------------------- \n")
  outfile.write("CHECKIMAGE_TYPE	APERTURES	# can be one of \"NONE\", \"BACKGROUND\", \n")
  outfile.write("				# \"MINIBACKGROUND\", \"-BACKGROUND\", \"OBJECTS\", \n")
  outfile.write("				# \"-OBJECTS\", \"SEGMENTATION\", \"APERTURES\", \n")
  outfile.write("				# or \"FILTERED\" (*) \n")
  outfile.write("CHECKIMAGE_NAME	check.fits	# Filename for the check-image (*) \n")
  
  outfile.write("#--------------------- Memory (change with caution!) ------------------------- \n")
  outfile.write("MEMORY_OBJSTACK	2000		# number of objects in stack \n")
  outfile.write("MEMORY_PIXSTACK	100000		# number of pixels in stack \n")
  outfile.write("MEMORY_BUFSIZE	512		# number of lines in buffer \n")
  
  outfile.write("#----------------------------- Miscellaneous --------------------------------- \n")
  outfile.write("VERBOSE_TYPE	NORMAL		# can be \"QUIET\", \"NORMAL\" or \"FULL\" (*) \n")
  
  outfile.write("#------------------------------- New Stuff ----------------------------------- \n")
  outfile.write("# Surprise!! \n")
  
  outfile.close()
