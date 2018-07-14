#!/home/abrahams/epd_free-7.3-2-rh5-x86/bin/python

import sys, os
import gt_apps as gaps
import numpy as np
from astLib import astCoords
import pyfits
from astLib import astWCS
import healpy as hp
import make3FGLxml
from scipy import ndimage

import pyLikelihood
from BinnedAnalysis import *
from SelectLoad import *
# ########################
# from ModelMake import *
from ModelMake_v2 import *
from LikeRun import *

from Type_count import Cloud


def CoordRead():
	f = open('CldNum.txt','r')
	dat = f.readlines()
	f.close()
	line = int(sys.argv[1])
	lon = float( dat[line].split(',')[0].split('(')[1] )
	lat = float( dat[line].split(',')[1].split(')')[0] )
	print "Great! Running cloud %s, at (l,b)=(%s,%s)"%(line-1,lon,lat)
	return (lon,lat)


# (250,30) is a clean, control sample.
###(glon,glat) = CoordFind()

emin = [250]
emax = [10000]
ebns = [25]
pre = '/home/abrahams/HICO_survey/SourceSearch/'

## t1 = 320377417
## t2 = 406585417

cnt = 0

(glon,glat) = CoordRead()
print "Making files for cloud at (%s,%s)"%(glon,glat)

f = open('CldNum.txt','r')
dat = f.readlines()
f.close()
line = int(sys.argv[1])
a = dat[line].split('\t')
mdd = []
#if ( len(a) <= 2 ):
#	sys.exit("Cloud Finished")
#else:
#	for m in range( len(a) ):
#		if ( 'EX' in a[m] ):
#			mdd.append( 'ptsrc' )
#		elif ( 'AGN' in a[m] ):
#			mdd.append( 'all_plus' )
#		elif ( 'CLD' in a[m] ):
#			mdd.append( 'nocld' )
#		elif ( 'CO' in a[m] ):
#			mdd.append( 'noco' )
#		elif ( 'DG' in a[m] ):
#			mdd.append( 'nodg' )

cld = Cloud(int(sys.argv[1]))

try:
	cld.getLike()
	for key in cld.like.keys():
		if np.isnan( cld.like[key] ):
			mdd.append( key )
except:
	mdd = ['all','nocodg','noco','nodg','ptsrc','all_plus']

ndir = pre +(('l%sb%s')%(int(glon),int(glat)))
print "\n\n\n\n Creating (l,b) = (%s,%s)"%(int(glon),int(glat))
try:
	os.mkdir(ndir)
except:
	print "Directory already exists ... continuing"
	pass

os.chdir(ndir)
for j in np.arange( len(emin) ):      # old tmax=392574813
	ins = (['glon=%s'%glon, 'glat=%s'%glat, 'rad=10', 'zcut=False', 
		'emin=%s'%emin[j], 'emax=%s'%emax[j], 'tmin=239557417',
                'tmax=455067824', 'nxpix=140', 'nypix=140', 'binsz=0.1',
                'enumbins=%s'%ebns[j]])
	Inputs = ReadParams(ins)
	velcut = 20	#km/s
	try:
		os.mkdir( ndir + '/ebin%s_%s'%(emin[j],emax[j]) )
	except:
		pass

	os.chdir( ndir + '/ebin%s_%s'%(emin[j],emax[j]) )
	Tspin = 125.
	Flname = FileNames(Inputs,pre,velcut,Tspin)
	cat = pyfits.open(Flname['cat'])
	Flname['cat_var'] = cat
	Flname['velct'] = velcut
	Flname['Tspin'] = Tspin
	Flname['time'] = False
	Select(Inputs,Flname)
	mktime(Inputs,Flname)

	# Cycle through the gtbin options, used for diagnostic
	# purposes (CMAP used in analysis)
	alg = ['CMAP','LC','CCUBE']
	for mode in alg:
		Flname['binalg'] = mode
		binning(Inputs,Flname)

	# Write "Config.txt" file necessary to reproduce analysis
	if ('Config.txt' not in os.listdir('.')):
		f = open('Config.txt','w')
		for elem in ins:
			f.write(elem + '\n')
	
		f.close()
   
	Flname['opacity'] = 'thick'
#	if ( (Flname['ltcube'] not in os.listdir(pre+'Analysis/')) 
#             and (cnt == 0) ):
#		expcube(Inputs,Flname,False)

#       Why check for dust? don't know
	if ( 'Dust.fits' not in os.listdir(ndir)):	
		HI_map(Inputs,Flname,1)			# Nearby
		HI_map(Inputs,Flname,2)			# far
		CO_map(Inputs,Flname)
		Avmap_make(Inputs,Flname)
		Dust_resid(Inputs,Flname)

	# Make a model -- with CO and dark gas and everything
	#              -- with no CO and DG
	#              -- replace CO/DG with a point source

	if ('bexpmap.fits' in os.listdir('.')):
		print "happy, already have binned exposure map."
	else:
		Bexpmap_mk(Inputs,Flname)

	if ('MODEL.xml' not in os.listdir('.')):
		mod = make3FGLxml.srcList(Flname['cat'],
                                          Flname['mkt_fl'],'MODEL.xml')
		mod.makeModel(Flname['gal'],'gll_iem_v06',
                              Flname['iso'],'iso_P8R2_SOURCE_V6_v06',
                              ExtraRad=3)

	os.chdir( ndir )
	Inputs = None
	Flname = None
	mod = None
	ins = None
	cat = None
	ndir = None
	del Inputs, Flname, mod, ins, cat, ndir

    

print "\n\n\n\nFinished creating files, starting analysis\n\n\n\n"

ndir = pre +(('l%sb%s')%(int(glon),int(glat)))
print "\n\n\n\n Analyzing (l,b) = (%s,%s)"%(int(glon),int(glat))
try:
	os.mkdir(ndir)
except:
	print "Directory already exists ... continuing"
	pass

os.chdir(ndir)
for j in np.arange( len(emin) ):      # old tmax=392574813
	ins = (['glon=%s'%glon, 'glat=%s'%glat, 'rad=10', 'zcut=False', 
		'emin=%s'%emin[j], 'emax=%s'%emax[j], 'tmin=239557417', 'tmax=455067824', 
		'nxpix=140', 'nypix=140', 'binsz=0.1', 'enumbins=%s'%ebns[j]])
	Inputs = ReadParams(ins)
	velcut = 20	#km/s
	try:
		os.mkdir( ndir + '/ebin%s_%s'%(emin[j],emax[j]) )
	except:
		pass

	os.chdir( ndir + '/ebin%s_%s'%(emin[j],emax[j]) )
	Tspin = 125.
	Flname = FileNames(Inputs,pre,velcut,Tspin)
	cat = pyfits.open(Flname['cat'])
	Flname['cat_var'] = cat
	Flname['velct'] = velcut
	Flname['Tspin'] = Tspin
	Flname['time'] = False
	Flname['opacity'] = 'thick'
	#####################################################################
	#####################################################################
	# Explanation of 'mode':					    #
	#####################################################################
	# orig		-- gal_2yearp7v6_v0: to fix point sources	    #
	# all		-- HI + CO + DG					    #
	# NoCODG	-- HI only					    #
	# nodg		-- HI + CO					    #
	# noco		-- HI + DG					    #
	# ptsrc		-- HI + point source at W_CO peak (no CO nor DG)    #
	# all_plus	-- HI + CO + DG + point source at W_CO peak	    #
	# power_law	-- all diffuse sources use single powerlaw spectrum #
	#####################################################################
	#####################################################################
#	Flname['mode'] = 'orig'
#	CatSrc_Fix(Inputs,Flname)
#	mode = ['all','nocld','nocodg','nodg','noco','ptsrc']#,'all_plus']#,'power_law']
#	mode = ['all','nocld','nocodg','nodg','noco']

	for Flname['mode'] in mdd:
		print "\nWorking on mode = %s"%Flname['mode']
		Mdl_Make(Inputs,Flname,(ndir + '/ebin%s_%s'%(emin[j],emax[j]) + '/'))
		Srcmap_mk(Inputs,Flname)
		try:
			BinLike(Inputs, Flname)
		except:
			print "BinLike failed somewhere ... check ...\n\n"
			m = open('LogLike.dat','a')
			m.write("BinLike failed somewhere during '%s'\n"%Flname['mode'])
			m.close()
			pass


	os.chdir( ndir )
	print "Done cloud (%s,%s)\n\n\n\n\n\n"%(glon,glat)
	Inputs = None
	Flname = None
	mod = None
	mode = None
	ins = None
	cat = None
	ndir = None
	del Inputs, Flname, mod, mode, ins, cat, ndir


print "Done! Tspin = %s"%Tspin
