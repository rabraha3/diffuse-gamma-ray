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



def CoordFind():
	f = open('CloudProperties_2sig.txt','r')
	dat = f.readlines()
	f.close()

	i = 0
	area   = np.array([])
	glon25 = np.array([])
	glat25 = np.array([])

	for i in np.arange( len(dat) ):
	#with open('CloudProperties_2sig.txt') as flObj:
	#    for line in flObj:
		if ( '(l,b)' in dat[i] ):
			l = float(dat[i].split(':(')[1].split(', ')[0])
			b = float( dat[i].split(':(')[1].split(', ')[1].split(')')[0])
			if ( (np.absolute(b) > 25) and (np.absolute(b) < 75) ):
				glon25 = np.append(glon25,l)
				glat25 = np.append(glat25,b)
				area   = np.append( area, float(dat[i+1].split('Area = ')[1]) )

	    

	return (glon25,glat25)



# (250,30) is a clean, control sample.
(glon,glat) = CoordFind()
# Done clouds 1,2

#glon = [250]#[ 94,  96, 104, 108, 138, 198, 221, 236, 254]
#glat = [ 30]#[-36, -51, -39, -53, -52,  32,  35,  39,  63]
#glon = [158]
#glat = [-33]
emin = [250]
emax = [10000]
ebns = [25]

#emin = [250,  350, 490, 685, 950, 1340, 1870, 2615, 3655, 5115, 7150]#[ 250]
#emax = [350,  490, 685, 950,1340, 1870, 2615, 3655, 5115, 7150, 10000]#[10000]
#ebns = [  1,    1,   1,   1,   1,    1,    1,    1,    1,    1,    1]#[ 25]
#emin = [375, 562,  845, 1265, 1900, 2850, 4275,  6410]#[ 250] 250
#emax = [562, 845, 1265, 1900, 2850, 4275, 6410, 10000]#[10000]375
#ebns = [  1,   1,    1,    1,    1,    1,    1,     1]#[ 25]    1

pre = '/home/abrahams/HICO_survey/SourceSearch/'

## t1 = 320377417
## t2 = 406585417
####
cnt = 0
####
# Inputs = ReadFile('Config.txt')
for i in np.arange( len(glon) ):
	# i = 2 incomplete, i = 3 error
	# 1 = 10 is the LMC
	# i = 14,15 done, 16. 17 NOT analyzed
	# Doing i = 17, 18, 19, 20, 21, 22 ... though 23 - 30 created.
#	if ( (i<17) or (i>22) ):
#	if ( (i<19) or (i>25) ): i = 19 sucks
	if ( (i<13) or (i>18) ):
		print "Skipping cloud at (%s,%s)"%(glon[i],glat[i])
		continue


	ndir = pre +(('l%sb%s')%(int(glon[i]),int(glat[i])))
	print "\n\n\n\n Analyzing (l,b) = (%s,%s)"%(int(glon[i]),int(glat[i]))
	try:
		os.mkdir(ndir)
	except:
		print "Directory already exists ... continuing"
		pass

	os.chdir(ndir)
	for j in np.arange( len(emin) ):      # old tmax=392574813
		ins = (['glon=%s'%glon[i], 'glat=%s'%glat[i], 'rad=10', 'zcut=False', 
			'emin=%s'%emin[j], 'emax=%s'%emax[j], 'tmin=239557417', 'tmax=455067824', 
			'nxpix=140', 'nypix=140', 'binsz=0.1', 'enumbins=%s'%ebns[j]])
#		ins = (['glon=%s'%glon[i], 'glat=%s'%glat[i], 'rad=15', 'zcut=True', 
#			'emin=%s'%emin[j], 'emax=%s'%emax[j], 'tmin=239557417', 'tmax=455067824', 
#			'nxpix=200', 'nypix=200', 'binsz=0.1', 'enumbins=%s'%ebns[j]])
#		ins = (['glon=%s'%glon[i], 'glat=%s'%glat[i], 'rad=7.5', 'zcut=True', 
#			'emin=%s'%emin[j], 'emax=%s'%emax[j], 'tmin=239557417', 'tmax=455067824',
#			'nxpix=100', 'nypix=100', 'binsz=0.1', 'enumbins=%s'%ebns[j]])
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
#		if ( (Flname['ltcube'] not in os.listdir(pre+'Analysis/')) and (cnt == 0) ):
#			expcube(Inputs,Flname,False)

		if ( 'Dust.fits' not in os.listdir(ndir)):	# Why check for dust? don't know
			HI_map(Inputs,Flname,1)			# Nearby
			HI_map(Inputs,Flname,2)			# far
			CO_map(Inputs,Flname)
			Avmap_make(Inputs,Flname)
			Dust_resid(Inputs,Flname)

		# Make a model -- with CO and dark gas and everything
		#              -- with no CO and DG
		#              -- replace CO/DG with a point source

		if ('bexpmap.fits' in os.listdir('.')):
##			os.remove('bexpmap.fits')
##			Bexpmap_mk(Inputs,Flname)
			print "happy, already have binned exposure map."
		else:
			Bexpmap_mk(Inputs,Flname)

		if ('MODEL.xml' not in os.listdir('.')):
			mod = make3FGLxml.srcList(Flname['cat'],Flname['mkt_fl'],'MODEL.xml')
			mod.makeModel(Flname['gal'],'gll_iem_v06',Flname['iso'],'iso_P8R2_SOURCE_V6_v06',ExtraRad=3)

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
#		Flname['mode'] = 'orig'
#		CatSrc_Fix(Inputs,Flname)
#		mode = ['all','nocld','nocodg','nodg','noco','ptsrc']#,'all_plus']#,'power_law']
		mode = ['all','nocld','nocodg','nodg','noco']
		for Flname['mode'] in mode:
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
		Inputs = None
		Flname = None
		mod = None
		mode = None
		ins = None
		cat = None
		ndir = None
		del Inputs, Flname, mod, mode, ins, cat, ndir

    

print "Done! Tspin = %s"%Tspin
