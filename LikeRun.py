#!/home/abrahams/epd_free-7.3-2-rh5-x86/bin/python

import sys, os
import gt_apps as gaps
import numpy as np
from astLib import astCoords
from astLib import astWCS
import pyLikelihood
from BinnedAnalysis import *
import pyfits
from scipy import ndimage
from UpperLimits import UpperLimits
from scipy import integrate

#########################################
##
##	Feb. 18, 2016:
##	From WeakSrc, cannot delete diffuse sources
##	anymore.
##
#########################################
##
##	Feb 24, 2016
##	Trying: remove upper limits calculations, and
##	writing Xml(#1) BEFORE writing results in
##	"Properties.txt"
##
#########################################
#
# Likelihood module:
# This module will define all the likelihood analysis routines
# deemed necessary including binned analysis, upper limit calculations,
# and probably others such as freezing weak sources (or removing them
# frmo the model file...). It is a work in progress and utilizes the 
# standard python Fermi tools, not the quickAnalysis routines.
## ##############################################
## ##############################################
def Bexpmap_mk(Input,fl):
	""" Create the binned exposure map. 

	Input: input data, ROI position, energy range, etc.
	fl   : list of file names from SelectLoad.FileNames """
    # need gaps.gtexpcube2 and it's parameters
	if (fl['bexpmp'] in os.listdir('.')):
		return "Binned exposure map made. Continue"

	(ra,de) = astCoords.convertCoords("GALACTIC","J2000",
		Input['glon'][1],Input['glat'][1],2000)
	if ( fl['time'] == False):
		gaps.gtexpcube2['infile'] = fl['ltcube']
	elif (fl['time'] == True):
		gaps.gtexpcube2['infile'] = fl['ltcube_t']
	else:
		print "Error in binned exposure map input."

	gaps.gtexpcube2['cmap'] = 'none'
	gaps.gtexpcube2['outfile'] = fl['bexpmp']
	gaps.gtexpcube2['irfs'] = fl['IRFS']
	gaps.gtexpcube2['nxpix'] = int(500)
	gaps.gtexpcube2['nypix'] = int(500)
	gaps.gtexpcube2['binsz'] = Input['binsz'][1]
	gaps.gtexpcube2['ebinalg'] = 'LOG'
	gaps.gtexpcube2['emin'] = Input['emin'][1]
	gaps.gtexpcube2['emax'] = Input['emax'][1]
	gaps.gtexpcube2['enumbins'] = int(Input['enumbins'][1])
	if (Input['glat'][1] > 70):
		gaps.gtexpcube2['coordsys'] = 'CEL'
		gaps.gtexpcube2['xref'] = ra
		gaps.gtexpcube2['yref'] = de
		gaps.gtexpcube2['proj'] = 'AIT'
		gaps.gtexpcube2.run()
		return 0
	else:
		gaps.gtexpcube2['coordsys'] = 'GAL'
		gaps.gtexpcube2['xref'] = Input['glon'][1]
		gaps.gtexpcube2['yref'] = Input['glat'][1]
		gaps.gtexpcube2['proj'] = 'CAR'
		gaps.gtexpcube2.run()
		return 0

## ##############################################
def Srcmap_mk(Input,fl):
	""" Make the source map for the given model XML file: convolves PSF and
	instrument responses and exposure maps together with the model 

	Input: input data, ROI position, energy range, etc.
	fl   : list of file names from SelectLoad.FileNames
	key  : all/NoCODG/Ptsrc -- CO/dark gas templates present or not, and
	       a point source replacing CO/DG templates """
	if ( fl['srcmp_%s'%fl['mode'].lower()] in os.listdir('.')):
		return "source map made for mode %s"%fl['mode']
	elif ( (fl['mode'] == 'noco') or (fl['mode'] == 'nodg') or (fl['mode'] == 'nocodg') or (fl['mode'] == 'ptsrc') or (fl['mode'] == 'all_plus') ):
		os.system('cp %s %s'%(fl['srcmp_all'],fl['srcmp_%s'%fl['mode'].lower()]))
		return "Copied original source map into directory for %s"%fl['mode']

	gaps.srcMaps['scfile'] = fl['SC_fl']
	if (fl['time'] == False):
		gaps.srcMaps['expcube']= fl['ltcube']
	elif (fl['time'] == True):
		gaps.srcMaps['expcube']= fl['ltcube_t']
	else:
		print "Error in source map input."

	gaps.srcMaps['cmap'] = fl['CCUBE']
#	gaps.srcMaps['ptsrc'] = 'no'

	gaps.srcMaps['evtype']= 3
	gaps.srcMaps['bexpmap'] = fl['bexpmp']
	gaps.srcMaps['irfs'] = fl['IRFS']
	gaps.srcMaps['srcmdl'] = fl['MODEL_%s'%fl['mode'].lower()]
	gaps.srcMaps['outfile'] = fl['srcmp_%s'%fl['mode'].lower()]

	gaps.srcMaps.run()
	return "Finished source map"

## ##############################################
def angsep(ra1,de1,ra2,de2):
	""" Find the angular separation between two (RA,dec) positions.
	Returns separation in degrees """
	d2r = np.pi/180.
	ra1 *= d2r
	ra2 *= d2r
	de1 *= d2r
	de2 *= d2r
	diffCosine=np.cos(de1)*np.cos(de2)*np.cos(ra1-ra2)+np.sin(de1)*np.sin(de2)
	dC=np.round(diffCosine,10)
	return np.arccos(dC)/d2r

## ##############################################
def FndSrcProp(Input,fl,name,cat):
	""" Find the radius from center of ROI and TS and returns
	the position in the list of sources from the model 

	Input: input data, ROI position, energy range, etc.
	fl   : list of file names from SelectLoad.FileNames
	name : name of source we want the properties of
	cat  : Fermi 2-year catalog file name """

	ind= np.where( cat[1].data['Source_Name'] == name )[0][0]

	TS = (cat[1].data['Sqrt_TS100_300'][ind]**2 + cat[1].data['Sqrt_TS300_1000'][ind]**2
		+ cat[1].data['Sqrt_TS1000_3000'][ind]**2 + cat[1].data['Sqrt_TS3000_10000'][ind]**2)
	rad= ( angsep(Input['glon'][1],Input['glat'][1],
			cat[1].data['GLON'][ind],cat[1].data['GLAT'][ind]) )

	return (ind,rad,TS)


## ##############################################
def FreezeSrc(Input,fl,like,cat,src):
	""" Freeze parameters of sources > 8.5 degrees from middle as well as
	weak sources (TS < 50). Remove sources > 15 degrees from middle 
	and sources with TS < 25. Freezes shape if > 5 degrees from middle
	Look for free parameters. If they exist check distance and TS.
	Freeze or remove as necessary 

	Input: input data, ROI position, energy range, etc.
	fl   : list of file names from SelectLoad.FileNames
	like : likelihood object for analysis
	cat  : Fermi 2-year catalog file name
	src  : name of source to be frozen"""
	wrt = open('Freezing.dat','a')
	try:
#	len( like.freePars(src)[0] )
		(ind, center_dist, TS) = FndSrcProp(Input,fl,src,cat)
		if ( TS > 75 ):		## Strong sources are free
			like = ThawParams(Input,fl,like,cat,src,ind)
			wrt.write('%s kept, TS = %s\n'%(src,TS))
			pass
		elif ( (center_dist > 4.5) ):
			like = FreezeSrc2(Input,fl,like,cat,src,ind)
			wrt.write('%s totally frozen......\n'%src)
		elif ( (center_dist > 2) and (TS > 250 ) ):
			like = FreezeShape(Input,fl,like,cat,src,ind)
			wrt.write('%s shape frozen......\n'%src)
		elif ( (center_dist > 2) and (TS < 250 ) ):
			like = FreezeSrc2(Input,fl,like,cat,src,ind)
			wrt.write('%s totally frozen......\n'%src)
		elif ( (center_dist < 1) and (TS < 500 ) ): 
			if ( (src[-1] != 'c') ):
				like.deleteSource( src )
				wrt.write('\t %s killed: too close,' +
					'CO dist = %s, TS = %s.. confused'+
					' source\n'%(src,CO_peak_dist,TS))
			elif ( (fl['mode'].lower() != 'all_plus') and 
					(fl['mode'].lower() != 'ptsrc') ):
				like.deleteSource( src )
				wrt.write('\t %s killed: too close, CO dist ' +
					'= %s, TS = %s......\n'%(src,CO_peak_dist,TS))
			elif ( (fl['mode'].lower() != 'all_plus') and (fl['mode'].lower() != 'ptsrc') ):
				wrt.write('\t %s kept: confused source at dist = %s\n'%(src,TS))
		else:
			wrt.write('\t %s not touched.\n'%src)
	except:
		wrt.write('%s error\n'%src)
		pass

	wrt.close()
	return like


## ##############################################
def FreezeShape(Input,fl,like,cat,src,ind):
	""" Freeze spectral shape.

	Input: input data, ROI position, energy range, etc.
	fl   : list of file names from SelectLoad.FileNames
	like : likelihood object for analysis
	cat  : Fermi 2-year catalog file name
	src  : name of source to be frozen"""
	wrt = open('Freezing.dat','a')
	try:
		if (cat[1].data['SpectrumType'][ind] == 'PowerLaw'):
			like.freeze(like.par_index(src,'Index'))
		elif (cat[1].data['SpectrumType'][ind] == 'PowerLaw2'):
			like.freeze(like.par_index(src,'Index'))
		elif (cat[1].data['SpectrumType'][ind] == 'LogParabola'):
			like.freeze(like.par_index(src,'alpha'))
			like.freeze(like.par_index(src,'beta'))
		elif (cat[1].data['SpectrumType'][ind] == 'BrokenPowerLaw'):
			like.freeze(like.par_index(src,'Index1'))
			like.freeze(like.par_index(src,'BreakValue'))
			like.freeze(like.par_index(src,'Index2'))
		elif (cat[1].data['SpectrumType'][ind] == 'BrokenPowerLaw2'):
			like.freeze(like.par_index(src,'Index1'))
			like.freeze(like.par_index(src,'BreakValue'))
			like.freeze(like.par_index(src,'Index2'))
		elif (cat[1].data['SpectrumType'][ind] == 'SmoothBrokenPowerLaw'):
			like.freeze(like.par_index(src,'Index1'))
			like.freeze(like.par_index(src,'BreakValue'))
			like.freeze(like.par_index(src,'Index2'))
			like.freeze(like.par_index(src,'Beta'))
		elif (cat[1].data['SpectrumType'][ind] == 'ExpCutoff'):
			like.freeze(like.par_index(src,'Index'))
			like.freeze(like.par_index(src,'Ebreak'))
			like.freeze(like.par_index(src,'P1'))
		elif (cat[1].data['SpectrumType'][ind] == 'BPLExpCutoff'):
			like.freeze(like.par_index(src,'Index1'))
			like.freeze(like.par_index(src,'Ebreak'))
			like.freeze(like.par_index(src,'P1'))
			like.freeze(like.par_index(src,'Index2'))
			like.freeze(like.par_index(src,'Eabs'))
		elif (cat[1].data['SpectrumType'][ind] == 'PLSuperExpCutoff'):
			like.freeze(like.par_index(src,'Index1'))
			like.freeze(like.par_index(src,'Cutoff'))
			like.freeze(like.par_index(src,'Index2'))
		elif (cat[1].data['SpectrumType'][ind] == 'Gaussian'):
			like.freeze(like.par_index(src,'Mean'))
			like.freeze(like.par_index(src,'Sigma'))
	except:
		wrt.write('%s error\n'%src)
		pass

	wrt.close()
	return like

## ##############################################
def FreezeSrc2(Input,fl,like,cat,src,ind):
	""" For the 2nd round through the likelihood analysis, we want to
	freeze all sources except those we care about. 

	Input: input data, ROI position, energy range, etc.
	fl   : list of file names from SelectLoad.FileNames
	like : likelihood object for analysis
	cat  : Fermi 2-year catalog file name
	src  : name of source to be frozen"""

	params = like.model[src].funcs['Spectrum'].paramNames
	num_params = len(params)
	for par in params:
		like.freeze(like.par_index(src,par))

	return like

## ##############################################
def ThawParams(Input,fl,like,cat,src,ind):
	""" For the 2nd round through the likelihood analysis, we want to
	thaw some sources. 

	Input: input data, ROI position, energy range, etc.
	fl   : list of file names from SelectLoad.FileNames
	like : likelihood object for analysis
	cat  : Fermi 2-year catalog file name
	src  : name of source to be frozen"""
	if ('Add' in src):
		like.thaw(like.par_index(src,'Integral'))
		like.thaw(like.par_index(src,'Index'))
	elif (cat[1].data['SpectrumType'][ind] == 'PowerLaw'):
		like.thaw(like.par_index(src,'Prefactor'))
#		like.thaw(like.par_index(src,'Index'))
	elif (cat[1].data['SpectrumType'][ind] == 'PowerLaw2'):
		like.thaw(like.par_index(src,'Integral'))
		like.thaw(like.par_index(src,'Index'))
	elif (cat[1].data['SpectrumType'][ind] == 'LogParabola'):
		like.thaw(like.par_index(src,'norm'))
#		like.thaw(like.par_index(src,'alpha'))
#		like.freeze(like.par_index(src,'beta'))
	elif (cat[1].data['SpectrumType'][ind] == 'BrokenPowerLaw'):
		like.thaw(like.par_index(src,'Prefactor'))
		like.thaw(like.par_index(src,'Index1'))
		like.thaw(like.par_index(src,'BreakValue'))
		like.thaw(like.par_index(src,'Index2'))
	elif (cat[1].data['SpectrumType'][ind] == 'BrokenPowerLaw2'):
		like.thaw(like.par_index(src,'Integral'))
		like.thaw(like.par_index(src,'Index1'))
		like.thaw(like.par_index(src,'BreakValue'))
		like.thaw(like.par_index(src,'Index2'))
	elif (cat[1].data['SpectrumType'][ind] == 'SmoothBrokenPowerLaw'):
		like.thaw(like.par_index(src,'Prefactor'))
		like.thaw(like.par_index(src,'Index1'))
		like.thaw(like.par_index(src,'BreakValue'))
		like.thaw(like.par_index(src,'Index2'))
		like.thaw(like.par_index(src,'Beta'))
	elif (cat[1].data['SpectrumType'][ind] == 'ExpCutoff'):
		like.thaw(like.par_index(src,'Prefactor'))
#		like.thaw(like.par_index(src,'Index'))
#		like.thaw(like.par_index(src,'Ebreak'))
#		like.thaw(like.par_index(src,'P1'))
	elif (cat[1].data['SpectrumType'][ind] == 'BPLExpCutoff'):
		like.thaw(like.par_index(src,'Prefactor'))
		like.thaw(like.par_index(src,'Index1'))
		like.thaw(like.par_index(src,'Ebreak'))
		like.thaw(like.par_index(src,'P1'))
		like.thaw(like.par_index(src,'Index2'))
		like.thaw(like.par_index(src,'Eabs'))
	elif (cat[1].data['SpectrumType'][ind] == 'PLSuperExpCutoff'):
		like.thaw(like.par_index(src,'Prefactor'))
#		like.thaw(like.par_index(src,'Index1'))
#		like.thaw(like.par_index(src,'Cutoff'))
#		like.thaw(like.par_index(src,'Index2'))
	elif (cat[1].data['SpectrumType'][ind] == 'Gaussian'):
		like.thaw(like.par_index(src,'Prefactor'))
		like.thaw(like.par_index(src,'Mean'))
		like.thaw(like.par_index(src,'Sigma'))

	return like



## ##############################################
def DiffFreeThaw(Input,fl,like,src,act):

	params = like.model[src].funcs['Spectrum'].paramNames
	num_params = len(params)
	if (act.lower() == 'freeze'):
		for par in params:
			like.freeze(like.par_index(src,par))
	elif (act.lower() == 'thaw'):
		for par in params:
			like.thaw(like.par_index(src,par))

	return like


## ##############################################
def ThawSrc(Input,fl,like,src):
	(ind, center_dist, TS) = FndSrcProp(Input,fl,src,fl['cat_var'])
	like = ThawParams(Input,fl,like,fl['cat_var'],src,ind)
	return like

## ##############################################
def MdlMap(Input,fl,lk_iter,flxrun):
	""" Create model map, because the likelihood routine is giving
	me a seg fault """
	gaps.model_map['srcmaps'] = fl['srcmp_%s'%fl['mode'].lower()]
	if (flxrun == 0):
		gaps.model_map['srcmdl']  = fl['outmdl%s_%s'%(str(lk_iter),fl['mode'].lower())]
		gaps.model_map['outfile'] = fl['outmap%s_%s'%(str(lk_iter),fl['mode'].lower())]
	elif(flxrun == 1):
		pass

	gaps.model_map['irfs']    = fl['IRFS']
	gaps.model_map['bexpmap'] = fl['bexpmp']
	if ( fl['time'] == False ):
		gaps.model_map['expcube'] = fl['ltcube']
	else:
		gaps.model_map['expcube'] = fl['ltcube_t']

	gaps.model_map.run()

## ##############################################
## ##############################################
def PrepAll(Input,fl,like,likerun):
	""" Do the likelihood analysis, here I will do the model prep. 
	Prepare model: freeze point sources & fit diffuse
	then freeze diffuse & thaw point sources & fit
	then unfreeze diffuse and do standard NewMinuit """

	if (fl['mode'] != 'all'):
		return (0,like)
	g = open('LogLike.dat','a')
    # Here: freeze point sources
	name = like.sourceNames()
	print "\n1st name = %s\n"%name[0]
	for nm in name:
		if ('3FGL' in nm):
			(ind, center_dist, TS) = FndSrcProp(Input,fl,nm,fl['cat_var'])
			like = FreezeSrc2(Input,fl,like,fl['cat_var'],nm,ind)
			print ("Freeeeeeeze %s"%nm)
#		elif ( (nm in fl['gal']) or (nm in fl['iso']) ):
#			like.deleteSource(nm)
		elif ( (fl['mode'].lower() == 'noco') and (nm == 'CO') ):
			like.deleteSource(nm)
		elif ( (fl['mode'].lower() == 'nodg') and (nm == 'E(B-V)') ):
			like.deleteSource(nm)
		elif ('Add' in nm):
			like = FreezeSrc2(Input,fl,like,fl['cat_var'],nm,ind)

	try:
#		if (like.optimizer == 'DRMNFB'):
#			lglike = like.fit()
#		elif (like.optimizer == 'Minuit'):
#			like.writeXml('MDL.xml')
#			like = BinnedAnalysis(obs,srcModel='MDL.xml',optimizer='Minuit')
#			likeobj = pyLike.Minuit(like.logLike)
#			lglike = like.fit(covar=True,optObject=likeobj)
#			g.write('Fit quality %s\t'%like1obj.getQuality())
		like1obj = pyLike.Minuit(like.logLike)
		lglike = like.fit(covar=True,optObject=like1obj,verbosity=0)
		g.write('Fit quality %s\t'%like1obj.getQuality())
		print ('Fit quality %s\t'%like1obj.getQuality())
	except:
		print "Boop boop"
		like.optimizer='minuit'
		print "Optimizer = 'minuit'"
		for nm in name:
			if ( ((nm == 'CO') or (nm == 'E(B-V)') or (nm == 'HI_far')) 
			     and (like.Ts(nm)< 3) ):
					like.deleteSource(nm)
			elif ( (nm == 'HI') and ( fl['gas_spec'].lower() != 'spec') ):
				like.model['HI'].funcs['Spectrum'].params['Index1'].setValue(-1.8)
				like.model['HI'].funcs['Spectrum'].params['Index1'].setScale(1)
				like.model['HI'].funcs['Spectrum'].params['Integral'].setValue(2)

		likeobj = pyLike.Minuit(like.logLike)
		lglike = like.fit(covar=True,optObject=likeobj,verbosity=0)

	# Delete weak diffuse sources and freeze strong diffuse sources and thaw good point sources
	print "Delete weak diffuse sources, freeze rest, thaw 'good' point sources."
	name = like.sourceNames()
	for nm in name:
		if ( ('3FGL' not in nm) and ('Add' not in nm) and ('iso' not in nm) ):
#			if ( like.Ts(nm) < 3 ):
#				like.deleteSource(nm)
#			else:
			like = DiffFreeThaw(Input,fl,like,nm,'freeze')
			print ("Frozen %s!"%nm)
		elif ( (('IC' in nm) or ( 'galprop' in nm )) and (like.Ts(nm) < 3) ):
			like.deleteSource(nm)
		elif ( '3FGL' in nm ):
			like = ThawSrc(Input,fl,like,nm)


	print ("Frozen Sources !! :). Fit point sources.\n")
	lgng = open('Logging.txt','a')
	lgng.write("Frozen Diffuse Sources !! :). Fit point sources.\n")
	lgng.close()
	# freeze some sources: fit only important ones
	like = MdlPrep(Input,fl,like,likerun)
	try:
		if (like.optimizer.lower() == 'drmnfb'):
			print "Starting 1nd fit of point sources with DRMNFB"
			lglike = like.fit(verbosity=0)
		elif (like.optimizer.lower() == 'minuit'):
			print "Starting 1nd fit of point sources with Minuit"
			likeobj2 = pyLike.Minuit(like.logLike)
			lglike2 = like.fit(covar=True,optObject=likeobj2,verbosity=0)
			print "Fit Quality = %s\n"%likeobj2.getQuality()
	except:
		print 'poop'
		lgng = open('Logging.txt','a')
		lgng.write("Poop, can't even accomplish step 1...\n")
		lgng.close()
		return ("poop")

	lgng = open('Logging.txt','a')
	lgng.write("Finished 1st fit of ponit sources\n")
	lgng.close()
	for nm in name:
		if ('3FGL' in nm):
			try:
				aa = len(like.freePars(nm))
				if ( like.Ts(nm) < 3 ):
					print "Deleting %s"%nm
					like.deleteSource(nm)
			except:
				pass

	# Need to freeze point sources and unfreeze diffuse sources...
	name = like.sourceNames()
	print "Freezing point sources, unfreezing diffuse"
	for nm in name:
		if ( ('3FGL' not in nm) and ('Add' not in nm) ):
			like = DiffFreeThaw(Input,fl,like,nm,'thaw')
		elif ( ('3FGL' in nm) ):
			like = DiffFreeThaw(Input,fl,like,nm,'freeze')
		elif ( 'Add' in nm ):
			like.freeze(like.par_index(nm,'Index'))
		elif ( (nm == 'CO') or (nm == 'E(B-V)') or (nm == 'HI') or (nm == 'HI_far') ):
			l_bnd=like.model[nm].funcs['Spectrum'].params['Normalization'].getBound()[0]
			nm_val=like.model[nm].funcs['Spectrum'].params['Normalization'].getValue()
			if ( ((nm_val-l_bnd)/l_bnd) < 1 ):
				like.deleteSource(nm)
				print "Deleted %s because normalization at lower limit"%nm

	try:
		print "\n\n2nd diffuse source fitting\n\n"
		lgng = open('Logging.txt','a')
		lgng.write("2nd diffuse source fitting.\n")
		lgng.close()
		likeobj = pyLike.Minuit(like.logLike)
		lglike = like.fit(covar=True,optObject=likeobj,verbosity=0)
		g.write('Fit quality %s\t'%like1obj.getQuality())
		print ('Fit quality %s\t'%like1obj.getQuality())
	except:
		lgng = open('Logging.txt','a')
		lgng.write("2nd diffuse fitting in 1st round = BAD\n")
		lgng.close()
		pass
	
	try:
		if ( like.model['IC'].funcs['Spectrum'].params['Value'].getValue() < 1e-4 ):
			like.deleteSource('IC')
			print "IC deleted, too weak"
	except:
		print "IC not found"
		pass

	try:
		if ( like.model['HI_far'].funcs['Spectrum'].params['Normalization'].getValue() < 0.1 ):
			like.deleteSource('HI_far')
			print "HI_far deleted, too weak"
	except:
		print "HI_far gone already"
		pass

	lgng = open('Logging.txt','a')
	lgng.write("Done round 1, writing model.\n")
	lgng.close()
	like.writeXml('mdl1_all.xml')
	print "sources fit, finding new sources.\n"


	return (lglike2,like)
	# Find new sources here
#
#	MdlMap(Input,fl,1,0)
#	nsrc_ = NewSrc_Detect(Input,fl)
#
#	if ("none" in nsrc_.lower()):
#		print "no new sources\n"
#		os.system('rm %s'%fl['outmap1_%s'%fl['mode'].lower()])
#		print "removed model map...\n"
#		lglike2 = lglike
#	else:
#		os.system('mv mdl1_all.xml mdl1_int1.xml')	# Not needed...?
#		print "new sources found!"
#		fl['mode'] = 'int1'
	#    No need to make new sourcemap, the likelihood procedure
	#    make one for us!
#
#		obs = BinnedObs(binnedExpMap=fl['bexpmp'],
#			srcMaps=fl['srcmp_all'],expCube=fl['ltcube'],irfs=fl['IRFS'])
#		like = BinnedAnalysis(obs, srcModel='mdl1_int1.xml',optimizer='minuit')
#		try:
#			if (like.optimizer.lower() == 'drmnfb'):
#				lglike2 = like.fit()
#			elif (like.optimizer.lower() == 'minuit'):
#				likeobj2 = pyLike.Minuit(like.logLike)
#				lglike2 = like.fit(covar=True,optObject=likeobj2)
#		except:
#			return ("poop")
#
#		fl['mode'] = 'all'
#	return (lglike,like)


## ##############################################
def BinLike(Input,fl):
	""" Perform the binned likelihood analysis using the standard 
	gt_apps: the python wrapper for the ballistic commands 
	Freeze sources completely if too far away and freeze the 
	shapes of sources at intermediate distances from ROI center
	Start likelihood with DRMNFB optimizer at coarse tolerance. 
	Write resulting fit to XML file. Then freeze everything 
	except gas templates and possibly added sources and do a 
	second likelihood analysis with NewMinuit optimizer."""

	srcmdl = fl['MODEL_%s'%fl['mode'].lower()]
	srcmap = fl['srcmp_%s'%fl['mode'].lower()]
	outmdl1 = fl['outmdl1_%s'%fl['mode'].lower()]
	outmdl2 = fl['outmdl2_%s'%fl['mode'].lower()]
	likerun = 1
	print "Start 1!"

	if ( fl['time'] == False):
		obs = BinnedObs(binnedExpMap=fl['bexpmp'],
			srcMaps=srcmap,expCube=fl['ltcube'],irfs=fl['IRFS'])
	elif ( fl['time'] == True):
		obs = BinnedObs(binnedExpMap=fl['bexpmp'],
			srcMaps=srcmap,expCube=fl['ltcube_t'],irfs=fl['IRFS'])
	else:
		print "Error in binned analysis input."

	print "Observations set up, creating likelihood object"
#	like = BinnedAnalysis(obs, srcModel=srcmdl,optimizer='DRMNFB')
	like = BinnedAnalysis(obs, srcModel=srcmdl,optimizer='Minuit')
	likeobj = pyLike.Minuit(like.logLike)

	g = open('LogLike.dat','a')
	print "Setting up likelihood analysis."
	optim = 1
	try:
		if (fl['mode'].lower() == 'all'):
			print "Fitting first round."
			lgng = open('Logging.txt','a')
			lgng.write("Starting round 1 ...\n")
			lgng.close()
			(lglike,like) = PrepAll(Input,fl,like,likerun)
			print "-log(likelihood) = %s"%lglike
			print "Done first round."
		else:
			print "mode =/= 'all', fitting..."
			lgng = open('Logging.txt','a')
			lgng.write("mode =/= all, fitting 1st round.\n")
			lgng.close()
			lglike = like.fit(covar=True,optObject=likeobj,verbosity=0)
	except:
		return "poops. Can't finish fitting for some reason. Check the 'PrepAll' routine, probably.\n\n\n Done \n\n\n :("

	print "Writing likelihoods.\n"
	g.write('-log(likelihood) of %s = %s\n'%(fl['mode'],lglike))

	# Check TS's for diffuse sources
	name = like.sourceNames()
#	for nm in name:
#		if ( ('3FGL' not in nm) and ('Add' not in nm) and ('iso' not in nm) ):
#			if ( like.Ts(nm) < 3 ):
#				like.deleteSource(nm)
#			print ("Deleted %s!"%nm)
#

	print "Weak diffuse sources gone, writing results\n."
	lgng = open('Logging.txt','a')
	lgng.write("Weak diffuse sources gone, writing results.\n")
	lgng.close()
	like.writeXml(outmdl1)
	try:
		WriteResults(Input,fl,like,optim)
	except:
		print "Results won't write, do it by hand."
		lgng = open('Logging.txt','a')
		lgng.write("Results won't write, do it by hand.\n")
		lgng.close()

	try:
		like = WeakSrc(Input,fl,like) # Get rid of weak sources, including diffuse
	except:
		lgng = open('Logging.txt','a')
		lgng.write("Can't kill weak sources ... let's soldier on.\n")
		lgng.close()

	# 2nd fitting ... setting up
	like = BinnedAnalysis(obs,srcModel=outmdl1,optimizer='NewMinuit')
	like2obj = pyLike.NewMinuit(like.logLike)
	MdlMap(Input,fl,1,0)

    	# 2nd likelihood anlaysis to fine tune. 
	like.tol=1e-4
	like = MdlPrep(Input,fl,like,2)

#	lging = open('Logging.txt','a')
#	lging.write('%s'%like.model)
#	lging.close()
	lgng = open('Logging.txt','a')
	lgng.write("2nd round set up, fitting.\n")
	lgng.close()
	try:
		lglike2 = like.fit(covar=True,optObject=like2obj,verbosity=0)
		g.write('likelihood return code = %s\t'%like2obj.getRetCode())
		g.write('Final -log(likelihood) of %s = %s\n'%(fl['mode'],lglike2))
		like.writeXml(outmdl2)
		MdlMap(Input,fl,2,0)
		error = 0
	except:
		g.write('ERROR IN NEWMINUIT FIT\n\n')
		lgng = open('Logging.txt','a')
		lgng.write("Problem in step 2 (New Minuit)...\n")
		lgng.close()
		error = 1
		pass


	g.close()
	try:
		if ( fl['mode'].lower() == 'all' ):
			covs= np.array(like.covariance)
			cov = open('Covariance.txt','a')
			cov.write('%s:'%len(covs))
			for i in np.arange( len(covs) ):
				cov.write('\n')
				for j in np.arange( len(covs) ):
					cov.write('%s,'%covs[i,j])

			cov.close()
	except:
		print "Error writing covariance matrix\n\n"

	WriteResults(Input,fl,like,error)
	obs = None
	like = None
	convs = None
	cov = None
	srcmdl = None
	srcmap = None
	outmdl1 = None
	outmdl2 = None
	likeobj = None
	like2obj = None
	del obs, like, convs, cov, srcmdl, srcmap, outmdl1, outmdl2, likeobj, like2obj
	lgng = open('Logging.txt','a')
	lgng.write("Finished fitting mode = %s\n\n\n\n\n\n\n\n\n\n\n\n\n\n"%fl['mode'])
	lgng.close()
	return "Done Likelihood. Run checks."


## ##############################################
def MdlPrep(Input,fl,like,likerun):
	""" We prepare the model in "like" for fitting. This means that
	we freeze sources deemed unworthy """

	like = ConfuSearch(Input,fl,like)
	nm = like.sourceNames()
	name = nm
	wrt = open('Freezing.dat','a')
	wrt.write('\n\n')
	try:
		for src in name:
			if ( (fl['mode'].lower() == 'noco') and (src == 'CO') ):
				like.deleteSource(src)
				wrt.write('Deleting %s'%src)
			elif ( (fl['mode'].lower() == 'nodg') and (src == 'E(B-V)') ):
				like.deleteSource(src)
				wrt.write('Deleting %s'%src)
	except:
		pass

	name = like.sourceNames()
	if ( likerun == 1 ):
		wrt.write('DRMNFB or Minuit step:\n')
	else:
		wrt.write('New Minuit step:\n')

# If we want to fix all the point sources after fitting to the 
# Fermi diffuse maps, then set likerun = 2, else, set 
# likerun=1
	if ( fl['mode'].lower() == 'fx_src' ):
		for src in name:
			if ( (src == 'CO') or (src == 'HI') or (src == 'E(B-V)') or (src == 'HI_far') ):
				like.freeze(like.par_index(src,'Integral'))
				like.freeze(like.par_index(src,'Index'))
	
		return like

	if ( likerun == 1 ):
		for src in name:
			if ('iso' in src ):
				#like.freeze(like.par_index(src,'Normalization'))
				print "Keeping %s\n"%src
			elif ( ('3FGL' in src) and ((src != 'SMC') or (src != 'LMC') or (src != " ")) ):
				try:
					like = FreezeSrc(Input,fl,like,fl['cat_var'],src)
					print ("%s frozen."%src)
				except:
					print "Passing; no freeze"
					pass
	elif( likerun == 2 ):
		for src in name:
			if ( 'iso' in src ):
				# like.freeze(like.par_index(src,'Normalization'))
				print "Keeping %s\n"%src
			elif ( 'Add' in src ):
				like.freeze(like.par_index(src,'Index'))
			elif ( ('3FGL' in src) and ((src != 'SMC') or 
					(src != 'LMC') or (src != " ") or
					(src[-1] != 'c')) ):
				try:
					(ind, center_dist, TS) = FndSrcProp(Input,fl,src,fl['cat_var'])
					like = FreezeSrc2(Input,fl,like,fl['cat_var'],src,ind)
					wrt.write('%s frozen.\n'%src)
				except:
					wrt.write('%s NOT FROZEN. ERROR!!!\n'%src)
					pass
			elif ( ('3FGL' not in src) ):
				like = DiffFreeThaw(Input,fl,like,src,'thaw')

	wrt.close()
	return like

## ##############################################
def ConfuSearch(Input,fl,like):
	""" I am searching through the model in order to find confused
	sources: source names that end in the character "c". I 
	should do something similar for extended sources sometime,
	those source names ending in "e".
	"""

	check_cut = 18
	check = []
	for src in like.sourceNames():
		if ( len(src) == check_cut ):
			check.append( src )

	for name in check:
		if ( (fl['mode'].lower() != 'ptsrc') and (name[check_cut-1] == 'c') ):
			like.deleteSource( name )
			print "Deleted %s."%name
		elif ( fl['mode'].lower() == 'ptsrc' ):
			like.deleteSource( 'Add_0' )

	return like



## ##############################################
def WriteResults(Input,fl,like,error):
	""" Write the fluxes and model files (for spectral indices).
	"""


	if ( fl['mode'].lower() == 'all' ):
		mm = open('Properties.txt','a')
		like.writeCountsSpectra('CountSpec_%s.fits'%fl['mode'])
	else:
		mm = open('%s/Properties.txt'%fl['mode'],'a')
		like.writeCountsSpectra('%s/CountSpec_%s.fits'%(fl['mode'],fl['mode']))

	mdls = ['HI','HI_far','CO','E(B-V)','IC','iso_P8R2_SOURCE_V6_v06']
	for i in np.arange( len(like.sourceNames())-1 ):
		mdls.append( 'Add_%s'%int(i) )

	nm = like.sourceNames()
	for name in nm:
		if ( name[len(name)-1]  == 'c' ):
			mdls.append( name )

	# (ulim == 0) means the UpperLimits routine failed. If so, then
	# 	upper limit is set to 3 times the found flux. (error==1) means
	# 	that the uncertainties could not be calculated.
	for name in nm:
		if ( name in mdls ):
			try:
				mm.write(' %s:\n'%name)
				if ( name == 'HI' ):
					mm.write(' Tspin = %s\n'%fl['Tspin'])
				mm.write(' TS = %s\n'%like.Ts( name ) )
			except:
				mm.write('Error in writing TS for %s\n'%name)
#			try:
#				if (like.Ts( name ) < 20):
#					ulim = UpperLim(Input,like, name)
#					if ((ulim != 0) and (error != 1)):
#						lim = float(str(ulim).split(' ph')[0])
#						mm.write(' Upper Limit = %s ph/cm^2/s \n'%lim)
#					elif ( (ulim == 0) and (error != 1) ):
#						lim = 3*like.flux(name,emin=Input['emin'][1],emax=Input['emax'][1])
#						limer=3*like.fluxError(name,emin=Input['emin'][1],emax=Input['emax'][1])
#						mm.write(' Upper Limit about (3*flux) %s +/- %s\n'%(lim,limer))
#					elif ( (ulim == 0) and (error == 1) ):
#						lim = 3*like.flux(name,emin=Input['emin'][1],emax=Input['emax'][1])
#						mm.write(' Upper Limit about (3*flux) %s \n'%lim)
#			except:
#				mm.write("Error in writing upper limit for source %s\n"%name)

			try:
				if ( error != 1 ):
					mm.write('energy flux = %s +/- %s \n'%(like.energyFlux(name,
						emin=Input['emin'][1],emax=Input['emax'][1]),
						like.energyFluxError('HI',emin=Input['emin'][1],
						emax=Input['emax'][1])))
					mm.write('flux = %s +/- %s'%(like.flux(name,
						emin=Input['emin'][1],emax=Input['emax'][1]),
						like.fluxError('HI',emin=Input['emin'][1],
						emax=Input['emax'][1])))
				else:
					mm.write('energy flux = %s \n'%(like.energyFlux(name,
						emin=Input['emin'][1],emax=Input['emax'][1])))
					mm.write('flux = %s'%(like.flux(name,
						emin=Input['emin'][1],emax=Input['emax'][1])))
			except:
				mm.write('Error in writing %s\n'%name)

			try:
				mm.write('\n')
				mm.write('%s'%like.model[name])
				mm.write('\n')
			except:
				mm.write('No model for %s\n'%name)

	mm.write('\n\n\n')
	mm.close()
	return "Done writing."



## ##############################################
def NewSrc_Detect(Input,fl):
	""" Rummage through the latest model map created too see if 
	there are any missed sources (large positive residuals)
	or if the model serverely overpredicted the emission, 
	usually the result of fixing a strong source.
	"""
	cmap = pyfits.open(fl['CMAP'])
	(nx,ny) = (cmap[0].header['NAXIS1'],cmap[0].header['NAXIS2'])
	cmap[0].data = np.transpose( cmap[0].data )
	try:
		mdl  = pyfits.open(fl['outmap1_%s'%fl['mode'].lower()])
		mdl[0].data = np.transpose( mdl[0].data )
	except:
#		print "Error loading model map. Refit."
		return "None"

	resi = cmap[0].data - mdl[0].data
	res = pyfits.PrimaryHDU( resi )
	res.header = cmap[0].header
	res.writeto('resid_%s.fits'%fl['mode'])
	resi_smooth = ndimage.gaussian_filter(resi, sigma=2)
	# Choose a 4-sigma cutoff. The smoothed residual map should be
	#  distributed like a normal curve if everything is modeled 
	#  well... What about holes? Those should fit out?
	resi_smooth[ resi_smooth < 5*resi_smooth.std() ] = 0
	if ( resi_smooth.max() == 0 ):	# If nothing stands out, exit
		return "None"

	f = open('mdl1_%s.xml'%fl['mode'],'r')	
	dat = f.readlines()
	f.close()
	del dat[-1]
	cnt = 1
	WCS = astWCS.WCS(fl['CMAP'])
	while (resi_smooth.std() > 0.0):
		maxx = np.where( resi_smooth == resi_smooth.max() )[0][0]
		maxy = np.where( resi_smooth == resi_smooth.max() )[1][0]
		(glon,glat) = WCS.pix2wcs( maxx, maxy )
		dat = Diffuse_ptsrc(Input,fl,dat,glon,glat,cnt)
		for ind in np.arange( nx*ny ):
			if ( np.linalg.norm( np.array( (maxx,maxy) ) - np.array( (ind%nx,ind/nx) ) ) < 10 ):
				resi_smooth[ ind%nx,ind/nx ] = 0

		cnt += 1
	

	dat.append('</source_library>')
#	os.system('mv mdl1_int1.xml MODEL_int1.xml')
	g = open('mdl1_int1.xml','w')
	g.writelines(dat)
	g.close()

	return "Done"

## ##############################################
def CatSrc_Fix(Input,fl):
	""" Here we perform a basic DRMNFB fit with the Fermi diffuse
	galactic emission: 'gal_2yearp7v6_v0'. We then write the
	model file and then *fix* all of the point sources while
	fitting the gas templates.
	"""
	os.mkdir('%s'%fl['mode'])
	Srcmap_mk(Input,fl)
	obs = BinnedObs(binnedExpMap=fl['bexpmp'],
			srcMaps=fl['srcmp_orig'], 
			expCube=fl['ltcube'],irfs=fl['IRFS'] )
	like = BinnedAnalysis(obs,srcModel=fl['MODEL_%s'%fl['mode'].lower()],
	  	optimizer='DRMNFB')
	like.setEnergyRange(emin=Input['emin'][1],emax=Input['emax'][1])
	g = open('LogLike.dat','a')
	g.write('\n\n')
	for src in like.sourceNames():
		if ('3FGL' in src):
			try:
				(ind,rad,TS) = FndSrcProp(Input,fl,src,fl['cat_var'])
				if (rad < Input['rad'][1]+2):
					like.deleteSource(src)
				else:
					FreezeShape(Input,fl,like,cat,src,ind)
			except:
				pass


	try:
		lglike = like.fit()
	except:
		g.write('Error in DRMNFB, going to Minuit.\n')
		like = BinnedAnalysis(obs, srcModel=fl['MODEL_orig'],optimizer='Minuit')
		like.setEnergyRange(emin=Input['emin'][1],emax=Input['emax'][1])
		like1obj = pyLike.Minuit(like.logLike)
		lglike = like.fit(covar=True,optObject=like1obj)

	g.write('-log(likelihood) of %s = %s\n'%(fl['mode'],lglike))
	g.write('\n\n')
	# Freeze some things before writing the model to a file
	likerun = 1
	like = WeakSrc(Input,fl,like)
	like.writeXml('MODEL_.xml')
	MdlMap(Input,fl,1,0)
	NewSrc_Detect(Input,fl)
	return "Finished fitting with 3FGL GDE map."


## ##############################################
def WeakSrc(Input,fl,like):
	""" After fitting the model with the gal_2yearp7v6_v0 diffuse
	file, scan through the point sources and remove those with
	low TS values (less than 50) """
#	nm = like.sourceNames()
	frz = open('TS.txt','a')
	if (fl['mode'].lower() == 'all'):
		mm = open('Properties.txt','a')
	else:
		mm = open('%s/Properties.txt'%fl['mode'],'a')

	frz.write("Number of sources = %s\n"%len(like.sourceNames()))
	for name in like.sourceNames():
		TS_ = like.Ts( name )
#		if (TS_ < 0):
#			if ((name == 'CO') or (name == 'HI') or (name=='E(B-V)') or (name == 'HI_far') or (name == 'IC')):
#				like.deleteSource( name )
#				frz.write('During %s fitting, %s killed with TS of %s\n'%(fl['mode'],name,TS_))
#				mm.write('During %s fitting, %s killed with TS of %s\n\n'%(fl['mode'],name,TS_))
#		elif (TS_ < 50):
		if (TS_ < 50):
			if ((name != 'CO') and (name != 'HI') and (name != 'E(B-V)') and (name != 'HI_far') and (name != 'IC') ):
				like.deleteSource( name )
				frz.write('During %s fitting, %s killed with TS of %s\n'%(fl['mode'],name,TS_))
		else:
			frz.write('During %s fitting, %s has TS of %s\n'%(fl['mode'],name,TS_))
	

	frz.close()
	mm.close()
	return like



## ##############################################
def Diffuse_ptsrc(Inputs,fl,dat,glon,glat,cnt):
	""" Modify model file to include additional point source 

	Inputs: standard analysis inputs
	fl    : file names for analysis
	dat   : model xml string list 
	glon  : galactic longitude of source
	glat  : galactic latitude of source"""

	(RA,dec) = astCoords.convertCoords("GALACTIC","J2000",glon,glat,2000)
	dat.append('<source name="Add_%s" type="PointSource">\n'%cnt)
	dat.append('\t<spectrum type="PowerLaw2">\n')
	dat.append('\t\t<parameter free="1" max="100000" min="1e-06" name="Integral" scale="1e-11" value="7" />\n')
	dat.append('\t\t<parameter free="1" max="1" min="-4" name="Index" scale="1" value="-2.2" />\n')
	dat.append('\t\t<parameter free="0" max="200000" min="20" name="LowerLimit" scale="1" value="%s" />\n'%Inputs['emin'][1])
	dat.append('\t\t<parameter free="0" max="200000" min="20" name="UpperLimit" scale="1" value="%s" />\n'%Inputs['emax'][1])
	dat.append('\t</spectrum>\n')
	dat.append('\t<spatialModel type="SkyDirFunction">\n')
	dat.append('\t\t<parameter free="0" max="360" min="-360" name="RA" scale="1" value="%s" />\n'%RA)
	dat.append('\t\t<parameter free="0" max="90" min="-90" name="DEC" scale="1" value="%s" />\n'%dec)
	dat.append('\t</spatialModel>\n')
	dat.append('</source>\n')

	return dat



# ###################################################### #
## #################################################### ##
### ################################################## ###
#### ################################################ ####
##### ############################################## #####
###### ############################################ ######
####### ########################################## #######
######## ######################################## ########

def Add_Srces(Input, fl):

	""" Take stock of confused sources, find their locations   """
	f = open('MODEL_all.xml','r')
	dat = f.readlines()
	f.close()
	nm = []
	for line in dat:
		if ('source name' in line):
			nm.append( line.split('name="')[1].split('"')[0] )

	check_cut = 18
	check = []
	for name in nm:
		if ( len(name) == check_cut ):
			check.append( name )

	mdl = []
	if ( len( check ) > 0 ):
		for nm in check:
			for i in np.arange( len(dat) ):
				if ( nm in dat[i] ):
					mdl = dat[i:i+11]

	GLON = []
	GLAT = []
	for j in np.arange( len(mdl) ):
		for line in mdl[j]:
			if ('RA' in line):
				RA = float(line.split('value="')[1].split('"')[0])
			elif ('DEC' in line):
				DE = float(line.split('value="')[1].split('"')[0])

		(lon,lat) = astCoords.convertCoords("J2000","GALACTIC",RA,DE,2000)
		GLON.append(lon)
		GLAT.append(lat)

	# Go through clouds. Check to see if CO clouds are close to any confused
	# sources.
	RAadd = []
	DEadd = []
	dirc = os.listdir('../Clouds_2sig')
	for nm in dirc:
		if (nm[:2] == 'CO'):
			close = 0
			cld = pyfits.open('../Clouds_2sig/' + nm )
			cldwcs = astWCS.WCS('../Clouds_2sig' + nm )
			if ( cld[0].data.max() > 3 ):
				(xx,yy) = np.where( cld[0].data == cld[0].data.max() )
				cld[0].data[ cld[0].data > 0 ] = 1
				for j in np.arange( len(GLON) ):
					(gln,glt) = cldwcs.wcs2pix(GLON[j],GLAT[j])
					if (np.sqrt( (xx-gln)**2 + (yy-glt)**2) < 15):
						close += 1

				if ( (close == 0) and (cld[0].data.sum() > 50) ):
					(l,b) = cldwcs.pix2wcs(xx,yy)
					(RA1,DE1) = astCoords.convertCoords("GALACTIC","J2000",l,b,2000)
					RAadd.append(RA1)
					DEadd.append(DE1)

#
#
#
#
# ################################################# #
def UpperLim(Input,like,nm):
	""" Determine if we need to calculate upper limit on flux and
	 	calculate it using the Upper Limits package. Use only after
	 	fitting with NewMinuit. """

	ul = UpperLimits(like)
	if (like.Ts(nm) < 20):
		try:
			ul[nm].compute(emin=Input['emin'][1],emax=Input['emax'][1])
		except:
			return 0


	return ul[nm].results[0]


# ################################################# #
def Pl2_MDL(like,Input):
	""" Here we change the spectrum of HI in the model to be
	power law 2. """
	integ = like.model['HI'].funcs['Spectrum'].getParam('Integral').value()
	like.setSpectrum('HI','PowerLaw2')
	like.model['HI'].funcs['Spectrum'].params['Integral'].setValue(integ)    
	like.model['HI'].funcs['Spectrum'].params['Integral'].setScale(1e-27)
	like.model['HI'].funcs['Spectrum'].params['Integral'].setBounds([1e-6,1e5])
	like.model['HI'].funcs['Spectrum'].params['Index'].setValue(-2.0)
	like.model['HI'].funcs['Spectrum'].params['Index'].setScale(1)
	like.model['HI'].funcs['Spectrum'].params['Index'].setBounds([-4,1])
	like.model['HI'].funcs['Spectrum'].params['LowerLimit'].setValue(Input['emin'][1])
	like.model['HI'].funcs['Spectrum'].params['UpperLimit'].setValue(Input['emax'][1])

	return like

# ################################################# #
def BPL2_MDL(like,Input,fl):
	""" Change spectrum of CO and E(B-V) in the model
	to be broken power law 2. This is run after an
	initial fitting where CO & E(B-V) are single
	power laws. """

	if ( (fl['mode'].lower() == 'all') or (fl['mode'].lower() == 'nodg') or (fl['mode'].lower() == 'all_plus')):
		integ = like.model['CO'].funcs['Spectrum'].getParam('Integral').value()
		like.setSpectrum('CO','BrokenPowerLaw2')
		like.model['CO'].funcs['Spectrum'].params['Integral'].setValue(integ)    
		like.model['CO'].funcs['Spectrum'].params['Integral'].setScale(1e-6)
		like.model['CO'].funcs['Spectrum'].params['Integral'].setBounds([1e-6,1e5])
		like.model['CO'].funcs['Spectrum'].params['Index1'].setValue(-2.0)
		like.model['CO'].funcs['Spectrum'].params['Index1'].setBounds([-4,1])
		like.model['CO'].funcs['Spectrum'].params['Index1'].setScale(1)
		like.model['CO'].funcs['Spectrum'].params['Index2'].setValue(-2.8)
		like.model['CO'].funcs['Spectrum'].params['Index2'].setScale(1)
		like.model['CO'].funcs['Spectrum'].params['Index2'].setBounds([-4,1])
		like.model['CO'].funcs['Spectrum'].params['BreakValue'].setValue(1100)
		like.model['CO'].funcs['Spectrum'].params['BreakValue'].setBounds([100,20000])
		like.model['CO'].funcs['Spectrum'].params['LowerLimit'].setValue(Input['emin'][1])
		like.model['CO'].funcs['Spectrum'].params['UpperLimit'].setValue(Input['emax'][1])
	# Now for E(B-V)

	if ( (fl['mode'].lower() == 'all') or (fl['mode'].lower() == 'noco') or (fl['mode'].lower() == 'all_plus') ):
		integ = like.model['E(B-V)'].funcs['Spectrum'].getParam('Integral').value()
		like.setSpectrum('E(B-V)','BrokenPowerLaw2')
		like.model['E(B-V)'].funcs['Spectrum'].params['Integral'].setValue(integ)    
		like.model['E(B-V)'].funcs['Spectrum'].params['Integral'].setScale(1e-6)
		like.model['E(B-V)'].funcs['Spectrum'].params['Integral'].setBounds([1e-6,1e5])
		like.model['E(B-V)'].funcs['Spectrum'].params['Index1'].setValue(-2.0)
		like.model['E(B-V)'].funcs['Spectrum'].params['Index1'].setScale(1)
		like.model['E(B-V)'].funcs['Spectrum'].params['Index1'].setBounds([-4,1])
		like.model['E(B-V)'].funcs['Spectrum'].params['Index2'].setValue(-2.8)
		like.model['E(B-V)'].funcs['Spectrum'].params['Index2'].setScale(1)
		like.model['E(B-V)'].funcs['Spectrum'].params['Index2'].setBounds([-4,1])
		like.model['E(B-V)'].funcs['Spectrum'].params['BreakValue'].setValue(1100)
		like.model['E(B-V)'].funcs['Spectrum'].params['BreakValue'].setBounds([100,20000])
		like.model['E(B-V)'].funcs['Spectrum'].params['LowerLimit'].setValue(Input['emin'][1])
		like.model['E(B-V)'].funcs['Spectrum'].params['UpperLimit'].setValue(Input['emax'][1])

	return like

# ############################################## #
# ###  Run this if we want to mask all point ### #
# ###  sources 				     ### #
# ############################################## #
def MaskSrc(Inputs,fl,like):
	""" We go through the likelihood model and mask the point
	sources in the CMAP, CCUBE, and bexpmap. 
	Input = input parameters, from data selection stage
	fl    = list of filenames used for analysis
	like  = likelihood object, used to get point sources
		to mask 
	-----------------------------------------------------
	Outputs:
	    Writes a masked binned exposure map and masked
	    counts map and counts cube to use to later analysis.
	    Also removes all point sources from the likelihood
	    objects and returns the object 'like'.
	"""
	bx = pyfits.open(fl['bexpmp'])
	cm = pyfits.open(fl['CMAP'])
	ccb= pyfits.open(fl['CCUBE'])
	co = pyfits.open(fl['COmpnon'])
	ebv= pyfits.open(fl['Ebv_%snon'%fl['opacity']])
	hi = pyfits.open(fl['HImp_%snon'%fl['opacity']])
	hi_= pyfits.open(fl['HImp_%snon_gt'%fl['opacity']])
	cat = fl['cat_var']
	cmwcs = astWCS.WCS( cm[0].header, mode="pyfits" )
	bxwcs = astWCS.WCS( bx[0].header, mode="pyfits" )
	cowcs = astWCS.WCS( co[0].header, mode="pyfits" )
	Eslp = (np.log10(Inputs['emax'][1])-np.log10(Inputs['emin'][1]))/ccb[0].data.shape[0]
	En = 10**(0.5*(np.log10(Inputs['emin'][1])+np.log10(Inputs['emax'][1])))

	# Prep for gas templates
	co1 = np.zeros( (ccb[0].header['NAXIS3'],co[0].header['NAXIS2'],
			 		co[0].header['NAXIS1']) )
	ebv1 = np.zeros( (ccb[0].header['NAXIS3'],co[0].header['NAXIS2'],
					co[0].header['NAXIS1']) )
	hi1 = np.zeros( (ccb[0].header['NAXIS3'],co[0].header['NAXIS2'],
					co[0].header['NAXIS1']) )
	hi2 = np.zeros( (ccb[0].header['NAXIS3'],co[0].header['NAXIS2'],
					co[0].header['NAXIS1']) )
	for i in np.arange( ccb[0].header['NAXIS3'] ):
		co1[i,:,:] = co[0].data
		ebv1[i,:,:]= ebv[0].data
		hi1[i,:,:] = hi[0].data
		hi2[i,:,:] = hi_[0].data

	# Cycle through names and remove....
	for name in like.sourceNames():
		if ('Add' in name ):
			ra = like.model[name].funcs['Position'].params['RA'].getValue()
			de = like.model[name].funcs['Position'].params['DEC'].getValue()

			(glon,glat) = astCoords.convertCoords("J2000","GALACTIC",ra,de,2000)
			(cm,ccb,bx) = Mask_RM(fl,glon,glat,-1,cm,ccb,
				  		bx,cmwcs,bxwcs,Eslp,En)
			(co1,ebv1,hi1,hi2) = TempMask_RM(fl,glon,glat,cowcs,co1,ebv1,hi1,hi2,-1)
			like.deleteSource(name)
		elif ( '3FGL' in name ):
			nm = name[1:5] + ' ' + name[5:]
			ind = np.where( cat[1].data['Source_Name'] == nm )[0][0] 
			(glon,glat) = ( cat[1].data['GLON'][ind],
					cat[1].data['GLAT'][ind] )
			(cm,ccb,bx) = Mask_RM(fl,glon,glat,ind,cm,ccb,
				  	bx,cmwcs,bxwcs,Eslp,En)
			(co1,ebv1,hi1,hi2) = TempMask_RM(fl,glon,glat,cowcs,co1,ebv1,hi1,hi2,ind)
			like.deleteSource(name)


	os.system('mv %s bexpmap_nomsk.fits'%fl['bexpmp'])
	os.system('mv %s CMAP_nomsk.fits'%fl['CMAP'])
	os.system('mv %s CCUBE_nomsk.fits'%fl['CCUBE'])
	bx.writeto(fl['bexpmp'])
	ccb.writeto(fl['CCUBE'])
	cm.writeto(fl['CMAP'])

	os.system('mv %s ../CO_temp_nomsk.fits'%fl['COmpnon'])
	os.system('mv %s ../Ebv_%s_nomsk.fits'%(fl['Ebv_%snon'%fl['opacity']],fl['opacity']))
	os.system('mv %s ../HI_%s_nomsk.fits'%(fl['HImp_%snon'%fl['opacity']],fl['opacity']))
	os.system('mv %s ../HI_%s_gt_nomsk.fits'%(fl['HImp_%snon_gt'%fl['opacity']],fl['opacity']))

	com= pyfits.PrimaryHDU( co1 )
	him= pyfits.PrimaryHDU( hi1 )
	hi_m=pyfits.PrimaryHDU( hi2 )
	ebvm=pyfits.PrimaryHDU( ebv1)
	com.header = co[0].header
	com.header['CTYPE3'] = ('ENERGY')
	com.header['CRPIX3'] = ( 1, 'Reference Pixel')
	com.header['CRVAL3'] = (Inputs['emin'][1], 'Energy at the reference pixel')
	com.header['CDELT3'] = (ccb[0].header['CDELT3'],'Z-axis incr per pixel of physical coord at posi')
	com.header['CUNIT3'] = ('MeV','Physical Units for z-axis')

	him.header = com.header
	ebvm.header= com.header
	hi_m.header= com.header

	return "Finished Masking"

####################################################################
def Mask_RM(fl, glon,glat,ind,cm,ccb,bx,cmwcs,bxwcs,Eslp,En):
	""" Remove ring around masked source ... just so I don't have to
	write this code twice... """
	cat = fl['cat_var']

	A = np.zeros( (ccb[0].data.shape[2],ccb[0].data.shape[1]) )
	# CMAP/CCUBE
	for i in np.arange( ccb[0].data.shape[2] ):
		for j in np.arange( ccb[0].data.shape[1] ):
			(l,b) = cmwcs.pix2wcs(i,j)
			A[j,i] = angsep(l,b,glon,glat)

	# bexpmap
	B = np.zeros( (bx[0].data.shape[2],bx[0].data.shape[1]) )
	for i in np.arange( bx[0].data.shape[2] ):
		for j in np.arange( bx[0].data.shape[1] ):
			(l,b) = bxwcs.pix2wcs(i,j)
			B[j,i] = angsep(l,b,glon,glat)

	#CMAP
	if ( ind == -1 ):
		rad = PSF_approx(En)
	else:
		rad = cat[1].data['Conf_95_SemiMajor'][ind] + PSF_approx(En)

	coor = np.where( A < rad )
	cm[0].data[coor] = 0

	#CCUBE/bexpmap
	for k in np.arange( ccb[0].data.shape[0] ):
		E = 10**( np.log10( ccb[0].header['CRVAL3'])+k*Eslp )
		if ( ind == -1 ):
			rad = PSF_approx(E)
		else:
			rad = cat[1].data['Conf_95_SemiMajor'][ind] + PSF_approx(E)

		coor = np.where( A < rad )
		ccb[0].data[k,coor[0],coor[1]] = 0
		coor = np.where( B < rad )
		bx[0].data[k,coor[0],coor[1]] = 0

	try:
		k += 1
		E = 10**( np.log10( ccb[0].header['CRVAL3'])+k*Eslp )
		if ( ind == -1 ):
			rad = PSF_approx(E)
		else:
			rad = cat[1].data['Conf_95_SemiMajor'][ind] + PSF_approx(E)

		coor = np.where( B < rad )
		bx[0].data[k,coor[0],coor[1]] = 0
	except:
		print "Bexpmap not that big!"
		pass

	return (cm,ccb,bx)

####################################################################
def TempMask_RM(fl,glon,glat,cowcs,co,ebv,hi,hi_,ind):


	A = np.ones( (co.shape[1],co.shape[2]) )
	for i in np.arange( co[0].data.shape[2] ):
		for j in np.arange( co[0].data.shape[1] ):
			(l,b) = cmwcs.pix2wcs(i,j)
			A[j,i] = angsep(l,b,glon,glat)

	for k in np.arange( co.shape[0] ):
		E = 10**( np.log10( ccb[0].header['CRVAL3'])+k*Eslp )
		if ( ind == -1 ):
			rad = PSF_approx(E)
		else:
			rad = cat[1].data['Conf_95_SemiMajor'][ind] + PSF_approx(E)

		coor = np.arange( A < rad )
		co[k,coor[0],coor[1]] = 0
		hi[k,coor[0],coor[1]] = 0
		ebv[k,coor[0],coor[1]]= 0
		hi_[k,coor[0],coor[1]]= 0


	return (co,ebv,hi,hi_)

######################################
def PSF_approx(En):
	E = np.array([30, 55, 95, 180, 300, 560, 950, 1080, 3000, 5600, 9500, 18000, 30000])
	E = np.log10(E)
	PSF_95=np.array([8.5, 5.5, 3.8, 2.3, 1.5, 1, 0.63, 0.4, 0.28, 0.2, 0.13, 0.1, 0.085])
	if (En < 30):
		return 10
	elif (En>30000):
		return 0.1
	else:
		ind = max(np.where(E<np.log10(En))[0])
		slope = (PSF_95[ind+1]-PSF_95[ind])/(E[ind+1]-E[ind])
	return ( slope*(np.log10(En)-E[ind]) + PSF_95[ind] )
