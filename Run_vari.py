#!/home/abrahams/epd_free-7.3-2-rh5-x86/bin/python

import sys, os
import gt_apps as gaps
import numpy as np
from astLib import astCoords
import pyfits
from astLib import astWCS
import healpy as hp
import make2FGLxml
from scipy import ndimage

import pyLikelihood
from BinnedAnalysis import *
from SelectLoad import *
# ########################
from ModelMake_v2 import *
from LikeRun import *


glon = [159]
glat = [-34]
emin = [250]
emax = [10000]
ebns = [25]
pre = '/home/abrahams/HICO_survey/SourceSearch/'
# Inputs = ReadFile('Config.txt')
for m in np.arange( len(glon) ):
    ndir = pre +(('l%sb%s')%(int(glon[m]),int(glat[m])))
    try:
	os.mkdir(ndir)
    except:
	print "Directory already exists ... continuing"
	pass

#    Time_Analysis(glon[m],glat[m],emin[m],emax[m],ebns[m],pre)
    os.chdir(ndir)
    ########################
    ## Generate time bins ##
    tmn = 239557417
    tmx = 392574813
#    dt  = 1.57785e7	# 3.15569e7 one year time bins (in seconds)
#			# this is 6 months
#    dt = 12096000      # 20 weeks in seconds ~4 months
#    dt = 6048000	# 10 weeks in seconds
    dt = 9072000	# 15 weeks in seconds
    tmin = []
    tmax = []
    i = 0
    while ( (tmn + i*dt) < tmx):
#    for i in np.arange( 5 ):
	tmin.append( tmn + i*dt )
	tmax.append( tmn + (i+1)*dt )
	i += 1

    tmax[len(tmax)-1] = tmx
    ##                    ##
    ########################
    # Brute force, I'm not very elegant...
    tmin_comb = [tmin[4],tmin[10]]
    tmax_comb = [tmin[5],tmax[11]]
    week  = [4,10]
    for mind in np.arange( len(tmin) ):
	# The following: check to see if good directories exist.
	# Skip is it does not...
#	if ('time%s'%mind in os.listdir('.') ):
#	    print "Starting on time%s"%mind
#	    ins = (['glon=%s'%glon[m], 'glat=%s'%glat[m], 'rad=10', 'zcut=True', 
#		'emin=%s'%emin[m], 'emax=%s'%emax[m], 'tmin=%s'%tmin[mind], 'tmax=%s'%tmax[mind], 
#		'nxpix=140', 'nypix=140', 'binsz=0.1', 'enumbins=%s'%ebns[m]])
#   	    Inputs = ReadParams(ins)
#  	    os.mkdir( ndir + '/time%s'%(mind) )
#    	    os.chdir( ndir + '/time%s'%(mind) )
#	elif ( 'time%s+%s'%(mind,mind+1) in os.listdir('.') ):
#	    print "Starting on time%s+%s"%(mind,mind+1)
#	    ins = (['glon=%s'%glon[m], 'glat=%s'%glat[m], 'rad=10', 'zcut=True', 
#		    'emin=%s'%emin[m], 'emax=%s'%emax[m], 'tmin=%s'%tmin[mind], 
#		    'tmax=%s'%tmax[mind+1], 'nxpix=140', 'nypix=140', 'binsz=0.1', 
#		    'enumbins=%s'%ebns[m]])
#	    Inputs = PreadParams(ins)
#	    os.mkdir( ndir + '/time%s+%s'%(mind,mind+1) )
#	    os.chdir( ndir + '/time%s+%s'%(mind,mind+1) )
#	else:
#	    print "Skipping bin number %s"%mind
#	    continue
	if ( mind not in week ):
	    print "Starting on time%s"%mind
	    ins = (['glon=%s'%glon[m], 'glat=%s'%glat[m], 'rad=10', 'zcut=True', 
		'emin=%s'%emin[m], 'emax=%s'%emax[m], 'tmin=%s'%tmin[mind], 'tmax=%s'%tmax[mind], 
		'nxpix=140', 'nypix=140', 'binsz=0.1', 'enumbins=%s'%ebns[m]])
    	    Inputs = ReadParams(ins)
    	    os.mkdir( ndir + '/time%s'%(mind) )
    	    os.chdir( ndir + '/time%s'%(mind) )
	elif ( mind in week ):
	    print "Starting on time%s+%s"%(mind,mind+1)
	    ins = (['glon=%s'%glon[m], 'glat=%s'%glat[m], 'rad=10', 'zcut=True', 
		    'emin=%s'%emin[m], 'emax=%s'%emax[m], 'tmin=%s'%tmin[mind], 
		    'tmax=%s'%tmax[mind+1], 'nxpix=140', 'nypix=140', 'binsz=0.1', 
		    'enumbins=%s'%ebns[m]])
	    Inputs = PreadParams(ins)
	    os.mkdir( ndir + '/time%s+%s'%(mind,mind+1) )
	    os.chdir( ndir + '/time%s+%s'%(mind,mind+1) )
	elif ( mind in (week + np.ones(len(A))) ):
	    print "Skipping bin number %s"%mind
	    continue
    	Flname = FileNames(Inputs,pre,10)
	Flname['Tspin'] = 125
	cat = pyfits.open(Flname['cat'])
	Flname['cat_var'] = cat
    	Select(Inputs,Flname)
    	mktime(Inputs,Flname)

    	# Cycle through the gtbin options, used for diagnostic
    	# purposes.
    	alg = ['CMAP','LC','CCUBE']
    	for mode in alg:
	    Flname['binalg'] = mode
	    binning(Inputs,Flname)

    	# Write "Config.txt" file necessary to reproduce analysis
    	f = open('Config.txt','w')
    	for elem in ins:
	    f.write(elem + '\n')
	
    	f.close()
    
	#   ### decide whether to make a livetime cube or not ###       #
	expcube(Inputs,Flname,True)
	Flname['opacity'] = 'thin'

#	HI_map(Inputs,Flname,mp)
#	CO_map(Inputs,Flname)
#	Avmap_make(Inputs,Flname)
#	os.system('mv CO_temp_norm.fits ../.')
#	os.system('mv HI_%s_norm.fits ../.'%mp)
#	os.system('mv CO_temp.fits ../.')
#	os.system('mv HI_%s.fits ../.'%mp)
#	os.system('mv CO_temp_error.fits ../.')
#	Dust_resid(Inputs,Flname,mp)
#	os.system('mv Ebv_%s.fits ../.'%mp)
#	os.system('mv Ebv_%s_norm.fits ../.'%mp)
    	# Make a model -- with CO and dark gas and everything
    	#              -- with no CO and DG
    	#              -- replace CO/DG with a point source

    	Bexpmap_mk(Inputs,Flname,True)
   	mod = make2FGLxml.srcList(Flname['cat'],Flname['mkt_fl'],'MODEL.xml')
    	mod.makeModel(Flname['gal'],'gal_2yearp7v6_v0',Flname['iso'],'iso_p7v6source')
	Flname['mode'] = 'orig'
	CatSrc_Fix(Inputs,Flname,True)
    	mode = ['all','NoCODG','Ptsrc','fx_src']
	for Flname['mode'] in mode:
	    Mdl_Make(Inputs,Flname,(ndir + '/ebin%s_%s'%(emin,emax) + '/'))
	    try:
		os.mkdir('%s'%Flname['mode'])
	    except:
		pass

	    Srcmap_mk(Inputs,Flname,True)
	    BinLike(Inputs, Flname,True)

    	os.chdir( ndir )

    print "Done variability analysis!"


print "Done ALL variability analyses!"
