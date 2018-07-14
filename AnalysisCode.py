#!/home/abrahams/epd_free-7.3-2-rh5-x86/bin/python

import sys, os
import gt_apps as gaps
import numpy as np
from astLib import astCoords
import pyfits
from astLib import astWCS
import healpy as hp
import scipy.optimize
from scipy import ndimage
import pyLikelihood
from BinnedAnalysis import *
from UpperLimits import UpperLimits
from scipy import integrate

# Here we are going to load in the parameters needed to select
# some of the data: we need the following input:
# glon, glat, rad, emin, emax, zcut, tmin, tmax
## ##############################################
# Only glon, glat, rad are necessary. Other parameters will
#     be set to defaults if not provided, as explained below
# Similarly, a Config.txt file can be passed to this routine.
## ##############################################
# Defaults:
# If emin or emax are empty, assume default: E = (250,30000)
# zcut is True/False: True means perform zenith cut during 
#     the gtmktime step, False: during the livetime cube
#     calculation. Default is True
# tmin and tmax will default to standard values if not given:
# tmin = 239557417, tmax = 455067824 (old:392574813)
## ##############################################
# Arguments will be passed as strings and then broken up into
# the necessary types. Example:
# ipython SelectLoad glon=90 glat=-40 rad=15
## ##############################################
def ReadParams(ins):
	""" Read in the inputs from command line. """
	arg=[]
	inp = {}
	for elem in ins:
		arg.append( elem.split('=') )

	try:
		inp['glon'] = [el for el in arg if el[0] == 'glon'][0][1]
		inp['glat'] = [el for el in arg if el[0] == 'glat'][0][1]
		inp['rad']  = [el for el in arg if el[0] == 'rad'][0][1]
	except:
		sys.exit("Invalid data selection parameters")

	try:
		inp['nxpix'] = [el for el in arg if el[0] == 'nxpix'][0][1]
		inp['nypix'] = [el for el in arg if el[0] == 'nypix'][0][1]
		inp['binsz'] = [el for el in arg if el[0] == 'binsz'][0][1]
		inp['enumbins'] = [el for el in arg if el[0] == 'enumbins'][0][1]
	except:
		sys.exit("Invalid binning parameters")

	try:
		inp['emin'] = [el for el in arg if el[0] == 'emin'][0][1]
	except:
		inp['emin'] = ['emin','250'][1]

	try:
		inp['emax'] = [el for el in arg if el[0] == 'emax'][0][1]
	except:
		inp['emax'] = ['emax','30000'][1]

	try:
		inp['zcut'] = [el for el in arg if el[0] == 'zcut'][0][1]
	except:
		inp['zcut'] = ['zcut','True'][1]

	try:
		inp['tmin'] = [el for el in arg if el[0] == 'tmin'][0][1]
	except:
		inp['tmin'] = ['tmin','239557417'][1]

	try:
		inp['tmax'] = [el for el in arg if el[0] == 'tmax'][0][1]
	except:
		inp['tmax'] = ['tmax','392574813'][1]

	return ValConvert(inp)

## ##############################################
## Read parameters from file, if option is taken
def ReadFile(fl):
	""" Read inputs from file 'Config.txt' """
	f = open(fl,'r')
	arg = f.readlines()
	f.close()
	inp = {}
	args = []
	for elem in arg:
		el = elem
		args.append( [el.split('=')[0], el.split('=')[1][:-1]] )

	inp['glon'] = [el for el in args if el[0] == 'glon'][0][1]
	inp['glat'] = [el for el in args if el[0] == 'glat'][0][1]
	inp['rad']  = [el for el in args if el[0] == 'rad'][0][1]
	inp['emin'] = [el for el in args if el[0] == 'emin'][0][1]
	inp['emax'] = [el for el in args if el[0] == 'emax'][0][1]
	inp['zcut'] = [el for el in args if el[0] == 'zcut'][0][1]
	inp['tmin'] = [el for el in args if el[0] == 'tmin'][0][1]
	inp['tmax'] = [el for el in args if el[0] == 'tmax'][0][1]
	inp['nxpix'] = [el for el in args if el[0] == 'nxpix'][0][1]
	inp['nypix'] = [el for el in args if el[0] == 'nypix'][0][1]
	inp['binsz'] = [el for el in args if el[0] == 'binsz'][0][1]
	inp['enumbins'] = [el for el in args if el[0] == 'enumbins'][0][1]
	return ValConvert(inp)

## ##############################################
## Convert the inputted strings to numbers/boolean
## as necessary for use later in the analysis.
def ValConvert(inp):
	""" Converts the inputted strings into floating point numbers. """
	inp['glon'] = float(inp['glon'])
	inp['glat'] = float(inp['glat'])
	inp['rad'] = float(inp['rad'])
	inp['emin'] = float(inp['emin'])
	inp['emax'] = float(inp['emax'])
	inp['tmin'] = float(inp['tmin'])
	inp['tmax'] = float(inp['tmax'])
	inp['zcut'] = str2bool(inp['zcut'])
	inp['nxpix'] = float(inp['nxpix'])
	inp['nypix'] = float(inp['nypix'])
	inp['binsz'] = float(inp['binsz'])
	inp['enumbins'] = float(inp['enumbins'])
	return inp



## ##############################################
## Convert string to boolean
def str2bool(v):
	""" Convert string to boolean, catches common strings for true/false. """
	if (v.lower() in ('yes','true','t','y','1')):
		return True
	elif (v.lower() in ('no','false','f','n','0')):
		return False

## ##############################################
def Select(Inputs,fl):
	""" Select data centered around specified coordinates. 
	
	Input: input data, ROI position, energy range, etc.
	fl   : list of file names from SelectLoad.FileNames"""
	(ra,de) = astCoords.convertCoords("GALACTIC","J2000",
		Inputs['glon'],Inputs['glat'],2000)
#	gaps.filter['evclass']=2		# P7REP and earlier
	gaps.filter['evclass'] = 128 # source = 128, ultracleanveto = 1024
	gaps.filter['evtype'] = 3
	gaps.filter['ra']=ra
	gaps.filter['dec']=de
	gaps.filter['rad']=Inputs['rad']
	gaps.filter['zmax']=90
	gaps.filter['tmin']=Inputs['tmin']
	gaps.filter['tmax']=Inputs['tmax']
	gaps.filter['infile']=fl['files']
	gaps.filter['outfile']=fl['sel_fl']
	gaps.filter['emin']=Inputs['emin']
	gaps.filter['emax']=Inputs['emax']

	if (fl['sel_fl'] in os.listdir('.')):
		return "Already ran gtselect, continue"
	else:
		gaps.filter.run()
		return "Done gtselect"

## ##############################################
def mktime(Inputs,fl):
	""" Filter selected data further: rock angle < 52 deg, etc. Follows
	the standard selection cuts recommended by LAT team. """
	gaps.maketime['scfile']=fl['SC_fl']
	gaps.maketime['filter']='DATA_QUAL==1 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52'
	if (Inputs['zcut'] == True):    
		gaps.maketime['roicut'] = 'yes'
	else:
		gaps.maketime['roicut'] = 'no'

	gaps.maketime['evfile']=fl['sel_fl']
	gaps.maketime['outfile']=fl['mkt_fl']
	if (fl['mkt_fl'] in os.listdir('.')):
		return "Already ran gtmktime, continue"
	else:
		gaps.maketime.run()
		return

## ##############################################
def FileNames(Inputs,pre,velcut,Tspin):
	""" Stores file names for the likelihood analysis """
	fl = {}
	# pre = '/home/abrahams/HICO_survey/SourceSearch/'
	fl['files'] = pre + 'all_sky_files.txt'
	fl['mkt_fl']='mkt.fits'
	fl['SC_fl']= pre + 'AllSky_p8/SC.fits'
#    fl['SC_fl']= pre + 'AllSky/SC.fits'
 #   fl['IRFS'] = 'P7SOURCE_V6'
 #   fl['IRFS'] = 'P7REP_SOURCE_V15'
#### Using Pass 8 ####
	fl['IRFS'] = 'P8R2_SOURCE_V6'
######################
	fl['sel_fl']='sel.fits'
	fl['CMAP'] = 'CMAP.fits'
	fl['CMAP_big'] = 'CMAP_big.fits'
	fl['CCUBE']= 'CCUBE.fits'
	fl['LC' ]  = 'LC.fits'
	fl['ltcube']= pre + ('Analysis/ltCube_%s' % (int(Inputs['rad']))) + 'deg.fits'
	fl['ltcube_t']= 'ltCube.fits'
	fl['bexpmp']='bexpmap.fits'
	fl['IC']   = '/home/abrahams/GALPROP/results_54_04d90003/ics_isotropic_mapcube_54_04d90003'
	fl['Dust_allsky'] = pre + 'DustModel_Planck.fits'
	fl['Dust'] = pre + 'l%sb%s/Dust.fits'%(int(Inputs['glon']),int(Inputs['glat']))
#    fl['iso']  = os.environ['FERMI_DIR']+ '/refdata/fermi/galdiffuse/' + 'iso_source_v05.txt'
#    fl['gal']  = os.environ['FERMI_DIR']+ '/refdata/fermi/galdiffuse/' + 'gll_iem_v05_rev1.fit'
	fl['iso'] = os.environ['FERMI_DIR'] + '/refdata/fermi/galdiffuse/' + 'iso_P8R2_SOURCE_V6_v06.txt'
	fl['gal'] = os.environ['FERMI_DIR'] + '/refdata/fermi/galdiffuse/' + 'gll_iem_v06.fits'

	fl['COplanck']= os.environ['FERMI_DIR']+ '/refdata/fermi/galdiffuse/' + 'CO_PlanckDR2_Type2.fits'
	fl['COplanck_type3'] = pre + 'CO_planck_type3.fits'
	fl['cat']  = os.environ['FERMI_DIR']+ '/refdata/fermi/' + 'gll_psc_v14.fit'

	fl['MODEL_newsrc'] = 'MODEL_newsrc.xml'
	fl['srcmp_int1'] = 'srcmap_int1.fits'
	fl['MODEL_int1'] = 'mdl1_int1.xml'
    # Models, source maps, and output models...
	# Generic

    # for default Fermi galactic diffuse map
	fl['MODEL_orig'] ='MODEL.xml'
	fl['srcmp_orig']= 'orig/srcmap_orig.fits'
	fl['outmdl1_orig']='MODEL_.xml'
	fl['outmap1_orig']='orig/mdlmap1_orig.fits'
    # HI + CO + DG
	fl['MODEL_all']= 'MODEL_all.xml'
	fl['srcmp_all']= 'srcmap_all.fits'
	fl['outmdl1_all']='mdl1_all.xml'
	fl['outmdl2_all']='mdl2_all.xml'
	fl['outmap1_all']='mdlmap1_all.fits'
	fl['outmap2_all']='mdlmap2_all.fits'
    # HI only
	fl['MODEL_nocodg']= 'nocodg/MODEL_nocodg.xml'
	fl['srcmp_nocodg']= 'nocodg/srcmap_srcmap.fits'
	fl['outmdl1_nocodg']='nocodg/mdl1_nocodg.xml'
	fl['outmdl2_nocodg']='nocodg/mdl2_nocodg.xml'
	fl['outmap1_nocodg']='nocodg/mdlmap1_nocodg.fits'
	fl['outmap2_nocodg']='nocodg/mdlmap2_nocodg.fits'
    # HI + point source located at W_CO peak
	fl['MODEL_ptsrc']= 'ptsrc/MODEL_ptsrc.xml'
	fl['srcmp_ptsrc']= 'ptsrc/srcmap_ptsrc.fits'
	fl['outmdl1_ptsrc']='ptsrc/mdl1_ptsrc.xml'
	fl['outmdl2_ptsrc']='ptsrc/mdl2_ptsrc.xml'
	fl['outmap1_ptsrc']='ptsrc/mdlmap1_ptsrc.fits'
	fl['outmap2_ptsrc']='ptsrc/mdlmap2_ptsrc.fits'
    # HI + CO + DG + point source at W_CO peak
	fl['MODEL_all_plus']= 'all_plus/MODEL_all_plus.xml'
	fl['srcmp_all_plus']= 'all_plus/srcmap_all_plus.fits'
	fl['outmdl1_all_plus']='all_plus/mdl1_all_plus.xml'
	fl['outmdl2_all_plus']='all_plus/mdl2_all_plus.xml'
	fl['outmap1_all_plus']='all_plus/mdlmap1_all_plus.fits'
	fl['outmap2_all_plus']='all_plus/mdlmap2_all_plus.fits'
    # HI + point source, where everything but the point source 
    # is fixed
	fl['MODEL_fx_src']='MODEL_all.xml'
	fl['srcmp_fx_src']= 'fx_src/srcmap_fx_src.fits'
	fl['outmdl1_fx_src']='fx_src/mdl1_fx_src.xml'
	fl['outmdl2_fx_src']='fx_src/mdl2_fx_src.xml'
	fl['outmap1_fx_src']='fx_src/mdlmap1_fx_src.fits'
	fl['outmap2_fx_src']='fx_src/mdlmap2_fx_src.fits'
    # HI + CO
	fl['MODEL_nodg']='nodg/MODEL_nodg.xml'
	fl['srcmp_nodg']= 'nodg/srcmap_nodg.fits'
	fl['outmdl1_nodg']='nodg/mdl1_nodg.xml'
	fl['outmdl2_nodg']='nodg/mdl2_nodg.xml'
	fl['outmap1_nodg']='nodg/mdlmap1_nodg.fits'
	fl['outmap2_nodg']='nodg/mdlmap2_nodg.fits'
    # HI + DG
	fl['MODEL_noco']='noco/MODEL_noco.xml'
	fl['srcmp_noco']= 'noco/srcmap_noco.fits'
	fl['outmdl1_noco']='noco/mdl1_noco.xml'
	fl['outmdl2_noco']='noco/mdl2_noco.xml'
	fl['outmap1_noco']='noco/mdlmap1_noco.fits'
	fl['outmap2_noco']='noco/mdlmap2_noco.fits'
    # nocld
	fl['MODEL_nocld']='nocld/MODEL_nocld.xml'
	fl['srcmp_nocld']= 'nocld/srcmap_nocld.fits'
	fl['outmdl1_nocld']='nocld/mdl1_nocld.xml'
	fl['outmdl2_nocld']='nocld/mdl2_nocld.xml'
	fl['outmap1_nocld']='nocld/mdlmap1_nocld.fits'
	fl['outmap2_nocld']='nocld/mdlmap2_nocld.fits'
    # Broken PL for HI
	fl['MODEL_power_law']='mdl2_all.xml'
	fl['srcmp_power_law']= 'power_law/srcmap_power_law.fits'
	fl['outmdl1_power_law']='power_law/mdl1_power_law.xml'
	fl['outmdl2_power_law']='power_law/mdl2_power_law.xml'
	fl['outmap1_power_law']='power_law/mdlmap1_power_law.fits'
	fl['outmap2_power_law']='power_law/mdlmap2_power_law.fits'

  	fl['COall'] = pre + 'l%sb%s/CO_temp_all.fits'%(int(Inputs['glon']),int(Inputs['glat']))
# velocity +/- 50 kms
	fl['HIthick'] = '/home/abrahams/HICO_survey/HImaps/%skms/HI_thick_Ts%s.fits'%(int(velcut),int(Tspin))
	fl['HIthick_gt']='/home/abrahams/HICO_survey/HImaps/%skms_gt/HI_thick_Ts%s.fits'%(int(velcut),int(Tspin))
	fl['HIthin'] = '/home/abrahams/HICO_survey/HImaps/%skms/HI_thin_Ts125.fits'%int(velcut)
	fl['HIthin_gt']='/home/abrahams/HICO_survey/HImaps/%skms_gt/HI_thin_Ts125.fits'%int(velcut)

	fl['HIthick_125']='/home/abrahams/HICO_survey/HImaps/%skms/HI_thick_Ts125.fits'%int(velcut)
	fl['HIthick_125_gt']='/home/abrahams/HICO_survey/HImaps/%skms_gt/HI_thick_Ts125.fits'%int(velcut)
	fl['HIthick_140']= '/home/abrahams/HICO_survey/HImaps/%skms/HI_thick_Ts140.fits'%int(velcut)
	fl['HIthick_140_gt']='/home/abrahams/HICO_survey/HImaps/%skms_gt/HI_thick_Ts140.fits'%int(velcut)
# Ts = 80K
	fl['HIthick_80']= '/home/abrahams/HICO_survey/HImaps/%skms/HI_thick_Ts80.fits'%int(velcut)
	fl['HIthick_80_gt']='/home/abrahams/HICO_survey/HImaps/%skms_gt/HI_thick_Ts80.fits'%int(velcut)
	fl['HIthin_125']='/home/abrahams/HICO_survey/HImaps/%skms/HI_thin_Ts125.fits'%int(velcut)
	fl['HIthin_125_gt']='/home/abrahams/HICO_survey/HImaps/%skms_gt/HI_thin_Ts125.fits'%int(velcut)
	fl['HIthin_140']='/home/abrahams/HICO_survey/HImaps/%skms/HI_thin_Ts140.fits'%int(velcut)
	fl['HIthin_140_gt']='/home/abrahams/HICO_survey/HImaps/%skms_gt/HI_thin_Ts140.fits'%int(velcut)
	fl['HIthin_80']='/home/abrahams/HICO_survey/HImaps/%skms/HI_thin_Ts140.fits'%int(velcut)
	fl['HIthin_80_gt']='/home/abrahams/HICO_survey/HImaps/%skms_gt/HI_thin_Ts140.fits'%int(velcut)
# HI spatial templates
	fl['HImp_thick'] = pre + 'l%sb%s/HI_thick_norm.fits'%(int(Inputs['glon']),int(Inputs['glat']))
	fl['HImp_thicknon'] = pre + 'l%sb%s/HI_thick.fits'%(int(Inputs['glon']),int(Inputs['glat']))
	fl['HImp_thin'] = pre + 'l%sb%s/HI_thin_norm.fits'%(int(Inputs['glon']),int(Inputs['glat']))
	fl['HImp_thinnon'] = pre + 'l%sb%s/HI_thin.fits'%(int(Inputs['glon']),int(Inputs['glat']))
	fl['HImp_thick_gt'] = pre + 'l%sb%s/HI_thick_gt_norm.fits'%(int(Inputs['glon']),int(Inputs['glat']))
	fl['HImp_thicknon_gt'] = pre + 'l%sb%s/HI_thick_gt.fits'%(int(Inputs['glon']),int(Inputs['glat']))
	fl['HImp_thin_gt'] = pre + 'l%sb%s/HI_thin_gt_norm.fits'%(int(Inputs['glon']),int(Inputs['glat']))
	fl['HImp_thinnon_gt'] = pre + 'l%sb%s/HI_thin_gt.fits'%(int(Inputs['glon']),int(Inputs['glat']))

	fl['COmp'] = pre + 'l%sb%s/CO_temp_norm.fits'%(int(Inputs['glon']),int(Inputs['glat']))
	fl['COmpnon'] = pre + 'l%sb%s/CO_temp.fits'%(int(Inputs['glon']),int(Inputs['glat']))
	fl['COer'] = pre + 'l%sb%s/CO_temp_error.fits'%(int(Inputs['glon']),int(Inputs['glat']))
	fl['COcldmsk'] = pre+'l%sb%s/CO_temp_nocld.fits'%(int(Inputs['glon']),int(Inputs['glat']))
	fl['Ebv_thick']  = pre + 'l%sb%s/Ebv_thick_norm.fits'%(int(Inputs['glon']),int(Inputs['glat']))
	fl['Ebvcldmsk'] = pre + 'l%sb%s/Ebv_thick_nocld.fits'%(int(Inputs['glon']),int(Inputs['glat']))
	fl['Ebv_thicknon']  = pre + 'l%sb%s/Ebv_thick.fits'%(int(Inputs['glon']),int(Inputs['glat']))
	fl['Ebv_thin']  = pre + 'l%sb%s/Ebv_thin_norm.fits'%(int(Inputs['glon']),int(Inputs['glat']))
	fl['Ebv_thinnon']  = pre + 'l%sb%s/Ebv_thin.fits'%(int(Inputs['glon']),int(Inputs['glat']))

	return fl

## ##############################################
def binning(Inputs,fl):
	""" Bin selected data in three ways. (1) Counts map: used primarily
	for visual inspection, not used in full likelihood analysis. 
	(2) Light curve: to be used to remove data during a blazar flare
	for model robustness. (3) Counts cube: data binned in spatial 
	position AND energy. """
	gaps.evtbin['algorithm'] = fl['binalg']
	gaps.evtbin['evfile'] = fl['mkt_fl']
	gaps.evtbin['scfile'] = fl['SC_fl']
	(ra,de) = astCoords.convertCoords("GALACTIC","J2000",
		Inputs['glon'],Inputs['glat'],2000)
	gaps.evtbin['axisrot'] = 0.0

	# Change other things in case we choose lightcurve, counts map,
	#        or counts cube
	if (fl['binalg'].lower() == 'cmap'):
		if (fl['CMAP'] in os.listdir('.')):
			return "CMAP made, continue"
		else:
			gaps.evtbin['nxpix'] = int(Inputs['nxpix'])
			gaps.evtbin['nypix'] = int(Inputs['nypix'])
			gaps.evtbin['binsz'] = Inputs['binsz']
			gaps.evtbin['outfile']=fl['CMAP']
			if (Inputs['glat'] > 70):
				gaps.evtbin['coordsys'] = 'CEL'
				gaps.evtbin['xref'] = ra
				gaps.evtbin['yref'] = de
				gaps.evtbin['proj'] = 'AIT'
				gaps.evtbin.run()
			else:
				gaps.evtbin['coordsys'] = 'GAL'
				gaps.evtbin['xref'] = Inputs['glon']
				gaps.evtbin['yref'] = Inputs['glat']
				gaps.evtbin['proj'] = 'CAR'
				gaps.evtbin.run()

			gaps.evtbin['nxpix'] = int(Inputs['nxpix'] + (10./Inputs['binsz']))
			gaps.evtbin['nypix'] = int(Inputs['nypix'] + (10./Inputs['binsz']))
			gaps.evtbin['outfile']=fl['CMAP_big']
			if (Inputs['glat'] > 70):
				gaps.evtbin['coordsys'] = 'CEL'
				gaps.evtbin['xref'] = ra
				gaps.evtbin['yref'] = de
				gaps.evtbin['proj'] = 'AIT'
				gaps.evtbin.run()
			else:
				gaps.evtbin['coordsys'] = 'GAL'
				gaps.evtbin['xref'] = Inputs['glon']
				gaps.evtbin['yref'] = Inputs['glat']
				gaps.evtbin['proj'] = 'CAR'
				gaps.evtbin.run()

			return 0
	elif (fl['binalg'].lower() == 'lc'):
		gaps.evtbin['outfile'] = fl['LC']
		gaps.evtbin['tbinalg'] = 'LIN'
		gaps.evtbin['tstart'] = Inputs['tmin']
		gaps.evtbin['tstop'] = Inputs['tmax']
		gaps.evtbin['dtime'] = 3592000		# 30 day time bins
		if (fl['LC'] in os.listdir('.')):
			return "Light Curve made, continue"
		else:
			gaps.evtbin.run()
			return 0
	elif (fl['binalg'].lower() == 'ccube'):
		if (fl['CCUBE'] in os.listdir('.')):
			return "CCUBE present, continue"
		else:
			gaps.evtbin['outfile']=fl['CCUBE']
			gaps.evtbin['nxpix'] = int(Inputs['nxpix'])
			gaps.evtbin['nypix'] = int(Inputs['nypix'])
			gaps.evtbin['binsz'] = Inputs['binsz']
			gaps.evtbin['ebinalg'] = 'LOG'
			gaps.evtbin['emin'] = Inputs['emin']
			gaps.evtbin['emax'] = Inputs['emax']
			gaps.evtbin['enumbins'] = int(Inputs['enumbins'])
			if (Inputs['glat'] > 70):
				gaps.evtbin['coordsys'] = 'CEL'
				gaps.evtbin['xref'] = ra
				gaps.evtbin['yref'] = de
				gaps.evtbin['proj'] = 'AIT'
				gaps.evtbin.run()
				return 0
			else:
				gaps.evtbin['coordsys'] = 'GAL'
				gaps.evtbin['xref'] = Inputs['glon']
				gaps.evtbin['yref'] = Inputs['glat']
				gaps.evtbin['proj'] = 'CAR'
				gaps.evtbin.run()
				return 0

## ##############################################
def expcube(Inputs,fl,time):
	""" Create live time cube for selected observations. """
	gaps.expCube['evfile'] = fl['mkt_fl']
	gaps.expCube['scfile'] = fl['SC_fl']
	gaps.expCube['dcostheta']=0.025
	gaps.expCube['binsz']  = 1
	if (time == False):
		gaps.expCube['outfile']= fl['ltcube']
	elif (time == True):
		gaps.expCube['outfile'] = fl['ltcube_t']
	else:
		print "Error in livetime cube input."

	if (Inputs['zcut'] == False):    
		gaps.expCube['zmax'] = 90

	di= '/home/abrahams/HICO_survey/SourceSearch/l%sb%s'%( int(Inputs['glon']),int(Inputs['glat']) )
	if ( (fl['ltcube'] not in os.listdir(di)) and (fl['ltcube_t'] not in os.listdir(di)) ):
		gaps.expCube.run()
		return 0
	else:
		return "Livetime cube already present."




# Here we make the model. We need to sift through the catalog and
# choose strong sources that are close to the center of the RoI.
# Then, we need to determine the spectral shape before writing then
# xml file. This routine also needs to create the spatial templates
# for HI, CO (from Planck), and dust.

## ##############################################
# weighted sum of nearest pixels
def Interpol(im, x, y):
	""" Weighted sum of nearest pixels given a fractional pixel position. """
	dx = x%1.0
	dy = y%1.0
	x1 = int(np.floor(x))
	y1 = int(np.floor(y))
	weight = np.array([(1-dx)*(1-dy),dx*(1-dy),(1-dx)*dy,dx*dy])
	val = np.zeros( 4 )
	val[0] = im[x,y]
	val[1] = im[(x+1),y]
	val[2] = im[x,(y+1)]
	val[3] = im[(x+1),(y+1)]
	return sum(val*weight)

## ##############################################
# Get the pixel value for the LAB survey, needs galactic longitude
# and the number of pixels on that axis.
def find_lon_hi(l, nl):
	if (l < 180):
		lon = 2*(180-l)
	else:
		lon = 2*(540-l)

	return lon

## ##############################################
def GenLon(Input):
	""" Longitude array for background maps """
	rad= ((Input['nxpix']*Input['binsz'])/2.0) + (6./2.0)
	lu = (Input['glon'] + rad)%360
	ld = (Input['glon'] - rad)%360
	df = int(lu/Input['binsz'])

	sz = int(2*rad/Input['binsz'])
	if ( (lu - ld) < 0 ):
		df = (lu - 0.)/Input['binsz']
		l = np.append( np.linspace(lu,0.,df),np.linspace((360-Input['binsz']),ld,int(sz-df)))
	else:
		l = np.linspace(lu,ld,int(2*rad/Input['binsz']))

	return l

## ##############################################
def GenLat(Input):
	""" Latitude array for background maps """
	rad= ((Input['nypix']*Input['binsz'])/2.0) + (6./2.0)
	bu = (Input['glat'] + rad)
	bd = (Input['glat'] - rad)

	return np.linspace(bd,bu,int(2*rad/Input['binsz']))

## ##############################################
def HI_map(Inputs, fl,dist):
	""" Make HI background map using LAB survey, integrated assuming 
	optically thick emission: N(HI) = -log( 1 - (Tb/(Ts - Tbg)) ). 
	Optically thin emission corresponds to Ts --> infinity

	Tb = antenna temperature (data)
	Ts = spin temperature
	Tbg= background temperature, assumed CMB temperature
	     (eventually use 1.4 GHz data, but small difference
	      at high Galactic latitudes)
	dist= 1 is near, 2 is far.
	"""
	import pyfits
	from astLib import astCoords

	if ( dist == 1 ):
		hi = pyfits.open(fl['HI%s'%(fl['opacity'])])
	elif ( dist == 2 ):
		hi = pyfits.open(fl['HI%s_gt'%(fl['opacity'])])
	else:
		print "Error -- no HI maps there."

	hi[0].data = np.transpose( hi[0].data )
	head = pyfits.getheader('CMAP.fits')
	cmap_w = astWCS.WCS(head, mode="pyfits")
	l = GenLon(Inputs)
	b = GenLat(Inputs)
	nl = head['NAXIS1']

	# cycle through coordinates, find pixel values and place it into array
	image = np.zeros( (l.size,b.size) )
	for i in np.arange( l.size ):
		for j in np.arange( b.size ):
			new_l = find_lon_hi(l[i],nl)
			new_b = 2*(90+b[j])
			image[i,j] = Interpol(hi[0].data, new_l, new_b)

	radx = ((Inputs['nxpix']*Inputs['binsz'])/2.0) + (6./2.0)
	rady = ((Inputs['nypix']*Inputs['binsz'])/2.0) + (6./2.0)
	pxcen = radx/Inputs['binsz']
	pycen = rady/Inputs['binsz']
	ima = pyfits.PrimaryHDU(np.transpose(image))
	ima.header['CTYPE1'] = ('GLON-CAR')
	ima.header['CRPIX1'] = ( pxcen, 'Reference Pixel')
	ima.header['CRVAL1'] = (Inputs['glon'], 'Galactic longitude at reference pixel')
	ima.header['CDELT1'] = (-Inputs['binsz'],'x-axis increment per pixel')
	ima.header['CUNIT1'] = ('deg','Physical Units for x-axis')
	ima.header['CTYPE2'] = ('GLAT-CAR')
	ima.header['CRPIX2'] = ( pycen, 'Reference pixel')
	ima.header['CRVAL2'] = (Inputs['glat'],'Galactic latitude at reference pixel')
	ima.header['CDELT2'] = (Inputs['binsz'], 'y-axis increment per pixel')
	ima.header['CUNIT2'] = ('deg','Physical units for y-axis')
	ima.header['Equinox'] = (2000.,'Equinox of coordinates')
	image[ image < 0.0 ] = 0.0
	# we need to normalize the map (making sure degrees are radians)
	if ( dist == 1 ):
		ima.writeto(fl['HImp_%snon'%fl['opacity']])
		hi = pyfits.open(fl['HImp_%snon'%fl['opacity']])
	elif ( dist == 2 ):
		ima.writeto(fl['HImp_%snon_gt'%fl['opacity']])
		hi = pyfits.open(fl['HImp_%snon_gt'%fl['opacity']])

	hi[0].data[ hi[0].data < 0.0 ] = 0.0
	norm = hi[0].data.sum() * (np.pi/180.)**2 * (Inputs['binsz']**2)
	hi[0].data /= norm
	if ( dist == 1 ):
		hi.writeto(fl['HImp_%s'%fl['opacity']])
	elif ( dist == 2 ):
		hi.writeto(fl['HImp_%s_gt'%fl['opacity']])


## ##############################################

def CO_map(Inputs, fl):
	""" Make CO background map from Planck CO - type 2 map. HealPy is used
	simply to find the index corresponding to the location of interest.
	We use the CO(J=1->0) line (1st column of the HealPix data). """
	import healpy as hp
	co = pyfits.open(fl['COplanck'])
	conv = (np.pi/180.)
	head = pyfits.getheader('CMAP.fits')
	l = GenLon(Inputs)
	b = GenLat(Inputs)
	mapp = np.zeros( (int(l.size),int(b.size)) )
	error = np.zeros( (int(l.size),int(b.size)) )
	for i in np.arange( int(l.size) ):
		for j in np.arange( int(b.size) ):
			index = hp.ang2pix(2048,(np.pi/2.)-(conv*b[j]),(conv*l[i]),nest=True)
			mapp[i,j] = co[1].data[index][0]
			error[i,j] = co[1].data[index][1]

	radx = ((Inputs['nxpix']*Inputs['binsz'])/2.0) + (6./2.0)
	rady = ((Inputs['nypix']*Inputs['binsz'])/2.0) + (6./2.0)
	pxcen = radx/Inputs['binsz']
	pycen = rady/Inputs['binsz']
# ADD IN THIS CORRECTION FACTOR, FOUND IN http://arxiv.org/abs/1501.03606
# :: Divide by 1.16 to remove 13CO contamination
	maps = pyfits.PrimaryHDU( np.transpose(mapp)/1.16 )
	maps.header['CTYPE1'] = ('GLON-CAR')
	maps.header['CRPIX1'] = ( pxcen, 'Reference Pixel')
	maps.header['CRVAL1'] = (Inputs['glon'], 'Galactic longitude at reference pixel')
	maps.header['CDELT1'] = (-Inputs['binsz'],'x-axis increment per pixel')
	maps.header['CUNIT1'] = ('deg','Physical Units for x-axis')
	maps.header['CTYPE2'] = ('GLAT-CAR')
	maps.header['CRPIX2'] = ( pycen, 'Reference pixel')
	maps.header['CRVAL2'] = (Inputs['glat'],'Galactic latitude at reference pixel')
	maps.header['CDELT2'] = (Inputs['binsz'], 'y-axis increment per pixel')
	maps.header['CUNIT2'] = ('deg','Physical units for y-axis')
	maps.header['Equinox'] = (2000.,'Equinox of coordinates')
	maps.writeto(fl['COall'])

	

	err  = pyfits.PrimaryHDU( np.transpose(error) )
	err.header = maps.header
	err.writeto(fl['COer'])
	COThresh(Inputs,fl)
	CO_Cldkill(Inputs,fl)

	if ( (np.absolute(Inputs['glon']) < 6) and (mdd == 'all')):
		# Make flat template if necessary
		co = pyfits.open(fl['COall'])
		co[0].data = np.ones( co[0].data.shape )
		co[0].data /= (co[0].data.sum() * (Inputs['binsz']**2)*(np.pi/180.)**2)
		co.writeto('FlatTemp_norm.fits')

	print "Done CO maps"


## ##############################################
def COThresh(Inputs,fl):
	""" We open the CO map and it's error map from Planck and 
	set about thresholding and normalizing the map for use. 

	Right now, it is set to a 2-sigma threshold.
	"""
	co = pyfits.open(fl['COall'])
	err= pyfits.open(fl['COer'])
	co[0].data[ co[0].data < 2*err[0].data ] = 0.0
	co[0].data[ co[0].data < 0 ] = 0.0
	co.writeto(fl['COmpnon'])
	norm = co[0].data.sum() * (np.pi/180.)**2 * (Inputs['binsz']**2)
	co[0].data /= norm
	co.writeto(fl['COmp'])
	return "bwoop."

###########################################################################
def CO_Cldkill(Inputs,fl):
	""" Take CO map and remove cloud of interest (center of ROI).
	"""
	l = Inputs['glon']
	b = Inputs['glat']
	co = pyfits.open(fl['COmpnon'])
	coWCS = astWCS.WCS(co[0].header,mode="pyfits")
	co[0].data = np.transpose( co[0].data )
	ny = co[0].header['NAXIS2']
	nx = co[0].header['NAXIS1']
	comsk = np.zeros( (ny,nx) )
	comsk[ co[0].data > 0 ] = 1
	co[0].data *= comsk
	(maxx,maxy) = coWCS.wcs2pix(l,b)
	maxx = np.where( co[0].data == co[0].data.max() )[0][0]
	maxy = np.where( co[0].data == co[0].data.max() )[1][0]
	cldmsk = np.zeros( (ny,nx) )
	cldmsk[maxx,maxy] = 1
	while True:
		diff_temp = comsk - cldmsk
		cldmsk = ndimage.binary_dilation(cldmsk)
		cldmsk *= comsk
		diff = comsk - cldmsk
		if ( (diff - diff_temp).sum() == 0 ):
			break

	co[0].data -= cldmsk*co[0].data
	co[0].data = np.transpose(co[0].data)
	co.writeto(fl['COcldmsk'])
	return "Removed cloud from CO map."



## ##############################################
def Avmap_make(Inputs,fl):
	""" Make color excess map from SFD98. Using older, non-HealPix
	version. Download of SFD website. Use multiply the E(B-V)
	by 3.1 to get Av (assuming Galactic average value of Rv).

	01 Dec 2014: using Planck-derived E(B-V)
	"""
	mp = pyfits.open(fl['Dust_allsky'])

	l = GenLon(Inputs)
	b = GenLat(Inputs)

	Av = np.zeros( (l.size,b.size) )
	Av_= np.zeros( (l.size,b.size) )

	conv = (np.pi/180.)
	for i in np.arange( int(l.size) ):
		for j in np.arange( int(b.size) ):
			index = hp.ang2pix(2048,(np.pi/2.)-(conv*b[j]),(conv*l[i]),nest=True)
			Av[i,j]= 3.1*mp[1].data[index][0]	# tau353
			Av_[i,j]=3.1*mp[1].data[index][1]	# tau353 error

	radx = ((Inputs['nxpix']*Inputs['binsz'])/2.0) + (6./2.0)
	rady = ((Inputs['nypix']*Inputs['binsz'])/2.0) + (6./2.0)
	pxcen = radx/Inputs['binsz']
	pycen = rady/Inputs['binsz']
#   notice the 3.1 above: need Av, not E(B-V)
#   and the 1.49e4 is the recommended conversion factor from
#   tau353 --> E(B-V)_{tau353}
	maps = pyfits.PrimaryHDU( 1.49e4*np.transpose(Av) )
	maps.header['CTYPE1'] = ('GLON-CAR')
	maps.header['CRPIX1'] = ( pxcen, 'Reference Pixel')
	maps.header['CRVAL1'] = (Inputs['glon'], 'Galactic longitude at reference pixel')
	maps.header['CDELT1'] = (-Inputs['binsz'],'x-axis increment per pixel')
	maps.header['CUNIT1'] = ('deg','Physical Units for x-axis')
	maps.header['CTYPE2'] = ('GLAT-CAR')
	maps.header['CRPIX2'] = ( pycen, 'Reference pixel')
	maps.header['CRVAL2'] = (Inputs['glat'],'Galactic latitude at reference pixel')
	maps.header['CDELT2'] = (Inputs['binsz'], 'y-axis increment per pixel')
	maps.header['CUNIT2'] = ('deg','Physical units for y-axis')
	maps.header['Equinox'] = (2000.,'Equinox of coordinates')
	maps.writeto(fl['Dust'])

	mapr = pyfits.PrimaryHDU(1.49e4*np.transpose(Av_))
	mapr.header = maps.header
	mapr.writeto('../Dust_err.fits')
#	ebvv=pyfits.PrimaryHDU(np.transpose(ebv))
#	ebvv.header = maps.header
#	ebvv.writeto(fl['Dust'])


## ##############################################
def Dust_resid(Inputs,fl):
	""" Make the standard dark gas template, E(B-V)res -- the color
	excess which is not explained by HI nor CO. Find the least-
	squares fit of HI+CO to the Av map. Remove it. This (generally
	positive) residual supposedly accounts for all gas/dust not
	traced by HI or CO. """
	var = np.array([5e21,1e20])
	if (fl['opacity'].lower() == 'thick'):
		(R,q) = scipy.optimize.fmin(fun_thick,var)
	elif (fl['opacity'].lower() == 'thin'):
		(R,q) = scipy.optimize.fmin(fun_thin,var)

	h = open('Properties.txt','a')
	h.write(' R = %s \n q = %s \n'%(R,q) )
	h.close()
	dust = pyfits.open(fl['Dust'])
	hi1  = pyfits.open(fl['HImp_%snon'%fl['opacity']])
	hi2  = pyfits.open(fl['HImp_%snon_gt'%fl['opacity']])
	co   = pyfits.open(fl['COmpnon'])
	co_msk=pyfits.open(fl['COcldmsk'])
	Eres = dust[0].data - ( (1.0/(R)) * ((hi1[0].data + hi2[0].data) + 
				(q*co[0].data)) )
	Eres[ Eres < 0.0 ] = 0.0
	resmap = pyfits.PrimaryHDU(Eres)
	resmap.header = hi1[0].header
	resmap.writeto(fl['Ebv_%snon'%fl['opacity']])
	norm = Eres.sum() * (np.pi/180)**2 * (Inputs['binsz']**2)

	# Now for the Av,res map without the CO-emitting area
	co[0].data -= co_msk[0].data
	Eres[ co[0].data > 0 ] = 0
	resmap = pyfits.PrimaryHDU(Eres)
	resmap.header = hi1[0].header
	resmap.writeto(fl['Ebvcldmsk'])

	Eres /= norm
	resmap = pyfits.PrimaryHDU(Eres)
	resmap.header = hi1[0].header
	resmap.writeto(fl['Ebv_%s'%fl['opacity']])

	
	



## ##############################################
def fun_thin(var):
	""" Used to make the E(B-V)residual map: optimize this linear 
	combination of E(B-V), N(HI), and W_CO to find the color excess
	that can't be explained by HI or CO. """
	dust = pyfits.open('../Dust.fits')
	hi1  = pyfits.open('../HI_thin.fits')
	hi2  = pyfits.open('../HI_thin_gt.fits')
	co   = pyfits.open('../CO_temp.fits')
	Eres = dust[0].data - ( (1.0/(var[0])) * ((hi1[0].data + hi2[0].data) + 
				(var[1]*co[0].data)) )
	return (Eres**2).sum()

## ##############################################
def fun_thick(var):
	""" Used to make the E(B-V)residual map: optimize this linear 
	combination of E(B-V), N(HI), and W_CO to find the color excess
	that can't be explained by HI or CO. """
	dust = pyfits.open('../Dust.fits')
	hi1  = pyfits.open('../HI_thick.fits')
	hi2  = pyfits.open('../HI_thick_gt.fits')
	co   = pyfits.open('../CO_temp.fits')
	Eres = dust[0].data - ( (1.0/(var[0])) * ((hi1[0].data + hi2[0].data) + 
				(var[1]*co[0].data)) )
	return (Eres**2).sum()



## ##############################################
def djs_angpos(xval):
	""" taken from Schlegel dust map routines Put angles into 0 <= angle < 360"""
	if (xval > 0 or xval == 0):
		return ( xval%360 )
	else:
		return ( (xval%360)+360 )


## ##############################################
def Diffuse_std(Inputs,fl,dat,fix):
	""" Modify model file to include the standard diffuse sources,
	ones that will not change: HI, inverse Compton, and 
	isotropic emission. HI assumes a broken power law 2 
	spectrum. Unlike some other Fermi papers, now (Nov.20, 2014)
	the isotropic component is left free. """

	dat.append('<source name="HI" type="DiffuseSource">\n')
	dat.append('<!-- point source units are cm^-2 s^-1 MeV^-1 -->\n')
	dat.append('\t<spectrum type="BrokenPowerLaw2">\n')
	dat.append('\t\t<parameter free="1" max="1000.0" min="0.001" name="Integral" scale="1e-27" value="1.0"/>\n')
	dat.append('\t\t<parameter free="0" max="-1.0" min="-5.0" name="Index1" scale="1.0" value="-1.8"/>\n')
	dat.append('\t\t<parameter free="1" max="-1.0" min="-5.0" name="Index2" scale="1.0" value="-2.8"/>\n')
	dat.append('\t\t<parameter free="1" max="5000.0" min="500.0" name="BreakValue" scale="1.0" value="1200.0"/>\n')
	dat.append('\t\t<parameter free="0" max="1000000.0" min="20.0" name="LowerLimit" scale="1.0" value="%s"/>\n'% Inputs['emin'])
	dat.append('\t\t<parameter free="0" max="1000000.0" min="20.0" name="UpperLimit" scale="1.0" value="%s"/>\n'% Inputs['emax'])
	dat.append('\t</spectrum>\n')
	dat.append('\t<spatialModel file="%s" type="SpatialMap">\n'%(fl['HImp_%snon'%fl['opacity']]))
	dat.append('\t\t<parameter free="0" max="1e3" min="1e-3" name="Prefactor" scale="1.0" value="1.0"/>\n')
	dat.append('\t</spatialModel>\n')
	dat.append('</source>\n')

    #### HI far #### HI_far scale 
	dat.append('<source name="HI_far" type="DiffuseSource">\n')
	dat.append('\t<spectrum type="PowerLaw2">\n')
	dat.append('\t\t<parameter free="1" max="1e5" min="1e-6" name="Integral" scale="1e-27" value="1"/>\n')
	dat.append('\t\t<parameter free="1" max="1" min="-4" name="Index" scale="1.0" value="-2"/>\n')
	dat.append('\t\t<parameter free="0" max="1000000.0" min="20.0" name="LowerLimit" scale="1.0" value="%s"/>\n'% Inputs['emin'])
	dat.append('\t\t<parameter free="0" max="1000000.0" min="20.0" name="UpperLimit" scale="1.0" value="%s"/>\n'% Inputs['emax'])
	dat.append('\t</spectrum>\n')
	dat.append('\t<spatialModel file="%s" type="SpatialMap">\n'%(fl['HImp_%snon_gt'%fl['opacity']]))
	dat.append('\t\t<parameter free="0" max="1e3" min="1e-3" name="Prefactor" scale="1.0" value="1.0"/>\n')
	dat.append('\t</spatialModel>\n')
	dat.append('</source>\n')
    #### IC ####
	dat.append('<source name="IC" type="DiffuseSource">\n')
	dat.append('\t<spectrum type="ConstantValue">\n')
	dat.append('\t\t<parameter free="1" max="100.0" min="1e-6" name="Value" scale="1.0" value="1.0"/>\n')
	dat.append('\t</spectrum>\n')
	dat.append('\t<spatialModel file="%s" type="MapCubeFunction">\n'% fl['IC'])
	dat.append('\t\t<parameter free="0" max="1e3" min="1e-3" name="Normalization" scale="1.0" value="1.0"/>\n')
	dat.append('\t</spatialModel>\n')
	dat.append('</source>\n')
    #### iso ####
	dat.append('<source name="iso_P8R2_SOURCE_V6_v06" type="DiffuseSource">\n')
	dat.append('\t<spectrum file="%s" type="FileFunction">\n'% fl['iso'])
	dat.append('\t\t<parameter free="%s" max="10" min="1e-2" name="Normalization" scale="1" value="1"/>\n'%int(fix))
	dat.append('\t</spectrum>\n')
	dat.append('\t<spatialModel type="ConstantValue">\n')
	dat.append('\t\t<parameter free="0" max="10.0" min="0.0" name="Value" scale="1.0" value="1.0"/>\n')
	dat.append('\t</spatialModel>\n')
	dat.append('</source>\n')
	return dat


## ##############################################
def Diffuse_iso(Inputs,fl,dat,fix):
	""" Modify model file to include the isotropic emission. """
    #### iso ####
	dat.append('<source name="iso_P8R2_SOURCE_V6_v06" type="DiffuseSource">\n')
	dat.append('\t<spectrum file="%s" type="FileFunction">\n'% fl['iso'])
	dat.append('\t\t<parameter free="%s" max="10" min="1e-2" name="Normalization" scale="1" value="1"/>\n'%int(fix))
	dat.append('\t</spectrum>\n')
	dat.append('\t<spatialModel type="ConstantValue">\n')
	dat.append('\t\t<parameter free="0" max="10.0" min="0.0" name="Value" scale="1.0" value="1.0"/>\n')
	dat.append('\t</spatialModel>\n')
	dat.append('</source>\n')
	return dat


## ##############################################
def Diffuse_CO(Inputs,fl,dat,fix,cld):
	""" Modify model file to include diffuse emission from CO 

	Inputs: standard analysis inputs
	fl    : file names for analysis
	dat   : model xml string list """
# CO scale was 1e-8, now 1e-7
	dat.append('<source name="CO" type="DiffuseSource">\n')
	dat.append('\t<spectrum type="PowerLaw2">\n')
	dat.append('\t\t<parameter free="%s" max="1e5" min="1e-6" name="Integral" scale="1e-6" value="1"/>\n'%int(fix))
	dat.append('\t\t<parameter free="%s" max="1" min="-4" name="Index" scale="1.0" value="-2"/>\n'%int(fix))
	dat.append('\t\t<parameter free="0" max="1000000.0" min="20.0" name="LowerLimit" scale="1.0" value="%s"/>\n'% Inputs['emin'])
	dat.append('\t\t<parameter free="0" max="1000000.0" min="20.0" name="UpperLimit" scale="1.0" value="%s"/>\n'% Inputs['emax'])
	dat.append('\t</spectrum>\n')
	if ( cld == True ):
		dat.append('\t<spatialModel file="%s" type="SpatialMap">\n'%(fl['COmpnon']))
	elif ( cld == False ):
		dat.append('\t<spatialModel file="%s" type="SpatialMap">\n'%(fl['COmcldmsk']))

	dat.append('\t\t<parameter free="0" max="1e3" min="1e-3" name="Prefactor" scale="1.0" value="1.0"/>\n')
	dat.append('\t</spatialModel>\n')
	dat.append('</source>\n')

	return dat


## ##############################################
def Diffuse_DG(Inputs,fl,dat,fix,cld):
	""" Modify model file to include diffuse emission from dark gas 

	Inputs: standard analysis inputs
	fl    : file names for analysis
	dat   : model xml string list """
# DG scale was 1e-8, now 1e-5
	dat.append('<source name="E(B-V)" type="DiffuseSource">\n')
	dat.append('\t<spectrum type="PowerLaw2">\n')
	dat.append('\t\t<parameter free="%s" max="1e5" min="1e-6" name="Integral" scale="1e-5" value="1"/>\n'%int(fix))
	dat.append('\t\t<parameter free="%s" max="1" min="-4" name="Index" scale="1.0" value="-2"/>\n'%int(fix))
	dat.append('\t\t<parameter free="0" max="1000000.0" min="20.0" name="LowerLimit" scale="1.0" value="%s"/>\n'% Inputs['emin'])
	dat.append('\t\t<parameter free="0" max="1000000.0" min="20.0" name="UpperLimit" scale="1.0" value="%s"/>\n'% Inputs['emax'])
	dat.append('\t</spectrum>\n')
	if ( cld ):
		dat.append('\t<spatialModel file="%s" type="SpatialMap">\n'%(fl['Ebv_%snon'%fl['opacity']]))
	else:
		dat.append('\t<spatialModel file="%s" type="SpatialMap">\n'%(fl['Ebvcldmsk']))

	dat.append('\t\t<parameter free="0" max="1e3" min="1e-3" name="Prefactor" scale="1.0" value="1.0"/>\n')
	dat.append('\t</spatialModel>\n')
	dat.append('</source>\n')

	return dat

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
	dat.append('\t\t<parameter free="1" max="100000" min="1e-06" name="Integral" scale="1e-9" value="7" />\n')
	dat.append('\t\t<parameter free="1" max="1" min="-4" name="Index" scale="1" value="-2.2" />\n')
	dat.append('\t\t<parameter free="0" max="1000000" min="20" name="LowerLimit" scale="1" value="%s" />\n'%Inputs['emin'])
	dat.append('\t\t<parameter free="0" max="1000000" min="20" name="UpperLimit" scale="1" value="%s" />\n'%Inputs['emax'])
	dat.append('\t</spectrum>\n')
	dat.append('\t<spatialModel type="SkyDirFunction">\n')
	dat.append('\t\t<parameter free="0" max="360" min="-360" name="RA" scale="1" value="%s" />\n'%RA)
	dat.append('\t\t<parameter free="0" max="90" min="-90" name="DEC" scale="1" value="%s" />\n'%dec)
	dat.append('\t</spatialModel>\n')
	dat.append('</source>\n')

	return dat

## ##############################################
def Mdl_Make(Inputs,fl,directory):
	""" Make a model XML file from "make3FGLxml" and then modifying it to include
	the HI, CO, Av_residual maps as backgrounds 

	Input: input data, ROI position, energy range, etc.
	fl   : list of file names from SelectLoad.FileNames
	key  : all/NoCODG/Ptsrc -- make model with or without CO/dark gas template
	       or with a point source replacing the CO and dark gas templates."""
	import numpy as np
	if ( fl['mode'] == 'orig' ):
		return "No need for now\n\n\n bwoop."

	if ( fl['mode'].lower() == 'all' ):
		try:
			f = open(fl['MODEL_newsrc'],'r')
		except:
			f = open('MODEL.xml','r')
	else:
		try:
			f = open(fl['outmdl2_all'],'r')
		except:
			f = open(fl['outmdl1_all'],'r')


	dat = f.readlines()
	f.close()
#    i = 0
#    while (i < 19 ):
#	del dat[-1]
#	i += 1
	del dat[-1]
	for i in np.arange( len(dat) ):
		try:
			if ( ('gll_iem' in dat[i]) ):
				del dat[i:i+10]
				break
		except:
			pass

	for i in np.arange( len(dat) ):
		try:
			if ('CO' in dat[i]):
				COint = float( dat[i+2].split('value="')[1].split('"')[0] )
			elif ('E(B-V)' in dat[i]):
				EBVint= float( dat[i+2].split('value="')[1].split('"')[0] )
			elif ('HI' in dat[i]):
				HIint = float( dat[i+2].split('value="')[1].split('"')[0] )
		except:
			pass

	for i in np.arange( len(dat) ):
		try:
			if (dat[i].split('name="')[1].split('" ')[0] == 'iso_P8R2_SOURCE_V6_v06'):
				del dat[i:i+8]
				break
		except:
			pass




	ind = np.where( np.array(dat) == '<source name="SMC" type="DiffuseSource">\n' )
#    try:
#	del dat[ind[0][0]:ind[0][0]+12]
#    except:
#	pass
	try:
		ind[0][0]
		smcfl = '/home/abrahams/Tools-April2012/ScienceTools-v9r27p1-fssc-20120410-i686-pc-linux-gnu-libc2.5/i686-pc-linux-gnu-libc2.5/refdata/fermi/genericSources/Templates/SMC.fits'
		i = 0
		while ( ("SMC.fits" in dat[i]) == False ):
			i += 1

		dat[i] = '\t<spatialModel file="%s" type="SpatialMap">\n'%(smcfl)
	except:
		pass

# Now to the separate cases
	mdd = fl['mode'].lower()
	if ( mdd == 'all' ):
		dat = Diffuse_std(Inputs,fl,dat,1)
		dat = Diffuse_CO(Inputs,fl,dat,1,True)
		dat = Diffuse_DG(Inputs,fl,dat,1,True)
	elif ( mdd == "nocodg" ):
		os.mkdir( '%s'%mdd )
		dat = Diffuse_iso(Inputs,fl,dat,1)
		for i in np.arange( len(dat) ):
			if ('CO' in dat[i]):
				if ( 'type="PowerLaw2"' in dat[i+1] ):
					del dat[i:i+11]
					break
				elif ('type="BrokenPowerLaw2"' in dat[i+1]):
					del dat[i:i+13]
					break

		for i in np.arange( len(dat) ):
			if ('E(B-V)' in dat[i]):
				if ('type="PowerLaw2"' in dat[i+1]):
					del dat[i:i+11]
					break
				elif ('type="BrokenPowerLaw2"' in dat[i+1]):
					del dat[i:i+13]
					break

#	dat = Diffuse_std(Inputs,fl,dat,1)
	elif ( mdd == "ptsrc" ):
		os.mkdir( '%s'%mdd )
		(glon,glat) = COPtsrc_coords(Inputs,fl)
#	dat = Diffuse_std(Inputs,fl,dat,1)
		dat = Diffuse_ptsrc(Inputs,fl,dat,glon,glat,0)
		dat = Diffuse_iso(Inputs,fl,dat,1)
		for i in np.arange( len(dat) ):
			if ('CO' in dat[i]):
				if ( 'type="PowerLaw2"' in dat[i+1] ):
					del dat[i:i+11]
					break
				elif ('type="BrokenPowerLaw2"' in dat[i+1]):
					del dat[i:i+13]
					break
		for i in np.arange( len(dat) ):
			if ('E(B-V)' in dat[i]):
				if ( 'type="PowerLaw2"' in dat[i+1] ):
					del dat[i:i+11]
					break
				elif ('type="BrokenPowerLaw2"' in dat[i+1]):
					del dat[i:i+13]
					break

	elif( mdd == 'nodg' ):
		os.mkdir('%s'%mdd)
#	dat = Diffuse_std(Inputs,fl,dat,1)
#	dat = Diffuse_CO(Inputs,fl,dat,1)
		dat = Diffuse_iso(Inputs,fl,dat,1)
		for i in np.arange( len(dat) ):
			if ('E(B-V)' in dat[i]):
				if ( 'type="PowerLaw2"' in dat[i+1] ):
					del dat[i:i+11]
					break
				elif ('type="BrokenPowerLaw2"' in dat[i+1]):
					del dat[i:i+13]
					break

	elif( mdd == 'noco' ):
		os.mkdir('%s'%mdd)
		dat = Diffuse_iso(Inputs,fl,dat,1)
#	dat = Diffuse_std(Inputs,fl,dat,1)
#	dat = Diffuse_DG(Inputs,fl,dat,1)
		for i in np.arange( len(dat) ):
			if ('CO' in dat[i]):
				if ( 'type="PowerLaw2"' in dat[i+1] ):
					del dat[i:i+11]
					break
				elif ('type="BrokenPowerLaw2"' in dat[i+1]):
					del dat[i:i+13]
					break

	elif ( mdd == "all_plus" ):
		os.mkdir( '%s'%mdd )
		(glon,glat) = COPtsrc_coords(Inputs,fl)
#	dat = Diffuse_std(Inputs,fl,dat,1)
#	dat = Diffuse_CO(Inputs,fl,dat,1)
#	dat = Diffuse_DG(Inputs,fl,dat,1)
		dat = Diffuse_ptsrc(Inputs,fl,dat,glon,glat,0)
		dat = Diffuse_iso(Inputs,fl,dat,0)
	elif ( (mdd == 'nocld') ):
		os.mkdir( '%s'%mdd )
		for i in np.arange( len(dat) ):
			if ( 'name="CO"' in dat[i] ):
				while ('</source>' not in dat[i]):
					i += 1
					if ( 'type="SpatialMap"' in dat[i] ):
						dat[i] = '\t<spatialModel file="%s" type="SpatialMap">\n'%fl['COcldmsk']
			elif ( 'name="E(B-V)"' in dat[i] ):
				while ('</source>' not in dat[i]):
					i += 1
					if ( 'type="SpatialMap"' in dat[i] ):
						dat[i] = '\t<spatialModel file="%s" type="SpatialMap">\n'%fl['Ebvcldmsk']
		
	elif ( (mdd =='fx_src') ):
		os.mkdir('%s'%mdd)
		dat = Diffuse_iso(Inputs,fl,dat,0)
#	dat = Diffuse_std(Inputs,fl,dat,0)
#	dat = Diffuse_CO(Inputs,fl,dat,0)
#	dat = Diffuse_DG(Inputs,fl,dat,0)
#	os.system('cp mdl2_all.xml %s/%s'%(mdd,fl['MODEL_power_law']))
	elif ( mdd == 'power_law' ):
		os.mkdir('%s'%mdd)
#	os.system('cp mdl2_all.xml %s/%s'%(mdd,fl['MODEL_power_law']))
		return 0

	if ( (np.absolute(Inputs['glon']) < 6) and (mdd == 'all')):
		dat.append('<source name="Flat" type="DiffuseSource">\n')
		dat.append('\t<spectrum type="PowerLaw2">\n')
		dat.append('\t\t<parameter free="%s" max="1e5" min="1e-6" name="Integral" scale="1e-7" value="1"/>\n'%int(1))
		dat.append('\t\t<parameter free="%s" max="1" min="-4" name="Index" scale="1.0" value="-2"/>\n'%int(1))
		dat.append('\t\t<parameter free="0" max="200000.0" min="20.0" name="LowerLimit" scale="1.0" value="%s"/>\n'% Inputs['emin'])
		dat.append('\t\t<parameter free="0" max="200000.0" min="20.0" name="UpperLimit" scale="1.0" value="%s"/>\n'% Inputs['emax'])
		dat.append('\t</spectrum>\n')
		dat.append('\t<spatialModel file="%s" type="SpatialMap">\n'%('/home/abrahams/HICO_survey/SourceSearch/l%sb%s/FlatTemp_norm.fits'%(int(Inputs['glon']),int(Inputs['glat']))))
		dat.append('\t\t<parameter free="0" max="1e3" min="1e-3" name="Prefactor" scale="1.0" value="1.0"/>\n')
		dat.append('\t</spatialModel>\n')
		dat.append('</source>\n')

	dat.append('</source_library>')
	g = open(fl['MODEL_%s'%mdd],'w')
	g.writelines(dat)
	g.close()
	return 0


## ##############################################
def Mdl_Gas(Inputs,fl,MDL):
	""" Make a model file of just the individual gas components

	Input: input data, ROI position, energy range, etc.
	fl   : list of file names from SelectLoad.FileNames
	MDL  : tell me which component I am looking at"""

	if ( fl['mode'] == 'orig' ):
		return "No need for now\n\n\n bwoop."


	dat = ['<?xml version="1.0" standalone="no"?>']
	dat.append('<source_library title="source library">')

# Now to the separate cases: CO+DG, no CO+DG, no CO+DG with point source #

	if ( MDL.lower() == 'hi' ):
		dat.append('<source name="HI" type="DiffuseSource">\n')
		dat.append('\t<spectrum type="PowerLaw2">\n')
		dat.append('\t\t<parameter free="1" max="1e5" min="1e-6" name="Integral" scale="2e-7" value="1"/>\n')
		dat.append('\t\t<parameter free="1" max="1" min="-4" name="Index" scale="1.0" value="-2"/>\n')
		dat.append('\t\t<parameter free="0" max="200000.0" min="20.0" name="LowerLimit" scale="1.0" value="%s"/>\n'% Inputs['emin'])
		dat.append('\t\t<parameter free="0" max="200000.0" min="20.0" name="UpperLimit" scale="1.0" value="%s"/>\n'% Inputs['emax'])
		dat.append('\t</spectrum>\n')
		dat.append('\t<spatialModel file="%s" type="SpatialMap">\n'%(fl['HImp_%s'%fl['opacity']]))
		dat.append('\t\t<parameter free="0" max="1e3" min="1e-3" name="Prefactor" scale="1.0" value="1.0"/>\n')
		dat.append('\t</spatialModel>\n')
		dat.append('</source>\n')
	elif( MDL.lower() == 'hi_far' ):
		dat.append('<source name="HI_far" type="DiffuseSource">\n')
		dat.append('\t<spectrum type="PowerLaw2">\n')
		dat.append('\t\t<parameter free="1" max="1e5" min="1e-6" name="Integral" scale="2e-8" value="1"/>\n')
		dat.append('\t\t<parameter free="1" max="1" min="-4" name="Index" scale="1.0" value="-2"/>\n')
		dat.append('\t\t<parameter free="0" max="200000.0" min="20.0" name="LowerLimit" scale="1.0" value="%s"/>\n'% Inputs['emin'])
		dat.append('\t\t<parameter free="0" max="200000.0" min="20.0" name="UpperLimit" scale="1.0" value="%s"/>\n'% Inputs['emax'])
		dat.append('\t</spectrum>\n')
		dat.append('\t<spatialModel file="%s" type="SpatialMap">\n'%(fl['HImp_%s_gt'%fl['opacity']]))
		dat.append('\t\t<parameter free="0" max="1e3" min="1e-3" name="Prefactor" scale="1.0" value="1.0"/>\n')
		dat.append('\t</spatialModel>\n')
		dat.append('</source>\n')
	elif( MDL.lower() == 'co' ):
		dat = Diffuse_CO(Inputs,fl,dat)
	elif( MDL.lower() == 'dg' ):
		dat = Diffuse_DG(Inputs,fl,dat)
	
	dat.append('</source_library>')
	g = open(fl['MODEL_%s'%MDL.lower()],'w')
	g.writelines(dat)
	g.close()
	return 0




## ##############################################
def COPtsrc_coords(Inputs,fl):
	""" We want to replace the CO/DG templates with point sources.
	Wherever a large clump of CO is found, replace with a 
	point source. Return coordinates of location of CO max
	close to the center of the ROI.
	We smooth the CO map with a Gaussian filter to average over
	noise and then threshold. """

	im = pyfits.open(fl['COmp'])
	im[0].data = np.transpose( im[0].data )
	im2= ndimage.gaussian_filter(im[0].data, sigma=4)

	dist = np.zeros( (im[0].header['CRPIX2'],im[0].header['CRPIX1']) )
	for i in np.arange( im[0].header['CRPIX1']):
		np.sqrt( (i - im[0].header['CRPIX1'])**2 + 
			(np.arange(im[0].header['CRPIX2']) - im[0].header['CRPIX2'])**2 )

	im2[dist > 50] = 0	# Kill CO outside 5 degrees or so
	max_ = np.where( im2 == np.max(im2) )

	wcs = astWCS.WCS(im[0].header, mode="pyfits")
	(glon,glat) = wcs.pix2wcs(max_[0][0],max_[1][0])
	h = open('Properties.txt','a')
	try:
		h.write('\n Coordinates of CO max = (%s, %s)\n'%(glon,glat))
		h.close()
		return (glon,glat)
	except:
		h.write('\n Coordinates of CO max = (%s, %s)\n'%(Inputs['glon'],Inputs['glat']))
		h.close()
		return (Inputs['glon'],Inputs['glat'])






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
		Input['glon'],Input['glat'],2000)
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
	gaps.gtexpcube2['binsz'] = Input['binsz']
	gaps.gtexpcube2['ebinalg'] = 'LOG'
	gaps.gtexpcube2['emin'] = Input['emin']
	gaps.gtexpcube2['emax'] = Input['emax']
	gaps.gtexpcube2['enumbins'] = int(Input['enumbins'])
	if (Input['glat'] > 70):
		gaps.gtexpcube2['coordsys'] = 'CEL'
		gaps.gtexpcube2['xref'] = ra
		gaps.gtexpcube2['yref'] = de
		gaps.gtexpcube2['proj'] = 'AIT'
		gaps.gtexpcube2.run()
		return 0
	else:
		gaps.gtexpcube2['coordsys'] = 'GAL'
		gaps.gtexpcube2['xref'] = Input['glon']
		gaps.gtexpcube2['yref'] = Input['glat']
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

	gaps.srcMaps['scfile'] = fl['SC_fl']
	if (fl['time'] == False):
		gaps.srcMaps['expcube']= fl['ltcube']
	elif (fl['time'] == True):
		gaps.srcMaps['expcube']= fl['ltcube_t']
	else:
		print "Error in source map input."

	gaps.srcMaps['cmap'] = fl['CCUBE']
    # The following ensures I do not have to remake after finding new sources
	gaps.srcMaps['ptsrc'] = 'no'
	gaps.srcMaps['evtype']= 3
	gaps.srcMaps['bexpmap'] = fl['bexpmp']
	gaps.srcMaps['irfs'] = fl['IRFS']
#    if (fl['mode'] == 'all'):
#	# Find new sources, add to MODEL_all.xml
#	NewSrc_Detect(Input,fl)

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

#    nm2 = ( name[1:5] + " " + name[5:] )
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
		if ( TS > 500 ):		## Super strong sources are free
			like = ThawParams(Input,fl,like,cat,src,ind)
			wrt.write('%s kept, TS = %s\n'%(src,TS))
			pass
#	elif ( (center_dist > 7) ):
#	    like.deleteSource( src )
#	    wrt.write('%s killed: too far, dist = %s, TS = %s......\n'%(src,center_dist,TS))
		elif ( (center_dist > 4.5) ):
			like = FreezeSrc2(Input,fl,like,cat,src,ind)
			wrt.write('%s totally frozen......\n'%src)
		elif ( (center_dist > 2.5) and (TS < 500 ) ):
			like = FreezeShape(Input,fl,like,cat,src,ind)
			wrt.write('%s shape frozen......\n'%src)
		elif ( (center_dist < 2) and (TS < 500 ) ): 
			if ( (src[-1] != 'c') ):
				like.deleteSource( src )
				wrt.write('\t %s killed: too close, CO dist = %s, TS = %s.. confused source\n'%(src,CO_peak_dist,TS))
			elif ( (fl['mode'].lower() != 'all_plus') and (fl['mode'].lower() != 'ptsrc') ):
				like.deleteSource( src )
				wrt.write('\t %s killed: too close, CO dist = %s, TS = %s......\n'%(src,CO_peak_dist,TS))
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
		like.thaw(like.par_index(src,'Index'))
	elif (cat[1].data['SpectrumType'][ind] == 'PowerLaw2'):
		like.thaw(like.par_index(src,'Integral'))
		like.thaw(like.par_index(src,'Index'))
	elif (cat[1].data['SpectrumType'][ind] == 'LogParabola'):
		like.thaw(like.par_index(src,'norm'))
		like.thaw(like.par_index(src,'alpha'))
		like.freeze(like.par_index(src,'beta'))
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
		like.thaw(like.par_index(src,'Index'))
		like.thaw(like.par_index(src,'Ebreak'))
		like.thaw(like.par_index(src,'P1'))
	elif (cat[1].data['SpectrumType'][ind] == 'BPLExpCutoff'):
		like.thaw(like.par_index(src,'Prefactor'))
		like.thaw(like.par_index(src,'Index1'))
		like.thaw(like.par_index(src,'Ebreak'))
		like.thaw(like.par_index(src,'P1'))
		like.thaw(like.par_index(src,'Index2'))
		like.thaw(like.par_index(src,'Eabs'))
	elif (cat[1].data['SpectrumType'][ind] == 'PLSuperExpCutoff'):
		like.thaw(like.par_index(src,'Prefactor'))
		like.thaw(like.par_index(src,'Index1'))
		like.thaw(like.par_index(src,'Cutoff'))
		like.thaw(like.par_index(src,'Index2'))
	elif (cat[1].data['SpectrumType'][ind] == 'Gaussian'):
		like.thaw(like.par_index(src,'Prefactor'))
		like.thaw(like.par_index(src,'Mean'))
		like.thaw(like.par_index(src,'Sigma'))

	return like



## ##############################################
def DiffFreeThaw(Input,fl,like,src,act):

	if (act.lower() == 'freeze'):
		params = like.model[src].funcs['Spectrum'].paramNames
		num_params = len(params)
		for par in params:
			like.freeze(like.par_index(src,par))
		
	elif (act.lower() == 'thaw'):
		if (src == 'HI'):
			like.thaw(like.par_index(src,'Integral'))
			like.thaw(like.par_index(src,'Index1'))
			like.thaw(like.par_index(src,'Index2'))
			like.thaw(like.par_index(src,'BreakValue'))
		elif ( (src == 'HI_far') or (src == 'CO') or (src == 'E(B-V)') ):
			like.thaw(like.par_index(src,'Integral'))
			like.thaw(like.par_index(src,'Index'))
		elif ( ('iso' in src) ):
			like.thaw(like.par_index(src,'Normalization'))
		elif ( src == 'IC' ):
			like.thaw(like.par_index(src,'Value'))
		elif ( '3FGL' in src ):
			(ind, center_dist, TS) = FndSrcProp(Input,fl,nm,fl['cat_var'])
			like = ThawParams(Input,fl,like,fl['cat_var'],src,ind)


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

	g = open('LogLike.dat','a')
    # Here: freeze point sources
	name = like.sourceNames()
	for nm in name:
		if ('3FGL' in nm):
			(ind, center_dist, TS) = FndSrcProp(Input,fl,nm,fl['cat_var'])
			like = FreezeSrc2(Input,fl,like,fl['cat_var'],nm,ind)
			print ("Freeeeeeeze %s"%nm)
		elif ( (nm in fl['gal']) or (nm in fl['iso']) ):
			like.deleteSource(nm)
		elif ('Add' in nm):
			like = FreezeSrc2(Input,fl,like,fl['cat_var'],nm,ind)

# The like.optimizer='Minuit' because: if parameter at limit, it stays there
# at that causes the whole thing to DIE
	try:
#		if (like.optimizer == 'DRMNFB'):
#			lglike = like.fit()
 #		elif (like.optimizer == 'Minuit'):
#			like.writeXml('MDL.xml')
#			like = BinnedAnalysis(obs,srcModel='MDL.xml',optimizer='Minuit')
#			likeobj = pyLike.Minuit(like.logLike)
#			lglike = like.fit(covar=True,optObject=likeobj)
##			g.write('Fit quality %s\t'%like1obj.getQuality())
		likeobj = pyLike.Minuit(like.logLike)
		lglike = like.fit(covar=True,optObject=likeobj)
		g.write('Fit quality %s\t'%like1obj.getQuality())
		print ('Fit quality %s\t'%like1obj.getQuality())
	except:
		like.optimizer='minuit'
		for nm in name:
			if ((nm == 'CO') or (nm == 'E(B-V)') or (nm == 'HI_far') ):
				if ( like.Ts(nm) < 3 ):
					like.deleteSource(nm)
			elif ( nm == 'HI' ):
				like.model['HI'].funcs['Spectrum'].params['Index1'].setValue(-1.8)
				like.model['HI'].funcs['Spectrum'].params['Integral'].setValue(2)

		likeobj = pyLike.Minuit(like.logLike)
		lglike = like.fit(covar=True,optObject=likeobj)

	# Delete weak diffuse sources and freeze strong diffuse sources and thaw good point sources
	name = like.sourceNames()
	for nm in name:
		if ( ('3FGL' not in nm) and ('Add' not in nm) ):
			if ( like.Ts(nm) < 3 ):
				like.deleteSource(nm)
			else:
				like = DiffFreeThaw(Input,fl,like,nm,'freeze')
				print ("Frozen %s!"%nm)
#		elif ( '3FGL' in nm ):				### DO NOT UNFREEZE PT SOURCES
#			like = ThawSrc(Input,fl,like,nm)
#		elif ( (nm == 'CO') or (nm == 'E(B-V)') or (nm == 'HI_far') ):
#			if ( like.Ts(nm) < 3 ):
#				like.deleteSource(nm)
#
#	print ("Frozen Sources !! :). Fit point sources.\n")
#	# freeze some sources: fit only important ones
#	like = MdlPrep(Input,fl,like,likerun)
#	try:
#		if (like.optimizer.lower() == 'drmnfb'):
#			lglike = like.fit()
 #		elif (like.optimizer.lower() == 'minuit'):
#			likeobj2 = pyLike.Minuit(like.logLike)
#			lglike2 = like.fit(covar=True,optObject=likeobj2)
#			print "Fit Quality = %s\n"%likeobj2.getQuality()
#	except:
#		print 'poop'
#		return ("poop")

	like.writeXml('mdl1_all.xml')
	print "sources fit, finding new sources.\n"
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

	# Need to freeze point sources and unfreeze diffuse sources...
	name = like.sourceNames()
	for nm in name:
		if ( ('3FGL' not in nm) and ('Add' not in nm) ):
			like = DiffFreeThaw(Input,fl,like,nm,'thaw')
		elif ( ('3FGL' in nm) ):
			like = DiffFreeThaw(Input,fl,like,nm,'freeze')
		elif ( 'Add' in nm ):
			like.freeze(like.par_index(nm,'Index'))


#	return (lglike2,like)
	return (lglike,like)


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

	if ( fl['time'] == False):
		obs = BinnedObs(binnedExpMap=fl['bexpmp'],
			srcMaps=srcmap,expCube=fl['ltcube'],irfs=fl['IRFS'])
	elif ( fl['time'] == True):
		obs = BinnedObs(binnedExpMap=fl['bexpmp'],
			srcMaps=srcmap,expCube=fl['ltcube_t'],irfs=fl['IRFS'])
	else:
		print "Error in binned analysis input."

#	like = BinnedAnalysis(obs, srcModel=srcmdl,optimizer='DRMNFB')
	like = BinnedAnalysis(obs, srcModel=srcmdl,optimizer='Minuit')
	likeobj = pyLike.Minuit(like.logLike)

	g = open('LogLike.dat','a')
	print "setting up likelihood analysis."
    # 1st likelihood analysis.
 #   try:
#	optim = 1
#	if (fl['mode'].lower() == 'all'):
#	    (lglike,like) = PrepAll(Input,fl,like,likerun)
#	else:
#	    lglike = like.fit()
#    except:
#	like.writeXml('mdl_int.xml')
#	g.write('Error in DRMNFB, going to Minuit.\n')
#	like.optimizer = 'Minuit'
#	like = BinnedAnalysis(obs, srcModel=srcmdl, optimizer='Minuit')
#	like1obj = pyLike.Minuit(like.logLike)
#	if ( fl['mode'].lower() == 'all' ):
#	    (lglike,like) = PrepAll(Input,fl,like,likerun)
#	    g.write('Fit quality %s\t'%like1obj.getQuality())
#	else:
#	    lglike = like.fit(covar=True,optObject=like1obj)

#	os.remove('mdl_int.xml')
#	optim = 0
	
	optim = 1
	try:
		if (fl['mode'].lower() == 'all'):
			print "Fitting first round."
			(lglike,like) = PrepAll(Input,fl,like,likerun)
			print "Done first round."
		else:
			lglike = like.fit(covar=True,optObject=likeobj)
	except:
		return "poops. Can't finish fitting for some reason. Check the 'PrepAll' routine, probably.\n\n\n Done \n\n\n :("

	print "Writing likelihoods.\n"
	g.write('-log(likelihood) of %s = %s\n'%(fl['mode'],lglike))
	print "Writing results\n."
	WriteResults(Input,fl,like,optim)
	like = WeakSrc(Input,fl,like) # Get rid of weak sources, including diffuse
	like.writeXml(outmdl1)
	like = BinnedAnalysis(obs,srcModel=outmdl1,optimizer='NewMinuit')
#	like.optimizer = 'NewMinuit'
	like2obj = pyLike.NewMinuit(like.logLike)
	MdlMap(Input,fl,1,0)

    # 2nd likelihood anlaysis to fine tune. 
	like.tol=1e-4
	like.setEnergyRange(emin=Input['emin'],emax=Input['emax'])
	like = MdlPrep(Input,fl,like,2)
	try:
		like.thaw(like.par_index('CO','Integral'))
	except:
		pass

	lging = open('Logging.txt','a')
	lging.write('%s'%like.model)
	lging.close()
	try:
		lglike2 = like.fit(covar=True,optObject=like2obj)
		g.write('likelihood return code = %s\t'%like2obj.getRetCode())
		g.write('Final -log(likelihood) of %s = %s\n'%(fl['mode'],lglike2))
		like.writeXml(outmdl2)
		MdlMap(Input,fl,2,0)
		error = 0
	except:
		g.write('ERROR IN NEWMINUIT FIT\n\n')
		error = 1
		pass


	if ( (like2obj.getRetCode() != 0) ):
		like.freeze(like.par_index('HI','Index1'))
		like.freeze(like.par_index('HI','Index2'))
		like.freeze(like.par_index('HI','BreakValue'))
		try:
			like.freeze(like.par_index('HI_far','Index'))
		except:
			pass 

		try:
			like.freeze(like.par_index('CO','Index'))
		except:
			pass

		try:
			like.freeze(like.par_index('E(B-V)','Index'))
		except:
			pass

		try:
			lglike2 = like.fit(covar=True,optObject=like2obj)
			g.write('likelihood return code = %s\t'%like2obj.getRetCode())
			g.write('Final -log(likelihood) of %s = %s\n'%(fl['mode'],lglike2))
			os.system('mv %s mdl2_int.xml'%outmdl2)
			like.writeXml(outmdl2)
			MdlMap(Input,fl,2,0)
			error = 0
		except:
			g.write('ERROR IN NEWMINUIT Re-FIT\n\n')
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
			elif ( (src != 'SMC') or (src != 'LMC') or (src != " ") ):
				try:
					like = FreezeSrc(Input,fl,like,fl['cat_var'],src)
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
			elif ( (src != 'SMC') or (src != 'LMC') or (src != " ") or (src[-1] != 'c')):
				try:
					(ind, center_dist, TS) = FndSrcProp(Input,fl,src,fl['cat_var'])
					like = FreezeSrc2(Input,fl,like,fl['cat_var'],src,ind)
					wrt.write('%s frozen.\n'%src)
				except:
					wrt.write('%s NOT FROZEN. ERROR!!!\n'%src)
					pass

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

	mdls = ['HI','HI_far','CO','E(B-V)','IC','iso_P8R2_SOURCE_V6_v06.txt']
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
			try:
				if (like.Ts( name ) < 20):
					ulim = UpperLim(Input,like, name)
					if ((ulim != 0) and (error != 1)):
						lim = float(str(ulim).split(' ph')[0])
						mm.write(' Upper Limit = %s ph/cm^2/s \n'%lim)
					elif ( (ulim == 0) and (error != 1) ):
						lim = 3*like.flux(name,emin=Input['emin'],emax=Input['emax'])
						limer=3*like.fluxError(name,emin=Input['emin'],emax=Input['emax'])
						mm.write(' Upper Limit about (3*flux) %s +/- %s\n'%(lim,limer))
					elif ( (ulim == 0) and (error == 1) ):
						lim = 3*like.flux(name,emin=Input['emin'],emax=Input['emax'])
						mm.write(' Upper Limit about (3*flux) %s \n'%lim)
			except:
				mm.write("Error in writing upper limit for source %s\n"%name)

			try:
				if ( error != 1 ):
					mm.write('energy flux = %s +/- %s \n'%(like.energyFlux(name,
						emin=Input['emin'],emax=Input['emax']),
						like.energyFluxError('HI',emin=Input['emin'],
						emax=Input['emax'])))
					mm.write('flux = %s +/- %s'%(like.flux(name,
						emin=Input['emin'],emax=Input['emax']),
						like.fluxError('HI',emin=Input['emin'],
						emax=Input['emax'])))
				else:
					mm.write('energy flux = %s \n'%(like.energyFlux(name,
						emin=Input['emin'],emax=Input['emax'])))
					mm.write('flux = %s'%(like.flux(name,
						emin=Input['emin'],emax=Input['emax'])))
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
	like.setEnergyRange(emin=Input['emin'],emax=Input['emax'])
	g = open('LogLike.dat','a')
	g.write('\n\n')
	for src in like.sourceNames():
		if ('3FGL' in src):
			try:
				(ind,rad,TS) = FndSrcProp(Input,fl,src,fl['cat_var'])
				if (rad < Input['rad']+2):
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
		like.setEnergyRange(emin=Input['emin'],emax=Input['emax'])
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
		if (TS_ < 1):
			if ((name == 'CO') or (name == 'HI') or (name=='E(B-V)') or (name == 'HI_far') or (name == 'IC')):
				like.deleteSource( name )
				frz.write('During %s fitting, %s killed with TS of %s\n'%(fl['mode'],name,TS_))
				mm.write('During %s fitting, %s killed with TS of %s\n\n'%(fl['mode'],name,TS_))
		elif (TS_ < 50):
			if ((name != 'CO') and (name != 'HI') and (name != 'E(B-V)') and (name != 'HI_far') and (name != 'IC') ):
				like.deleteSource( name )
				frz.write('During %s fitting, %s killed with TS of %s\n'%(fl['mode'],name,TS_))
#		elif (TS_ < 50):
#			if ((name != 'CO') and (name != 'HI') and (name != 'E(B-V)') and (name != 'HI_far') and (name != 'IC') ):
#				try:
#					(ind, center_dist, TS) = FndSrcProp(Input,fl,src,fl['cat_var'])
#					FreezeShape(Input,fl,like,cat,src,ind)
#					frz.write('During %s fitting, %s shape frozen with TS = %s'%(fl['mode'],name,TS_))
#				except:
#					frz.write("During %s fitting, %s NOT killed due to catalog error\n"%(fl['mode'],name))
		else:
			frz.write('During %s fitting, %s has TS of %s\n'%(fl['mode'],name,TS_))
	

	frz.close()
	mm.close()
	return like


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
			ul[nm].compute(emin=Input['emin'],emax=Input['emax'])
		except:
			return 0


	return ul[nm].results[0]



