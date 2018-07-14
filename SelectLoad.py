#!/home/abrahams/epd_free-7.3-2-rh5-x86/bin/python

import sys, os
import gt_apps as gaps
import numpy as np
from astLib import astCoords
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
		inp['glon'] = [el for el in arg if el[0] == 'glon'][0]
		inp['glat'] = [el for el in arg if el[0] == 'glat'][0]
		inp['rad']  = [el for el in arg if el[0] == 'rad'][0]
	except:
		sys.exit("Invalid data selection parameters")

	try:
		inp['nxpix'] = [el for el in arg if el[0] == 'nxpix'][0]
		inp['nypix'] = [el for el in arg if el[0] == 'nypix'][0]
		inp['binsz'] = [el for el in arg if el[0] == 'binsz'][0]
		inp['enumbins'] = [el for el in arg if el[0] == 'enumbins'][0]
	except:
		sys.exit("Invalid binning parameters")

	try:
		inp['emin'] = [el for el in arg if el[0] == 'emin'][0]
	except:
		inp['emin'] = ['emin','250']

	try:
		inp['emax'] = [el for el in arg if el[0] == 'emax'][0]
	except:
		inp['emax'] = ['emax','30000']

	try:
		inp['zcut'] = [el for el in arg if el[0] == 'zcut'][0]
	except:
		inp['zcut'] = ['zcut','True']

	try:
		inp['tmin'] = [el for el in arg if el[0] == 'tmin'][0]
	except:
		inp['tmin'] = ['tmin','239557417']

	try:
		inp['tmax'] = [el for el in arg if el[0] == 'tmax'][0]
	except:
		inp['tmax'] = ['tmax','392574813']

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

	inp['glon'] = [el for el in args if el[0] == 'glon'][0]
	inp['glat'] = [el for el in args if el[0] == 'glat'][0]
	inp['rad']  = [el for el in args if el[0] == 'rad'][0]
	inp['emin'] = [el for el in args if el[0] == 'emin'][0]
	inp['emax'] = [el for el in args if el[0] == 'emax'][0]
	inp['zcut'] = [el for el in args if el[0] == 'zcut'][0]
	inp['tmin'] = [el for el in args if el[0] == 'tmin'][0]
	inp['tmax'] = [el for el in args if el[0] == 'tmax'][0]
	inp['nxpix'] = [el for el in args if el[0] == 'nxpix'][0]
	inp['nypix'] = [el for el in args if el[0] == 'nypix'][0]
	inp['binsz'] = [el for el in args if el[0] == 'binsz'][0]
	inp['enumbins'] = [el for el in args if el[0] == 'enumbins'][0]
	return ValConvert(inp)

## ##############################################
## Convert the inputted strings to numbers/boolean
## as necessary for use later in the analysis.
def ValConvert(inp):
	""" Converts the inputted strings into floating point numbers. """
	inp['glon'][1] = float(inp['glon'][1])
	inp['glat'][1] = float(inp['glat'][1])
	inp['rad'][1] = float(inp['rad'][1])
	inp['emin'][1] = float(inp['emin'][1])
	inp['emax'][1] = float(inp['emax'][1])
	inp['tmin'][1] = float(inp['tmin'][1])
	inp['tmax'][1] = float(inp['tmax'][1])
	inp['zcut'][1] = str2bool(inp['zcut'][1])
	inp['nxpix'][1] = float(inp['nxpix'][1])
	inp['nypix'][1] = float(inp['nypix'][1])
	inp['binsz'][1] = float(inp['binsz'][1])
	inp['enumbins'][1] = float(inp['enumbins'][1])
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
	fl   : list of file nafl['srcmp_nocodg']= 'NoCODG/srcmap.fits'
	fl['outmdl1_nocodg']='NoCODG/mdl1_nocodg.xml'
	fl['outmdl2_nocodg']='NoCODG/mdl2_nocodg.xml'
	fl['outmap1_nocodg']='NoCODG/mdlmap1_nocodg.fits'
	fl['outmap2_nocodg']='NoCODG/mdlmap2_nocodg.fits'mes from SelectLoad.FileNames"""
	(ra,de) = astCoords.convertCoords("GALACTIC","J2000",
		Inputs['glon'][1],Inputs['glat'][1],2000)
#	gaps.filter['evclass']=2		# P7REP and earlier
	gaps.filter['evclass'] = 128 # source = 128, ultracleanveto = 1024
	gaps.filter['evtype'] = 3
	gaps.filter['ra']=ra
	gaps.filter['dec']=de
	gaps.filter['rad']=Inputs['rad'][1]
	gaps.filter['zmax']=90
	gaps.filter['tmin']=Inputs['tmin'][1]
	gaps.filter['tmax']=Inputs['tmax'][1]
	gaps.filter['infile']=fl['files']
	gaps.filter['outfile']=fl['sel_fl']
	gaps.filter['emin']=Inputs['emin'][1]
	gaps.filter['emax']=Inputs['emax'][1]

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
	ndir = pre + 'l%sb%s/ebin250_10000/'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
	fl['files'] = pre + 'all_sky_files.txt'
	fl['mkt_fl']='mkt.fits'
	fl['SC_fl']= pre + 'AllSky_p8/SC.fits'
	fl['gas_spec'] = 'spec'
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
	fl['ltcube']= pre + ('Analysis/ltCube_%s' % (int(Inputs['rad'][1]))) + 'deg.fits'
	fl['ltcube_t']= 'ltCube.fits'
	fl['bexpmp']='bexpmap.fits'
	fl['MODEL_newsrc'] = 'MODEL_newsrc.xml'
	fl['srcmp_int1'] = 'srcmap_int1.fits'
	fl['MODEL_int1'] = 'mdl1_int1.xml'
    # Models, source maps, and output models...
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
	fl['srcmp_nocodg']= 'nocodg/srcmap_nocodg.fits'
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

  	fl['COall'] = pre + 'l%sb%s/CO_temp_all.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
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
	fl['HImp_thick'] = pre + 'l%sb%s/HI_thick_norm.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
	fl['HImp_thicknon'] = pre + 'l%sb%s/HI_thick.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
	fl['HImp_thin'] = pre + 'l%sb%s/HI_thin_norm.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
	fl['HImp_thinnon'] = pre + 'l%sb%s/HI_thin.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
	fl['HImp_thick_gt'] = pre + 'l%sb%s/HI_thick_gt_norm.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
	fl['HImp_thicknon_gt'] = pre + 'l%sb%s/HI_thick_gt.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
	fl['HImp_thin_gt'] = pre + 'l%sb%s/HI_thin_gt_norm.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
	fl['HImp_thinnon_gt'] = pre + 'l%sb%s/HI_thin_gt.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))

	fl['COmp'] = pre + 'l%sb%s/CO_temp_norm.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
	fl['COmpnon'] = pre + 'l%sb%s/CO_temp.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
	fl['COer'] = pre + 'l%sb%s/CO_temp_error.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
	fl['COcldmsk'] = pre+'l%sb%s/CO_temp_nocld.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
	fl['Ebv_thick']  = pre + 'l%sb%s/Ebv_thick_norm.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
	fl['Ebvcldmsk'] = pre + 'l%sb%s/Ebv_thick_nocld.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
	fl['Ebv_thicknon']  = pre + 'l%sb%s/Ebv_thick.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
	fl['Ebv_thin']  = pre + 'l%sb%s/Ebv_thin_norm.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
	fl['Ebv_thinnon']  = pre + 'l%sb%s/Ebv_thin.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
#	Old IC (from 77varhX... whatever)
	fl['IC']   = '/home/abrahams/GALPROP/results_54_04d90003/ics_isotropic_mapcube_54_04d90003'
#	New IC (parameters taken from FermiPlanck Chamaeleon paper.
#	fl['IC']   = '/home/abrahams/GALPROP/Chamrun_new/ics_isotropic_mapcube_54_04d90003'

	fl['Dust_allsky'] = pre + 'DustModel_Planck.fits'
	fl['Dust'] = pre + 'l%sb%s/Dust.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))
#    fl['iso']  = os.environ['FERMI_DIR']+ '/refdata/fermi/galdiffuse/' + 'iso_source_v05.txt'
#    fl['gal']  = os.environ['FERMI_DIR']+ '/refdata/fermi/galdiffuse/' + 'gll_iem_v05_rev1.fit'
	fl['iso'] = os.environ['FERMI_DIR'] + '/refdata/fermi/galdiffuse/' + 'iso_P8R2_SOURCE_V6_v06.txt'
	fl['gal'] = os.environ['FERMI_DIR'] + '/refdata/fermi/galdiffuse/' + 'gll_iem_v06.fits'
#   fl['ndust']= pre + 'SFD_dust_4096_ngp.fits'
#    fl['sdust']= pre + 'SFD_dust_4096_sgp.fits'
	fl['COplanck']= os.environ['FERMI_DIR']+ '/refdata/fermi/galdiffuse/' + 'CO_PlanckDR2_Type2.fits'
	fl['COplanck_type3'] = pre + 'CO_planck_type3.fits'
	fl['cat']  = os.environ['FERMI_DIR']+ '/refdata/fermi/' + 'gll_psc_v14.fit'

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
		Inputs['glon'][1],Inputs['glat'][1],2000)
	gaps.evtbin['axisrot'] = 0.0

	# Change other things in case we choose lightcurve, counts map,
	#        or counts cube
	if (fl['binalg'].lower() == 'cmap'):
		if (fl['CMAP'] in os.listdir('.')):
			return "CMAP made, continue"
		else:
			gaps.evtbin['nxpix'] = int(Inputs['nxpix'][1])
			gaps.evtbin['nypix'] = int(Inputs['nypix'][1])
			gaps.evtbin['binsz'] = Inputs['binsz'][1]
			gaps.evtbin['outfile']=fl['CMAP']
			if (Inputs['glat'][1] > 70):
				gaps.evtbin['coordsys'] = 'CEL'
				gaps.evtbin['xref'] = ra
				gaps.evtbin['yref'] = de
				gaps.evtbin['proj'] = 'AIT'
				gaps.evtbin.run()
			else:
				gaps.evtbin['coordsys'] = 'GAL'
				gaps.evtbin['xref'] = Inputs['glon'][1]
				gaps.evtbin['yref'] = Inputs['glat'][1]
				gaps.evtbin['proj'] = 'CAR'
				gaps.evtbin.run()

			gaps.evtbin['nxpix'] = int(Inputs['nxpix'][1] + (10./Inputs['binsz'][1]))
			gaps.evtbin['nypix'] = int(Inputs['nypix'][1] + (10./Inputs['binsz'][1]))
			gaps.evtbin['outfile']=fl['CMAP_big']
			if (Inputs['glat'][1] > 70):
				gaps.evtbin['coordsys'] = 'CEL'
				gaps.evtbin['xref'] = ra
				gaps.evtbin['yref'] = de
				gaps.evtbin['proj'] = 'AIT'
				gaps.evtbin.run()
			else:
				gaps.evtbin['coordsys'] = 'GAL'
				gaps.evtbin['xref'] = Inputs['glon'][1]
				gaps.evtbin['yref'] = Inputs['glat'][1]
				gaps.evtbin['proj'] = 'CAR'
				gaps.evtbin.run()

			return 0
	elif (fl['binalg'].lower() == 'lc'):
		gaps.evtbin['outfile'] = fl['LC']
		gaps.evtbin['tbinalg'] = 'LIN'
		gaps.evtbin['tstart'] = Inputs['tmin'][1]
		gaps.evtbin['tstop'] = Inputs['tmax'][1]
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
			gaps.evtbin['nxpix'] = int(Inputs['nxpix'][1])
			gaps.evtbin['nypix'] = int(Inputs['nypix'][1])
			gaps.evtbin['binsz'] = Inputs['binsz'][1]
			gaps.evtbin['ebinalg'] = 'LOG'
			gaps.evtbin['emin'] = Inputs['emin'][1]
			gaps.evtbin['emax'] = Inputs['emax'][1]
			gaps.evtbin['enumbins'] = int(Inputs['enumbins'][1])
			if (Inputs['glat'][1] > 70):
				gaps.evtbin['coordsys'] = 'CEL'
				gaps.evtbin['xref'] = ra
				gaps.evtbin['yref'] = de
				gaps.evtbin['proj'] = 'AIT'
				gaps.evtbin.run()
				return 0
			else:
				gaps.evtbin['coordsys'] = 'GAL'
				gaps.evtbin['xref'] = Inputs['glon'][1]
				gaps.evtbin['yref'] = Inputs['glat'][1]
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

	di= '/home/abrahams/HICO_survey/SourceSearch/l%sb%s'%( int(Inputs['glon'][1]),int(Inputs['glat'][1]) )
	if ( (fl['ltcube'] not in os.listdir(di)) and (fl['ltcube_t'] not in os.listdir(di)) ):
		gaps.expCube.run()
		return 0
	else:
		return "Livetime cube already present."


## ##############################################
## ##############################################
#
#if ((len(sys.argv) == 1)):
#    Input = ReadFile('Config.txt')
#elif ((sys.argv[1].lower() == 'config.txt')):
#    Input = ReadFile('Config.txt')
#else:
#    Input = ReadParams(sys.argv)
#
#print Input
### When this stuff is imported...
# The way to use this:
# ins = ['working_dir','glon=10','glat=20','rad=10',...]
# 
# Flname = FileNames(Input)
# Select(Input,Flname)
# mktime(Input,Flname)
# 
# Flname['binalg'] = 'CMAP'
# need to make a counts map
#
# Flname['binalg'] = 'LC'
# need to make light curve
#
# flname['binalg'] = 'CCUBE'
# need to make counts cube
