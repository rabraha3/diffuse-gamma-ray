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
#
# 5August2014: Fix fx_src for the lightcurve!!! The model is not created:
#	       To fix, just copy mdl2_all.xml and in LikeRun freeze the
#	       sources. Easy enough.
#
#
#
#
# Must be have changed into analysis directory!
#
# Here we make the model. We need to sift through the catalog and
# choose strong sources that are close to the center of the RoI.
# Then, we need to determine the spectral shape before writing then
# xml file. This routine also needs to create the spatial templates
# for HI, CO (from Planck), and dust.
#
# For many of these functions, we'll need to import a bunch of things
# that are present in the "SelectLoad.py" section.
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

	# previously else: lon = nl - (l-180)*2
	return lon

## ##############################################
def GenLon(Input):
	""" Longitude array for background maps """
	rad= ((Input['nxpix'][1]*Input['binsz'][1])/2.0) + (6./2.0)
	lu = (Input['glon'][1] + rad)%360
	ld = (Input['glon'][1] - rad)%360
	df = int(lu/Input['binsz'][1])

	sz = int(2*rad/Input['binsz'][1])
	if ( (lu - ld) < 0 ):
		df = (lu - 0.)/Input['binsz'][1]
		l = np.append( np.linspace(lu,0.,df),np.linspace((360-Input['binsz'][1]),ld,int(sz-df)))
	else:
		l = np.linspace(lu,ld,int(2*rad/Input['binsz'][1]))

	return l

## ##############################################
def GenLat(Input):
	""" Latitude array for background maps """
	rad= ((Input['nypix'][1]*Input['binsz'][1])/2.0) + (6./2.0)
	bu = (Input['glat'][1] + rad)
	bd = (Input['glat'][1] - rad)

#	return np.linspace(bd,bu,int(Input['nxpix'][1]+ (10./Input['binsz'][1])))
	return np.linspace(bd,bu,int(2*rad/Input['binsz'][1]))

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

	radx = ((Inputs['nxpix'][1]*Inputs['binsz'][1])/2.0) + (6./2.0)
	rady = ((Inputs['nypix'][1]*Inputs['binsz'][1])/2.0) + (6./2.0)
	pxcen = radx/Inputs['binsz'][1]
	pycen = rady/Inputs['binsz'][1]
	ima = pyfits.PrimaryHDU(np.transpose(image))
	ima.header['CTYPE1'] = ('GLON-CAR')
	ima.header['CRPIX1'] = ( pxcen, 'Reference Pixel')
	ima.header['CRVAL1'] = (Inputs['glon'][1], 'Galactic longitude at reference pixel')
	ima.header['CDELT1'] = (-Inputs['binsz'][1],'x-axis increment per pixel')
	ima.header['CUNIT1'] = ('deg','Physical Units for x-axis')
	ima.header['CTYPE2'] = ('GLAT-CAR')
	ima.header['CRPIX2'] = ( pycen, 'Reference pixel')
	ima.header['CRVAL2'] = (Inputs['glat'][1],'Galactic latitude at reference pixel')
	ima.header['CDELT2'] = (Inputs['binsz'][1], 'y-axis increment per pixel')
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
	norm = hi[0].data.sum() * (np.pi/180.)**2 * (Inputs['binsz'][1]**2)
	hi[0].data /= norm
	if ( dist == 1 ):
		hi.writeto(fl['HImp_%s'%fl['opacity']])
	elif ( dist == 2 ):
		hi.writeto(fl['HImp_%s_gt'%fl['opacity']])

	hi = None
	ima= None
	image = None
	cmap_w = None
	head = None
	del hi, ima, image, cmap_w, head


## ##############################################

def CO_map(Inputs, fl):
	""" Make CO background map from Planck CO - type 2 map. HealPy is used
	simply to find the index corresponding to the location of interest.
	We use the CO(J=1->0) line (1st column of the HealPix data). """
	import healpy as hp
	mmm = open('COmake.txt','a')
	mmm.write('Starting CO\n')
	co = pyfits.open(fl['COplanck'])
	conv = (np.pi/180.)
	head = pyfits.getheader('CMAP.fits')
	l = GenLon(Inputs)
	b = GenLat(Inputs)
	mmm.write("Opened and generated lat/lon\n")
	mapp = np.zeros( (int(l.size),int(b.size)) )
	error = np.zeros( (int(l.size),int(b.size)) )
	mmm.write("Starting map create\n")
	for i in np.arange( int(l.size) ):
		for j in np.arange( int(b.size) ):
			index = hp.ang2pix(2048,(np.pi/2.)-(conv*b[j]),(conv*l[i]),nest=True)
			mapp[i,j] = co[1].data[index][0]
			error[i,j] = co[1].data[index][1]

	mmm.write('Finished map\n')
	co = None
	del co
	radx = ((Inputs['nxpix'][1]*Inputs['binsz'][1])/2.0) + (6./2.0)
	rady = ((Inputs['nypix'][1]*Inputs['binsz'][1])/2.0) + (6./2.0)
	pxcen = radx/Inputs['binsz'][1]
	pycen = rady/Inputs['binsz'][1]
# ADD IN THIS CORRECTION FACTOR, FOUND IN http://arxiv.org/abs/1501.03606
# :: Divide by 1.16 to remove 13CO contamination
	mmm.write('Writing CO(J=1-0) map\n')
	maps = pyfits.PrimaryHDU( np.transpose(mapp)/1.16 )
	maps.header['CTYPE1'] = ('GLON-CAR')
	maps.header['CRPIX1'] = ( pxcen, 'Reference Pixel')
	maps.header['CRVAL1'] = (Inputs['glon'][1], 'Galactic longitude at reference pixel')
	maps.header['CDELT1'] = (-Inputs['binsz'][1],'x-axis increment per pixel')
	maps.header['CUNIT1'] = ('deg','Physical Units for x-axis')
	maps.header['CTYPE2'] = ('GLAT-CAR')
	maps.header['CRPIX2'] = ( pycen, 'Reference pixel')
	maps.header['CRVAL2'] = (Inputs['glat'][1],'Galactic latitude at reference pixel')
	maps.header['CDELT2'] = (Inputs['binsz'][1], 'y-axis increment per pixel')
	maps.header['CUNIT2'] = ('deg','Physical units for y-axis')
	maps.header['Equinox'] = (2000.,'Equinox of coordinates')
	maps.writeto(fl['COall'])

	
	mmm.write('Writing error map\n')
	err  = pyfits.PrimaryHDU( np.transpose(error) )
	err.header = maps.header
	err.writeto(fl['COer'])
	COThresh(Inputs,fl)
	CO_Cldkill(Inputs,fl)

	if ( (np.absolute(Inputs['glon'][1]) < 15) ):
		# Make flat template if necessary
		co = pyfits.open(fl['COall'])
		co[0].data = np.ones( co[0].data.shape )
		co[0].data /= (co[0].data.sum() * (Inputs['binsz'][1]**2)*(np.pi/180.)**2)
		co.writeto('../FlatTemp_norm.fits')

	mapp = None
	error = None
	co = None
	maps = None
	head = None
	del mapp, error, co, maps, head

	mmm.write('Done CO maps\n')
	mmm.close()
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
	norm = co[0].data.sum() * (np.pi/180.)**2 * (Inputs['binsz'][1]**2)
	co[0].data /= norm
	co.writeto(fl['COmp'])
	co = None
	err = None
	norm = None
	del co, err, norm
	return "bwoop."

###########################################################################
def CO_Cldkill(Inputs,fl):
	""" Take CO map and remove cloud of interest (center of ROI).
	"""
	l = Inputs['glon'][1]
	b = Inputs['glat'][1]
	co = pyfits.open(fl['COmpnon'])
	xcen = co[0].header['CRPIX1']
	ycen = co[0].header['CRPIX2']
	co[0].data = np.transpose( co[0].data )
	ny = co[0].header['NAXIS2']
	nx = co[0].header['NAXIS1']
	comsk = np.zeros( (ny,nx) )
	comsk[ co[0].data > 0 ] = 1
	co[0].data *= comsk
	cldmsk = np.zeros( (ny,nx) )
	cldmsk[xcen,ycen] = 1
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
	
	co = None
	comsk = None
	cldmsk = None
	diff = None
	diff_temp = None
	del co, comsk, cldmsk, diff, diff_temp
	return "Removed cloud from CO map."


####################################################################


#def fun_ebv(var):
#    ebv = pyfits.open('../EBV_pl.fits')
#    tau = pyfits.open('../Dust_.fits')
#    return ((ebv[0].data[ebv[0].data<0.2] - var*tau[0].data[ebv[0].data<0.2])**2).sum()

## ##############################################
def Avmap_make(Inputs,fl):
	""" Make color excess map from SFD98. Using older, non-HealPix
	version. Download of SFD website. Use multiply the E(B-V)
	by 3.1 to get Av (assuming Galactic average value of Rv).

	01 Dec 2014: using Planck-derived E(B-V)
	"""
#    Used for SFD map
#	if (Inputs['glat'][1] > 0):
#		mp = pyfits.open(fl['ndust'])
#	else:
#	mp = pyfits.open(fl['sdust'])
	mp = pyfits.open(fl['Dust_allsky'])

	l = GenLon(Inputs)
	b = GenLat(Inputs)

	Av = np.zeros( (l.size,b.size) )
	Av_= np.zeros( (l.size,b.size) )

	conv = (np.pi/180.)
#	ebv = np.zeros( (l.size,b.size) )
	for i in np.arange( int(l.size) ):
		for j in np.arange( int(b.size) ):
			index = hp.ang2pix(2048,(np.pi/2.)-(conv*b[j]),(conv*l[i]),nest=True)
			Av[i,j]= 3.1*mp[1].data[index][0]	# tau353
			Av_[i,j]=3.1*mp[1].data[index][1]	# tau353 error
#			ebv[i,j]=3.1*mp[1].data[index][2]	# E(B-V)

##			Av[i,j] = 3.55*nicest[1].data[int(index)/1024][0][int(index)%1024]
##			For the strange ordering . . . 

# The above is A_V(tau353)
# Ok, now, least squares fit: tau353 to E(B-V): fit to low color excess?? Haven't done yet
# as of September 1, 2015

	radx = ((Inputs['nxpix'][1]*Inputs['binsz'][1])/2.0) + (6./2.0)
	rady = ((Inputs['nypix'][1]*Inputs['binsz'][1])/2.0) + (6./2.0)
	pxcen = radx/Inputs['binsz'][1]
	pycen = rady/Inputs['binsz'][1]
#   notice the 3.1 when writing the fits file: need Av, not E(B-V)
#	maps = pyfits.PrimaryHDU( 3.1*np.transpose(Av) )
	maps = pyfits.PrimaryHDU( 1.49e4*np.transpose(Av) )
	maps.header['CTYPE1'] = ('GLON-CAR')
	maps.header['CRPIX1'] = ( pxcen, 'Reference Pixel')
	maps.header['CRVAL1'] = (Inputs['glon'][1], 'Galactic longitude at reference pixel')
	maps.header['CDELT1'] = (-Inputs['binsz'][1],'x-axis increment per pixel')
	maps.header['CUNIT1'] = ('deg','Physical Units for x-axis')
	maps.header['CTYPE2'] = ('GLAT-CAR')
	maps.header['CRPIX2'] = ( pycen, 'Reference pixel')
	maps.header['CRVAL2'] = (Inputs['glat'][1],'Galactic latitude at reference pixel')
	maps.header['CDELT2'] = (Inputs['binsz'][1], 'y-axis increment per pixel')
	maps.header['CUNIT2'] = ('deg','Physical units for y-axis')
	maps.header['Equinox'] = (2000.,'Equinox of coordinates')
	maps.writeto(fl['Dust'])
#
	mapr = pyfits.PrimaryHDU(1.49e4*np.transpose(Av_))
	mapr.header = maps.header
	mapr.writeto('../Dust_err.fits')
	mp = None
	Av = None
	Av_= None
	maps = None
	del mp, Av, Av_, maps


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
	norm = Eres.sum() * (np.pi/180)**2 * (Inputs['binsz'][1]**2)

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

	dust = None
	hi1 = None
	hi2 = None
	co = None
	co_msk = None
	Eres = None
	resmap = None
	del dust, hi1, hi2, co, co_msk, Eres, resmap
	



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

	dust = None
	hi1 = None
	hi2 = None
	co = None
	del dust, hi1, hi2, co
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
	dust = None
	hi1 = None
	hi2 = None
	co = None
	del dust, hi1, hi2, co
	
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
    # HI scale was 1e-8, now 1e-27
#    dat.append('<source name="HI" type="DiffuseSource">\n')
#    dat.append('\t<spectrum type="PowerLaw2">\n')
#    dat.append('\t\t<parameter free="%s" max="1e5" min="1e-6" name="Integral" scale="1e-27" value="1"/>\n'%int(fix))
#    dat.append('\t\t<parameter free="%s" max="1" min="-4" name="Index" scale="1.0" value="-2"/>\n'%int(fix))
#    dat.append('\t\t<parameter free="0" max="200000.0" min="20.0" name="LowerLimit" scale="1.0" value="%s"/>\n'% Inputs['emin'][1])
#    dat.append('\t\t<parameter free="0" max="200000.0" min="20.0" name="UpperLimit" scale="1.0" value="%s"/>\n'% Inputs['emax'][1])
#    dat.append('\t</spectrum>\n')
#    dat.append('\t<spatialModel file="%s" type="SpatialMap">\n'%(fl['HImp_%snon'%fl['opacity']]))
#    dat.append('\t\t<parameter free="0" max="1e3" min="1e-3" name="Prefactor" scale="1.0" value="1.0"/>\n')
#    dat.append('\t</spatialModel>\n')
#    dat.append('</source>\n')
	dat.append('<source name="HI" type="DiffuseSource">\n')
	dat.append('<!-- point source units are cm^-2 s^-1 MeV^-1 -->\n')
	if ( fl['gas_spec'].lower() == 'BPL2' ):
		dat.append('\t<spectrum type="BrokenPowerLaw2">\n')
		dat.append('\t\t<parameter free="1" max="1000.0" min="0.001" name="Integral" scale="1e-27" value="1.0"/>\n')
		dat.append('\t\t<parameter free="0" max="-1.0" min="-5.0" name="Index1" scale="1.0" value="-1.8"/>\n')
		dat.append('\t\t<parameter free="1" max="-1.0" min="-5.0" name="Index2" scale="1.0" value="-2.8"/>\n')
		dat.append('\t\t<parameter free="1" max="5000.0" min="500.0" name="BreakValue" scale="1.0" value="1200.0"/>\n')
		dat.append('\t\t<parameter free="0" max="1000000.0" min="20.0" name="LowerLimit" scale="1.0" value="%s"/>\n'% Inputs['emin'][1])
		dat.append('\t\t<parameter free="0" max="1000000.0" min="20.0" name="UpperLimit" scale="1.0" value="%s"/>\n'% Inputs['emax'][1])
	elif ( fl['gas_spec'].lower() == 'spec' ):
		dat.append('\t<spectrum file="/home/abrahams/HICO_survey/SourceSearch/Analysis/LIS_cas.txt" type="FileFunction">\n')
		dat.append('\t\t<parameter free="1" max="10" min="0.01" name="Normalization" scale="1" value="1" />\n')

	dat.append('\t</spectrum>\n')
	dat.append('\t<spatialModel file="%s" type="SpatialMap">\n'%(fl['HImp_%snon'%fl['opacity']]))
	dat.append('\t\t<parameter free="0" max="1e3" min="1e-3" name="Prefactor" scale="1.0" value="1.0"/>\n')
	dat.append('\t</spatialModel>\n')
	dat.append('</source>\n')

    #### HI far #### HI_far scale 
	dat.append('<source name="HI_far" type="DiffuseSource">\n')
	if ( fl['gas_spec'].lower() == 'BPL2' ):
		dat.append('\t<spectrum type="PowerLaw2">\n')
		dat.append('\t\t<parameter free="1" max="1e5" min="1e-6" name="Integral" scale="1e-27" value="1"/>\n')
		dat.append('\t\t<parameter free="1" max="1" min="-4" name="Index" scale="1.0" value="-2"/>\n')
		dat.append('\t\t<parameter free="0" max="1000000.0" min="20.0" name="LowerLimit" scale="1.0" value="%s"/>\n'% Inputs['emin'][1])
		dat.append('\t\t<parameter free="0" max="1000000.0" min="20.0" name="UpperLimit" scale="1.0" value="%s"/>\n'% Inputs['emax'][1])
	elif ( fl['gas_spec'].lower() == 'spec' ):
		dat.append('\t<spectrum file="/home/abrahams/HICO_survey/SourceSearch/Analysis/LIS_cas.txt" type="FileFunction">\n')
		dat.append('\t\t<parameter free="1" max="10" min="0.01" name="Normalization" scale="1" value="1" />\n')

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
	if ( fl['gas_spec'].lower() == 'spec' ):
		dat.append('\t<spectrum file="/home/abrahams/HICO_survey/SourceSearch/Analysis/qCO_cas.txt" type="FileFunction">\n')
		dat.append('\t\t<parameter free="1" max="10" min="0.01" name="Normalization" scale="1" value="1" />\n')
		dat.append('\t</spectrum>\n')
	else:
		dat.append('\t<spectrum type="PowerLaw2">\n')
		dat.append('\t\t<parameter free="%s" max="1e5" min="1e-6" name="Integral" scale="1e-6" value="1"/>\n'%int(fix))
		dat.append('\t\t<parameter free="%s" max="1" min="-4" name="Index" scale="1.0" value="-2"/>\n'%int(fix))
		dat.append('\t\t<parameter free="0" max="1000000.0" min="20.0" name="LowerLimit" scale="1.0" value="%s"/>\n'% Inputs['emin'][1])
		dat.append('\t\t<parameter free="0" max="1000000.0" min="20.0" name="UpperLimit" scale="1.0" value="%s"/>\n'% Inputs['emax'][1])
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
	if ( fl['gas_spec'].lower() == 'spec' ):
		dat.append('\t<spectrum file="/home/abrahams/HICO_survey/SourceSearch/Analysis/qdg_cas.txt" type="FileFunction">\n')
		dat.append('\t\t<parameter free="1" max="10" min="0.01" name="Normalization" scale="1" value="1" />\n')
		dat.append('\t</spectrum>\n')
	else:
		dat.append('\t<spectrum type="PowerLaw2">\n')
		dat.append('\t\t<parameter free="%s" max="1e5" min="1e-6" name="Integral" scale="1e-5" value="1"/>\n'%int(fix))
		dat.append('\t\t<parameter free="%s" max="1" min="-4" name="Index" scale="1.0" value="-2"/>\n'%int(fix))
		dat.append('\t\t<parameter free="0" max="1000000.0" min="20.0" name="LowerLimit" scale="1.0" value="%s"/>\n'% Inputs['emin'][1])
		dat.append('\t\t<parameter free="0" max="1000000.0" min="20.0" name="UpperLimit" scale="1.0" value="%s"/>\n'% Inputs['emax'][1])
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
	dat.append('\t\t<parameter free="0" max="1000000" min="20" name="LowerLimit" scale="1" value="%s" />\n'%Inputs['emin'][1])
	dat.append('\t\t<parameter free="0" max="1000000" min="20" name="UpperLimit" scale="1" value="%s" />\n'%Inputs['emax'][1])
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
	import os
	mdd = fl['mode'].lower()
	if ( fl['MODEL_%s'%mdd] in os.listdir('.') ):
		print "Model already made. what?"
		return "Model already made. what?"

	if ( fl['mode'] == 'orig' ):
		return "No need for now\n\n\n bwoop."

	if ( fl['mode'].lower() == 'all' ):
		try:
			f = open('MODEL.xml','r')
		except:
			print "No model!"
	else:
		try:
			f = open('/home/abrahams/HICO_survey/SourceSearch/l%sb%s/ebin250_10000/'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))+fl['outmdl2_all'],'r')
			print "Yes NewMinuit fit."
			print "In %s"%os.listdir('.')
		except:
			try:
				f = open(fl['outmdl1_all'],'r')
				print "Pooping, no %s"%fl['outmdl1_all']
			except:
				f = open('MODEL.xml','r')
				print "Nothing fit, nothing to see here."


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
			if ( 'iso_P8R2' in dat[i] ):
				del dat[i:i+8]
				break
		except:
			pass


	# Added Dec.2, 2015:
	# Give proper LMC template
#	ind = np.where( np.array(dat) == '<source name="LMC" type="DiffuseSource">\n' )
	smcdir = '/home/abrahams/Tools-Old/Tools-April2012/ScienceTools-v9r27p1-fssc-20120410-i686-pc-linux-gnu-libc2.5/i686-pc-linux-gnu-libc2.5/refdata/fermi/genericSources/Templates/'
	for i in range( len(dat) ):
		if ( 'LMC.fits' in dat[i] ):
			smcnm = smcdir + 'LMC.fits'
			dat[i] = '\t<spatialModel file="%s" type="SpatialMap">\n'%(smcnm)
		elif ( 'SMC.fits' in dat[i] ):
			smcnm = smcdir + 'SMC.fits'
			dat[i] = '\t<spatialModel file="%s" type="SpatialMap">\n'%(smcnm)

# Now to the separate cases
	if ( mdd == 'all' ):
		dat = Diffuse_std(Inputs,fl,dat,1)
		dat = Diffuse_CO(Inputs,fl,dat,1,True)
		dat = Diffuse_DG(Inputs,fl,dat,1,True)
	elif ( mdd == "nocodg" ):
		try:
			os.mkdir( '%s'%mdd )
		except:
			pass

		dat = Diffuse_iso(Inputs,fl,dat,1)
		for i in np.arange( len(dat) ):
			if ('name="CO"' in dat[i]):
				if ( 'type="PowerLaw2"' in dat[i+1] ):
					del dat[i:i+11]
					break
				elif ('type="BrokenPowerLaw2"' in dat[i+1]):
					del dat[i:i+13]
					break
				elif ( 'type="FileFunction"' in dat[i+1]):
					del dat[i:i+8]
					break

		for i in np.arange( len(dat) ):
			if ('E(B-V)' in dat[i]):
				if ('type="PowerLaw2"' in dat[i+1]):
					del dat[i:i+11]
					break
				elif ('type="BrokenPowerLaw2"' in dat[i+1]):
					del dat[i:i+13]
					break
				elif ( 'type="FileFunction"' in dat[i+1]):
					del dat[i:i+8]
					break

#	dat = Diffuse_std(Inputs,fl,dat,1)
	elif ( mdd == "ptsrc" ):
		try:
			os.mkdir( '%s'%mdd )
		except:
			pass

#		(glon,glat) = COPtsrc_coords(Inputs,fl)
#	dat = Diffuse_std(Inputs,fl,dat,1)
		glon = Inputs['glon'][1]
		glat = Inputs['glat'][1]
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
				elif ( 'type="FileFunction"' in dat[i+1]):
					del dat[i:i+8]
					break
		for i in np.arange( len(dat) ):
			if ('E(B-V)' in dat[i]):
				if ( 'type="PowerLaw2"' in dat[i+1] ):
					del dat[i:i+11]
					break
				elif ('type="BrokenPowerLaw2"' in dat[i+1]):
					del dat[i:i+13]
					break
				elif ( 'type="FileFunction"' in dat[i+1]):
					del dat[i:i+8]
					break

	elif( mdd == 'nodg' ):
		try:
			os.mkdir( '%s'%mdd )
		except:
			pass

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
				elif ( 'type="FileFunction"' in dat[i+1]):
					del dat[i:i+8]
					break

	elif( mdd == 'noco' ):
		try:
			os.mkdir( '%s'%mdd )
		except:
			pass

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
				elif ( 'type="FileFunction"' in dat[i+1]):
					del dat[i:i+8]
					break

	elif ( mdd == "all_plus" ):
		try:
			os.mkdir( '%s'%mdd )
		except:
			pass

#		(glon,glat) = COPtsrc_coords(Inputs,fl)
#	dat = Diffuse_std(Inputs,fl,dat,1)
#	dat = Diffuse_CO(Inputs,fl,dat,1)
#	dat = Diffuse_DG(Inputs,fl,dat,1)
		glon = Inputs['glon'][1]
		glat = Inputs['glat'][1]
#		dat = Diffuse_std(Inputs,fl,dat,1)
#		dat = Diffuse_CO(Inputs,fl,dat,1,True)
#		dat = Diffuse_DG(Inputs,fl,dat,1,True)
		dat = Diffuse_ptsrc(Inputs,fl,dat,glon,glat,0)
		dat = Diffuse_iso(Inputs,fl,dat,0)
	elif ( (mdd == 'nocld') ):
		## Was failing because "directory already exists"
		try:
			os.mkdir( '%s'%mdd )
		except:
			pass

		dat = Diffuse_iso(Inputs,fl,dat,1)
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
		try:
			os.mkdir( '%s'%mdd )
		except:
			pass

		dat = Diffuse_iso(Inputs,fl,dat,0)
#	dat = Diffuse_std(Inputs,fl,dat,0)
#	dat = Diffuse_CO(Inputs,fl,dat,0)
#	dat = Diffuse_DG(Inputs,fl,dat,0)
#	os.system('cp mdl2_all.xml %s/%s'%(mdd,fl['MODEL_power_law']))
	elif ( mdd == 'power_law' ):
		os.mkdir('%s'%mdd)
#	os.system('cp mdl2_all.xml %s/%s'%(mdd,fl['MODEL_power_law']))
		return 0

	if ( (np.absolute(Inputs['glon'][1]) < 15) and (mdd == 'all')):
		dat.append('<source name="Flat" type="DiffuseSource">\n')
		dat.append('\t<spectrum file="/home/abrahams/HICO_survey/SourceSearch/Analysis/fermi_bubbles_spectra.txt" type="FileFunction">\n')
		dat.append('\t\t<parameter free="1" max="10" min="0.01" name="Normalization" scale="1" value="1" />\n')
		dat.append('\t</spectrum>\n')


		dat.append('\t<spectrum type="PowerLaw2">\n')
		dat.append('\t\t<parameter free="%s" max="1e5" min="1e-6" name="Integral" scale="1e-7" value="1"/>\n'%int(1))
		dat.append('\t\t<parameter free="%s" max="1" min="-4" name="Index" scale="1.0" value="-2"/>\n'%int(1))
		dat.append('\t\t<parameter free="0" max="200000.0" min="20.0" name="LowerLimit" scale="1.0" value="%s"/>\n'% Inputs['emin'][1])
		dat.append('\t\t<parameter free="0" max="200000.0" min="20.0" name="UpperLimit" scale="1.0" value="%s"/>\n'% Inputs['emax'][1])
		dat.append('\t</spectrum>\n')
		dat.append('\t<spatialModel file="%s" type="SpatialMap">\n'%('/home/abrahams/HICO_survey/SourceSearch/l%sb%s/FlatTemp_norm.fits'%(int(Inputs['glon'][1]),int(Inputs['glat'][1]))))
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
		dat.append('\t\t<parameter free="0" max="200000.0" min="20.0" name="LowerLimit" scale="1.0" value="%s"/>\n'% Inputs['emin'][1])
		dat.append('\t\t<parameter free="0" max="200000.0" min="20.0" name="UpperLimit" scale="1.0" value="%s"/>\n'% Inputs['emax'][1])
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
		dat.append('\t\t<parameter free="0" max="200000.0" min="20.0" name="LowerLimit" scale="1.0" value="%s"/>\n'% Inputs['emin'][1])
		dat.append('\t\t<parameter free="0" max="200000.0" min="20.0" name="UpperLimit" scale="1.0" value="%s"/>\n'% Inputs['emax'][1])
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
		h.write('\n Coordinates of CO max = (%s, %s)\n'%(Inputs['glon'][1],Inputs['glat'][1]))
		h.close()
		return (Inputs['glon'][1],Inputs['glat'][1])



