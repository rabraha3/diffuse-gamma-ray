
import numpy as np
import pyfits,os
from astLib import astWCS
from scipy import ndimage

def CO_Cldkill():
	""" Take CO map and remove cloud of interest (center of ROI).
	"""
	co = pyfits.open('CO_temp.fits')
	coWCS = astWCS.WCS(co[0].header,mode="pyfits")
	co[0].data = np.transpose( co[0].data )
	ny = co[0].header['NAXIS2']
	nx = co[0].header['NAXIS1']
	comsk = np.zeros( (ny,nx) )
	comsk[ co[0].data > 0 ] = 1
	co[0].data *= comsk
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
	co.writeto('CO_temp_nocld.fits')




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


(glon,glat) = CoordFind()
pre = '/home/abrahams/HICO_survey/SourceSearch/'

for i in np.arange( len(glon) ):
	ndir = pre +(('l%sb%s')%(int(glon[i]),int(glat[i])))
	os.chdir(ndir)

	co = pyfits.open('CO_temp.fits')
	try:
		co_msk = pyfits.open('CO_temp_nocld.fits')
	except:
		CO_Cldkill()
		co_msk = pyfits.open('CO_temp_nocld.fits')

	Ebv = pyfits.open('Ebv_thick.fits')
	try:
		m = pyfits.open('Ebv_thick_nocld.fits')
		continue
	except:
		pass

	if ( co[0].data.shape != co_msk[0].data.shape ):
		co_msk[0].data = np.transpose(co_msk[0].data)

	co[0].data -= co_msk[0].data
	Ebv[0].data[ co[0].data > 0 ] = 0
	Ebv.writeto('Ebv_thick_nocld.fits')
	
	

print "Done!"



