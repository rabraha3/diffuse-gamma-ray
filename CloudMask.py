
import pyfits
from astLib import astWCS
from scipy import ndimage

f = open('CloudProperties_2sig.txt','r')
dat = f.readlines()
f.close()

## Need (l,b) of the ROI center ##

co = pyfits.open('CO_temp.fits')
coWCS = astWCS.WCS(co[0].header,mode='pyfits')
co[0].data = np.transpose( co[0].data )
l = 159
b = -34
(x,y) = coWCS.wcs2pix(l,b)
co[0].data = np.transpose( co[0].data )


ny = co[0].header['NAXIS2']
nx = co[0].header['NAXIS1']
comsk = np.zeros( (ny,nx) )
comsk[ co[0].data > 0 ] = 1
co[0].data *= comsk

cldmsk = np.zeros( (ny,nx) )
cldmsk[int(x),int(y)] = 1
while True:
	diff_temp = comsk - cldmsk
	cldmsk = ndimage.binary_dilation(cldmsk)
	cldmsk *= comsk
	diff = comsk - cldmsk
	if ( (diff - diff_temp).sum() == 0 ):
		break

co[0].data -= cldmsk*co[0].data
co.writeto('CO_temp_nocld.fits')
