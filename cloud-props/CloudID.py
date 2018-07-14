
import pyfits, os
import numpy as np
from scipy import ndimage
from astLib import astWCS
import sys

def CldResults(cldmsk,hi,co,dg,cldnum,coWCS,f):
    """ Take a cloud (defined by CO). Write down some properties of
	the cloud. """

    cocld = np.zeros( (co[0].header['NAXIS1'],co[0].header['NAXIS2']) )
    dgcld = np.zeros( (co[0].header['NAXIS1'],co[0].header['NAXIS2']) )
    hicld = np.zeros( (co[0].header['NAXIS1'],co[0].header['NAXIS2']) )
    cocld[cldmsk > 0] = co[0].data[ cldmsk > 0]
    dgcld[cldmsk > 0] = dg[0].data[ cldmsk > 0]
    hicld[cldmsk > 0] = hi[0].data[ cldmsk > 0]

    cm_x = 0
    cm_y = 0
    nx = co[0].header['NAXIS1']
    for i in np.arange( cldmsk.size ):
	if (cldmsk[i%nx,i/nx] > 0):
	    cm_x += i%nx
	    cm_y += i/nx

    cm_x /= cldmsk.sum()    
    cm_y /= cldmsk.sum()
    (glon,glat) = coWCS.pix2wcs(cm_y,cm_x)
    f.write("Cloud %s:\n"%cldnum)
    f.write("Pixel coordinates (%s, %s) = (l,b):(%s, %s)\n"%(cm_y,cm_x,glon,glat))
    f.write("Area = %s\n"%cldmsk.sum())
    f.write("\tWCO peak(K km s-1) = %s\n"%cocld[cocld!=0].max() )
    f.write("\tWCO minimum (K km s-1) = %s\n"%cocld[cocld!=0].min() )
    f.write("\tWCO mean (K km s-1) = %s\n"%cocld[cocld!=0].mean() )
    f.write("\tWCO standard deviation (K km s-1) = %s\n"%cocld[cocld!=0].std() )
    f.write("\tWCO sum (K km s-1) = %s\n"%cocld[cocld!=0].sum() )
    f.write("\n")
    f.write("\tDG peak (mag) = %s\n"%dgcld[dgcld!=0].max() )
    f.write("\tDG minimum (mag) = %s\n"%dgcld[dgcld!=0].min() )
    f.write("\tDG mean (mag) = %s\n"%dgcld[dgcld!=0].mean() )
    f.write("\tDG standard deviation (mag) = %s\n"%dgcld[dgcld!=0].std() )
    f.write("\tDG sum (mag) = %s\n"%dgcld[dgcld!=0].sum() )
    f.write("\n")
    f.write("\tN(HI) peak (cm-2) = %s\n"%hicld[hicld!=0].max() )
    f.write("\tN(HI) minimum (cm-2) = %s\n"%hicld[hicld!=0].min() )
    f.write("\tN(HI) mean (cm-2) = %s\n"%hicld[hicld!=0].mean() )
    f.write("\tN(HI) standard deviation (cm-2) = %s\n"%hicld[hicld!=0].std() )
    f.write("\tN(HI) sum (cm-2) = %s\n"%hicld[hicld!=0].sum() )
    f.write("\n\n")
    coco = pyfits.PrimaryHDU( cocld )
    coco.header = co[0].header
    hihi = pyfits.PrimaryHDU( hicld )
    hihi.header = hi[0].header
    dgdg = pyfits.PrimaryHDU( dgcld )
    dgdg.header = dg[0].header
    coco.writeto('Clouds_%ssig/CO_cloud%s.fits'%(sigma,cldnum))
    hihi.writeto('Clouds_%ssig/HI_cloud%s.fits'%(sigma,cldnum))
    dgdg.writeto('Clouds_%ssig/DG_cloud%s.fits'%(sigma,cldnum))
    







glon = 315
glat = -30
sigma = 2

#os.chdir('../l%sb%s'%(int(glon),int(glat)))
#os.mkdir('Clouds_%ssig'%sigma)
os.chdir('/home/abrahams/HICO_survey/Planck/PlanckDR1/COMaps/')

# go into appropriate directory, which contains the spatial templates
# Cut out some of the fluff elsewhere:
#	perform a combination of erosion, and dilation. Erode 
#	the image to further kill unwanted things. Then dilate 
#	again to recover the edges killed. Sure, we lose a bit 
#	of information, but we get rid of a LOT of fluff.
co = pyfits.open('COsky_tenthdeg_gt25.fits')
coer=pyfits.open('COsky_tenthdeg_error.fits')
ratio = co[0].data/coer[0].data

co[0].data[ ratio < sigma ] = 0
ratio = None
ny = co[0].header['NAXIS2']
nx = co[0].header['NAXIS1']
#co.writeto('CO_temp_1sig.fits')

#co = pyfits.open('CO_temp_1sig.fits')
#opacity = 'thick'
#dg = pyfits.open('Ebv_%s.fits'%opacity)
#hi = pyfits.open('HI_%s.fits'%opacity)
#coWCS = astWCS.WCS(co[0].header,mode="pyfits")

glon = co[0].header['CDELT1']*( np.arange( co[0].header['NAXIS1'] ) - (co[0].header['CRPIX1'])) + co[0].header['CRVAL1']
glat = co[0].header['CDELT2']*( np.arange( co[0].header['NAXIS2'] ) - (co[0].header['CRPIX2'])) + co[0].header['CRVAL2']

co[0].data = np.transpose( co[0].data )
comsk = np.zeros( (nx,ny) )
comsk[ co[0].data > 0 ] = 1


#comsk = ndimage.binary_erosion( comsk )
#for x in np.linspace( 1,nx-2,nx-2 ):
#    for y in np.linspace( 1,ny-2,ny-2 ):
#    	if ( (comsk[ x+1,y] == 0) and (comsk[ x-1,y] == 0) and (comsk[ x,y+1] == 0) and (comsk[ x,y-1] == 0) ):
#	    comsk[x,y] = 0

#comsk = ndimage.binary_dilation( comsk )

co[0].data *= comsk
cldnum = 0
f = open('CloudProperties_%ssig.txt'%sigma,'a')
# Iterate: find CO max, puff outward
while (co[0].data.sum() > 0):
    cldnum += 1
    maxx = np.where( co[0].data == co[0].data.max() )[0][0]
    maxy = np.where( co[0].data == co[0].data.max() )[1][0]
    cldmsk = np.zeros( (nx,ny) )
    cldmsk[maxx,maxy] = 1
    while True:
	diff_temp = comsk - cldmsk
	cldmsk = ndimage.binary_dilation(cldmsk)
	cldmsk *= comsk
	diff = comsk - cldmsk
	if ( (diff - diff_temp).sum() == 0 ):
	    break

# So ... justification for this cutoff? My routine will also combine 
# some clouds. All right?
    if ( cldmsk.sum() > 40):
#	CldResults(cldmsk,hi,co,dg,cldnum,coWCS,f)
	# f: cloud number || x || y || glon || glat || area || COpeak
	f.write('%s %s %s %s %s %s %s\n'%(cldnum,maxx+1,maxy+1,glon[maxx],glat[maxy],cldmsk.sum(),co[0].data.max()))

    co[0].data -= cldmsk*co[0].data


f.close()


