
import numpy as np
import pyfits
import os

def distance(x,y,x0,y0):
    return np.sqrt( (x-x0)**2 + (y-y0)**2 )


glon = 159
glat = -34
cld  = 1
# Cloud number (cld) is 1 for most ... but ... 
# (72, -43) = cloud 4
os.chdir('../l%sb%s'%(int(glon),int(glat)))


co = pyfits.open('CO_temp_all.fits')
co_= pyfits.open('CO_temp_error.fits')
cocld = pyfits.open('Clouds_2sig/CO_cloud%s.fits'%cld)
cocld1= pyfits.open('Clouds_1sig/CO_cloud%s.fits'%cld)

for i in [1,2,3]:
    co[0].data[ co[0].data < i*co_[0].data ] = 0
    fl = 'CO_temp_%ssig.fits'%i
    if ( fl in os.listdir('.') ):
	os.remove(fl)

    co.writeto(fl)
    

co1 = pyfits.open('CO_temp_1sig.fits')
co2 = pyfits.open('CO_temp_2sig.fits')
co3 = pyfits.open('CO_temp_3sig.fits')


co3[0].data[ cocld[0].data == 0 ] = 0
co2[0].data[ cocld[0].data == 0 ] = 0


dst = 0
(nx,ny) = (co2[0].header['NAXIS1'],co2[0].header['NAXIS2'])
dst = np.zeros( (nx,ny) )
for i in np.arange( nx ):
    for j in np.arange( ny ):
	dst[i,j] = distance(i,j,100,100)

dst_mx = dst[ co2[0].data > 0 ].max()


area3 = co3[0].data[ co3[0].data > 0 ].size/100.
area2 = co2[0].data[ co2[0].data > 0 ].size/100.

cocld1[0].data[ cocld1[0].data > 0 ] = 1
cocld1[0].data[ dst*cocld1[0].data > dst_mx ] = 0

area1 = cocld1[0].data[ cocld1[0].data > 0 ].size/100.

print "Angular size = %s - %s + %s"%(area2, (area2-area3), (area1-area2))


