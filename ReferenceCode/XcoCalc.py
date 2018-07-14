
import os, pyfits
import numpy as np

glon = 159
glat = -34
cloud = 1
opacity = 'thick'

xco = 1.23e19
xav = 8.17e20

os.chdir('/home/abrahams/HICO_survey/SourceSearch/l%sb%s/'%(glon,glat))

co = pyfits.open('CO_temp.fits')
dg = pyfits.open('Ebv_%s.fits'%opacity)
dg[0].data[ dg[0].data < 0 ] = 0.0
comsk = pyfits.open('Clouds_2sig/CO_cloud%s.fits'%cloud)
comsk[0].data[ comsk[0].data != 0 ] = 1.0
co[0].data *= comsk[0].data
dg[0].data *= comsk[0].data

XCO_ = xco + xav*(dg[0].data.sum()/co[0].data.sum())
print "Cloud at (%s,%s)"%(glon,glat)
print "Xco' = %s"%XCO_

