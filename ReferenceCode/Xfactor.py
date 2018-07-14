
import numpy as np
import pyfits, os

# make sure we're in coorect directory
glon = 315
glat = -30
cld = 1
# Cloud 4 for ( 72,-43)
# Cloud 2 for (210,-36)
os.chdir('../Pilot/l%sb%s'%(int(glon),int(glat)))

Xco = 11.1500244339e19
Xco_= 3.23745389426e19
Xav = 35.2744746701e20
Xav_= 9.14214938931e20
#Nhi = 9.38
#Nco = 0.19
#Nav = 1.15


co = pyfits.open('Clouds_2sig/CO_cloud%s.fits'%cld)
#dg = pyfits.open('Clouds_2sig/DG_cloud%s.fits'%cld)
#hi = pyfits.open('Clouds_2sig/HI_cloud%s.fits'%cld)

#f_up = ( (Nco + Nav) / (Nhi + Nco + Nav) )
#f_dn = ( Nco / (Nhi + Nco + Nav) )
#f = 0.5*(f_up + f_dn)

#Xcomp = Xco + (f*dg[0].data[ dg[0].data > 0 ].mean()/co[0].data[ co[0].data > 0 ].mean() ) *Xav

#Xcomp_up = Xco + (f_up*dg[0].data[ dg[0].data > 0 ].mean()/co[0].data[ co[0].data > 0 ].mean() ) *Xav
#Xcomp_dn = Xco + (f_dn*dg[0].data[ dg[0].data > 0 ].mean()/co[0].data[ co[0].data > 0 ].mean() ) *Xav

#print Xcomp_dn
#print Xcomp_up
#print Xcomp

co[0].data[ co[0].data > 0 ] = 1
co_ = pyfits.open('CO_temp_all.fits')
co_[0].data *= co[0].data
dg = pyfits.open('Ebv_thick.fits')
dg[0].data *= co[0].data


#### FOR ALL
# co_[0].data[ co_[0].data < 0.0 ] = 0
#### #######

Xcomp_up = Xco
Xcomp_dn = Xco + (dg[0].data[ dg[0].data > 0 ].mean()/co[0].data[ co[0].data > 0 ].mean() ) *Xav

Xc = Xco + (dg[0].data/co_[0].data ) *Xav
#Xc[ np.isnan(Xc) ] = 0

sig_Xc = np.sqrt( Xco_**2 + (Xav_*dg[0].data/co_[0].data)**2 + (Xav*.16*dg[0].data/co_[0].data)**2 + (Xav*dg[0].data/co_[0].data**2)**2 )


#### For ALL
Xc[ np.isinf(Xc) ] = 'nan'
sig_Xc[ np.isinf(sig_Xc) ] = 'nan'



print "Max = %s +/- %s"%(Xc[ np.isnan(Xc) == False ].max(),sig_Xc[ np.isnan(Xc) == False ].max())
print "Mean= %s +/- %s"%(Xc[ np.isnan(Xc) == False ].mean(),sig_Xc[ np.isnan(Xc) == False ].mean())
print "Min = %s +/- %s"%(Xc[ np.isnan(Xc) == False ].min(),sig_Xc[ np.isnan(Xc) == False ].min())

# mean - min
minus = Xc[ np.isnan(Xc) == False ].mean() - Xc[ np.isnan(Xc) == False ].min()
# max - mean
maxim = Xc[ np.isnan(Xc) == False ].max() - Xc[ np.isnan(Xc) == False ].mean()

print "Xc = %s + %s - %s"%(Xc[np.isnan(Xc)==False].mean(),maxim,minus)

##mm = pyfits.PrimaryHDU(Xc)
##mm.header = dg[0].header
##mm.writeto('Xfactor.fits')
#
##mm_ = pyfits.PrimaryHDU(sig_Xc)
##mm_.header = dg[0].header
##mm_.writeto('Xfactor_err.fits')
# Here are some aplpy commands

# fig = plt.figure()
# im1 = aplpy.FITSFigure('Xco_comp.fits',figure=fig, subplot=[0.1, 0.5, 0.38, 0.38], convention='calabretta')
# im1.show_colorscale(vmin=1e18,vmax=1e20)
