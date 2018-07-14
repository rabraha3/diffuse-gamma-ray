
import numpy as np
import pyfits

m = pyfits.open('CountSpec_all.fits')
n_ebn = m[1].header['NAXIS2']

emn = np.zeros( n_ebn )
emx = np.zeros( n_ebn )
obs = np.zeros( n_ebn )
CO  = np.zeros( n_ebn )
DG  = np.zeros( n_ebn )
HI  = np.zeros( n_ebn )
iso = np.zeros( n_ebn )
mdl = np.zeros( n_ebn )

n_hdr = len(m[1].header)
n_src = len(m[1].data[0])
for i in np.arange( n_ebn ):
    obs[i] = m[1].data[i][0]
    CO[i]  = m[1].data[i][-4]
    DG[i]  = m[1].data[i][-3]
    HI[i]  = m[1].data[i][-2]
    iso[i] = m[1].data[i][-1]
    emn[i] = m[3].data[i][0]
    emx[i] = m[3].data[i][1]
    mdl[i] = np.sum( m[1].data[i][1:] )

ebn = np.sqrt(emn*emx)
e_er= 0.5*(emx-emn)


plt.errorbar(ebn,obs,xerr=e_er,yerr=np.sqrt(obs),fmt='o',color='black',label='Data')
plt.plot(ebn,mdl,label='Model')
plt.plot(ebn,HI,label='HI')
plt.plot(ebn,iso,label='isotropic')


