
import numpy as np

HI_int = 6.139e-27
HI_int_= 0.3112e-27
HI_int_= np.sqrt( HI_int_**2 + (0.17*HI_int)**2 + (0.15*HI_int)**2 )

CO_int = 1.369e-6
CO_int_= 0.1e-6
CO_int_= np.sqrt( CO_int_**2 + (0.05*CO_int)**2 + (0.15*CO_int)**2 )

Av_int = 4.331e-5
Av_int_= 0.4976e-5
AV_int_= np.sqrt( Av_int_**2 + (0.08*Av_int)**2 + (0.15*Av_int)**2 )


Xco = CO_int/(2*HI_int)
sig_xco = np.sqrt( (CO_int_/(2*HI_int))**2 + (CO_int*HI_int_/(2*HI_int**2))**2 )
Xav = Av_int/HI_int
sig_xav = np.sqrt( (Av_int_/(HI_int))**2 + (Av_int*HI_int_/(HI_int**2))**2 )




print "Xco = %s +/- %s"%(Xco,sig_xco)
print "Xav/2 = %s +/- %s"%(Xav/2,sig_xav/2)




