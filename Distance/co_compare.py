
import numpy as np
import pyfits
from astLib import astWCS,astCoords
import json, requests
import os

def query(lon,lat,coordsys='gal',mode='full'):
	"""
	send a line-of-sigh reddening query to the Argonaut web server.
	Inputs:
	 lon, lat: longitude and latitude, in degrees.
	 coordsys: 'gal' for Galactic, 'equ' for Equatorial (J2000).
	 mode: 'full', 'lite', or 'sfd'

	In 'full' mode, outputs a dictionary containing, among other things:
	 'distmod':	The distance moduli that define the distance bins.
	 'best':	The best-fit (maximum probability density) 
			line-of-sight reddening, in units of SFD-equivalent
			E(B-V), to each distance modulus in 'distmod.' See
			Schlafly & Finkberiner (2011) for a definition of the 
			reddening vector (use R_V = 3.1).
	 'samples':	Samples of the line-of-sight reddening, drawn from
			the probability density on reddening profiles.
	 'success':	1 if query succeeded, and 0 otherwise.
	 'converged':	1 if the line-of-sight reddening fit converged, and
			0 otherwise.
	 'n_stars':	@ of stars used to fit the line-of-sight reddening.
	 'DM_reliable_min':	Minimum reliable distance modulus in pixel.
	 'DM_reliable_max':	Maximum reliable distance modulus in pixel.

	Less infomation is returned in 'lite' mode, while in 'sfd' mode,
	the Schlegel, Finkbeiner & Davis (1998) E(B-V) is returned.
	"""
	url = 'http://argonaut.skymaps.info/gal-lb-query-light'
	payload = {'mode':mode}
	if coordsys.lower() in ['gal','g']:
		payload['l'] = lon
		payload['b'] = lat
	elif coordsys.lower() in ['equ','e']:
		payload['ra'] = lon
		payload['dec'] = lat
	else:
		raise ValueError("coordsys '{0}' not understood.".format(coordsys))

	headers = {'content-type':'application/json'}
	r = requests.post(url, data=json.dumps(payload), headers=headers)
	
	try:
		r.raise_for_status()
	except requests.exceptions.HTTPError as e:
		print('Response received from Argonaut:')
		print(r.text)
		raise e

	return json.loads(r.text)


##############################################################################
def distWrite(l,b,flname):
	glon = l + np.arange( -2, 2.25, 0.25 )
	glat = b + np.arange( -2, 2.25, 0.25 )
	lon = []
	lat = []
	for i in range( len(glon) ):
		for k in range( len(glat) ):
			lon.append( glon[i] )
			lat.append( glat[k] )

	qresult = query(lon,lat)
	dmod = np.array(qresult['distmod'])
	f = open(flname,'a')
	for i in range( len(lon) ):
		if ( qresult['converged'][i] == 1 ):
			best = np.array(qresult['best'][i])
			dav = best[1:] - best[:-1]
			dav = np.append(dav,0)

			dmod_max = dmod[ dav == dav.max() ]
			dmod_range = dmod[[np.where(dav==dav.max())[0][0]-1,np.where(dav==dav.max())[0][0]+1]]
			dist = 10**( (dmod_max/5.)+1 )
			f.write('(%s,%s),converge = %s\t\t%s max:%s min:%s\n'%(lon[i],lat[i], 1, dist,10**( dmod_range[0]/5.+1 ), 10**(dmod_range[1]/5.+1)))
		else:
			f.write('(%s,%s),converge = %s\t\t[ 000.0000] max:000.00000 min:000.00000\n'%(lon[i],lat[i],0))



	f.close()


l = 353
b = 17
flname = 'Dist_2deg.txt'
if ( flname in os.listdir('.') ):
	os.remove(flname)

distWrite(l,b,flname)


f = open(flname,'r')
dat = f.readlines()
f.close()
#coim = pyfits.open('/home/abrahams/HICO_survey/SourceSearch/l%sb%s/CO_temp.fits'%(l,b))
coim = pyfits.open('/home/abrahams/HICO_survey/SourceSearch/RhoOph/l353b17/ro_pass8_files/Aug7_Ryan_files/CO_temp.fits')
co_wcs = astWCS.WCS(coim[0].header,mode='pyfits')
coim[0].data = np.transpose( coim[0].data )

d = np.array([])
co = np.array([])

for i in range( len(dat) ):
	glon = float( dat[i].split(',')[0][1:] )
	glat = float( dat[i].split(',')[1][:-1] )
	(glon,glat) = astCoords.convertCoords("GALACTIC","J2000",glon,glat,2000)
	try:
		converge = int( dat[i].split('converge = ')[1] )
		d = np.append( d, 0 )
		continue
	except:
		try:
			d = np.append(d, float(dat[i].split(',')[2].split('[ ')[1].split(']')[0]))
		except:
			d = np.append(d,0)


	(x,y) = co_wcs.wcs2pix(glon,glat)
	co = np.append( co, coim[0].data[x,y] )

if ( len(co) == len(d)-1 ):
	co = np.append(co,0)


wmean = (d*co).sum()/co.sum()
