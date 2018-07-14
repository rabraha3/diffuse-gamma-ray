

import json, requests, pyfits
from astLib import astWCS
import numpy as np

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


##########################################################
###### Find distance to enumerated/specified clouds
#f = open('DistNeed.txt','r')
#data = f.readlines()
#f.close()
#glon = []
#glat = []
#num_clds = len(data)
#for i in range( num_clds ):
#	glon.append( float(data[i].split(',')[0][1:]) )
#	glat.append( float(data[i].split(',')[1][:-2]))
##########################################################

glon = 353 + np.arange( -2, 2.25, 0.25 )
glat = 17 + np.arange( -2, 2.25, 0.25 )
lon = []
lat = []
for i in range( len(glon) ):
	for k in range( len(glat) ):
		lon.append( glon[i] )
		lat.append( glat[k] )

qresult = query(lon,lat)
dmod = np.array(qresult['distmod'])
f = open('Distance/Dist_mbm12_2deg.txt','a')
for i in range( len(lon) ):
	#dist = 10**( (np.array(qresult['distmod'][i])/5.)+1 )
	if ( qresult['converged'][i] == 1 ):
		best = np.array(qresult['best'][i])
		dav = best[1:] - best[:-1]
		dav = np.append(dav,0)

		dmod_max = dmod[ dav == dav.max() ]
		dmod_range = dmod[[np.where(dav==dav.max())[0][0]-1,np.where(dav==dav.max())[0][0]+1]]
		dist = 10**( (dmod_max/5.)+1 )
#		dmax = dist[dav == dav.max()]
#		dmax_= dist[ dav > dav.max()*0.75 ]		# If len(dmax_) == 1, ignore
		f.write('(%s,%s),converge = %s\t\t%s max:%s min:%s\n'%(lon[i],lat[i], 1, dist,10**( dmod_range[0]/5.+1 ), 10**(dmod_range[1]/5.+1)))
	else:
		f.write('(%s,%s),converge = %s\n'%(lon[i],lat[i],0))



f.close()




