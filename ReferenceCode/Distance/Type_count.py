
import numpy as np
##################################################
## 	Need to grab emissivity errors next	##





class Cloud(object):
	""" Grab some properties about a cloud, listed in "CldNum.txt"
	    cld = Cloud(number) ## will create a Cloud object and look at cloud #number,
				   as listed in "CldNum.txt"
	    The following methods are available:
		cld.getLike()		# Retrieve likelihoods from "LogLike.dat".
		cld.getTs()		# Calculate TS values from likelihoods.
		cld.getDistance()	# Find distance from Schlafly et al. (2014).
					# Call: cld.getDistance(mode='lallement'). Without specifying mode,
					# method defaults to Schlafly distance.
		cld.getCoords()		# Print coordinates of ROI center.
		cld.getEmissivity()	# Retrieve emissivity from "Properties.txt".
					# and component TS.
		cld.peakComponents()	# Find peak HI, CO, and/or E(B-V) emission IN cloud.
	    The following attributes are available:
		cld.like		# Likelihoods dictionary, keys = analysis modes.
		cld.ts			# Overall ROI TS dictionary, keys found in cld._allTs.
		cld.dist		# Distance to cloud, as found in .getDistance().
		cld.distError		# Standard deviation on the distance measurement.
		cld.emissivity		# Emissivity dictionary, keys = diffuse components.
		cld.emissivityTs	# Component TS dictionary, keys = diffuse components.
		cld.home		# Home directory of analysis.
		cld.path		# Path to directory of cloud.
		cld.peak		# Peak HI, CO, E(B-V) emission.
	    Hidden methods:
		cld._coordRead()	## Read "CldNum.txt" to find GLON, GLAT.
		cld._angularDist()	## Find angular distance between two points on sphere.
		cld._schlaflyDistance()	## Query Schlafly E(B-V) cube and extract largest E(B-V)
					   jump for distance.
		cld._query(l,b)		## Interface with Harvard website (argonaut) to grab data.
					   for schlafly distance.
		cld._getEmissivityScale	## Emissivity scale, from integrating Casandjian's spectra from Emin-Emax
	    Hidden attributes:
		cld._allModes		## All possible modes of likelihood analysis.
		cld._component		## Diffuse components in initial model.
		cld._allTs		## Calculation of TS's.
	"""
	def __init__(self, cldnum):
		""" Initialize: which cloud number of we looking at? Grab coordinates and 
		    set up the properties we care about.
		    _allModes tells us all POSSIBLE modes of analysis
		    _allTs tells us all POSSIBLE TS values
		"""
		self.cloudNum = cldnum
		self.home = '/home/abrahams/HICO_survey/SourceSearch/'
		(self.glon,self.glat) = self._coordRead(cldnum)
		self.path = self.home + 'l%sb%s/ebin250_10000/'%(int(self.glon), int(self.glat))
		self.like = {}
		self.ts = {}
		self.emissivity = {}
		self.emissivityTs={}
		self.emissivityError = {}
		self.emissivityScale = {'HI':1,'HI_far':1,'CO':1,'E(B-V)':1}
		self._allModes = ['all','nocodg','noco','nodg','nocld','ptsrc','all_plus']
		self._component= ['HI','HI_far','CO','E(B-V)','IC','iso_P8R2_SOURCE_V6_v06']
		self._allTs = {'TS_H2':'nocodg-all', 'TS_CO':'noco-all', 'TS_DG':'nodg-all', 'TS_CLD':'nocld-all', 'TS_EX':'ptsrc-all', 'TS_AGN':'all_plus-all'}


	#################################################################
	def _coordRead(self,line):
		""" Read coordinates of cloud "line" from a file. Use:
		"""
		f = open(self.home+'Analysis/CldNum.txt','r')
		dat = f.readlines()
		f.close()
		lon = float( dat[line].strip('\n').split(',')[0].split('(')[1] )
		lat = float( dat[line].strip('\n').split(',')[1].split(')')[0] )
		if ( lon > 180 ):
			lon = lon-360
		#print "Great! Running cloud %s, at (l,b)=(%s,%s)"%(line,lon,lat)
		return (lon,lat)


	#################################################################
	def getCoords(self,RA=False):
		""" Print out longitude, latitude in Galactic (default) or ecliptic
		    cld.getCoords(), or cld.getCoords(RA=False) if in Galactic
		    cld.getCoords(RA=True) if J2000 ecliptic coordinates. 
		    No support for other coordinate systems.
		"""
		if ( RA ):
			from astLib import astCoords
			return astCoords.convertCoords("GALACTIC", "J2000", self.glon, self.glat, 2000)
		else:
			return (self.glon,self.glat)


	#################################################################
	def _findStart(self,x,mode='all'):
		""" Find the line where we should start counting likelihoods in
		    LogLike.dat . . . the analysis may have been run a number of
		    times...
		"""
		lookingFor = " " + mode + " "
		for i in xrange(len(x)-1,-1,-1):
			if ( ( (lookingFor in x[i]) and ('all_plus' not in x[i])) or (i == 0) ):
				return i
			

	#################################################################
	def getLike(self,mode = None):
		""" Grab the log(likelihoods) 
		    note: this does NOT keep track of repeated
		    analyses. It ONLY saves the latest analysis of MODE.
		    If mode =/= None [[ .like(mode='all') ]], we search for
		    the likelihood of ONLY that mode. Otherwise, ALL modes.
		    Use: for ALL modes available, 
			 cld.getLike(), or cld.getLike(mode=None)
			 for only a single mode, example:
			 cld.getLike(mode='all')
		"""
		f = open(self.path+'LogLike.dat','r')
		dat = f.readlines()
		f.close()

		istart = self._findStart(dat,mode='all')
		
		for i in range( istart,len(dat) ):
			if ( mode != None ):
				a = mode
			elif ( ' of ' in dat[i] ):
				a = dat[i].split(' of ')[1].split(' = ')[0]
			else:
				continue

			self.like[a] = float( dat[i].split('of')[1].split(' = ')[1] )

		## last pass through all modes: if mode hasn't been done, set likelihood to NaN
		for md in self._allModes:
			try:
				self.like[ md ] > 0
			except KeyError:
				self.like[ md ] = np.nan


	#################################################################
	def getTs(self,mode = None):
		""" Grab the OVERALL TS value from likelihood ratio """
		from numpy import fabs
		if ( (len(self.like) == 0) and (mode == None) ):
			self.getLike()
		elif ( (len(self.like) == 0) and (mode != None) ):
			self.getLike(mode='all')
			self.getLike(mode=mode)

		if (mode != None):
			a = [mode,' ']
		else:
			a = self.like.keys()

		a.sort()

		for i in range( 1,len(a) ):
			tmp = a[i].lower()
			if ( tmp == 'nocodg' ):
				mode_ts = 'TS_H2'
			elif ( tmp == 'noco' ):
				mode_ts = 'TS_CO'
			elif ( tmp == 'nodg' ):
				mode_ts = 'TS_DG'
			elif ( tmp == 'nocld' ):
				mode_ts = 'TS_CLD'
			elif ( tmp == 'ptsrc' ):
				mode_ts = 'TS_EX'
			elif ( tmp == 'all_plus' ):
				mode_ts = 'TS_AGN'
			else:
				return None


			self.ts[mode_ts] = round(fabs( 2*(self.like[a[i]] - self.like['all']) ),4)


	#################################################################
	def getEmissivity(self,mode=None):
		""" Get emissivity of each diffuse component (HI,HI_far, CO, E(B-V), IC, iso)"""
		f = open(self.path+'Properties.txt','r')
		dat = f.readlines()
		f.close()

		
		for i in range( len(dat) ):
			for k in range( len(self._component) ):
				comp = self._component[k]
				if ( comp.lower() == 'hi_far' ):
					self.emissivityScale[comp] = self._getEmissivityScale('hi')
				elif ( comp.lower() in ['hi','co','e(b-v)'] ):
					self.emissivityScale[comp] = self._getEmissivityScale(comp)
				if ( i >= len(dat)-2 ):
					continue

				if ( (comp in dat[i]) and ('Spectrum' in dat[i+1]) ):
					self.emissivity[comp] = float(dat[i+2].split(':  ')[1].split('  ')[0])
					self.emissivityError[comp] = float(dat[i+2].split(':  ')[1].split('  ')[1])

				elif ( (comp in dat[i]) and ('TS' in dat[i+1]) ):
					self.emissivityTs[comp] = float(dat[i+1].split(' = ')[1])

		for line in dat[:3]:
			if ( 'R = ' in line ):
				self.dustToGas = float(line.split(' = ')[1])
			elif ( 'q = ' in line ):
				self.xco = float( line.split(' = ')[1] )

		
		

	#################################################################
	def _getEmissivityScale(self,component):
		""" Energy range: 250 MeV -- 10 GeV 
		    Read out Casandjian's LIS for HI, CO, or DG. Integrate
		    via trapezoid rule over our energy range to calculate
		    emissivity scaling.
		    use: self._getEmissivityScale(comp),
			 comp = 'HI', 'CO', or 'E(B-V)' to get the scale
			 for each component, respectively.
		    If the emissivity spectrum is not available, return 
		    error.
		"""
		from numpy import array,append
		fl = {'hi':'LIS_cas.txt', 'co':'qCO_cas.txt', 'e(b-v)':'qdg_cas.txt'}
		if ( component.lower() not in fl.keys() ):
			return "Not valid component"

		try:
			f = open(self.home+'Analysis/%s'%fl[component.lower()],'r')
		except:
			return "No emissivity spectrum for %s"%component

		lis = f.readlines()
		f.close()
		ee = array([])
		emis=array([])
		for i in range( len(lis) ):
			E = float(lis[i].split(' ')[0])
			if ( (E >= 240) and (E <= 10100) ):
				# For rounding purposes ...
				ee = append(ee,E)
				emis=append( emis, float(lis[i].split(' ')[1]) )

		return (0.5*(emis[1:]+emis[:-1])*(ee[1:]-ee[:-1])).sum()
		

	#################################################################
	def getDistance(self,mode='schlafly'):
		""" Find if Schlafly et al. (2014) found the distances to the cloud.
		    Return [min,max] distances in (pc).
		"""
		import pyfits, os
		from numpy import mean
		self._numClds = 0
		if ( mode.lower() not in ['schlafly','schlegel'] ):
			f = open(self.home + 'Analysis/CldsDist_lal.dat','r')
			dat = f.readlines()
			f.close()
			a = dat[self.cloudNum].strip('\n').strip(' ').split('    ')
			b = []
			for num in a:
				try:
					b.append(float(num))
				except:
					continue

			self.dist = b[2]
			self.distError = b[-1]
			return

		try:
			if ( 'mbmcloud.fits' in os.listdir('.') ):
				mbm = pyfits.open('mbmcloud.fits')
			else:
				mbm = pyfits.open('Distance/mbmcloud.fits')
		except:
			self.dist = None
			raise ValueError("File mbmcloud.fits not found.")

		dist_possible = []
		for i in range( len(mbm[1].data) ):
			if ( self._angularDist(mbm[1].data['l'][i],mbm[1].data['b'][i],self.glon,self.glat) < 1.0 ):
				a = 10**( (mbm[1].data['m16'][i]/5.)+1 )
				b = 10**( (mbm[1].data['m84'][i]/5.)+1 )
				dist_possible.append( [a,b] )
				self._numClds += 1
		
		del mbm
		try:
			if ( 'bigcloud.fits' in os.listdir('.') ):
				mbm = pyfits.open('bigcloud.fits')
			else:
				mbm = pyfits.open('Distance/bigcloud.fits')
		except:
			self.dist = None
			raise ValueError("File bigcloud.fits not found.")

		for i in range( len(mbm[1].data) ):
			if ( self._angularDist(mbm[1].data['l'][i],mbm[1].data['b'][i],self.glon,self.glat) < 1.0 ):
				a = 10**( (mbm[1].data['m16'][i]/5.)+1 )
				b = 10**( (mbm[1].data['m84'][i]/5.)+1 )
				dist_possible.append( [a,b] )
				self._numClds += 1

		if ( len(dist_possible) == 0 ):
			(dst,std) = self._schlaflyDistance()
			self.dist = [dst,dst]
			self.distError = std
		else:
			self.dist = mean( dist_possible, axis=0 )
			self.distError = 0.5*(dist_possible[0][1]-dist_possible[0][0])
		

	#################################################################
	def _schlaflyDistance(self):
		""" Query Schlafly et al. (2014) data cube and calculate
		    the distance from the differential E(B-V) data. 
		    Take a 4x4 degree box around the ROI center and find
		    the distances, weighted by CO emission. 
		"""
		import numpy as np
		import pyfits
		from astLib import astWCS
		self.coim = pyfits.open(self.path+'../CO_temp.fits')
		co_wcs = astWCS.WCS(self.coim[0].header,mode='pyfits')
		self.coim[0].data = np.transpose( self.coim[0].data )

		glon = self.glon + np.arange( -2, 2.1, 0.1 )
		glat = self.glat + np.arange( -2, 2.1, 0.1 )
		lon = []
		lat = []
		for i in range( len(glon) ):
			for k in range( len(glat) ):
				lon.append( glon[i] )
				lat.append( glat[k] )

		del glon,glat	##################
		qresult = self._query(lon,lat)
		dmod = np.array(qresult['distmod'])
		dist = np.array([])
		dist_= np.array([])
		co = np.array([])
		for i in range( len(lon) ):
			if ( qresult['converged'][i] == 1 ):
				best = np.array(qresult['best'][i])
				dav = best[1:] - best[:-1]
				dav = np.append(dav,0)

				dmod_max = dmod[ dav == dav.max() ]

				d1 = 10**( (dmod[np.where(dav==dav.max())[0][0]-1]/5.) + 1 )
				d2 = 10**( (dmod[np.where(dav==dav.max())[0][0]+1]/5.) + 1 )
				dist_ = np.append( dist_, np.mean( [d1,d2] )/2. )
				dist = np.append(dist, 10**( (dmod_max/5.)+1 ) )
				(x,y) = co_wcs.wcs2pix(lon[i],lat[i])
				co = np.append(co, self.coim[0].data[x,y] )
			elif ( qresult['converged'][i] == 0 ):
				continue

		if ( len( dist ) > 0 ):
			wmean = (dist*co).sum()/co.sum()
			wstd  = np.sqrt( ((dist_**2*co**2).sum())/(co.sum()**2) )
			return (wmean,wstd)
		else:
			return (0,0)
		


	#################################################################
	def _query(self,lon,lat,coordsys='gal',mode='full'):
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
		import json,requests
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

	#################################################################
	def _angularDist(self,l,b,l0,b0):
		""" Calculate the angular distance between two points on the sky.
		    (l,b) used because things will be given in Galactic coordinates
		    (degrees).
		"""
		from math import cos,acos,sin,fabs,pi
		conv = (pi/180.)
		l *= conv
		b *= conv
		l0*= conv
		b0*= conv
		angle = acos( sin(b)*sin(b0) + cos(b)*cos(b0)*cos(fabs(l-l0)) )
		return (1/conv)*angle

	#################################################################
	def peakComponents(self,comp=['hi','co','ebv']):
		""" Extract out peak Wco, NHI, and E(B-V), any or all. To do this, it
		    opens the FITS template files, considers only the inner 2x2 degree
		    region, and get the Value of the maximum.
		    Use:
		    (hi,co,dg) = __.peakComponents(comp = ['hi','co'])
				If you only want the HI and CO components
		    The "comp" option can be 'hi','co',or 'ebv' or 'HI','WCO','E(B-V)'
		    Output = [#,#,#], where the length of list depends on length of "comp",
			     i.e., the number of components we want in the order we specify.
		"""
		import pyfits
		hydro = ['hi','hydrogen']
		cmono = ['co','wco','carbon monoxide']
		dgas  = ['ebv','e(b-v)']
		out = [0]*len(comp)
		error = [0]*len(comp)
		comp = map(str.lower,comp)
		for cc in range(len(comp)):
			if ( (comp[cc] in hydro) ):
				im = pyfits.open(self.path+'../HI_thick.fits')
			elif ( comp[cc] in cmono ):
				im = pyfits.open(self.path+'../CO_temp.fits')
				err = pyfits.open(self.path+'../CO_temp_error.fits')
			elif ( comp[cc] in dgas ):
				im = pyfits.open(self.path+'../Ebv_thick.fits')
				err = pyfits.open(self.path+'../Dust_err.fits')
			else:
				continue
		
			im_size = [im[0].header['NAXIS1'],im[0].header['NAXIS2']]
			if ( comp[cc] in cmono ):
				out[cc] = im[0].data[int(im_size[1]/2.-10):int(im_size[1]/2.+10), 
					    int(im_size[0]/2.-10):int(im_size[0]/2.+10)].max()
				error[cc] = err[0].data[int(im_size[1]/2.-10):int(im_size[1]/2.+10), 
					    int(im_size[0]/2.-10):int(im_size[0]/2.+10)].mean()
			elif ( comp[cc] in dgas ):
				out[cc] = im[0].data[int(im_size[1]/2.-10):int(im_size[1]/2.+10), 
					    int(im_size[0]/2.-10):int(im_size[0]/2.+10)].max()
				error[cc] = err[0].data[int(im_size[1]/2.-10):int(im_size[1]/2.+10), 
					    int(im_size[0]/2.-10):int(im_size[0]/2.+10)].mean()
			elif ( comp[cc] in hydro ):
				out[cc] = im[0].data[int(im_size[1]/2.-10):int(im_size[1]/2.+10), 
					    int(im_size[0]/2.-10):int(im_size[0]/2.+10)].max()

		self.peak = out
		self.peakError = error

