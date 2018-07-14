
from Type_count import Cloud
import numpy as np

def quer(cld):
	""" If Schlafly et al. (2014) tables do not contain the distance
	    information, then we should query their data sube 
	"""
	



f = open('DistNeed.txt','a')


HI = np.array([])
CO = np.array([])
DG = np.array([])
#dist = np.array([])
glon = np.array([])
glat = np.array([])
for i in range( 106 ):
	try:
		cld = Cloud(i)
		cld.getTs()
		cld.getEmissivity()
#		cld.getDistance()
	except:
		continue

	glon = np.append(glon, cld.glon)
	glat = np.append(glat, cld.glat)
	try:
		HI = np.append(HI, float(cld.emissivity['HI']))
	except:
		HI = np.append(HI, 0)

	try:
		CO = np.append(CO, float(cld.emissivity['CO']))
	except:
		CO = np.append(CO, 0)

	try:
		DG = np.append(DG, float(cld.emissivity['E(B-V)']))
	except:
		DG = np.append(DG, 0)


#	if ( cld.dist == None):
#		f.write('(%s,%s)\n'%(cld.glon,cld.glat))

#	try:
#		dist = np.append(dist, np.mean(cld.dist))
#	except:
#		dist = np.append(dist, 0)


f.close()
