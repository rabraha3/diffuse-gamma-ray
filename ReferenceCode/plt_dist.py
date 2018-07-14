
import numpy as np

f = open('Clds_Distance_1.0deg.txt','r')
g = open('CldNum.txt','r')
h = open('Clds_Emissivity.txt','r')

datDist = f.readlines()
datCoor = g.readlines()
datEmis = h.readlines()
f.close()
g.close()
h.close()

HI = np.array([])
HI_= np.array([])
CO = np.array([])
CO_= np.array([])
dist = np.array([])
DG = np.array([])
DG_= np.array([])
glon = np.array([])
glat = np.array([])
## Right now, HI_, etc. are emissivity TS values, which we do not care about.
for i in range( len(datCoor) ):
	## Coordinates ##
	glon = np.append(glon, float(datCoor[i].split(',')[0].split('(')[1]) )
	glat = np.append(glat, float(datCoor[i].split(',')[1].split(')')[0]) )

	## Distance ##
	if ( len(datDist[i].split('\t')) != 4 ):
		dist = np.append(dist, 0 )
	else:
		dista = float(datDist[i].split('\t')[-2])
		distb = float(datDist[i].split('\t')[-1])
		if ( distb == 0 ):
			dist = np.append(dist, dista )
		else:
			dist = np.append(dist, 0.5*(dista+distb) )

	## Emissivity ##
	a = datEmis[i].split('\t')
	if ( len(a) <= 1 ):
		HI = np.append( HI, 0 )
		HI_ = np.append(HI_,0 )
		CO = np.append( CO, 0 )
		CO_ = np.append(CO_,0 )
		DG = np.append( DG, 0 )
		DG_ = np.append(DG_,0 )
		continue
	else:
		aa = a[2].split(' ')
		for m in range(len(aa)):
			if ( ('HI' in aa[m]) and ('HI_far' not in aa[m]) ):
				HI = np.append( HI, float(aa[m+1][1:-2]) )
				continue
			elif ( ('CO' in aa[m]) ):
				CO = np.append( CO, float(aa[m+1][1:-2]) )
				continue
			elif ( 'E(B-V)' in aa[m] ):
				DG = np.append( DG, float(aa[m+1][1:-2]) )
				continue
		try:
			aa = a[3].split(' ')
			for m in range(len(aa)):
				if ( ('HI' in aa[m]) and ('HI_far' not in aa[m]) ):
					HI_ = np.append(HI_, float(aa[m+1][1:-2]) )
					continue
				elif ( ('CO' in aa[m]) ):
					CO_ = np.append(CO_, float(aa[m+1][1:-2]) )
					continue
				elif ( 'E(B-V)' in aa[m] ):
					DG_ = np.append(DG_, float(aa[m+1][1:-2]) )
					continue
		except:
			CO_ = np.append(CO_,0)
			HI_ = np.append(HI_,0)
			DG_ = np.append(DG_,0)

	
x = dist * np.cos( np.deg2rad(glon) ) * np.cos( np.deg2rad(glat) )
x = 8.33-(x/1000.)		## Because l=0 means towards Galactic center, means SMALLER R

