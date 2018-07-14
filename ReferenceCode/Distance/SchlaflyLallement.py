
import numpy as np

f = open('Clds_Distance.txt','r')
g = open('CldsDist_lal.dat','r')

dats = f.readlines()
datl = g.readlines()

f.close()
g.close()

if ( len(dats) != len(datl) ):
	print "Lists not same length"
	

dist_schlafly = np.array([])
dist_serror = np.array([])
dist_lallement= np.array([])
dist_lerror = np.array([])
for i in range( len(dats) ):
	a = dats[i].split('\t')
	b = datl[i].split('     ')
	if ( len(a) == 3 ):
		dist_schlafly = np.append(dist_schlafly,0)
		dist_serror = np.append( dist_serror,0 )
	elif ( len(a) > 3 ):
		dist_schlafly = np.append(dist_schlafly,float(a[2]))
		dist_serror = np.append( dist_serror,float(a[-1]))

	dist_lallement = np.append(dist_lallement,float(b[3]))
	dist_lerror = np.append( dist_lerror, float( b[-1] ))


