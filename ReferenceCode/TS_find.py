
import numpy as np
import os


############################################################################
############################################################################
def FindTS(num,glon,glat,home_dir):
	try:
		f = open('LogLike.dat','r')
	except:
		print "No likelihood data ... perform/redo fit"
		return "Nothing"
    

	dat = f.readlines()
	f.close()
	
	mdls = ['all','nocodg','nodg','noco','nocld']# haven't yet done ptsrc, all_plus
	NewMin = []
	MDL = []
	
	for i in np.arange(len(dat)):
		if ( len(dat[i].split('Final ')) > 1 ):
			NewMin.append( float(dat[i].split('Final ')[1].split(' = ')[1]) )
			MDL.append( dat[i].split('Final ')[1].split(' ')[2] )
		elif ('ERROR IN NEWMINUIT FIT' in dat[i]):
			NewMin.append( float(dat[i-1].split(' = ')[1]) )
    			MDL.append( dat[i-1].split(' ')[2] )
	
	#for i in np.arange( len(dat) ):
	#    if ( 'Final' in dat[i] ):
	#	
	
	TS_H2 = 2*(NewMin[1] - NewMin[0])
	
	ind = []
	for i in np.arange( len(mdls) ):
		for j in np.arange( len(MDL) ):
			if (MDL[j].lower() == mdls[i].lower()):
				ind.append( j )
	
	TS = []
	for i in np.arange( len(ind) ):
		TS.append( 2*(NewMin[i] - NewMin[0]))

	# Result: TS = [0.0, TS_H2, TS_DG, TS_CO, TS_PS, TS_AG]
	# Result: TS = [0.0, TS_H2, TS_GC, TS_CO, TS_CLD]
	#print "TS_CLD= %s"%TS[4]
	MDs = ['TS_H2', 'TS_DG', 'TS_CO', 'TS_CLD']
	aa = open(home_dir + '/TS.dat','a')
	aa.write('Cloud %s at \t (%s,%s)\n'%(num,glon,glat))
	aa.write("%s = %s\n"%(MDs[0],TS[1]))
	aa.write("%s = %s\n"%(MDs[1],TS[2]))
	aa.write("%s = %s\n"%(MDs[2],TS[3]))
	aa.write("%s= %s\n\n\n"%(MDs[3],TS[4]))
	aa.close()
	return 0



############################################################################

def CoordRead(line):
	f = open('CldNum.txt','r')
	dat = f.readlines()
	f.close()
	###   line = int(sys.argv[1])
	lon = float( dat[line].split(',')[0].split('(')[1] )
	lat = float( dat[line].split(',')[1].split(')')[0] )
	#print "Great! Running cloud %s, at (l,b)=(%s,%s)"%(line-1,lon,lat)
	return (lon,lat)

############################################################################
############################################################################

home_dir = os.getcwd()
print home_dir
print home_dir+'TS.dat'

for i in range( 37,106,1 ):
	(glon,glat) = CoordRead(i)
	os.chdir('../l%sb%s/ebin250_10000/'%(int(glon),int(glat)))
	try:
		FindTS(i,glon,glat,home_dir)
	except:
		print "Nope, missing something at cloud %s"%i

	os.chdir(home_dir)

print "Done calculating TS values"



