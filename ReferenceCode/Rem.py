
import os



m = os.listdir('../')
keep = ['bexpmap.fits','CMAP_big.fits','CMAP.fits','Config.txt','LC.fits', 'MODEL.xml', 'sel.fits', 'mkt.fits']

for i in range( len(m) ):
	if ( (m[i][0] == 'l') and (m[i] != 'l158b-33') ):
		try:
			os.chdir('../'+m[i]+'/ebin250_10000')
			dr = os.listdir('.')
			#for i in range( len(dr) ):
			#	if (dr[i] not in keep):
			#		try:
			#			os.remove(dr[i])
			#		except OSError:
			#			os.system('rm -r %s'%dr[i])
			os.remove('bexpmap.fits')
			os.chdir('/home/abrahams/HICO_survey/SourceSearch/Analysis')
		except:
			print "not a folder!"
