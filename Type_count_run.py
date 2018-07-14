
from Type_count import Cloud
import numpy as np
a = open('Clds_TS.txt','a')

TS = ['TS_H2', 'TS_CO', 'TS_DG', 'TS_CLD', 'TS_EX','TS_AGN']
Comp = ['HI','HI_far','CO','E(B-V)','IC','iso_P8R2_SOURCE_V6_v06']
def spaces(x):
	""" Writing coordinates in a human-readable format, without \t.
		Not necessary for latex output.
	"""
	sp = ''
	if np.abs(x)<10:
		sp += "  "
	elif ( np.abs(x)<100 ):
		sp += " "
	
	if x>0:
		sp += " "

	return sp



### Let's write a table in latex. This will just be the "data" part of the table . . . maybe.
# Table: cloud properties. Need:
#	 coordinates (l,b)
#	 Peak WCO, NHI, E(B-V)res
#	 Distance
#	 Area??

#A = open('CldProperty_table.txt','a')
#table = ''
#table += 'l & b & Peak \hi & Peak \WCO & Peak \Avres & Distance \\\\ \n'
#table += 'deg & deg & $10^{21}$ \cmtwo & K \kms & mag & pc \\\\ \hline\n'
#for i in xrange( 93 ):
#	try:
#		cld = Cloud(i)
#		cld.getTs()
#		cld.peakComponents()
#		cld.getDistance(mode='lallement')
#	except:
#		print '%.1f, %.1f'%(cld.glon,cld.glat)
#		continue
#
#	
#	table += '%.1f & %.1f & $%.2f^{%.2f}_{%.2f}$ & $%.2f \pm %.2f$ & $%.2f \pm %.2f$ & $%d \pm %d$ \\\\ \n'%(cld.glon%360, cld.glat, cld.peak[0]/1e21, (0.07*cld.peak[0])/1e21, (0.15*cld.peak[0])/1e21, cld.peak[1], cld.peakError[1], cld.peak[2],  cld.peakError[2], cld.dist,cld.distError)
#
#	try:
#		cld.getDistance(
#
#
#table+='\hline'
#
#A.write(table)
#A.close()


#c = open('CldDistance.txt','a')
a = open('Clds_TS_table.txt','w')
for i in np.arange( 93 ):
	try:
		cld = Cloud(i)
		cld.getTs()
		cld.getEmissivity()
#		cld.getDistance()
	except:
		a.write('%s\n'%i)
		continue
#
#	## Write TS ##
#	a_ = spaces(cld.glon)+'%s,'%int(cld.glon)+spaces(cld.glat)+ '%s\t'%(int(cld.glat))
	a_ = '%.1f & %.1f & '%(cld.glon%360,cld.glat)
	for ts_i in TS:
		try:
#			a_ += '%s\t\t'%round(cld.ts[ts_i],2)
			a_ += '%d'%int( cld.ts[ts_i] )
			a_ += ' & '
		except:
			a_ += '    \t\t'

#	a_ += '\n'
	a_ += '\\\\\n'
	a.write(a_)
#
#	## Write Emissivity ##
#	b_ = '%s\t\t'%i
#	for component in Comp:
#		try:
#			b_ += '%s\t%s\t\t'%(cld.emissivity[component],cld.emissivityTs[component])
#		except:
#			b_ += '        \t    \t\t'
#	try:
#		b_ += '%s\t%s'%(cld.emissivity,cld.emissivityTs)
#		b_ += '%s\t%s\t%s'%(cld.emissivity,cld.emissivityError,cld.emissivityScale)
#	except:
#		pass
#
#	b_ += '\n'
#	b.write(b_)

	## Write distance ##
#	c_ = '%s (%s)\t\t'%(i,cld._numClds)
#	try:
#		c_ += '%s\t%s\n'%(round(np.mean(cld.dist)),round(cld.distError))
#	except:
#		c_ += '\n'
#
#	c.write(c_)


a.close()
#b.close()
#c.close()
###### 
## Next time: run this loop, storing EACH thing in own array:
## HI = [cld1, cld2, ... ]
## HI_TS = ...
## TS_H2 = [cld1, cld2, ... ]
## ...
