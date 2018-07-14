
import numpy as np
import os

glon = 5
glat = 38

os.chdir('../l%sb%s/ebin250_10000/'%(glon,glat))

try:
    f = open('LogLike.dat','r')
except:
    print "No likelihood data ... perform/redo fit"
    

dat = f.readlines()
f.close()

mdls = ['all','NoCODG','nodg','noco','Ptsrc','all_plus']
like = []

for i in np.arange( len(dat) ):
    if ( 'Final' in dat[i] ):
	line = dat[i].split('Final -log(likelihood) of ')[1].split(' = ')
	like.append( [line[0],float(line[1])] )

