
import aplpy
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size':16})
fig = plt.figure(figsize=(9,9))




im = []
labels = ['CO','E(B-V)','HI','IC','ISO','J0322.0+2335','J0245.4+2410','J0325.5+2223','J0213.0+2245']
# Slices = 15 means 2.5 GeV
## TOP ROW
im.append( aplpy.FITSFigure('srcmap_all_norm_0.001cut.fits',convention='calabretta',hdu=33,slices=[15],subplot=[0.1,0.65,0.28,0.28],figure=fig) )
im.append( aplpy.FITSFigure('srcmap_all_norm_0.001cut.fits',convention='calabretta',hdu=34,slices=[15],subplot=[0.375,0.65,0.28,0.28],figure=fig) )
im.append( aplpy.FITSFigure('srcmap_all_norm_0.001cut.fits',convention='calabretta',hdu=35,slices=[15],subplot=[0.65,0.65,0.28,0.28],figure=fig) )

## MIDDLE ROW
im.append( aplpy.FITSFigure('srcmap_all_norm_0.001cut.fits',convention='calabretta',hdu=37,slices=[15],subplot=[0.1,0.375,0.28,0.28],figure=fig) )
im.append( aplpy.FITSFigure('srcmap_all_norm_0.001cut.fits',convention='calabretta',hdu=38,slices=[15],subplot=[0.375,0.375,0.28,0.28],figure=fig) )
im.append( aplpy.FITSFigure('srcmap_all_norm_0.001cut.fits',convention='calabretta',hdu=29,slices=[15],subplot=[0.65,0.375,0.28,0.28],figure=fig) )

## BOTTOM ROW
im.append( aplpy.FITSFigure('srcmap_all_norm_0.001cut.fits',convention='calabretta',hdu=18,slices=[15],subplot=[0.1,0.1,0.28,0.28],figure=fig) )
im.append( aplpy.FITSFigure('srcmap_all_norm_0.001cut.fits',convention='calabretta',hdu=30,slices=[15],subplot=[0.375,0.1,0.28,0.28],figure=fig) )
im.append( aplpy.FITSFigure('srcmap_all_norm_0.001cut.fits',convention='calabretta',hdu=5,slices=[15],subplot=[0.65,0.1,0.28,0.28],figure=fig) )



for i in range(9):
	im[i].ticks.set_xspacing(10)
	im[i].ticks.set_yspacing(10)
	im[i].tick_labels.hide()
	im[i].show_colorscale(vmin=0.01,vmax=1,cmap=plt.cm.hot_r,stretch='log')
	im[i].show_colorbar()
	im[i].colorbar.hide()
	im[i].axis_labels.hide()
	if ( i in [2,4]):
		im[i].add_label(0.5,0.05,labels[i],relative=True,color='white',fontsize=14)
	else:
		im[i].add_label(0.5,0.05,labels[i],relative=True,color='black',fontsize=14)



im[0].tick_labels.show_y()
im[0].axis_labels.show_y()
im[0].tick_labels.set_yformat('dd.')
im[3].tick_labels.show_y()
im[3].axis_labels.show_y()
im[3].tick_labels.set_yformat('dd.')
im[6].tick_labels.show()
im[6].axis_labels.show()
im[6].tick_labels.set_xformat('dd.')
im[6].tick_labels.set_yformat('dd.')
im[7].tick_labels.show_x()
im[7].axis_labels.show_x()
im[7].tick_labels.set_xformat('dd.')
im[8].tick_labels.show_x()
im[8].axis_labels.show_x()
im[8].tick_labels.set_xformat('dd.')
im[8].colorbar.set_ticks([0.01,0.1,0.5,1])


fig.savefig('srcmaps_3x3.png')


