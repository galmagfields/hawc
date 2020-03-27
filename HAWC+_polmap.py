import astropy.io.fits as fits
from aplpy import FITSFigure 
import astropy.units as u 
import os as os
import numpy as np 
import matplotlib.pyplot as plt 

def get_data(infile):
	#get data
	hawc = fits.open(infile)
	p    = hawc['DEBIASED PERCENT POL']
	perr = hawc['ERROR PERCENT POL']
	pa   = hawc['ROTATED POL ANGLE']
	stkI = hawc['STOKES I']
	stkIerr = hawc['ERROR I']
	pi   = hawc['DEBIASED POL FLUX']
	pierr   = hawc['ERROR POL FLUX']

	#Jy/px to Jy/sqarcsec
	pxscale = stkI.header['CDELT2']*3600
	stkI.data /= pxscale**2
	pi.data /= pxscale**2
	stkIerr.data /= pxscale**2

	#
	pi = fits.PrimaryHDU(pi.data,header=stkI.header)
	pierr = fits.PrimaryHDU(pierr.data,header=stkI.header)

	return p,perr,pa,stkI,stkIerr,pi,pierr,pxscale

def quality_cuts(stkI,stkIerr,p,perr,SNRp_cut,p_cut,SNRi_cut):
	#snr in P
	SNRp = p.data/perr.data
	mask_snrp = np.where(SNRp < SNRp_cut)
	p.data[mask_snrp] = np.nan
	#p_cut
	maskp = np.where(p.data > p_cut)
	p.data[maskp] = np.nan
	#snr in P
	SNRi = stkI.data/stkIerr.data
	mask_snri = np.where(SNRi < SNRi_cut)
	p.data[mask_snri] = np.nan
	return p


#### USER INPUT ####
## Files
DIR = 'Nominal/'
file = 'NGC1068_C_quarterbeam.fits'

# figure
title = 'NGC 1068 (89 $\mu$m; B-vectors)'
width  = 150
height = 150
vmax = 150
vmin = 0
figname = 'NGC1068_C.png'
cmap = 'jet'

title_size = 20
tick_labels = 15
label_plot = 15
label_colorbar = 15
tick_colorbar = 15
label_fontsize = 20

# pol map cuts
SNRp_cut = 3.0
p_cut = 15
SNRi_cut = 10
scalevec = 0.4 #1px = scalevec * 1% pol 
vec_legend = 5
#######


#### SCRIPT
fig = plt.figure(figsize=(13,10))

#Data
p,perr,pa,stkI,stkIerr,pi,pierr,pxscale = get_data(DIR+file)

stkI.data *= 1000.
stkIerr.data *= 1000.


#central coordinates
RA = (stkI.header['OBSRA']*u.hourangle).to(u.deg)
DEC = stkI.header['OBSDEC']*u.deg

###Figure
gc = FITSFigure(stkI,figure=fig)
#gc = FITSFigure(stkI,figure=fig)

##STOKES I
#colorscale
gc.show_colorscale(cmap=cmap,vmin=vmin,vmax=vmax,smooth=1,kernel='gauss')
#colorbae
gc.add_colorbar()
gc.colorbar.set_axis_label_text('I (mJy/sqarcsec)')
gc.colorbar.set_axis_label_font(size=label_colorbar)
gc.colorbar.set_font(size=tick_colorbar)
#recenter
gc.recenter(RA,DEC,width=width/3600,height=height/3600)
#contours
sigmaI = np.nanstd(stkIerr.data)
levelsI = sigmaI * 1.8**(np.arange(2,12,0.5))
#levelsI = np.arange(3,155,3)*sigmaI
gc.show_contour(colors='black',levels=levelsI,linewidth=0.1,\
				filled=False,smooth=1,kernel='box',alpha=0.4)
#gc.show_contour(levels=levelsI,\
#				filled=True,smooth=1,kernel='box',cmap='viridis')
#beam
gc.add_beam(color='red')
#title
gc.set_title(title,fontsize=title_size)


##Pol map
p = quality_cuts(stkI,stkIerr,p,perr,SNRp_cut,p_cut,SNRi_cut)

#2. polmap
gc.show_vectors(p,pa,scale=scalevec,\
				step=2,color='black',linewidth=4.0)
gc.show_vectors(p,pa,scale=scalevec,\
				step=2,color='grey',linewidth=2.0)

#show label
gc.add_label(0.1,0.95,'Total Flux',fontsize=label_fontsize,\
	color='white',weight='bold',relative=True)

#legend vector
vecscale = scalevec * pxscale/3600
gc.add_scalebar(vec_legend*vecscale,'P ='+np.str(vec_legend)+'%',\
				corner='bottom right',frame=True,color='black',facecolor='blue')
#Figure parameters
gc.tick_labels.set_font(size=tick_labels)
gc.axis_labels.set_font(size=label_plot)

fig.savefig(figname,dpi=300)
os.system('open '+figname)










