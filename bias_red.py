import numpy as np
from astropy import units as u
from astropy.nddata import CCDData
import ccdproc
import matplotlib.pyplot as pl

bias_list = [CCDData.read(f"data/imager-00{'{:02d}'.format(i)}bias.fit",unit='electron') for i in range(1,21)]
dark_list = [CCDData.read(f"data/imager-00{'{:02d}'.format(i)}dark_final.fit",unit='electron') for i in range(1,11)]

master_bias = ccdproc.combine(bias_list,method='median')
master_bias_average = ccdproc.combine(bias_list,method='average')

pl.figure(figsize=(16, 16))
pl.title(f"Master bias")
pl.imshow(master_bias,cmap='gray')
pl.colorbar()

pl.figure(figsize=(16, 16))
pl.title(f"Master bias average")
pl.imshow(master_bias_average,cmap='gray')
pl.colorbar()


dark_list = [CCDData.read(f"data/imager-00{'{:02d}'.format(i)}dark_final.fit",unit='electron') for i in range(1,11)]
master_dark = ccdproc.combine(dark_list,method='median')
#master_dark.data = np.clip(master_dark, 0, 1500)
master_dark_average = ccdproc.combine(dark_list,method='average')
#master_dark_average.data = np.clip(master_dark_average, 0, 1500)
dark_without_bias = ccdproc.subtract_bias(master_dark,master_bias)

master_dark.data = np.clip(master_dark,0,1500)
dark_without_bias.data = np.clip(dark_without_bias,0,1500)

pl.figure(figsize=(16, 16))
pl.title(f"Master dark")
pl.imshow(master_dark,cmap='gray')
pl.colorbar()

pl.figure(figsize=(16, 16))
pl.title(f"Master dark average")
pl.imshow(master_dark_average,cmap='gray')
pl.colorbar()


pl.figure(figsize=(16, 16))
pl.title(f"Master dark without bias")
pl.imshow(dark_without_bias,cmap='gray')
pl.colorbar()
pl.show()