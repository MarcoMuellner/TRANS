# Set the WCS information manually by setting properties of the WCS
# object.

import numpy as np
from astropy import wcs
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle, Latitude, Longitude
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from astroquery.skyview import SkyView
from copy import deepcopy
from astropy.wcs import DistortionLookupTable
from astropy.table import Table



# Create a new WCS object.  The number of axes must be set
# from the start

hdulist = fits.open("16_01_out/imager-0001V_150.fit")

"""
490.2965 634.2495 40.61360560457 42.75804396956
583.4383 452.2284 024206.4057787989 +424556.787813359
456.8745 438.5215 024210.1035930438 +424511.284196129
"""
x_1,y_1 = 490.20833,634.63308
x_2,y_2 = 256.92361,489.85185
x_3,y_3 = 886.29707,458.55262
x_4,y_4 = 492.69707,574.95262

#x_2,y_2 = 493.76234,575.24076

alpha_1,delta_1 = 40.61360560457,42.75804396956

alpha_2,delta_2 = 40.55472965405, 42.69919233792
alpha_3,delta_3 = 40.56015662019,42.86645226731
alpha_4,delta_4 = 40.59228414797,42.76013370404

x,y,alpha,delta = [x_1,x_2,x_3,x_4],[y_1,y_2,y_3,y_4],[alpha_1,alpha_2,alpha_3,alpha_4],[delta_1,delta_2,delta_3,delta_4]

t = Table([x,y,alpha,delta],names=("x","y","alpha","delta"))



delta_alpha = (alpha_2 - alpha_1)/(x_2-x_1)
delta_delta = (delta_2 - delta_1)/(y_2-y_1)

print(f"{delta_alpha},{delta_delta}")


w = wcs.WCS(naxis = 2)
"""
w.wcs.crpix = [-234.75, 8.3393]
w.wcs.cdelt = np.array([-0.066667, 0.066667])
w.wcs.crval = [0, -90]
w.wcs.ctype = ["RA---AIR", "DEC--AIR"]
w.wcs.set_pv([(2, 1, 45.0)])
"""
coords = SkyCoord(40.61360560457,42.75804396956,unit=(u.deg,u.deg))

w.wcs.crpix = np.array((x_1,y_1))
w.wcs.crval = np.array((alpha_1,delta_1))
w.wcs.cdelt = [delta_alpha,delta_delta]
w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
w_plot = deepcopy(w)
w.wcs.crota = [0,90]
#w.wcs.set_pv([(0, 0, 0)])

#DistortionLookupTable(t,[x_1,y_1],[alpha_1,alpha_2],[delta_alpha,delta_delta])

pixcrd = np.array([[x_2, y_2]], np.float_)
pixcrd_2 = np.array([[x_3, y_3]], np.float_)
pixcrd_3 = np.array([[x_4, y_4]], np.float_)

world = w.all_pix2world(pixcrd,1)
print(f"{world},{alpha_2}:{delta_2})")
world = w.all_pix2world(pixcrd_2,1)
print(f"{world},{alpha_3}:{delta_3})")
world = w.all_pix2world(pixcrd_3,1)
print(f"{world},{alpha_4}:{delta_4})")

header = w.to_header()

hdulist[0].header = header

fig : Figure= plt.figure()
ax : Axes= fig.add_subplot(111, projection=w_plot)
ax.imshow(np.rot90(hdulist[0].data,3), origin='lower', cmap=plt.cm.viridis)
ax.coords.grid(color='white', alpha=0.5, linestyle='solid')
ax.coords['ra'].set_axislabel('Right Ascension')
ax.coords['dec'].set_axislabel('Declination')

#ax.invert_xaxis()
#plt.show()