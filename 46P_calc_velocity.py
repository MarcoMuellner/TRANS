import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pandas import DataFrame
#from matplotlib.figure import Figure
#from matplotlib.axes import Axes
from subprocess import call

comet_R = "comet_exp_out/p46_johnr_t_90_001.fit" 
comet_I = "comet_exp_out/p46_johni_t_90_001.fit"
comet_V = "comet_exp_out/p46_johnv_t_90_001.fit"

call(['sex',comet_R+"[0]",'-c','sex_config_comet_R.txt'])
call(['sex',comet_I+"[0]",'-c','sex_config_comet_I.txt'])
call(['sex',comet_V+"[0]",'-c','sex_config_comet_V.txt'])

band_r = np.genfromtxt("sex_out_comet_R.cat",comments="#",usecols=[0,2,3])
band_i = np.genfromtxt("sex_out_comet_I.cat",comments="#",usecols=[0,2,3])
band_v = np.genfromtxt("sex_out_comet_V.cat",comments="#",usecols=[0,2,3])

print(band_r)
print(band_i)

delta_x = band_r[0][1]-band_i[0][1]
delta_y = band_r[0][2]-band_i[0][2]

deltaXY=[delta_x,delta_y]
print("delta x in px" , delta_x)
print("delta y in px" , delta_y)

fits_R = fits.open(comet_R)
fits_I = fits.open(comet_I)
header_R = fits_R[0].header
header_I = fits_I[0].header
print(header_R['DATE-OBS'],header_I['DATE-OBS']) # --> 4:11 minutes difference
delta_t = 251/60  # in minutes

velocity = [delta_x/delta_t, delta_y/delta_t]

print("velocity in pixels/minute", velocity)


# star 1 ID Gaia : Gaia DR2 37452806312572672
# star 1 pos Gaia : RA: 54.29164785983 DEC: +11.78326216157 	
# star 1 pos R: 615.6107, 367.22
# star 1 pos I: 615.5116, 368.0657
# star 2 ID Gaia : Gaia DR2 37449095460828288  
# star 2 pos Gaia: RA: 54.38409078240 DEC:+11.74309025637
# star 2 pos R: 490.2602, 718.3314
# star 2 pos I: 490.1379, 719.535

x_1_R = 615.6107
y_1_R = 367.22

x_2_R = 490.2602
y_2_R = 718.3314

x_1_I = 615.5116
y_1_I = 368.0657

x_2_I = 490.1379
y_2_I = 719.535

x_1_G = 54.29164785983
y_1_G = 11.78326216157

x_2_G = 54.38409078240
y_2_G = 11.74309025637

deltaG = [x_1_G-x_2_G,y_1_G-y_2_G]
deltaR = [x_1_R-x_2_R,y_1_R-y_2_R]
deltaI = [x_1_I-x_2_I,y_1_I-y_2_I]
shiftRI1 = [x_1_R-x_1_I,y_1_R-y_1_I]
shiftRI2 = [x_2_R-x_2_I,y_2_R-y_2_I]

print(deltaG, deltaR,deltaI,shiftRI1,shiftRI2)

def vecabs(deltalist):
	return(np.sqrt(deltalist[0]**2+deltalist[1]**2))

print(vecabs(deltaG),vecabs(deltaR),vecabs(deltaI),vecabs(deltaXY),vecabs(shiftRI1),vecabs(shiftRI2))

print(3600*vecabs(velocity)*vecabs(deltaG)/np.mean([vecabs(deltaR),vecabs(deltaI)]))
