import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pandas import DataFrame
#from matplotlib.figure import Figure
#from matplotlib.axes import Axes
from subprocess import call
import astropy.units as u

comet_R = "astrometry_out/wcs_46P_01_R_090.fits" 
comet_I = "astrometry_out/wcs_46P_01_I_090.fits"
comet_V = "astrometry_out/wcs_46P_01_V_090.fits"

comet_R_2 = "astrometry_out/wcs_46P_02_R_090.fits" 
comet_I_2 = "astrometry_out/wcs_46P_02_I_090.fits"

call(['sex',comet_R+"[0]",'-c','sex_config_comet_R.txt'])
call(['sex',comet_I+"[0]",'-c','sex_config_comet_I.txt'])
call(['sex',comet_V+"[0]",'-c','sex_config_comet_V.txt'])


band_r = np.genfromtxt("sex_out_comet_R.cat",comments="#",usecols=[2,4,5,8,9])
band_i = np.genfromtxt("sex_out_comet_I.cat",comments="#",usecols=[2,4,5,8,9])
band_v = np.genfromtxt("sex_out_comet_V.cat",comments="#",usecols=[2,4,5,8,9])

ra_v = band_v[0][3] *u.degree
dec_v = band_v[0][4] *u.degree
ra_r = band_r[0][3] *u.degree
dec_r = band_r[0][4] *u.degree
ra_i = band_i[0][3] *u.degree
dec_i = band_i[0][4] *u.degree

fits_R = fits.open(comet_R)
fits_I = fits.open(comet_I)
fits_V = fits.open(comet_V)

obsT_R = fits_R[0].header['JD-HELIO'] *u.day -2458466 *u.day
obsT_I = fits_I[0].header['JD-HELIO'] *u.day -2458466 *u.day
obsT_V = fits_V[0].header['JD-HELIO'] *u.day -2458466 *u.day

call(['sex',comet_R_2+"[0]",'-c','sex_config_comet_R.txt'])
call(['sex',comet_I_2+"[0]",'-c','sex_config_comet_I.txt'])

band_r_2 = np.genfromtxt("sex_out_comet_R.cat",comments="#",usecols=[2,4,5,8,9])
band_i_2 = np.genfromtxt("sex_out_comet_I.cat",comments="#",usecols=[2,4,5,8,9])

ra_r2 = band_r_2[0][3] *u.degree
dec_r2 = band_r_2[0][4] *u.degree
ra_i2 = band_i_2[0][3] *u.degree
dec_i2 = band_i_2[0][4] *u.degree

fits_R_2 = fits.open(comet_R_2)
fits_I_2 = fits.open(comet_I_2)

obsT_R2 = fits_R_2[0].header['JD-HELIO'] *u.day -2458466 *u.day
obsT_I2 = fits_I_2[0].header['JD-HELIO'] *u.day -2458466 *u.day

print(obsT_V,obsT_R,obsT_R2,obsT_I,obsT_I2)
# ORDER : VR, RR2, R2I, II2

delT_VR = obsT_V - obsT_R
delT_RR2 = obsT_R - obsT_R2
delT_R2I = obsT_R2 - obsT_I
delT_II2 = obsT_I - obsT_I2

print(delT_VR,delT_RR2,delT_R2I,delT_II2)

pmra_VR = ((ra_v - ra_r)/delT_VR).to(u.arcsec/u.minute)
pmdec_VR = ((dec_v - dec_r)/delT_VR).to(u.arcsec/u.minute)
pmabs_VR = np.sqrt(pmra_VR**2*np.cos(np.mean([dec_v.value,dec_r.value]))**2+pmdec_VR**2)

pmra_RR2 = ((ra_r - ra_r2)/delT_RR2).to(u.arcsec/u.minute)
pmdec_RR2 = ((dec_r - dec_r2)/delT_RR2).to(u.arcsec/u.minute)
pmabs_RR2 = np.sqrt(pmra_RR2**2*np.cos(np.mean([dec_r.value,dec_r2.value]))**2+pmdec_RR2**2)

pmra_R2I = ((ra_r2 - ra_i)/delT_R2I).to(u.arcsec/u.minute)
pmdec_R2I = ((dec_r2 - dec_i)/delT_R2I).to(u.arcsec/u.minute)
pmabs_R2I = np.sqrt(pmra_R2I**2*np.cos(np.mean([dec_r2.value,dec_i.value]))**2+pmdec_R2I**2)

pmra_II2 = ((ra_i - ra_i2)/delT_II2).to(u.arcsec/u.minute)
pmdec_II2 = ((dec_i - dec_i2)/delT_II2).to(u.arcsec/u.minute)
pmabs_II2 = np.sqrt(pmra_II2**2*np.cos(np.mean([dec_i.value,dec_i2.value]))**2+pmdec_II2**2)

print(pmra_VR,pmdec_VR,pmabs_VR)
print(pmra_RR2,pmdec_RR2,pmabs_RR2)
print(pmra_R2I,pmdec_R2I,pmabs_R2I)
print(pmra_II2,pmdec_II2,pmabs_II2)
print(np.mean([pmabs_VR.value,pmabs_RR2.value,pmabs_R2I.value,pmabs_II2.value])*u.arcsec/u.minute,np.std([pmabs_VR.value,pmabs_RR2.value,pmabs_R2I.value,pmabs_II2.value])*u.arcsec/u.minute)

# RESULT: 2.232657404811999 arcsec / min +/-  1.3051232974035107 arcsec / min

