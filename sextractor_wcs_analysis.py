import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from subprocess import call
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

# location of input images

# M34 #

#image_r = "astrometry_out/wcs_M34_01_R_150.fits" 
#image_i = "astrometry_out/wcs_M34_01_I_250.fits"
#image_v = "astrometry_out/wcs_M34_01_V_150.fits"
#image_b = "astrometry_out/wcs_M34_01_B_300.fits"

# C0001+557  # does not work as there is one star (index 11) without Gaia entry

#image_r = "astrometry_out/wcs_C0001+557_02_R_200.fits" #01 doesn't work for some reason
#image_i = "astrometry_out/wcs_C0001+557_01_I_100.fits"
#image_v = "astrometry_out/wcs_C0001+557_01_V_200.fits"
#image_b = "astrometry_out/wcs_C0001+557_01_B_500.fits"

# Teutsch55

image_r = "astrometry_out/wcs_Teutsch55_03_R_400.fits" 
image_i = "astrometry_out/wcs_Teutsch55_03_I_400.fits"
image_v = "astrometry_out/wcs_Teutsch55_03_V_400.fits"
image_b = "astrometry_out/wcs_Teutsch55_02_B_300.fits"


# call sextractor with custom config for each filter
# the "[0]" is needed as sextractor otherwise uses all extensions of the file and adds all detections to the catalogue!

call(['sex',image_r+"[0]",'-c','sex_config_R.txt'])
call(['sex',image_i+"[0]",'-c','sex_config_I.txt'])
call(['sex',image_v+"[0]",'-c','sex_config_V.txt'])
call(['sex',image_b+"[0]",'-c','sex_config_B.txt'])

# load sextractor output (apertures of detected sources are also available)
# column 2 is mag, 4 is x in pixel, 5 is y in pixel, 8,9 are RA/DEC see custom.param and .cat files

band_r = np.genfromtxt("sex_out_R.cat",comments="#",usecols=[2,4,5,8,9])
band_i = np.genfromtxt("sex_out_I.cat",comments="#",usecols=[2,4,5,8,9])
band_v = np.genfromtxt("sex_out_V.cat",comments="#",usecols=[2,4,5,8,9])
band_b = np.genfromtxt("sex_out_B.cat",comments="#",usecols=[2,4,5,8,9])

max_length = max(len(band_r),len(band_i),len(band_v),len(band_b))

def add_empty(band):
	if len(band) != max_length:
		diff = max_length - len(band)
		add = np.empty((diff,5))
		add[:] = np.nan
		band = np.vstack((band,add))
	return band

band_r = add_empty(band_r)
band_i = add_empty(band_i)
band_v = add_empty(band_v)
band_b = add_empty(band_b)

sigma = 3

def get_item(x_1,y_1,x_2,y_2,m):
	if np.abs(x_1 - x_2) < sigma and np.abs(y_1 - y_2) < sigma:
		return m
	else:
		return None

res_dict = {"r": {"mag": [],"x": [],"y": [],"RA": [],"DEC": []}, "i": {"mag": [],"x": [],"y": [],"RA": [],"DEC": []}, "v": {"mag": [],"x": [],"y": [],"RA": [],"DEC": []}, "b": {"mag": [],"x": [],"y": [],"RA": [],"DEC": []}}

for (m_r,x_r,y_r,r_r,d_r),(m_i,x_i,y_i,r_i,d_i),(m_v,x_v,y_v,r_v,d_v),(m_b,x_b,y_b,r_b,d_b) in zip(band_r,band_i,band_v,band_b):
	res_dict["r"]["mag"].append(m_r)
	res_dict["r"]["x"].append(x_r)
	res_dict["r"]["y"].append(y_r)
	res_dict["r"]["RA"].append(r_r)
	res_dict["r"]["DEC"].append(d_r)
	for band,n in zip([band_i,band_v,band_b],["i","v","b"]):
		value_set = False
		for m,x,y,r,d in band:
			if get_item(x_r,y_r,x,y,m) is not None:
				res_dict[n]['mag'].append(m)
				res_dict[n]['x'].append(x)
				res_dict[n]['y'].append(y)
				res_dict[n]['RA'].append(r)
				res_dict[n]['DEC'].append(d)
				value_set = True
				break
		if not value_set:
			res_dict[n]['mag'].append(None)
			res_dict[n]['x'].append(None)
			res_dict[n]['y'].append(None)
			res_dict[n]['RA'].append(None)
			res_dict[n]['DEC'].append(None)

#print(res_dict['v']['x'])

final_dic= {'r':{'mag':[],'x':[],'y':[],'RA':[],'DEC':[]},'i':{'mag':[],'x':[],'y':[],'RA':[],'DEC':[]},'v':{'mag':[],'x':[],'y':[],'RA':[],'DEC':[]},'b':{'mag':[],'x':[],'y':[],'RA':[],'DEC':[]}}

for ind in range(len(res_dict['r']['mag'])):
	if res_dict['i']['mag'][ind] != None and res_dict['v']['mag'][ind] != None and res_dict['b']['mag'][ind] != None:
		final_dic['r']['mag'].append(res_dict['r']['mag'][ind])
		final_dic['r']['x'].append(res_dict['r']['x'][ind])
		final_dic['r']['y'].append(res_dict['r']['y'][ind])
		final_dic['r']['RA'].append(res_dict['r']['RA'][ind])
		final_dic['r']['DEC'].append(res_dict['r']['DEC'][ind])
		final_dic['i']['mag'].append(res_dict['i']['mag'][ind])
		final_dic['i']['x'].append(res_dict['i']['x'][ind])
		final_dic['i']['y'].append(res_dict['i']['y'][ind])
		final_dic['i']['RA'].append(res_dict['i']['RA'][ind])
		final_dic['i']['DEC'].append(res_dict['i']['DEC'][ind])
		final_dic['v']['mag'].append(res_dict['v']['mag'][ind])
		final_dic['v']['x'].append(res_dict['v']['x'][ind])
		final_dic['v']['y'].append(res_dict['v']['y'][ind])
		final_dic['v']['RA'].append(res_dict['v']['RA'][ind])
		final_dic['v']['DEC'].append(res_dict['v']['DEC'][ind])
		final_dic['b']['mag'].append(res_dict['b']['mag'][ind])
		final_dic['b']['x'].append(res_dict['b']['x'][ind])
		final_dic['b']['y'].append(res_dict['b']['y'][ind])
		final_dic['b']['RA'].append(res_dict['b']['RA'][ind])
		final_dic['b']['DEC'].append(res_dict['b']['DEC'][ind])

print("max sources:",len(res_dict['r']['mag']))
print("matching sources:",len(final_dic['r']['mag']))

df = DataFrame(data=final_dic)
#print(df.r.mag) # very nice

fig : Figure = plt.figure(figsize=(16,4))
ax : Axes = fig.add_subplot(131)
ax.plot(np.asarray(df.b.mag)-np.asarray(df.v.mag),df.v.mag,'o',color='k',markersize=2)
ax.invert_yaxis()
ax.set_xlabel("B-V")
ax.set_ylabel("V")

ax : Axes = fig.add_subplot(132)
ax.plot(np.asarray(df.b.mag)-np.asarray(df.i.mag),df.i.mag,'o',color='k',markersize=2)
ax.invert_yaxis()
ax.set_xlabel("B-I")
ax.set_ylabel("I")

ax : Axes = fig.add_subplot(133)
ax.plot(np.asarray(df.b.mag)-np.asarray(df.r.mag),df.r.mag,'o',color='k',markersize=2)
ax.invert_yaxis()
ax.set_xlabel("B-R")
ax.set_ylabel("R")
#pl.savefig("HR_Teutsch55_22.pdf",bbox_inches='tight')
plt.show()

Gaia_dic = {"g":{"mag":[],"RA":[],"DEC":[],"pmra":[],"pmra_err":[],"pmdec":[],"pmdec_err":[],"parallax":[],"parallax_err":[],"rv":[],"rv_err":[],"dist":[],"x":[],"y":[],"ID":[]}}

for jnd in range(len(final_dic['r']['RA'])):
	coord = SkyCoord(ra=final_dic['r']['RA'][jnd], dec=final_dic['r']['DEC'][jnd], unit=(u.degree, u.degree), frame='icrs')
	width = u.Quantity(0.001, u.deg)
	height = u.Quantity(0.001, u.deg)
	g = Gaia.query_object_async(coordinate=coord, width=width, height=height)
	#print(g)	
	#print(g["phot_g_mean_mag"])
	#print(jnd)	
	if g["dist"].quantity.value == 0 or not g["dist"].quantity.value:
		print("FAILURE TO FIND GAIA MATCH FOR INDEX",jnd)
	else:	
		Gaia_dic['g']['mag'].append(g["phot_g_mean_mag"].quantity.value[0])
		Gaia_dic['g']['RA'].append(g["ra"].quantity.value[0])
		Gaia_dic['g']['DEC'].append(g["dec"].quantity.value[0])
		Gaia_dic['g']['pmra'].append(g["pmra"].quantity.value[0])
		Gaia_dic['g']['pmra_err'].append(g["pmra_error"].quantity.value[0])
		Gaia_dic['g']['pmdec'].append(g["pmdec"].quantity.value[0])
		Gaia_dic['g']['pmdec_err'].append(g["pmdec_error"].quantity.value[0])
		Gaia_dic['g']['parallax'].append(g["parallax"].quantity.value[0])
		Gaia_dic['g']['parallax_err'].append(g["parallax_error"].quantity.value[0])
		Gaia_dic['g']['rv'].append(g["radial_velocity"].quantity.value[0])
		Gaia_dic['g']['rv_err'].append(g["radial_velocity_error"].quantity.value[0])
		Gaia_dic['g']['dist'].append(g["dist"].quantity.value[0])
		Gaia_dic['g']['x'].append(final_dic['r']['x'][jnd])
		Gaia_dic['g']['y'].append(final_dic['r']['y'][jnd])
		Gaia_dic['g']['ID'].append(g["source_id"].quantity.value[0])
#print(Gaia_dic)

GaiaDF = DataFrame(data=Gaia_dic)
print(GaiaDF)
#GaiaDF.to_csv('gaiatest.csv')
