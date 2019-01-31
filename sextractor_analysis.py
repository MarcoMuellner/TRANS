import numpy as np
import matplotlib.pyplot as pl
from pandas import DataFrame
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from subprocess import call

# location of input images
image_r = "16_01_out/imager-0001R_150.fit" 
image_i = "16_01_out/imager-0001I_250.fit"
image_v = "16_01_out/imager-0001V_150.fit"
image_b = "16_01_out/imager-0001B_300.fit"

# call sextractor with custom config for each filter
# the "[0]" is needed as sextractor otherwise uses all extensions of the file and adds all detections to the catalogue!
call(['sex',image_r+"[0]",'-c','sex_config_R.txt'])
call(['sex',image_i+"[0]",'-c','sex_config_I.txt'])
call(['sex',image_v+"[0]",'-c','sex_config_V.txt'])
call(['sex',image_b+"[0]",'-c','sex_config_B.txt'])

# load sextractor output (apertures of detected sources are also available)
band_r = np.genfromtxt("sex_out_R.cat",comments="#",usecols=[0,2,3])
band_i = np.genfromtxt("sex_out_I.cat",comments="#",usecols=[0,2,3])
band_v = np.genfromtxt("sex_out_V.cat",comments="#",usecols=[0,2,3])
band_b = np.genfromtxt("sex_out_B.cat",comments="#",usecols=[0,2,3])

max_length = max(len(band_r),len(band_i),len(band_v),len(band_b))

def add_empty(band):
    if len(band) != max_length:
        diff = max_length - len(band)
        add = np.empty((diff,3))
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
        #print(f"{x_1},{y_1} - {x_2},{y_2} = {np.abs(x_1 - x_2)},{np.abs(y_1 - y_2)}")
        return m
    else:
        return None

res_dict = {"r": [],
        "i": [],
        "v": [],
        "b":[]}

for (m_r,x_r,y_r),(m_i,x_i,y_i),(m_v,x_v,y_v),(m_b,x_b,y_b) in zip(band_r,band_i,band_v,band_b):
    res_dict["r"].append(m_r)
    for band,n in zip([band_i,band_v,band_b],["i","v","b"]):
        value_set = False
        for m,x,y in band:
            if get_item(x_r,y_r,x,y,m) is not None:
                value_set = True
                res_dict[n].append(m)
                break

        if not value_set:
            res_dict[n].append(None)

df = DataFrame(data=res_dict)

fig : Figure = pl.figure(figsize=(16,4))
ax : Axes = fig.add_subplot(131)
ax.plot(df.b-df.v,df.v,'o',color='k',markersize=2)
ax.invert_yaxis()
ax.set_xlabel("B-V")
ax.set_ylabel("V")

ax : Axes = fig.add_subplot(132)
ax.plot(df.b-df.i,df.i,'o',color='k',markersize=2)
ax.invert_yaxis()
ax.set_xlabel("B-I")
ax.set_ylabel("I")

ax : Axes = fig.add_subplot(133)
ax.plot(df.b-df.r,df.r,'o',color='k',markersize=2)
ax.invert_yaxis()
ax.set_xlabel("B-R")
ax.set_ylabel("R")
pl.show()
