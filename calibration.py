import numpy as np
from astropy import units as u
from astropy.nddata import CCDData
import ccdproc
import matplotlib.pyplot as pl
import argparse
from directoryManager import cd
import os
from astropy.wcs import WCS
from copy import deepcopy
from os import makedirs
from shutil import rmtree
from typing import List, Dict
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from collections import OrderedDict

imstats = lambda dat: (dat.min(), dat.max(), dat.mean(), dat.std())

parser = argparse.ArgumentParser()
parser.add_argument("folder", help="The folder containing the data", type=str)

str_dark = "Dark Frame"
str_dark_master = "Dark Frame Master"
str_bias = 'Bias Frame'
str_bias_master = 'Bias Frame Master'

args = parser.parse_args()
in_path = args.folder
out_path = args.folder + "_out"

try:
    rmtree(out_path)
except:
    pass

try:
    makedirs(out_path)
except:
    pass


def sort_files(path: str):
    result_files = {}

    with cd(path):
        files = [f for f in os.listdir(".") if os.path.isfile(os.path.join(".", f))]

        for file in files:
            if file.endswith(".fit") or file.endswith(".fits"):
                file_data = CCDData.read(file, unit='electron')
                temp = file_data.header["SET-TEMP"]

                if temp not in result_files.keys():
                    result_files[temp] = {str_bias: [], str_dark: []}

                if str_dark == file_data.header['IMAGETYP']:
                    result_files[temp][str_dark].append(file_data)
                elif str_bias == file_data.header['IMAGETYP']:
                    result_files[temp][str_bias].append(file_data)

    return result_files


def overscan_trim_and_sigma_clip_median(image_list):
    """
    Combine a list of images using median

    This function does several steps:

    1. Subtract overscan
    2. Trim image
    3. sigma clip image using a median of the unclipped stack as the baseline
    4. combine the images on the list using median

    ** It modifies the images in the input list. **
    """
    combo = ccdproc.Combiner(image_list)
    combo.sigma_clipping(func=np.ma.median)
    return combo


def get_master_bias(bias: List[CCDData]) -> CCDData:
    bias_combine = ccdproc.Combiner(bias)
    master_bias = bias_combine.average_combine()
    #master_bias = ccdproc.combine(bias, method='median')
    master_bias.header = bias[0].header

    return master_bias


def get_master_dark(dark: List[CCDData], master_bias: CCDData) -> CCDData:
    # make a combiner for sigma clipping and median combine
    for i in range(0, len(dark)):
        dark[i] = ccdproc.subtract_bias(dark[i], master_bias)
    a_combiner = overscan_trim_and_sigma_clip_median(dark)
    master_dark = a_combiner.median_combine(median_func=np.ma.median)
    master_dark.data[master_dark.data < 0] = 0
    master_dark.data = master_dark.data / dark[0].header["EXPOSURE"]
    master_dark.header = dark[0].header
    master_dark.header["EXPOSURE"] = 1
    return master_dark


res_files = sort_files(in_path)

for temp, values in res_files.items():
    master_bias = get_master_bias(values[str_bias])
    master_dark = get_master_dark(values[str_dark], master_bias)

    CCDData.write(master_bias, out_path + f"/master_bias_m_{temp}.fit")
    CCDData.write(master_dark, out_path + f"/master_dark_m_{temp}.fit")

    res_files[temp][str_bias_master] = master_bias
    res_files[temp][str_dark_master] = master_dark


res_files = OrderedDict(reversed(sorted(res_files.items())))
n = len(res_files)

y = 3
x = n - y - 1
fig_bias: Figure = pl.figure(figsize=(16, 10))
fig_bias_hist: Figure = pl.figure(figsize=(16, 10))
fig_dark: Figure = pl.figure(figsize=(16, 10))
fig_dark_hist: Figure = pl.figure(figsize=(16, 10))
fig_dark_progression :Figure = pl.figure(figsize=(16,10))
ax_progression : Axes = fig_dark_progression.add_subplot(111)
ax_progression.set_xlabel("Temperature")
ax_progression.set_ylabel("Mean dark current")
ax_progression.invert_xaxis()

for (temp, values), i in zip(res_files.items(), range(1, n + 1)):

    for im_type in [str_bias_master, str_dark_master]:
        d_min, d_max, d_mean, d_std = imstats(np.asarray(values[im_type]))
        d_median = np.median(values[im_type])
        flat = values[im_type].data.flatten()

        if im_type == str_bias_master:
            ax: Axes = fig_bias.add_subplot(x, y, i)
            ax_hist : Axes = fig_bias_hist.add_subplot(x, y, i)
            im = ax.imshow(values[im_type], vmax=d_mean + 4 * d_std, vmin=d_mean - 4 * d_std, cmap='gray')
            range_hist = (d_mean - 4 * d_std, d_mean + 4 * d_std)
        else:
            ax: Axes = fig_dark.add_subplot(x, y, i)
            ax_hist: Axes = fig_dark_hist.add_subplot(x, y, i)
            im = ax.imshow(values[im_type], vmax=d_median + 0.5 * d_std, vmin=0, cmap='gray')
            range_hist = (0, d_mean + 0.5 * d_std)
            #ax_progression.plot(temp, d_mean,'x', label=f'{temp}')
            ax_progression.errorbar(temp,d_mean,d_std/10,fmt='x',label=f"{temp}")
        ax_hist.hist(flat, bins=30, density=True, range=range_hist, histtype='step',color='k')

        pl.colorbar(im, ax=ax)
        ax.set_title(f"{im_type} {values[im_type].header['SET-TEMP']}")
        ax_hist.set_title(f"{im_type} {values[im_type].header['SET-TEMP']}")

ax_progression.legend()
fig_dark_progression.savefig(out_path+f"/dark_mean_progression.pdf")
fig_dark.savefig(out_path+f"/dark.pdf")
fig_bias.savefig(out_path+f"/bias.pdf")
fig_bias_hist.savefig(out_path+f"/bias_hist.pdf")
fig_dark_hist.savefig(out_path+f"/dark_hist.pdf")
pl.close()
#pl.show()

# print(res_files)
