import numpy as np
from astropy import units as u
from astropy.nddata import CCDData
import ccdproc
import matplotlib.pyplot as pl
import argparse
from directoryManager import cd
import os
from copy import deepcopy
from os import makedirs
from shutil import rmtree
from typing import List,Dict

parser = argparse.ArgumentParser()
parser.add_argument("folder", help="The folder containing the data", type=str)
parser.add_argument("output", help="Output folder", type=str)

str_flat = "Flat Field"
str_dark = "Dark Frame"
str_light = "Light Frame"
str_bias = 'Bias Frame'

args = parser.parse_args()

try:
    rmtree(args.output)
except:
    pass

try:
    makedirs(args.output)
except:
    pass


def sort_files(path: str):
    bias = []
    dark = []
    images = {}
    flat = {}

    with cd(args.folder):
        files = [f for f in os.listdir(".") if os.path.isfile(os.path.join(".", f))]

        for file in files:
            if file.endswith(".fit") or file.endswith(".fits"):
                file_data = CCDData.read(file, unit='electron')
                if str_flat == file_data.header['IMAGETYP']:
                    try:
                        flat[file_data.header['FILTER']].append(file_data)
                    except KeyError:
                        flat[file_data.header['FILTER']] = [file_data]
                elif str_dark == file_data.header['IMAGETYP']:
                    dark.append(file_data)
                elif str_light == file_data.header['IMAGETYP']:
                    images[file] = file_data
                elif str_bias == file_data.header['IMAGETYP']:
                    bias.append(file_data)

    return bias, dark, flat, images


def get_master_bias(bias: List[CCDData]) -> CCDData:
    master_bias = ccdproc.combine(bias, method='median')
    pl.figure(figsize=(8, 8))
    pl.title(f"Master bias")
    pl.imshow(master_bias, cmap='gray')
    pl.colorbar()
    pl.savefig(args.output + "master_bias.pdf")
    pl.close()

    return master_bias


def get_master_dark(dark: List[CCDData], master_bias : CCDData) -> CCDData:
    master_dark = ccdproc.combine(dark,method='median')
    master_dark = ccdproc.subtract_bias(master_dark,master_bias)
    master_dark.data = master_dark.data/master_dark.header["EXPOSURE"]
    master_dark.header["EXPOSURE"] = 1
    pl.figure(figsize=(8, 8))
    pl.title(f"Master dark")
    pl.imshow(master_dark,
              cmap='gray',vmax= np.median(master_dark.data) + 6*np.std(master_dark.data))
    pl.colorbar()
    pl.savefig(args.output + "master_dark.pdf")
    pl.close()
    return master_dark

def get_master_flat(flat : Dict[str,List[CCDData]],master_bias : CCDData,master_dark : CCDData)->Dict[str,CCDData]:
    master_flat = {}
    median_cutoff = 4000
    for key, value in flat.items():
        flat_list = []
        for data in value:
            print(f"Mean value {np.median(data.data)}")
            if np.median(data.data) < median_cutoff:
                continue
            data = ccdproc.subtract_bias(data, master_bias)
            data.data = data.data / np.amax(data.data)
            flat_list.append(data)
        combined_flat = ccdproc.combine(flat_list, method='median')
        median = np.median(combined_flat.data)
        std = np.std(combined_flat.data)
        flat_mask = np.logical_or(combined_flat.data < median - 3 * std, combined_flat.data > median + 3 * std)
        print(median)
        print(std)
        combined_flat.data[flat_mask] = median
        combined_flat.data = combined_flat.data / np.amax(combined_flat.data)
        pl.figure(figsize=(16, 16))
        pl.title(f"Master flat {key}")
        pl.imshow(combined_flat, cmap='gray')
        pl.colorbar()
        pl.savefig(args.output + f"master_flat{key}.pdf")
        pl.close()
        master_flat[key] = combined_flat

    return master_flat

bias, dark, flat, images = sort_files(args.folder)

master_bias = get_master_bias(bias)
master_dark = get_master_dark(dark,master_bias)


CCDData.write(master_bias, args.output + f"master_bias.fit")
CCDData.write(master_dark, args.output + f"master_dark.fit")

# master_dark.data[inv_mask] = master_dark.data[inv_mask] / master_dark.header["EXPOSURE"]


median_cutoff = 4000
master_flat = get_master_flat(flat,master_bias,master_dark)

for key, image in images.items():
    fig = pl.figure(figsize=(10, 7))
    ax = fig.add_subplot(121)
    ax.set_title(f"Image before")
    im = ax.imshow(image, cmap='gray', vmin=0)  # ,vmax=np.median(image.data)*2)
    pl.colorbar(im, ax=ax)
    image = ccdproc.subtract_bias(image, master_bias)
    CCDData.write(image, args.output + f"after_bias_{key}")
    image = ccdproc.subtract_dark(image, master_dark,
                                  dark_exposure=master_dark.header["EXPOSURE"] * u.second,
                                  data_exposure=image.header["EXPOSURE"] * 100*u.second,
                                  scale=True,
                                  exposure_unit=u.second)
    CCDData.write(image, args.output + f"after_dark{key}")
    image = ccdproc.flat_correct(image, master_flat[image.header["FILTER"]])
    CCDData.write(image, args.output + f"after_flat{key}")
    # image.data = (ccdproc.subtract_bias(image,master_bias).data - tmp_dark)/master_flat[image.header["FILTER"]]
    image.data[image.data < 0] = 0

    ax = fig.add_subplot(122)
    ax.set_title(f"Image after")
    im = ax.imshow(image, cmap='gray', vmin=0)  # ,vmax=np.median(image.data)*2)
    pl.colorbar(im, ax=ax)
    pl.savefig(args.output + f"{key}.pdf")
    pl.close()
    CCDData.write(image, args.output + f"{key}")
