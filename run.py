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


imstats = lambda dat: (dat.min(), dat.max(), dat.mean(), dat.std())

def oscan_and_trim(image_list):
    """
    Remove overscan and trim a list of images. The original list is replaced by a list of images
    with the changes applied.
    """
    for idx, img in enumerate(image_list):
        oscan = ccdproc.subtract_overscan(img, img[:, 3075:3079], add_keyword={'oscan_sub': True, 'calstat': 'O'}, model=models.Polynomial1D(1))
        image_list[idx] = ccdproc.trim_image(oscan[:, :3073], add_keyword={'trimmed': True, 'calstat': 'OT'})

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
    bias_combine = ccdproc.Combiner(bias)
    master_bias = bias_combine.average_combine()
    #master_bias = ccdproc.combine(bias, method='median')
    pl.figure(figsize=(8, 8))
    pl.title(f"Master bias")
    bias_min, bias_max, bias_mean, bias_std = imstats(np.asarray(master_bias))
    pl.imshow(master_bias, vmax=bias_mean + 4 * bias_std, vmin=bias_mean - 4 * bias_std,cmap='gray')
    pl.colorbar()
    pl.savefig(args.output + "master_bias.pdf")
    pl.close()

    return master_bias


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


def get_master_dark(dark: List[CCDData], master_bias : CCDData) -> CCDData:
    # make a combiner for sigma clipping and median combine
    for i in range(0,len(dark)):
        dark[i] = ccdproc.subtract_bias(dark[i],master_bias)
    a_combiner = overscan_trim_and_sigma_clip_median(dark)
    master_dark = a_combiner.median_combine(median_func=np.ma.median)
    master_dark.data[master_dark.data < 0] = 0
    master_dark.data = master_dark.data/dark[0].header["EXPOSURE"]
    master_dark.header["EXPOSURE"] = 1
    d_min, d_max, d_mean, d_std = imstats(np.asarray(master_dark))
    pl.figure(figsize=(8, 8))
    pl.imshow(master_dark, vmax=d_mean + 4 * d_std, vmin=0,cmap='gray')
    pl.colorbar()
    pl.savefig(args.output + "master_dark.pdf")
    pl.close()
    return master_dark

def get_master_flat(flat : Dict[str,List[CCDData]],master_bias : CCDData,master_dark : CCDData)->Dict[str,CCDData]:
    master_flat = {}
    median_cutoff = 4000
    for key, value in flat.items():
        flat_list = []
        try:
            makedirs(args.output+f"{key}")
        except:
            pass
        for data,i in zip(value,range(0,len(value))):
            print(f"Mean value {np.median(data.data)}")
            if np.median(data.data) < median_cutoff:
                continue
            data = ccdproc.subtract_bias(data, master_bias)
            data = ccdproc.subtract_dark(data,master_dark,exposure_time="EXPOSURE",scale=True, exposure_unit=u.second)
            data.data = data.data / np.median(data.data)
            flat_list.append(data)

            pl.figure(figsize=(16, 16))
            pl.title(f"Master flat {key}")
            d_min, d_max, d_mean, d_std = imstats(np.asarray(data))
            pl.imshow(data, vmax=d_mean + 4 * d_std, vmin=d_mean - 4 * d_std, cmap='gray')
            pl.colorbar()
            pl.savefig(args.output + f"/{key}/{key}_flat_{i}.pdf")
            pl.close()

        combined_flat = ccdproc.combine(flat_list, method='median')
        median = np.median(combined_flat.data)
        std = np.std(combined_flat.data)
        #flat_mask = np.logical_or(combined_flat.data < median - 3 * std, combined_flat.data > median + 3 * std)
        print(median)
        print(std)
        #combined_flat.data[flat_mask] = median
        d_min, d_max, d_mean, d_std = imstats(np.asarray(combined_flat))
        combined_flat.data = combined_flat.data / np.amax(combined_flat.data)
        pl.figure(figsize=(16, 16))
        pl.title(f"Master flat {key}")
        pl.imshow(combined_flat, vmax=d_mean + 4 * d_std, vmin=d_mean - 4 * d_std, cmap='gray')
        pl.colorbar()
        pl.savefig(args.output + f"master_flat{key}.pdf")
        pl.close()
        master_flat[key] = combined_flat
        CCDData.write(combined_flat, args.output + f"master_flat_{key}.fits")

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
    image = ccdproc.subtract_dark(image, master_dark,
                                  dark_exposure=master_dark.header["EXPOSURE"] * u.second,
                                  data_exposure=image.header["EXPOSURE"] * 100*u.second,
                                  scale=True,
                                  exposure_unit=u.second)
    image = ccdproc.flat_correct(image, master_flat[image.header["FILTER"]])
    # image.data = (ccdproc.subtract_bias(image,master_bias).data - tmp_dark)/master_flat[image.header["FILTER"]]
    image.data[image.data < 0] = 0

    ax = fig.add_subplot(122)
    ax.set_title(f"Image after")
    im = ax.imshow(image, cmap='gray', vmin=0)  # ,vmax=np.median(image.data)*2)
    pl.colorbar(im, ax=ax)
    pl.savefig(args.output + f"{key}.pdf")
    #pl.show()
    pl.close()
    CCDData.write(image, args.output + f"{key}")
