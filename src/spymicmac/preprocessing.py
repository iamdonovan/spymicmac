"""
spymicmac.preprocessing is a collection of tools for preprocessing images
"""
import os
import argparse
import shutil
import tarfile
import numpy as np
from glob import glob
from skimage import io, filters
from spymicmac.image import join_hexagon
from spymicmac.matching import find_reseau_grid, remove_crosses
from spymicmac.resample import resample_hex
from spymicmac import micmac


def extract(tar_ext='.tgz'):
    """
    Extract image parts from tar files in the current directory.

    :param str tar_ext: the extension of the files to look for (default: .tgz)
    """
    # make a directory to move the tarballs to
    os.makedirs('tarballs', exist_ok=True)

    tarlist = glob('*' + tar_ext)
    tarlist.sort()

    if len(tarlist) > 0:
        print('Extracting images from tar files.')
        for tarball in tarlist:
            print(tarball)
            with tarfile.open(tarball, 'r') as tfile:
                tfile.extractall('.')
            shutil.move(tarball, 'tarballs')
    else:
        print('No tar files found, skipping.')
        # we still want to make sure that we have a list of images to work with.
        tarlist = list(set([os.path.splitext(fn)[0].split('_')[0] for fn in glob('DZB*.tif')]))
        tarlist.sort()

    return tarlist


def check_reseau():
    """
    Check whether all KH-9 Mapping Camera images in the current directory have a corresponding Measures{fn_img}.xml
    file in Ori-InterneScan.

    :return:
    """
    imlist = glob('DZB*.tif')
    measlist = glob('MeasuresIm*.tif.xml', dir_fd='Ori-InterneScan')

    # if there are no DZB*.tif, but there are MeasuresIm files, we've probably done this step
    if len(imlist) == 0 and len(measlist) > 0:
        return True
    else:
        # check that each image file has a matching measuresim file
        return all([any([im in meas for meas in measlist]) for im in imlist])

