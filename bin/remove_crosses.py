#!/usr/bin/env python
import os
import shutil
import argparse
from glob import glob
import pandas as pd 
import xml.etree.ElementTree as ET 
from skimage.io import imsave, imread
from skimage.morphology import disk
import numpy as np 
import matplotlib.pyplot as plt
from pybob.image_tools import nanmedian_filter
from pybob.bob_tools import mkdir_p
import sPyMicMac.image as imtools


def fix_cross(subimg):
    subimg = subimg.astype(np.float32)

    cross = imtools.cross_template(subimg.shape[0], width=5)
    cross[:, :16] = 0
    cross[:, -16:] = 0
    cross[:16, :] = 0
    cross[-16:, :] = 0
    if subimg.shape[0] != cross.shape[0]:
        cross = cross[:subimg.shape[0], :]

    if subimg.shape[1] != cross.shape[1]:
        cross = cross[:, :subimg.shape[1]]

    # print(cross.shape)
    subimg[cross != 0] = np.nan
    fixed = nanmedian_filter(subimg, footprint=disk(7))
    subimg[np.isnan(subimg)] = fixed[np.isnan(subimg)]
    return subimg.astype(np.uint8)


def parse_im_meas(fn_meas): 
    gcp_df = pd.DataFrame() 
    root = ET.parse(fn_meas).getroot() 
    measures = root.findall('MesureAppuiFlottant1Im')[0] 
    for i, mes in enumerate(measures.findall('OneMesureAF1I')): 
        gcp_df.loc[i, 'name'] = mes.find('NamePt').text 
        pt = mes.find('PtIm').text.split() 
        gcp_df.loc[i, 'i'] = float(pt[1]) 
        gcp_df.loc[i, 'j'] = float(pt[0]) 
    return gcp_df


def _argparser():
    parser = argparse.ArgumentParser(description="Remove Reseau marks from KH-9 image(s).",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('img', action='store', type=str, nargs='+', help='Image(s) to remove crosses from.')
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    mkdir_p('original')

    for im in args.img:
        print(im)

        fn_meas = os.path.join('Ori-InterneScan', 'MeasuresIm-{}.xml'.format(im))
        img = imread(im)
        gcps = parse_im_meas(fn_meas)

        pt = np.round([gcps.i[0], gcps.j[0]]).astype(int)
        subim, row_, col_ = imtools.make_template(img, pt, 200)
        cross = imtools.cross_template(subim.shape[0], width=5)
        cross[:, 16:22] = 1
        cross[:, -24:-16] = 1
        cross[16:22, :] = 1
        cross[-24:-16, :] = 1

        subim = subim.astype(np.float32)
        subim[cross != 0] = np.nan
        fixed = nanmedian_filter(subim, footprint=disk(7))
        subim[np.isnan(subim)] = fixed[np.isnan(subim)]
        img[int(pt[0])-row_[0]:int(pt[0])+row_[1]+1, int(pt[1])-col_[0]:int(pt[1])+col_[1]+1] = subim.astype(np.uint8)

        for i, row in gcps.loc[1:].iterrows():
            pt = np.round([row.i, row.j]).astype(int)
            subim, row_, col_ = imtools.make_template(img, pt, 200)
            img[pt[0]-row_[0]:pt[0]+row_[1]+1, pt[1]-col_[0]:pt[1]+col_[1]+1] = fix_cross(subim)

        shutil.move(im, 'original')
        imsave(im, img.astype(np.uint8))


if __name__ == "__main__":
    main()
