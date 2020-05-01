#!/usr/bin/env python
import argparse
import os
from glob import glob
from skimage.io import imread, imsave
from skimage.measure import ransac
from skimage.transform import ProjectiveTransform, AffineTransform, EuclideanTransform, warp
import gdal
# import pyvips
import numpy as np
import sPyMicMac.image_tools as imtools


def _argparser():
    parser = argparse.ArgumentParser(description="", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-o', '--overlap', action='store', type=int, default=2000,
                        help='overlap search width between two images [2000]')
    return parser


parser = _argparser()
args = parser.parse_args()

imlist = [f.split('_a.tif')[0] for f in glob('DZB*a.tif')]
for im in imlist:
    print(im)

    left = imread('{}_a.tif'.format(im))
    right = imread('{}_b.tif'.format(im))

    left_gd = gdal.Open('{}_a.tif'.format(im))
    right_gd = gdal.Open('{}_b.tif'.format(im))

    left_meta = left_gd.GetMetadata_Dict()
    right_meta = right_gd.GetMetadata_Dict()

    xres = np.mean(np.array([left_meta['TIFFTAG_XRESOLUTION'],
                             right_meta['TIFFTAG_XRESOLUTION']]).astype(np.float32))
    yres = np.mean(np.array([left_meta['TIFFTAG_YRESOLUTION'],
                             right_meta['TIFFTAG_YRESOLUTION']]).astype(np.float32))

    src_pts = []
    dst_pts = []

    row_inds = list(range(0, left.shape[0]+1, args.overlap))
    if row_inds[-1] != left.shape[0]:
        row_inds.append(-1)
    row_inds = np.array(row_inds)

    for i, ind in enumerate(row_inds[:-1]):
        try:
            l_sub = left[ind:row_inds[i+1], -args.overlap:]
            r_sub = right[ind:row_inds[i+1], :args.overlap]

            kp, des, matches = imtools.get_matches(l_sub, r_sub)
            src_pts.extend([np.array(kp[0][m.queryIdx].pt) + np.array([left.shape[1]-args.overlap, ind]) for m in matches])
            dst_pts.extend([np.array(kp[1][m.trainIdx].pt) + np.array([0, ind]) for m in matches])
        except:
            continue

    M, inliers = ransac((np.array(src_pts), np.array(dst_pts)), ProjectiveTransform,
                        min_samples=25, residual_threshold=1, max_trials=1000)
    print('{} tie points found'.format(np.count_nonzero(inliers)))

    out_shape = (left.shape[0], left.shape[1] + right.shape[1])

    combined_right = warp(right, M, output_shape=out_shape, preserve_range=True, order=5)

    combined_left = np.zeros(out_shape)
    combined_left[:, :left.shape[1]] = left

    first = np.where(np.sum(combined_right, axis=0) > 0)[0][0]
    last = left.shape[1]

    m = 1 / (first-last)
    alpha = np.ones(out_shape)
    alpha[:, last:] = 0
    for i, ind in enumerate(np.arange(first, last)):
        alpha[:, ind] = 1 + m * (ind-first)
    # combined_left[combined_left == 0] = np.nan
    # combined_right[combined_right == 0] = np.nan

    # combined = np.nansum(np.array([alpha * combined_left,
    #                               (1-alpha) * combined_right]), axis=0)
    # combined[np.isnan(combined)] = 0
    combined = alpha * combined_left + (1-alpha) * combined_right

    last_ind = np.where(np.sum(combined, axis=0) > 0)[0][-1]
    combined = combined[:, :last_ind+1]

    imsave('{}.tif'.format(im), combined.astype(np.uint8))

    out_gd = gdal.Open('{}.tif'.format(im), gdal.GA_Update)
    out_meta = out_gd.GetMetadata_Dict()

    out_meta['TIFFTAG_XRESOLUTION'] = str(xres)
    out_meta['TIFFTAG_YRESOLUTION'] = str(yres)

    out_gd.SetMetadata(out_meta)

    # img = pyvips.Image.new_from_file('tmp.tif', memory=True)
    # img_bal = img.hist_equal()
    # img_bal.write_to_file('{}.tif'.format(im))
    # os.remove('tmp.tif')
