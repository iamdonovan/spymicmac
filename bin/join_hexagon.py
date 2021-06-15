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
import sPyMicMac.image as imtools


def blend_halves(_left, _right, left_shape):
    first = np.where(np.sum(_right, axis=0) > 0)[0][0]
    last = left_shape[1]

    m = 1 / (first - last)
    alpha = np.ones(_left.shape, dtype=np.float32)
    alpha[:, last:] = 0
    for i, ind in enumerate(np.arange(first, last)):
        alpha[:, ind] = 1 + m * (ind - first)

    return alpha * _left + (1 - alpha) * _right


def find_transform(left, right, args):
    src_pts = []
    dst_pts = []

    row_inds = list(range(0, left.shape[0] + 1, args.overlap))
    if row_inds[-1] != left.shape[0]:
        row_inds.append(-1)
    row_inds = np.array(row_inds)

    for i, ind in enumerate(row_inds[:-1]):
        try:
            l_sub = left[ind:row_inds[i + 1], -args.overlap:]
            r_sub = right[ind:row_inds[i + 1], :args.overlap]

            kp, des, matches = imtools.get_matches(l_sub, r_sub)
            src_pts.extend(
                [np.array(kp[0][m.queryIdx].pt) + np.array([left.shape[1] - args.overlap, ind]) for m in matches])
            dst_pts.extend([np.array(kp[1][m.trainIdx].pt) + np.array([0, ind]) for m in matches])
        except:
            continue

    M, inliers = ransac((np.array(src_pts), np.array(dst_pts)), EuclideanTransform,
                        min_samples=25, residual_threshold=1, max_trials=1000)
    print('{} tie points found'.format(np.count_nonzero(inliers)))
    return M


def _argparser():
    parser = argparse.ArgumentParser(description="", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-o', '--overlap', action='store', type=int, default=2000,
                        help='overlap search width between two images [2000]')
    parser.add_argument('-p', '--pattern', action='store', type=str, default='DZB',
                       help='Match pattern for images [DZB]')
    parser.add_argument('-b', '--blend', action='store_true',
                        help='Blend across image halves to prevent a sharp line at edge.')
    parser.add_argument('-c', '--corona', action='store_true',
                        help='Images are Corona KH-4/KH-4A scans (i.e., there are 4 parts).')
    return parser


parser = _argparser()
args = parser.parse_args()

imlist = [f.split('_a.tif')[0] for f in glob('{}*a.tif'.format(args.pattern))]
for im in imlist:
    print(im)

    if not args.corona:
        left = imread('{}_a.tif'.format(im))
        right = imread('{}_b.tif'.format(im))

        left_gd = gdal.Open('{}_a.tif'.format(im))
        right_gd = gdal.Open('{}_b.tif'.format(im))

        M = find_transform(left, right, args)

        out_shape = (left.shape[0], left.shape[1] + right.shape[1])

        combined_right = warp(right, M, output_shape=out_shape, preserve_range=True, order=3)

        combined_left = np.zeros(out_shape)
        combined_left[:, :left.shape[1]] = left

        if args.blend:
            combined = blend_halves(combined_left, combined_right, left.shape)
        else:
            combined = combined_left + combined_right

        last_ind = np.where(np.sum(combined, axis=0) > 0)[0][-1]
        combined = combined[:, :last_ind+1]

        imsave('{}.tif'.format(im), combined.astype(np.uint8))
    else:
        left = imread('{}_d.tif'.format(im))
        for right_img in ['c', 'b', 'a']:
            right = imread('{}_{}.tif'.format(im, right_img))

            M = find_transform(left, right, args)
            out_shape = (left.shape[0], left.shape[1] + right.shape[1])

            combined_right = warp(right, M, output_shape=out_shape, preserve_range=True, order=3)

            combined_left = np.zeros(out_shape)
            combined_left[:, :left.shape[1]] = left

            if args.blend:
                combined = blend_halves(combined_left, combined_right, left.shape)
            else:
                combined = combined_left + combined_right

            last_ind = np.where(np.sum(combined, axis=0) > 0)[0][-1]
            left = combined[:, :last_ind + 1]

        imsave('{}.tif'.format(im), combined.astype(np.uint8))

    # out_gd = gdal.Open('{}.tif'.format(im), gdal.GA_Update)
    # out_meta = out_gd.GetMetadata_Dict()

    # out_meta['TIFFTAG_XRESOLUTION'] = str(xres)
    # out_meta['TIFFTAG_YRESOLUTION'] = str(yres)

    # out_gd.SetMetadata(out_meta)

    # img = pyvips.Image.new_from_file('tmp.tif', memory=True)
    # img_bal = img.hist_equal()
    # img_bal.write_to_file('{}.tif'.format(im))
    # os.remove('tmp.tif')

# if __name__ == "__main__":
#     main()
