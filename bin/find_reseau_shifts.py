#!/usr/bin/env python
from __future__ import print_function, division
import os
import errno
import argparse
import cv2
import multiprocessing as mp
from functools import partial
from scipy.ndimage.filters import median_filter
# from skimage.io import imsave
from skimage.morphology import disk
from skimage.filters import rank
from scipy.interpolate import RectBivariateSpline as RBS
# import skimage.transform as tf
from scipy import ndimage
import numpy as np
import matplotlib.pyplot as plt
import gdal
# import pyvips
from numba import jit
from llc import jit_filter_function
import pandas as pd
import lxml.etree as etree
import lxml.builder as builder
from pybob.bob_tools import mkdir_p
import sPyMicMac.image_tools as imtools
import sPyMicMac.micmac_tools as mmt


def get_im_meas(gcps, E):
    pt_els = []
    for ind, row in gcps.iterrows():
        this_mes = E.OneMesureAF1I(
                        E.NamePt(row['gcp']),
                        E.PtIm('{} {}'.format(row['im_col'], row['im_row']))
                        )
        pt_els.append(this_mes)
    return pt_els


def make_grid():
    pass


def get_grid_matches():
    pass


def warp_image():
    pass


def _argparser():
    parser = argparse.ArgumentParser(description="Find Reseau marks and warp Hexagon image.",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('img', action='store', type=str, help='Image to warp.')
    # parser.add_argument('right_img', action='store', type=str, help='Right image to warp.')
    parser.add_argument('-orig', action='store', nargs=2, type=float, default=None,
                        help='Location of origin to search from in image [row, column]')
    # parser.add_argument('-r_orig', action='store', nargs=2, type=int, default=None,
    #                    help='Location of origin to search from in right image [row, column]')
    parser.add_argument('-csize', action='store', type=int, default=361, help='cross size [361 pixels]')
    parser.add_argument('-tsize', action='store', type=int, default=300, help='half-size of search window [300 pixels]')
    parser.add_argument('-scanres', action='store', type=float, default=None,
                        help='Scanning resolution in pixels per cm [None]')
    parser.add_argument('-scanres_y', action='store', type=float, default=None,
                        help='Scanning resolution in pixels per cm in y direction [None]')
    parser.add_argument('-filt_size', action='store', type=int, default=5, help='Size of median filter to use [5x5]')
    parser.add_argument('-nproc', action='store', type=int, default=1,
                        help='Number of processors to use [1].')
    parser.add_argument('-joined', action='store_true',
                        help='Use if img represents an already re-joined image.')
    parser.add_argument('-s', '--save_outputs', action='store_true',
                        help='Save image showing location of match for each reseau mark.')
    return parser


def main():
    np.seterr(divide='ignore', invalid='ignore')
    parser = _argparser()
    args = parser.parse_args()

    # outname = os.path.splitext(args.left_img)[0].split('_')[0] + '_shift' + os.path.splitext(args.left_img)[-1]

    print('Reading {}'.format(args.img))

    gd = gdal.Open(args.img)
    metadata = gd.GetMetadata_Dict()
    img = gd.ReadAsArray()
    print('Image read.')

    if args.orig is None:
        raise ValueError('Origin must be provided for image.')

    tmp_cross = imtools.cross_template(args.csize, width=3)
    # cross = np.random.randint(0, 255, size=(361, 361)).astype(np.uint8)
    cross = np.ones((args.csize, args.csize)).astype(np.uint8)
    cross[tmp_cross == 1] = 255

    orig_y, orig_x = args.orig
    subimg, _, _ = imtools.make_template(img, args.orig, args.tsize)
    # find the origin in the given image
    res, orig_i, orig_j = imtools.find_match(subimg, cross)
    oi, oj = orig_i-args.tsize+orig_y, orig_j-args.tsize+orig_x

    this_res_list = [res]

    print('Found origin location: {}, {}'.format(oi, oj))
    if args.save_outputs:
        mkdir_p('gcp_imgs')
        plt.figure(figsize=(8,8))
        plt.imshow(subimg, cmap='gray')
        plt.plot(orig_j, orig_i, 'r+', ms=8, linewidth=2)
        plt.savefig('gcp_imgs/{}_match_0_0.png'.format(os.path.splitext(args.img)[0]), bbox_inches='tight', dpi=200)
        plt.close()

    if args.scanres is None:
        npix_x = np.float(metadata['TIFFTAG_XRESOLUTION'])
        npix_y = np.float(metadata['TIFFTAG_YRESOLUTION'])
    else:
        npix_x = args.scanres
        if args.scanres_y is None:
            npix_y = args.scanres
        else:
            npix_y = args.scanres_y

    # now, go every 10mm in x, y direction until we cover the image.
    # should have 47 marks across, 23 up
    i_list = -npix_y * np.arange(22, -1, -1) + oi
    # make j variable based on origin location, size of image
    # if img.shape[1] > 40000:
    #    j_list = npix * np.arange(0, 47) + oj
    # else:
    if not args.joined:
        full_j = 23 * npix_x + oj < img.shape[1]

        if full_j:
            j_list = npix_x * np.arange(0, 24) + oj
        else:
            j_list = npix_x * np.arange(0, 23) + oj
    else:
        j_list = npix_x * np.arange(0,47) + oj

    nj = j_list.size

    J, I = np.meshgrid(np.arange(0, j_list.size), np.arange(0, i_list.size))
    gcp_names = list(zip(I[0, :], J[0, :]))
    for i in range(1, i_list.size):
        gcp_names.extend(list(zip(I[i, :], J[i, :])))

    JJ, II = np.meshgrid(np.round(j_list).astype(int), np.round(i_list).astype(int))
    ij = list(zip(II[-1, :], JJ[-1, :]))
    for i in range(21, -1, -1):
        ij.extend(list(zip(II[i, :], JJ[i, :])))

    print('Finding grid points in image...')
    # now that we have a list of coordinates, we can loop through and find the real grid locations
    warped_ij = [(oi, oj)] # using oi, oj to preserve float values
    if args.nproc > 1:
        subimgs = []
        for loc in ij[1:]:
            subimg, _, _ = imtools.make_template(img, loc, args.tsize)
            subimgs.append(subimg)
        pool = mp.Pool(args.nproc)
        outputs = pool.map(partial(imtools.find_match, template=cross), subimgs)
        pool.close()
        # have to map outputs to warped_ij
        for n, output in enumerate(outputs):
            res, this_i, this_j = output
            this_res_list.append(res)
            rel_i = this_i - args.tsize
            rel_j = this_j - args.tsize
            i, j = ij[n+1]
            warped_ij.append((rel_i+i, rel_j+j))
            if args.save_outputs:
                plt.figure(figsize=(8,8))
                plt.imshow(subimgs[n], cmap='gray')
                plt.plot(this_j, this_i, 'r+', ms=8, linewidth=2)
                plt.savefig('gcp_imgs/{}_match_{}_{}.png'.format(os.path.splitext(args.img)[0], gcp_names[n+1][0],
                                                                 gcp_names[n+1][1]), bbox_inches='tight', dpi=200)
                plt.close()
    else:
        for loc in ij[1:]:
            i, j = loc
            subimg, _, _ = imtools.make_template(img, loc, args.tsize)
            # img_eq = rank.equalize(subimg, selem=selem)
            res, this_i, this_j = imtools.find_match(subimg, cross)
            warped_ij.append((this_i+i, this_j+j))

    # res_list.append(this_res_list)
    ij = np.array(ij)
    warped_ij = np.array(warped_ij)
    warp = warped_ij - ij

    ux = np.reshape(warp[:,1], (23, nj))
    uy = np.reshape(warp[:,0], (23, nj))

    if args.filt_size > 0:
        ux_ = median_filter(ux, size=args.filt_size, mode='nearest')
        uy_ = median_filter(uy, size=args.filt_size, mode='nearest')
        # see how much the vectors changed after the median filter
        dx_ = ux - ux_
        dy_ = uy - uy_
        # unless there are huge deflections, keep the original values
        ux_[np.abs(dx_) < 5] = ux[np.abs(dx_) < 5]
        uy_[np.abs(dy_) < 5] = uy[np.abs(dy_) < 5]
    else:
        ux_ = ux
        uy_ = uy

    map_x = ij[:,1]
    map_y = ij[:,0]
    image_x = map_x + ux_.reshape(-1)
    image_y = map_y + uy_.reshape(-1)
    # image_x = warped_ij[:,1]
    # image_y = warped_ij[:,0]

    # make sure that the first point is the same - no shift for the origin point.
    image_x[0] = map_x[0]
    image_y[0] = map_y[0]

    # warp_list.append(warp)
    # ij_list.append(ij)
    # img_xy_list.append(np.column_stack((image_x,image_y)))

    print('Grid points found.')
    # print('Warping {} left image.'.format(name_list[ind]))
    # generate a transform using the reseau field
    # t = tf.PolynomialTransform()
    # t.estimate(warped_ij, ij)

    # apply the transformation with bicubic resampling
    # warped_img = tf.warp(img_list[0], t, order=3)

    # warped_list.append(warped_img)
    # save a results image so we can see if there's anything suspicious
    plt.figure(figsize=(12,12))
    plt.imshow(img[::10, ::10], extent=[0, img.shape[1], img.shape[0], 0], cmap='gray')
    plt.quiver(image_x, image_y, -ux_.reshape(-1), -uy_.reshape(-1), color='r')
    # plt.quiver(map_x, map_y, warp[:,1], warp[:,0], color='r')
    #
    this_out = os.path.splitext(args.img)[0]
    plt.savefig(this_out + '_matches.png', bbox_inches='tight', dpi=200)
    plt.close()
    #
    # gcp_list = []
    gcp_df = pd.DataFrame()
    for i, ind in enumerate(gcp_names):
        row, col = ind
        gcp_df.loc[i, 'gcp'] = 'GCP_{}_{}'.format(row, col)
        gcp_df.loc[i, 'im_row'] = image_y[i]
        gcp_df.loc[i, 'im_col'] = image_x[i]

    E = builder.ElementMaker()
    ImMes = E.MesureAppuiFlottant1Im(E.NameIm(args.img))

    pt_els = mmt.get_im_meas(gcp_df, E)
    for p in pt_els:
        ImMes.append(p)
    mkdir_p('Ori-InterneScan')

    outxml = E.SetOfMesureAppuisFlottants(ImMes)
    tree = etree.ElementTree(outxml)
    tree.write('Ori-InterneScan/MeasuresIm-' + args.img + '.xml', pretty_print=True,
               xml_declaration=True, encoding="utf-8")


if __name__ == "__main__":
    main()
