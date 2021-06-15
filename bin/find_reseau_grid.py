#!/usr/bin/env python
import os
import argparse
import cv2
import multiprocessing as mp
from functools import partial
from PIL import Image
from shapely.geometry import LineString
from skimage.feature import peak_local_max
from skimage.morphology import disk
from skimage import filters
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
import gdal
import pandas as pd
import lxml.etree as etree
import lxml.builder as builder
from pybob.bob_tools import mkdir_p
from pybob.ddem_tools import nmad
import sPyMicMac.image as imtools
import sPyMicMac.micmac as mmt


def moving_average(a, n=5):
    ret = np.cumsum(a)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


def lsq_fit(x, y, z):
    tmp_A = []
    tmp_b = []

    for kk in range(z.size):
        _ind = np.unravel_index(kk, z.shape)
        tmp_A.append([x[_ind], y[_ind], 1])
        tmp_b.append(z[_ind])

    A = np.matrix(tmp_A)
    b = np.matrix(tmp_b).T

    A = A[np.isfinite(z.reshape(-1, 1)).any(axis=1)]
    b = b[np.isfinite(b)]

    fit = (A.T * A).I * A.T * b.T
    out_lsq = np.zeros(z.shape)
    for kk in range(z.size):
        _ind = np.unravel_index(kk, z.shape)
        out_lsq[_ind] = fit[0] * x[_ind] + fit[1] * y[_ind] + fit[2]
    return out_lsq


def find_cross(_img, _pt, _cross, _args):
    _subimg, _, _ = imtools.make_template(_img, _pt, _args.tsize)
    res, this_i, this_j = imtools.find_match(_subimg, _cross)
    inv_res = res.max() - res

    pks = peak_local_max(inv_res, min_distance=5, num_peaks=2)
    this_z_corrs = []
    for pk in pks:
        max_ = inv_res[pk[0], pk[1]]
        this_z_corrs.append((max_ - inv_res.mean()) / inv_res.std())
    if max(this_z_corrs) > 4 and max(this_z_corrs)/min(this_z_corrs) > 1.15:
        return this_i-_args.tsize+_pt[0], this_j-_args.tsize+_pt[1]
    else:
        return np.nan, np.nan


def get_im_meas(gcps, E):
    pt_els = []
    for ind, row in gcps.iterrows():
        this_mes = E.OneMesureAF1I(
                        E.NamePt(row['name']),
                        E.PtIm('{} {}'.format(row['match_j'], row['match_i']))
                        )
        pt_els.append(this_mes)
    return pt_els


def make_grid(ledge, redge):
    i_grid = []
    j_grid = []
    search_pts = []

    for ii in range(0, 23):
        this_left = ledge.interpolate(ii * left_edge.length / 22)
        this_right = redge.interpolate(ii * right_edge.length / 22)
        this_line = LineString([this_left, this_right])
        for jj in range(0, 24):
            this_pt = this_line.interpolate(jj * this_line.length / 23)
            i_grid.append(ii)
            j_grid.append(jj)
            search_pts.append((this_pt.y, this_pt.x))

    return i_grid, j_grid, search_pts


def get_perp(corners):
    a = np.array([corners[1][0] - corners[0][0],
                  corners[1][1] - corners[0][1]])
    a = a / np.linalg.norm(a)
    return np.array([-a[1], a[0]])


def _argparser():
    _parser = argparse.ArgumentParser(description="Find Reseau marks and warp Hexagon image.",
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    _parser.add_argument('img', action='store', type=str, help='Image to warp.')
    _parser.add_argument('-orig', action='store', nargs=2, type=float, default=None,
                         help='Location of origin to search from in image [row, column]')
    _parser.add_argument('-csize', action='store', type=int, default=361, help='cross size [361 pixels]')
    _parser.add_argument('-tsize', action='store', type=int, default=300, help='half-size of search window [300 pixels]')
    _parser.add_argument('-scanres', action='store', type=float, default=None,
                         help='Scanning resolution in pixels per cm [None]')
    _parser.add_argument('-scanres_y', action='store', type=float, default=None,
                         help='Scanning resolution in pixels per cm in y direction [None]')
    _parser.add_argument('-nproc', action='store', type=int, default=1,
                         help='Number of processors to use [1].')
    return _parser


def main():
    np.seterr(divide='ignore', invalid='ignore')
    parser = _argparser()
    args = parser.parse_args()

    print('Reading {}'.format(args.img))
    gd = gdal.Open(args.img)
    metadata = gd.GetMetadata_Dict()
    img = gd.ReadAsArray()
    print('Image read.')

    img_ = Image.fromarray(img)
    img_lowres = np.array(img_.resize((np.array(img_.size)/10).astype(int), Image.LANCZOS))
    del img_

    tmp_cross = imtools.cross_template(args.csize, width=3)
    cross = np.ones((args.csize, args.csize)).astype(np.uint8)
    cross[tmp_cross == 1] = 255

    fig = plt.figure(figsize=(7, 12))
    ax = fig.add_subplot(111)
    ax.imshow(img_lowres, cmap='gray', extent=[0, img.shape[1], img.shape[0], 0])

    if args.orig is None:
        img_lowres_mask = np.zeros(img_lowres.shape)
        img_lowres_mask[np.logical_or(img_lowres < 25,
                                      filters.median(img_lowres, selem=disk(16)) < 25)] = 1

        vert = np.count_nonzero(filters.sobel_v(img_lowres_mask)**2 > 0.8, axis=0)
        hori = np.count_nonzero(filters.sobel_h(img_lowres_mask)**2 > 0.8, axis=1)

        vpeaks = peak_local_max(moving_average(vert), min_distance=6000, num_peaks=1, exclude_border=10)
        hpeaks = peak_local_max(moving_average(hori), min_distance=3000, num_peaks=2, exclude_border=10)

        top, bot = 10 * hpeaks.min(), 10 * hpeaks.max()
        edge = 10 * vpeaks.min()
        if edge < 5000:
            search_corners = [(bot-625, edge+180), (top+625, edge+180)]
        else:
            search_corners = [(bot - 625, edge-180), (top + 625, edge-180)]

        # corners = [(bot, lft), (bot, rgt), (top, rgt), (top, lft)]
        # search_corners = [(bot-625, lft+180), (bot-625, rgt-180),
        #                   (top+625, rgt-180), (top+625, lft+180)] # average distance of cross centers to img corners
        grid_corners = []
        for c in search_corners:
            _i, _j = find_cross(img, c, cross, args)
            if any(np.isnan([_j, _i])):
                grid_corners.append((c[1], c[0]))
            else:
                grid_corners.append((_j, _i))
        pixres = (grid_corners[0][1] - grid_corners[1][1]) / 22
        npix = 23 * pixres

        perp = get_perp(grid_corners)

        if edge < 5000:
            left_edge = LineString([grid_corners[0], grid_corners[1]])
            right_edge = LineString([grid_corners[0] + npix * perp,
                                     grid_corners[1] + npix * perp])
        else:
            right_edge = LineString([grid_corners[0], grid_corners[1]])
            left_edge = LineString([grid_corners[0] - npix * perp,
                                    grid_corners[1] - npix * perp])

    else:
        orig_y, orig_x = args.orig
        oj, oi = find_cross(img, args.orig, cross, args)
        if any(np.isnan([oj, oi])):
            print('Failed to find valid origin. Falling back to the one provided.')
            oi, oj = orig_y, orig_x
        else:
            print('Found origin location: {}, {}'.format(oi, oj))

        if args.scanres is None:
            npix_x = np.float(metadata['TIFFTAG_XRESOLUTION'])
            npix_y = np.float(metadata['TIFFTAG_YRESOLUTION'])
        else:
            npix_x = args.scanres
            if args.scanres_y is None:
                npix_y = args.scanres
            else:
                npix_y = args.scanres_y

        left_edge = LineString([(oj, oi), (oj, oi+22*npix_y)])
        right_edge = LineString([(oj+23*npix_x, oi), (oj+23*npix_x, oi+23*npix_y)])

    II, JJ, search_grid = make_grid(left_edge, right_edge)
    matched_grid = []

    print('Finding grid points in image...')
    # now that we have a list of coordinates, we can loop through and find the real grid locations
    if args.nproc > 1:
        subimgs = []
        std_devs = []
        means = []
        for loc in search_grid:
            subimg, _, _ = imtools.make_template(img, loc, args.tsize)
            subimgs.append(subimg)
        pool = mp.Pool(args.nproc)
        outputs = pool.map(partial(imtools.find_match, template=cross), subimgs)
        pool.close()
        # have to map outputs to warped_ij
        for n, output in enumerate(outputs):
            res, this_i, this_j = output
            maxj, maxi = cv2.minMaxLoc(res)[2]  # 'peak' is actually the minimum location, remember.
            inv_res = res.max() - res
            std_devs.append(inv_res.std())
            means.append(inv_res.mean())

            pks = peak_local_max(inv_res, min_distance=5, num_peaks=2)
            if pks.size > 0:
                this_z_corrs = []
                for pk in pks:
                    max_ = inv_res[pk[0], pk[1]]
                    this_z_corrs.append((max_ - inv_res.mean()) / inv_res.std())

                if max(this_z_corrs) > 5 and max(this_z_corrs)/min(this_z_corrs) > 1.15:
                    rel_i = this_i - args.tsize
                    rel_j = this_j - args.tsize
                    i, j = search_grid[n]
                    matched_grid.append((rel_i+i, rel_j+j))
                else:
                    matched_grid.append((np.nan, np.nan))
            else:
                matched_grid.append((np.nan, np.nan))
    else:
        for loc in search_grid:
            _j, _i = find_cross(img, loc, cross, args)
            matched_grid.append((_j, _i))

    matched_grid = np.array(matched_grid)
    search_grid = np.array(search_grid)

    gcps_df = pd.DataFrame()
    for i, pr in enumerate(list(zip(II, JJ))):
        gcps_df.loc[i, 'gcp'] = 'GCP_{}_{}'.format(pr[0], pr[1])

    gcps_df['search_j'] = search_grid[:, 1]
    gcps_df['search_i'] = search_grid[:, 0]
    gcps_df['match_j'] = matched_grid[:, 1]
    gcps_df['match_i'] = matched_grid[:, 0]

    dx = gcps_df['match_j'] - gcps_df['search_j']
    dy = gcps_df['match_i'] - gcps_df['search_i']

    nomatch = np.isnan(gcps_df.match_j)

    # ux = griddata(gcps_df.loc[~nomatch, ['search_j', 'search_i']], gcps_df.loc[~nomatch, 'dj'],
    #               gcps_df[['search_j', 'search_i']])
    # uy = griddata(gcps_df.loc[~nomatch, ['search_j', 'search_i']], gcps_df.loc[~nomatch, 'di'],
    #               gcps_df[['search_j', 'search_i']])

    ux = lsq_fit(gcps_df.search_j.values, gcps_df.search_i.values, gcps_df.match_j.values)
    uy = lsq_fit(gcps_df.search_j.values, gcps_df.search_i.values, gcps_df.match_i.values)

    gcps_df.loc[nomatch, 'match_j'] = ux[nomatch]
    gcps_df.loc[nomatch, 'match_i'] = uy[nomatch]

    gcps_df['dj'] = gcps_df['match_j'] - gcps_df['search_j']
    gcps_df['di'] = gcps_df['match_i'] - gcps_df['search_i']

    gcps_df['im_row'] = gcps_df['match_i']
    gcps_df['im_col'] = gcps_df['match_j']

    print('Grid points found.')
    mkdir_p('match_imgs')

    ax.quiver(gcps_df.search_j, gcps_df.search_i, gcps_df.dj, gcps_df.di, color='r')
    ax.plot(gcps_df.search_j[nomatch], gcps_df.search_i[nomatch], 'b+')
    this_out = os.path.splitext(args.img)[0]
    fig.savefig('match_imgs/' + this_out + '_matches.png', bbox_inches='tight', dpi=200)
    # fig.close()

    E = builder.ElementMaker()
    ImMes = E.MesureAppuiFlottant1Im(E.NameIm(args.img))

    pt_els = mmt.get_im_meas(gcps_df, E)
    for p in pt_els:
        ImMes.append(p)
    mkdir_p('Ori-InterneScan')

    outxml = E.SetOfMesureAppuisFlottants(ImMes)
    tree = etree.ElementTree(outxml)
    tree.write('Ori-InterneScan/MeasuresIm-' + args.img + '.xml', pretty_print=True,
               xml_declaration=True, encoding="utf-8")


if __name__ == "__main__":
    main()
