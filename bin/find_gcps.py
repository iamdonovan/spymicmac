#!/usr/bin/env python
import os
import shutil
import subprocess
import argparse
import difflib
from glob import glob
from itertools import chain
import cv2
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import skimage
from skimage.io import imsave, imread
from skimage.measure import ransac
from skimage.feature import peak_local_max
from skimage.transform import match_histograms, warp, AffineTransform, EuclideanTransform
from skimage import exposure
from shapely.geometry.point import Point
import geopandas as gpd
import pandas as pd
import lxml.etree as etree
import xml.etree.ElementTree as ET
import lxml.builder as builder
from pybob.GeoImg import GeoImg
from pybob.image_tools import create_mask_from_shapefile, reshape_geoimg, match_hist
from pybob.plot_tools import plot_geoimg_sidebyside
from pybob.bob_tools import mkdir_p
from pybob.ddem_tools import nmad
import sPyMicMac.image_tools as imtools
import sPyMicMac.micmac_tools as mmt
from pymmaster.mmaster_tools import orient_footprint


def get_match_pattern(imlist):
    matches = []
    for i, this_im in enumerate(imlist[:-1]):
        for im in imlist[i + 1:]:
            matches.extend(list(difflib.SequenceMatcher(None, this_im, im).get_matching_blocks()))

    good_matches = set([m for m in matches if m.size > 0 and m.a == m.b])
    start_inds = set([m.a for m in good_matches])
    ind_lengths = [(ind, min([m.size for m in good_matches if m.a == ind])) for ind in start_inds]
    ind_lengths.sort()

    first, last = ind_lengths[0][1], ind_lengths[1][0]
    return imlist[0][:first] + '(' + '|'.join([im[first:last] for im in imlist]) + ')' + imlist[0][last:]


def write_auto_mesures(gcp_df, subscript, out_dir):
    with open(os.path.join(out_dir, 'AutoMeasures{}.txt'.format(subscript)), 'w') as f:
        for i, row in gcp_df.iterrows():
            print('{} {} {}'.format(row.rel_x, row.rel_y, row.el_rel), file=f)


def write_auto_gcps(gcp_df, subscript, out_dir, utm_zone):
    with open(os.path.join(out_dir, 'AutoGCPs{}.txt'.format(subscript)), 'w') as f:
        # print('#F= N X Y Z Ix Iy Iz', file=f)
        print('#F= N X Y Z', file=f)
        print('#Here the coordinates are in UTM {} X=Easting Y=Northing Z=Altitude'.format(utm_zone), file=f)
        for i, row in gcp_df.iterrows():
            # print('{} {} {} {} {} {} {}'.format(row.id, row.geometry.x, row.geometry.y, row.elevation,
            #                                        5/row.z_corr, 5/row.z_corr, 1), file=f)
            print('{} {} {} {}'.format(row.id, row.geometry.x, row.geometry.y, row.elevation), file=f)


def get_matches_mask(matches):
    matchesMask = [[0,0] for ii in range(len(matches))]
    for i, p in enumerate(matches):
        try:
            m, n = p
            if m.distance < 0.75*n.distance:
                matchesMask[i] = [1, 0]
        except:
            continue
    return matchesMask


def splitter(img, nblocks, overlap=0):
    split1 = np.array_split(img, nblocks[0], axis=0)
    split2 = [np.array_split(im, nblocks[1], axis=1) for im in split1]
    olist = [np.copy(a) for a in list(chain.from_iterable(split2))]
    return olist


def get_subimg_offsets(split, shape):
    ims_x = np.array([s.shape[1] for s in split])
    ims_y = np.array([s.shape[0] for s in split])

    rel_x = np.cumsum(ims_x.reshape(shape), axis=1)
    rel_y = np.cumsum(ims_y.reshape(shape), axis=0)

    rel_x = np.concatenate((np.zeros((shape[0], 1)), rel_x[:, :-1]), axis=1)
    rel_y = np.concatenate((np.zeros((1, shape[1])), rel_y[:-1, :]), axis=0)

    return rel_x.astype(int), rel_y.astype(int)


def get_bascule_residuals(fn_basc, gcp_df):
    root = ET.parse(fn_basc).getroot()
    gcp_res = root.findall('Residus')
    gcp_names = np.array([res.find('Name').text for res in gcp_res])
    # residuals = np.array([float(res.find('Dist').text) for res in gcp_res])
    x_res = np.array([float(res.find('Offset').text.split()[0]) for res in gcp_res])
    y_res = np.array([float(res.find('Offset').text.split()[1]) for res in gcp_res])
    for data_ in zip(gcp_names, x_res):
        gcp_df.loc[gcp_df.id == data_[0], 'xres'] = data_[1]

    for data_ in zip(gcp_names, y_res):
        gcp_df.loc[gcp_df.id == data_[0], 'yres'] = data_[1]

    gcp_df['residual'] = np.sqrt(gcp_df['xres'].values**2 + gcp_df['yres'].values**2)

    return gcp_df


def get_campari_residuals(fn_resids, gcp_df):
    camp_root = ET.parse(fn_resids).getroot()
    camp_gcp_names = [a.find('Name').text for a in camp_root.findall('Iters')[-1].findall('OneAppui')]
    err_max = [float(a.find('EcartImMax').text) for a in camp_root.findall('Iters')[-1].findall('OneAppui')]

    for data_ in zip(camp_gcp_names, err_max):
        gcp_df.loc[gcp_df.id == data_[0], 'camp_res'] = data_[1]

    return gcp_df


def check_mask(pts, M, mask):
    isin = np.zeros((pts.shape[0], 1), dtype=np.uint8)
    pts = np.concatenate((pts, np.ones((pts.shape[0], 1))), axis=1)
    t_pts = np.dot(M, pts.T).T.astype(int)

    for i, p in enumerate(t_pts):
        if 0 < p[0] < mask.shape[1] and 0 < p[1] < mask.shape[0]:
            if mask[p[1], p[0]] > 0:
                isin[i] = 1
    return isin


def sliding_window_filter(img_shape, pts_df, winsize, stepsize=None, mindist=2000):
    if stepsize is None:
        stepsize = winsize / 2

    out_inds = []
    out_pts = []

    for x_ind in np.arange(stepsize, img_shape[1], winsize):
        for y_ind in np.arange(stepsize, img_shape[0], winsize):
            min_x = x_ind - winsize / 2
            max_x = x_ind + winsize / 2
            min_y = y_ind - winsize / 2
            max_y = y_ind + winsize / 2
            samp_ = pts_df.loc[np.logical_and.reduce([pts_df.i > min_x,
                                                      pts_df.i < max_x,
                                                      pts_df.j > min_y,
                                                      pts_df.j < max_y])]
            if samp_.shape[0] == 0:
                continue
            samp_.sort_values('residual', ascending=True, inplace=True)
            if len(out_inds) == 0:
                best_ind = samp_.index[0]
                best_pt = Point(samp_.loc[best_ind, ['j', 'i']].values)

                out_inds.append(best_ind)
                out_pts.append(best_pt)
            else:
                for ind, row in samp_.iterrows():
                    this_pt = Point(row[['j', 'i']].values)
                    this_min_dist = np.array([this_pt.distance(pt) for pt in out_pts]).min()
                    if this_min_dist > mindist:
                        out_inds.append(ind)
                        out_pts.append(this_pt)

    return np.array(out_inds)


def get_subpixel_matches(inv_lsat, ortho, inv_mask, src_pts, dst_pts, args, M=None):
    match_pts = np.array(src_pts)
    z_corrs = []
    peak_corrs = []
    matching_results = []

    for i, src_ in enumerate(src_pts):
        try:
            # dst_pt = np.dot(invM, np.array([dst_[0], dst_[1], 1]).T)
            if M is not None:
                dst_ = dst_pts[i]
                dst_pt = Mout.inverse(dst_).flatten()
            else:
                dst_pt = dst_pts[i]
            print(src_, dst_pt)

            testchip, _, _ = imtools.make_template(inv_lsat, (dst_pt[1], dst_pt[0]), args.rsize)
            dst_chip, _, _ = imtools.make_template(ortho, (src_[1], src_[0]), args.tsize)

            test = np.ma.masked_values(imtools.highpass_filter(testchip), 0)
            dest = np.ma.masked_values(imtools.highpass_filter(dst_chip), 0)

            corr_res, this_i, this_j = imtools.find_gcp_match(dest.astype(np.float32), test.astype(np.float32))
            # find the top two peaks in the correlation image, and determine how unique they are.
            peak_corr = cv2.minMaxLoc(corr_res)[1]

            pks = peak_local_max(corr_res, min_distance=5, num_peaks=2)
            this_z_corrs = []
            for pk in pks:
                max_ = corr_res[pk[0], pk[1]]
                this_z_corrs.append((max_ - corr_res.mean()) / corr_res.std())
            dz_corr = max(this_z_corrs) / min(this_z_corrs)
            z_corr = max(this_z_corrs)

            # if the correlation peak is very high, or very unique, add it as a match
            if 0.2 < peak_corr or (z_corr > 5 and dz_corr > 1.5):
                out_i, out_j = this_i - args.tsize + src_[1], this_j - args.tsize + src_[0]
                z_corrs.append(z_corr)
                peak_corrs.append(peak_corr)
                match_pts[i] = [out_j, out_i]
                matching_results.append((i, testchip, dst_chip, corr_res))
            else:
                # src_pts[i] = (-100, -100)
                z_corrs.append(-1)
                peak_corrs.append(-1)
        except:
            # src_pts[i] = (-100, -100)
            z_corrs.append(-1)
            peak_corrs.append(-1)

    z_corrs = np.array(z_corrs)
    peak_corrs = np.array(peak_corrs)

    return match_pts, z_corrs, peak_corrs, matching_results


def _argparser():
    parser = argparse.ArgumentParser(description="", formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('ortho', action='store', type=str, help='non-referenced orthophoto mosaic')
    parser.add_argument('lsat', action='store', type=str, help='georeferenced satellite image')
    parser.add_argument('dem', action='store', type=str, help='dem')
    parser.add_argument('rel_dem', action='store', type=str, help='relative dem corresponding to ortho')
    parser.add_argument('-glacmask', action='store', type=str, default=None,
                        help='path to shapefile of glacier outlines')
    parser.add_argument('-landmask', action='store', type=str, default=None,
                        help='path to shapefile of land outlines')
    parser.add_argument('-footprint', action='store', type=str, default=None,
                        help='approximate outline of hexagon images to crop satellite image to')
    parser.add_argument('-im_subset', action='store', type=str, default=None, nargs='+',
                        help='subset of raw images to work with (default all)')
    parser.add_argument('-corr_thresh', action='store', type=float, default=0.5,
                        help='minimum correlation value to use for accepting a match.')
    parser.add_argument('-tfm_pts', action='store', type=str, default=None,
                        help='CSV containing set of 4-5 points to estimate a rough transform, in the form I,J,X,Y.')
    parser.add_argument('-gcps', action='store', type=str, default=None,
                        help='Shapefile (or CSV) containing pre-selected control points.')
    parser.add_argument('-b', '--block', action='store', type=str, default=None,
                        help='Block number to use if multiple image blocks exist in directory.')
    parser.add_argument('-rsize', action='store', type=int, default=50,
                        help='half-size of reference search window [50 pixels]')
    parser.add_argument('-tsize', action='store', type=int, default=600,
                        help='half-size of search window [600 pixels]')
    return parser


# def main():
parser = _argparser()
args = parser.parse_args()

plt.ion()
print('start.')

out_dir = 'auto_gcps'
mkdir_p(out_dir)
if args.block is not None:
    subscript = '_block{}'.format(args.block)
else:
    subscript = ''
# mkdir_p('')
# mkdir_p('')

ortho = imread(args.ortho)
ortho_ = Image.fromarray(ortho)
ortho_lowres = np.array(ortho_.resize((np.array(ortho_.size)/100).astype(int), Image.LANCZOS))

if args.im_subset is None:
    imlist = glob('OIS*.tif')
else:
    imlist = args.im_subset

if args.tfm_pts is not None:
    tfm_pts = pd.read_csv(args.tfm_pts)
    tmp_lsat = GeoImg(args.lsat)

    dst_pts = np.array([(tmp_lsat.xy2ij(pt)[1], tmp_lsat.xy2ij(pt)[0])
                        for pt in zip(tfm_pts['mapX'], tfm_pts['mapY'])]).reshape(-1, 2) * (tmp_lsat.dx/800)
    src_pts = np.array(list(zip(tfm_pts['pixelX'], np.abs(tfm_pts['pixelY'])))).reshape(-1, 2) / 100

    M, _ = cv2.estimateAffine2D(src_pts.reshape(-1, 1, 2), dst_pts.reshape(-1, 1, 2))

    lsat_lowres = reshape_geoimg(args.lsat, 800, 800)
    lowres_size = lsat_lowres.img.shape

else:
    print('finding initial transformation between {}, {}'.format(args.ortho, args.lsat))
    M, success, lowres_size = imtools.get_initial_transformation(ortho_lowres, args.lsat,
                                                            landmask=args.landmask,
                                                            footmask=args.footprint,
                                                            imlist=imlist)
    if not success:
        if args.footprint is not None:
            print('failed to find valid transformation. Trying with provided image footprints.')
            lsat_lowres = reshape_geoimg(args.lsat, 800, 800)

            fp_mask, fprint = imtools.get_footprint_mask(args.footprint, lsat_lowres, imlist, fprint_out=True)
            oprint = orient_footprint(fprint)
            h, w = ortho_lowres.shape
            src_pts = np.array([[0, 0], [0, h-1], [w-1, h-1], [w-1, 0]]).reshape(-1, 2)
            x, y = oprint.boundary.xy
            ij = [lsat_lowres.xy2ij(pt) for pt in zip(x,y)]
            i = [pt[0] for pt in ij]
            j = [pt[1] for pt in ij]
            inds = np.array([1, 0, 3, 2]).reshape(-1, 1)
            dst_pts = np.array(list(zip(np.array(j)[inds], np.array(i)[inds]))).reshape(-1, 2)
            # M, _ = cv2.estimateAffine2D(src_pts.reshape(-1, 1, 2), dst_pts.reshape(-1, 1, 2))
            M = skimage.transform.estimate_transform('euclidean', dst_pts, src_pts)
        else:
            raise RuntimeError("Failed to find valid transformation between {}, {}. "
                               "Try running again with footprint option.".format(args.ortho, args.lsat))
print('initial transformation found.')
print('loading full-resolution images.')

lsat_fullres = GeoImg(args.lsat)

epsg_str = str(lsat_fullres.epsg)
hemi_dict = {'6': 'N', '7': 'S'}
utm_str = epsg_str[-2:] + hemi_dict[epsg_str[2]]

bitdepth = np.ceil(np.log2(np.nanmax(lsat_fullres.img)))

lsat_fullres.img = (lsat_fullres.img/(2**bitdepth-1)).astype(np.float32)
# lsat_eq = match_hist(lsat_fullres.img, np.array(ortho_lowres))
# lsat_eq = lsat_fullres.img
lsat_eq = (255 * exposure.equalize_adapthist(lsat_fullres.img, kernel_size=32)).astype(np.uint8)
lsat_eq[np.isnan(lsat_fullres.img)] = 0
lsat_h, lsat_w = lsat_eq.shape

# ortho_tfm = cv2.warpAffine(ortho_lowres, M, (lowres_size[1], lowres_size[0]))
ortho_tfm = warp(ortho_lowres, M, output_shape=lowres_size, preserve_range=True)

plt.figure(figsize=(7,5))
plt.subplot(121); plt.imshow(ortho_tfm, cmap='gray')
plt.subplot(122); plt.imshow(lsat_eq[::10, ::10], cmap='gray')
plt.savefig(os.path.join(out_dir,
                         'initial_transformation{}.png'.format(subscript)),
            dpi=200, bbox_inches='tight')
plt.close(plt.gcf())

print('creating image masks')
lsat_mask = 255 * np.ones(lsat_fullres.img.shape, dtype=np.uint8)
lsat_mask[lsat_eq == 0] = 0

if args.landmask is not None:
    lmask = create_mask_from_shapefile(lsat_fullres, args.landmask, buffer=-200)
    lsat_mask[~lmask] = 0
if args.glacmask is not None:
    gmask = create_mask_from_shapefile(lsat_fullres, args.glacmask)
    lsat_mask[gmask] = 0

ortho_mask = 255 * np.ones(ortho.shape, dtype=np.uint8)
ortho_mask[ortho == 0] = 0

# break image up into subimgs
x_tiles = np.floor(ortho.shape[1] / 1000).astype(int)
y_tiles = np.floor(ortho.shape[0] / 1000).astype(int)

split_ortho = splitter(ortho, (y_tiles, x_tiles))
nblocks = len(split_ortho)

rel_x, rel_y = get_subimg_offsets(split_ortho, (y_tiles, x_tiles))

src_pts = []
dst_pts = []

print('finding keypoints in full-resolution images.')
for i, subim in enumerate(split_ortho):
    print('do one bloc {} of {} ----------------------------'.format(i, nblocks-1))
    try:
        iy, ix = np.unravel_index(i, (y_tiles, x_tiles))
        ox = rel_x[iy, ix]
        oy = rel_y[iy, ix]
        im_y, im_x = subim.shape
        corners = np.array([[ox, oy], [ox + im_x, oy],
                            [ox + im_x, oy + im_y], [ox, oy + im_y]]) / 100
        # corners = np.array([[ox, oy, 100], [ox+im_x, oy, 100],
        #                     [ox+im_x, oy+im_y, 100], [ox, oy+im_y, 100]]) / 100
        # t_x, t_y = np.dot(M, corners.T) * (800/lsat_fullres.dx) # downscaled landsat corners in initial transformation
        corners_tfm = M.inverse(corners) * (800 / lsat_fullres.dx)

        t_x, t_y = corners_tfm[:, 0], corners_tfm[:, 1]
        lind = int(max(0, min(t_x)-400))
        rind = int(min(lsat_w, max(t_x)+400))
        tind = int(max(0, min(t_y)-400))
        bind = int(min(lsat_h, max(t_y)+400))

        lsat_sub = 255 * exposure.equalize_hist(lsat_eq[tind:bind, lind:rind])
        lmask_sub = lsat_mask[tind:bind, lind:rind]

        sub_mask = ortho_mask[oy:oy+im_y, ox:ox+im_x]
        if np.any(lmask_sub > 0) and np.any(sub_mask > 0):
            kp, des, matches = imtools.get_matches(subim, lsat_sub, mask1=sub_mask, mask2=lmask_sub)
            src = np.array([kp[0][m.queryIdx].pt for m in matches])
            dst = np.array([kp[1][m.trainIdx].pt for m in matches])
            # subM, good = cv2.estimateAffine2D(src.reshape(-1, 1, 2), dst.reshape(-1, 1, 2), ransacReprojThreshold=10)
            # notmasked = check_mask(src, subM, lmask_sub)
            # good_ = np.where(np.logical_and(notmasked, good))[0]

            dist = np.array([m.distance for m in matches])
            lastind = min(10, len(dist))
            best = np.argsort(dist)[:lastind+1]  # return only the top 10 best matches
            # ind = best[np.argsort(dist[best])[0]] # return the best match
            # for each of the best matches
            for ind in best:
                src_pts.append(np.array(kp[0][matches[ind].queryIdx].pt) + np.array([ox, oy]))
                dst_pts.append(np.array(kp[1][matches[ind].trainIdx].pt) + np.array([lind, tind]))
    except:
        print('end bloc {} -------------------------------------'.format(i))
        continue
    print('end bloc {} -------------------------------------'.format(i))

print('keypoints found.')
print('finding sub-pixel keypoint locations')

# Mout, _ = cv2.estimateAffine2D(np.array(src_pts).reshape(-1, 1, 2),
#                               np.array(dst_pts).reshape(-1, 1, 2), ransacReprojThreshold=2)
Mout, inliers = ransac((np.array(src_pts), np.array(dst_pts)), AffineTransform,
                       min_samples=10, residual_threshold=25, max_trials=1000)

print('{} points used to estimate refined transformation'.format(np.count_nonzero(inliers)))
# inv_lsat = cv2.warpAffine(lsat_eq, invM, (ortho.shape[1], ortho.shape[0]))
inv_lsat = warp(lsat_eq, Mout, output_shape=ortho.shape, preserve_range=True)
inv_mask = warp(lsat_mask, Mout, output_shape=ortho.shape, preserve_range=True, order=0).astype(np.uint8)

if args.gcps is not None:
    in_gcps = gpd.read_file(args.gcps)

    dst_ = np.array([lsat_fullres.xy2ij([(row.geometry.x, row.geometry.y)]) for i, row in in_gcps.iterrows()])
    dst_pts = Mout.inverse(np.array([dst_[:,1], dst_[:,0]]).T)
    src_pts = dst_pts.copy()
    match_pts, z_corrs, peak_corrs, matching_results = get_subpixel_matches(inv_lsat, ortho, inv_mask,
                                                                            src_pts, dst_pts, M=None)
else:
    src_pts = np.array(src_pts)  # [inliers]
    dst_pts = np.array(dst_pts)  # [inliers]

    plt.figure(figsize=(8, 10))
    plt.subplot(211)
    plt.imshow(ortho[::10, ::10], cmap='gray', extent=[0, ortho.shape[1], ortho.shape[0], 0])
    plt.subplot(212)
    plt.imshow(inv_lsat[::10, ::10], cmap='gray', extent=[0, inv_lsat.shape[1], inv_lsat.shape[0], 0])
    plt.savefig(os.path.join(out_dir,
                'refined_transformation{}.png'.format(subscript)),
                bbox_inches='tight', dpi=200)
    plt.close(plt.gcf())

    fig1 = plt.figure(figsize=(8, 10))
    plt.imshow(ortho[::5, ::5], cmap='gray', extent=[0, ortho.shape[1], ortho.shape[0], 0])
    plt.plot(np.array(src_pts)[:, 0], np.array(src_pts)[:, 1], 'r^')

    match_pts, z_corrs, peak_corrs, matching_results = get_subpixel_matches(inv_lsat, ortho, inv_mask,
                                                                            src_pts, dst_pts, args, Mout)

    src_pts = np.array(src_pts)
    dst_pts = np.array(dst_pts)

    plt.plot(src_pts[:, 0], src_pts[:, 1], 'k+')
    plt.axis([0, ortho.shape[1], ortho.shape[0], 0])

    # print(src_pts.shape, dst_pts.shape)
# now, get the elevations from the dem, and write out the valid points
print('loading dems')
dem = GeoImg(args.dem)
rel_dem = GeoImg(args.rel_dem)

xy = np.array([lsat_fullres.ij2xy((pt[1], pt[0])) for pt in dst_pts]).reshape(-1, 2)

fn_tfw = args.ortho.replace('.tif', '.tfw')
with open(fn_tfw, 'r') as f:
    gt = [float(l.strip()) for l in f.readlines()]

gcps = gpd.GeoDataFrame(columns=['elevation', 'el_rel', 'rel_x', 'rel_y',
                                 'mi', 'mj', 'i', 'j', 'orig_ind', 'geometry'])
gcps['geometry'] = [Point(pt) for pt in xy]
gcps['j'] = match_pts[:, 0]
gcps['i'] = match_pts[:, 1]
gcps['oj'] = src_pts[:, 0]
gcps['oi'] = src_pts[:, 1]
gcps['mj'] = dst_pts[:, 0]
gcps['mi'] = dst_pts[:, 1]
gcps['rel_x'] = gt[4] + match_pts[:, 0] * gt[0]
gcps['rel_y'] = gt[5] + match_pts[:, 1] * gt[3]
gcps['orig_ind'] = range(match_pts.shape[0])
gcps['z_corr'] = z_corrs
gcps['z_max'] = peak_corrs

gcps.crs = {'init': 'epsg:{}'.format(lsat_fullres.epsg)}
for i, row in gcps.to_crs(crs=dem.proj4).iterrows():
    gcps.loc[i, 'elevation'] = dem.raster_points([(row.geometry.x, row.geometry.y)], nsize=9, mode='cubic')
    gcps.loc[i, 'el_rel'] = rel_dem.raster_points([(row.rel_x, row.rel_y)], nsize=9, mode='cubic')

# drop any gcps where we don't have a DEM value or a good match
gcps.loc[gcps.z_max < 0, 'z_max'] = np.nan
gcps.dropna(inplace=True)

ref_src = gcps[['j', 'i']].values
ref_dst = gcps[['mj', 'mi']].values

gcps.index = range(gcps.shape[0])  # make sure index corresponds to row we're writing out
gcps['id'] = ['GCP{}'.format(i) for i in range(gcps.shape[0])]
gcps.to_file(os.path.join(out_dir, 'AutoGCPs{}.shp'.format(subscript)))

print('writing AutoGCPs.txt')
write_auto_gcps(gcps, subscript, out_dir, utm_str)

print('converting AutoGCPs.txt to AutoGCPs.xml')
subprocess.Popen(['mm3d', 'GCPConvert', 'AppInFile',
                  os.path.join(out_dir, 'AutoGCPs{}.txt'.format(subscript))]).wait()

print('writing AutoMeasures.txt')
write_auto_mesures(gcps, subscript, out_dir)

print('running get_autogcp_locations.sh to get rough image locations for each point')
subprocess.Popen(['get_autogcp_locations.sh', 'Ori-Relative',
                  os.path.join(out_dir, 'AutoMeasures{}.txt'.format(subscript))] + imlist).wait()

print('searching for points in non-orthorectified images')
E = builder.ElementMaker()
MesureSet = E.SetOfMesureAppuisFlottants()
for im in imlist:
    print(im)
    img = imread(im)
    maxi, maxj = img.shape

    impts = pd.read_csv('Auto-{}.txt'.format(im), sep=' ', names=['j', 'i'])
    impts_nodist = pd.read_csv('NoDist-{}.txt'.format(im), sep=' ', names=['j', 'i'])

    this_im_mes = E.MesureAppuiFlottant1Im(E.NameIm(im))
    for ind, row in impts.iterrows():
        in_im = 0 < row.j < maxj and 0 < row.i < maxi
        in_nd = -200 < impts_nodist.j[ind] < maxj+200 and -200 < impts_nodist.i[ind] < maxi+200
        if in_im and in_nd:
            try:
                testchip, _, _ = imtools.make_template(ortho, (gcps.i[ind], gcps.j[ind]), 20)
                dst_chip, _, _ = imtools.make_template(img, (row.i, row.j), 200)

                # test = imtools.highpass_filter(testchip)
                # dest = imtools.highpass_filter(dst_chip)
                test = np.ma.masked_values(testchip, 0)
                dest = np.ma.masked_values(dst_chip, 0)
                corr_res, this_i, this_j = imtools.find_gcp_match(dest.astype(np.float32), test.astype(np.float32))
                # z_corr = (corr_res.max() - corr_res.min()) / corr_res.std()
                peak_corr = cv2.minMaxLoc(corr_res)[1]

                # pks = peak_local_max(corr_res, min_distance=5, num_peaks=2)
                # this_z_corrs = []
                # for pk in pks:
                #     max_ = corr_res[pk[0], pk[1]]
                #     this_z_corrs.append((max_ - corr_res.mean()) / corr_res.std())
                # dz_corr = max(this_z_corrs) / min(this_z_corrs)
                # z_corr = max(this_z_corrs)

                # if args.corr_thresh < peak_corr or (z_corr > 5 and dz_corr > 1.5):
                if args.corr_thresh < peak_corr:
                    out_i, out_j = this_i - 200 + int(row.i), this_j - 200 + int(row.j)
                    this_mes = E.OneMesureAF1I(E.NamePt('GCP{}'.format(ind)),
                                                E.PtIm('{} {}'.format(out_j, out_i)))
                    this_im_mes.append(this_mes)
            except:
                continue
    MesureSet.append(this_im_mes)

tree = etree.ElementTree(MesureSet)
tree.write(os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript)),
           pretty_print=True, xml_declaration=True, encoding="utf-8")

print('running mm3d GCPBascule to estimate terrain errors')
subprocess.Popen(['mm3d', 'GCPBascule', get_match_pattern(imlist), 'Relative', 'TerrainRelAuto{}'.format(subscript),
                  os.path.join(out_dir, 'AutoGCPs{}.xml'.format(subscript)),
                  os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript))]).wait()

gcps = get_bascule_residuals(os.path.join('Ori-TerrainRelAuto{}'.format(subscript), 'Result-GCP-Bascule.xml'), gcps)

# inds = gcps.residual < 5 * nmad(gcps.residual)
# gcps = gcps.loc[inds]

# gcps.to_file(os.path.join(out_dir, 'AutoGCPs{}.shp'.format(subscript)))
# write_auto_gcps(gcps, subscript, out_dir, utm_str)
# subprocess.Popen(['mm3d', 'GCPConvert', 'AppInFile',
#                   os.path.join(out_dir, 'AutoGCPs{}.txt'.format(subscript))]).wait()
#
# auto_root = ET.parse(os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript))).getroot()
# for im in auto_root.findall('MesureAppuiFlottant1Im'):
#     for pt in im.findall('OneMesureAF1I'):
#         if pt.find('NamePt').text not in gcps.id.values:
#             im.remove(pt)
# # # # save AutoMeasures
# out_xml = ET.ElementTree(auto_root)
# out_xml.write(os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript)),
#               encoding="utf-8", xml_declaration=True)
#
# subprocess.Popen(['mm3d', 'GCPBascule', get_match_pattern(imlist), 'Relative',
#                   'TerrainRelAuto{}'.format(subscript),
#                   os.path.join(out_dir, 'AutoGCPs{}.xml'.format(subscript)),
#                   os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript))]).wait()
#
# gcps = get_bascule_residuals(os.path.join('Ori-TerrainRelAuto{}'.format(subscript), 'Result-GCP-Bascule.xml'), gcps)
#
# out_inds = sliding_window_filter([ortho.shape[1], ortho.shape[0]], gcps, 1500, mindist=800)
# gcps = gcps.loc[out_inds]
# residuals = gcps.residual
# gcp_names = gcps.id
# #
# # out = np.abs(gcps.residual) > 5 * nmad(gcps.residual)
# # gcps = gcps[gcps.id.isin(gcp_names)]
# # out_list = gcp_names[out]
# #
# # gcps.loc[gcps.id.isin(out_list), 'z_corr'] = np.nan
# # gcps.dropna(inplace=True)
# #
# gcps.to_file(os.path.join(out_dir, 'AutoGCPs{}.shp'.format(subscript)))
# write_auto_gcps(gcps, subscript, out_dir, utm_str)
# subprocess.Popen(['mm3d', 'GCPConvert', 'AppInFile',
#                   os.path.join(out_dir, 'AutoGCPs{}.txt'.format(subscript))]).wait()
#
# auto_root = ET.parse(os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript))).getroot()
# for im in auto_root.findall('MesureAppuiFlottant1Im'):
#     for pt in im.findall('OneMesureAF1I'):
#         if pt.find('NamePt').text not in gcps.id.values:
#             im.remove(pt)
# # # save AutoMeasures
# out_xml = ET.ElementTree(auto_root)
# out_xml.write(os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript)),
#               encoding="utf-8", xml_declaration=True)
#
# subprocess.Popen(['mm3d', 'GCPBascule', get_match_pattern(imlist), 'Relative',
#                   'TerrainRelAuto{}'.format(subscript),
#                   os.path.join(out_dir, 'AutoGCPs{}.xml'.format(subscript)),
#                   os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript))]).wait()
#
# subprocess.Popen(['mm3d', 'Campari', get_match_pattern(imlist),
#                   'TerrainRelAuto{}'.format(subscript),
#                   'TerrainFinal{}'.format(subscript),
#                   'GCP=[{},5,{},2]'.format(os.path.join(out_dir, 'AutoGCPs{}.xml'.format(subscript)),
#                                            os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript))),
#                   'SH=Homol', 'AllFree=1']).wait()
#
# gcps = get_bascule_residuals(os.path.join('Ori-TerrainRelAuto{}'.format(subscript), 'Result-GCP-Bascule.xml'), gcps)
# gcps = get_campari_residuals('Ori-TerrainFinal_block0/Residus.xml', gcps)
#
# # gcps = gcps[gcps.camp_res < 5*nmad(gcps.camp_res)]
#
# fig1 = plt.figure(figsize=(7, 5))
# plt.imshow(ortho[::5, ::5], cmap='gray', extent=[0, ortho.shape[1], ortho.shape[0], 0])
# plt.plot(gcps.j, gcps.i, 'r+')
# plt.quiver(gcps.j, gcps.i, gcps.xres, gcps.yres, color='r', scale=1000)
#
# plt.savefig(os.path.join(out_dir, 'relative_gcps{}.png'.format(subscript)), bbox_inches='tight', dpi=200)
# plt.close(fig1)
#
# fig2 = lsat_fullres.display(sfact=10, fig=plt.figure(figsize=(7,5)))
# plt.plot(gcps.geometry.x, gcps.geometry.y, 'r+')
# plt.savefig(os.path.join(out_dir, 'world_gcps{}.png'.format(subscript)), bbox_inches='tight', dpi=200)
# plt.close(fig2)
#
# print('cleaning up.')
# # remove Auto-im.tif.txt, NoDist-im.tif.txt, etc. Ori-Relative-NoDist
# shutil.rmtree('Ori-Relative-NoDist/')
#
# for txtfile in glob('Auto-OIS*.tif.txt') + \
#                glob('NoDist-OIS*.tif.txt'):
#     os.remove(txtfile)
print('end.')