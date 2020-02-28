#!/usr/bin/env python
import argparse
import os
import shutil
from itertools import chain
import numpy as np
import cv2
import gdal
import pandas as pd
import geopandas as gpd
from shapely.geometry.point import Point
from skimage.feature import peak_local_max
from skimage.measure import ransac
from skimage.transform import EuclideanTransform, AffineTransform, warp
from pybob.bob_tools import mkdir_p
from pybob.image_tools import create_mask_from_shapefile
from pybob.GeoImg import GeoImg
import sPyMicMac.image_tools as imtools


def splitter(img, nblocks):
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


def get_dense_keypoints(img, mask, npix=200, return_des=False):
    orb = cv2.ORB_create(50)
    keypts = []
    if return_des:
        descriptors = []

    x_tiles = np.floor(img.shape[1] / npix).astype(int)
    y_tiles = np.floor(img.shape[0] / npix).astype(int)

    split_img = splitter(img, (y_tiles, x_tiles))
    split_msk = splitter(mask, (y_tiles, x_tiles))

    rel_x, rel_y = get_subimg_offsets(split_img, (y_tiles, x_tiles))

    for i, img_ in enumerate(split_img):
        iy, ix = np.unravel_index(i, (y_tiles, x_tiles))

        ox = rel_x[iy, ix]
        oy = rel_y[iy, ix]

        kp, des = orb.detectAndCompute(img_, mask=split_msk[i])
        if return_des:
            for ds in des:
                descriptors.append(ds)

        for p in kp:
            p.pt = p.pt[0] + ox, p.pt[1] + oy
            keypts.append(p)

    if return_des:
        return keypts, descriptors
    else:
        return keypts


def get_initial_transformation(mst, slv, landmask=None):
    mst_lowres = mst.resample(400, method=gdal.GRA_NearestNeighbour)
    slv_lowres = slv.resample(400, method=gdal.GRA_NearestNeighbour)

    slv_lowres_hp = imtools.highpass_filter(slv_lowres.img)
    mst_lowres_hp = imtools.highpass_filter(mst_lowres.img)

    mst_lowres_hp[np.isnan(mst_lowres_hp)] = 0
    slv_lowres_hp[np.isnan(slv_lowres_hp)] = 0
    slv_lowres_hp = np.ma.masked_values(slv_lowres_hp, 0)

    if landmask is not None:
        lmask = create_mask_from_shapefile(mst_lowres, landmask)

    mst_rescale = (256 * (mst_lowres.img - np.nanmin(mst_lowres.img) /
                          (np.nanmax(mst_lowres.img) - np.nanmin(mst_lowres.img))))

    _mask = 255 * np.ones(mst_lowres.img.shape, dtype=np.uint8)

    if landmask is not None:
        _mask[np.logical_or(mst_rescale == 0, ~lmask)] = 0
    else:
        _mask[np.logical_or(np.isnan(mst_rescale), mst_rescale == 0)] = 0

    kpts = get_dense_keypoints(mst_rescale.astype(np.uint8), _mask)

    src_ = np.array([p.pt for p in kpts])

    z_corrs = []
    peak_corrs = []
    match_pts = []
    res_imgs = []

    for pt in src_:
        try:
            chip, _, _ = imtools.make_template(mst_lowres_hp, (pt[1], pt[0]), 25)
            chip = np.ma.masked_values(chip, 0)

            corr_res, this_i, this_j = imtools.find_gcp_match(slv_lowres_hp.astype(np.float32), chip.astype(np.float32))

            peak_corr = cv2.minMaxLoc(corr_res)[1]

            pks = peak_local_max(corr_res, min_distance=5, num_peaks=2)
            this_z_corrs = []
            for pk in pks:
                max_ = corr_res[pk[0], pk[1]]
                this_z_corrs.append((max_ - corr_res.mean()) / corr_res.std())

            z_corrs.append(max(this_z_corrs))
            peak_corrs.append(peak_corr)
            match_pts.append([this_j, this_i])
            res_imgs.append(corr_res)
        except:
            z_corrs.append(-1)
            peak_corrs.append(-1)
            match_pts.append([-1, -1])
            res_imgs.append(-1)
            continue

    dst_pts = np.array(match_pts)

    lowres_gcps = pd.DataFrame()
    lowres_gcps['pk_corr'] = peak_corrs
    lowres_gcps['z_corr'] = z_corrs
    lowres_gcps['match_j'] = dst_pts[:, 0]
    lowres_gcps['match_i'] = dst_pts[:, 1]
    lowres_gcps['src_j'] = src_[:, 0]
    lowres_gcps['src_i'] = src_[:, 1]

    best_lowres = lowres_gcps[lowres_gcps['pk_corr'] > lowres_gcps['pk_corr'].quantile(0.75)]

    dst_pts = best_lowres[['match_j', 'match_i']].values
    src_pts = best_lowres[['src_j', 'src_i']].values

    dst_scale = (slv_lowres.dx / slv.dx)
    src_scale = (mst_lowres.dx / mst.dx)

    Minit, inliers = ransac((dst_scale * dst_pts, src_scale * src_pts), EuclideanTransform, min_samples=5,
                            residual_threshold=3, max_trials=1000)
    print('{} points used to find initial transformation'.format(np.count_nonzero(inliers)))

    return Minit


def qa_mask(img):
    return np.logical_not(np.logical_or.reduce((img == 672, img == 676, img == 680, img == 684,
                                                img == 1696, img == 1700, img == 1704, img == 1708)))


def get_mask(mst, slv, tfm, _args):
    mask = 255 * np.ones(slv.img.shape, dtype=np.uint8)

    if _args.quality_mask is not None:
        slv_qa = GeoImg(_args.quality_mask)
        quality_mask = qa_mask(slv_qa.img)
        mask[quality_mask] = 0

    if _args.landmask is not None:
        lm = create_mask_from_shapefile(mst, _args.landmask)
        lm_tfm = warp(lm, tfm, output_shape=slv.img.shape, order=0)
        mask[lm_tfm == 0] = 0

    if _args.glacmask is not None:
        gm = create_mask_from_shapefile(mst, _args.glacmask, buffer=200)
        gm_tfm = warp(gm, tfm, output_shape=slv.img.shape, order=0)
        mask[gm_tfm > 0] = 0

    return mask


def _argparser():
    parser = argparse.ArgumentParser(description="Automatically register a Landsat L1GS image to a Landsat L1TP image",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('master', action='store', type=str, help='non-referenced orthophoto mosaic')
    parser.add_argument('slave', action='store', type=str, help='georeferenced satellite image')
    parser.add_argument('-q', '--quality_mask', action='store', type=str, default=None,
                        help='BQA band for the L1GS image (recommended)')
    parser.add_argument('-dem', action='store', type=str, default=None,
                        help='digital elevation model')
    parser.add_argument('-glacmask', action='store', type=str, default=None,
                        help='path to shapefile of glacier outlines')
    parser.add_argument('-landmask', action='store', type=str, default=None,
                        help='path to shapefile of land outlines')

    return parser


def main():
    parser_ = _argparser()
    args = parser_.parse_args()

    mst_fullres = GeoImg(args.master)
    slv_fullres = GeoImg(args.slave)

    Minit = get_initial_transformation(mst_fullres, slv_fullres, landmask=args.landmask)

    rough_tfm = warp(mst_fullres.img, Minit, output_shape=slv_fullres.img.shape, preserve_range=True)
    rough_tfm[np.isnan(rough_tfm)] = 0

    mask = get_mask(mst_fullres, slv_fullres, Minit, args)
    mask[rough_tfm == 0] = 0

    # for each of these pairs (src, dst), find the precise subpixel match (or not...)
    match_pts = []
    z_corrs = []
    peak_corrs = []
    res_imgs = []

    jj = np.arange(0, slv_fullres.img.shape[1], 50)
    ii = np.arange(0, slv_fullres.img.shape[0], 50)

    I, J = np.meshgrid(ii, jj)
    search_pts = np.hstack([J.reshape(-1, 1), I.reshape(-1, 1)])

    for pt in search_pts:
        if mask[pt[1], pt[0]] == 0:
            match_pts.append((-1, -1))
            z_corrs.append(np.nan)
            peak_corrs.append(np.nan)
            res_imgs.append(np.nan)
            continue

        try:
            testchip, _, _ = imtools.make_template(rough_tfm, (pt[1], pt[0]), 40)
            dst_chip, _, _ = imtools.make_template(slv_fullres.img, (pt[1], pt[0]), 80)
            dst_chip[np.isnan(dst_chip)] = 0

            test = np.ma.masked_values(imtools.highpass_filter(testchip), 0)
            dest = np.ma.masked_values(imtools.highpass_filter(dst_chip), 0)

            corr_res, this_i, this_j = imtools.find_gcp_match(dest.astype(np.float32), test.astype(np.float32))
            peak_corr = cv2.minMaxLoc(corr_res)[1]

            pks = peak_local_max(corr_res, min_distance=5, num_peaks=2)
            this_z_corrs = []
            for pk in pks:
                max_ = corr_res[pk[0], pk[1]]
                this_z_corrs.append((max_ - corr_res.mean()) / corr_res.std())
            dz_corr = max(this_z_corrs) / min(this_z_corrs)
            z_corr = max(this_z_corrs)

            # if the correlation peak is very high, or very unique, add it as a match
            out_i, out_j = this_i - 80 + pt[1], this_j - 80 + pt[0]
            z_corrs.append(z_corr)
            peak_corrs.append(peak_corr)
            match_pts.append([out_j, out_i])
            res_imgs.append(corr_res)
        except:
            match_pts.append((-1, -1))
            z_corrs.append(np.nan)
            peak_corrs.append(np.nan)
            res_imgs.append(np.nan)

    _src = np.dot(Minit.params, np.hstack([search_pts,
                                           np.ones(search_pts[:, 0].shape).reshape(-1, 1)]).T).T[:, :2]
    _dst = np.array(match_pts)

    xy = np.array([mst_fullres.ij2xy((pt[1], pt[0])) for pt in _src]).reshape(-1, 2)

    gcps = gpd.GeoDataFrame()
    gcps['geometry'] = [Point(pt) for pt in xy]
    gcps['pk_corr'] = peak_corrs
    gcps['z_corr'] = z_corrs
    gcps['match_j'] = _dst[:, 0]
    gcps['match_i'] = _dst[:, 1]
    gcps['src_j'] = search_pts[:, 0]
    gcps['src_i'] = search_pts[:, 1]
    gcps['dj'] = gcps['src_j'] - gcps['match_j']
    gcps['di'] = gcps['src_i'] - gcps['match_i']
    gcps['elevation'] = 0
    gcps.crs = mst_fullres.proj4

    gcps.dropna(inplace=True)

    if args.dem is not None:
        dem = GeoImg(args.dem)
        for i, row in gcps.to_crs(crs=dem.proj4).iterrows():
            gcps.loc[i, 'elevation'] = dem.raster_points([(row.geometry.x, row.geometry.y)], nsize=9, mode='cubic')

    best = gcps[np.logical_and(gcps.pk_corr > 0.2, gcps.z_corr > 8)]

    Mfin, inliers_fin = ransac((best[['match_j', 'match_i']].values, best[['src_j', 'src_i']].values), AffineTransform,
                               min_samples=10, residual_threshold=1, max_trials=1000)
    best = best[inliers_fin]

    gcp_list = []
    outname = os.path.splitext(os.path.basename(args.slave))[0]
    with open('{}_gcps.txt'.format(outname), 'w') as f:
        for i, row in best.iterrows():
            gcp_list.append(gdal.GCP(row.geometry.x, row.geometry.y, row.elevation, row.match_j, row.match_i))
            print(row.geometry.x, row.geometry.y, row.elevation, row.match_j, row.match_i, file=f)

    shutil.copy(args.slave, 'tmp.tif')
    in_ds = gdal.Open('tmp.tif', gdal.GA_Update)
    gcp_wkt = mst_fullres.proj_wkt
    in_ds.SetGCPs(gcp_list, gcp_wkt)

    print('warping image to new geometry')
    mkdir_p('warped')
    gdal.Warp(os.path.join('warped', os.path.basename(args.slave)),
              in_ds,
              xRes=slv_fullres.dx,
              yRes=slv_fullres.dy,
              resampleAlg=gdal.GRA_Bilinear)

    print('cleaning up.')
    os.remove('tmp.tif')


if __name__ == "__main__":
    main()
