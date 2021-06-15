"""
sPyMicMac.image is a collection of tools for working with aerial images.
"""
import os
import shutil
from glob import glob
import multiprocessing as mp
from functools import partial
from itertools import chain
import cv2
import gdal
import matplotlib.pyplot as plt
import pandas as pd
import lxml.etree as etree
import lxml.builder as builder
from rtree import index
from skimage.filters import rank, median, sobel_v, sobel_h
from skimage import exposure, transform, morphology, io
from skimage.morphology import binary_dilation, disk
from skimage.measure import ransac
from skimage.feature import peak_local_max
from skimage.transform import match_histograms, warp, AffineTransform, EuclideanTransform
from scipy.interpolate import RectBivariateSpline as RBS
from scipy import ndimage
import numpy as np
from shapely.ops import cascaded_union
from shapely.geometry import LineString
import geopandas as gpd
# import pyvips
from llc import jit_filter_function
from pybob.image_tools import match_hist, reshape_geoimg, create_mask_from_shapefile, nanmedian_filter
from pybob.bob_tools import mkdir_p
from pymmaster.mmaster_tools import orient_footprint
from sPyMicMac.micmac import get_im_meas, parse_im_meas


######################################################################################################################
# image filtering tools
######################################################################################################################
@jit_filter_function
def nanstd(a):
    return np.nanstd(a)


def cross_template(shape, width=3):
    if isinstance(shape, int):
        rows = shape
        cols = shape
    else:
        rows, cols = shape
    half_r = int((rows - 1) / 2)
    half_c = int((cols - 1) / 2)
    half_w = int((width - 1) / 2)

    cross = np.zeros((rows, cols))
    cross[half_r - half_w - 1:half_r + half_w + 2:width + 1, :] = 2
    cross[:, half_c - half_w - 1:half_c + half_w + 2:width + 1] = 2

    cross[half_r - half_w:half_r + half_w + 1, :] = 1
    cross[:, half_c - half_w:half_c + half_w + 1] = 1
    return cross


def make_template(img, pt, half_size):
    nrows, ncols = img.shape
    row, col = np.round(pt).astype(int)
    left_col = max(col - half_size, 0)
    right_col = min(col + half_size, ncols)
    top_row = max(row - half_size, 0)
    bot_row = min(row + half_size, nrows)
    row_inds = [row - top_row, bot_row - row]
    col_inds = [col - left_col, right_col - col]
    template = img[top_row:bot_row + 1, left_col:right_col + 1].copy()
    return template, row_inds, col_inds


def find_match(img, template):
    img_eq = rank.equalize(img, selem=morphology.disk(100))
    # res = cross_filter(img_eq, template)
    res = cv2.matchTemplate(img_eq, template, cv2.TM_CCORR_NORMED)
    i_off = (img.shape[0] - res.shape[0]) / 2
    j_off = (img.shape[1] - res.shape[1]) / 2
    minval, _, minloc, _ = cv2.minMaxLoc(res)
    # maxj, maxi = maxloc
    minj, mini = minloc
    sp_delx, sp_dely = get_subpixel(res)
    # sp_delx, sp_dely = 0, 0
    return res, mini + i_off + sp_dely, minj + j_off + sp_delx


def get_subpixel(res, how='min'):
    assert how in ['min', 'max'], "have to choose min or max"

    mgx, mgy = np.meshgrid(np.arange(-1, 1.01, 0.1), np.arange(-1, 1.01, 0.1), indexing='xy')  # sub-pixel mesh

    if how == 'min':
        peakval, _, peakloc, _ = cv2.minMaxLoc(res)
        mml_ind = 2
    else:
        _, peakval, _, peakloc = cv2.minMaxLoc(res)
        mml_ind = 3

    rbs_halfsize = 3  # size of peak area used for spline for subpixel peak loc
    rbs_order = 4  # polynomial order for subpixel rbs interpolation of peak location

    if ((np.array([n - rbs_halfsize for n in peakloc]) >= np.array([0, 0])).all()
            & (np.array([(n + rbs_halfsize) for n in peakloc]) < np.array(list(res.shape))).all()):
        rbs_p = RBS(range(-rbs_halfsize, rbs_halfsize + 1), range(-rbs_halfsize, rbs_halfsize + 1),
                    res[(peakloc[1] - rbs_halfsize):(peakloc[1] + rbs_halfsize + 1),
                    (peakloc[0] - rbs_halfsize):(peakloc[0] + rbs_halfsize + 1)],
                    kx=rbs_order, ky=rbs_order)

        b = rbs_p.ev(mgx.flatten(), mgy.flatten())
        mml = cv2.minMaxLoc(b.reshape(21, 21))
        # mgx,mgy: meshgrid x,y of common area
        # sp_delx,sp_dely: subpixel delx,dely
        sp_delx = mgx[mml[mml_ind][0], mml[mml_ind][1]]
        sp_dely = mgy[mml[mml_ind][0], mml[mml_ind][1]]
    else:
        sp_delx = 0.0
        sp_dely = 0.0
    return sp_delx, sp_dely


def highpass_filter(img):
    """

    :param img:
    :return:
    """
    v = img.copy()
    v[np.isnan(img)] = 0
    vv = ndimage.gaussian_filter(v, 3)

    w = 0 * img.copy() + 1
    w[np.isnan(img)] = 0
    ww = ndimage.gaussian_filter(w, 3)

    tmplow = vv / ww
    tmphi = img - tmplow
    return tmphi


def splitter(img, nblocks, overlap=0):
    new_width = int(np.floor(img.shape[1]/nblocks[1]))
    new_height = int(np.floor(img.shape[0]/nblocks[0]))

    simages = []
    top_inds = []
    left_inds = []
    for i in range(nblocks[0]):
        this_col = []
        this_top = []
        this_left = []

        for j in range(nblocks[1]):
            lind = max(0, j*new_width-overlap)
            rind = min(img.shape[1], (j+1)*new_width + overlap)
            tind = max(0, i*new_height - overlap)
            bind = min(img.shape[0], (i+1)*new_height + overlap)
            this_col.append(img[tind:bind, lind:rind].copy())
            this_top.append(tind)
            this_left.append(lind)
        simages.append(this_col)
        top_inds.append(this_top)
        left_inds.append(this_left)

    return list(chain.from_iterable(simages)), list(chain.from_iterable(top_inds)), list(chain.from_iterable(left_inds))


def get_subimg_offsets(split, shape, overlap=0):
    ims_x = np.array([s.shape[1] - 2 * overlap for s in split])
    ims_y = np.array([s.shape[0] - 2 * overlap for s in split])

    rel_x = np.cumsum(ims_x.reshape(shape), axis=1)
    rel_y = np.cumsum(ims_y.reshape(shape), axis=0)

    rel_x = np.concatenate((np.zeros((shape[0], 1)), rel_x[:, :-1]), axis=1)
    rel_y = np.concatenate((np.zeros((1, shape[1])), rel_y[:-1, :]), axis=0)

    return rel_x.astype(int), rel_y.astype(int)


def stretch_image(img, scale=(0,1), mult=255, outtype=np.uint8, mask=None):
    """

    :param img:
    :param scale:
    :param mult:
    :param outtype:
    :param mask:
    :return:
    """
    if mask is None:
        maxval = np.nanquantile(img, max(scale))
        minval = np.nanquantile(img, min(scale))
    else:
        maxval = np.nanquantile(img[mask], max(scale))
        minval = np.nanquantile(img[mask], min(scale))

    img[img > maxval] = maxval
    img[img < minval] = minval

    return (mult * (img - minval) / (maxval - minval)).astype(outtype)


def contrast_enhance(fn_img, mask_value=None, qmin=0.02, qmax=0.98, gamma=1.25):
    """

    :param fn_img:
    :param mask_value:
    :param qmin:
    :param qmax:
    :param gamma:
    :return:
    """
    img = io.imread(fn_img)
    if mask_value is not None:
        img = img.astype(np.float32)
        img[img == mask_value] = np.nan

    filt = nanmedian_filter(img, footprint=disk(3))

    stretch = stretch_image(filt, scale=(qmin, qmax))
    gamma = exposure.adjust_gamma(stretch, gamma=gamma)
    return gamma


def make_binary_mask(img, erode=0, mask_value=0):
    """

    :param img:
    :param erode:
    :param mask_value:
    :return:
    """
    _mask = 255 * np.ones(img.shape, dtype=np.uint8)
    if np.isfinite(mask_value):
        _mask[img == mask_value] = 0
    else:
        _mask[np.isnan(img)] = 0

    if erode > 0:
        erode_mask = morphology.binary_erosion(_mask, selem=morphology.disk(erode))
        _mask[~erode_mask] = 0

    return _mask


def balance_image(img):
    """

    :param img:
    :return:
    """
    img_eq = (255 * exposure.equalize_adapthist(img)).astype(np.uint8)
    img_filt = median(img_eq, selem=disk(1))
    return img_filt


######################################################################################################################
# GCP matching tools
######################################################################################################################
def keypoint_grid(img, spacing=25, size=10):
    _j = np.arange(0, img.shape[1], spacing)
    _i = np.arange(0, img.shape[0], spacing)
    _gridj, _gridi = np.meshgrid(_j, _i)
    _grid = np.concatenate((_gridj.reshape(-1, 1), _gridi.reshape(-1, 1)), axis=1)

    return [cv2.KeyPoint(pt[0], pt[1], size) for pt in _grid]


def get_dense_keypoints(img, mask, npix=100, nblocks=None, return_des=False):
    """

    :param img:
    :param mask:
    :param npix:
    :param nblocks:
    :param return_des:
    :return:
    """
    orb = cv2.ORB_create()
    keypts = []
    if return_des:
        descriptors = []

    if nblocks is None:
        x_tiles = np.floor(img.shape[1] / npix).astype(int)
        y_tiles = np.floor(img.shape[0] / npix).astype(int)
    else:
        x_tiles = nblocks
        y_tiles = nblocks

    olap = int(max(0.25 * img.shape[1]/x_tiles, 0.25 * img.shape[0]/y_tiles))

    split_img, oy, ox = splitter(img, (y_tiles, x_tiles), overlap=olap)
    split_msk, _, _ = splitter(mask, (y_tiles, x_tiles), overlap=olap)

    # rel_x, rel_y = get_subimg_offsets(split_img, (y_tiles, x_tiles), overlap=olap)

    for i, img_ in enumerate(split_img):

        kp, des = orb.detectAndCompute(img_, mask=split_msk[i])
        if return_des:
            if des is not None:
                for ds in des:
                    descriptors.append(ds)

        for p in kp:
            p.pt = p.pt[0] + ox[i], p.pt[1] + oy[i]
            keypts.append(p)

    if return_des:
        return keypts, np.array(descriptors)
    else:
        return keypts


def get_footprint_overlap(fprints):
    """

    :param fprints:
    :return:
    """
    if fprints.shape[0] == 1:
        return fprints.geometry.values[0]

    idx = index.Index()

    for pos, row in fprints.iterrows():
        idx.insert(pos, row['geometry'].bounds)

    intersections = []
    for poly in fprints['geometry']:
        merged = cascaded_union([fprints.loc[pos, 'geometry'] for pos in idx.intersection(poly.bounds) if fprints.loc[pos, 'geometry'] != poly])
        intersections.append(poly.intersection(merged))

    intersection = cascaded_union(intersections)

    # return intersection.minimum_rotated_rectangle
    return intersection


def get_footprint_mask(shpfile, geoimg, filelist, fprint_out=False):
    """

    :param shpfile:
    :param geoimg:
    :param filelist:
    :param fprint_out:
    :return:
    """
    imlist = [im.split('OIS-Reech_')[-1].split('.tif')[0] for im in filelist]
    if isinstance(shpfile, str):
        footprints_shp = gpd.read_file(shpfile)
        fp = footprints_shp[footprints_shp.ID.isin(imlist)].copy()
    else:
        fp = shpfile[shpfile.ID.isin(imlist)].copy()

    fprint = get_footprint_overlap(fp.to_crs(geoimg.proj4))

    tmp_gdf = gpd.GeoDataFrame(columns=['geometry'])
    tmp_gdf.loc[0, 'geometry'] = fprint
    tmp_gdf.crs = geoimg.proj4
    tmp_gdf.to_file('tmp_fprint.shp')

    maskout = create_mask_from_shapefile(geoimg, 'tmp_fprint.shp')

    for f in glob('tmp_fprint.*'):
        os.remove(f)
    if fprint_out:
        return maskout, fprint
    else:
        return maskout


def get_rough_geotransform(img1, img2, landmask=None, footmask=None,
                           stretch=(0.05, 0.95), equalize=False):
    """

    :param img1:
    :param img2:
    :param landmask:
    :param footmask:
    :param stretch:
    :param equalize:
    :return:
    """
    if equalize:
        img2_eq = (255 * exposure.equalize_adapthist(img2.img.astype(np.uint16),
                                                     clip_limit=0.03)).astype(np.uint8)
    else:
        img2_eq = stretch_image(img2.img, scale=stretch)

    img1_mask = make_binary_mask(img1, erode=3, mask_value=0)

    img2_mask = make_binary_mask(img2.img, erode=3, mask_value=np.nan)

    if landmask is not None:
        # lm = create_mask_from_shapefile(img2, landmask)
        img2_mask[~landmask] = 0

    if footmask is not None:
        # mask_ = get_footprint_mask(footmask, img2, imlist)
        img2_mask[~footmask] = 0

    kp, des, matches = get_matches(img1.astype(np.uint8), img2_eq,
                                   mask1=img1_mask, mask2=img2_mask, dense=True)
    print('{} matches found.'.format(len(matches)))

    src_pts = np.array([kp[0][m.queryIdx].pt for m in matches])
    dst_pts = np.array([kp[1][m.trainIdx].pt for m in matches])

    Minit, inliers = ransac((dst_pts, src_pts), AffineTransform,
                            min_samples=5, residual_threshold=10, max_trials=1000)
    print('{} matches used for initial transformation'.format(np.count_nonzero(inliers)))

    img1_tfm = warp(img1, Minit, output_shape=img2_eq.shape, preserve_range=True)

    return img1_tfm, Minit, (dst_pts, src_pts, inliers)


def get_initial_transformation(img1, img2, pRes=800, landmask=None, footmask=None, imlist=None):
    """

    :param img1:
    :param img2:
    :param pRes:
    :param landmask:
    :param footmask:
    :param imlist:
    :return:
    """
    im2_lowres = reshape_geoimg(img2, pRes, pRes)

    im2_eq = match_hist(im2_lowres.img, np.array(img1))
    im1_mask = 255 * np.ones(img1.shape, dtype=np.uint8)
    im1_mask[img1 == 0] = 0  # nodata from ortho

    im2_mask = 255 * np.ones(im2_eq.shape, dtype=np.uint8)
    im2_mask[im2_eq == 0] = 0

    if imlist is None:
        imlist = glob('OIS*.tif')

    if landmask is not None:
        mask_ = create_mask_from_shapefile(im2_lowres, landmask)
        im2_mask[~mask_] = 0
    if footmask is not None:
        mask_ = get_footprint_mask(footmask, im2_lowres, imlist)
        im2_mask[~mask_] = 0

    kp, des, matches = get_matches(img1, im2_eq, mask1=im1_mask, mask2=im2_mask)
    src_pts = np.array([kp[0][m.queryIdx].pt for m in matches])
    dst_pts = np.array([kp[1][m.trainIdx].pt for m in matches])
    # aff_matrix, good_mask = cv2.estimateAffine2D(src_pts, dst_pts, ransacReprojThreshold=25)
    Mout, inliers = ransac((dst_pts, src_pts), EuclideanTransform,
                           min_samples=5, residual_threshold=2, max_trials=1000)
    # check that the transformation was successful by correlating the two images.
    # im1_tfm = cv2.warpAffine(img1, aff_matrix, (im2_lowres.img.shape[1], im2_lowres.img.shape[0]))
    im1_tfm = warp(img1, Mout, output_shape=im2_lowres.img.shape, preserve_range=True)
    im1_pad = np.zeros(np.array(im1_tfm.shape) + 2, dtype=np.uint8)
    im1_pad[1:-1, 1:-1] = im1_tfm
    res = cv2.matchTemplate(np.ma.masked_values(im2_eq, 0),
                            np.ma.masked_values(im1_pad, 0),
                            cv2.TM_CCORR_NORMED)
    print(res[1, 1])
    success = res[1, 1] > 0.5

    return Mout, success, im2_eq.shape


def get_matches(img1, img2, mask1=None, mask2=None, dense=False):

    if dense:
        if np.any(np.array(img1.shape) < 100) or np.any(np.array(img2.shape) < 100):
            kp1, des1 = get_dense_keypoints(img1.astype(np.uint8), mask1, nblocks=1, return_des=True)
            kp2, des2 = get_dense_keypoints(img2.astype(np.uint8), mask2, nblocks=1, return_des=True)
        else:
            kp1, des1 = get_dense_keypoints(img1.astype(np.uint8), mask1, return_des=True)
            kp2, des2 = get_dense_keypoints(img2.astype(np.uint8), mask2, return_des=True)
    else:
        orb = cv2.ORB_create()
        kp1, des1 = orb.detectAndCompute(img1.astype(np.uint8), mask=mask1)
        kp2, des2 = orb.detectAndCompute(img2.astype(np.uint8), mask=mask2)

    flann_idx = 6
    index_params = dict(algorithm=flann_idx, table_number=12, key_size=20, multi_probe_level=2)
    search_params = dict(checks=10000)

    flann = cv2.FlannBasedMatcher(index_params, search_params)
    raw_matches = flann.knnMatch(des1, des2, k=2)
    matches = []
    for m in raw_matches:
        if len(m) == 2 and m[0].distance < m[1].distance * 0.75:
            matches.append(m[0])

    return (kp1, kp2), (des1, des2), matches


def find_gcp_match(img, template, method=cv2.TM_CCORR_NORMED):
    res = cv2.matchTemplate(img, template, method)
    i_off = (img.shape[0] - res.shape[0]) / 2
    j_off = (img.shape[1] - res.shape[1]) / 2
    _, maxval, _, maxloc = cv2.minMaxLoc(res)
    maxj, maxi = maxloc
    try:
        sp_delx, sp_dely = get_subpixel(res, how='max')
    except ValueError as e:
        sp_delx = 0
        sp_dely = 0

    return res, maxi + i_off + sp_dely, maxj + j_off + sp_delx


def find_grid_matches(tfm_img, refgeo, mask, initM=None, spacing=200, srcwin=60, dstwin=600):
    """

    :param tfm_img:
    :param refgeo:
    :param mask:
    :param initM:
    :param spacing:
    :param srcwin:
    :param dstwin:
    :return:
    """
    match_pts = []
    z_corrs = []
    peak_corrs = []
    res_imgs = []

    jj = np.arange(srcwin, spacing * np.ceil((refgeo.img.shape[1]-srcwin) / spacing) + 1, spacing).astype(int)
    ii = np.arange(srcwin, spacing * np.ceil((refgeo.img.shape[0]-srcwin) / spacing) + 1, spacing).astype(int)

    search_pts = []

    for _i in ii:
        for _j in jj:
            search_pts.append((_j, _i))
            # for pt in search_pts:
            # if mask[pt[1], pt[0]] == 0:
            submask, _, _ = make_template(mask, (_i, _j), srcwin)
            if np.count_nonzero(submask) / submask.size < 0.05:
                match_pts.append((-1, -1))
                z_corrs.append(np.nan)
                peak_corrs.append(np.nan)
                res_imgs.append(np.nan)
                continue
            try:
                testchip, _, _ = make_template(refgeo.img, (_i, _j), srcwin)
                dst_chip, _, _ = make_template(tfm_img, (_i, _j), dstwin)

                testchip[np.isnan(testchip)] = 0
                dst_chip[np.isnan(dst_chip)] = 0

                test = highpass_filter(testchip)
                dest = highpass_filter(dst_chip)

                testmask = binary_dilation(testchip == 0, selem=disk(8))
                destmask = binary_dilation(dst_chip == 0, selem=disk(8))

                test[testmask] = np.random.rand(test.shape[0], test.shape[1])[testmask]
                dest[destmask] = np.random.rand(dest.shape[0], dest.shape[1])[destmask]

                corr_res, this_i, this_j = find_gcp_match(dest.astype(np.float32), test.astype(np.float32))
                peak_corr = cv2.minMaxLoc(corr_res)[1]

                pks = peak_local_max(corr_res, min_distance=5, num_peaks=2)
                this_z_corrs = []
                for pk in pks:
                    max_ = corr_res[pk[0], pk[1]]
                    this_z_corrs.append((max_ - corr_res.mean()) / corr_res.std())
                dz_corr = max(this_z_corrs) / min(this_z_corrs)
                z_corr = max(this_z_corrs)

                # if the correlation peak is very high, or very unique, add it as a match
                out_i, out_j = this_i - dstwin + _i, this_j - dstwin + _j
                z_corrs.append(z_corr)
                peak_corrs.append(peak_corr)
                match_pts.append([out_j, out_i])
                res_imgs.append(corr_res)
            except:
                match_pts.append((-1, -1))
                z_corrs.append(np.nan)
                peak_corrs.append(np.nan)
                res_imgs.append(np.nan)

    search_pts = np.array(search_pts)
    _dst = np.array(match_pts)

    gcps = gpd.GeoDataFrame()
    gcps['pk_corr'] = peak_corrs
    gcps['z_corr'] = z_corrs
    gcps['match_j'] = _dst[:, 0]  # points matched in transformed image
    gcps['match_i'] = _dst[:, 1]

    if initM is not None:
        # have to find the initial, which means back-transforming match_pts with Minit_full
        # _src = np.dot(initM.params, np.hstack([_dst, np.ones(_dst[:, 0].shape).reshape(-1, 1)]).T).T[:, :2]
        _src = initM(_dst)
        gcps['orig_j'] = _src[:, 0]  # this should be the back-transformed match_pts
        gcps['orig_i'] = _src[:, 1]

    gcps['search_j'] = search_pts[:, 0]
    gcps['search_i'] = search_pts[:, 1]
    gcps['dj'] = gcps['search_j'] - gcps['match_j']
    gcps['di'] = gcps['search_i'] - gcps['match_i']

    gcps.dropna(inplace=True)

    return gcps


def transform_from_fprint(img, geo, fprint):
    """

    :param img:
    :param geo:
    :param fprint:
    :return:
    """
    oprint = orient_footprint(fprint)
    h, w = img.shape
    src_pts = np.array([[0, 0], [0, h-1], [w-1, h-1], [w-1, 0]]).reshape(-1, 2)
    x, y = oprint.boundary.xy
    ij = [geo.xy2ij(pt) for pt in zip(x, y)]
    i = [pt[0] for pt in ij]
    j = [pt[1] for pt in ij]
    inds = np.array([1, 0, 3, 2]).reshape(-1, 1)
    dst_pts = np.array(list(zip(np.array(j)[inds], np.array(i)[inds]))).reshape(-1, 2)

    return transform.estimate_transform('euclidean', dst_pts, src_pts)


def find_cross(img, pt, cross, tsize=300):
    subimg, _, _ = make_template(img, pt, tsize)
    res, this_i, this_j = find_match(subimg, cross)
    inv_res = res.max() - res

    pks = peak_local_max(inv_res, min_distance=5, num_peaks=2)
    this_z_corrs = []
    for pk in pks:
        max_ = inv_res[pk[0], pk[1]]
        this_z_corrs.append((max_ - inv_res.mean()) / inv_res.std())
    if max(this_z_corrs) > 4 and max(this_z_corrs)/min(this_z_corrs) > 1.15:
        return this_i-tsize+pt[0], this_j-tsize+pt[1]
    else:
        return np.nan, np.nan


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


def _moving_average(a, n=5):
    ret = np.cumsum(a)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


def get_perp(corners):
    a = np.array([corners[1][0] - corners[0][0],
                  corners[1][1] - corners[0][1]])
    a = a / np.linalg.norm(a)
    return np.array([-a[1], a[0]])


def find_border(img, cross, tsize=300):
    img_mask = np.zeros(img.shape)
    img_mask[np.logical_or(img < 25, median(img, selem=disk(16)) < 25)] = 1

    vert = np.count_nonzero(sobel_v(img_mask)**2 > 0.8, axis=0)
    hori = np.count_nonzero(sobel_h(img_mask)**2 > 0.8, axis=1)

    vpeaks = peak_local_max(_moving_average(vert), min_distance=6000, num_peaks=1, exclude_border=10)
    hpeaks = peak_local_max(_moving_average(hori), min_distance=3000, num_peaks=2, exclude_border=10)

    top, bot = 10 * hpeaks.min(), 10 * hpeaks.max()
    edge = 10 * vpeaks.min()
    if edge < 5000:
        search_corners = [(bot-625, edge+180), (top+625, edge+180)]
    else:
        search_corners = [(bot - 625, edge-180), (top + 625, edge-180)]

    grid_corners = []
    for c in search_corners:
        _i, _j = find_cross(img, c, cross, tsize)
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

    return left_edge, right_edge, top, bot


def _search_grid(ledge, redge, left, right):
    i_grid = []
    j_grid = []
    search_pts = []

    for ii in range(0, 23):
        this_left = ledge.interpolate(ii * left.length / 22)
        this_right = redge.interpolate(ii * right.length / 22)
        this_line = LineString([this_left, this_right])
        for jj in range(0, 24):
            this_pt = this_line.interpolate(jj * this_line.length / 23)
            i_grid.append(ii)
            j_grid.append(jj)
            search_pts.append((this_pt.y, this_pt.x))

    return i_grid, j_grid, search_pts


def find_reseau_grid(fn_img, csize=361, tsize=300, nproc=1, return_val=False):
    """

    :param fn_img:
    :param csize:
    :param tsize:
    :param scanres:
    :param scanres_y:
    :param nproc:
    :param return_val:
    :return:
    """
    print('Reading {}'.format(fn_img))
    img = io.imread(fn_img)

    print('Image read.')
    tmp_cross = cross_template(csize, width=3)
    cross = np.ones((csize, csize)).astype(np.uint8)
    cross[tmp_cross == 1] = 255

    fig = plt.figure(figsize=(7, 12))
    ax = fig.add_subplot(111)
    ax.imshow(img[::10, ::10], cmap='gray', extent=[0, img.shape[1], img.shape[0], 0])

    left_edge, right_edge, top, bot = find_border(img)

    II, JJ, search_grid = _search_grid(left_edge, right_edge)
    matched_grid = []

    print('Finding grid points in {}...'.format(fn_img))
    if nproc > 1:
        subimgs = []
        std_devs = []
        means = []
        for loc in search_grid:
            subimg, _, _ = make_template(img, loc, tsize)
            subimgs.append(subimg)
        pool = mp.Pool(nproc)
        outputs = pool.map(partial(find_match, template=cross), subimgs)
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

                if max(this_z_corrs) > 5 and max(this_z_corrs) / min(this_z_corrs) > 1.15:
                    rel_i = this_i - tsize
                    rel_j = this_j - tsize
                    i, j = search_grid[n]
                    matched_grid.append((rel_i + i, rel_j + j))
                else:
                    matched_grid.append((np.nan, np.nan))
            else:
                matched_grid.append((np.nan, np.nan))
    else:
        for loc in search_grid:
            _j, _i = find_cross(img, loc, cross, tsize=tsize)
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

    model = AffineTransform()
    model.estimate(gcps_df[['search_j', 'search_i']].values, gcps_df[['match_j', 'match_i']].values)
    dst = model(gcps_df[['search_j', 'search_i']].values)

    nomatch = np.isnan(gcps_df.match_j)

    ux = lsq_fit(gcps_df.search_j.values, gcps_df.search_i.values, gcps_df.match_j.values)
    uy = lsq_fit(gcps_df.search_j.values, gcps_df.search_i.values, gcps_df.match_i.values)

    gcps_df.loc[nomatch, 'match_j'] = ux[nomatch]
    gcps_df.loc[nomatch, 'match_i'] = uy[nomatch]

    gcps_df['dj'] = gcps_df['match_j'] - dst[:, 0]
    gcps_df['di'] = gcps_df['match_i'] - dst[:, 1]

    gcps_df['im_row'] = gcps_df['match_i']
    gcps_df['im_col'] = gcps_df['match_j']

    print('Grid points found.')
    mkdir_p('match_imgs')

    ax.quiver(gcps_df.search_j, gcps_df.search_i, gcps_df.dj, gcps_df.di, color='r')
    ax.plot(gcps_df.search_j[nomatch], gcps_df.search_i[nomatch], 'b+')
    this_out = os.path.splitext(fn_img)[0]
    fig.savefig(os.path.join('match_imgs', this_out + '_matches.png'), bbox_inches='tight', dpi=200)

    E = builder.ElementMaker()
    ImMes = E.MesureAppuiFlottant1Im(E.NameIm(fn_img))

    pt_els = get_im_meas(gcps_df, E)
    for p in pt_els:
        ImMes.append(p)
    mkdir_p('Ori-InterneScan')

    outxml = E.SetOfMesureAppuisFlottants(ImMes)
    tree = etree.ElementTree(outxml)
    tree.write(os.path.join('Ori-InterneScan', 'MeasuresIm-' + fn_img + '.xml'), pretty_print=True,
               xml_declaration=True, encoding="utf-8")

    if return_val:
        return gcps_df


def _fix_cross(subimg):
    subimg = subimg.astype(np.float32)

    cross = cross_template(subimg.shape[0], width=5)
    cross[:, :16] = 0
    cross[:, -16:] = 0
    cross[:16, :] = 0
    cross[-16:, :] = 0
    if subimg.shape[0] != cross.shape[0]:
        cross = cross[:subimg.shape[0], :]

    if subimg.shape[1] != cross.shape[1]:
        cross = cross[:, :subimg.shape[1]]

    subimg[cross != 0] = np.nan
    fixed = nanmedian_filter(subimg, footprint=disk(7))
    subimg[np.isnan(subimg)] = fixed[np.isnan(subimg)]
    return subimg.astype(np.uint8)


def remove_crosses(fn_img):
    """

    :param fn_img:
    :return:
    """
    fn_meas = os.path.join('Ori-InterneScan', 'MeasuresIm-{}.xml'.format(fn_img))
    img = io.imread(fn_img)
    gcps = parse_im_meas(fn_meas)

    pt = np.round([gcps.i[0], gcps.j[0]]).astype(int)
    subim, row_, col_ = make_template(img, pt, 200)
    cross = cross_template(subim.shape[0], width=5)
    cross[:, 16:22] = 1
    cross[:, -24:-16] = 1
    cross[16:22, :] = 1
    cross[-24:-16, :] = 1

    subim = subim.astype(np.float32)
    subim[cross != 0] = np.nan
    fixed = nanmedian_filter(subim, footprint=disk(7))
    subim[np.isnan(subim)] = fixed[np.isnan(subim)]
    img[int(pt[0]) - row_[0]:int(pt[0]) + row_[1] + 1, int(pt[1]) - col_[0]:int(pt[1]) + col_[1] + 1] = subim.astype(
        np.uint8)

    for i, row in gcps.loc[1:].iterrows():
        pt = np.round([row.i, row.j]).astype(int)
        subim, row_, col_ = make_template(img, pt, 200)
        img[pt[0] - row_[0]:pt[0] + row_[1] + 1, pt[1] - col_[0]:pt[1] + col_[1] + 1] = _fix_cross(subim)

    mkdir_p('original')
    shutil.move(fn_img, 'original')
    io.imsave(fn_img, img.astype(np.uint8))


def join_hexagon(im_pattern, overlap=2000, blend=True, corona=False):
    """

    :param im_pattern:
    :param overlap:
    :param blend:
    :param corona:
    :return:
    """
    if not corona:
        left = io.imread('{}_a.tif'.format(im_pattern))
        right = io.imread('{}_b.tif'.format(im_pattern))

        left_gd = gdal.Open('{}_a.tif'.format(im_pattern))
        right_gd = gdal.Open('{}_b.tif'.format(im_pattern))

        M = match_halves(left, right, overlap=overlap)

        out_shape = (left.shape[0], left.shape[1] + right.shape[1])

        combined_right = warp(right, M, output_shape=out_shape, preserve_range=True, order=3)

        combined_left = np.zeros(out_shape)
        combined_left[:, :left.shape[1]] = left

        if blend:
            combined = _blend(combined_left, combined_right, left.shape)
        else:
            combined = combined_left + combined_right

        last_ind = np.where(np.sum(combined, axis=0) > 0)[0][-1]
        combined = combined[:, :last_ind+1]

        io.imsave('{}.tif'.format(im_pattern), combined.astype(np.uint8))
    else:
        left = io.imread('{}_d.tif'.format(im_pattern))
        for right_img in ['c', 'b', 'a']:
            right = io.imread('{}_{}.tif'.format(im_pattern, right_img))

            M = match_halves(left, right, overlap=overlap)
            out_shape = (left.shape[0], left.shape[1] + right.shape[1])

            combined_right = warp(right, M, output_shape=out_shape, preserve_range=True, order=3)

            combined_left = np.zeros(out_shape)
            combined_left[:, :left.shape[1]] = left

            if blend:
                combined = _blend(combined_left, combined_right, left.shape)
            else:
                combined = combined_left + combined_right

            last_ind = np.where(np.sum(combined, axis=0) > 0)[0][-1]
            left = combined[:, :last_ind + 1]

        io.imsave('{}.tif'.format(im_pattern), combined.astype(np.uint8))


def _blend(_left, _right, left_shape):
    first = np.where(np.sum(_right, axis=0) > 0)[0][0]
    last = left_shape[1]

    m = 1 / (first - last)
    alpha = np.ones(_left.shape, dtype=np.float32)
    alpha[:, last:] = 0
    for i, ind in enumerate(np.arange(first, last)):
        alpha[:, ind] = 1 + m * (ind - first)

    return alpha * _left + (1 - alpha) * _right


def match_halves(left, right, overlap):
    """

    :param left:
    :param right:
    :param overlap:
    :return:
    """
    src_pts = []
    dst_pts = []

    row_inds = list(range(0, left.shape[0] + 1, overlap))
    if row_inds[-1] != left.shape[0]:
        row_inds.append(-1)
    row_inds = np.array(row_inds)

    for i, ind in enumerate(row_inds[:-1]):
        try:
            l_sub = left[ind:row_inds[i + 1], -overlap:]
            r_sub = right[ind:row_inds[i + 1], :overlap]

            kp, des, matches = get_matches(l_sub, r_sub)
            src_pts.extend(
                [np.array(kp[0][m.queryIdx].pt) + np.array([left.shape[1] - overlap, ind]) for m in matches])
            dst_pts.extend([np.array(kp[1][m.trainIdx].pt) + np.array([0, ind]) for m in matches])
        except:
            continue

    model, inliers = ransac((np.array(src_pts), np.array(dst_pts)), EuclideanTransform,
                        min_samples=25, residual_threshold=1, max_trials=1000)
    print('{} tie points found'.format(np.count_nonzero(inliers)))
    return model
