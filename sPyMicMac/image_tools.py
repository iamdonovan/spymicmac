"""
sPyMicMac.image_tools is a collection of tools for working with KH-9 Hexagon imagery.
"""
import os
from glob import glob
import cv2
from itertools import chain
import gdal
from rtree import index
from skimage import morphology
from skimage.filters import rank
from skimage import exposure, transform
from skimage.measure import ransac
from skimage.feature import peak_local_max
from skimage.transform import match_histograms, warp, AffineTransform, EuclideanTransform
from scipy.interpolate import RectBivariateSpline as RBS
from scipy import ndimage
import numpy as np
from shapely.ops import cascaded_union
import geopandas as gpd
# import pyvips
from llc import jit_filter_function
from pybob.image_tools import match_hist, reshape_geoimg, create_mask_from_shapefile
from pymmaster.mmaster_tools import orient_footprint


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


def cross_filter(img, cross):
    cross_edge = cross == 2
    cross_cent = cross == 1
    edge_std = ndimage.filters.generic_filter(highpass_filter(img), nanstd, footprint=cross_edge)
    cent_std = ndimage.filters.generic_filter(highpass_filter(img), nanstd, footprint=cross_cent)
    return np.where(np.logical_and(edge_std != 0, cent_std != 0), cent_std / edge_std, 2)


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
    v = img.copy()
    v[np.isnan(img)] = 0
    vv = ndimage.gaussian_filter(v, 3)

    w = 0 * img.copy() + 1
    w[np.isnan(img)] = 0
    ww = ndimage.gaussian_filter(w, 3)

    tmplow = vv / ww
    tmphi = img - tmplow
    return tmphi


# def splitter(img, nblocks, overlap=0):
#     split1 = np.array_split(img, nblocks[0], axis=0)
#     split2 = [np.array_split(im, nblocks[1], axis=1) for im in split1]
#     olist = [np.copy(a) for a in list(chain.from_iterable(split2))]
#     return olist
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


def stretch_image(img, scale=(0,1), mult=255, outtype=np.uint8):
    maxval = np.nanquantile(img, max(scale))
    minval = np.nanquantile(img, min(scale))

    img[img > maxval] = maxval
    img[img < minval] = minval

    return (mult * (img - minval) / (maxval - minval)).astype(outtype)


def make_binary_mask(img, erode=0, mask_value=0):
    _mask = 255 * np.ones(img.shape, dtype=np.uint8)
    if np.isfinite(mask_value):
        _mask[img == mask_value] = 0
    else:
        _mask[np.isnan(img)] = 0

    if erode > 0:
        erode_mask = morphology.binary_erosion(_mask, selem=morphology.disk(erode))
        _mask[~erode_mask] = 0

    return _mask


######################################################################################################################
# GCP matching tools
######################################################################################################################
def keypoint_grid(img, spacing=25, size=10):
    _j = np.arange(0, img.shape[1], spacing)
    _i = np.arange(0, img.shape[0], spacing)
    _gridj, _gridi = np.meshgrid(_j, _i)
    _grid = np.concatenate((_gridj.reshape(-1, 1), _gridi.reshape(-1, 1)), axis=1)

    return [cv2.KeyPoint(pt[0], pt[1], size) for pt in _grid]


def get_dense_keypoints(img, mask, npix=200, nblocks=None, return_des=False):
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
    idx = index.Index()

    for pos, row in fprints.iterrows():
        idx.insert(pos, row['geometry'].bounds)

    intersections = []
    for poly in fprints['geometry']:
        merged = cascaded_union([fprints.loc[pos, 'geometry'] for pos in idx.intersection(poly.bounds) if fprints.loc[pos, 'geometry'] != poly])
        intersections.append(poly.intersection(merged))
    intersection = cascaded_union(intersections)

    return intersection.minimum_rotated_rectangle


def get_footprint_mask(shpfile, geoimg, filelist, fprint_out=False):
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
        if np.any(np.array(img1.shape) < 200) or np.any(np.array(img2.shape) < 200):
            kp1, des1 = get_dense_keypoints(img1.astype(np.uint8), mask1, nblocks=2, return_des=True)
            kp2, des2 = get_dense_keypoints(img2.astype(np.uint8), mask2, nblocks=2, return_des=True)
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

def find_grid_matches(tfm_img, refgeo, mask, initM=None, spacing=200, srcwin=40, dstwin=600):
    match_pts = []
    z_corrs = []
    peak_corrs = []
    res_imgs = []

    jj = np.arange(0, tfm_img.shape[1], spacing)
    ii = np.arange(0, tfm_img.shape[0], spacing)

    search_pts = []

    for _i in ii:
        for _j in jj:
            search_pts.append((_j, _i))
            # for pt in search_pts:
            # if mask[pt[1], pt[0]] == 0:
            if mask[_i, _j] == 0:
                match_pts.append((-1, -1))
                z_corrs.append(np.nan)
                peak_corrs.append(np.nan)
                res_imgs.append(np.nan)
                continue
            try:
                testchip, _, _ = make_template(tfm_img, (_i, _j), srcwin)
                dst_chip, _, _ = make_template(refgeo.img, (_i, _j), dstwin)

                dst_chip[np.isnan(dst_chip)] = 0

                test = np.ma.masked_values(highpass_filter(testchip), 0)
                dest = np.ma.masked_values(highpass_filter(dst_chip), 0)

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
    gcps['match_j'] = _dst[:, 0]  # points matched in master image
    gcps['match_i'] = _dst[:, 1]

    if initM is not None:
        # have to find the initial, which means back-transforming search_pts with Minit_full
        _src = np.dot(initM.params, np.hstack([search_pts,
                                               np.ones(search_pts[:, 0].shape).reshape(-1, 1)]).T).T[:, :2]
        gcps['orig_j'] = _src[:, 0]  # this should be the back-transformed search_pts
        gcps['orig_i'] = _src[:, 1]

    gcps['search_j'] = search_pts[:, 0]
    gcps['search_i'] = search_pts[:, 1]
    gcps['dj'] = gcps['search_j'] - gcps['match_j']
    gcps['di'] = gcps['search_i'] - gcps['match_i']

    gcps.dropna(inplace=True)

    return gcps


def transform_from_fprint(img, geo, fprint):
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


######################################################################################################################
# image writing
######################################################################################################################
# def join_halves(img, overlap, indir='.', outdir='.', color_balance=True):
#     """
#     Join scanned halves of KH-9 image into one, given a common overlap point.
#
#     :param img: KH-9 image name (i.e., DZB1215-500454L001001) to join. The function will look for open image halves
#         img_a.tif and img_b.tif, assuming 'a' is the left-hand image and 'b' is the right-hand image.
#     :param overlap: Image coordinates for a common overlap point, in the form [x1, y1, x2, y2]. Best results tend to be
#         overlaps toward the middle of the y range. YMMV.
#     :param indir: Directory containing images to join ['.']
#     :param outdir: Directory to write joined image to ['.']
#     :param color_balance: Attempt to color balance the two image halves before joining [True].
#
#     :type img: str
#     :type overlap: array-like
#     :type indir: str
#     :type outdir: str
#     :type color_balance: bool
#     """
#
#     left = pyvips.Image.new_from_file(os.path.sep.join([indir, '{}_a.tif'.format(img)]), memory=True)
#     right = pyvips.Image.new_from_file(os.path.sep.join([indir, '{}_b.tif'.format(img)]), memory=True)
#     outfile = os.path.sep.join([outdir, '{}.tif'.format(img)])
#
#     if len(overlap) < 4:
#         x1, y1 = overlap
#         if x1 < 0:
#             join = left.merge(right, 'horizontal', x1, y1)
#         else:
#             join = right.merge(left, 'horizontal', x1, y1)
#
#         join.write_to_file(outfile)
#     else:
#         x1, y1, x2, y2 = overlap
#
#         join = left.mosaic(right, 'horizontal', x1, y1, x2, y2, mblend=0)
#         if color_balance:
#             balance = join.globalbalance(int_output=True)
#             balance.write_to_file(outfile)
#         else:
#             join.write_to_file(outfile)
#
#     return
