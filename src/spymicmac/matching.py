"""
spymicmac.matching is a collection of tools for matching templates in images
"""
import os
import shutil
import multiprocessing as mp
import cv2
import matplotlib.pyplot as plt
import pandas as pd
import lxml.etree as etree
import lxml.builder as builder
from skimage import morphology, io, filters
from skimage.morphology import binary_dilation, disk
from skimage.measure import ransac
from skimage.feature import peak_local_max
from skimage.transform import AffineTransform, EuclideanTransform
from scipy.interpolate import RectBivariateSpline as RBS
import numpy as np
from shapely.ops import nearest_points
from shapely.geometry import LineString, MultiPoint, Point
import geopandas as gpd
from pybob.image_tools import nanmedian_filter
from pybob.ddem_tools import nmad
from spymicmac import image, micmac, resample


######################################################################################################################
# tools for matching fiducial markers (or things like fiducial markers)
######################################################################################################################
def cross_template(shape, width=3):
    """
    Create a cross-shaped template for matching reseau or fiducial marks.

    :param int shape: the output shape of the template
    :param int width: the width of the cross at the center of the template (default: 3 pixels).
    :return: **cross** (*array-like*) -- the cross template
    """
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


def find_crosses(img, cross):
    """
    Find all cross markers in an image.

    :param array-like img: the image
    :param array-like cross: the cross template to use
    :return: **grid_df** (*pandas.DataFrame*) -- a dataframe of marker locations and offsets
    """

    sub_coords = []
    simgs, top_inds, left_inds = image.splitter(img, (4, 8), overlap=4*cross.shape[0])

    for ind, simg in enumerate(simgs):
        img_inv = img.max() - simg
        res = cv2.matchTemplate(img_inv.astype(np.uint8), cross.astype(np.uint8), cv2.TM_CCORR_NORMED)

        these_coords = peak_local_max(res, min_distance=int(1.5*cross.shape[0]),
                                      threshold_abs=np.median(res)).astype(np.float64)

        these_coords += cross.shape[0] / 2 - 0.5

        these_coords[:, 0] += top_inds[ind]
        these_coords[:, 1] += left_inds[ind]

        sub_coords.append(these_coords)

    coords = np.concatenate(sub_coords, axis=0)

    grid_df = match_reseau_grid(img, coords, cross)
    # if we have missing values, we find matches individually
    if grid_df.dist.isnull().sum() > 0:
        tsize = int(cross.shape[0] * 1.5)
        for ind, row in grid_df[grid_df.isnull().any(axis=1)].iterrows():
            subimg, _, _ = make_template(img, row[['search_i', 'search_j']].astype(int), tsize)
            res, this_i, this_j = find_match(subimg.astype(np.uint8), cross.astype(np.uint8), how='min')

            grid_df.loc[ind, 'match_i'] = this_i - tsize + row['search_i']
            grid_df.loc[ind, 'match_j'] = this_j - tsize + row['search_j']
            grid_df.loc[ind, 'dist'] = np.sqrt((grid_df.loc[ind, 'match_i'] - grid_df.loc[ind, 'grid_i'])**2 +
                                               (grid_df.loc[ind, 'match_j'] - grid_df.loc[ind, 'grid_j'])**2)

    return grid_df


def match_reseau_grid(img, coords, cross):
    """
    Find the best match for each KH-9 mapping camera reseau grid point, given a list of potential matches.

    :param array-like img: the image to use
    :param array-like coords: the coordinates of the potential matches
    :param array-like cross: the cross template to use.
    :return: **grid_df** (*pandas.DataFrame*) -- a DataFrame of grid locations and match points
    """
    matchpoints = MultiPoint(coords[:, ::-1])
    # top, bot, left, right = find_grid_border(coords)
    left, right, top, bot = image.get_rough_frame(img)

    lr_valid = left > 0 and right < 1e10
    tb_valid = top > 0 and bot < 1e10

    if lr_valid and tb_valid:
        scale = np.mean([(bot - top - 3 * cross.shape[0]) / 22, (right - left - cross.shape[0]) / 46])

    elif lr_valid and not tb_valid:
        # if the left/right border is okay, use it to guess the top/bottom border
        scale = (right - left - cross.shape[0]) / 46

        if top > 0:
            bot = int(top + scale * 22 + 1.5 * cross.shape[0])
        elif bot < 1e10:
            top = int(bot - scale * 22 - 1.5 * cross.shape[0])
        else: # if both or wrong, use very approximate average values
            top = 1200
            bot = 34000

    elif not lr_valid and tb_valid:
        # if the top/bottom border is okay, use it to guess the left/right border
        scale = (bot - top - 3 * cross.shape[0]) / 22
        if left > 0:
            right = int(left + scale * 46 + 0.5 * cross.shape[0])
        elif right < 1e10:
            left = int(right - scale * 46 - 0.5 * cross.shape[0])
        else:
            left = 2000
            right = 68500
    else:
        # if we can't find the image border, try using the mark coordinates.
        top, bot, left, right = _grid_border(coords)

    left_edge = LineString([(left + 0.5 * cross.shape[0], bot - 1.5 * cross.shape[0]),
                            (left + 0.5 * cross.shape[0], top + 1.5 * cross.shape[0])])
    right_edge = LineString([(right - 0.5 * cross.shape[0], bot - 1.5 * cross.shape[0]),
                             (right - 0.5 * cross.shape[0], top + 1.5 * cross.shape[0])])

    II, JJ, search_grid = _search_grid(left_edge, right_edge, True)

    grid_df = pd.DataFrame()
    for ii, pr in enumerate(list(zip(II, JJ))):
        grid_df.loc[ii, 'gcp'] = 'GCP_{}_{}'.format(pr[0], pr[1])
        grid_df.loc[ii, 'grid_j'] = left + scale * pr[1] + 0.5 * cross.shape[0]
        grid_df.loc[ii, 'grid_i'] = bot - scale * pr[0] - 1.5 * cross.shape[0]

    grid_df['search_j'] = np.array(search_grid)[:, 1]
    grid_df['search_i'] = np.array(search_grid)[:, 0]

    for ind, row in grid_df.iterrows():
        gridpt = Point(row.grid_j, row.grid_i)
        _, matchpt = nearest_points(gridpt, matchpoints)
        grid_df.loc[ind, 'match_i'] = matchpt.y
        grid_df.loc[ind, 'match_j'] = matchpt.x
        grid_df.loc[ind, 'dist'] = gridpt.distance(matchpt)

    model, inliers = ransac((grid_df[['grid_j', 'grid_i']].values, grid_df[['match_j', 'match_i']].values),
                            AffineTransform, min_samples=10, residual_threshold=grid_df.dist.median(), max_trials=5000)
    grid_df.loc[~inliers, 'dist'] = np.nan
    grid_df.loc[~inliers, 'match_j'] = np.nan
    grid_df.loc[~inliers, 'match_i'] = np.nan

    return grid_df


def _grid_border(coords):
    top = np.percentile(coords[:, 0], 1)
    bot = np.percentile(coords[:, 0], 99)
    left = np.percentile(coords[:, 1], 1)
    right = np.percentile(coords[:, 1], 99)

    return top, bot, left, right


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


def remove_crosses(fn_img, nproc=1):
    """
    Remove the Reseau marks from a KH-9 image before re-sampling.

    :param str fn_img: the image filename.
    :param int nproc: the number of subprocesses to use (default: 1).
    """
    fn_meas = os.path.join('Ori-InterneScan', 'MeasuresIm-{}.xml'.format(fn_img))
    img = io.imread(fn_img)
    gcps = micmac.parse_im_meas(fn_meas)

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

    if nproc == 1:
        for i, row in gcps.loc[1:].iterrows():
            pt = np.round([row.i, row.j]).astype(int)
            subim, row_, col_ = make_template(img, pt, 200)
            img[pt[0] - row_[0]:pt[0] + row_[1] + 1, pt[1] - col_[0]:pt[1] + col_[1] + 1] = _fix_cross(subim)
    else:
        pool = mp.Pool(nproc, maxtasksperchild=1)

        subimgs = []
        rows = []
        cols = []
        points = []

        for ind, row in gcps.loc[1:].iterrows():
            pt = np.round([row.i, row.j]).astype(int)
            subim, row_, col_ = make_template(img, pt, 200)
            subimgs.append(subim)
            rows.append(row_)
            cols.append(col_)
            points.append(pt)

        outputs = pool.map(_fix_cross, subimgs, chunksize=1)
        pool.close()
        pool.join()

        for subim, row, col, pt in zip(outputs, rows, cols, points):
            img[pt[0] - row_[0]:pt[0] + row_[1] + 1, pt[1] - col_[0]:pt[1] + col_[1] + 1] = subim

    os.makedirs('original', exist_ok=True)
    shutil.move(fn_img, 'original')
    io.imsave(fn_img, img.astype(np.uint8))


def _search_grid(left, right, joined=False):
    i_grid = []
    j_grid = []
    search_pts = []

    if joined:
        nj = 47
    else:
        nj = 24
    for ii in range(0, 23):
        this_left = left.interpolate(ii * left.length / 22)
        this_right = right.interpolate(ii * right.length / 22)
        this_line = LineString([this_left, this_right])
        for jj in range(0, nj):
            this_pt = this_line.interpolate(jj * this_line.length / (nj - 1))
            i_grid.append(ii)
            j_grid.append(jj)
            search_pts.append((this_pt.y, this_pt.x))

    return i_grid, j_grid, search_pts


def _outlier_filter(vals, n=3):
    return np.abs(vals - np.nanmean(vals)) > n * np.nanstd(vals)


def find_reseau_grid(fn_img, csize=361, return_val=False):
    """
    Find the locations of the Reseau marks in a scanned KH-9 image. Locations are saved
    to Ori-InterneScan/MeasuresIm-:fn_img:.xml.

    :param str fn_img: the image filename.
    :param int csize: the size of the cross template (default: 361 -> 361x361)
    :param bool return_val: return a pandas DataFrame of the Reseau mark locations (default: False).
    :return: **gcps_df** (*pandas.DataFrame*) -- a DataFrame of the Reseau mark locations (if return_val=True).
    """
    print('Reading {}'.format(fn_img))
    img = io.imread(fn_img)
    img_lowres = resample.downsample(img, fact=10)

    print('Image read.')
    cross = cross_template(csize, width=3)
    cross[cross > 1] = 0
    cross *= 255

    fig = plt.figure(figsize=(7, 12))
    ax = fig.add_subplot(111)
    ax.imshow(img_lowres, cmap='gray', extent=[0, img.shape[1], img.shape[0], 0])
    ax.set_xticks([])
    ax.set_yticks([])

    print('Finding grid points in {}...'.format(fn_img))
    grid_df = find_crosses(img, cross)

    model, inliers = ransac((grid_df[['grid_j', 'grid_i']].values, grid_df[['match_j', 'match_i']].values),
                            AffineTransform, min_samples=10, residual_threshold=10, max_trials=5000)
    grid_df['resid'] = model.residuals(grid_df.dropna()[['grid_j', 'grid_i']].values,
                                       grid_df.dropna()[['match_j', 'match_i']].values)

    grid_df.loc[~inliers, ['match_j', 'match_i']] = np.nan

    dst = model(grid_df[['grid_j', 'grid_i']].values)

    x_res = grid_df['match_j'] - dst[:, 0]
    y_res = grid_df['match_i'] - dst[:, 1]

    ux = x_res.values.reshape(23, 47)
    uy = y_res.values.reshape(23, 47)

    xdiff = ux - nanmedian_filter(ux, footprint=disk(3))
    ydiff = uy - nanmedian_filter(uy, footprint=disk(3))

    xout = _outlier_filter(xdiff, n=5)
    yout = _outlier_filter(ydiff, n=5)

    outliers = np.logical_or.reduce([~inliers, xout.flatten(), yout.flatten()])

    ux[outliers.reshape(23, 47)] = np.nan
    uy[outliers.reshape(23, 47)] = np.nan

    ux[outliers.reshape(23, 47)] = nanmedian_filter(ux, footprint=disk(1))[outliers.reshape(23, 47)]
    uy[outliers.reshape(23, 47)] = nanmedian_filter(uy, footprint=disk(1))[outliers.reshape(23, 47)]

    grid_df.loc[outliers, 'match_j'] = dst[outliers, 0] + ux.flatten()[outliers]
    grid_df.loc[outliers, 'match_i'] = dst[outliers, 1] + uy.flatten()[outliers]

    grid_df.loc[np.isnan(grid_df['match_j']), 'match_j'] = dst[np.isnan(grid_df['match_j']), 0]
    grid_df.loc[np.isnan(grid_df['match_i']), 'match_i'] = dst[np.isnan(grid_df['match_i']), 1]

    grid_df['im_row'] = grid_df['match_i']
    grid_df['im_col'] = grid_df['match_j']

    grid_df['dj'] = grid_df['match_j'] - dst[:, 0]
    grid_df['di'] = grid_df['match_i'] - dst[:, 1]

    print('Grid points found.')
    os.makedirs('match_imgs', exist_ok=True)

    print('Mean x residual: {:.2f} pixels'.format(grid_df.loc[~outliers, 'dj'].abs().mean()))
    print('Mean y residual: {:.2f} pixels'.format(grid_df.loc[~outliers, 'di'].abs().mean()))
    print('Mean residual: {:.2f} pixels'.format(grid_df.loc[~outliers, 'resid'].mean()))

    ax.quiver(grid_df.match_j, grid_df.match_i, grid_df.dj, grid_df.di, color='r')
    ax.plot(grid_df.match_j[outliers], grid_df.match_i[outliers], 'b+')

    this_out = os.path.splitext(fn_img)[0]
    fig.savefig(os.path.join('match_imgs', this_out + '_matches.png'), bbox_inches='tight', dpi=200)

    E = builder.ElementMaker()
    ImMes = E.MesureAppuiFlottant1Im(E.NameIm(fn_img))

    pt_els = micmac.get_im_meas(grid_df, E)
    for p in pt_els:
        ImMes.append(p)
    os.makedirs('Ori-InterneScan', exist_ok=True)

    outxml = E.SetOfMesureAppuisFlottants(ImMes)
    tree = etree.ElementTree(outxml)
    tree.write(os.path.join('Ori-InterneScan', 'MeasuresIm-' + fn_img + '.xml'), pretty_print=True,
               xml_declaration=True, encoding="utf-8")

    if return_val:
        return grid_df


def wagon_wheel(size, width=3, mult=255):
    """
    Creates a template in the shape of a "wagon wheel" (a cross inscribed in a ring).

    :param int size: the width (and height) of the template, in pixels
    :param int width: the width/thickness of the cross, in pixels
    :param mult: a multiplier to use for the template [default: 255]

    :return: **template** (*array-like*) the wagon wheel template
    """
    cross = cross_template(size, width)
    cross[cross > 1] = 0

    templ = disk(int(size / 2))
    padded = np.zeros(templ.shape, dtype=templ.dtype)
    padded[width:-width, width:-width] = disk(int((size - 2 * width) / 2))

    templ -= padded
    templ += cross.astype(templ.dtype)
    templ[templ > 1] = 1

    return mult * templ


# all credit to joe kennedy for the name of this function.
def ocm_show_wagon_wheels(img, size, width=3, img_border=None):
    """
    Find all "wagon wheel" markers in an image.

    :param array-like img: the image
    :param int size: the size of the marker (in pixels)
    :param int width: the width/thickness of the cross, in pixels (default: 3)
    :param img_border: the approximate top and bottom rows of the image frame. If not set,
        calls get_rough_frame() on the image.
    :return: **coords** an Nx2 array of the location of the detected markers.
    """
    if img_border is None:
        _, _, top, bot = image.get_rough_frame(img)
        if top < 0 or bot > 1e10:
            raise RuntimeError("Unable to find image border. Try running again with approximate values.")
    else:
        top, bot = img_border

    templ = wagon_wheel(size, width=width)

    img_top = np.zeros(img[:top, :].shape, dtype=np.uint8)
    img_bot = np.zeros(img[bot:, :].shape, dtype=np.uint8)

    img_top[img[:top, :] > filters.threshold_local(img[:top, :], 101, method='mean')] = 1
    img_bot[img[bot:, :] > filters.threshold_local(img[bot:, :], 101, method='mean')] = 1

    res_top = cv2.matchTemplate(img_top.astype(np.uint8), templ.astype(np.uint8), cv2.TM_CCORR_NORMED)
    res_bot = cv2.matchTemplate(img_bot.astype(np.uint8), templ.astype(np.uint8), cv2.TM_CCORR_NORMED)

    coords_top = peak_local_max(res_top, threshold_abs=0.7).astype(np.float64)
    coords_bot = peak_local_max(res_bot, threshold_abs=0.7).astype(np.float64)

    coords_bot[:, 0] += bot

    coords_top += size / 2 - 0.5
    coords_bot += size / 2 - 0.5

    return np.concatenate((coords_top, coords_bot), axis=0)


def find_rail_marks(img):
    """
    Find all rail marks along the bottom edge of a KH-4 style image.

    :param array-like img: the image to find the rail marks in.
    :return: **coords** (*array-like*) -- an Nx2 array of the location of the detected markers.
    """
    left, right, top, bot = image.get_rough_frame(img)
    img_lowres = resample.downsample(img, fact=10)

    templ = np.zeros((21, 21), dtype=np.uint8)
    templ[5:-5, 5:-5] = 255 * disk(5)  # rail marks are approximately 100 x 100 pixels, so lowres is 10 x 10
    res = cv2.matchTemplate(img_lowres.astype(np.uint8), templ.astype(np.uint8), cv2.TM_CCORR_NORMED)

    norm_res = res - res.mean()

    coords = peak_local_max(norm_res, threshold_abs=2.5*norm_res.std(), min_distance=10).astype(np.float64)
    coords += templ.shape[0] / 2 - 0.5
    coords *= 10

    bottom_rail = np.logical_and.reduce([coords[:, 1] > left,
                                         coords[:, 1] < right,
                                         np.abs(coords[:, 0] - bot) < 400])
    out_coords = coords[bottom_rail]

    valid = _refine_rail(coords[bottom_rail])

    return out_coords[valid]


def _refine_rail(coords):
    valid = np.isfinite(coords[:, 1])
    prev_valid = np.count_nonzero(valid)

    nout = 1

    while nout > 0:
        p = np.polyfit(coords[valid, 1], coords[valid, 0], 1)
        fit = np.polyval(p, coords[:, 1])

        diff = coords[:, 0] - fit
        valid = np.abs(diff - np.median(diff)) < 4 * nmad(diff)

        nout = prev_valid - np.count_nonzero(valid)
        prev_valid = np.count_nonzero(valid)

    return valid


def notch_template(size):
    """
    Create a notch-shaped ("^") template.

    :param int size: the size of the template, in pixels
    :return: **template** (*array-like*) -- the notch template
    """
    template = np.zeros((size, size), dtype=np.uint8)
    template[-1, :] = 1
    for ind in range(1, int(size/2) + 1):
        template[-ind-1, ind:-ind] = 1
    return 255 * template


def find_kh4_notches(img, size=101):
    """
    Find all 4 notches along the top of a KH-4 style image.

    :param array-like img: the image.
    :param int size: the size of the notch template to use.
    :return: **coords** (*array-like*) -- a 4x2 array of notch locations
    """
    left, right, top, bot = image.get_rough_frame(img)

    templ = notch_template(size)

    pcts = np.array([0.03, 0.5, 0.6, 0.97])

    search_j = (left + pcts * (right - left)).astype(int)
    search_i = top * np.ones(4).astype(int)

    matches = []

    for jj, ii in zip(search_j, search_i):
        subimg, _, _ = make_template(img, [ii, jj], 400)
        # res, this_i, this_j = find_match(subimg, templ, how='max', eq=False)
        res = cv2.matchTemplate(subimg.astype(np.uint8), templ.astype(np.uint8), cv2.TM_CCORR_NORMED)

        try:
            this_i, this_j = peak_local_max(res, min_distance=50, num_peaks=1)[0]
        except IndexError:
            this_i = this_j = np.nan

        matches.append((this_i - 400 + ii, this_j - 400 + jj))

    return np.array(matches)


######################################################################################################################
# tools for matching gcps
######################################################################################################################
def make_template(img, pt, half_size):
    """
    Return a sub-section of an image to use for matching.

    :param array-like img: the image from which to create the template
    :param tuple pt: the (row, column) center of the template
    :param int half_size: the half-size of the template; template size will be 2 * half_size + 1
    :return:
        - **template** (*array-like*) -- the template
        - **row_inds** (*list*) -- the number of rows above/below the center of the template
        - **col_inds** (*list*) -- the number of columns left/right of the center of the template
    """
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


def find_match(img, template, how='min', eq=True):
    """
    Given an image and a template, find a match using openCV's normed cross-correlation.

    :param array-like img: the image to find a match in
    :param array-like template: the template to use for matching
    :param str how: determines whether the match is the minimum or maximum correlation (default: min)
    :param bool eq: use a rank equalization filter before matching the templates (default: True)
    :return:
        - **res** (*array-like*) -- the correlation image
        - **match_i** (*float*) -- the row location of the match
        - **match_j** (*float*) -- the column location of the match
    """
    assert how in ['min', 'max'], "have to choose min or max"

    if eq:
        img_eq = filters.rank.equalize(img, footprint=morphology.disk(100))
        res = cv2.matchTemplate(img_eq, template, cv2.TM_CCORR_NORMED)
    else:
        res = cv2.matchTemplate(img, template, cv2.TM_CCORR_NORMED)

    i_off = (img.shape[0] - res.shape[0]) / 2
    j_off = (img.shape[1] - res.shape[1]) / 2
    if how == 'min':
        val, _, loc, _ = cv2.minMaxLoc(res)
    elif how == 'max':
        _, val, _, loc = cv2.minMaxLoc(res)

    this_j, this_i = loc

    sp_delx, sp_dely = _subpixel(res, how=how)
    # sp_delx, sp_dely = 0, 0
    return res, this_i + i_off + sp_dely, this_j + j_off + sp_delx


def _subpixel(res, how='min'):
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


def find_gcp_match(img, template, method=cv2.TM_CCORR_NORMED):
    res = cv2.matchTemplate(img, template, method)
    i_off = (img.shape[0] - res.shape[0]) / 2
    j_off = (img.shape[1] - res.shape[1]) / 2
    _, maxval, _, maxloc = cv2.minMaxLoc(res)
    maxj, maxi = maxloc
    try:
        sp_delx, sp_dely = _subpixel(res, how='max')
    except ValueError as e:
        sp_delx = 0
        sp_dely = 0

    return res, maxi + i_off + sp_dely, maxj + j_off + sp_delx


def find_grid_matches(tfm_img, refgeo, mask, initM=None, spacing=200, srcwin=60, dstwin=600):
    """
    Find matches between two images on a grid using normalized cross-correlation template matching.

    :param array-like tfm_img: the image to use for matching.
    :param GeoImg refgeo: the reference image to use for matching.
    :param array-like mask: a mask indicating areas that should be used for matching.
    :param initM: the model used for transforming the initial, non-georeferenced image.
    :param int spacing: the grid spacing, in pixels (default: 200 pixels)
    :param int srcwin: the half-size of the template window.
    :param int dstwin: the half-size of the search window.
    :return: **gcps** (*pandas.DataFrame*) -- a DataFrame with GCP locations, match strength, and other information.
    """
    match_pts = []
    z_corrs = []
    peak_corrs = []

    jj = np.arange(srcwin, spacing * np.ceil((refgeo.img.shape[1]-srcwin) / spacing) + 1, spacing).astype(int)
    ii = np.arange(srcwin, spacing * np.ceil((refgeo.img.shape[0]-srcwin) / spacing) + 1, spacing).astype(int)

    search_pts = []

    for _i in ii:
        for _j in jj:
            search_pts.append((_j, _i))
            match, z_corr, peak_corr = do_match(tfm_img, refgeo.img, mask, (_i, _j), srcwin, dstwin)
            match_pts.append(match)
            z_corrs.append(z_corr)
            peak_corrs.append(peak_corr)

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


def do_match(dest_img, ref_img, mask, pt, srcwin, dstwin):
    """
    Find a match between two images using normalized cross-correlation template matching.

    :param array-like dest_img: the image to search for the matching point in.
    :param array-like ref_img: the reference image to use for matching.
    :param array-like mask: a mask indicating areas that should be used for matching.
    :param array-like pt: the index (i, j) to search for a match for.
    :param int srcwin: the half-size of the template window.
    :param int dstwin: the half-size of the search window.
    :return:
        - **match_pt** (*tuple*) -- the matching point (j, i) found in dest_img
        - **z_corr** (*float*) -- number of standard deviations (z-score) above other potential matches
        - **peak_corr** (*float*) -- the correlation value of the matched point
    """
    _i, _j = pt
    submask, _, _ = make_template(mask, pt, srcwin)

    if np.count_nonzero(submask) / submask.size < 0.05:
        return (np.nan, np.nan), np.nan, np.nan

    try:
        testchip, _, _ = make_template(ref_img, pt, srcwin)
        dst_chip, _, _ = make_template(dest_img, pt, dstwin)

        testchip[np.isnan(testchip)] = 0
        dst_chip[np.isnan(dst_chip)] = 0

        test = image.highpass_filter(testchip)
        dest = image.highpass_filter(dst_chip)

        testmask = binary_dilation(testchip == 0, footprint=disk(8))
        destmask = binary_dilation(dst_chip == 0, footprint=disk(8))

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

    except Exception as e:
        return (np.nan, np.nan), np.nan, np.nan

    return (out_j, out_i), z_corr, peak_corr


def get_matches(img1, img2, mask1=None, mask2=None, dense=False):
    """
    Return keypoint matches found using openCV's ORB implementation.

    :param array-like img1: the first image to match
    :param array-like img2: the second image to match
    :param array-like mask1: a mask to use for the first image. (default: no mask)
    :param array-like mask2: a mask to use for the second image. (default: no mask)
    :param bool dense: compute matches over sub-blocks (True) or the entire image (False). (default: False)
    :return:
        - **keypoints** (*tuple*) -- the keypoint locations for the first and second image.
        - **descriptors** (*tuple*) -- the descriptors for the first and second image.
        - **matches** (*list*) -- a list of matching keypoints between the first and second image
    """
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


def keypoint_grid(img, spacing=25, size=10):
    _j = np.arange(0, img.shape[1], spacing)
    _i = np.arange(0, img.shape[0], spacing)
    _gridj, _gridi = np.meshgrid(_j, _i)
    _grid = np.concatenate((_gridj.reshape(-1, 1), _gridi.reshape(-1, 1)), axis=1)

    return [cv2.KeyPoint(pt[0], pt[1], size) for pt in _grid]


def get_dense_keypoints(img, mask, npix=100, nblocks=None, return_des=False):
    """
    Find ORB keypoints by dividing an image into smaller parts.

    :param array-like img: the image to use.
    :param array-like mask: a mask to use for keypoint generation.
    :param int npix: the block size (in pixels) to divide the image into.
    :param int nblocks: the number of blocks to divide the image into. If set, overrides value given by npix.
    :param bool return_des: return the keypoint descriptors, as well
    :return:
        - **keypoints** (*list*) -- a list of keypoint locations
        - **descriptors** (*list*) -- if requested, a list of keypoint descriptors.
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

    split_img, oy, ox = image.splitter(img, (y_tiles, x_tiles), overlap=olap)
    split_msk, _, _ = image.splitter(mask, (y_tiles, x_tiles), overlap=olap)

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


def match_halves(left, right, overlap, block_size=None):
    """
    Find a transformation to join the left and right halves of an image scan.

    :param array-like left: the left-hand image scan.
    :param array-like right: the right-hand image scan.
    :param int overlap: the estimated overlap between the two images, in pixels.
    :param int block_size: the number of rows each sub-block should cover. Defaults to overlap.
    :return: **model** (*EuclideanTransform*) -- the estimated Euclidean transformation between the two image halves.
    """
    src_pts = []
    dst_pts = []

    if block_size is None:
        block_size = overlap

    row_inds = list(range(0, left.shape[0] + 1, block_size))
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
                        min_samples=10, residual_threshold=2, max_trials=25000)

    print('{} tie points found'.format(np.count_nonzero(inliers)))
    return model, np.count_nonzero(inliers)
