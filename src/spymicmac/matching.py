"""
spymicmac.matching is a collection of tools for matching templates in images
"""
import os
from pathlib import Path
import shutil
import itertools
import multiprocessing as mp
import cv2
import matplotlib.pyplot as plt
import pandas as pd
from itertools import combinations
from skimage import morphology, io, filters
from skimage.morphology import binary_dilation, disk
from skimage.measure import ransac
from skimage.feature import peak_local_max, ORB
from skimage.transform import ProjectiveTransform, AffineTransform, EuclideanTransform, SimilarityTransform
from scipy.interpolate import RectBivariateSpline as RBS
from scipy import ndimage
import numpy as np
from shapely.ops import nearest_points
from shapely.geometry import LineString, MultiPoint, Point
import geopandas as gpd
import geoutils as gu
from . import image, micmac, resample, register
from numpy.typing import NDArray
from typing import Union


######################################################################################################################
# tools for matching fiducial markers (or things like fiducial markers)
######################################################################################################################
def find_fiducials(fn_img: str, templates: dict, fn_cam: Union[str, Path, None] = None, thresh_tol: float = 0.9,
                   npeaks: int = 5, min_dist: int = 1, angle: Union[float, None] = None,
                   use_frame: bool = True, tsize: Union[int, None] = None, threshold: bool =True,
                   dual: bool = False) -> float:
    """
    Match the location of fiducial markers for a scanned aerial photo.

    :param fn_img: the filename of the image to find fiducial markers in.
    :param templates: a dict of (name, template) pairs corresponding to each fiducial marker.
    :param fn_cam: the filename of the MeasuresIm file for the template image, if templates are created using
        templates_from_meas().
    :param thresh_tol: the minimum relative peak intensity to use for detecting matches
    :param npeaks: maximum number of potential matches to accept for each fiducial marker template
    :param min_dist: the minimum distance allowed between potential peaks
    :param angle: the angle by which to rotate the points in MeasuresCam
    :param use_frame: use the rough image frame to try to find fiducial markers
    :param tsize: target half-size to use for matching (default: calculated based on image size)
    :param threshold: use a local threshold to help find matches
    :param dual: match using both thresholding and not thresholding to help find matches
    :return: the mean residual between the matches and the fiducial marker locations
    """
    img = io.imread(fn_img)

    measures_cam = micmac.parse_im_meas(Path('Ori-InterneScan', 'MeasuresCamera.xml'))
    measures_cam.set_index('name', inplace=True)

    if angle is not None:
        measures_cam = _rotate_meas(measures_cam, angle)

    if fn_cam is not None:
        measures_img = micmac.parse_im_meas(fn_cam)
        measures_img.set_index('name', inplace=True)
        measures_img.rename(columns={'i': 'rough_i', 'j': 'rough_j'}, inplace=True)

        measures_cam = measures_cam.join(measures_img, lsuffix='_cam')
    else:
        # now, get the fractional locations in the image of each marker
        if not use_frame:
            measures_cam = _get_rough_locs(measures_cam)
            measures_cam['rough_j'] *= img.shape[1]
            measures_cam['rough_i'] *= img.shape[0]
        else:
            measures_cam = _get_rough_locs(measures_cam, img)

    # get all potential matches based on our input parameters
    coords_all = _get_all_fid_matches(img, templates, measures_cam, thresh_tol, min_dist,
                                      npeaks, use_frame, tsize, threshold)
    if dual:
        coords_all = pd.concat([coords_all,
                                _get_all_fid_matches(img, templates, measures_cam, thresh_tol, min_dist,
                                                     npeaks, use_frame, tsize, ~threshold)],
                               ignore_index=True)

    # filter based on the best match
    coords_all = _filter_fid_matches(coords_all, measures_cam)

    # now, drop any duplicated values - if we have these, we need to replace/estimate
    # coords_all = coords_all.sort_values('resid').drop_duplicates(subset=['im_col', 'im_row']).sort_values('gcp')
    if len(coords_all) < len(measures_cam):
        print('One or more markers could not be found.')
        # coords_all = _fix_fiducials(coords_all, measures_cam)
        residuals = np.array(len(coords_all) * [np.nan])
    else:
        _, residuals = _get_residuals(coords_all, measures_cam)

    print(f"Mean residual: {residuals.mean():.2f} pixels")

    # write the measures
    micmac.write_measures_im(coords_all, fn_img)

    return residuals.mean()


def _get_all_fid_matches(img, templates, measures_cam, thresh_tol=0.9, min_dist=1, npeaks=5,
                         use_frame=False, tsize=None, threshold=True):

    coords_all = []

    if tsize is None:
        if not use_frame:
            tsize = int(min(0.075 * np.array([measures_cam.rough_j.max() - measures_cam.rough_j.min(),
                                              measures_cam.rough_i.max() - measures_cam.rough_i.min()])))
        else:
            tsize = int(min(0.04 * np.array([measures_cam.rough_j.max() - measures_cam.rough_j.min(),
                                             measures_cam.rough_i.max() - measures_cam.rough_i.min()])))

    bsize = _odd(int(tsize/4))

    for fid, row in measures_cam.iterrows():
        templ = templates[fid]

        subimg, isize, jsize = make_template(img, (row['rough_i'], row['rough_j']), half_size=tsize)
        if threshold:
            thresh = subimg > filters.threshold_local(subimg, block_size=bsize)
            res = cv2.matchTemplate(thresh.astype(np.uint8), templ.astype(np.uint8), cv2.TM_CCORR_NORMED)
        else:
            res = cv2.matchTemplate(subimg.astype(np.uint8), templ.astype(np.uint8), cv2.TM_CCORR_NORMED)

        coords = peak_local_max(res, threshold_rel=thresh_tol, min_distance=min_dist, num_peaks=npeaks).astype(float)

        # get subpixel by looking in a small window around each "peak"
        for ind, coord in enumerate(coords):
            ii, jj = coord.astype(int)
            subres, _, _ = make_template(res, (ii, jj), half_size=3)

            sub_x, sub_y = _subpixel(subres, how='max')
            coords[ind, 1] += sub_x
            coords[ind, 0] += sub_y

        coords += templ.shape[0] / 2 - 0.5

        coords[:, 1] += np.round(row['rough_j']) - jsize[0]
        coords[:, 0] += np.round(row['rough_i']) - isize[0]

        these_coords = pd.DataFrame()
        these_coords['im_col'] = coords[:, 1]
        these_coords['im_row'] = coords[:, 0]
        these_coords['gcp'] = fid

        coords_all.append(these_coords)

    return pd.concat(coords_all, ignore_index=True)


def _filter_fid_matches(coords_all, measures_cam):
    nfids = len(coords_all.gcp.unique())
    if nfids < len(measures_cam) - 1:
        print(f'Unable to find a transformation with only {nfids} points.')
        inds = []
        for fid in coords_all.gcp.unique():
            inds.append((coords_all['gcp'] == fid).argmax())

        return coords_all.loc[inds]

    combs = list(itertools.combinations(coords_all.index, nfids))

    filtered_combs = [list(c) for c in combs if len(set(coords_all.loc[list(c), 'gcp'].to_list())) == nfids]

    resids = []
    for c in filtered_combs:
        these_meas = coords_all.loc[c].set_index('gcp').join(measures_cam)

        model, inliers = ransac((these_meas[['im_col', 'im_row']].values, these_meas[['j', 'i']].values),
                                AffineTransform, min_samples=3, residual_threshold=10, max_trials=20)
        try:
            resids.append(model.residuals(these_meas[['im_col', 'im_row']].values,
                                          these_meas[['j', 'i']].values).mean())
        except AttributeError:
            resids.append(np.nan)

    return coords_all.loc[filtered_combs[np.nanargmin(resids)]]


def _get_scale(scale, units):
    if units == 'dpi':
        scale = 1 / ((1 / scale) * 25.4)
    else:
        scale = 1 / scale * 1000

    return scale


def _odd(num):
    if num & 1:
        return num
    else:
        return num + 1


def _fix_fiducials(coords, measures_cam):

    model, residuals = _get_residuals(coords, measures_cam)

    missing = ~measures_cam.index.isin(coords['gcp'])
    coords.set_index('gcp', inplace=True)

    for ind, row in measures_cam.loc[missing].iterrows():
        print(f'Predicting location of {ind}')
        x, y = model.inverse(row[['j', 'i']].values).flatten()
        coords.loc[ind, 'im_col'] = x
        coords.loc[ind, 'im_row'] = y

    return coords.reset_index()


def fix_measures_xml(fn_img: str, fn_cam: Union[str, Path, None] = None):
    """
    Use an affine transformation to estimate locations of any missing fiducial markers in an image. Output written to
        Ori-InterneScan/MeasuresIm-{fn_img}.xml

    :param fn_img: the filename of the image.
    :param fn_cam: the Measures file to use. Defaults to Ori-InterneScan/MeasuresCamera.xml.
    """
    if fn_cam is None:
        fn_cam = Path('Ori-InterneScan', 'MeasuresCamera.xml')

    measures_cam = micmac.parse_im_meas(fn_cam).set_index('name')
    measures_img = micmac.parse_im_meas(Path('Ori-InterneScan', f'MeasuresIm-{fn_img}.xml')).set_index('name')

    meas = measures_cam.join(measures_img, lsuffix='_cam', rsuffix='_img').dropna()

    model = AffineTransform()
    model.estimate(meas[['j_img', 'i_img']].values,
                   meas[['j_cam', 'i_cam']].values)

    missing = ~measures_cam.index.isin(measures_img.index)

    for ind, row in measures_cam.loc[missing].iterrows():
        print(f"Predicting location of {ind}")
        x, y = model.inverse(row[['j', 'i']].values).flatten()
        measures_img.loc[ind, 'j'] = x
        measures_img.loc[ind, 'i'] = y

    model.estimate(meas[['j_img', 'i_img']].values,
                   meas[['j_cam', 'i_cam']].values)
    residuals = model.residuals(meas[['j_img', 'i_img']].values,
                                meas[['j_cam', 'i_cam']].values)
    print(f"Mean residual: {residuals.mean():.2f} pixels")

    measures_img.reset_index(inplace=True)
    measures_img.rename(columns={'name': 'gcp', 'j': 'im_col', 'i': 'im_row'}, inplace=True)
    measures_img.dropna(subset=['im_col', 'im_row'], inplace=True)

    micmac.write_measures_im(measures_img, fn_img)


def _get_residuals(meas_img, meas_cam):
    joined = meas_cam.join(meas_img.set_index('gcp'))

    model = AffineTransform()
    est = model.estimate(joined[['im_col', 'im_row']].values, joined[['j', 'i']].values)
    if est:
        return model, model.residuals(joined[['im_col', 'im_row']].values, joined[['j', 'i']].values)
    else:
        raise RuntimeError("Unable to estimate an affine transformation")


def _rotate_meas(meas, angle, pp=None):
    M = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])

    rot = meas.copy()

    if pp is not None:
        shift_j, shift_i = pp.x, pp.y
    else:
        shift_j, shift_i = meas.j.mean(), meas.i.mean()

    rot.j -= shift_j
    rot.i -= shift_i

    rot[['j', 'i']] = rot[['j', 'i']].values.dot(M)

    rot.j += shift_j
    rot.i += shift_j

    return rot


def tfm_measures(meas: pd.DataFrame, pairs: list[tuple],
                 angles: Union[dict, pd.Series], inverted: bool = True) -> tuple[pd.DataFrame, float]:
    """
    Rotate a set of fiducial marker measures so that the angle made by each fiducial marker and the principal point
    is as close to the expected angle as possible.

    :param meas: a DataFrame of fiducial marker locations (as read by micmac.parse_im_meas)
    :param pairs: a list of pairs of co-linear fiducial markers; e.g., [(P1, P2), (P3, P4)]
    :param angles: a Series or dict of fiducial marker names and their angle with respect to the principal point.
        i.e., a mid-side marker on the right side of the frame should have an angle of 0, the fiducial marker in the
        upper right-hand corner should have an angle of 45° (pi / 4), a mid-side marker on the top of the frame should
         have an angle of 90° (pi / 2), and so on.
    :param inverted: the y-axis is inverted.
    :return:
        - **rotated** - a DataFrame of the rotated fiducial marker locations
        - **angle** - the angle by which the markers were rotated
    """
    collinear = [LineString(meas.loc[p, ['j', 'i']].values) for p in pairs]

    for ind, pair in enumerate(pairs):
        meas.loc[pair, ['collim_dist']] = collinear[ind].length

    ppx, ppy = _meas_center(meas, pairs)

    meas['j'] -= ppx
    meas['i'] -= ppy

    if inverted:
        meas['i'] *= -1

    meas['angle'] = np.arctan2(meas['i'], meas['j'])
    meas.loc[meas['angle'] < 0, 'angle'] += 2 * np.pi

    if isinstance(angles, dict):
        if max(angles.values()) > 2 * np.pi:
            angles = pd.Series(angles).apply(np.deg2rad)
        else:
            angles = pd.Series(angles)

    rot_angle = np.arctan2(np.mean(np.sin(meas['angle'] - angles)),
                           np.mean(np.cos(meas['angle'] - angles)))

    rotated = _rotate_meas(meas, rot_angle, pp=Point(0, 0))

    rotated['angle'] = np.arctan2(rotated['i'], rotated['j'])

    return rotated, rot_angle


def _meas_center(meas: pd.DataFrame, pairs: list[tuple]) -> tuple[float, float]:
    collims = [LineString(meas.loc[p, ['j', 'i']].values) for p in pairs]
    pp = MultiPoint([a.intersection(b) for a, b in list(combinations(collims, 2))]).centroid

    return pp.x, pp.y


def _get_rough_locs(meas, img=None):
    # get the rough locations of the corners and mid-side fiducial markers in an image
    scaled = meas.copy()
    scaled['j'] /= scaled.j.max()
    scaled['i'] /= scaled.i.max()

    if img is None:
        rough_x, rough_y = np.meshgrid(np.array([0.075, 0.5, 0.925]), np.array([0.075, 0.5, 0.925]))
        rough_pts = [Point(x, y) for x, y in zip(rough_x.flatten(), rough_y.flatten())]

    else:
        left, right, top, bot = _clean_frame(image.get_rough_frame(img, fact=4), img)

        lr = (right - left)
        tb = (bot - top)

        x_mid = left + lr / 2
        y_mid = top + tb / 2

        left += 0.025 * lr
        right -= 0.025 * lr
        top += 0.025 * tb
        bot -= 0.025 * tb

        scaled['j'] = scaled['j'] * (right - left) + left
        scaled['i'] = scaled['i'] * (bot - top) + top

        rough_x = np.array([left, x_mid, right, left, right, left, x_mid, right])
        rough_y = np.array([top, top, top, y_mid, y_mid, bot, bot, bot])

        rough_pts = [Point(x, y) for x, y in zip(rough_x, rough_y)]

    for ind, row in scaled.iterrows():
        pt = Point(row['j'], row['i'])
        dists = [pt.distance(_pt) for _pt in rough_pts]
        nind = np.argmin(dists)

        meas.loc[ind, 'rough_j'] = rough_x.flatten()[nind]
        meas.loc[ind, 'rough_i'] = rough_y.flatten()[nind]

    return meas


def _clean_frame(frame, img):
    left, right, top, bot = frame

    left = max(0, left)
    right = min(img.shape[1], right)
    top = max(0, top)
    bot = min(img.shape[0], bot)
    return left, right, top, bot


def _corner(size):
    templ = np.zeros((size, size), dtype=np.uint8)
    templ[:int(size/2)+1, int(size/2)+1:] = 255
    return templ


def _box(size):
    templ = np.zeros((size, size), dtype=np.uint8)
    templ[:int(size / 2) + 1, int(size / 2) + 1:] = 255
    templ[int(size / 2) + 1:, :int(size / 2) + 1] = 255

    # now, remove the inner 4 pixels along the vertical axis
    # and the inner 3 pixels along the horizontal axis
    templ[:, int(size/2)-2:-(int(size/2)-1)] = 0
    templ[int(size / 2) - 1:-(int(size / 2) - 1), :] = 0

    return templ


def _inscribe(outer, inner):
    padded = np.zeros(outer.shape)
    pad = int((outer.shape[0] - inner.shape[0]) / 2)
    padded[pad:-pad, pad:-pad] = inner
    return outer - padded


def padded_dot(size: int, disk_size: int) -> NDArray:
    """
    Pad a disk-shaped marker with zeros. Works for, e.g., Zeiss RMK mid-side fiducials.

    :param size: the size of the padded template
    :param disk_size: the half-size of the disk to use
    :return: **padded** -- the disk with a padding of zeros around it
    """
    template = 255 * np.ones((size, size))
    dot = 255 * disk(disk_size)

    return 255 - _inscribe(template, dot)


def inscribed_cross(size: int, cross_size: int, width: int = 3, angle: Union[float, None] = 45) -> NDArray:
    """
    Create a cross-shaped template inscribed inside of a circle for matching fiducial marks.

    :param size: the half-size of the template. Final size will be (2 * size + 1, 2 * size + 1).
    :param cross_size: the size of the cross template to create
    :param width: the width of the cross at the center of the template.
    :param angle: the angle to rotate the template by.
    :return: **template** -- the output template
    """

    circle = 255 * disk(size)
    cross = cross_template(cross_size, width=width, angle=angle, no_border=True)
    cross[cross > 0.8] = 255

    pad = int((circle.shape[0] - cross.shape[0]) / 2)
    padded = np.zeros(circle.shape)
    padded[pad:-pad, pad:-pad] = cross

    return circle - padded


def templates_from_meas(fn_img: Union[str, Path], half_size: int = 100, threshold: bool = False) -> dict:
    """
    Create fiducial templates from points in a MeasuresIm file.

    :param str fn_img: the filename of the image to use. Points for templates will be taken from
        Ori-InterneScan-Measuresim{fn-img}.xml.
    :param int half_size: the half-size of the template to create, in pixels
    :param bool threshold: return binary templates based on otsu thresholding
    :return: **templates** -- a dict of (name, template) pairs for each fiducial marker.
    """
    dir_img = os.path.dirname(fn_img)
    if dir_img == '':
        dir_img = '.'

    bn_img = os.path.basename(fn_img)

    fn_meas = Path(dir_img, 'Ori-InterneScan', f'MeasuresIm-{bn_img}.xml')
    meas_im = micmac.parse_im_meas(fn_meas)

    img = io.imread(fn_img)

    templates = []
    for ind, row in meas_im.iterrows():
        subimg, _, _ = make_template(img, (row.i, row.j), half_size=half_size)
        if threshold:
            subimg = (subimg > filters.threshold_otsu(subimg)).astype(np.uint8)

        templates.append(subimg)

    return dict(zip(meas_im.name.values, templates))


def match_fairchild(fn_img: Union[str, Path], size: int, model: str, data_strip: str, fn_cam: Union[str, Path] = None,
                    dot_size: int = 4, **kwargs) -> float:
    """
    Match the fiducial locations for a Fairchild-style camera (4 fiducial markers markers on the side).

    :param fn_img: the filename of the image to match
    :param size: the size of the marker to match
    :param model: the type of fiducial marker: T11 style with either checkerboard-style markers (T11S) or dot style
        markers (T11D), side + corner dot style markers (T12), or K17 style ("wing" style markers). Must be one of
        [K17, T11S, T11D, T12].
    :param data_strip: the location of the data strip in the image (left, right, top, bot). For T11 style cameras,
        the data strip should be along the left-hand side; for K17 style cameras, the "data strip" (focal length
        indicator) should be on the right-hand side. Be sure to check your images, as the scanned images may be rotated
        relative to this expectation.
    :param fn_cam: the filename of the MeasuresCamera.xml file corresponding to the image
    :param dot_size: the half-size of the dot to use for T11D style fiducial markers (default: 4 -> 9x9)
    :param kwargs: additional keyword arguments to pass to matching.find_fiducials()
    :return: the mean residual between the matches and the fiducial marker locations
    """
    assert model.upper() in ['K17', 'T11S', 'T11D', 'T12'], "model must be one of [K17, T11S, T11D, T12]"
    assert data_strip in ['left', 'right', 'top', 'bot'], "data_strip must be one of [left, right, top, bot]"

    if model.upper in ['K17', 'T11S', 'T11D']:
        fids = [f'P{n}' for n in range(1, 5)]
    else:
        fids = [f'P{n}' for n in range(1, 9)]

    if model.upper() == 'K17':
        templ = _corner(size)
        templates = [templ, np.fliplr(templ), templ.T, np.fliplr(templ).T]

        locs = ['right', 'top', 'left', 'bot']
        angles = [None, np.deg2rad(-90), np.deg2rad(180), np.deg2rad(90)]
        ldict = dict(zip(locs, angles))
    else:
        if model.upper() == 'T11S':
            templ = _box(size)
            templates = [templ, templ, np.fliplr(templ), np.fliplr(templ)]
        elif model.upper() == 'T11D':
            templ = padded_dot(size, dot_size)
            templates = [templ, templ, np.fliplr(templ), np.fliplr(templ)]
        elif model.upper() == 'T12':
            templates = 8 * [padded_dot(size, dot_size)]

        locs = ['left', 'top', 'right', 'bot']
        angles = [None, np.deg2rad(-90), np.deg2rad(180), np.deg2rad(90)]
        ldict = dict(zip(locs, angles))

    tdict = dict(zip(fids, templates))
    angle = ldict[data_strip]

    return find_fiducials(fn_img, tdict, fn_cam=fn_cam, angle=angle, **kwargs)


def _zeiss_corner(size):
    return 4 * [cross_template(size, no_border=True)]


def _zeiss_midside(size, dot_size):
    templ = padded_dot(size, dot_size)
    return 4 * [templ]


def match_zeiss_rmk(fn_img: Union[str, Path], size: int, dot_size: int,
                    data_strip: str = 'left', fn_cam: Union[str, Path, None] = None,
                    corner_size: Union[int, None] = None, **kwargs) -> float:
    """
    Match the fiducial locations for a Zeiss RMK-style camera (4 dot-shaped markers on the side, possibly 4 cross-shaped
    markers in the corners).

    :param str fn_img: the filename of the image to match
    :param int size: the size of the marker to match
    :param int dot_size: the size of the dot marker to match
    :param str data_strip: the location of the data strip in the image (left, right, top, bot). Most calibration reports
        assume the data strip is along the left-hand side, but scanned images may be rotated relative to this.
    :param str fn_cam: the filename of the MeasuresCamera.xml file corresponding to the image
    :param int corner_size: the size of the corner markers (default: do not find corner markers)
    :param kwargs: additional keyword arguments to pass to matching.find_fiducials()
    :return: the mean residual between the matches and the fiducial marker locations
    """
    assert data_strip in ['left', 'right', 'top', 'bot'], "data_strip must be one of [left, right, top, bot]"

    if corner_size is not None:
        fids = [f'P{n}' for n in range(1, 9)]
        ctempl = _zeiss_corner(corner_size)
        stempl = _zeiss_midside(size, dot_size)
        templates = ctempl + stempl
    else:
        fids = [f'P{n}' for n in range(1, 5)]
        templates = _zeiss_midside(size, dot_size)

    if data_strip == 'left':
        angle = None
    elif data_strip == 'top':
        angle = np.deg2rad(-90)
    elif data_strip == 'right':
        angle = np.deg2rad(180)
    else:
        angle = np.deg2rad(90)

    tdict = dict(zip(fids, templates))
    return find_fiducials(fn_img, tdict, fn_cam=fn_cam, angle=angle, **kwargs)


def _wild_corner(size, model, circle_size=None, ring_width=7, width=3, gap=None, vgap=None, dot_size=None, pad=10):

    target_angle = 45
    if model.upper() in ['RC5']:
        if circle_size is None:
            circle_size = _odd(int(1.2 * size))  # approximate but probably good enough
        template = inscribed_cross(circle_size, size, angle=45)
        template = np.pad(template, 20)

    elif model.upper() in ['RC5A', 'RC8']:
        template = cross_template(size, width=width, angle=45, no_border=True)

        rows, cols = template.shape

        if gap is not None:
            half_r = int((rows - 1) / 2)
            half_c = int((rows - 1) / 2)

            if vgap is None:
                half_h = int((gap - 1) / 2)
                half_v = int((gap - 1) / 2)
            else:
                half_h = int((gap - 1) / 2)
                half_v = int((vgap - 1) / 2)

            template[half_r - half_v:half_r + half_v + 1, :] = 0
            template[:, half_c - half_h:half_c + half_h + 1] = 0

        if dot_size is not None:
            template += padded_dot(size, dot_size)

        template[template > 0.8] = 255
    else:
        template = wagon_wheel(size, width=width, circle_size=circle_size, circle_width=ring_width, angle=target_angle)

    return np.pad(template, pad)


def _wild_midside(size, model, circle_size, ring_width):

    if model.upper() in ['RC5', 'RC8']:
        target_angle = 45
    else:
        target_angle = None

    template = wagon_wheel(size, width=3, circle_size=circle_size, circle_width=ring_width, angle=target_angle)

    return template


def match_wild_rc(fn_img: Union[str, Path], size: int, model: str, data_strip: str = 'left',
                  fn_cam: Union[str, Path] = None, width: int = 3, circle_size: Union[int, None] = None,
                  ring_width: int = 7, gap: int = 9, vgap: Union[int, None] = None, dot_size: Union[int, None] = None,
                  pad: int = 10, **kwargs) -> float:
    """
    Match the fiducial locations for a Wild RC-style camera (4 cross/bulls-eye markers in the corner, possibly
    4 bulls-eye markers along the sides).

    :param fn_img: the filename of the image to match
    :param size: the size of the marker to match
    :param model: whether the camera is an RC5/RC8 (4 corner markers) or RC10-style (corner + midside markers)
    :param data_strip: the location of the data strip in the image (left, right, top, bot). Most calibration reports
        assume the data strip is along the left-hand side, but scanned images may be rotated relative to this.
    :param fn_cam: the filename of the MeasuresCamera.xml file corresponding to the
        image (default: Ori-InterneScan/MeasuresCamera.xml)
    :param width: the thickness of the cross template, in pixels
    :param circle_size: the size of the circle in which to inscribe the cross-shaped marker (default: no circle)
    :param ring_width: the width of the ring if the marker(s) are a cross inscribed with a ring. Only used for RC10
        models.
    :param gap: the width, in pixels, of the gap in the middle of the cross
    :param vgap: the height, in pixels, of the gap in the middle of the cross (default: same as gap)
    :param dot_size: the half-size, in pixels, of the dot in the middle of the cross (default: no dot)
    :param pad: the size of the padding around the outside of the cross to include
    :param kwargs: additional keyword arguments to pass to matching.find_fiducials()
    :return: the mean residual betweent he matches and the fiducial marker locations
    """
    assert model.upper() in ['RC5', 'RC5A', 'RC8', 'RC10'], "model must be one of [RC5, RC5A, RC8, RC10]"
    assert data_strip in ['left', 'right', 'top', 'bot'], "data_strip must be one of [left, right, top, bot]"

    if model.upper() in ['RC5', 'RC5A', 'RC8']:
        fids = [f'P{n}' for n in range(1, 5)]
        templates = 4 * [_wild_corner(size, model, circle_size, ring_width, width=width,
                                      gap=gap, vgap=vgap, dot_size=dot_size, pad=pad)]
    else:
        fids = [f'P{n}' for n in range(1, 9)]
        stempl = _wild_midside(size, model, circle_size, ring_width)
        ctempl = _wild_corner(size, model, circle_size, ring_width, width=width,
                              gap=gap, vgap=vgap, dot_size=dot_size, pad=pad)
        templates = 4 * [ctempl] + 4 * [stempl]

    locs = ['left', 'top', 'right', 'bot']
    angles = [None, np.deg2rad(-90), np.deg2rad(180), np.deg2rad(90)]
    ldict = dict(zip(locs, angles))

    tdict = dict(zip(fids, templates))
    return find_fiducials(fn_img, tdict, fn_cam=fn_cam, angle=ldict[data_strip], **kwargs)


def cross_template(shape: int, width: int = 3,
                   angle: Union[float, None] = None, no_border: bool = False) -> NDArray:
    """
    Create a cross-shaped template for matching Réseau or fiducial marks.

    :param int shape: the output shape of the template
    :param int width: the width of the cross at the center of the template
    :param float angle: the angle to rotate the template by (default: no rotation)
    :param bool no_border: do not include a border around the cross
    :return: **cross** -- the cross template
    """
    if isinstance(shape, int):
        rows = shape
        cols = shape
    else:
        rows, cols = shape

    if angle is not None:
        rows *= (np.sin(np.deg2rad(angle)) + np.cos(np.deg2rad(angle)))
        rows = int(np.round(rows))
        if rows % 2 == 0:
            rows += 1

        cols *= (np.sin(np.deg2rad(angle)) + np.cos(np.deg2rad(angle)))
        cols = int(np.round(cols))
        if cols % 2 == 0:
            cols += 1

    half_r = int((rows - 1) / 2)
    half_c = int((cols - 1) / 2)
    half_w = int((width - 1) / 2)

    cross = np.zeros((rows, cols))
    if not no_border:
        cross[half_r - half_w - 1:half_r + half_w + 2:width + 1, :] = 2
        cross[:, half_c - half_w - 1:half_c + half_w + 2:width + 1] = 2

    cross[half_r - half_w:half_r + half_w + 1, :] = 1
    cross[:, half_c - half_w:half_c + half_w + 1] = 1

    if angle is None:
        return cross
    else:
        cross = np.round(ndimage.rotate(cross, angle, reshape=False))
        rs, = np.where(cross.sum(axis=1) > 0)
        cs, = np.where(cross.sum(axis=0) > 0)

        return cross[rs[0]+1:rs[-1], cs[0]+1:cs[-1]]


def find_crosses(img: NDArray, cross: NDArray) -> pd.DataFrame:
    """
    Find all cross markers in an image.

    :param img: the image
    :param cross: the cross template to use
    :return: **grid_df** -- a dataframe of marker locations and offsets
    """

    sub_coords = []
    simgs, top_inds, left_inds = image.splitter(img, (4, 8), overlap=4*cross.shape[0])

    for ind, simg in enumerate(simgs):
        img_inv = img.max() - simg
        res = cv2.matchTemplate(img_inv.astype(np.uint8), cross.astype(np.uint8), cv2.TM_CCORR_NORMED)

        these_coords = peak_local_max(res, min_distance=int(1.5*cross.shape[0]),
                                      threshold_abs=np.quantile(res, 0.6)).astype(np.float64)

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


def match_reseau_grid(img: NDArray, coords: NDArray, cross: NDArray) -> pd.DataFrame:
    """
    Find the best match for each KH-9 mapping camera Réseau grid point, given a list of potential matches.

    :param img: the image to use
    :param coords: the coordinates of the potential matches
    :param cross: the cross template to use.
    :return: **grid_df** -- a DataFrame of grid locations and match points
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
        grid_df.loc[ii, 'gcp'] = f"GCP_{pr[0]}_{pr[1]}"
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
    fixed = image.nanmedian_filter(subimg, footprint=disk(7))
    subimg[np.isnan(subimg)] = fixed[np.isnan(subimg)]
    return subimg.astype(np.uint8)


def remove_crosses(fn_img: Union[str, Path], nproc: int = 1) -> None:
    """
    Remove the Réseau marks from a KH-9 image before re-sampling.

    :param fn_img: the image filename
    :param nproc: the number of subprocesses to use
    """
    fn_meas = Path('Ori-InterneScan', f"MeasuresIm-{fn_img}.xml")
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
    fixed = image.nanmedian_filter(subim, footprint=disk(7))
    subim[np.isnan(subim)] = fixed[np.isnan(subim)]
    img[int(pt[0]) - row_[0]:int(pt[0]) + row_[1] + 1, int(pt[1]) - col_[0]:int(pt[1]) + col_[1] + 1] = subim.astype(
        np.uint8)

    if nproc == 1:
        for ind, row in gcps.loc[1:].iterrows():
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


def find_reseau_grid(fn_img: Union[str, Path],
                     csize: int = 361, return_val: bool = False) -> Union[pd.DataFrame, None]:
    """
    Find the locations of the Réseau marks in a scanned KH-9 image. Locations are saved
    to Ori-InterneScan/MeasuresIm-{fn_img}.xml.

    :param fn_img: the image filename.
    :param csize: the size of the cross template (default: 361 -> 361x361)
    :param return_val: return a pandas DataFrame of the Réseau mark locations
    :return: **gcps_df** (*pandas.DataFrame*) -- a DataFrame of the Réseau mark locations (if return_val=True).
    """
    img = io.imread(fn_img)
    img_lowres = resample.downsample(img, fact=10)

    cross = cross_template(csize, width=3)
    cross[cross > 1] = 0
    cross *= 255

    fig = plt.figure(figsize=(7, 12))
    ax = fig.add_subplot(111)
    ax.imshow(img_lowres, cmap='gray', extent=(0, img.shape[1], img.shape[0], 0))
    ax.set_xticks([])
    ax.set_yticks([])

    print(f"Finding grid points in {fn_img}...")
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

    xdiff = ux - image.nanmedian_filter(ux, footprint=disk(3))
    ydiff = uy - image.nanmedian_filter(uy, footprint=disk(3))

    xout = _outlier_filter(xdiff, n=5)
    yout = _outlier_filter(ydiff, n=5)

    outliers = np.logical_or.reduce([~inliers, xout.flatten(), yout.flatten()])

    ux[outliers.reshape(23, 47)] = np.nan
    uy[outliers.reshape(23, 47)] = np.nan

    ux[outliers.reshape(23, 47)] = image.nanmedian_filter(ux, footprint=disk(1))[outliers.reshape(23, 47)]
    uy[outliers.reshape(23, 47)] = image.nanmedian_filter(uy, footprint=disk(1))[outliers.reshape(23, 47)]

    grid_df.loc[outliers, 'match_j'] = dst[outliers, 0] + ux.flatten()[outliers]
    grid_df.loc[outliers, 'match_i'] = dst[outliers, 1] + uy.flatten()[outliers]

    grid_df.loc[np.isnan(grid_df['match_j']), 'match_j'] = dst[np.isnan(grid_df['match_j']), 0]
    grid_df.loc[np.isnan(grid_df['match_i']), 'match_i'] = dst[np.isnan(grid_df['match_i']), 1]

    grid_df['im_row'] = grid_df['match_i']
    grid_df['im_col'] = grid_df['match_j']

    grid_df['dj'] = grid_df['match_j'] - dst[:, 0]
    grid_df['di'] = grid_df['match_i'] - dst[:, 1]

    os.makedirs('match_imgs', exist_ok=True)
    print(f"Grid points found for {fn_img}.")
    print(f"Mean x residual: {grid_df.loc[~outliers, 'dj'].abs().mean():.2f} pixels")
    print(f"Mean y residual: {grid_df.loc[~outliers, 'di'].abs().mean():.2f} pixels")
    print(f"Mean residual: {grid_df.loc[~outliers, 'resid'].mean():.2f} pixels")

    ax.quiver(grid_df.match_j, grid_df.match_i, grid_df.dj, grid_df.di, color='r')
    ax.plot(grid_df.match_j[outliers], grid_df.match_i[outliers], 'b+')

    this_out = os.path.splitext(fn_img)[0]
    fig.savefig(Path('match_imgs', this_out + '_matches.png'), bbox_inches='tight', dpi=200)

    micmac.write_measures_im(grid_df, fn_img)

    if return_val:
        return grid_df
    else:
        return None


def wagon_wheel(size: int, width: int = 3, mult: Union[int, float] = 255,
                circle_size: Union[int, None] = None, circle_width: Union[int, None] = None,
                angle: Union[float, None] = None) -> NDArray:
    """
    Creates a template in the shape of a "wagon wheel" (a cross inscribed in a ring).

    :param size: the width (and height) of the template, in pixels
    :param width: the width/thickness of the cross, in pixels
    :param mult: a multiplier to use for the template
    :param circle_size: the size of the circle to inscribe the cross into (default: same as cross size)
    :param circle_width: the width of the ring to inscribe the cross into (default: same as cross width)
    :param angle: the angle by which to rotate the cross (default: do not rotate)
    :return: **template** the wagon wheel template
    """
    cross = cross_template(size, width, angle=angle)
    cross[cross > 0.8] = 1

    if circle_size is None:
        templ = disk(int(size / 2))
        padded = np.zeros(templ.shape, dtype=templ.dtype)
        padded[width:-width, width:-width] = disk(int((size - 2 * width) / 2))
    else:
        templ = disk(int(circle_size / 2))
        padded = np.zeros(templ.shape, dtype=templ.dtype)
        if circle_width is None:
            padded[width:-width, width:-width] = disk(int((circle_size - 2 * width) / 2))
        else:
            padded[circle_width:-circle_width,
                   circle_width:-circle_width] = disk(int((circle_size - 2 * circle_width) / 2))

    templ -= padded

    if circle_size is None:
        templ += cross.astype(templ.dtype)
    else:
        padded = np.zeros(cross.shape)
        pad = int((cross.shape[0] - circle_size) / 2)
        padded[pad:-pad, pad:-pad] = templ
        padded += cross

        templ = padded

    templ[templ > 1] = 1

    return mult * templ


# all credit to joe kennedy for the name of this function.
def ocm_show_wagon_wheels(img: NDArray, size: int, width: int = 3,
                          img_border: Union[tuple[int, int], None] = None) -> NDArray:
    """
    Find all "wagon wheel" markers in an image.

    :param img: the image
    :param size: the size of the marker (in pixels)
    :param width: the width/thickness of the cross, in pixels
    :param img_border: the approximate top and bottom rows of the image frame. If not set,
        calls get_rough_frame() on the image.
    :return: **coords** -- an Nx2 array of the location of the detected markers.
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


def find_rail_marks(img: NDArray, marker: NDArray) -> NDArray:
    """
    Find all rail marks along the bottom edge of a KH-4 style image.

    :param img: the image to find the rail marks in.
    :param marker: the marker template to use for matching
    :return: **coords** -- Nx2 array of the location (row, col) of the detected markers.
    """
    left, right, top, bot = image.get_rough_frame(img)
    img_lowres = resample.downsample(img, fact=10)

    res = cv2.matchTemplate(img_lowres.astype(np.uint8), marker.astype(np.uint8), cv2.TM_CCORR_NORMED)

    coords = peak_local_max(res, threshold_rel=np.percentile(res, 99.9) / res.max(),
                            min_distance=10).astype(np.float64)
    coords += marker.shape[0] / 2 - 0.5
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
        valid = np.abs(diff - np.median(diff)) < 4 * register.nmad(diff)

        nout = prev_valid - np.count_nonzero(valid)
        prev_valid = np.count_nonzero(valid)

    return valid


def notch_template(size: int) -> NDArray:
    """
    Create a notch-shaped ("^") template.

    :param size: the size of the template, in pixels
    :return: **template** -- the notch template
    """
    template = np.zeros((size, size), dtype=np.uint8)
    template[-1, :] = 1
    for ind in range(1, int(size/2) + 1):
        template[-ind-1, ind:-ind] = 1
    return 255 * template


def find_kh4_notches(img: NDArray, size: int = 101) -> NDArray:
    """
    Find all 4 notches along the top of a KH-4 style image.

    :param img: the image.
    :param size: the size of the notch template to use.
    :return: **coords** -- a 4x2 array of notch locations
    """
    left, right, top, bot = image.get_rough_frame(img)

    templ = notch_template(size)

    pcts = np.array([0.03, 0.5, 0.6, 0.97])

    search_j = (left + pcts * (right - left)).astype(int)
    search_i = top * np.ones(4).astype(int)

    matches = []

    for jj, ii in zip(search_j, search_i):
        subimg, _, _ = make_template(img, [ii, jj], 4 * size)
        # res, this_i, this_j = find_match(subimg, templ, how='max', eq=False)
        res = cv2.matchTemplate(subimg.astype(np.uint8), templ.astype(np.uint8), cv2.TM_CCORR_NORMED)

        try:
            this_i, this_j = peak_local_max(res, min_distance=50, num_peaks=1)[0]
        except IndexError:
            this_i = this_j = np.nan

        this_i += (subimg.shape[0] - res.shape[0]) / 2
        this_j += (subimg.shape[1] - res.shape[1]) / 2

        matches.append((this_i - 4*size + ii, this_j - 4*size + jj))

    return np.array(matches)


######################################################################################################################
# tools for matching gcps
######################################################################################################################
def make_template(img: NDArray, pt: tuple[int, int], half_size: int) -> tuple[NDArray, list, list]:
    """
    Return a sub-section of an image to use for matching.

    :param img: the image from which to create the template
    :param pt: the (row, column) center of the template
    :param half_size: the half-size of the template; template size will be 2 * half_size + 1
    :return:
        - **template** -- the template
        - **row_inds** -- the number of rows above/below the center of the template
        - **col_inds** -- the number of columns left/right of the center of the template
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


def find_match(img: NDArray, template: NDArray, how: str = 'min', eq: bool = True) -> tuple[NDArray, float, float]:
    """
    Given an image and a template, find a match using openCV's normed cross-correlation.

    :param img: the image to find a match in
    :param template: the template to use for matching
    :param how: determines whether the match is the minimum or maximum correlation
    :param eq: use a rank equalization filter before matching the templates
    :return:
        - **res** -- the correlation image
        - **match_i** -- the row location of the match
        - **match_j** -- the column location of the match
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


def _random_points(mask, npts):
    pctmask = np.count_nonzero(mask) / mask.size

    xx = np.random.uniform(0, mask.shape[1], int(npts / pctmask)).astype(int)
    yy = np.random.uniform(0, mask.shape[0], int(npts / pctmask)).astype(int)

    masked = mask[yy, xx] == 0

    return xx[~masked], yy[~masked]


def _match_grid(refgeo, spacing, srcwin):

    jj, ii = np.meshgrid(np.arange(srcwin, spacing * np.ceil((refgeo.shape[1]-srcwin) / spacing) + 1, spacing),
                         np.arange(srcwin, spacing * np.ceil((refgeo.shape[0]-srcwin) / spacing) + 1, spacing))

    return jj.astype(int).flatten(), ii.astype(int).flatten()


def find_matches(tfm_img: NDArray, refgeo: gu.Raster, mask: NDArray, points: Union[gpd.GeoDataFrame, None] = None,
                 initM: Union[ProjectiveTransform, None] = None, strategy: str = 'grid',
                 spacing: int = 200, srcwin: int = 60, dstwin: int = 600) -> pd.DataFrame:
    """
    Find matches between two images using normalized cross-correlation template matching. If point locations are not
    given, generates a two-dimensional grid of evenly spaced points.

    :param tfm_img: the image to use for matching.
    :param refgeo: the reference image to use for matching.
    :param mask: a mask indicating areas that should be used for matching.
    :param points: a GeoDataFrame of point locations
    :param initM: the model used for transforming the initial, non-georeferenced image.
    :param strategy: strategy for generating points. Must be one of 'grid' or 'random'. Note that if
        'random' is used, density is the approximate number of points, rather than the distance between
        grid points
    :param spacing: the grid spacing, in pixels
    :param srcwin: the half-size of the template window.
    :param dstwin: the half-size of the search window.
    :return: **gcps** -- a DataFrame with GCP locations, match strength, and other information.
    """
    assert strategy in ['grid', 'random', 'chebyshev'], f"{strategy} must be one of [grid, random]"

    match_pts = []
    z_corrs = []
    peak_corrs = []

    if points is None:
        if strategy == 'grid':
            jj, ii = np.array(_match_grid(refgeo, spacing, srcwin))
        elif strategy == 'random':
            jj, ii = _random_points(mask, spacing)
        elif strategy == 'chebyshev':
            pass
    else:
        jj, ii = points.search_j, points.search_i

    search_pts = []

    for _i, _j in zip(ii, jj):
        search_pts.append((_j, _i))
        match, z_corr, peak_corr = do_match(tfm_img, refgeo.data, mask, (int(_i), int(_j)), srcwin, dstwin)
        match_pts.append(match)
        z_corrs.append(z_corr)
        peak_corrs.append(peak_corr)

    search_pts = np.array(search_pts)
    _dst = np.array(match_pts)

    if points is None:
        gcps = gpd.GeoDataFrame()
    else:
        gcps = points.copy()

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


def do_match(dest_img: NDArray, ref_img: NDArray, mask: NDArray, pt: tuple[int, int],
             srcwin: int, dstwin: int) -> tuple[tuple, float, float]:
    """
    Find a match between two images using normalized cross-correlation template matching.

    :param dest_img: the image to search for the matching point in.
    :param ref_img: the reference image to use for matching.
    :param mask: a mask indicating areas that should be used for matching.
    :param pt: the index (i, j) to search for a match for.
    :param srcwin: the half-size of the template window.
    :param dstwin: the half-size of the search window.
    :return:
        - **match_pt** (*tuple*) -- the matching point (j, i) found in dest_img
        - **z_corr** (*float*) -- number of standard deviations (z-score) above other potential matches
        - **peak_corr** (*float*) -- the correlation value of the matched point
    """
    _i, _j = pt
    submask, _, _ = make_template(mask, pt, srcwin)

    if submask.size == 0:
        return (np.nan, np.nan), np.nan, np.nan
    elif np.count_nonzero(submask) / submask.size < 0.05:
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


def get_matches(img1: NDArray, img2: NDArray, mask1: Union[NDArray, None] = None, mask2: Union[NDArray, None] = None,
                dense: bool = False, npix: int = 100, nblocks: Union[int, None] = None) -> tuple[tuple, tuple, list]:
    """
    Return keypoint matches found using openCV's ORB implementation.

    :param img1: the first image to match
    :param img2: the second image to match
    :param mask1: a mask to use for the first image. (default: no mask)
    :param mask2: a mask to use for the second image. (default: no mask)
    :param dense: compute matches over sub-blocks (True) or the entire image (False).
    :param npix: the block size (in pixels) to divide the image into, if doing dense matching.
    :param nblocks: the number of blocks to divide the image into. If set, overrides value given by npix.
    :return:
        - **keypoints** -- the keypoint locations for the first and second image.
        - **descriptors** -- the descriptors for the first and second image.
        - **matches** -- a list of matching keypoints between the first and second image
    """
    if dense:
        if np.any(np.array(img1.shape) < 100) or np.any(np.array(img2.shape) < 100):
            kp1, des1 = get_dense_keypoints(img1.astype(np.uint8), mask1, nblocks=1, return_des=True)
            kp2, des2 = get_dense_keypoints(img2.astype(np.uint8), mask2, nblocks=1, return_des=True)
        else:
            kp1, des1 = get_dense_keypoints(img1.astype(np.uint8), mask1, npix=npix, nblocks=nblocks, return_des=True)
            kp2, des2 = get_dense_keypoints(img2.astype(np.uint8), mask2, npix=npix, nblocks=nblocks, return_des=True)
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


def _dense_skimage(split_img: list, oy: list, ox: list, return_des: bool, detector_kwargs: dict) -> tuple[list, ...]:
    orb = ORB(**detector_kwargs)

    keypoints = []
    descriptors = []

    for ind, img in enumerate(split_img):
        try:
            orb.detect_and_extract(img)
            kp, des = orb.keypoints, orb.descriptors

            descriptors.append(des)

            kp[:, 1] += ox[ind]
            kp[:, 0] += oy[ind]

            keypoints.append(kp)
        except RuntimeError as e:
            if "ORB found no features." in e.args[0]:
                continue
            else:
                raise e

    keypoints = np.concatenate(keypoints, axis=0)
    descriptors = np.concatenate(descriptors, axis=0)

    if return_des:
        return keypoints, descriptors
    else:
        return keypoints


def _dense_opencv(split_img: list, oy: list, ox: list,
                  split_msk: list, return_des: bool, detector_kwargs: dict) -> tuple[list, ...]:
    orb = cv2.ORB_create(**detector_kwargs)

    keypoints = []
    descriptors = []

    for ind, img in enumerate(split_img):
        kp, des = orb.detectAndCompute(img, mask=split_msk[ind])

        for ds in des:
            descriptors.append(ds)

        for p in kp:
            p.pt = p.pt[0] + ox[ind], p.pt[1] + oy[ind]
            keypoints.append(p)

    if return_des:
        return keypoints, np.array(descriptors)
    else:
        return keypoints


def get_dense_keypoints(img: NDArray, mask: NDArray, npix: int = 100, nblocks: int = None, return_des: bool = False,
                        use_skimage: bool = False, detector_kwargs: dict = {}) -> tuple[list, ...]:
    """
    Find ORB keypoints by dividing an image into smaller parts.

    :param img: the image to use.
    :param mask: a mask to use for keypoint generation.
    :param npix: the block size (in pixels) to divide the image into.
    :param nblocks: the number of blocks to divide the image into. If set, overrides value given by npix.
    :param return_des: return the keypoint descriptors, as well
    :param use_skimage: use the scikit-image implementation of ORB rather than OpenCV
    :param detector_kwargs: additional keyword arguments to pass when creating the ORB detector. For details,
        see the documentation for cv2.ORB_create or skimage.feature.ORB.
    :return:
        - **keypoints** -- a list of keypoint locations
        - **descriptors** -- if requested, a list of keypoint descriptors.
    """

    if nblocks is None:
        x_tiles = max(1, np.floor(img.shape[1] / npix).astype(int))
        y_tiles = max(1, np.floor(img.shape[0] / npix).astype(int))
    else:
        x_tiles = nblocks
        y_tiles = nblocks

    olap = int(max(0.25 * img.shape[1]/x_tiles, 0.25 * img.shape[0]/y_tiles))

    split_img, oy, ox = image.splitter(img, (y_tiles, x_tiles), overlap=olap)
    if mask is not None:
        split_msk, _, _ = image.splitter(mask, (y_tiles, x_tiles), overlap=olap)
    else:
        split_msk = [None] * len(split_img)

    if use_skimage:
        return _dense_skimage(split_img, oy, ox, return_des, detector_kwargs)
    else:
        return _dense_opencv(split_img, oy, ox, split_msk, return_des, detector_kwargs)


def match_halves(left: NDArray, right: NDArray, overlap: int, block_size: int = None) -> EuclideanTransform:
    """
    Find a transformation to join the left and right halves of an image scan.

    :param left: the left-hand image scan.
    :param right: the right-hand image scan.
    :param overlap: the estimated overlap between the two images, in pixels.
    :param block_size: the number of rows each sub-block should cover. Defaults to overlap.
    :return: **model** -- the estimated Euclidean transformation between the two image halves.
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

    models = []
    inliers = []

    for ii in range(20):
        mod, inl = ransac((np.array(src_pts), np.array(dst_pts)), EuclideanTransform,
                          min_samples=10, residual_threshold=0.2, max_trials=2500)
        models.append(mod)
        inliers.append(inl)

    num_inliers = [np.count_nonzero(inl) for inl in inliers]
    best_ind = np.argmax(num_inliers)

    print(f"{num_inliers[best_ind]} tie points found")
    return models[best_ind], num_inliers[best_ind]
