"""
spymicmac.image is a collection of tools for working with images.
"""
import os
from glob import glob
from itertools import chain
from skimage import exposure, morphology, io, filters
from skimage.morphology import disk
from skimage.feature import peak_local_max
from skimage.transform import warp
from scipy import ndimage
import numpy as np
from numba import jit
from pybob.image_tools import nanmedian_filter
from spymicmac import matching, resample


######################################################################################################################
# image filtering tools
######################################################################################################################
@jit
def nanstd(a):
    return np.nanstd(a)


def highpass_filter(img):
    """
    Subtract a low-pass from an image, to return a highpass filter.

    :param array-like img: the image to filter.
    :return: **highpass** (*array-like*) -- the highpass-filtered image
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
    """
    Split an image into (m, n) blocks with a given overlap.

    :param array-like img: the image to split
    :param tuple nblocks: the number of blocks to create along each axis (m, n)
    :param int overlap: the number of pixels to overlap each block. (default: 0)

    :return:
        - **blocks** (*list*) -- a list of the image blocks created
        - **top_inds** (*list*) -- a list of the original row index of the top edge of each block
        - **left_inds** (*list*) -- a list of the original column index of the left edge of each block
    """
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


def _subimg_offsets(split, shape, overlap=0):
    ims_x = np.array([s.shape[1] - 2 * overlap for s in split])
    ims_y = np.array([s.shape[0] - 2 * overlap for s in split])

    rel_x = np.cumsum(ims_x.reshape(shape), axis=1)
    rel_y = np.cumsum(ims_y.reshape(shape), axis=0)

    rel_x = np.concatenate((np.zeros((shape[0], 1)), rel_x[:, :-1]), axis=1)
    rel_y = np.concatenate((np.zeros((1, shape[1])), rel_y[:-1, :]), axis=0)

    return rel_x.astype(int), rel_y.astype(int)


def stretch_image(img, scale=(0.0, 1.0), mult=255, imgmin=0, outtype=np.uint8, mask=None):
    """
    Apply a linear stretch to an image by clipping and stretching to quantiles.

    :param array-like img: the image to stretch.
    :param tuple scale: the minimum and maximum quantile to stretch to. (default: (0, 1) - the minimum/maximum values
        of the image)
    :param int|float mult: a multiplier to scale the result to. (default: 255)
    :param int|float imgmin: the minimum value in the output image (default: 0)
    :param numpy.dtype outtype: the numpy datatype to return the stretched image as. (default: np.uint8)
    :param array-like mask: a mask of pixels to ignore when calculating quantiles.
    :return: **stretched** (*array-like*) -- the stretched image.
    """
    if mask is None:
        maxval = np.nanquantile(img, max(scale))
        minval = np.nanquantile(img, min(scale))
    else:
        maxval = np.nanquantile(img[mask], max(scale))
        minval = np.nanquantile(img[mask], min(scale))

    img[img > maxval] = maxval
    img[img < minval] = minval

    return (mult * (img - minval) / (maxval - minval + imgmin)).astype(outtype)


def remove_scanner_stripes(img, dtype=np.uint8, scan_axis=1):
    """
    Remove horizontal (or vertical) stripes from an image.

    :param array-like img: the image to remove stripes from.
    :param numpy.dtype dtype: the original datatype of the image.
    :param int scan_axis: the axis corresponding to the direction of the stripes. A scan_axis of 1 corresponds to
        horizontal stripes (the default), while a scan_axis of 0 corresponds to vertical stripes.

    :return: **destriped**: the original image with the stripes (mostly) removed.
    """
    assert scan_axis in [0, 1], "scan_axis corresponds to image axis of scan direction [0, 1]"
    if scan_axis == 0:
        outimg = img.astype(np.float32) + \
                 ((np.median(img.flatten()) - np.median(img, axis=scan_axis)).reshape(1, -1) * np.ones(img.shape))
    else:
        outimg = img.astype(np.float32) + \
                 (np.ones(img.shape) * (np.median(img.flatten()) - np.median(img, axis=scan_axis)).reshape(-1, 1))

    return stretch_image(outimg, outtype=dtype)


def contrast_enhance(fn_img, mask_value=None, qmin=0.02, qmax=0.98, gamma=1.25, disksize=3, imgmin=0):
    """
    Enhance image contrast in a three-step process. First, the image is processed with a median filter to reduce
    noise. Next, a linear contrast stretch is applied, and finally, a gamma adjustment is applied.

    :param str fn_img: the image filename.
    :param int|float mask_value: a mask value to use when filtering the image.
    :param float qmin: the minimum quantile to use for the linear contrast stretch (default: 0.02)
    :param float qmax: the maximum quantile to use for the linear contrast stretch (default: 0.98)
    :param float gamma: the value to use for the gamma adjustment
    :param int disksize: the filter disk size (input to skimage.morphology.disk; default: 3)
    :param int|float imgmin: the minimum value in the output image (default: 0)
    :return: **enhanced** (*array-like*) -- the contrast-enhanced image.
    """
    img = io.imread(fn_img)
    if mask_value is not None:
        img = img.astype(np.float32)
        img[img == mask_value] = np.nan

    filt = nanmedian_filter(img, footprint=disk(disksize))

    stretch = stretch_image(filt, scale=(qmin, qmax), imgmin=imgmin)
    gamma = exposure.adjust_gamma(stretch, gamma=gamma)

    if mask_value is not None:
        gamma[img == mask_value] = mask_value

    return gamma


def make_binary_mask(img, mult_value=255, erode=0, mask_value=0):
    """
    Create a binary mask for an image based on a given mask value. Values equal to mask_value will be given a value
    of 0 in the mask, while all other values will be set equal to mult_value.

    :param array-like img: the image to create a mask for.
    :param int|float mult_value: the value indicating a non-masked value (default: 255).
    :param int erode: the size of the erosion operation to apply (default: 0).
    :param int mask_value: the value to mask in the image (default: 0).
    :return: **mask** (*array-like*) the binary mask.
    """
    _mask = mult_value * np.ones(img.shape, dtype=np.uint8)
    if np.isfinite(mask_value):
        _mask[img == mask_value] = 0
    else:
        _mask[np.isnan(img)] = 0

    if erode > 0:
        erode_mask = morphology.binary_erosion(_mask, footprint=morphology.disk(erode))
        _mask[~erode_mask] = 0

    return _mask


def balance_image(img):
    """
    Apply contrast-limited adaptive histogram equalization (CLAHE) on an image, then apply a de-noising filter.

    :param array-like img: the image to balance.
    :return: **img_filt** (*array-like*) -- the balanced, filtered image.
    """
    img_eq = (255 * exposure.equalize_adapthist(img)).astype(np.uint8)
    img_filt = filters.median(img_eq, footprint=disk(1))
    return img_filt


# thanks to SO user Jamie for this answer
# https://stackoverflow.com/a/14314054
def _moving_average(a, n=5):
    ret = np.cumsum(a)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


def get_rough_frame(img):
    """
    Find the rough location of an image frame/border.

    :param array-like img: the image to find a border for
    :return: **xmin**, **xmax**, **ymin**, **ymax** (*float*) -- the left, right, top, and bottom indices for the rough border.
    """
    img_lowres = resample.downsample(img, fact=10)

    rowmean = img_lowres.mean(axis=0)
    smooth_row = _moving_average(rowmean, n=5)

    colmean = img_lowres.mean(axis=1)
    smooth_col = _moving_average(colmean, n=5)

    # xmin = 10 * np.where(rowmean > np.percentile(rowmean, 10))[0][0]
    # xmin = 10 * (np.argmax(np.diff(smooth_row)) + 1)

    # if xmin / img.shape[1] < 0.001:
    #     xmin = np.nan

    # xmax = 10 * np.where(rowmean > np.percentile(rowmean, 10))[0][-1]
    # xmax = 10 * (np.argmin(np.diff(smooth_row)) + 1)

    # if xmax / img.shape[1] > 0.999:
    #     xmax = np.nan
    # sorted_row = np.argsort(np.diff(smooth_row))

    # get the location in the sorted array that corresponds to the minimum and maximum
    # of the difference, that's also in the right half of the image
    # min_ind = np.where(sorted_row < 0.2 * sorted_row.size)[0][-1]
    # max_ind = np.where(sorted_row > 0.8 * sorted_row.size)[0][0]
    col_peaks = peak_local_max(np.diff(smooth_col), min_distance=20, threshold_rel=0.1, num_peaks=2).flatten()
    col_troughs = peak_local_max(-np.diff(smooth_col), min_distance=20, threshold_rel=0.1, num_peaks=2).flatten()

    row_peaks = peak_local_max(np.diff(smooth_row), min_distance=20, threshold_rel=0.1, num_peaks=2).flatten()
    row_troughs = peak_local_max(-np.diff(smooth_row), min_distance=20, threshold_rel=0.1, num_peaks=2).flatten()

    left_ind = np.max(row_peaks[np.where(row_peaks < 0.2 * rowmean.size)[0]], initial=-1e10)
    right_ind = np.min(row_troughs[np.where(row_troughs > 0.8 * rowmean.size)[0]], initial=1e10)

    top_ind = np.max(col_peaks[np.where(col_peaks < 0.2 * colmean.size)[0]], initial=-1e10)
    bot_ind = np.min(col_troughs[np.where(col_troughs > 0.8 * colmean.size)[0]], initial=1e10)

    # xmin = 10 * (sorted_row[min_ind] + 1)
    # xmax = 10 * (sorted_row[max_ind] + 1)
    xmin = 10 * (left_ind + 1)
    xmax = 10 * (right_ind + 1)

    # ymin = 10 * np.where(colmean > np.percentile(colmean, 10))[0][0]
    # ymax = 10 * np.where(colmean > np.percentile(colmean, 10))[0][-1]
    # sorted_col = np.argsort(np.diff(smooth_col))

    # get the location in the sorted array that corresponds to the minimum and maximum
    # of the difference, that's also in the right half of the image
    # min_ind = np.where(sorted_col < 0.2 * sorted_col.size)[0][-1]
    # max_ind = np.where(sorted_col > 0.8 * sorted_col.size)[0][0]
    ymin = 10 * (top_ind + 1)
    ymax = 10 * (bot_ind + 1)

    # ymin = 10 * (sorted_col[min_ind] + 1)
    # ymax = 10 * (sorted_col[max_ind] + 1)

    return xmin, xmax, ymin, ymax


def get_parts_list(im_pattern):
    """
    Find all of the parts of a scanned image that match a given filename pattern

    :param str im_pattern: the image pattern to match
    :return: **parts_list** (*list*) -- a list of all parts of the image that match the pattern.
    """
    imlist = glob(im_pattern + '*.tif')
    imlist.sort()

    return [os.path.splitext(fn_img.split('_')[-1])[0] for fn_img in imlist]


def join_hexagon(im_pattern, overlap=2000, block_size=None, blend=True, is_reversed=False):
    """
    Join multiple parts of a scanned image.

    :param str im_pattern: the base name of the image to use (e.g., DZB1216-500280L002001).
    :param int overlap: the overlap, in pixels, between the image parts.
    :param int block_size: the number of rows each sub-block should cover. Defaults to overlap.
    :param bool blend: apply a linear blend between the two scanned halves (default: True).
    :param bool is_reversed: parts are in reversed order (i.e., part b is the left part, part a is the right part)
    """
    parts = get_parts_list(im_pattern)

    if is_reversed:
        parts.reverse()

    left = io.imread('{}_{}.tif'.format(im_pattern, parts[0]))
    for part in parts[1:]:
        right = io.imread('{}_{}.tif'.format(im_pattern, part))

        left = join_halves(left, right, overlap, block_size=block_size, blend=blend)
        if len(parts) > 2:
            io.imsave('tmp_left.tif', left.astype(np.uint8))

    io.imsave('{}.tif'.format(im_pattern), left.astype(np.uint8))


def _blend(_left, _right, left_shape):
    first = np.where(np.sum(_right, axis=0) > 0)[0][0]
    last = left_shape[1]

    m = 1 / (first - last)
    alpha = np.ones(_left.shape, dtype=np.float32)
    alpha[:, last:] = 0
    for ind in np.arange(first, last):
        alpha[:, ind] = 1 + m * (ind - first)

    return alpha * _left + (1 - alpha) * _right


def join_halves(left, right, overlap, block_size=None, blend=True, trim=None):
    """
    Join two halves of a scanned image together.

    :param array-like left: the left half of the image
    :param array-like right: the right half of the image
    :param int overlap: the amount of overlap, in pixels, between the two halves.
    :param int block_size: the number of rows each sub-block should cover. Defaults to overlap.
    :param bool blend: apply a linear blend between the two scanned halves (default: True).
    :param int trim: the amount to trim the right side of the image by. (default: None).
    """
    M, num_inliers = matching.match_halves(left, right, overlap=overlap, block_size=block_size)

    if num_inliers < 10:
        print('Not enough tie points found. Re-trying with a larger overlap.')
        M, num_inliers = matching.match_halves(left, right, overlap=2*overlap, block_size=block_size)

        if num_inliers < 10:
            raise RuntimeError("Unable to find a reliable transformation between left and right halves.")

    out_shape = (left.shape[0], left.shape[1] + right.shape[1])

    combined_right = warp(right, M, output_shape=out_shape, preserve_range=True, order=3)

    combined_left = np.zeros(out_shape, dtype=np.uint8)
    combined_left[:, :left.shape[1]] = left

    if blend:
        combined = _blend(combined_left, combined_right, left.shape)
    else:
        combined_right[:, :left.shape[1]] = 0
        combined = combined_left + combined_right

    # last_ind = np.where(np.sum(combined, axis=0) > 0)[0][-1]
    right_edge = M.inverse(np.array([[right.shape[1], 0], [right.shape[1], right.shape[0]]]))
    last_ind = int(min(right_edge[:, 0]))

    if trim is not None:  # there has to be a way to get the "correct" end here
        return combined[:, :last_ind - int(trim)]
    else:
        return combined[:, :last_ind]

