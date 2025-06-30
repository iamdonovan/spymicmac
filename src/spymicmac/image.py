"""
spymicmac.image is a collection of tools for working with images.
"""
import os
from glob import glob
from itertools import chain, product
from skimage import exposure, morphology, io, filters
from skimage.morphology import disk
from skimage.feature import peak_local_max
from skimage.transform import warp
from scipy import ndimage
from scipy.ndimage.filters import generic_filter
import numpy as np
from numba import jit
from . import matching, resample
from numpy.typing import NDArray, DTypeLike
from typing import Union


np.seterr(divide='ignore', invalid='ignore')
######################################################################################################################
# image filtering tools
######################################################################################################################
@jit(nopython=True)
def nanstd(a: NDArray) -> NDArray:
    return np.nanstd(a)


def nanmedian_filter(img: NDArray, **kwargs: tuple) -> NDArray:
    """
    Calculate a multi-dimensional median filter that respects NaN values
    and masked arrays.

    :param img: image on which to calculate the median filter
    :param kwargs: additional arguments to ndimage.generic_filter
        Note that either size or footprint must be defined. size gives the shape
        that is taken from the input array, at every element position, to define
        the input to the filter function. footprint is a boolean array that
        specifies (implicitly) a shape, but also which of the elements within
        this shape will get passed to the filter function. Thus size=(n,m) is
        equivalent to footprint=np.ones((n,m)). We adjust size to the number
        of dimensions of the input array, so that, if the input array is
        shape (10,10,10), and size is 2, then the actual size used is (2,2,2).
    :type img: array-like

    :returns filtered: Filtered array of same shape as input.
    """
    # set up the wrapper function to call generic filter
    @jit(nopython=True)
    def nanmed(a):
        return np.nanmedian(a)

    return generic_filter(img, nanmed, **kwargs)


def highpass_filter(img: NDArray) -> NDArray:
    """
    Subtract a low-pass from an image, to return a highpass filter.

    :param img: the image to filter.
    :return: the highpass-filtered image
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


def splitter(img: NDArray, nblocks: tuple[int, int], overlap: int = 0) -> tuple[list, list, list]:
    """
    Split an image into (m, n) blocks with a given overlap.

    :param img: the image to split
    :param nblocks: the number of blocks to create along each axis (m, n)
    :param overlap: the number of pixels to overlap each block.

    :return:
        - **blocks** -- a list of the image blocks created
        - **top_inds** -- a list of the original row index of the top edge of each block
        - **left_inds** -- a list of the original column index of the left edge of each block
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
            lind = max(0, j*new_width - overlap - 1)
            rind = min(img.shape[1], (j+1)*new_width + overlap)
            tind = max(0, i*new_height - overlap - 1)
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


def stretch_image(img: NDArray, scale: tuple[float, float] = (0.0, 1.0),
                  mult: Union[int, float] = 255, imgmin: Union[int, float] = 0,
                  outtype: DTypeLike = np.uint8, mask: Union[NDArray, None] = None) -> NDArray:
    """
    Apply a linear stretch to an image by clipping and stretching to quantiles.

    :param img: the image to stretch.
    :param scale: the minimum and maximum quantile to stretch to. Defaults to the minimum/maximum values
        of the image.
    :param mult: a multiplier to scale the result to.
    :param imgmin: the minimum value in the output image
    :param outtype: the numpy datatype to return the stretched image as.
    :param mask: a mask of pixels to ignore when calculating quantiles.
    :return: the stretched image.
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


def remove_scanner_stripes(img: NDArray, dtype: DTypeLike = np.uint8, scan_axis: int = 1) -> NDArray:
    """
    Remove horizontal (or vertical) stripes from an image.

    :param img: the image to remove stripes from.
    :param dtype: the original datatype of the image.
    :param scan_axis: the axis corresponding to the direction of the stripes. A scan_axis of 1 corresponds to
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


def contrast_enhance(fn_img: str, mask_value: Union[int, float, None] = None,
                     qmin: float = 0.02, qmax: float = 0.98, gamma: float = 1.25,
                     disksize: int = 3, imgmin: Union[int, float] = 0) -> NDArray:
    """
    Enhance image contrast in a three-step process. First, the image is processed with a median filter to reduce
    noise. Next, a linear contrast stretch is applied, and finally, a gamma adjustment is applied.

    :param fn_img: the image filename.
    :param mask_value: a mask value to use when filtering the image.
    :param qmin: the minimum quantile to use for the linear contrast stretch
    :param qmax: the maximum quantile to use for the linear contrast stretch
    :param gamma: the value to use for the gamma adjustment
    :param disksize: the filter disk size (input to skimage.morphology.disk)
    :param imgmin: the minimum value in the output image
    :return: **enhanced** -- the contrast-enhanced image.
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


def high_low_subtract(img: NDArray) -> NDArray:
    """
    Remove the column and row median values from an image, then stretch/scale to the new min/max values.

    :param img: the image to adjust
    :returns: the adjusted image
    """

    row_med = np.median(img, axis=1)
    col_med = np.median(img, axis=0)

    subtract = img - col_med * np.ones_like(img)
    subtract -= row_med.reshape(-1, 1) * np.ones_like(img)

    return stretch_image(subtract)


def make_binary_mask(img: NDArray, mult_value: Union[int, float] = 255,
                     erode: int = 0, mask_value: int = 0) -> NDArray:
    """
    Create a binary mask for an image based on a given mask value. Values equal to mask_value will be given a value
    of 0 in the mask, while all other values will be set equal to mult_value.

    :param img: the image to create a mask for.
    :param mult_value: the value indicating a non-masked value
    :param erode: the size of the erosion operation to apply
    :param mask_value: the value to mask in the image
    :return: the binary mask.
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


def balance_image(img: NDArray, clip_limit: float = 0.01) -> NDArray:
    """
    Apply contrast-limited adaptive histogram equalization (CLAHE) on an image.

    :param img: the image to balance.
    :param clip_limit: Clipping limit, normalized between 0 and 1 (higher values give
        more contrast).
    :return: **img_filt** -- the balanced, filtered image.
    """
    return (255 * exposure.equalize_adapthist(img, clip_limit=clip_limit)).astype(np.uint8)


# thanks to SO user Jamie for this answer
# https://stackoverflow.com/a/14314054
def _moving_average(a: NDArray, n: int = 5) -> NDArray:
    ret = np.cumsum(a)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


# because sometimes, usgs adds a border with a bright white line on it
def _spike_filter(img: NDArray, axis: int) -> NDArray:
    img_mean = img.mean(axis=axis)

    otsu = img > filters.threshold_otsu(img)
    pct_bright = np.count_nonzero(otsu, axis=axis) / otsu.shape[axis]
    z = (pct_bright - pct_bright.mean()) / pct_bright.std()

    inds, = np.where(z > 3)
    if inds.size > 0:
        window = []
        for ind in inds:
            window.extend(list(range(ind-10, ind+11)))
        _inds = np.array(list(set(window)))
        _inds = _inds[np.logical_and(_inds >= 0, _inds < img_mean.size)]
        img_mean[_inds] = np.mean([img_mean[min(_inds)], img_mean[max(_inds)]])

    return img_mean


def get_rough_frame(img: NDArray, fact: int = 10) -> tuple[float, float, float, float]:
    """
    Find the rough location of an image frame/border.

    :param img: the image to find a border for
    :param fact: the scaling factor for the low-resolution image
    :return: **xmin**, **xmax**, **ymin**, **ymax** -- the left, right, top, and bottom indices for the rough border.
    """
    img_lowres = filters.gaussian(resample.downsample(img, fact=fact), 4)
    aspect = min(img.shape) / max(img.shape)
    if aspect > 0.75:
        lower, upper = 0.1, 0.9
    else:
        lower, upper = 0.15, 0.85

    rowmean = _spike_filter(img_lowres, axis=0)
    smooth_row = _moving_average(rowmean, n=20)

    colmean = _spike_filter(img_lowres, axis=1)
    smooth_col = _moving_average(colmean, n=20)

    lr = int(lower * smooth_row.size)
    ur = int(upper * smooth_row.size)

    lc = int(lower * smooth_col.size)
    uc = int(upper * smooth_col.size)
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
    col_peaks = peak_local_max(np.diff(smooth_col[:lc]), min_distance=int(0.001 * smooth_col.size),
                               threshold_rel=0.25).flatten()
    col_troughs = peak_local_max(-np.diff(smooth_col[uc:]), min_distance=int(0.001 * smooth_col.size),
                                 threshold_rel=0.25).flatten() + uc

    row_peaks = peak_local_max(np.diff(smooth_row[:lr]), min_distance=int(0.001 * smooth_row.size),
                               threshold_rel=0.25).flatten()
    row_troughs = peak_local_max(-np.diff(smooth_row[ur:]), min_distance=int(0.001 * smooth_row.size),
                                 threshold_rel=0.25).flatten() + ur

    # left_ind = np.max(row_peaks[np.where(row_peaks < lower * rowmean.size)[0]], initial=-1e10)
    # right_ind = np.min(row_troughs[np.where(row_troughs > upper * rowmean.size)[0]], initial=1e10)
    # left_ind, right_ind = _maximum_sep(row_peaks, row_troughs)
    left_ind = _best_sep(row_peaks, smooth_row, False)
    right_ind = _best_sep(row_troughs, smooth_row, True)

    if (right_ind - left_ind) / smooth_row.size < (upper - lower):
        left_ind = np.max(row_peaks[np.where(row_peaks < lower * rowmean.size)[0]], initial=-1e10)
        right_ind = np.min(row_troughs[np.where(row_troughs > upper * rowmean.size)[0]], initial=1e10)

    # top_ind = np.max(col_peaks[np.where(col_peaks < lower * colmean.size)[0]], initial=-1e10)
    # bot_ind = np.min(col_troughs[np.where(col_troughs > upper * colmean.size)[0]], initial=1e10)
    # top_ind, bot_ind = _maximum_sep(col_peaks, col_troughs)
    top_ind = _best_sep(col_peaks, smooth_col, False)
    bot_ind = _best_sep(col_troughs, smooth_col, True)

    if (bot_ind - top_ind) / smooth_col.size < (upper - lower):
        top_ind = np.max(col_peaks[np.where(col_peaks < lower * colmean.size)[0]], initial=-1e10)
        bot_ind = np.min(col_troughs[np.where(col_troughs > upper * colmean.size)[0]], initial=1e10)

    # xmin = 10 * (sorted_row[min_ind] + 1)
    # xmax = 10 * (sorted_row[max_ind] + 1)
    xmin = fact * (left_ind + 1)
    xmax = fact * (right_ind + 1)

    # ymin = 10 * np.where(colmean > np.percentile(colmean, 10))[0][0]
    # ymax = 10 * np.where(colmean > np.percentile(colmean, 10))[0][-1]
    # sorted_col = np.argsort(np.diff(smooth_col))

    # get the location in the sorted array that corresponds to the minimum and maximum
    # of the difference, that's also in the right half of the image
    # min_ind = np.where(sorted_col < 0.2 * sorted_col.size)[0][-1]
    # max_ind = np.where(sorted_col > 0.8 * sorted_col.size)[0][0]
    ymin = fact * (top_ind + 1)
    ymax = fact * (bot_ind + 1)

    # ymin = 10 * (sorted_col[min_ind] + 1)
    # ymax = 10 * (sorted_col[max_ind] + 1)

    return xmin, xmax, ymin, ymax


def _best_sep(pks, means, trough):
    davg = min([min(pks), 100, means.size - max(pks)])
    seps = [means[pk-davg:pk].mean() - means[pk:pk+davg].mean() for pk in pks]
    if trough:
        return pks[np.argmax(seps)]
    else:
        return pks[np.argmin(seps)]


def _maximum_sep(peaks, troughs):
    combs = list(product(peaks, troughs))
    dists = [max(pair) - min(pair) for pair in combs]

    best = combs[np.argmax(dists)]

    return min(best), max(best)


def get_parts_list(im_pattern: str) -> list[str]:
    """
    Find all of the parts of a scanned image that match a given filename pattern

    :param im_pattern: the image pattern to match
    :return: **parts_list** -- a list of all parts of the image that match the pattern.
    """
    imlist = glob(im_pattern + '*.tif')
    imlist.sort()

    return [os.path.splitext(fn_img.split('_')[-1])[0] for fn_img in imlist]


def join_hexagon(im_pattern: str, overlap: int = 2000, block_size: Union[int, None] = None,
                 blend: bool = True, is_reversed: bool = False) -> None:
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

    left = io.imread(f"{im_pattern}_{parts[0]}.tif")
    for part in parts[1:]:
        right = io.imread(f"{im_pattern}_{part}.tif")

        left = join_halves(left, right, overlap, block_size=block_size, blend=blend)
        if len(parts) > 2:
            io.imsave('tmp_left.tif', left.astype(np.uint8))

    io.imsave(f"{im_pattern}.tif", left.astype(np.uint8))

    if len(parts) > 2:
        os.remove('tmp_left.tif') # clean up after we're done


def _blend(_left, _right, left_shape):
    first = np.where(np.sum(_right, axis=0) > 0)[0][0]
    last = left_shape[1]

    m = 1 / (first - last)
    alpha = np.ones(_left.shape, dtype=np.float32)
    alpha[:, last:] = 0
    for ind in np.arange(first, last):
        alpha[:, ind] = 1 + m * (ind - first)

    return alpha * _left + (1 - alpha) * _right


def join_halves(left: NDArray, right: NDArray, overlap: int,
                block_size: Union[int, None] = None, blend: bool = True, trim: Union[int, None] = None) -> NDArray:
    """
    Join two halves of a scanned image together.

    :param left: the left half of the image
    :param right: the right half of the image
    :param overlap: the amount of overlap, in pixels, between the two halves.
    :param block_size: the number of rows each sub-block should cover. Defaults to same value as overlap.
    :param blend: apply a linear blend between the two scanned halves.
    :param trim: the amount to trim the right side of the image by. Default is no trimming.
    :return: the joined image.
    """
    M, num_inliers = matching.match_halves(left, right, overlap=overlap, block_size=block_size)

    if num_inliers < 10:
        print('Not enough tie points found. Re-trying with a larger overlap.')
        M, num_inliers = matching.match_halves(left, right, overlap=2*overlap, block_size=block_size)

        if num_inliers < 10:
            raise RuntimeError("Unable to find a reliable transformation between left and right halves.")

    out_shape = (left.shape[0], left.shape[1] + right.shape[1])

    combined_right = warp(right, M, output_shape=out_shape, preserve_range=True, order=3)

    # io.imsave('tmp_right.tif', combined_right.astype(np.uint8))

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

