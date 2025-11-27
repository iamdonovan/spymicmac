"""
spymicmac.resample is a collection of tools for resampling images
"""
import os
from pathlib import Path
import multiprocessing as mp
import PIL.Image
import pandas as pd
from osgeo import gdal
from skimage import io, filters
from skimage.measure import ransac, LineModelND
from skimage.transform import AffineTransform, warp
from shapely.geometry import Point, Polygon, LineString, MultiPoint
from scipy import ndimage
from scipy.stats import rankdata
import numpy as np
from . import image, micmac, matching
from numpy.typing import NDArray
from typing import Union


def downsample(img: NDArray, fact: Union[int, float] = 4) -> NDArray:
    """
    Rescale an image using Lanczos resampling

    :param img: the image to rescale
    :param fact: the number by which to divide the image width and height
    :return: **rescaled** -- the rescaled image
    """
    _img = PIL.Image.fromarray(img)
    return np.array(_img.resize((np.array(_img.size) / fact).astype(int), PIL.Image.Resampling.LANCZOS))


def resample_hex(fn_img: Union[str, Path], scale: int, ori: str = 'InterneScan',
                 alg=gdal.GRA_Bilinear, tps: bool = True, order: Union[int, None] = None) -> None:
    """
    Resample a KH-9 Mapping Camera image based on the reseau grid, using gdal.Warp

    :param fn_img: the filename of the image to resample
    :param scale: the number of pixels per mm of the scanned image
    :param ori: the Ori directory that contains both MeasuresCamera.xml and MeasuresIm
    :param alg: the gdal resampling algorithm to use (default: gdal.GRA_Bilinear)
    :param tps: use a thin plate spline transformer to transform based on reseau grid
    :param order: the order (1-3) of polynomial GCP interpolation (default: not used)
    """
    cam_meas = micmac.parse_im_meas(Path(f"Ori-{ori}", 'MeasuresCamera.xml'))
    img_meas = micmac.parse_im_meas(Path(f"Ori-{ori}", f"MeasuresIm-{fn_img}.xml"))

    all_meas = img_meas.set_index('name').join(cam_meas.set_index('name'), lsuffix='_img', rsuffix='_cam')
    all_meas['i_cam'] *= scale
    all_meas['j_cam'] *= scale

    gcp_list = [gdal.GCP(row.j_cam, row.i_cam, 1, row.j_img, row.i_img) for _, row in all_meas.iterrows()]

    ds = gdal.Open(fn_img)
    ds.SetGCPs(gcp_list, '')
    ds.FlushCache()
    ds = None

    options = {'xRes': 1, 'yRes': 1,
               'outputBounds': [0, 0, all_meas.j_cam.max(), all_meas.i_cam.max()],
               'resampleAlg': alg,
               'tps': tps}

    if order is not None:
        options['polynomialOrder'] = order

    out_ds = gdal.Warp(f"tmp_{fn_img}", fn_img, **options)

    meta_shp = '{"shape": ' + f"[{out_ds.RasterYSize}, {out_ds.RasterXSize}]" + '}'
    out_ds.SetMetadata({'TIFFTAG_IMAGEDESCRIPTION': meta_shp})
    out_ds.FlushCache()
    out_ds = None

    img = io.imread(f"tmp_{fn_img}")
    io.imsave(f"OIS-Reech_{fn_img}", np.flipud(img).astype(np.uint8))

    os.remove(f"tmp_{fn_img}")
    os.remove(f"{fn_img}.aux.xml")


def rotate_from_rails(img: NDArray, rails: NDArray) -> tuple[NDArray, float]:
    """
    Use the rail marks or other horizontal points in an image to rotate the image.

    :param img: the image to rotate.
    :param rails: an Nx2 array of (row, col) points
    :return:
        - **rotated** -- the rotated image
        - **angle** -- the calculated angle of rotation, in degrees
    """
    slope, intercept = np.polyfit(rails[:, 1], rails[:, 0], 1)
    angle = np.rad2deg(np.arctan(slope))
    print(f"Calculated angle of rotation: {angle:.4f}")

    return ndimage.rotate(img, angle), angle


def crop_panoramic(fn_img: Union[str, Path], fact: Union[int, None] = None,
                   fn_coords: Union[None, str, Path] = None) -> None:
    """
    Crop a declassified panoramic (KH4 or KH9) image by finding the frame.

    :param fn_img: the filename of the image to rotate and crop
    :param fact: the number by which to divide the image width and height to scale the image (default: do not scale)
    :param fn_coords:
    """
    img = io.imread(fn_img)
    img_lowres = downsample(img, fact=10)

    if fn_coords is not None:
        box = _load_box(fn_coords)
    else:
        box = _find_rectangle(img_lowres)

    cropped, coords = box_transform(img, box)

    src, dst = coords
    coords = pd.DataFrame(data={'orig_x': src[:, 0], 'orig_y': src[:, 1],
                                'new_x': dst[:, 0], 'new_y': dst[:, 1]})
    coords.to_csv(f"{os.path.splitext(fn_img)[0]}_coords.csv", index=False)

    if fact is not None:
        cropped = downsample(cropped, fact=fact)

    # if the image is from the aft camera, flip it before saving
    # USGS naming convention is *[A|F]00N, so to figure out if it's an aft image,
    # remove the extension, check the character 4 places from the end
    is_aft = os.path.splitext(fn_img)[0][-4] == 'A'
    if is_aft:
        cropped = np.fliplr(np.flipud(cropped))

    # save the resampled image
    io.imsave('OIS-Reech_' + fn_img, cropped.astype(np.uint8))


def _load_box(fn_coords):
    coords = pd.read_csv(fn_coords)
    return Polygon(zip(coords.orig_x, coords.orig_y))


def _find_line(img, axis):
    peaks = np.argmax(img, axis=axis)

    if axis == 0:
        xs = np.arange(0, img.shape[1])
    else:
        xs = np.arange(0, img.shape[0])

    model, inliers = ransac(np.column_stack([xs, peaks]), LineModelND, min_samples=2, residual_threshold=5)
    p = np.polyfit(xs[inliers], peaks[inliers], 1)

    return p, xs[inliers], peaks[inliers]


def _find_rectangle(img):
    filtered = filters.sobel(filters.gaussian(img, 4))

    cols = np.arange(0, filtered.shape[1], 100)
    rows = np.arange(0, filtered.shape[0], 100)

    row_width = int(0.3 * filtered.shape[0])
    col_width = int(0.1 * filtered.shape[1])

    p, xtop, ptop = _find_line(filtered[:row_width, :], 0)
    top_line = LineString(list(zip(cols, np.polyval(p, cols))))

    p, xbot, pbot = _find_line(filtered[-row_width:-100, :], 0)
    bot_line = LineString(list(zip(cols, np.polyval(p, cols) + (filtered.shape[0] - row_width))))

    p, xleft, pleft = _find_line(filtered[:, :col_width], 1)
    left_line = LineString(list(zip(np.polyval(p, rows), rows)))

    p, xright, pright = _find_line(filtered[:, -col_width:], 1)
    right_line = LineString(list(zip(np.polyval(p, rows) + (filtered.shape[1] - col_width), rows)))

    ul = top_line.intersection(left_line)
    ur = top_line.intersection(right_line)
    ll = bot_line.intersection(left_line)
    lr = bot_line.intersection(right_line)

    x_all = np.concatenate([xtop, xbot, pleft, pright + (filtered.shape[1] - col_width)])
    p_all = np.concatenate([ptop, pbot + (filtered.shape[0] - row_width), xleft, xright])

    points = np.array([Point(x, y) for x, y in zip(x_all, p_all)])
    init_box = Polygon([ul, ur, lr, ll]).buffer(10, cap_style='square', join_style='mitre')

    points = points[[pt.within(init_box) for pt in points]]

    box = MultiPoint(points).minimum_rotated_rectangle.buffer(-13, join_style='mitre', cap_style='square')

    bx, by = box.exterior.coords.xy

    bx = np.array(bx) * 10
    by = np.array(by) * 10

    return Polygon(zip(bx, by))


def box_transform(img: NDArray, rect: Polygon):
    """
    Crop an image to a rectangle while rotating so the image is horizontal, using an affine transformation.

    :param img: the image to crop + rotate
    :param rect: the rectangle representing the bounds of the new image
    """
    # get width/height of the rectangle
    x, y = rect.exterior.coords.xy
    edges = [Point(x[ind], y[ind]).distance(Point(x[ind + 1], y[ind + 1])) for ind in range(4)]

    height = int(min(edges))
    width = int(max(edges))

    # rank the x, y values to find the lowest combined rank (should be UL corner)
    # if the rectangle is given by minimum_rotated_rectangle, it is oriented - so we should be able
    # to cycle around the corners counter-clockwise, starting from the upper left corner
    comb_ranks = rankdata(x[:-1]) + rankdata(y[:-1])
    ul_ind = np.argmin(comb_ranks)
    ll_ind = (ul_ind + 1) % 4
    lr_ind = (ll_ind + 1) % 4
    ur_ind = (lr_ind + 1) % 4

    ordered_x = [x[ul_ind], x[ll_ind], x[lr_ind], x[ur_ind]]
    ordered_y = [y[ul_ind], y[ll_ind], y[lr_ind], y[ur_ind]]

    new_x = [0, 0, width, width]
    new_y = [0, height, height, 0]

    src = np.column_stack([ordered_x, ordered_y])
    dst = np.column_stack([new_x, new_y])

    tfm = AffineTransform()
    tfm.estimate(src, dst)

    return warp(img, tfm.inverse, output_shape=(height, width), preserve_range=True, order=3), (src, dst)


def crop_from_extent(fn_img: Union[str, Path], border: NDArray,
                     angle: Union[float, None] = None, fact: Union[int, None] = None) -> None:
    """
    Crop an image given the coordinates of the image border.

    :param fn_img: the filename of the image to rotate and crop
    :param border: the estimated image border coordinates (left, right, top, bot)
    :param angle: the angle by which to rotate the image (default: None)
    :param fact: the number by which to divide the image width and height to scale the image (default: do not scale)
    """
    img = io.imread(fn_img)

    if angle is not None:
        img = ndimage.rotate(img, angle)

    left, right, top, bot = border

    print(f'Cropping to window (left, right, top, bot): {int(left)}, {int(right)}, {int(top)}, {int(bot)}')
    cropped = img[int(top):int(bot), int(left):int(right)]

    if fact is not None:
        cropped = downsample(cropped, fact=fact)

    # if the image is from the aft camera, flip it before saving
    # USGS naming convention is *[A|F]00N, so to figure out if it's an aft image,
    # remove the extension, check the character 4 places from the end
    is_aft = os.path.splitext(fn_img)[0][-4] == 'A'
    if is_aft:
        cropped = np.fliplr(np.flipud(cropped))

    io.imsave('OIS-Reech_' + fn_img, cropped.astype(np.uint8))


def _border_mask(img, border):
    mask = 255 * np.ones(img.shape, dtype=np.uint8)
    mask[border:-border, border:-border] = 0
    return mask


def align_image_borders(fn_left, fn_right, border):

    left = io.imread(fn_left)
    right = io.imread(fn_right)

    mask_left = _border_mask(left, border)
    mask_right = _border_mask(right, border)


def resample_fiducials(fn_img: Union[str, Path, list[str], list[Path]],
                       scale: float, transform=AffineTransform(),
                       fn_cam: Union[str, Path, None] = None, nproc: int = 1) -> None:
    """
    Resample image(s) using fiducial markers.

    :param fn_img: the filename, or a list of filenames, of the image(s)
    :param scale: the image scale (in mm/pixel) to use
    :param transform: the type of transformation to use. Should be an instance of skimage.transform (default: AffineTransform)
    :param fn_cam: the filename for the MeasuresCamera.xml file (default: Ori-InterneScan/MeasuresCamera.xml)
    :param nproc: the number of processors to use (default: 1)
    """
    if type(fn_img) is str:
        _fiducials(fn_img=fn_img, scale=scale, fn_cam=fn_cam, transform=transform)
    else:
        if nproc > 1:
            pool = mp.Pool(nproc, maxtasksperchild=1)
            arg_dict = {'scale': scale, 'fn_cam': fn_cam, 'transform': transform}
            pool_args = [{'fn_img': fn} for fn in fn_img]
            for d in pool_args:
                d.update(arg_dict)

            pool.map(_fiducials_wrapper, pool_args, chunksize=1)
            pool.close()
            pool.join()
        else:
            for fn in fn_img:
                _fiducials(fn_img=fn, scale=scale, fn_cam=fn_cam, transform=transform)


def _fiducials_wrapper(args):
    _fiducials(**args)


def _fiducials(fn_img=None, scale=None, fn_cam=None, transform=AffineTransform):
    print(fn_img)
    meas = micmac.parse_im_meas(Path('Ori-InterneScan', f'MeasuresIm-{fn_img}.xml'))
    if fn_cam is None:
        measures_cam = micmac.parse_im_meas(Path('Ori-InterneScan', 'MeasuresCamera.xml'))
    else:
        measures_cam = micmac.parse_im_meas(fn_cam)

    measures_cam['i'] /= scale
    measures_cam['j'] /= scale

    joined = meas.set_index('name').join(measures_cam.set_index('name'), lsuffix='_img', rsuffix='_cam')

    transform.estimate(joined[['j_img', 'i_img']].values,
                       joined[['j_cam', 'i_cam']].values)

    outshape = (int(joined['i_cam'].max() - joined['i_cam'].min()),
                int(joined['j_cam'].max() - joined['j_cam'].min()))

    img = io.imread(fn_img)
    img_tfm = warp(img, transform.inverse, output_shape=outshape, preserve_range=True, order=5)

    io.imsave('OIS-Reech_' + fn_img, img_tfm.astype(np.uint8))
