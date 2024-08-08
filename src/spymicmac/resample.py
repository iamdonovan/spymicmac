"""
spymicmac.resample is a collection of tools for resampling images
"""
import os
import multiprocessing as mp
import PIL.Image
from osgeo import gdal
from skimage import io
from skimage.transform import AffineTransform, SimilarityTransform, warp
from scipy import ndimage
import numpy as np
from spymicmac import image, micmac, matching


def downsample(img, fact=4):
    """
    Rescale an image using Lanczos resampling

    :param array-like img: the image to rescale
    :param numeric fact: the number by which to divide the image width and height (default: 4)
    :return: **rescaled** (*array-like*) -- the rescaled image
    """
    _img = PIL.Image.fromarray(img)
    return np.array(_img.resize((np.array(_img.size) / fact).astype(int), PIL.Image.Resampling.LANCZOS))


def resample_hex(fn_img, scale, ori='InterneScan', alg=gdal.GRA_Bilinear, tps=True, order=None):
    """
    Resample a KH-9 Mapping Camera image based on the reseau grid, using gdal.Warp

    :param str fn_img: the filename of the image to resample
    :param int scale: the number of pixels per mm of the scanned image
    :param str ori: the Ori directory that contains both MeasuresCamera.xml and MeasuresIm (default: InterneScan)
    :param alg: the gdal resampling algorithm to use (default: gdal.GRA_Bilinear)
    :param bool tps: use a thin plate spline transformer to transform based on reseau grid (default: False)
    :param int order: the order (1-3) of polynomial GCP interpolation (default: not used)
    """
    cam_meas = micmac.parse_im_meas(os.path.join('Ori-{}'.format(ori), 'MeasuresCamera.xml'))
    img_meas = micmac.parse_im_meas(os.path.join('Ori-{}'.format(ori), 'MeasuresIm-{}.xml'.format(fn_img)))

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

    out_ds = gdal.Warp('tmp_{}'.format(fn_img), fn_img, **options)

    meta_shp = '{"shape": ' + '[{}, {}]'.format(out_ds.RasterYSize, out_ds.RasterXSize) + '}'
    out_ds.SetMetadata({'TIFFTAG_IMAGEDESCRIPTION': meta_shp})
    out_ds.FlushCache()
    out_ds = None

    img = io.imread('tmp_{}'.format(fn_img))
    io.imsave('OIS-Reech_{}'.format(fn_img), np.flipud(img).astype(np.uint8))

    os.remove('tmp_{}'.format(fn_img))
    os.remove('{}.aux.xml'.format(fn_img))


def rotate_from_rails(img, rails):
    """
    Use the rail marks or other horizontal points in an image to rotate the image.

    :param array-like img: the image to rotate.
    :param array-like rails: an Nx2 array of (row, col) points
    :return:
        - **rotated** (*array-like*) -- the rotated image
        - **angle** (*float*) -- the calculated angle of rotation, in degrees
    """
    slope, intercept = np.polyfit(rails[:, 1], rails[:, 0], 1)
    angle = np.rad2deg(np.arctan(slope))
    print('Calculated angle of rotation: {:.4f}'.format(angle))

    return ndimage.rotate(img, angle), angle


def crop_panoramic(fn_img, flavor, marker_size=31, fact=None, return_vals=False):
    """
    Crop a declassified panoramic (KH4 or KH9) image, after rotating based on horizontal rail markers or "wagon wheel"
    fiducial markers.

    :param str fn_img: the filename of the image to rotate and crop
    :param str flavor: the camera type (KH4 or KH9)
    :param int marker_size: The approximate size of the wagon wheels to identify in the image (default: 31 pixels)
    :param int fact: the number by which to divide the image width and height to scale the image (default: do not scale)
    :param bool return_vals: Return estimated image border and rotation angle (default: False)
    :returns:
        - **border** (*tuple*) -- the estimated image border (left, right, top, bot)
        - **angle**  (*float*) -- the estimated rotation angle.
        Only returned if **return_vals** is True.
    """
    assert flavor in ['KH4', 'KH9'], "flavor must be one of [KH4, KH9]"

    img = io.imread(fn_img)
    if flavor == 'KH4':
        rails = matching.find_rail_marks(img)
        rotated, angle = rotate_from_rails(img, rails)

        # crop the lower part of the image to avoid introducing a bright line at the bottom
        rotated = rotated[:-int(0.1 * rails[:, 0].mean()), :]
    else:
        rails = matching.ocm_show_wagon_wheels(img, size=marker_size)

        # restrict ourselves to the top rail
        rails = rails[rails[:, 0] < 0.1 * img.shape[0]]

        # refine the choice to ensure the points are on the same line
        valid = matching._refine_rail(rails)

        rotated, angle = rotate_from_rails(img, rails[valid])

    # get a rough idea of where the image frame should be
    left, right, top, bot = image.get_rough_frame(rotated)

    # buffer by 0.05% (left/right) or 0.5% (top/bottom) to ensure we have a clean image
    left += (right - left) * 0.0005
    right -= (right - left) * 0.0005

    top += (bot - top) * 0.005
    bot -= (bot - top) * 0.005

    # crop the image to the buffered window
    print(f'Cropping to window (left, right, top, bot): {int(left)}, {int(right)}, {int(top)}, {int(bot)}')
    cropped = rotated[int(top):int(bot), int(left):int(right)]

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

    if return_vals:
        return (left, right, top, bot), angle


def crop_from_extent(fn_img, border, angle=None, fact=None):
    """
    Crop an image given the coordinates of the image border.

    :param str fn_img: the filename of the image to rotate and crop
    :param array-like border: the estimated image border coordinates (left, right, top, bot)
    :param float angle: the angle by which to rotate the image (default: None)
    :param int fact: the number by which to divide the image width and height to scale the image (default: do not scale)
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


def resample_fiducials(fn_img, scale, transform=AffineTransform(), fn_cam=None, nproc=1):
    """
    Resample image(s) using fiducial markers.

    :param str|list fn_img: the filename, or a list of filenames, of the image(s)
    :param float scale: the image scale (in mm/pixel) to use
    :param transform: the type of transformation to use. Should be an instance of skimage.transform (default: AffineTransform)
    :param str fn_cam: the filename for the MeasuresCamera.xml file (default: Ori-InterneScan/MeasuresCamera.xml)
    :param int nproc: the number of processors to use (default: 1)
    :return:
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
    meas = micmac.parse_im_meas(os.path.join('Ori-InterneScan', f'MeasuresIm-{fn_img}.xml'))
    if fn_cam is None:
        measures_cam = micmac.parse_im_meas(os.path.join('Ori-InterneScan', 'MeasuresCamera.xml'))
    else:
        measures_cam = micmac.parse_im_meas(fn_cam)

    measures_cam['i'] /= scale
    measures_cam['j'] /= scale

    joined = meas.set_index('name').join(measures_cam.set_index('name'), lsuffix='_img', rsuffix='_cam')

    transform.estimate(joined[['j_img', 'i_img']].values,
                       joined[['j_cam', 'i_cam']].values)

    outshape = ((joined['i_cam'].max() - joined['i_cam'].min()),
                (joined['j_cam'].max() - joined['j_cam'].min()))

    img = io.imread(fn_img)
    img_tfm = warp(img, transform.inverse, output_shape=outshape, preserve_range=True, order=5)

    io.imsave('OIS-Reech_' + fn_img, img_tfm.astype(np.uint8))
