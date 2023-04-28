"""
spymicmac.resample is a collection of tools for resampling images
"""
import os
import PIL.Image
from osgeo import gdal
from skimage import io
from scipy import ndimage
import numpy as np
from spymicmac import image, micmac


def downsample(img, fact=4):
    """
    Rescale an image using Lanczos resampling

    :param array-like img: the image to rescale
    :param numeric fact: the number by which to divide the image width and height (default: 4)
    :return: **rescaled** (*array-like*) -- the rescaled image
    """
    _img = PIL.Image.fromarray(img)
    return np.array(_img.resize((np.array(_img.size) / fact).astype(int), PIL.Image.Resampling.LANCZOS))


def resample_hex(fn_img, scale, ori='InterneScan'):
    """
    Resample a KH-9 Mapping Camera image based on the reseau grid, using gdal.Warp

    :param str fn_img: the filename of the image to resample
    :param int scale: the number of pixels per mm of the scanned image
    :param str ori: the Ori directory that contains both MeasuresCamera.xml and MeasuresIm (default: InterneScan)
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

    out_ds = gdal.Warp('tmp_{}'.format(fn_img), fn_img, xRes=1, yRes=1,
                       outputBounds=[0, 0, all_meas.j_cam.max(), all_meas.i_cam.max()],
                       resampleAlg=gdal.GRA_Lanczos)
    meta_shp = '{"shape": ' + '[{}, {}]'.format(out_ds.RasterYSize, out_ds.RasterXSize) + '}'
    out_ds.SetMetadata({'TIFFTAG_IMAGEDESCRIPTION': meta_shp})
    out_ds.FlushCache()
    out_ds = None

    img = io.imread('tmp_{}'.format(fn_img))
    io.imsave('OIS-Reech_{}'.format(fn_img), np.flipud(img).astype(np.uint8))

    os.remove('tmp_{}'.format(fn_img))
    os.remove('{}.aux.xml'.format(fn_img))


def rotate_kh4(img):
    """
    Use the rail marks in a KH-4 image to rotate the image.

    :param array-like img: the image to rotate.
    :return: **rotated** (*array-like*) -- the rotated image
    """
    rails = image.find_rail_marks(img)
    slope, intercept = np.polyfit(rails[:, 1], rails[:, 0], 1)
    angle = np.rad2deg(np.arctan(slope))
    print('Calculated angle of rotation: {:.4f}'.format(angle))

    return ndimage.rotate(img, angle)


def resample_kh4(img):
    """
    **INCOMPLETE**

    :param array-like img: the image to resample
    :return:
    """
    rotated = rotate_kh4(img)
    left, right, top, bot = image.get_rough_frame(rotated)
    rails = image.find_rail_marks(rotated)

