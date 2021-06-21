"""
spymicmac.register is a collection of tools for registering images and finding GCPs.
"""
import os
import re
import subprocess
import shutil
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
import pandas as pd
import geopandas as gpd
from PIL import Image
from glob import glob
from shapely.geometry.point import Point
from skimage.io import imread
from skimage.measure import ransac
from skimage.filters import median
from skimage.morphology import disk, binary_dilation
from skimage.transform import EuclideanTransform, AffineTransform, warp
from pybob.bob_tools import mkdir_p
from pybob.ddem_tools import nmad
from pybob.image_tools import create_mask_from_shapefile
from pybob.GeoImg import GeoImg
import spymicmac.image as imtools
import spymicmac.micmac as mmtools
from spymicmac.usgs import get_usgs_footprints


def sliding_window_filter(img_shape, pts_df, winsize, stepsize=None, mindist=2000, how='residual', is_ascending=True):
    """
    Given a DataFrame of indices representing points, use a sliding window filter to keep only the 'best' points within
    a given window and separation distance.

    :param array-like img_shape: The shape of the image to filter indices from.
    :param DataFrame pts_df: A DataFrame that contains the pixel locations of points (column: orig_j, row: orig_i)
    :param int winsize: the size of the window to use for the filter.
    :param int stepsize: how large of a step size to use for the sliding window (default: winsize/2)
    :param int mindist: the minimum distance (in pixels) to use between points (default: 2000)
    :param str how: how to sort pts_df to determine the 'best' points to keep (default: 'residual')
    :param bool is_ascending: whether the column named in 'how' should be sorted ascending or descending (default: True)
    :return:
        - **out_inds** (*array-like*) -- an array with the filtered indices
    """
    if stepsize is None:
        stepsize = winsize / 2

    _out_inds = []
    _out_pts = []

    for x_ind in np.arange(stepsize, img_shape[1], winsize):
        for y_ind in np.arange(stepsize, img_shape[0], winsize):
            min_x = x_ind - winsize / 2
            max_x = x_ind + winsize / 2
            min_y = y_ind - winsize / 2
            max_y = y_ind + winsize / 2
            samp_ = pts_df.loc[np.logical_and.reduce([pts_df.orig_i > min_x,
                                                      pts_df.orig_i < max_x,
                                                      pts_df.orig_j > min_y,
                                                      pts_df.orig_j < max_y])].copy()
            if samp_.shape[0] == 0:
                continue
            # only take the above-average z_corr values
            samp_ = samp_[samp_.z_corr >= samp_.z_corr.quantile(0.5)]
            # make sure we get the best residual
            samp_.sort_values(how, ascending=is_ascending, inplace=True)
            if len(_out_inds) == 0:
                best_ind = samp_.index[0]
                best_pt = Point(samp_.loc[best_ind, ['orig_j', 'orig_i']].values)

                _out_inds.append(best_ind)
                _out_pts.append(best_pt)
            else:
                for _ind, _row in samp_.iterrows():
                    this_pt = Point(_row[['orig_j', 'orig_i']].values)
                    this_min_dist = np.array([this_pt.distance(pt) for pt in _out_pts]).min()
                    if this_min_dist > mindist:
                        _out_inds.append(_ind)
                        _out_pts.append(this_pt)

    return np.array(_out_inds)


def get_imlist(im_subset):
    """
    Given either a list of filenames or a regex match pattern, return a list of filenames and a pattern to provide
    to MicMac.

    :param list im_subset: a list of filenames or a match pattern (e.g., ['OIS.*tif']) representing filenames
    :return:
        - **imlist** (*list*) -- the list of filenames matching the provided pattern.
        - **match_pattern** (*str*) -- the match pattern to be provided to MicMac.
    """
    if im_subset is None:
        imlist = glob('OIS*.tif')
        match_pattern = 'OIS.*.tif'
    else:
        if len(im_subset) > 1:
            imlist = im_subset
            # match_pattern = mmtools.get_match_pattern(imlist)
            match_pattern = '|'.join(imlist)
        else:
            match_pattern = im_subset[0] + '.*tif'
            imlist = [f for f in glob('OIS*.tif') if re.search(match_pattern, f)]

    return imlist, match_pattern


def get_lowres_transform(lowres, ref, fprint, ref_mask, imgsource='DECLASSII'):

    if imgsource == 'DECLASSII':
        M = imtools.transform_from_fprint(lowres, ref, fprint, ref_mask)
    else:
        try:
            ortho_mask = 255 * np.ones(lowres.shape, dtype=np.uint8)
            ortho_mask[binary_dilation(lowres == 0, selem=disk(5))] = 0

            kp, des, matches = imtools.get_matches(lowres.astype(np.uint8),
                                                   ref.img.astype(np.uint8),
                                                   mask1=ortho_mask,
                                                   mask2=ref_mask,
                                                   dense=True)
            print('{} matches found.'.format(len(matches)))

            src_pts = np.array([kp[0][m.queryIdx].pt for m in matches])
            dst_pts = np.array([kp[1][m.trainIdx].pt for m in matches])

            M, inliers = ransac((dst_pts, src_pts), EuclideanTransform,
                                min_samples=10, residual_threshold=20, max_trials=10000)
            resids = M.residuals(dst_pts, src_pts)[inliers]
            print('{} matches used for initial transformation'.format(np.count_nonzero(inliers)))

        except ValueError:
            raise RuntimeError("Unable to find initial transformation with the parameters provided.")

    return M


def get_utm_str(epsg):
    """
    Given a GeoImg, read the associated EPSG code and return a UTM zone.

    Examples:
        - get_utm_str(32608) -> 'UTM Zone 8N'
        - get_utm_str(32708) -> 'UTM Zone 8S'

    :param int|str epsg: a str or int representing an EPSG Code
    :return:
        - **utm_str** (*str*) -- the UTM string representation
    """
    epsg_str = str(epsg)
    hemi_dict = {'6': 'N', '7': 'S'}

    if epsg_str[:2] == '32':
        utm_str = epsg_str[-2:] + hemi_dict[epsg_str[2]]
    else:
        utm_str = 'not utm'

    return utm_str


def transform_centers(img_gt, ref, imlist, footprints, ori):
    """
    Use the camera centers in relative space provided by MicMac Orientation files, along with camera footprints,
    to estimate a transformation between the relative coordinate system and the absolute coordinate system.

    :param array-like img_gt: the image GeoTransform (as provided by gdal.Dataset.GetGeoTransform)
    :param GeoImg ref: the reference image to use to determine the output image shape
    :param list imlist: a list of of the images that were used for the relative orthophoto
    :param GeoDataFrame footprints: the (approximate) image footprints - the centroid will be used for the absolute camera positions.
    :param str ori: name of orientation directory (after Ori-)
    :return:
        - **model** (*AffineTransform*) -- the estimated Affine Transformation between relative and absolute space
        - **inliers** (*array-like*) -- a list of the inliers returned by skimage.measure.ransac
        - **join** (*GeoDataFrame*) -- the joined image footprints and relative orientation files
    """

    rel_ori = mmtools.load_all_orientation(ori, imlist=imlist)

    footprints = footprints.to_crs(epsg=ref.epsg).copy()

    for i, row in footprints.iterrows():
        footprints.loc[i, 'x'] = row['geometry'].centroid.x
        footprints.loc[i, 'y'] = row['geometry'].centroid.y
        footprints.loc[i, 'name'] = 'OIS-Reech_' + row['ID'] + '.tif'

    join = footprints.set_index('name').join(rel_ori.set_index('name'), lsuffix='abs', rsuffix='rel')
    join.dropna(inplace=True)

    ref_ij = np.array([ref.xy2ij((row.xabs, row.yabs)) for i, row in join.iterrows()])
    rel_ij = np.array([((row.xrel - img_gt[4]) / img_gt[0],
                        (row.yrel - img_gt[5]) / img_gt[3]) for i, row in join.iterrows()])

    model, inliers = ransac((ref_ij[:, ::-1], rel_ij), AffineTransform, min_samples=10, residual_threshold=100, max_trials=5000)
    print('{} inliers for center transformation'.format(np.count_nonzero(inliers)))
    return model, inliers, join


def refine_lowres_tfm(rough_tfm, ref, mask, Minit):
    rough_gcps = imtools.find_grid_matches(rough_tfm, ref, mask, Minit, spacing=20, srcwin=40, dstwin=100)
    max_d = 100 - 40

    rough_gcps = rough_gcps[np.logical_and.reduce([rough_gcps.di.abs() < max_d,
                                                   rough_gcps.di < rough_gcps.di.median() + 2 * nmad(rough_gcps.di),
                                                   rough_gcps.di > rough_gcps.di.median() - 2 * nmad(rough_gcps.di),
                                                   rough_gcps.dj.abs() < max_d,
                                                   rough_gcps.dj < rough_gcps.dj.median() + 2 * nmad(rough_gcps.dj),
                                                   rough_gcps.dj > rough_gcps.dj.median() - 2 * nmad(rough_gcps.dj)])]

    Mref, inliers = ransac((rough_gcps[['search_j', 'search_i']].values,
                            rough_gcps[['orig_j', 'orig_i']].values),
                           AffineTransform, min_samples=10, residual_threshold=10, max_trials=5000)
    return Mref, inliers


def scale_factor(centers, img_gt, ref):

    ref_ij = np.array([ref.xy2ij((row.xabs, row.yabs)) for i, row in centers.iterrows()])
    rel_ij = np.array([((row.xrel - img_gt[4]) / img_gt[0],
                        (row.yrel - img_gt[5]) / img_gt[3]) for i, row in centers.iterrows()])
    ref_pt = Point(ref_ij[0])
    rel_pt = Point(rel_ij[0])

    ref_dist = []
    rel_dist = []

    for i in range(1, ref_ij.shape[0]):
        ref_dist.append(ref_pt.distance(Point(ref_ij[i])))
        rel_dist.append(rel_pt.distance(Point(rel_ij[i])))

    return np.nanmean(np.array(ref_dist) / np.array(rel_dist))


def get_mask(footprints, img, imlist, landmask=None, glacmask=None):
    """
    Create a mask for an image from different sources.

    :param GeoDataFrame footprints: vector data representing image footprints
    :param GeoImg img: the GeoImg to create a mask for
    :param array-like imlist: a list of image names
    :param str landmask: path to file of land outlines (i.e., an inclusion mask)
    :param str glacmask: path to file of glacier outlines (i.e., an exclusion mask)

    :returns:
        - **mask** (*array-like*) -- the mask
        - **fmask** (*GeoImg*) -- the georeferenced footprint mask
        - **img** (*GeoImg*) -- the GeoImg, cropped to a 10 pixel buffer around the image footprints
    """
    fmask, fprint = imtools.get_footprint_mask(footprints, img, imlist, fprint_out=True)
    fmask_geo = img.copy(new_raster=fmask)

    xmin, ymin, xmax, ymax = fprint.buffer(img.dx * 10).bounds

    img = img.crop_to_extent([xmin, xmax, ymin, ymax], pixel_size=img.dx)

    fmask = fmask_geo.crop_to_extent([xmin, xmax, ymin, ymax], pixel_size=img.dx).img == 1

    mask = imtools.make_binary_mask(img.img, erode=3, mask_value=np.nan)
    if landmask is not None:
        lmask = create_mask_from_shapefile(img, landmask)
        mask[~lmask] = 0
    if glacmask is not None:
        gmask = create_mask_from_shapefile(img, glacmask)
        mask[gmask] = 0

    mask[~fmask] = 0

    return mask, fmask, img


def register_ortho_old(fn_ortho, fn_ref, fn_reldem, fn_dem, glacmask=None, landmask=None, footprints=None,
                   im_subset=None, tfm_points=None, block_num=None, ori='Relative',
                   init_res=400, ortho_res=8, imgsource='DECLASSII', density=200, out_dir=None):
    
    print('start.')

    if out_dir is None:
        out_dir = 'auto_gcps'

    mkdir_p(out_dir)

    if block_num is not None:
        subscript = '_block{}'.format(block_num)
    else:
        subscript = ''

    ort_dir = os.path.dirname(fn_ortho)

    imlist, match_pattern = get_imlist(im_subset)

    ref_img = GeoImg(fn_ref)

    utm_str = get_utm_str(ref_img.epsg)

    ortho = imread(fn_ortho)

    ortho = median(ortho, selem=disk(1))

    ortho_ = Image.fromarray(ortho)

    lowres_ = int(init_res / ortho_res)

    ortho_lowres = np.array(ortho_.resize((np.array(ortho_.size) / lowres_).astype(int), Image.LANCZOS))

    ref_lowres = ref_img.resample(init_res)

    if footprints is not None:
        fmask, fprint = imtools.get_footprint_mask(footprints, ref_lowres, imlist, fprint_out=True)
    else:
        clean_imlist = [im.split('OIS-Reech_')[-1].split('.tif')[0] for im in imlist]
        print('Attempting to get image footprints from USGS EarthExplorer.')
        _fprints = get_usgs_footprints(clean_imlist, dataset=imgsource)
        fmask, fprint = imtools.get_footprint_mask(_fprints, ref_lowres, imlist, fprint_out=True)

    fmask_geo = ref_lowres.copy(new_raster=fmask)

    xmin, ymin, xmax, ymax = fprint.buffer(init_res * 10).bounds

    ref_img = ref_img.crop_to_extent([xmin, xmax, ymin, ymax], pixel_size=ref_img.dx)

    if isinstance(ref_lowres.img, np.ma.masked_array):
        ref_eq = imtools.stretch_image(ref_lowres.img.data.copy(), scale=(0.02, 0.98), mask=fmask)
    else:
        ref_eq = imtools.stretch_image(ref_lowres.img.copy(), scale=(0.02, 0.98), mask=fmask)    

    ref_lowres = ref_lowres.copy(new_raster=ref_eq, datatype=gdal.GDT_Byte)

    ref_lowres = ref_lowres.crop_to_extent([xmin, xmax, ymin, ymax], pixel_size=init_res)
    fmask = fmask_geo.crop_to_extent([xmin, xmax, ymin, ymax], pixel_size=init_res).img == 1

    lowres_mask = imtools.make_binary_mask(ref_lowres.img, erode=3, mask_value=np.nan)
    if landmask is not None:
        lmask = create_mask_from_shapefile(ref_lowres, landmask)
        lowres_mask[~lmask] = 0
    if glacmask is not None:
        gmask = create_mask_from_shapefile(ref_lowres, glacmask)
        lowres_mask[gmask] = 0

    lowres_mask[~fmask] = 0

    if tfm_points is not None:
        pts = pd.read_csv(tfm_points)
        src = pts[['J', 'I']].values / lowres_
        dst = np.array([ref_lowres.xy2ij((row.X, row.Y)) for i, row in pts.iterrows()])
        M = EuclideanTransform()
        M.estimate(src, dst[:, ::-1])
    else:
        M = get_lowres_transform(ortho_lowres, ref_lowres, fprint, lowres_mask, imgsource=imgsource)

    init_tfm = warp(ortho_lowres, M, output_shape=ref_lowres.img.shape, preserve_range=True)

    # if glacmask is not None:
    #   lowres_mask[gmask] = 255
    #    lowres_mask[~fmask] = 0

    rough_gcps = imtools.find_grid_matches(init_tfm, ref_lowres, lowres_mask, M, spacing=20, srcwin=40, dstwin=100)
    max_d = 100 - 40

    rough_gcps = rough_gcps[fmask[rough_gcps.search_i, rough_gcps.search_j] > 0]

    rough_gcps = rough_gcps[np.logical_and.reduce([rough_gcps.di.abs() < max_d,
                                                   rough_gcps.di < rough_gcps.di.median() + 2 * nmad(rough_gcps.di),
                                                   rough_gcps.di > rough_gcps.di.median() - 2 * nmad(rough_gcps.di),
                                                   rough_gcps.dj.abs() < max_d,
                                                   rough_gcps.dj < rough_gcps.dj.median() + 2 * nmad(rough_gcps.dj),
                                                   rough_gcps.dj > rough_gcps.dj.median() - 2 * nmad(rough_gcps.dj)])]

    Minit, inliers = ransac((rough_gcps[['search_j', 'search_i']].values,
                             rough_gcps[['orig_j', 'orig_i']].values),
                            AffineTransform, min_samples=10, residual_threshold=10, max_trials=5000)

    rough_gcps['residuals'] = Minit.residuals(rough_gcps[['search_j', 'search_i']].values,
                                              rough_gcps[['orig_j', 'orig_i']].values)
    # inliers = rough_gcps.residuals < min(10 * nmad(rough_gcps.residuals), 25)
    print('{} valid matches used for initial transformation'.format(np.count_nonzero(inliers)))

    ortho_tfm = warp(ortho_lowres, Minit, output_shape=ref_lowres.img.shape, preserve_range=True, order=5)

    dst_scale = ref_lowres.dx / ref_img.dx

    fig, ax = plt.subplots(1, 2, figsize=(7, 5))
    ax[0].imshow(ortho_tfm, cmap='gray')
    ax[1].imshow(ref_lowres.img, cmap='gray')

    plt.savefig('initial_transformation{}.png'.format(subscript), dpi=200, bbox_inches='tight')
    plt.close(fig)

    i_ = np.arange(0, ortho_lowres.shape[0], 10)
    j_ = np.arange(0, ortho_lowres.shape[1], 10)

    I, J = np.meshgrid(i_, j_)

    src_grd = np.array(list(zip(J.reshape(-1, 1), I.reshape(-1, 1)))).reshape(-1, 2)

    dst_tfm = []
    for pt in src_grd:
        dst_tfm.append(Minit.inverse(pt) * dst_scale)
    dst_tfm = np.array(dst_tfm).reshape(-1, 2)

    Minit_full, _ = ransac((dst_tfm, lowres_ * src_grd), AffineTransform, min_samples=3,
                           residual_threshold=1, max_trials=1000)
    # Minit_full = AffineTransform()
    # Minit_full.estimate((dst_tfm, lowres_ * src_grd))

    rough_tfm = warp(ortho, Minit_full, output_shape=ref_img.img.shape, preserve_range=True)

    mask_full = imtools.make_binary_mask(ref_img.img, mask_value=np.nan)
    if landmask is not None:
        lm = create_mask_from_shapefile(ref_img, landmask)
        mask_full[~lm] = 0
    if glacmask is not None:
        gm = create_mask_from_shapefile(ref_img, glacmask)
        mask_full[gm] = 0
    mask_full[rough_tfm == 0] = 0

    # for each of these pairs (src, dst), find the precise subpixel match (or not...)
    gcps = imtools.find_grid_matches(rough_tfm, ref_img, mask_full, Minit_full, spacing=density, dstwin=600)

    xy = np.array([ref_img.ij2xy((pt[1], pt[0])) for pt in gcps[['search_j', 'search_i']].values]).reshape(-1, 2)
    gcps['geometry'] = [Point(pt) for pt in xy]

    gcps = gcps[mask_full[gcps.search_i, gcps.search_j] == 255]

    max_d = 400 - 60

    # gcps = gcps[np.logical_and.reduce([gcps.di.abs() < max_d,
    #                                    gcps.di.abs() < 2 * nmad(gcps.di),
    #                                    gcps.dj.abs() < max_d,
    #                                    gcps.dj.abs() < 2 * nmad(gcps.dj)])]

    # gcps = gcps[gcps.z_corr > gcps.z_corr.quantile(0.5)]

    if os.path.exists(fn_ortho.replace('.tif', '.tfw')):
        fn_tfw = fn_ortho.replace('.tif', '.tfw')
        with open(fn_tfw, 'r') as f:
            gt = [float(l.strip()) for l in f.readlines()]

        gcps['rel_x'] = gt[4] + gcps['orig_j'].values * gt[0]  # need the original image coordinates
        gcps['rel_y'] = gt[5] + gcps['orig_i'].values * gt[3]

    gcps['elevation'] = 0
    gcps.crs = ref_img.proj4

    gcps.dropna(inplace=True)

    print('loading dems')
    dem = GeoImg(fn_dem, dtype=np.float32)
    rel_dem = GeoImg(fn_reldem)

    # gcps.crs = {'init': 'epsg:{}'.format(ref_img.epsg)}
    for i, row in gcps.to_crs(crs=dem.proj4).iterrows():
        gcps.loc[i, 'elevation'] = dem.raster_points([(row.geometry.x, row.geometry.y)], nsize=3, mode='linear')
        if os.path.exists(fn_ortho.replace('.tif', '.tfw')):
            gcps.loc[i, 'el_rel'] = rel_dem.raster_points([(row.rel_x, row.rel_y)], nsize=3, mode='linear')

    # drop any gcps where we don't have a DEM value or a valid match
    gcps.dropna(inplace=True)

    # run ransac to find the matches between the transformed image and the master image make a coherent transformation
    # residual_threshold is 10 pixels to allow for some local distortions, but get rid of the big blunders
    Mref, inliers_ref = ransac((gcps[['search_j', 'search_i']].values, gcps[['match_j', 'match_i']].values),
                               AffineTransform, min_samples=6, residual_threshold=20, max_trials=5000)

    gcps['aff_resid'] = Mref.residuals(gcps[['search_j', 'search_i']].values,
                                       gcps[['match_j', 'match_i']].values)

    gcps_orig = gcps.copy()

    # gcps = gcps[inliers_ref]

    out = sliding_window_filter([ortho.shape[1], ortho.shape[0]], gcps,
                                min(1000, ortho.shape[1] / 4, ortho.shape[0] / 4),
                                mindist=500, how='aff_resid', is_ascending=False)
    gcps = gcps.loc[out]

    # Mfin, inliers_fin = ransac((gcps[['search_j', 'search_i']].values, gcps[['match_j', 'match_i']].values),
    #                            AffineTransform, min_samples=6, residual_threshold=20, max_trials=5000)
    # gcps = gcps[inliers_fin]

    print('{} valid matches found'.format(gcps.shape[0]))

    gcps.index = range(gcps.shape[0])  # make sure index corresponds to row we're writing out
    gcps['id'] = ['GCP{}'.format(i) for i in range(gcps.shape[0])]
    gcps.to_file(os.path.join(out_dir, 'AutoGCPs{}.shp'.format(subscript)))

    print('writing AutoGCPs.txt')
    mmtools.write_auto_gcps(gcps, subscript, out_dir, utm_str)

    print('converting AutoGCPs.txt to AutoGCPs.xml')
    subprocess.Popen(['mm3d', 'GCPConvert', 'AppInFile',
                      os.path.join(out_dir, 'AutoGCPs{}.txt'.format(subscript))]).wait()

    print('writing AutoMeasures.txt')
    mmtools.write_auto_mesures(gcps, subscript, out_dir)

    print('running get_autogcp_locations.sh to get rough image locations for each point')
    subprocess.Popen(['get_autogcp_locations.sh', 'Ori-{}'.format(ori),
                      os.path.join(out_dir, 'AutoMeasures{}.txt'.format(subscript))] + imlist).wait()

    # print('searching for points in orthorectified images')
    print('finding image measures')
    mmtools.write_image_mesures(imlist, gcps, out_dir, subscript, ort_dir=ort_dir)

    print('running mm3d GCPBascule to estimate terrain errors')
    gcps = mmtools.run_bascule(gcps, out_dir, match_pattern, subscript, ori)
    gcps['res_dist'] = np.sqrt(gcps.xres**2 + gcps.yres**2)

    gcps = gcps[gcps.res_dist < 4 * nmad(gcps.res_dist)]

    mmtools.save_gcps(gcps, out_dir, utm_str, subscript)
    gcps = mmtools.run_bascule(gcps, out_dir, match_pattern, subscript, ori)
    gcps['res_dist'] = np.sqrt(gcps.xres**2 + gcps.yres**2)

    gcps = gcps[gcps.res_dist < 4 * nmad(gcps.res_dist)]

    mmtools.save_gcps(gcps, out_dir, utm_str, subscript)

    gcps = mmtools.run_bascule(gcps, out_dir, match_pattern, subscript, ori)
    gcps = mmtools.run_campari(gcps, out_dir, match_pattern, subscript, ref_img.dx, ortho_res)
    gcps['camp_dist'] = np.sqrt(gcps.camp_xres**2 + gcps.camp_yres**2)

    # if args.imgsource != 'DECLASSII':
    #     out = sliding_window_filter([ortho.shape[1], ortho.shape[0]], gcps, 2000, mindist=1000, how='camp_dist')
    #    gcps = gcps.loc[out]
    niter = 0
    while any([np.any(gcps.camp_res > 4 * nmad(gcps.camp_res)),
               np.any(gcps.camp_dist > 4 * nmad(gcps.camp_dist)),
               gcps.camp_res.max() > 2]) and niter <= 5:
        valid_inds = np.logical_and.reduce((gcps.camp_res < 4 * nmad(gcps.camp_res),
                                           gcps.camp_res < gcps.camp_res.max(),
                                           gcps.z_corr > gcps.z_corr.min()))
        if np.count_nonzero(valid_inds) < 10:
            break

        gcps = gcps.loc[valid_inds]
        mmtools.save_gcps(gcps, out_dir, utm_str, subscript)
        gcps = mmtools.run_bascule(gcps, out_dir, match_pattern, subscript, ori)
        gcps['res_dist'] = np.sqrt(gcps.xres ** 2 + gcps.yres ** 2)

        gcps = mmtools.run_campari(gcps, out_dir, match_pattern, subscript, ref_img.dx, ortho_res)
        gcps['camp_dist'] = np.sqrt(gcps.camp_xres**2 + gcps.camp_yres**2)
        niter += 1

    # final write of gcps to disk.
    gcps.to_file(os.path.join(out_dir, 'AutoGCPs{}.shp'.format(subscript)))

    fig1 = plt.figure(figsize=(7, 5))
    plt.imshow(ortho[::5, ::5], cmap='gray', extent=[0, ortho.shape[1], ortho.shape[0], 0])
    plt.plot(gcps.orig_j, gcps.orig_i, 'r+')
    plt.quiver(gcps.orig_j, gcps.orig_i, gcps.camp_xres, gcps.camp_yres, color='r')

    plt.savefig(os.path.join(out_dir, 'relative_gcps{}.png'.format(subscript)), bbox_inches='tight', dpi=200)
    plt.close(fig1)

    fig2 = ref_img.display(sfact=10, fig=plt.figure(figsize=(7, 5)))
    plt.plot(gcps.geometry.x, gcps.geometry.y, 'r+')
    plt.savefig(os.path.join(out_dir, 'world_gcps{}.png'.format(subscript)), bbox_inches='tight', dpi=200)
    plt.close(fig2)

    print('cleaning up.')
    # remove Auto-im.tif.txt, NoDist-im.tif.txt, etc. Ori-Relative-NoDist
    shutil.rmtree('Ori-{}-NoDist/'.format(ori))

    for txtfile in glob('Auto-OIS*.tif.txt') + \
                   glob('NoDist-OIS*.tif.txt'):
        os.remove(txtfile)
    print('end.')
    # embed()


def register_ortho(fn_ortho, fn_ref, fn_reldem, fn_dem, glacmask=None, landmask=None, footprints=None,
                   im_subset=None, block_num=None, ori='Relative', ortho_res=8.,
                   imgsource='DECLASSII', density=200, out_dir=None, allfree=True):
    """
    Register a relative orthoimage and DEM to a reference orthorectified image and DEM.

    :param str fn_ortho: path to relative orthoimage
    :param str fn_ref: path to reference orthorectified image
    :param str fn_reldem: path to relative DEM
    :param str fn_dem: path to reference DEM
    :param str glacmask: path to file of glacier outlines (i.e., an exclusion mask)
    :param str landmask: path to file of land outlines (i.e., an inclusion mask)
    :param str footprints: path to shapefile of image outlines. If not set, will download from USGS.
    :param str im_subset: subset of raw images to work with
    :param str block_num: block number to use if processing multiple image blocks
    :param str ori: name of orientation directory (after Ori-) (default: Relative)
    :param float ortho_res: approx. ground sampling distance (pixel resolution) of ortho image (default: 8 m)
    :param str imgsource: USGS dataset name for images (default: DECLASSII)
    :param int density: pixel spacing to look for GCPs (default: 200)
    :param str out_dir: output directory to save auto GCP files to (default: auto_gcps)
    :param bool allfree: run Campari setting all parameters free (default: True)
    """
    print('start.')

    if out_dir is None:
        out_dir = 'auto_gcps'

    mkdir_p(out_dir)

    if block_num is not None:
        subscript = '_block{}'.format(block_num)
    else:
        subscript = ''

    ort_dir = os.path.dirname(fn_ortho)

    imlist, match_pattern = get_imlist(im_subset)

    ref_img = GeoImg(fn_ref)
    # ref_lowres = ref_img.resample(init_res)
    # resamp_fact = init_res / ref_img.dx

    utm_str = get_utm_str(ref_img.epsg)

    ortho = imread(fn_ortho)
    # ortho_ = Image.fromarray(ortho)

    fn_tfw = fn_ortho.replace('.tif', '.tfw')
    with open(fn_tfw, 'r') as f:
        ortho_gt = [float(l.strip()) for l in f.readlines()]

    if footprints is None:
        clean_imlist = [im.split('OIS-Reech_')[-1].split('.tif')[0] for im in imlist]
        print('Attempting to get image footprints from USGS EarthExplorer.')
        footprints = get_usgs_footprints(clean_imlist, dataset=imgsource)
    else:
        footprints = gpd.read_file(footprints)

    mask_full, _, ref_img = get_mask(footprints, ref_img, imlist, landmask, glacmask)

    Minit, _, centers = transform_centers(ortho_gt, ref_img, imlist, footprints, 'Ori-{}'.format(ori))
    rough_tfm = warp(ortho, Minit, output_shape=ref_img.shape, preserve_range=True)

    fig, ax = plt.subplots(1, 2, figsize=(7, 5))
    ax[0].imshow(rough_tfm[::10, ::10], extent=[0, rough_tfm.shape[1], rough_tfm.shape[0], 0], cmap='gray')
    ax[1].imshow(ref_img.img[::10, ::10], extent=[0, ref_img.shape[1], ref_img.shape[0], 0], cmap='gray')

    plt.savefig('initial_transformation{}.png'.format(subscript), dpi=200, bbox_inches='tight')
    plt.close(fig)

    # for each of these pairs (src, dst), find the precise subpixel match (or not...)
    gcps = imtools.find_grid_matches(rough_tfm, ref_img, mask_full, Minit, spacing=density, dstwin=400)

    xy = np.array([ref_img.ij2xy((pt[1], pt[0])) for pt in gcps[['search_j', 'search_i']].values]).reshape(-1, 2)
    gcps['geometry'] = [Point(pt) for pt in xy]

    gcps = gcps[mask_full[gcps.search_i, gcps.search_j] == 255]

    gcps['rel_x'] = ortho_gt[4] + gcps['orig_j'].values * ortho_gt[0]  # need the original image coordinates
    gcps['rel_y'] = ortho_gt[5] + gcps['orig_i'].values * ortho_gt[3]

    gcps['elevation'] = 0
    gcps.crs = ref_img.proj4

    gcps.dropna(inplace=True)

    print('loading dems')
    dem = GeoImg(fn_dem, dtype=np.float32)
    rel_dem = GeoImg(fn_reldem)

    # gcps.crs = {'init': 'epsg:{}'.format(ref_img.epsg)}
    for i, row in gcps.to_crs(crs=dem.proj4).iterrows():
        gcps.loc[i, 'elevation'] = dem.raster_points([(row.geometry.x, row.geometry.y)], nsize=3, mode='linear')
        if os.path.exists(fn_ortho.replace('.tif', '.tfw')):
            gcps.loc[i, 'el_rel'] = rel_dem.raster_points([(row.rel_x, row.rel_y)], nsize=3, mode='linear')

    # drop any gcps where we don't have a DEM value or a valid match
    gcps.dropna(inplace=True)

    # run ransac to find the matches between the transformed image and the master image make a coherent transformation
    # residual_threshold is 10 pixels to allow for some local distortions, but get rid of the big blunders
    Mref, inliers_ref = ransac((gcps[['search_j', 'search_i']].values, gcps[['match_j', 'match_i']].values),
                               AffineTransform, min_samples=6, residual_threshold=10, max_trials=5000)

    gcps['aff_resid'] = Mref.residuals(gcps[['search_j', 'search_i']].values,
                                       gcps[['match_j', 'match_i']].values)

    out = sliding_window_filter([ortho.shape[1], ortho.shape[0]], gcps,
                                min(1000, ortho.shape[1] / 4, ortho.shape[0] / 4),
                                mindist=500, how='aff_resid', is_ascending=False)
    gcps = gcps.loc[out]

    print('{} valid matches found'.format(gcps.shape[0]))

    gcps.index = range(gcps.shape[0])  # make sure index corresponds to row we're writing out
    gcps['id'] = ['GCP{}'.format(i) for i in range(gcps.shape[0])]
    gcps.to_file(os.path.join(out_dir, 'AutoGCPs{}.shp'.format(subscript)))

    print('writing AutoGCPs.txt')
    mmtools.write_auto_gcps(gcps, subscript, out_dir, utm_str)

    print('converting AutoGCPs.txt to AutoGCPs.xml')
    subprocess.Popen(['mm3d', 'GCPConvert', 'AppInFile',
                      os.path.join(out_dir, 'AutoGCPs{}.txt'.format(subscript))]).wait()

    print('writing AutoMeasures.txt')
    mmtools.write_auto_mesures(gcps, subscript, out_dir)

    print('running get_autogcp_locations.sh to get rough image locations for each point')
    subprocess.Popen(['get_autogcp_locations.sh', 'Ori-{}'.format(ori),
                      os.path.join(out_dir, 'AutoMeasures{}.txt'.format(subscript))] + imlist).wait()

    # print('searching for points in orthorectified images')
    print('finding image measures')
    mmtools.write_image_mesures(imlist, gcps, out_dir, subscript, ort_dir=ort_dir)

    print('running mm3d GCPBascule to estimate terrain errors')
    gcps = mmtools.run_bascule(gcps, out_dir, match_pattern, subscript, ori)
    gcps['res_dist'] = np.sqrt(gcps.xres ** 2 + gcps.yres ** 2)

    gcps = gcps[gcps.res_dist < 4 * nmad(gcps.res_dist)]

    mmtools.save_gcps(gcps, out_dir, utm_str, subscript)
    gcps = mmtools.run_bascule(gcps, out_dir, match_pattern, subscript, ori)
    gcps['res_dist'] = np.sqrt(gcps.xres ** 2 + gcps.yres ** 2)

    gcps = gcps[gcps.res_dist < 4 * nmad(gcps.res_dist)]

    mmtools.save_gcps(gcps, out_dir, utm_str, subscript)

    gcps = mmtools.run_bascule(gcps, out_dir, match_pattern, subscript, ori)
    gcps = mmtools.run_campari(gcps, out_dir, match_pattern, subscript, ref_img.dx, ortho_res, allfree=allfree)
    gcps['camp_dist'] = np.sqrt(gcps.camp_xres ** 2 + gcps.camp_yres ** 2)

    niter = 0
    while any([np.any(gcps.camp_res > 4 * nmad(gcps.camp_res)),
               np.any(gcps.camp_dist > 4 * nmad(gcps.camp_dist)),
               gcps.camp_res.max() > 2]) and niter <= 5:
        valid_inds = np.logical_and.reduce((gcps.camp_res < 4 * nmad(gcps.camp_res),
                                            gcps.camp_res < gcps.camp_res.max(),
                                            gcps.z_corr > gcps.z_corr.min()))
        if np.count_nonzero(valid_inds) < 10:
            break

        gcps = gcps.loc[valid_inds]
        mmtools.save_gcps(gcps, out_dir, utm_str, subscript)
        gcps = mmtools.run_bascule(gcps, out_dir, match_pattern, subscript, ori)
        gcps['res_dist'] = np.sqrt(gcps.xres ** 2 + gcps.yres ** 2)

        gcps = mmtools.run_campari(gcps, out_dir, match_pattern, subscript, ref_img.dx, ortho_res, allfree=allfree)
        gcps['camp_dist'] = np.sqrt(gcps.camp_xres ** 2 + gcps.camp_yres ** 2)
        niter += 1

    # final write of gcps to disk.
    gcps.to_file(os.path.join(out_dir, 'AutoGCPs{}.shp'.format(subscript)))

    fig1 = plt.figure(figsize=(7, 5))
    plt.imshow(ortho[::5, ::5], cmap='gray', extent=[0, ortho.shape[1], ortho.shape[0], 0])
    plt.plot(gcps.orig_j, gcps.orig_i, 'r+')
    plt.quiver(gcps.orig_j, gcps.orig_i, gcps.camp_xres, gcps.camp_yres, color='r')

    plt.savefig(os.path.join(out_dir, 'relative_gcps{}.png'.format(subscript)), bbox_inches='tight', dpi=200)
    plt.close(fig1)

    fig2 = ref_img.display(sfact=10, fig=plt.figure(figsize=(7, 5)))
    plt.plot(gcps.geometry.x, gcps.geometry.y, 'r+')
    plt.savefig(os.path.join(out_dir, 'world_gcps{}.png'.format(subscript)), bbox_inches='tight', dpi=200)
    plt.close(fig2)

    print('cleaning up.')
    # remove Auto-im.tif.txt, NoDist-im.tif.txt, etc. Ori-Relative-NoDist
    shutil.rmtree('Ori-{}-NoDist/'.format(ori))

    for txtfile in glob('Auto-OIS*.tif.txt') + \
                   glob('NoDist-OIS*.tif.txt'):
        os.remove(txtfile)
    print('end.')
    # embed()
