"""
spymicmac.register is a collection of tools for registering images and finding GCPs.
"""
import os
import re
import subprocess
import shutil
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import lxml.etree as etree
import lxml.builder as builder
from glob import glob
from rtree import index
from shapely.ops import unary_union
from shapely.geometry.point import Point
from skimage.io import imread
from skimage.measure import ransac
from skimage.transform import AffineTransform, warp
from pybob.ddem_tools import nmad
from pybob.image_tools import create_mask_from_shapefile
from pybob.GeoImg import GeoImg
from spymicmac import data, image, matching, micmac, orientation


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
    :return: **out_inds** (*array-like*) -- an array with the filtered indices
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


def get_imlist(im_subset, dirname='.', strip_text=None):
    """
    Given either a list of filenames or a regex match pattern, return a list of filenames and a pattern to provide
    to MicMac.

    :param list|str im_subset: a list of filenames or a match pattern (e.g., ['OIS.*tif']) representing filenames
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
            match_pattern = micmac.get_match_pattern(imlist)
            # match_pattern = '|'.join(imlist)
        else:
            match_pattern = im_subset[0] + '.*tif'
            imlist = [f for f in glob('OIS*.tif') if re.search(match_pattern, f)]

    imlist.sort()

    return imlist, match_pattern


def get_utm_str(epsg):
    """
    Given a GeoImg, read the associated EPSG code and return a UTM zone.

    Examples:
        - get_utm_str(32608) -> 'UTM Zone 8N'
        - get_utm_str(32708) -> 'UTM Zone 8S'

    :param int|str epsg: a str or int representing an EPSG Code
    :return: **utm_str** (*str*) -- the UTM string representation
    """
    epsg_str = str(epsg)
    hemi_dict = {'6': 'N', '7': 'S'}

    if epsg_str[:2] == '32':
        utm_str = epsg_str[-2:] + hemi_dict[epsg_str[2]]
    else:
        utm_str = 'not utm'

    return utm_str


def warp_image(model, ref, img):
    """
    Given a transformation model between two coordinate systems, warp an image to a reference image

    :param GeometricTransform model: the transformation model between the coordinate systems
    :param GeoImg ref: the reference GeoImg
    :param GeoImg img: the GeoImg to be transformed
    :return:
        - **tfm_img** (*np.array*) -- the input image transformed to the same extent as the reference image
        - **this_model** (*AffineTransform*) -- the estimated Affine Transformation between the two images
        - **inliers** (*array-like*) -- a list of the inliers returned by skimage.measure.ransac
    """

    # get image points in rel coords
    img_x, img_y = img.xy()
    img_x = img_x[::100, ::100].flatten()
    img_y = img_y[::100, ::100].flatten()

    rel_pts = np.hstack((img_x.reshape(-1, 1), img_y.reshape(-1, 1)))
    # get image corners, center in abs coords
    ref_pts = model.inverse(rel_pts)

    # get the new transformation - have to go from gdal geotransform to tfw geotransform first
    this_model, inliers = orientation.transform_points(ref, ref_pts, _to_tfw(img.gt), rel_pts)

    # transform the image
    tfm_img = warp(img.img, this_model, output_shape=ref.img.shape, preserve_range=True)

    return tfm_img, this_model, inliers


def _to_tfw(gt):
    # convert from gdal GeoTransform to TFW format
    return [gt[1], gt[2], gt[4], gt[5], gt[0], gt[3]]


def get_footprint_overlap(fprints):
    """
    Return the area where image footprints overlap.

    :param GeoDataFrame fprints: a GeoDataFrame of image footprints
    :return: **intersection** (*shapely.Polygon*) -- the overlapping area (unary union) of the images.
    """
    if fprints.shape[0] == 1:
        return fprints.geometry.values[0]

    idx = index.Index()

    for pos, row in fprints.iterrows():
        idx.insert(pos, row['geometry'].bounds)

    intersections = []
    for poly in fprints['geometry']:
        merged = unary_union([fprints.loc[pos, 'geometry'] for pos in idx.intersection(poly.bounds)
                              if fprints.loc[pos, 'geometry'] != poly])
        intersections.append(poly.intersection(merged))

    intersection = unary_union(intersections)

    # return intersection.minimum_rotated_rectangle
    return intersection


def get_footprint_mask(shpfile, geoimg, filelist, fprint_out=False):
    """
    Return a footprint mask for an image.

    :param str|GeoDataFrame shpfile: a filename or a GeoDataFrame representation of the footprints.
    :param GeoImg geoimg: the image to create a mask for.
    :param list filelist: a list of image names to use to create the footprint mask.
    :param bool fprint_out: return the polygon representation of the footprint mask.
    :return:
        - **mask** (*array-like*) -- the footprint mask
        - **fprint** (*shapely.Polygon*) -- the footprint polygon, if requested.
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


# copied from pymmaster.mmaster_tools
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
    fmask, fprint = get_footprint_mask(footprints, img, imlist, fprint_out=True)
    fmask_geo = img.copy(new_raster=fmask)

    xmin, ymin, xmax, ymax = fprint.buffer(img.dx * 10).bounds

    img = img.crop_to_extent([xmin, xmax, ymin, ymax], pixel_size=img.dx)

    fmask = fmask_geo.crop_to_extent([xmin, xmax, ymin, ymax], pixel_size=img.dx).img == 1

    mask = image.make_binary_mask(img.img, erode=3, mask_value=np.nan)
    if landmask is not None:
        lmask = create_mask_from_shapefile(img, landmask)
        mask[~lmask] = 0
    if glacmask is not None:
        gmask = create_mask_from_shapefile(img, glacmask)
        mask[gmask] = 0

    mask[~fmask] = 0

    return mask, fmask, img


def _search_size(imshape):
    min_dst = 100
    dst = min(400, imshape[0] / 10, imshape[1] / 10)
    return int(max(min_dst, dst))


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

    os.makedirs(out_dir, exist_ok=True)

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
        footprints = data.get_usgs_footprints(clean_imlist, dataset=imgsource)
    else:
        footprints = gpd.read_file(footprints)

    mask_full, _, ref_img = get_mask(footprints, ref_img, imlist, landmask, glacmask)

    Minit, _, centers = orientation.transform_centers(ortho_gt, ref_img, imlist, footprints, 'Ori-{}'.format(ori))
    rough_tfm = warp(ortho, Minit, output_shape=ref_img.shape, preserve_range=True)

    rough_spacing = np.round(max(ref_img.shape) / 20 / 1000) * 1000

    rough_gcps = matching.find_grid_matches(rough_tfm, ref_img, mask_full, Minit,
                                            spacing=int(rough_spacing), dstwin=int(rough_spacing))

    model, inliers = ransac((rough_gcps[['search_j', 'search_i']].values,
                             rough_gcps[['orig_j', 'orig_i']].values), AffineTransform,
                            min_samples=6, residual_threshold=200, max_trials=5000)

    rough_tfm = warp(ortho, model, output_shape=ref_img.shape, preserve_range=True)

    rough_geo = ref_img.copy(new_raster=rough_tfm)
    rough_geo.write('Orthophoto{}_geo.tif'.format(subscript), dtype=np.uint8)

    fig, ax = plt.subplots(1, 2, figsize=(7, 5))
    ax[0].imshow(rough_tfm[::10, ::10], extent=[0, rough_tfm.shape[1], rough_tfm.shape[0], 0], cmap='gray')
    ax[1].imshow(ref_img.img[::10, ::10], extent=[0, ref_img.shape[1], ref_img.shape[0], 0], cmap='gray')

    fig.savefig('initial_transformation{}.png'.format(subscript), dpi=200, bbox_inches='tight')
    plt.close(fig)

    # for each of these pairs (src, dst), find the precise subpixel match (or not...)
    gcps = matching.find_grid_matches(rough_tfm, ref_img, mask_full, model,
                                      spacing=density, dstwin=_search_size(rough_tfm.shape))

    xy = np.array([ref_img.ij2xy((pt[1], pt[0])) for pt in gcps[['search_j', 'search_i']].values]).reshape(-1, 2)
    gcps['geometry'] = [Point(pt) for pt in xy]

    gcps = gcps.loc[mask_full[gcps.search_i, gcps.search_j] == 255]

    gcps['rel_x'] = ortho_gt[4] + gcps['orig_j'].values * ortho_gt[0]  # need the original image coordinates
    gcps['rel_y'] = ortho_gt[5] + gcps['orig_i'].values * ortho_gt[3]

    gcps['elevation'] = 0
    gcps.set_crs(crs=ref_img.proj4, inplace=True)

    gcps.dropna(inplace=True)
    print('{} potential matches found'.format(gcps.shape[0]))

    print('loading dems')
    dem = GeoImg(fn_dem, dtype=np.float32)
    rel_dem = GeoImg(fn_reldem)

    # gcps.crs = {'init': 'epsg:{}'.format(ref_img.epsg)}
    for i, row in gcps.to_crs(crs=dem.proj4).iterrows():
        gcps.loc[i, 'elevation'] = dem.raster_points([(row.geometry.x, row.geometry.y)], nsize=3, mode='linear')
        if os.path.exists(fn_ortho.replace('.tif', '.tfw')):
            gcps.loc[i, 'el_rel'] = rel_dem.raster_points([(row.rel_x, row.rel_y)], nsize=3, mode='linear')

    # drop any gcps where we don't have a DEM value or a valid match
    gcps.loc[np.abs(gcps.elevation - dem.NDV) < 1, 'elevation'] = np.nan
    gcps.dropna(inplace=True)
    print('{} matches with valid elevations'.format(gcps.shape[0]))

    # run ransac to find the matches between the transformed image and the master image make a coherent transformation
    # residual_threshold is 20 pixels to allow for some local distortions, but get rid of the big blunders
    Mref, inliers_ref = ransac((gcps[['search_j', 'search_i']].values, gcps[['match_j', 'match_i']].values),
                               AffineTransform, min_samples=6, residual_threshold=20, max_trials=5000)
    gcps['aff_resid'] = Mref.residuals(gcps[['search_j', 'search_i']].values,
                                       gcps[['match_j', 'match_i']].values)

    # valid = np.abs(gcps.aff_resid - gcps.aff_resid.median()) < nmad(gcps.aff_resid)
    # gcps = gcps.loc[valid]
    gcps = gcps.loc[inliers_ref]

    # out = sliding_window_filter([ortho.shape[1], ortho.shape[0]], gcps,
    #                             min(1000, ortho.shape[1] / 4, ortho.shape[0] / 4),
    #                             mindist=500, how='aff_resid', is_ascending=False)
    # gcps = gcps.loc[out]

    print('{} valid matches found after estimating transformation'.format(gcps.shape[0]))

    gcps.index = range(gcps.shape[0])  # make sure index corresponds to row we're writing out
    gcps['id'] = ['GCP{}'.format(i) for i in range(gcps.shape[0])]
    gcps.to_file(os.path.join(out_dir, 'AutoGCPs{}.shp'.format(subscript)))

    print('writing AutoGCPs.txt')
    micmac.write_auto_gcps(gcps, subscript, out_dir, utm_str)

    print('converting AutoGCPs.txt to AutoGCPs.xml')
    subprocess.Popen(['mm3d', 'GCPConvert', 'AppInFile',
                      os.path.join(out_dir, 'AutoGCPs{}.txt'.format(subscript))]).wait()

    print('writing AutoMeasures.txt')
    micmac.write_auto_mesures(gcps, subscript, out_dir)

    print('running get_autogcp_locations to get rough image locations for each point')
    micmac.get_autogcp_locations(f'Ori-{ori}', os.path.join(out_dir, f'AutoMeasures{subscript}.txt'), imlist)

    # print('searching for points in orthorectified images')
    print('finding image measures')
    micmac.write_image_mesures(imlist, gcps, out_dir, subscript, ort_dir=ort_dir)

    print('running mm3d GCPBascule to estimate terrain errors')
    gcps = micmac.bascule(gcps, out_dir, match_pattern, subscript, ori)
    gcps['res_dist'] = np.sqrt(gcps.xres ** 2 + gcps.yres ** 2)

    gcps = gcps.loc[np.abs(gcps.res_dist - gcps.res_dist.median()) < 2 * nmad(gcps.res_dist)]
    # gcps = gcps[np.logical_and(np.abs(gcps.xres - gcps.xres.median()) < 2 * nmad(gcps.xres),
    #                            np.abs(gcps.yres - gcps.yres.median()) < 2 * nmad(gcps.yres))]

    micmac.save_gcps(gcps, out_dir, utm_str, subscript)
    gcps = micmac.bascule(gcps, out_dir, match_pattern, subscript, ori)
    gcps['res_dist'] = np.sqrt(gcps.xres ** 2 + gcps.yres ** 2)
    gcps = gcps.loc[np.abs(gcps.res_dist - gcps.res_dist.median()) < 2 * nmad(gcps.res_dist)]
    # gcps = gcps[np.logical_and(np.abs(gcps.xres - gcps.xres.median()) < 2 * nmad(gcps.xres),
    #                            np.abs(gcps.yres - gcps.yres.median()) < 2 * nmad(gcps.yres))]

    micmac.save_gcps(gcps, out_dir, utm_str, subscript)

    # now, iterate campari to refine the orientation
    gcps = micmac.iterate_campari(gcps, out_dir, match_pattern, subscript, ref_img.dx, ortho_res, allfree=allfree)

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


def register_individual(dir_ortho, fn_ref, fn_reldem, fn_dem, glacmask=None, landmask=None, footprints=None,
                        im_subset=None, block_num=None, ori='Relative', ortho_res=8., imgsource='DECLASSII', density=200,
                        out_dir='auto_gcps', allfree=True, is_geo=False):
    """
    IN DEVELOPMENT: Register individual orthoimages to a reference orthorectified image and DEM.

    :param str dir_ortho: directory containing orthoimages (e.g., Ortho-MEC-Relative)
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
    :param bool is_geo: True if orthoimages are already in absolute coordinates (default: False)
    """

    if block_num is not None:
        subscript = '_block{}'.format(block_num)
    else:
        subscript = ''

    # TODO: fix this
    if im_subset is None:
        imlist = [os.path.basename(fn).split('Ort_')[1] for fn in glob(os.path.join(dir_ortho, 'Ort_*.tif'))]
        match_pattern = 'OIS.*tif'
    else:
        imlist = [os.path.basename(fn).split('Ort_')[1] for fn in glob(os.path.join(dir_ortho, 'Ort_*.tif')) if fn in im_subset]

    os.makedirs('TmpMeasures{}'.format(subscript), exist_ok=True)
    if out_dir is not None:
        os.makedirs(out_dir, exist_ok=True)
    else:
        out_dir = '.'

    ref_img = GeoImg(fn_ref)
    ref_dem = GeoImg(fn_dem)

    utm_str = get_utm_str(ref_img.epsg)

    rel_dem = GeoImg(fn_reldem)

    if footprints is None:
        clean_imlist = [im.split('OIS-Reech_')[-1].split('.tif')[0] for im in imlist]
        print('Attempting to get image footprints from USGS EarthExplorer.')
        footprints = data.get_usgs_footprints(clean_imlist, dataset=imgsource)
    else:
        footprints = gpd.read_file(footprints)

    mask_full, _, ref_img = get_mask(footprints, ref_img, imlist, landmask, glacmask)
    mask_geo = ref_img.copy(new_raster=mask_full)

    if not is_geo:
        with open(os.path.join(dir_ortho, 'Orthophotomosaic.tfw'), 'r') as f:
            ortho_gt = [float(l.strip()) for l in f.readlines()]
        model, _, centers = orientation.transform_centers(ortho_gt, ref_img, imlist,
                                                          footprints, 'Ori-{}'.format(ori), imgeom=False)

    # get a grid of points from the reference image, spaced by $density
    x, y = ref_img.xy()
    x = x[::density, ::density].flatten()
    y = y[::density, ::density].flatten()
    mask_pts = mask_full[::density, ::density].flatten()

    x = x[mask_pts != 0]
    y = y[mask_pts != 0]

    all_pts = np.hstack((x.reshape(-1, 1), y.reshape(-1, 1)))
    all_geom = [Point(pt) for pt in all_pts]
    pt_names = np.array(['GCP{}'.format(ii) for ii in range(len(all_pts))])

    gcps = gpd.GeoDataFrame()
    gcps['id'] = pt_names
    gcps['geometry'] = all_geom
    gcps.set_crs(epsg=ref_img.epsg, inplace=True)

    for ii, row in gcps.to_crs(crs=ref_dem.proj4).iterrows():
        gcps.loc[ii, 'elevation'] = ref_dem.raster_points([(row.geometry.x, row.geometry.y)], nsize=3, mode='linear')

    gcps.to_file(os.path.join(out_dir, 'AutoGCPs{}.shp'.format(subscript)))

    print('writing AutoGCPs.txt')
    micmac.write_auto_gcps(gcps, subscript, out_dir, utm_str)

    print('converting AutoGCPs.txt to AutoGCPs.xml')
    subprocess.Popen(['mm3d', 'GCPConvert', 'AppInFile',
                      os.path.join(out_dir, 'AutoGCPs{}.txt'.format(subscript))]).wait()

    E = builder.ElementMaker()
    MesureSet = E.SetOfMesureAppuisFlottants()

    for im in imlist:
        print(im)
        # have to reproject image if needed
        this_img = GeoImg('{}/Ort_{}'.format(dir_ortho, im))

        if is_geo:
            this_ref = ref_img.reproject(this_img)
            this_mask = mask_geo.reproject(this_img)
            this_mask.img[this_img.img == 0] = 0
        else:
            this_ext = np.array([(this_img.xmin, this_img.ymin),
                                 (this_img.xmax, this_img.ymin),
                                 (this_img.xmax, this_img.ymax),
                                 (this_img.xmin, this_img.ymax)])
            warp_ext = model.inverse(this_ext)
            this_x = warp_ext[:, 0]
            this_y = warp_ext[:, 1]

            this_ref = ref_img.crop_to_extent([min(this_x), max(this_x), min(this_y), max(this_y)])
            tfm_img, this_model, _ = warp_image(model, this_ref, this_img)

            orig_img = this_img.copy()
            this_img = this_ref.copy(new_raster=tfm_img)
            this_mask = mask_geo.reproject(this_img)
            this_mask.img[this_img.img == 0] = 0

        # only use grid points that are within this image
        isin = np.array([not this_img.outside_image(pt, index=False) for pt in all_pts])
        these_pts = all_pts[isin]
        these_names = pt_names[isin]

        these_meas = pd.DataFrame()
        these_meas['id'] = these_names
        these_meas['x'] = these_pts[:, 0]
        these_meas['y'] = these_pts[:, 1]

        for ii, row in these_meas.iterrows():
            pt = np.array(this_ref.xy2ij((row.x, row.y))).astype(int)
            these_meas.loc[ii, 'search_j'] = pt[1]
            these_meas.loc[ii, 'search_i'] = pt[0]

            # TODO: choose src, dstwin size based on image size, "goodness" of transform
            match_pt, this_z_corr, this_peak_corr = matching.do_match(this_img.img, this_ref.img,
                                                                      this_mask.img, pt, srcwin=20, dstwin=25)
            these_meas.loc[ii, 'match_j'] = match_pt[0]
            these_meas.loc[ii, 'match_i'] = match_pt[1]

            if is_geo:
                ix, iy = this_img.ij2xy(match_pt)
            else:
                # _x, _y = this_img.ij2xy(match_pt)
                j, i = this_model(match_pt)[0]
                ix, iy = orig_img.ij2xy((i, j))

            # will need to use xyz2im to transform this point to the original image
            these_meas.loc[ii, 'img_x'] = ix
            these_meas.loc[ii, 'img_y'] = iy
            if np.isfinite(ix) and np.isfinite(iy):
                these_meas.loc[ii, 'img_z'] = rel_dem.raster_points([(ix, iy)], nsize=3, mode='linear')
            else:
                these_meas.loc[ii, 'img_z'] = np.nan
            these_meas.loc[ii, 'z_corr'] = this_z_corr
            these_meas.loc[ii, 'pk_corr'] = this_peak_corr

        if these_meas.shape[0] > 0:
            these_meas['dj'] = these_meas['match_j'] - these_meas['search_j']
            these_meas['di'] = these_meas['match_i'] - these_meas['search_i']

            these_meas.dropna(how='any', inplace=True)
            these_meas.index = range(these_meas.shape[0])
            # valid = np.logical_and.reduce([these_meas['match_j'] > 0,
            #                                these_meas['match_j'] < this_img.npix_x,
            #                                these_meas['match_i'] > 0,
            #                                these_meas['match_i'] < this_img.npix_y])
            # these_meas = these_meas.loc[valid]

            with open('TmpAutoMeasures.txt', 'w') as f:
                for ii, row in these_meas.iterrows():
                    print('{} {} {}'.format(row.img_x, row.img_y, row.img_z), file=f)

            subprocess.Popen(['mm3d', 'XYZ2Im', 'Ori-{}/Orientation-{}.xml'.format(ori, im),
                              'TmpAutoMeasures.txt', 'Auto-{}.txt'.format(im)]).wait()

            # now have to read the output
            impts = pd.read_csv('Auto-{}.txt'.format(im), sep=' ', names=['j', 'i'])
            these_meas['j'] = impts['j']
            these_meas['i'] = impts['i']

            this_im_mes = E.MesureAppuiFlottant1Im(E.NameIm(im))
            raw_img = imread(im)
            maxi, maxj = raw_img.shape

            for ii, row in these_meas.iterrows():
                if np.logical_and.reduce([np.isfinite(row.j), np.isfinite(row.i),
                                          0 < row.j < maxj, 0 < row.i < maxi]):
                    this_mes = E.OneMesureAF1I(E.NamePt(these_meas.loc[ii, 'id']),
                                               E.PtIm('{} {}'.format(row.j, row.i)))
                    this_im_mes.append(this_mes)

            these_meas.to_csv('TmpMeasures{}/{}.csv'.format(subscript, im), index=False)
            MesureSet.append(this_im_mes)

            os.remove('Auto-{}.txt'.format(im))
            os.remove('TmpAutoMeasures.txt')

    tree = etree.ElementTree(MesureSet)
    tree.write(os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript)),
               pretty_print=True, xml_declaration=True, encoding="utf-8")

    for ii in range(2):
        gcps = micmac.bascule(gcps, out_dir, match_pattern, subscript, ori)

        valid = np.abs(gcps.residual - gcps.residual.median()) < nmad(gcps.residual)
        gcps = gcps.loc[valid]

        micmac.save_gcps(gcps, out_dir, utm_str, subscript)
