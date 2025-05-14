"""
spymicmac.register is a collection of tools for registering images and finding GCPs.
"""
import os
from pathlib import Path
import re
import subprocess
import shutil
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
import geoutils as gu
from glob import glob
import pandas as pd
from rtree import index
from shapely.ops import unary_union
from shapely.geometry.point import Point
from skimage.measure import ransac
from skimage.transform import AffineTransform, warp
from spymicmac import data, image, matching, micmac, orientation


def nmad(values, nfact=1.4826):
    """
    Calculate the normalized median absolute deviation (NMAD) of an array.

    :param array-like values: input data
    :param float nfact: normalization factor for the data; default is 1.4826

    :returns nmad: (normalized) median absolute deviation of data.
    """
    m = np.nanmedian(values)
    return nfact * np.nanmedian(np.abs(values - m))


def _sliding_window_filter(img_shape, pts_df, winsize, stepsize=None, mindist=2000, how='residual', is_ascending=True):
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
            # samp_ = samp_[samp_.z_corr >= samp_.z_corr.quantile(0.5)]
            # make sure we get the best residual
            samp_.sort_values(how, ascending=is_ascending, inplace=True)
            if len(_out_inds) == 0:
                best_ind = samp_.index[0]
                best_pt = Point(samp_.loc[best_ind, ['orig_j', 'orig_i']].values.astype(float))

                _out_inds.append(best_ind)
                _out_pts.append(best_pt)
            else:
                for _ind, _row in samp_.iterrows():
                    this_pt = Point(_row[['orig_j', 'orig_i']].values.astype(float))
                    this_min_dist = np.array([this_pt.distance(pt) for pt in _out_pts]).min()
                    if this_min_dist > mindist:
                        _out_inds.append(_ind)
                        _out_pts.append(this_pt)

    return np.array(_out_inds)


def _get_imlist(im_subset, dirname='.', strip_text=None):
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


def _get_utm_str(epsg):
    """
    Return a UTM zone name based on a 5-digit EPSG code.

    Examples:
        - get_utm_str(32608) -> 'UTM Zone 8N'
        - get_utm_str(32708) -> 'UTM Zone 8S'

    :param int|str epsg: a str or int representing an EPSG Code
    :return: **utm_str** (*str*) -- the UTM string representation
    """
    epsg_str = str(epsg)
    hemi_dict = {'6': 'N', '7': 'S'}

    if epsg is None:
        return 'not utm'

    if epsg_str[:2] == '32':
        return epsg_str[-2:] + hemi_dict[epsg_str[2]]
    else:
        return 'not utm'


def _split_cps_gcps(gcps, frac):
    cps = gcps.sample(frac=frac)
    gcps.drop(cps.index, inplace=True)
    return gcps.sort_index(), cps.sort_index()


def warp_image(model, ref, img):
    """
    Given a transformation model between two coordinate systems, warp an image to a reference image

    :param GeometricTransform model: the transformation model between the coordinate systems
    :param Raster ref: the reference Raster
    :param Raster img: the Raster to be transformed
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
    tfm_img = warp(img.data, this_model, output_shape=ref.shape, preserve_range=True)

    return tfm_img, this_model, inliers


def _to_tfw(gt):
    # convert from gdal GeoTransform to TFW format
    return [gt[1], gt[2], gt[4], gt[5], gt[0], gt[3]]


def _get_footprint_overlap(fprints):
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


def _get_footprint_mask(shpfile, rast, filelist, fprint_out=False):
    """
    Return a footprint mask for an image.

    :param str|GeoDataFrame shpfile: a filename or a GeoDataFrame representation of the footprints.
    :param Raster rast: the image to create a mask for.
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

    fprint = _get_footprint_overlap(fp.to_crs(rast.crs))

    tmp_gdf = gpd.GeoDataFrame(columns=['geometry'])
    tmp_gdf.loc[0, 'geometry'] = fprint
    tmp_gdf.crs = rast.crs
    tmp_gdf.to_file('tmp_fprint.shp')

    maskout = gu.Vector('tmp_fprint.shp').create_mask(rast)

    for f in glob('tmp_fprint.*'):
        os.remove(f)
    if fprint_out:
        return maskout, fprint
    else:
        return maskout


def _get_mask(footprints, img, imlist, landmask=None, glacmask=None):
    """
    Create a mask for an image from different sources.

    :param GeoDataFrame footprints: vector data representing image footprints
    :param Raster img: the Raster to create a mask for
    :param array-like imlist: a list of image names
    :param str landmask: path to file of land outlines (i.e., an inclusion mask)
    :param str glacmask: path to file of glacier outlines (i.e., an exclusion mask)

    :returns:
        - **mask** (*array-like*) -- the mask
        - **fmask** (*Raster*) -- the georeferenced footprint mask
        - **img** (*Raster*) -- the input Raster, cropped to a buffer around the image footprints
    """
    fmask, fprint = _get_footprint_mask(footprints, img, imlist, fprint_out=True)

    xmin, ymin, xmax, ymax = fprint.bounds
    buff_dist = 0.05 * min(abs(xmax - xmin), abs(ymax - ymin))

    img.crop(fprint.buffer(buff_dist).bounds, mode='match_pixel', inplace=True)
    fmask.crop(fprint.buffer(buff_dist).bounds, mode='match_pixel', inplace=True)

    mask = fmask.copy(new_array=np.ones(img.shape, dtype=np.uint8))
    mask.astype(np.uint8, inplace=True)
    mask.data[~img.data.mask] = 255
    mask.data[mask.data == 1] = 0

    if landmask is not None:
        lmask = gu.Vector(landmask).create_mask(img)
        mask[~lmask] = 0
    if glacmask is not None:
        gmask = gu.Vector(glacmask).create_mask(img)
        mask[gmask] = 0

    mask[~fmask] = 0

    return mask, fmask, img


def _search_size(imshape):
    min_dst = 100
    dst = min(400, imshape[0] / 10, imshape[1] / 10)
    return int(max(min_dst, dst))


def _get_last_malt(dirmec):
    dem_list = sorted(glob('Z_Num*STD-MALT.tif', root_dir=dirmec))
    level = int(re.findall(r'\d+', dem_list[-1].split('_')[1])[0])
    zoomf = int(re.findall(r'\d+', dem_list[-1].split('_')[2])[0])

    return f'Z_Num{level}_DeZoom{zoomf}_STD-MALT.tif'


def _load_auto_mask(dirmec):
    last_dem = _get_last_malt(dirmec)
    num = int(last_dem.split('_')[1].strip('Num')) - 1

    shutil.copy(Path(dirmec, last_dem.replace('.tif', '.tfw')),
                Path(dirmec, f"AutoMask_STD-MALT_Num_{num}.tfw"))

    return gu.Raster(Path(dirmec, f"AutoMask_STD-MALT_Num_{num}.tif")) == 1


def _read_gcps(fn_gcp, ref_img):
    gcps = gpd.read_file(fn_gcp)
    gcps = gcps.rename(columns={'z': 'elevation', 'name': 'id'})

    if isinstance(gcps, gpd.GeoDataFrame):
        gcps = gcps.to_crs(ref_img.crs)
    else:
        gcps = gpd.GeoDataFrame(gcps, geometry=gpd.points_from_xy(gcps.x, gcps.y, crs=ref_img.crs))

    gcps['search_i'], gcps['search_j'] = ref_img.xy2ij(gcps.geometry.x, gcps.geometry.y)

    return gcps


def register_relative(dirmec, fn_dem, fn_ref=None, fn_ortho=None, glacmask=None, landmask=None, footprints=None,
                      im_subset=None, block_num=None, subscript=None, ori='Relative', ortho_res=8.,
                      imgsource='DECLASSII', density=200, out_dir=None, allfree=True, useortho=False, max_iter=5,
                      use_cps=False, cp_frac=0.2, use_orb=False, fn_gcps=None):
    """
    Register a relative DEM or orthoimage to a reference DEM and/or orthorectified image.

    :param str dirmec: the name of the MEC directory to read the relative DEM from (e.g., MEC-Relative)
    :param str fn_dem: path to reference DEM
    :param str fn_ref: path to reference orthorectified image (optional)
    :param str fn_ortho: path to relative orthoimage (optional)
    :param str glacmask: path to file of glacier outlines (i.e., an exclusion mask)
    :param str landmask: path to file of land outlines (i.e., an inclusion mask)
    :param str footprints: path to shapefile of image outlines. If not set, will download from USGS.
    :param str im_subset: subset of raw images to work with
    :param str block_num: block number to use if processing multiple image blocks
    :param str subscript: optional subscript to use for output filenames (default: None)
    :param str ori: name of orientation directory (after Ori-) (default: Relative)
    :param float ortho_res: approx. ground sampling distance (pixel resolution) of ortho image (default: 8 m)
    :param str imgsource: USGS dataset name for images (default: DECLASSII)
    :param int density: pixel spacing to look for GCPs (default: 200)
    :param str out_dir: output directory to save auto GCP files to (default: auto_gcps)
    :param bool allfree: run Campari setting all parameters free (default: True)
    :param bool useortho: use the orthomosaic in Ortho-{dirmec} rather than the DEM (default: False). If fn_ortho is
        set, uses that file instead.
    :param int max_iter: the maximum number of Campari iterations to run. (default: 5)
    :param bool use_cps: split the GCPs into GCPs and CPs, to quantify the uncertainty of the camera model (default: False)
    :param float cp_frac: the fraction of GCPs to use when splitting into GCPs and CPs (default: 0.2)
    :param bool use_orb: use skimage.feature.ORB to identify GCP locations in the reference image
        (default: use regular grid for matching)
    :param str fn_gcps: (optional) shapefile or CSV of GCP coordinates to use. Column names should be [(name | id),
        (z | elevation), x, y]. If CSV is used, x,y should have the same CRS as the reference image.
    """
    print('start.')

    if fn_ortho is not None or useortho:
        assert fn_ref is not None, "If using ortho image, fn_ref must be set."

    if out_dir is None:
        out_dir = 'auto_gcps'

    os.makedirs(out_dir, exist_ok=True)

    if subscript is not None:
        subscript = '_' + subscript
    elif block_num is not None:
        subscript = f"_block{block_num}"
    else:
        subscript = ''

    # use the supplied dirmec to get the ortho directory
    ort_dir = 'Ortho-' + dirmec

    imlist, match_pattern = _get_imlist(im_subset)

    # if we're using the ortho image, load the reference ortho
    if useortho:
        ref_img = gu.Raster(fn_ref)
    # otherwise, load the reference dem
    else:
        ref_img = gu.Raster(fn_dem)

    utm_str = _get_utm_str(ref_img.crs.to_epsg())

    rel_dem = gu.Raster(os.path.join(dirmec, _get_last_malt(dirmec)))
    rel_mask = _load_auto_mask(dirmec)
    rel_dem[~rel_mask] = np.nan

    if useortho:
        if fn_ortho is None:
            fn_reg = os.path.join(ort_dir, 'Orthophotomosaic.tif')
        else:
            fn_reg = fn_ortho
        reg_img = gu.Raster(fn_reg)
        reg_img.set_nodata(0)
    else:
        fn_reg = os.path.join(dirmec, _get_last_malt(dirmec))
        reg_img = gu.Raster(fn_reg)
        reg_img[~rel_mask] = np.nan
    print(f"Loaded relative image {fn_reg}.")

    # set the fill value for the rough_tfm image (0 or np.nan)
    if np.issubdtype(reg_img.data.dtype, np.floating):
        tfm_fill = np.nan
    else:
        tfm_fill = 0

    if footprints is None:
        if os.path.isfile('Footprints.gpkg'):
            print('Using existing Footprints.gpkg file in current directory.')
            footprints = gpd.read_file('Footprints.gpkg')
        else:
            clean_imlist = [im.split('OIS-Reech_')[-1].split('.tif')[0] for im in imlist]
            print('Attempting to get image footprints from USGS EarthExplorer.')
            footprints = data.get_usgs_footprints(clean_imlist, dataset=imgsource)

            print('Saving footprints to current directory.')
            footprints.to_file('Footprints.gpkg')
    else:
        footprints = gpd.read_file(footprints)

    mask_full, _, ref_img = _get_mask(footprints, ref_img, imlist, landmask, glacmask)

    Minit, _, centers = orientation.transform_centers(reg_img, ref_img, imlist, footprints, f"Ori-{ori}")
    rough_tfm = warp(reg_img.data, Minit, output_shape=ref_img.shape, preserve_range=True, cval=tfm_fill)

    rough_spacing = max(1000, np.round(max(ref_img.shape) / 20 / 1000) * 1000)

    rough_gcps = matching.find_matches(rough_tfm, ref_img, mask_full.data.data, initM=Minit,
                                       spacing=int(rough_spacing), dstwin=int(rough_spacing))

    try:
        model, inliers = ransac((rough_gcps[['search_j', 'search_i']].values,
                                 rough_gcps[['orig_j', 'orig_i']].values), AffineTransform,
                                min_samples=6, residual_threshold=20, max_trials=5000)
        if model is None or np.count_nonzero(inliers) < 6:
            raise ValueError()

        rough_tfm = warp(reg_img.data, model, output_shape=ref_img.shape, preserve_range=True, cval=tfm_fill)

    except ValueError as e:
        print('Unable to refine transformation with rough GCPs. Using transform estimated from footprints.')
        model = Minit

    rough_geo = ref_img.copy(new_array=rough_tfm)
    rough_geo.save(f"Register{subscript}_rough_geo.tif")

    fig, axs = plt.subplots(1, 2, figsize=(7, 5))
    axs[0].imshow(rough_tfm[::10, ::10], extent=[0, rough_tfm.shape[1], rough_tfm.shape[0], 0], cmap='gray',
                 vmin=np.nanpercentile(rough_tfm, 0.5), vmax=np.nanpercentile(rough_tfm, 99.5))
    axs[1].imshow(ref_img.data[::10, ::10], extent=[0, ref_img.shape[1], ref_img.shape[0], 0], cmap='gray')

    for ax in axs:
        ax.set_xticks([])
        ax.set_yticks([])

    fig.savefig(f"initial_transformation{subscript}.png", dpi=200, bbox_inches='tight')
    plt.close(fig)

    if fn_gcps is not None:
        gcps = _read_gcps(fn_gcps, ref_img)
        gcps = matching.find_matches(rough_tfm, ref_img, mask_full.data.data, points=gcps, initM=model,
                                     spacing=density, dstwin=_search_size(rough_tfm.shape))
    else:
        if use_orb:
            ref_hp = image.highpass_filter(ref_img.data)
            ref_hp[np.isnan(ref_hp)] = 0

            keypoints = matching.get_dense_keypoints(ref_hp, mask=None, npix=int(density/2), use_skimage=True,
                                                     detector_kwargs={'n_keypoints': 5})

            gcps = pd.DataFrame(data=keypoints, columns=['search_i', 'search_j'])
            gcps = matching.find_matches(rough_tfm, ref_img, mask_full.data.data, points=gcps, initM=model,
                                         dstwin=_search_size(rough_tfm.shape))

        else:
            gcps = matching.find_matches(rough_tfm, ref_img, mask_full.data.data, initM=model,
                                         spacing=density, dstwin=_search_size(rough_tfm.shape))

    x, y = ref_img.ij2xy(gcps['search_i'], gcps['search_j'])
    gcps = gpd.GeoDataFrame(gcps, geometry=gpd.points_from_xy(x, y, crs=ref_img.crs))

    gcps = gcps.loc[mask_full.data.data[gcps.search_i.astype(int), gcps.search_j.astype(int)] == 255]

    gcps['rel_x'], gcps['rel_y'] = reg_img.ij2xy(gcps.orig_i, gcps.orig_j)

    gcps.dropna(inplace=True)
    print(f"{gcps.shape[0]} potential matches found")

    if useortho:
        print('loading dems')
        dem = gu.Raster(fn_dem)
    else:
        dem = ref_img

    if 'elevation' not in gcps.columns:
        gcps['elevation'] = dem.interp_points((gcps.geometry.x, gcps.geometry.y))
    gcps['el_rel'] = rel_dem.interp_points((gcps.rel_x, gcps.rel_y))

    # drop any gcps where we don't have a DEM value or a valid match
    if dem.nodata is not None:
        gcps.loc[np.abs(gcps.elevation - dem.nodata) < 1, 'elevation'] = np.nan
    gcps.dropna(inplace=True)
    print(f"{gcps.shape[0]} matches with valid elevations")

    # run ransac to find the matches between the transformed image and the master image make a coherent transformation
    # residual_threshold is 20 pixels to allow for some local distortions, but get rid of the big blunders
    gcps['offset'] = np.sqrt(gcps['dj'] ** 2 + gcps['di'] ** 2)
    thresh = gcps['offset'].median() + 4 * nmad(gcps['offset'])

    Mref, inliers_ref = ransac((gcps[['search_j', 'search_i']].values, gcps[['match_j', 'match_i']].values),
                               AffineTransform, min_samples=6, residual_threshold=thresh, max_trials=5000)
    gcps['aff_resid'] = Mref.residuals(gcps[['search_j', 'search_i']].values,
                                       gcps[['match_j', 'match_i']].values)

    # valid = np.abs(gcps.aff_resid - gcps.aff_resid.median()) < nmad(gcps.aff_resid)
    # gcps = gcps.loc[valid]
    gcps = gcps.loc[inliers_ref]

    # out = _sliding_window_filter([reg_img.shape[1], reg_img.shape[0]], gcps,
    #                              min(500, reg_img.shape[1] / 4, reg_img.shape[0] / 4),
    #                              mindist=500, how='pk_corr', is_ascending=True)
    # gcps = gcps.loc[out]

    print(f"{gcps.shape[0]} valid matches found after estimating transformation")

    gcps.index = range(gcps.shape[0])  # make sure index corresponds to row we're writing out
    if 'id' not in gcps.columns:
        gcps['id'] = [f'GCP{ind}' for ind in range(gcps.shape[0])]

    gcps.to_file(Path(out_dir, f"AutoGCPs{subscript}.shp"))

    print('writing AutoGCPs.txt')
    micmac.write_auto_gcps(gcps, subscript, out_dir, utm_str)

    print('converting AutoGCPs.txt to AutoGCPs.xml')
    subprocess.Popen(['mm3d', 'GCPConvert', 'AppInFile',
                      Path(out_dir, f"AutoGCPs{subscript}.txt")]).wait()

    print('writing AutoMeasures.txt')
    micmac.write_auto_mesures(gcps, subscript, out_dir)

    print('running get_autogcp_locations to get rough image locations for each point')
    micmac.get_autogcp_locations(f"Ori-{ori}", os.path.join(out_dir, f"AutoMeasures{subscript}.txt"), imlist)

    # print('searching for points in orthorectified images')
    print('finding image measures')
    micmac.write_image_mesures(imlist, gcps, out_dir, subscript, ort_dir=ort_dir)

    print('running mm3d GCPBascule to estimate terrain errors')
    gcps = micmac.bascule(gcps, out_dir, match_pattern, subscript, ori)
    gcps['res_dist'] = np.sqrt(gcps.xres ** 2 + gcps.yres ** 2)

    gcps = gcps.loc[np.abs(gcps.res_dist - gcps.res_dist.median()) < 3 * nmad(gcps.res_dist)]
    # gcps = gcps[np.logical_and(np.abs(gcps.xres - gcps.xres.median()) < 2 * nmad(gcps.xres),
    #                            np.abs(gcps.yres - gcps.yres.median()) < 2 * nmad(gcps.yres))]

    micmac.save_gcps(gcps, out_dir, utm_str, subscript)
    gcps = micmac.bascule(gcps, out_dir, match_pattern, subscript, ori)
    gcps['res_dist'] = np.sqrt(gcps.xres ** 2 + gcps.yres ** 2)
    gcps = gcps.loc[np.abs(gcps.res_dist - gcps.res_dist.median()) < 3 * nmad(gcps.res_dist)]
    # gcps = gcps[np.logical_and(np.abs(gcps.xres - gcps.xres.median()) < 2 * nmad(gcps.xres),
    #                            np.abs(gcps.yres - gcps.yres.median()) < 2 * nmad(gcps.yres))]

    micmac.save_gcps(gcps, out_dir, utm_str, subscript)

    # if we're getting cps, split them now
    if use_cps:
        gcps, cps = _split_cps_gcps(gcps, cp_frac)

        micmac.save_gcps(gcps, out_dir, utm_str, subscript)

        micmac.write_auto_gcps(cps, subscript, out_dir, utm_str, outname='AutoCPs')
        subprocess.Popen(['mm3d', 'GCPConvert', 'AppInFile',
                          Path(out_dir, f"AutoCPs{subscript}.txt")]).wait()

        micmac.write_auto_mesures(cps, subscript, out_dir, outname='AutoCPMeasures')
        micmac.get_autogcp_locations(f"Ori-{ori}", Path(out_dir, f"AutoCPMeasures{subscript}.txt"), imlist)
        micmac.write_image_mesures(imlist, cps, out_dir, subscript, ort_dir=ort_dir, outname='AutoCPMeasures')

    # now, iterate campari to refine the orientation
    gcps = micmac.iterate_campari(gcps, out_dir, match_pattern, subscript, ref_img.res[0], ortho_res,
                                  rel_ori=ori, allfree=allfree, max_iter=max_iter)

    if use_cps:
        cp_resids = micmac.checkpoints(match_pattern, f"Ori-TerrainFinal{subscript}",
                                       fn_cp=Path(out_dir, f"AutoCPs{subscript}.xml"),
                                       fn_meas=Path(out_dir, f"AutoCPMeasures{subscript}-S2D.xml"),
                                       fn_resids=Path(out_dir, f"AutoCPs{subscript}", ret_df=True))

        # merge cps and cp_resids
        cps.drop(['xres', 'yres', 'zres', 'residual', 'res_dist'], axis=1, inplace=True)

        cps = cps.merge(cp_resids, left_on='id', right_on='id')
        cps['res_dist'] = np.sqrt(cps.xres ** 2 + cps.yres ** 2)
        cps['residual'] = np.sqrt(cps.xres ** 2 + cps.yres ** 2 + cps.zres ** 2)

        cps.to_file(Path(out_dir, f"AutoCPs{subscript}.shp"))

    # final write of gcps to disk.
    gcps.to_file(os.path.join(out_dir, f"AutoGCPs{subscript}.shp"))

    fig1, ax1 = plt.subplots(1, 1, figsize=(7, 5))
    ax1.imshow(reg_img[::5, ::5], cmap='gray', extent=(0, reg_img.shape[1], reg_img.shape[0], 0))
    ax1.plot(gcps.orig_j, gcps.orig_i, 'r+')
    ax1.quiver(gcps.orig_j, gcps.orig_i, gcps.camp_xres, gcps.camp_yres, color='r')
    if use_cps:
        ax1.plot(cps.orig_j, cps.orig_i, 'bs')
        ax1.quiver(cps.orig_j, cps.orig_i, cps.xres, cps.yres, color='b')

    fig1.savefig(os.path.join(out_dir, f"relative_gcps{subscript}.png"), bbox_inches='tight', dpi=200)
    plt.close(fig1)

    fig2, ax2 = plt.subplots(1, 1, figsize=(7, 5))
    xmin, ymin, xmax, ymax = ref_img.bounds
    ax2.imshow(ref_img.data[::10, ::10], cmap='gray', extent=(xmin, xmax, ymin, ymax))
    ax2.plot(gcps.geometry.x, gcps.geometry.y, 'r+')
    if use_cps:
        ax2.plot(cps.geometry.x, cps.geometry.y, 'bs')

    fig2.savefig(os.path.join(out_dir, f"world_gcps{subscript}.png"), bbox_inches='tight', dpi=200)
    plt.close(fig2)

    print('cleaning up.')

    for txtfile in glob('Auto-OIS*.tif.txt') + \
                   glob('NoDist-OIS*.tif.txt'):
        os.remove(txtfile)
    print('end.')
    # embed()
