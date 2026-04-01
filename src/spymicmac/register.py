"""
spymicmac.register is a collection of tools for registering images and finding GCPs.
"""
import os
import warnings
from pathlib import Path
import re
import subprocess
import shutil
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
import geoutils as gu
import xdem
from glob import glob
import pandas as pd
from rtree import index
from scipy import stats
from shapely.ops import unary_union
from shapely.geometry import Point, Polygon, LineString
from skimage.filters import gaussian
from skimage.measure import ransac
from skimage.transform import AffineTransform, warp
from . import data, image, matching, micmac, orientation
from pyogrio.errors import DataSourceError
from numpy.typing import NDArray
from typing import Union
from skimage.transform._geometric import _GeometricTransform


def nmad(values: NDArray, nfact: Union[float, int] = 1.4826) -> NDArray:
    """
    Calculate the normalized median absolute deviation (NMAD) of an array.

    :param values: input data
    :param nfact: normalization factor for the data
    :returns: **nmad** -- (normalized) median absolute deviation of data.
    """
    m = np.nanmedian(values)
    return nfact * np.nanmedian(np.abs(values - m))


def rms(values: NDArray) -> NDArray:
    """
    Calculate the root mean square (rms) of an array of values.

    :param values: the values to use
    :returns: rms of the input values
    """
    return np.sqrt(np.nanmean(np.asarray(values) ** 2))


def _sliding_window_filter(img_shape: NDArray, pts_df: pd.DataFrame, winsize: int, stepsize: Union[int, None] = None,
                           mindist: int = 2000, how: str = 'residual', is_ascending: bool = True) -> NDArray:
    """
    Given a DataFrame of indices representing points, use a sliding window filter to keep only the 'best' points within
    a given window and separation distance.

    :param img_shape: The shape of the image to filter indices from.
    :param pts_df: A DataFrame that contains the pixel locations of points (column: orig_j, row: orig_i)
    :param winsize: the size of the window to use for the filter.
    :param stepsize: how large of a step size to use for the sliding window (default: winsize/2)
    :param mindist: the minimum distance (in pixels) to use between points
    :param how: how to sort pts_df to determine the 'best' points to keep
    :param is_ascending: whether the column named in 'how' should be sorted ascending or descending
    :return: **out_inds** -- an array with the filtered indices
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


def _get_imlist(im_subset: Union[list, str, None], globstr: str = 'OIS*.tif') -> tuple[list, str]:
    """
    Given either a list of filenames or a regex match pattern, return a list of filenames and a pattern to provide
    to MicMac.

    :param im_subset: a list of filenames or a match pattern (e.g., ['OIS.*tif']) representing filenames. If None,
        returns all images matching "OIS.*tif"
    :param globstr: the search string to use to find images in the current directory
    :return:
        - **imlist** -- the list of filenames matching the provided pattern.
        - **match_pattern** -- the match pattern to be provided to MicMac.
    """
    if im_subset is None:
        imlist = glob(globstr)
        if '*.' in globstr:
            match_pattern = globstr.replace('*.', '.*')
        else:
            match_pattern = globstr
    else:
        if len(im_subset) > 1:
            imlist = im_subset
            match_pattern = micmac.get_match_pattern(imlist)
            # match_pattern = '|'.join(imlist)
        else:
            match_pattern = im_subset[0] + '.*tif'
            imlist = [f for f in glob(globstr) if re.search(match_pattern, f)]

    imlist.sort()

    return imlist, match_pattern


def _get_utm_str(epsg: Union[int, str]) -> str:
    """
    Return a UTM zone name based on a 5-digit EPSG code.

    Examples:
        - get_utm_str(32608) -> 'UTM Zone 8N'
        - get_utm_str(32708) -> 'UTM Zone 8S'

    :param epsg: a str or int representing an EPSG Code
    :return: **utm_str** -- the UTM string representation
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


def warp_image(model: _GeometricTransform, ref: gu.Raster, img: gu.Raster) -> tuple[NDArray, AffineTransform, NDArray]:
    """
    Given a transformation model between two coordinate systems, warp an image to a reference image

    :param model: the transformation model between the coordinate systems
    :param ref: the reference Raster
    :param img: the Raster to be transformed
    :return:
        - **tfm_img** -- the input image transformed to the same extent as the reference image
        - **this_model** -- the estimated Affine Transformation between the two images
        - **inliers** -- an array of the inliers returned by skimage.measure.ransac
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

def _tfm_from_corners(rel_img, ref_img, model):
    rel_corners = np.array([(0, 0),
                            (rel_img.shape[1], 0),
                            (rel_img.shape[1], rel_img.shape[0]),
                            (0, rel_img.shape[0])])
    ref_corners = model.inverse(rel_corners)

    rel_xy = np.array(rel_img.ij2xy(rel_corners[:,1], rel_corners[:,0])).T
    ref_xy = np.array(ref_img.ij2xy(ref_corners[:,1], ref_corners[:,0])).T

    new_model = AffineTransform()
    new_model.estimate(rel_xy, ref_xy)

    return new_model

def _corners_from_tfm(rel_img, ref_img, model):
    rel_corners = np.array([(0, 0),
                            (rel_img.shape[1], 0),
                            (rel_img.shape[1], rel_img.shape[0]),
                            (0, rel_img.shape[0])])
    rel_xy = np.array(rel_img.ij2xy(rel_corners[:,1], rel_corners[:,0])).T

    ref_xy = model(rel_xy)
    ref_corners = np.array(ref_img.xy2ij(ref_xy[:,0], ref_xy[:,1])).T[:, ::-1]

    new_model = AffineTransform()
    new_model.estimate(ref_corners, rel_corners)

    return new_model


def _reproject_footprint(poly: Polygon, img: gu.Raster, model: AffineTransform, invert: bool = True):
    x, y = poly.exterior.coords.xy

    src = np.array(img.xy2ij(x, y)).T[:, ::-1]

    if invert:
        return model.inverse(src)
    else:
        return model(src)


def _reproject_footprints(fprints: gpd.GeoDataFrame, rel_img: gu.Raster,
                          ref_img: gu.Raster, model: AffineTransform) -> gpd.GeoDataFrame:

    geoms = []
    for poly in fprints['geometry']:
        coords = _reproject_footprint(poly, rel_img, model, invert=True)
        x, y = ref_img.ij2xy(coords[:, 1], coords[:, 0])
        geoms.append(Polygon(zip(x, y)))

    ids = [im.split('OIS-Reech_')[-1].split('.tif')[0] for im in fprints.index]

    return gpd.GeoDataFrame(data={'filename': fprints.index, 'ID': ids}, geometry=geoms, crs=ref_img.crs)


def _get_footprint_overlap(fprints: gpd.GeoDataFrame) -> Polygon:
    """
    Return the area where image footprints overlap.

    :param fprints: a GeoDataFrame of image footprints
    :return: **intersection** -- the overlapping area (unary union) of the images.
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


def _get_footprint_mask(shpfile: Union[gpd.GeoDataFrame, str], rast: gu.Raster,
                        filelist: list, fprint_out: bool = False) -> Union[gu.Mask, tuple]:
    """
    Return a footprint mask for an image.

    :param shpfile: a filename or a GeoDataFrame representation of the footprints.
    :param rast: the image to create a mask for.
    :param filelist: a list of image names to use to create the footprint mask.
    :param fprint_out: return the polygon representation of the footprint mask.
    :return:
        - **mask** -- the footprint mask
        - **fprint** -- the footprint polygon, if requested.
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
    tmp_gdf.to_file('tmp_fprint.gpkg')

    maskout = gu.Vector('tmp_fprint.gpkg').create_mask(rast)

    for f in glob('tmp_fprint.*'):
        os.remove(f)
    if fprint_out:
        return maskout, fprint
    else:
        return maskout


def _get_mask(footprints: gpd.GeoDataFrame, img: gu.Raster, imlist: list,
              landmask: Union[str, Path, None] = None,
              glacmask: Union[str, Path, None] = None) -> tuple[gu.Mask, gu.Raster, gu.Raster]:
    """
    Create a mask for an image from different sources.

    :param footprints: vector data representing image footprints
    :param img: the Raster to create a mask for
    :param imlist: a list of image names
    :param landmask: path to file of land outlines (i.e., an inclusion mask)
    :param glacmask: path to file of glacier outlines (i.e., an exclusion mask)

    :returns:
        - **mask** -- the mask
        - **fmask** -- the georeferenced footprint mask
        - **img** -- the input Raster, cropped to a buffer around the image footprints
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
        lmask = _safe_rasterize(landmask, img, True)
        mask[~lmask] = 0
    if glacmask is not None:
        gmask = _safe_rasterize(glacmask, img, False)
        mask[gmask] = 0

    mask[~fmask] = 0

    return mask, fmask, img


def _safe_rasterize(shape, img, inclusive=True):
    try:
        _mask = _load_mask_from_file(shape, img)
    except ValueError as e:
        warnings.warn("No valid geometry objects found to rasterize. Be sure to double-check your input files.", UserWarning)
        _mask = img.copy(new_array=np.ones_like(img.data))
        if inclusive:
            _mask = _mask == 1
        else:
            _mask = _mask == 0

    return _mask


def _load_mask_from_file(fn_mask, img):
    try:
        _mask = gu.Vector(fn_mask).create_mask(img)
    except DataSourceError as e:
        _mask = gu.Raster(fn_mask)
        _mask = _mask.reproject(img) != 0

    return _mask


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
    gcps[['x', 'y', 'elevation']] = gcps[['x', 'y', 'elevation']].astype(float)

    if isinstance(gcps, gpd.GeoDataFrame):
        gcps = gcps.to_crs(ref_img.crs)
    else:
        gcps = gpd.GeoDataFrame(gcps, geometry=gpd.points_from_xy(gcps.x, gcps.y, crs=ref_img.crs))

    gcps['search_i'], gcps['search_j'] = ref_img.xy2ij(gcps.geometry.x, gcps.geometry.y)

    return gcps

def _resid_subplot(ax, gdf, cols, names, colors, styles):
    lower = np.floor(min(gdf[cols].min()))
    upper = np.ceil(max(gdf[cols].max()))

    bins = np.linspace(lower, upper, 100)

    for res, dd, cc, ls in zip(cols, names, colors, styles):
        density = stats.gaussian_kde(gdf[res])
        ax.plot(bins, density(bins), label=dd, linewidth=2, color=cc, linestyle=ls)

    ax.legend()
    ax.set_xlabel('residual (m)')
    ax.set_ylabel('frequency')

    return ax


def _plot_residuals(gdf, res='camp'):
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharey=True)

    if res != '':
        dcols = [f"{res}_{dd}" for dd in ['xres', 'yres', 'zres']]
        ecols = ['camp_xy', 'camp_dist']
    else:
        dcols = ['xres', 'yres', 'zres']
        ecols = ['res_xy', 'res_dist']

    _resid_subplot(axs[0], gdf, dcols, ['x', 'y', 'z'], ['0.2', '0.6', 'k'], ['-', '--', '-.'])
    _resid_subplot(axs[1], gdf, ecols, ['plani', '3d'], ['k', '0.5'], ['--', '-'])

    locs = [(0.05, 0.92), (0.05, 0.85), (0.05, 0.78)]
    xyz_stats = '\n'.join([f"rms$_{name}$: {rms(gdf[res]):.1f}" for res, name in zip(dcols, 'xyz')])

    axs[0].text(0.05, 0.95, xyz_stats, transform=axs[0].transAxes, fontsize=10, va='top', family='monospace',
                bbox={'boxstyle': 'round', 'facecolor': 'w', 'alpha': 0.8})


    stats_3d = '\n'.join([f"3d mean: {gdf[ecols[-1]].mean():.1f}",
                          f"3d med : {gdf[ecols[-1]].median():.1f}",
                          f"3d rms : {rms(gdf[ecols[-1]]):.1f}"])

    axs[1].text(0.05, 0.95, stats_3d, transform=axs[1].transAxes, fontsize=10, va='top', family='monospace',
                bbox={'boxstyle': 'round', 'facecolor': 'w', 'alpha': 0.8})


    axs[1].set_ylabel('')

    fig.suptitle(f"$n={len(gdf)}$ gcps", x=0.5, y=0.9, va='bottom')

    return fig


def _get_rectangle_rotation(box):

    coords = box.exterior.coords
    sides = [LineString([pt1, pt2]) for pt1, pt2 in zip(coords, coords[1:])]
    dists = [s.length for s in sides]

    long_ax = max(sides, key=lambda x: x.length)
    short_ax = min(sides, key=lambda x: x.length)

    if np.allclose(dists[0], long_ax.length):
        rect_coords = np.array([(0, 0),
                                (0, short_ax.length),
                                (long_ax.length, short_ax.length),
                                (long_ax.length, 0)])
    else:
        rect_coords = np.array([(0, 0),
                                (long_ax.length, 0),
                                (long_ax.length, short_ax.length),
                                (0, short_ax.length)])

    img_coords = np.array(box.exterior.coords)[:-1]

    rect_tfm = AffineTransform()
    rect_tfm.estimate(rect_coords, img_coords)

    return rect_tfm, long_ax, short_ax


def _chebyshev_grid(img, spacing, tfm, mask = None):

    if mask is None:
        rady, radx = img.shape

        nptsx = int(radx / spacing)
        nptsy = int(rady / spacing)

    else:
        dx, dy = img.res

        rect_tfm, long_ax, short_ax = _get_rectangle_rotation(mask)

        radx = long_ax.length
        rady = short_ax.length

        nptsx = int(radx / dx / spacing)
        nptsy = int(rady / dy / spacing)

    xx = np.cos(np.arange(0, nptsx) * np.pi / nptsx)
    yy = np.cos(np.arange(0, nptsy) * np.pi / nptsy)

    xx += 1
    yy += 1

    xx *= radx / 2
    yy *= rady / 2

    XX, YY = np.meshgrid(xx, yy)

    if mask is None:
        return tfm.inverse(np.array([XX.flatten(), YY.flatten()]).T)
    else:
        img_pts = rect_tfm(np.array([XX.flatten(), YY.flatten()]).T)
        yi, xi = img.xy2ij(img_pts[:, 0], img_pts[:, 1])
        return tfm.inverse(np.array([xi, yi]).T)


def _prepare_hillshades(dem, tfm_img, **hillshade_kwargs):

    # stretch tfm_img to approximate dem
    min_el = dem[np.isfinite(tfm_img)].min()
    max_el = dem[np.isfinite(tfm_img)].max()

    stretched = image.stretch_image(tfm_img, (0.005, 0.995), 1, 0, np.float32)

    dem_hs = xdem.DEM(dem).hillshade(**hillshade_kwargs)

    if 'z_factor' in hillshade_kwargs.keys():
        hillshade_kwargs.update({'z_factor': hillshade_kwargs['z_factor'] * (max_el - min_el)})
    else:
        hillshade_kwargs.update({'z_factor': 4 * (max_el - min_el)})

    tfm_hs = xdem.DEM(dem_hs.copy(new_array=stretched)).hillshade(**hillshade_kwargs)

    return dem_hs, tfm_hs


def register_relative(dirmec: str, fn_dem: Union[str, Path], fn_ref: Union[str, Path, None] = None,
                      fn_ortho: Union[str, Path, None] = None, glacmask: Union[str, Path, None] = None,
                      landmask: Union[str, Path, None] = None, footprints: Union[str, Path, None] = None,
                      globstr: str = 'OIS*.tif', im_subset: Union[str, None] = None, block_num: Union[str, None] = None,
                      subscript: Union[str, None] = None, ori: str ='Relative', ortho_res: Union[float, int] = 8.,
                      imgsource: str = 'DECLASSII', strategy: str = 'grid', density: int = 200,
                      out_dir: Union[str, Path, None] = None, allfree: bool = True, dir_homol: str = 'Homol',
                      useortho: bool = False, max_iter: int = 5, use_cps: bool = False, cp_frac: float = 0.2,
                      use_orb: bool = False, fn_gcps: Union[str, Path, None] = None,
                      use_blur: bool = False, use_highpass: bool = True,
                      use_hillshade: bool = False, hillshade_kwargs: dict = {}) -> None:
    """
    Register a relative DEM or orthoimage to a reference DEM and/or orthorectified image.

    :param dirmec: the name of the MEC directory to read the relative DEM from (e.g., MEC-Relative)
    :param fn_dem: path to reference DEM
    :param fn_ref: path to reference orthorectified image (optional)
    :param fn_ortho: path to relative orthoimage (optional)
    :param glacmask: path to file of glacier outlines (i.e., an exclusion mask)
    :param landmask: path to file of land outlines (i.e., an inclusion mask)
    :param footprints: path to shapefile of image outlines. If not set, will look for Footprints.gpkg in the current
        directory. If this file does not exist, will attempt to download from USGS using imgsource.
    :param globstr: the search string to use to find images in the current directory
    :param im_subset: subset of raw images to work with
    :param block_num: block number to use if processing multiple image blocks
    :param subscript: optional subscript to use for output filenames
    :param ori: name of orientation directory (after Ori-)
    :param imgsource: USGS dataset name for images
    :param strategy: strategy for generating GCPs. Must be one of 'grid', 'random', 'chebyshev', or 'peaks'. Note that
        if 'random' is used, density is the approximate number of points, rather than the distance between grid points.
        If 'peaks' is used, peaks will be spaced by at least (density/2) pixels.
    :param density: pixel spacing to look for GCPs
    :param out_dir: output directory to save auto GCP files to
    :param allfree: run Campari setting all parameters free
    :param dir_homol: the name of the Homol directory to use for Campari (default: Homol)
    :param useortho: use the orthomosaic in Ortho-{dirmec} rather than the DEM. If fn_ortho is
        set, uses that file instead.
    :param max_iter: the maximum number of Campari iterations to run.
    :param use_cps: split the GCPs into GCPs and checkpoints (CPs), to quantify the uncertainty of the camera model
    :param cp_frac: the fraction of GCPs to use as CPs when splitting into GCPs and CPs
    :param use_orb: use skimage.feature.ORB to identify GCP locations in the reference image
        (default: use regular grid for matching)
    :param fn_gcps: (optional) shapefile or CSV of GCP coordinates to use. Column names should be [(name | id),
        (z | elevation), x, y]. If CSV is used, x,y should have the same CRS as the reference image.
    :param use_blur: use a gaussian blur on the relative image before matching
    :param use_highpass: match templates using a highpass filter
    :param use_hillshade: match templates using DEM hillshade rather than elevation
    :param hillshade_kwargs: kwargs to pass to xdem.DEM.hillshade()
    """
    assert strategy in ['grid', 'random', 'chebyshev', 'peaks'], \
        f"{strategy} must be one of [grid, random, chebyshev, peaks]"
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

    imlist, match_pattern = _get_imlist(im_subset, globstr)

    # if we're using the ortho image, load the reference ortho
    if useortho:
        ref_img = gu.Raster(fn_ref)
    # otherwise, load the reference dem
    else:
        ref_img = gu.Raster(fn_dem)

    utm_str = _get_utm_str(ref_img.crs.to_epsg())

    rel_dem = gu.Raster(Path(dirmec, _get_last_malt(dirmec)))
    rel_mask = _load_auto_mask(dirmec)
    rel_dem[~rel_mask] = np.nan

    if useortho:
        if fn_ortho is None:
            fn_reg = Path(ort_dir, 'Orthophotomosaic.tif')
        else:
            fn_reg = fn_ortho
        reg_img = gu.Raster(fn_reg)
        reg_img.set_nodata(0)
    else:
        fn_reg = Path(dirmec, _get_last_malt(dirmec))
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

    mask_full, _, ref_img_crop = _get_mask(footprints, ref_img, imlist, landmask, glacmask)

    Minit, _, centers = orientation.transform_centers(reg_img, ref_img_crop, imlist, footprints, f"Ori-{ori}")
    rough_tfm = warp(reg_img.data, Minit, output_shape=ref_img_crop.shape, preserve_range=True, cval=tfm_fill)
    rough_tfm = rough_tfm.astype(reg_img.data.dtype)

    rough_spacing = max(1000, np.round(max(ref_img.shape) / 20 / 1000) * 1000)

    rough_gcps = matching.find_matches(rough_tfm, ref_img_crop, mask_full.data.data, initM=Minit,
                                       spacing=int(rough_spacing), dstwin=int(rough_spacing))

    try:
        model, inliers = ransac((rough_gcps[['search_j', 'search_i']].values,
                                 rough_gcps[['orig_j', 'orig_i']].values), AffineTransform,
                                min_samples=6, residual_threshold=20, max_trials=5000)
        if model is None or np.count_nonzero(inliers) < 6:
            raise ValueError()

    except ValueError as e:
        print('Unable to refine transformation with rough GCPs. Using transform estimated from footprints.')
        model = Minit

    # re-do the mask and the footprints with the re-projected relative footprints
    micmac.drone_footprint(f"({'|'.join(imlist)})", ori)
    rel_footprints = gpd.read_file(f"{ori}_footprints.gpkg").set_index('filename')
    new_footprints = _reproject_footprints(rel_footprints, reg_img, ref_img_crop, model)

    # transformation from relative coords to absolute coords
    rel2abs = _tfm_from_corners(reg_img, ref_img_crop, model)

    # re-do the mask with the new footprints
    mask_full, _, ref_img = _get_mask(new_footprints, ref_img, imlist, landmask, glacmask)

    model = _corners_from_tfm(reg_img, ref_img, rel2abs)
    rough_tfm = warp(reg_img.data, model, output_shape=ref_img.shape, preserve_range=True, cval=tfm_fill)
    rough_tfm = rough_tfm.astype(reg_img.data.dtype)

    if use_blur:
        print("Smoothing relative image with a Gaussian blur.")
        rough_tfm = gaussian(rough_tfm, 2)

    if use_hillshade:
        print("Using DEM hillshades for matching.")
        ref_img, rough_tfm = _prepare_hillshades(ref_img, rough_tfm, **hillshade_kwargs)

    if not use_hillshade:
        rough_geo = ref_img.copy(new_array=rough_tfm)
    else:
        rough_geo = rough_tfm

    rough_geo.save(f"Register{subscript}_rough_geo.tif", co_opts={'BIGTIFF': 'YES'})

    np.savetxt(f"Register{subscript}_rough_tfm.csv", model.params, delimiter=',')

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
    else:
        if use_orb:
            ref_hp = image.highpass_filter(ref_img.data)
            ref_hp[np.isnan(ref_hp)] = 0

            keypoints = matching.get_dense_keypoints(ref_hp, mask=None, npix=int(density/2), use_skimage=True,
                                                     detector_kwargs={'n_keypoints': 5})

            gcps = pd.DataFrame(data=keypoints, columns=['search_i', 'search_j'])

        elif strategy == 'chebyshev':
            rel_mask = warp(mask_full.data.data, model.inverse,
                            output_shape=reg_img.shape,
                            preserve_range=True, cval=0).astype(np.uint8)

            rel_tfm = list(reg_img.transform)[:6]
            rel_tfm[0] *= 10
            rel_tfm[4] *= 10

            rel_mask = gu.Raster.from_array(rel_mask[::10, ::10],
                                            transform=tuple(rel_tfm[:6]),
                                            crs=None)

            mask_bounds = (rel_mask > 0).polygonize().ds.union_all().minimum_rotated_rectangle

            gridpts = _chebyshev_grid(reg_img, density, model, mask_bounds)
            gcps = pd.DataFrame(data=gridpts, columns=['search_j', 'search_i'])
            interior = np.logical_and.reduce([0 <= gcps.search_j, gcps.search_j <= ref_img.shape[1],
                                              0 <= gcps.search_i, gcps.search_i <= ref_img.shape[0]])
            gcps = gcps.loc[interior]

            masked = mask_full.data.data[gcps.search_i.astype(int), gcps.search_j.astype(int)] == 0
            gcps = gcps.loc[~masked]

        elif strategy == 'peaks':
            ref_hp = image.highpass_filter(ref_img.data)
            ref_hp[np.isnan(ref_hp)] = 0

            gridpts = matching.peaks(ref_hp, int(density / 2), mask_full.data.data)
            gcps = pd.DataFrame(data=gridpts, columns=['search_i', 'search_j'])

        else:
            gcps = None

    gcps = matching.find_matches(rough_tfm, ref_img, mask_full.data.data, points=gcps, initM=model, strategy=strategy,
                                 spacing=density, srcwin=100, dstwin=_search_size(rough_tfm.shape),
                                 use_highpass=use_highpass)

    x, y = ref_img.ij2xy(gcps['search_i'], gcps['search_j'])
    gcps = gpd.GeoDataFrame(gcps, geometry=gpd.points_from_xy(x, y, crs=ref_img.crs))

    gcps = gcps.loc[mask_full.data.data[gcps.search_i.astype(int), gcps.search_j.astype(int)] == 255]

    gcps['rel_x'], gcps['rel_y'] = reg_img.ij2xy(gcps.orig_i, gcps.orig_j)

    gcps.dropna(inplace=True)
    print(f"{gcps.shape[0]} potential matches found")

    if useortho or use_hillshade:
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
    thresh = np.ceil(min(20, gcps['offset'].median() + 4 * nmad(gcps['offset'])))

    models = []
    inliers = []

    for ii in range(20):
        mod, inl = ransac((gcps[['search_j', 'search_i']].values, gcps[['match_j', 'match_i']].values),
                          AffineTransform, min_samples=6, residual_threshold=thresh, max_trials=5000)
        models.append(mod)
        inliers.append(inl)

    num_inliers = [np.count_nonzero(inl) for inl in inliers]
    best_ind = np.argmax(num_inliers)

    Mref = models[best_ind]
    inliers_ref = inliers[best_ind]

    if num_inliers[best_ind] < 10:
        raise ValueError("Unable to estimate valid transformation using matches found.")

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
    micmac.get_autogcp_locations(f"Ori-{ori}", Path(out_dir, f"AutoMeasures{subscript}.txt"), imlist)

    # print('searching for points in orthorectified images')
    print('finding image measures')
    if Path(ort_dir).exists():
        micmac.write_image_mesures(imlist, gcps, out_dir, subscript, ort_dir=ort_dir)
    else:
        micmac.write_image_mesures(imlist, gcps, out_dir, subscript, ort_dir=None, ori=ori)

    print('running mm3d GCPBascule to estimate terrain errors')
    gcps = micmac.bascule(gcps, out_dir, match_pattern, subscript, ori)
    gcps['res_dist'] = np.sqrt(gcps.xres ** 2 + gcps.yres ** 2)

    gcps = gcps.loc[np.abs(gcps.res_dist - gcps.res_dist.median()) < 2 * nmad(gcps.res_dist)]
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

        if Path(ort_dir).exists():
            micmac.write_image_mesures(imlist, cps, out_dir, subscript, ort_dir=ort_dir, outname='AutoCPMeasures')
        else:
            micmac.write_image_mesures(imlist, cps, out_dir, subscript, ort_dir=None, outname='AutoCPMeasures', ori=ori)

    # now, iterate campari to refine the orientation
    gcps = micmac.iterate_campari(gcps, out_dir, match_pattern, subscript, dx=ref_img.res[0],
                                  rel_ori=ori, allfree=allfree, max_iter=max_iter, homol=dir_homol)

    if use_cps:
        cp_resids = micmac.checkpoints(match_pattern, f"Ori-TerrainFinal{subscript}",
                                       fn_cp=Path(out_dir, f"AutoCPs{subscript}.xml"),
                                       fn_meas=Path(out_dir, f"AutoCPMeasures{subscript}-S2D.xml"),
                                       fn_resids=Path(out_dir, f"AutoCPs{subscript}", ret_df=True))

        # merge cps and cp_resids
        cps.drop(['xres', 'yres', 'zres', 'residual', 'res_dist'], axis=1, inplace=True)

        cps = cps.merge(cp_resids, left_on='id', right_on='id')
        cps['res_xy'] = np.sqrt(cps.xres ** 2 + cps.yres ** 2)
        cps['res_dist'] = np.sqrt(cps.xres ** 2 + cps.yres ** 2 + cps.zres ** 2)

        cps.to_file(Path(out_dir, f"AutoCPs{subscript}.shp"))

    # final write of gcps to disk.
    gcps.to_file(Path(out_dir, f"AutoGCPs{subscript}.shp"))

    fig1, ax1 = plt.subplots(1, 1, figsize=(7, 5))
    ax1.imshow(reg_img[::5, ::5], cmap='gray', extent=(0, reg_img.shape[1], reg_img.shape[0], 0))
    ax1.plot(gcps.orig_j, gcps.orig_i, 'r+')
    ax1.quiver(gcps.orig_j, gcps.orig_i, gcps.camp_xres, gcps.camp_yres, color='r')
    if use_cps:
        ax1.plot(cps.orig_j, cps.orig_i, 'bs')
        ax1.quiver(cps.orig_j, cps.orig_i, cps.xres, cps.yres, color='b')
    ax1.set_xticks([])
    ax1.set_yticks([])

    fig1.savefig(Path(out_dir, f"relative_gcps{subscript}.png"), bbox_inches='tight', dpi=200)
    plt.close(fig1)

    fig2, ax2 = plt.subplots(1, 1, figsize=(7, 5))
    xmin, ymin, xmax, ymax = ref_img.bounds
    ax2.imshow(ref_img.data[::10, ::10], cmap='gray', extent=(xmin, xmax, ymin, ymax))
    ax2.plot(gcps.geometry.x, gcps.geometry.y, 'r+')
    footprints.to_crs(ref_img.crs).boundary.plot(ax=ax2, color='b')
    if use_cps:
        ax2.plot(cps.geometry.x, cps.geometry.y, 'bs')

    ax2.set_xlim(xmin, xmax)
    ax2.set_ylim(ymin, ymax)

    ax2.set_xticks([])
    ax2.set_yticks([])

    fig2.savefig(Path(out_dir, f"world_gcps{subscript}.png"), bbox_inches='tight', dpi=200)
    plt.close(fig2)

    fig = _plot_residuals(gcps, res='camp')
    fig.savefig(Path(out_dir, f"gcp_residuals_plot{subscript}.png"), bbox_inches='tight', dpi=200)
    plt.close(fig)

    if use_cps:
        fig = _plot_residuals(cps, res='')
        fig.savefig(Path(out_dir, f"cp_residuals_plot{subscript}.png"), bbox_inches='tight', dpi=200)
        plt.close(fig)

    print('cleaning up.')

    for txtfile in glob('Auto*.tif.txt') + \
                   glob('NoDist*.tif.txt'):
        os.remove(txtfile)
    print('end.')
    # embed()
