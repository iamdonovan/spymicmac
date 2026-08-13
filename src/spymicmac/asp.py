"""
spymicmac.asp is a collection of tools for interfacing with Ames Stereo Pipeline
"""
import os
import io
from pathlib import Path
from glob import glob
import subprocess
import numpy as np
import pyproj
import geoutils as gu
import pandas as pd
import geopandas as gpd
from osgeo import gdal
from rasterio.crs import CRS
from shapely.ops import split, orient
from shapely.geometry import LineString, Point, Polygon
from sklearn.neighbors import BallTree
from . import data, declass, micmac, register
from typing import Union


gdal.UseExceptions()


def _isaft(fn_img: str) -> bool:
    return os.path.splitext(fn_img)[0][-4] == 'A'


def sorted_framelist(globstr: str, split_ext: bool = False, root_dir: Union[str, Path] = '.') -> list:
    """
    Get a list of images in the current directory, sorted by forward/aft camera and then frame type.

    :param globstr: the glob pattern get a list of filenames (e.g., "D3C*.tif" for KH-9 PC images)
    :param split_ext: remove extension from filenames
    :param root_dir: the root directory to search for image files
    """
    fn_imgs = [fn for fn in glob(globstr, root_dir=root_dir) if 'map' not in fn]

    if split_ext:
        fn_imgs = [os.path.splitext(fn)[0] for fn in fn_imgs]

    imgs = pd.DataFrame(data={'filename': fn_imgs})
    imgs['is_aft'] = [_isaft(fn) for fn in imgs['filename']]
    imgs['frame'] = [os.path.splitext(fn)[0][-3:] for fn in imgs['filename']]

    return list(imgs.sort_values(['is_aft', 'frame'])['filename'])


def _parse_cam(fn_cam: str) -> dict:
    with open(fn_cam, 'r') as f:
        cam_lines = [l.strip() for l in f.readlines()]

    cam = dict()
    cam['version'] = cam_lines[0]
    cam['type'] = cam_lines[1]
    cam['image_size'] = tuple([int(p) for p in cam_lines[2].split(' = ')[-1].split()])

    for ll in cam_lines[2:]:
        name, val = ll.split(' = ')
        if len(val.split()) < 2:
            try:
                cam[name] = float(val)
            except ValueError as e:
                cam[name] = val
        else:
            try:
                cam[name] = [float(p) for p in val.split()]
            except ValueError as e:
                cam[name] = val

    return cam


def _init_center(fprint: Polygon) -> tuple[float, float, float]:
    cx, cy = fprint.centroid.x, fprint.centroid.y
    alt = 180000  # very rough estimated altitude of 180 km

    geocent = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    geodet = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

    transformer = pyproj.Transformer.from_proj(geodet, geocent)
    x, y, z = transformer.transform(cx, cy, alt)

    return x, y, z


def add_motion_comp(cam: str, params: dict) -> dict:
    """
    Add a default motion compensation factor value to write to an ASP camera.

    :param cam: the panoramic camera flavor to use. Must be
    :param params: the dict describing the camera attributes
    :returns: the updated parameter dict
    """
    # values based on reported results from Ghuffar et al. 2022
    imc_params = {'KH4': 0.014, 'KH4A': 0.014, 'KH4B': 1e-4}

    assert cam in imc_params.keys(), f"{cam} not recognized as a valid camera [{imc_params.keys()}]"

    params['motion_comp'] = imc_params[cam]
    return params


def optical_bar_cam(fn_img: str, flavor: str, out_name: str,
                    fprint: Union[Polygon, None] = None,
                    scan_res: float = 7e-6,
                    mean_el: Union[float, int] = 1000,
                    use_3d_vel: bool = True) -> None:
    """
    Generate a sample ASP camera file for a KH-4 Optical Bar camera.

    :param fn_img: the filename of the image. Used to read the image size, and determine whether the image is from
        the aft or forward camera.
    :param flavor: what type of camera the image came from - currently either KH4 or KH9
    :param out_name: the filename to write the camera file to.
    :param fprint: an optional image, footprint used to estimate the initial camera position
    :param scan_res: the image scanning resolution, in m per pixel (e.g., 7 microns -> 7.0e-6)
    :param mean_el: the mean elevation covered by the image
    :param use_3d_vel: use a 3D velocity vector, rather than a 1D speed. Requires ASP 3.6.0 or greater.
    """
    assert flavor in declass.sample_params.keys(), f"flavor must be one of {list(declass.sample_params.keys())}"
    ds = gdal.Open(fn_img)
    width, height = ds.RasterXSize, ds.RasterYSize
    cx, cy = width / 2, height / 2

    ds = None  # close the image

    params = declass.sample_params[flavor]

    if flavor == 'KH4':
        cam = declass.get_declass_camera(fn_img)
        params = add_motion_comp(cam, params)
    else:
        params['motion_comp'] = 0.014

        frame_width_cm = scan_res * width * 100 # frame width in cm
        scan_time, scan_angle = declass.match_scan_angle(frame_width_cm)
        print(f"Estimate frame width is {frame_width_cm} cm, corresponding to a {scan_angle}° scan.")
        print(f"Using initial scan time guess of {scan_time:.4f} s.")

        params['scan_time'] = scan_time

    with open(out_name, 'w') as f:
        print('VERSION_4', file=f)
        print('OPTICAL_BAR', file=f)

        print(f'image_size = {width} {height}', file=f)
        print(f'image_center = {cx} {cy}', file=f)

        print(f'pitch = {scan_res}', file=f)
        print(f'f = {params["f"]}', file=f)

        if use_3d_vel:
            print(f'scan_time = 1', file=f)
            print(f'forward_tilt = 0', file=f)
        else:
            print(f'scan_time = {params["scan_time"]}', file=f)
            if _isaft(fn_img):
                print(f'forward_tilt = {-params["tilt"]}', file=f)
            else:
                print(f'forward_tilt = {params["tilt"]}', file=f)

        if fprint is not None and not use_3d_vel:
            icx, icy, icz = _init_center(fprint)
        else:
            icx, icy, icz = 0, 0, 0
        print(f'iC = {icx} {icy} {icz}', file=f)
        print('iR = 1 0 0 0 1 0 0 0 1', file=f)

        if use_3d_vel:
            print('speed = 0', file=f)
        else:
            print(f'speed = {params["speed"]}', file=f)

        print('mean_earth_radius = 6371000', file=f)
        # need a better value than this
        print(f"mean_surface_elevation = {mean_el}", file=f)

        if use_3d_vel:
            print(f"motion_compensation_factor = 0", file=f)
            print('scan_dir = right', file=f)
            print('velocity = 0 0 0', file=f)

        else:
            print(f"motion_compensation_factor = {params['motion_comp']}", file=f)
            if _isaft(fn_img):
                print('scan_dir = left', file=f)
            else:
                print('scan_dir = right', file=f)


def cam_from_footprint(fn_img: str, flavor: str, scan_res: float, fn_dem: Union[str, Path],
                       north_up: bool = True, footprints: gpd.GeoDataFrame = None, datum: Union[str, None] = None,
                       mean_el: Union[float, int, None] = 1000, use_3d_vel: bool = True):
    """
    Generate a camera (.tsai) file from an image footprint.

    :param fn_img: the filename of the image to generate a camera for.
    :param flavor: what type of camera the image came from - currently either KH4 or KH9
    :param scan_res: the scanning resolution of the image
    :param fn_dem: the filename of the reference DEM
    :param north_up: whether the top of the image corresponds to North or not
    :param footprints: a GeoDataFrame containing image footprints and an ID field with image names. If not
        provided, will attempt to download from USGS.
    :param datum: the geodetic datum to use. If not set, will try to guess from the DEM (or default to WGS84). See
        ASP docs for list of options.
    :param mean_el: the mean surface elevation covered by the image. If None, uses DEM and footprint to
        calculate the value.
    :param use_3d_vel: use a 3D velocity vector, rather than a 1D speed. Requires ASP 3.6.0 or greater.
    :return:
    """
    clean_name = fn_img.split('OIS-Reech_')[-1].split('.tif')[0]

    # now, get the image footprint, and use ul_corner to get the ul, ur, lr, ll coordinates
    if footprints is None:
        footprints = data.get_usgs_footprints([clean_name], dataset=declass.usgs_datasets[flavor])
        fprint = footprints.loc[0, 'geometry']
    else:
        fprint = footprints.loc[footprints['ID'] == clean_name, 'geometry'].values[0]

    if mean_el is None:
        dem = gu.Raster(fn_dem)
        mask = gu.Vector(footprints.loc[footprints['ID'] == clean_name]).create_mask(dem)
        mean_el = dem[mask].mean()

    if _isaft(fn_img):
        optical_bar_cam(fn_img, flavor, 'samp_aft.tsai', fprint,
                        scan_res=scan_res, mean_el=mean_el, use_3d_vel=use_3d_vel)
        fn_samp = 'samp_aft.tsai'
    else:
        optical_bar_cam(fn_img, flavor, 'samp_for.tsai', fprint,
                        scan_res=scan_res, mean_el=mean_el, use_3d_vel=use_3d_vel)
        fn_samp = 'samp_for.tsai'

    coords = _stanrogers(fprint, north_up)

    cl_args = ['cam_gen', '--sample-file', fn_samp, '--camera-type', 'opticalbar',
               '--lon-lat-values', '  '.join([f'{x} {y}' for x, y in coords]), fn_img,
               '--reference-dem', fn_dem, '--gcp-file', fn_img.replace('.tif', '-cam.gcp'),
               '--refine-camera', '-o', fn_img.replace('.tif', '.tsai')]

    if datum is not None:
        cl_args.append(['--datum', datum])

    print(cl_args)

    p = subprocess.Popen(cl_args)
    p.wait()

    os.remove(fn_samp)


# helper functions to help sort polygons from north to south and east to west
def _cenlat(poly: Polygon) -> float:
    return poly.centroid.y


def _cenlon(poly: Polygon) -> float:
    return poly.centroid.x


# return the upper left, upper right, lower right, lower left coordinates for an image
def _stanrogers(fprint: Polygon, north_up: bool) -> tuple[tuple[float, float], ...]:

    # oriented_envelope (mrr) goes lr, ur, ul, ll
    # use orient to ensure that it is properly oriented - for some reason this isn't always the case with mrr?
    outer = orient(fprint.buffer(0.05).minimum_rotated_rectangle)
    inner = orient(fprint.buffer(0.01).minimum_rotated_rectangle)
    x, y = outer.exterior.coords.xy
    coords = np.array(list(zip(x, y)))

    # get the right, top, left, bottomrig sides of the envelope
    right = LineString(coords[0:2])
    top = LineString(coords[1:3])
    left = LineString(coords[2:4])
    bot = LineString(coords[3:])

    horizontal = LineString([left.centroid, right.centroid])
    vertical = LineString([top.centroid, bot.centroid])

    # split the envelope into upper and lower halves
    # split the geometry, sort by centroid latitude, in descending order if the top of the image is north
    lower, upper = sorted(split(inner, horizontal).geoms, key=_cenlat, reverse=(not north_up))

    upper_left, upper_right = sorted(split(upper, vertical).geoms, key=_cenlon, reverse=(not north_up))
    lower_left, lower_right = sorted(split(lower, vertical).geoms, key=_cenlon, reverse=(not north_up))

    fx, fy = fprint.exterior.xy
    vertices = np.array([Point(x, y) for x, y in zip(fx[:-1], fy[:-1])])

    # upper left, upper right, lower right, lower left
    ul = vertices[[upper_left.contains(pt) for pt in vertices]][0]
    ur = vertices[[upper_right.contains(pt) for pt in vertices]][0]
    lr = vertices[[lower_right.contains(pt) for pt in vertices]][0]
    ll = vertices[[lower_left.contains(pt) for pt in vertices]][0]

    return (ul.x, ul.y), (ur.x, ur.y), (lr.x, lr.y), (ll.x, ll.y)


def bundle_adjust(fn_imgs: Union[list[Union[str, Path]], str],
                  out_prefix: str,
                  cam_prefix: Union[str, Path] = '',
                  map_suffix: Union[str, None] = None,
                  session_type: Union[str, None] = None,
                  gcp_patt: Union[str, None] = None,
                  num_iter: int = 20,
                  num_pass: int = 2,
                  ba_kwargs: dict = {},
                  ba_flags: list = []) -> None:
    """
    Run bundle_adjust on a given set of images. For more information about bundle_adjust, see the ASP documentation.

    :param fn_imgs: the filename of the images to run bundle_adjust on.
    :param out_prefix: the output prefix to use for the files produced by bundle_adjust.
    :param cam_prefix: the prefix/path to use for the camera (e.g., .tsai) files.
    :param map_suffix: the suffix used for the map-projected images, if map-projected images are being used. One of
        session_type or map_suffix must be specified.
    :param session_type: the stereo session type to use for processing. One of session_type or map_suffix must be
        specified.
    :param gcp_patt: the matching pattern for the GCP file(s) to use in the bundle adjustment.
    :param num_iter: the maximum number of iterations.
    :param num_pass: how many passes of bundle adjustment to do, with given number of iterations in each pass.
        For more than one pass, outliers will be removed between passes using --remove-outliers-params, and
        re-optimization will take place. Residual files and a copy of the match files with the outliers removed
        (*-clean.match) will be written to disk.
    :param ba_kwargs: additional kwargs to pass to bundle_adjust, for any arguments/flags that take a value.
        Keys should not include the '--' prefix - for example, use 'intrinsics-to-float' to define which
        intrinsics should be floated, rather than '--intrinsics-to-float'.
    :param ba_flags: additional flags to pass to bundle adjust, for any arguments/flags that do not take a value.
        Flags should not include the '--' prefix - for example, use 'fix-gcp-xyz' rather than '--fix-gcp-xyz'.
    """
    if map_suffix is None and session_type is None:
        raise KeyError('One of map_suffix or session_type must be specified.')

    if isinstance(fn_imgs, str):
        fn_imgs = sorted_framelist(fn_imgs)

    # get a list of image names without extensions
    clean_names = [os.path.splitext(fn)[0] for fn in fn_imgs]

    cl_args = ['bundle_adjust']

    # add the images
    cl_args.extend(fn_imgs)

    # add the cam files
    if cam_prefix == '':
        fn_cams = [fn + '.tsai' for fn in clean_names]
    else:
        fn_cams = [f"{cam_prefix}-{fn}.tsai" for fn in clean_names]
    cl_args.extend(fn_cams)

    if gcp_patt is not None:
        fn_gcp = glob(gcp_patt)
        cl_args.extend(fn_gcp)

    if map_suffix is not None:
        map_args = ['--mapprojected-data']
        fn_map = ['.'.join([fn, map_suffix]) for fn in clean_names]
        map_args.append(' '.join(fn_map))

        cl_args.extend(map_args)

    if session_type is not None:
        cl_args.extend(['-t', session_type])

    cl_args.extend(['-o', out_prefix])
    cl_args.extend(['--num-iterations', str(num_iter)])
    cl_args.extend(['--num-passes', str(num_pass)])

    for arg in ba_flags:
        cl_args.append('--' + arg)

    for kwarg in ba_kwargs:
        cl_args.extend(['--' + kwarg, str(ba_kwargs[kwarg])])

    print(cl_args)

    p = subprocess.Popen(cl_args)
    p.wait()


def meas_to_asp_gcp(fn_gcp: Union[str, Path], fn_meas: Union[str, Path], imlist: list,
                    outname: Union[str, None] = None, scale: int = 1, singles: bool = False) -> None:
    """
    Convert image measures stored in a micmac xml file to an ASP .gcp file format.

    :param str fn_gcp: the filename of the shapefile with the GCP coordinates
    :param str fn_meas: the filename of the xml file with the image measures
    :param list imlist: the image(s) to write point locations for
    :param str outname: the name of the output filename to create (default: {fn_meas}.gcp)
    :param int scale: the factor by which to scale the image point locations (default: 1)
    :param bool singles: write gcps present in only a single image (default: False)
    """
    if outname is None:
        outname = fn_meas.replace('.xml', '.gcp')

    gcps = gpd.read_file(fn_gcp).to_crs(crs='epsg:4326').set_index('id')
    meas = micmac.parse_im_meas(fn_meas)

    meas = meas.loc[meas['image'].isin(imlist)]

    gcp_list = sorted(meas.name.unique())

    write_asp_gcp(outname, gcps, gcp_list, imlist=imlist, scale=scale, singles=singles, meas=meas)

def mapproject(fn_dem: Union[str, Path], fn_img: Union[str, Path], fn_cam: Union[str, Path],
               res: Union[None, float, int]=None, fn_out: Union[str, Path]=None):
    """
    Run mapproject to project an image using a DEM.

    :param fn_dem: the filename of the DEM to use
    :param fn_img: the filename of the image
    :param fn_cam: the filename of the camera corresponding to the image
    :param res: the resolution (ground sample distance) of the output file
    :param fn_out: the filename of the output image
    """

    if fn_out is None:
        fn_out = '.'.join([os.path.splitext(fn_img)[0], 'map', os.path.splitext(fn_img)[-1]])

    cl_args = ['mapproject', fn_dem, fn_img, fn_cam, fn_out]

    if res is not None:
        cl_args.extend(['--tr', str(res)])

    p = subprocess.Popen(cl_args)
    p.wait()


def camera_footprint(fn_img: Union[str, Path],
                     fn_cam: Union[str, Path],
                     fn_dem: Union[str, Path],
                     fn_out: Union[str, Path],
                     quick: bool = True,
                     crs = None) -> None:
    """
    Project camera footprint to a shapefile using a camera file and a DEM.

    :param fn_img: the image filename
    :param fn_cam: the camera filename
    :param fn_dem: the DEM filename
    :param fn_out: the output filename
    :param quick: use a faster, less accurate projection.
    :param crs: the CRS to project the footprint to.
    """
    cl_args = ['camera_footprint', fn_img, fn_cam]

    cl_args.extend(['--output-kml', 'tmp.kml'])
    cl_args.extend(['--dem-file', fn_dem])

    if quick:
        cl_args.append('--quick')

    p = subprocess.Popen(cl_args)
    p.wait()

    tmp = gpd.read_file('tmp.kml')
    fp = Polygon(tmp.loc[0, 'geometry'])

    if crs is not None:
        out = gpd.GeoDataFrame(data={'filename': [os.path.splitext(fn_img)[0]], 'geometry': [fp]}, crs=tmp.crs).to_crs(crs)
    else:
        out = gpd.GeoDataFrame(data={'filename': [os.path.splitext(fn_img)[0]], 'geometry': [fp]}, crs=tmp.crs)

    out.to_file(fn_out)

    os.remove('tmp.kml')


def mapprojected_footprint(fn_img, out_crs=None, nodata=None) -> gpd.GeoDataFrame:
    """
    Get the footprint of the valid (not nodata) areas of a raster.

    :param fn_img: the image filename
    :param out_crs: the output CRS of the footprint
    :param nodata: the nodata value to use. If not provided, uses the image nodata value if it is set; otherwise, 0.
    """
    img = gu.Raster(fn_img)

    if nodata is None:
        if img.nodata is None:
            nodata = 0
        else:
            nodata = img.nodata

    tmp = (img != nodata).polygonize()

    fp = tmp.union_all().ds.loc[0, 'geometry']

    if out_crs is None:
        return gpd.GeoDataFrame(data={'filename': fn_img, 'geometry': fp}, index=[0], crs=img.crs)
    else:
        return gpd.GeoDataFrame(data={'filename': fn_img, 'geometry': fp}, index=[0], crs=img.crs).to_crs(out_crs)

def _gcp_row(_gcp, name, scale, meas):
    first = _gcp.iloc[0]

    if {'lon', 'lat'} <= set(_gcp.columns):
        lon, lat = first.lon, first.lat
    else:
        lon, lat = first.geometry.x, first.geometry.y

    out_gcp = ','.join([str(name).strip('GCP'), str(lat), str(lon), str(first.elevation),
                        str(first['sigma_x']), str(first['sigma_y']), str(first['sigma_z'])])

    for num in range(len(_gcp)):
        _meas = _gcp.iloc[num]

        if meas is not None:
            img = _meas['image_name']
            row, col = meas.loc[(meas.image == img) & (meas.name == name), ['i', 'j']].values[0]

            out_gcp += ',' + ','.join([img,
                                       f"{col / scale}",
                                       f"{row / scale}",
                                       f"{_meas['sx_px']}",
                                       f"{_meas['sy_px']}"])
        else:
            out_gcp += ',' + ','.join([_meas['image_name'],
                                       f"{_meas['pixel_x'] / scale}",
                                       f"{_meas['pixel_y'] / scale}",
                                       f"{_meas['sx_px']}",
                                       f"{_meas['sy_px']}"])

    return out_gcp


def write_asp_gcp(fn_gcp: Union[str, Path], gcp_df: gpd.GeoDataFrame,
                  gcp_list: Union[None, list] = None, imlist: Union[None, list] = None,
                  scale: int = 1, singles: bool = True, meas: Union[None, pd.DataFrame] = None,
                  headers: Union[None, list[str]] = None,
                  gcp_sig: Union[int, float, list[int], list[float], None] = None,
                  px_sig: Union[int, float, list[int], list[float], None] = None) -> None:
    """
    Write GCPs in ASP format.

    :param fn_gcp: the filename to write the GCPs to
    :param gcp_df: a GeoDataFrame of GCP locations
    :param gcp_list: a list of which GCPs to write
    :param imlist: a list of what images to write GCPs for
    :param scale: the scale to use for scaling image coordinates
    :param singles: whether to write GCPs that are only found in one image
    :param meas: a DataFrame of image measurements, as created by mm3d SaisieAppuis and read by micmac.parse_im_meas()
    :param headers: the header rows from the original .gcp file
    :param gcp_sig: the uncertainty of the GCP position, either as a list of values for (x, y, z), or a single value.
        If not set, looks for (sigma_x, sigma_y, sigma_z) in the GCP dataframe. If those columns are not present,
        defaults to 1.0
    :param px_sig: the uncertainty of the GCP location in the image, either as a list of values for (j, i), or a
        single value. If not set, looks for columns like (sx_px, sy_px) in the GCP dataframe. If those columns are
        not present, defaults to 1.0.
    """

    if gcp_sig is None:
        for dd in 'xyz':
            if not f"sigma_{dd}" in gcp_df.columns:
                gcp_df[f"sigma_{dd}"] = 1

    elif isinstance(gcp_sig, (int, float)):
        for dd in 'xyz':
            gcp_df[f"sigma_{dd}"] = gcp_sig

    elif isinstance(gcp_sig, (list, tuple)):
        assert len(gcp_sig) == 3, "must provide 3 values for GCP uncertainty (sx, sy , sz)"
        for ind, dd in enumerate('xyz'):
            gcp_df[f"sigma_{dd}"] = gcp_sig[ind]

    if 'id' not in gcp_df.columns:
        gcp_df['id'] = gcp_df.index

    if 'height_above_datum' in gcp_df.columns:
        gcp_df['elevation'] = gcp_df['height_above_datum']

    pivoted = pd.wide_to_long(gcp_df,
                              stubnames=['image_name', 'pixel_x', 'pixel_y', 'sx_px', 'sy_px'],
                              i='id', j='num',
                              sep='.', suffix=r'\d+')
    pivoted.dropna(inplace=True)

    if px_sig is None:
        for dd in 'xy':
            if not f"s{dd}_px" in pivoted.columns:
                gcp_df[f"s{dd}_px"] = 1

    elif isinstance(px_sig, (int, float)):
        for dd in 'xy':
            gcp_df[f"s{dd}_px"] = gcp_sig

    elif isinstance(px_sig, (list, tuple)):
        assert len(px_sig) == 2, "must provide 3 values for GCP uncertainty (sx, sy , sz)"
        for ind, dd in enumerate('xy'):
            gcp_df[f"s{dd}_px"] = px_sig[ind]

    with open(fn_gcp, 'w') as f:
        if headers is not None:
            for header in headers:
                print(header, file=f)

        if gcp_list is None:
            gcp_list = sorted(set(pivoted.index.get_level_values(0).tolist()))

        for gcp in gcp_list:
            _gcp = pivoted.loc[gcp]
            if not singles and len(_gcp) == 1:
                continue

            print(_gcp_row(_gcp, gcp, scale, meas), file=f)


def _parse_gcp(fn_gcp, delim=None):
    cols = ['id', 'lat', 'lon', 'height_above_datum', 'sigma_x', 'sigma_y', 'sigma_z']
    img_headers = ['image_name', 'pixel_x', 'pixel_y', 'sx_px', 'sy_px']


    with open(fn_gcp, 'r') as f:
        all_lines = [l.strip() for l in f.readlines()]

        crs_wkt = all_lines[0].strip().replace('# ', '')

        if delim is None:
            delim = ' '
            ncols = max([len(l.split(delim)) for l in all_lines[2:]])

            if ncols < 2:
                delim = ','
                ncols = max([len(l.split(delim)) for l in all_lines[2:]])
        else:
            ncols = max([len(l.split(delim)) for l in all_lines[2:]])

    nimg = 0
    while len(cols) < ncols:
        nimg += 1
        cols += ['.'.join([c, str(nimg)]) for c in img_headers]

    padded = [l + delim * (ncols - len(l.split(delim))) + '\n' for l in all_lines[2:]]

    cols += ['blank']
    df = pd.read_csv(io.StringIO(''.join(padded)), delimiter=delim, names=cols).set_index('id')
    del df['blank']

    return gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat, crs='epsg:4326'))


def mask_asp_dem_gcps(fn_gcp: Union[str, Path],
                      fn_inc: Union[None, str, Path, list[str], list[Path]] = None,
                      fn_exc: Union[None, str, Path, list[str], list[Path]] = None,
                      overwrite: bool = True) -> None:
    """
    Filter/mask GCPs from a .gcp file, using inclusion/exclusion masks. Inclusion masks indicate areas where points
    inside should be included as stable terrain; exclusion masks indicate areas where points inside should be excluded.
    At least one of fn_inc, fn_exc must be set in order to filter.

    TODO: implement multiple masks (WIP)
    :param fn_gcp: the .gcp file to filter
    :param fn_inc: the inclusion mask(s) to use, pointing to a vector file format
    :param fn_exc: the exclusion mask(s) to use, pointing to a vector file format
    :param overwrite: whether to overwrite existing file. If False, appends '_filt' to existing filename.
    """
    assert not all(el is None for el in [fn_inc, fn_exc]), "at least one mask has to be set!"

    gcps = _parse_gcp(fn_gcp)
    with open(fn_gcp, 'r') as f:
        headers = [l.strip() for l in f.readlines() if l.strip()[0] == '#']

    if fn_inc is not None:
        inc_mask = gpd.read_file(fn_inc).to_crs(gcps.crs)
        within_inc = inc_mask.sindex.query(gcps.geometry, predicate='intersects')[0]
    else:
        within_inc = gcps.index

    if fn_exc is not None:
        exc_mask = gpd.read_file(fn_exc).to_crs(gcps.crs)
        within_exc = exc_mask.sindex.query(gcps.geometry, predicate='intersects')[0]
    else:
        within_exc = gcps.index[~gcps.index.isin(gcps.index)]

    valid = (gcps.index.isin(within_inc)) & (~gcps.index.isin(within_exc))

    if overwrite:
        fn_out = fn_gcp
    else:
        fn_out = fn_gcp.replace('.gcp', '_filt.gcp')

    write_asp_gcp(fn_out, gcps.loc[valid], headers=headers)


def gcps_from_dem(img_pair: tuple[str, str],
                  fn_dem: Union[str, Path],
                  fn_ref: Union[str, Path],
                  camera_prefix: str,
                  fn_gcp: str,
                  fn_landmask: Union[str, Path, None] = None,
                  fn_glacmask: Union[str, Path, None] = None,
                  warp_prefix: str = 'warp/run',
                  match_prefix: Union[str, None] = None,
                  use_clean: bool = True,
                  ps_kwargs: dict = {},
                  ps_args: list = [],
                  gcp_kwargs: dict = {},
                  gcp_args: list = []) -> None:
    """
    Use ASP's parallel_stereo and dem2gcp to generate a GCP file for a pair of images.

    Runs the following steps, in order:
        - uses dem_mosaic to blur the input (warped) dem
        - uses gdaldem hillshade to generate a hillshade of the blurred warped dem
        - masks the reference DEM to the valid areas of the blurred warped DEM hillshade
        - uses gdaldem hillshade to generate a hillshade of the cropped, blurred reference dem
        - runs parallel_stereo in correlator-mode to find the disparity between the two DEM hillshade (if do_corr=True)
        - runs dem2gcp to generate a .gcp file with GCPs for the two input images.

    :param img_pair: the names of the two images to use. Should be ordered as (left_img, right_img), the same way that
        parallel_stereo was used to generate the warped DEM.
    :param fn_dem: the filename of the input (warped) DEM.
    :param fn_ref: the filename of the reference DEM.
    :param camera_prefix: the prefix of the camera files for the two input images (e.g., ba/run)
    :param fn_gcp: the filename of the output GCP file. Should use the extension .gcp
    :param fn_landmask: the (optional) filename for an inclusion (i.e., land/stable terrain) mask for matching.
    :param fn_glacmask: the (optional) filename for an exclusion (i.e., unstable terrain) mask for matching.
    :param warp_prefix: the prefix to save the output of parallel_stereo to. Defaults to warp/run
    :param match_prefix: the prefix to use for the match files. Defaults to the same as camera_prefix.
    :param use_clean: use the 'clean' match file, rather than the full match file
    :param ps_kwargs: optional keyword parameters for parallel_stereo, for any arguments/flags that take a value.
        Keys should not include the '--' prefix - for example, use 'stereo-algorithm' to define which stereo algorithm
        to use, rather than '--stereo-algorithm'.
    :param ps_args: optional flags to pass to parallel_stereo, for any arguments/flags that do not take a value.
        Flags should not include the '--' prefix.
    :param gcp_kwargs: optional keyword parameters for dem2gcp, for any arguments/flags that take a value.
        Keys should not include the '--' prefix - for example, use 'gcp-sigma' to set the GCP uncertainty
        to use, rather than '--gcp-sigma'.
    :param gcp_args: optional flags to pass to dem2gcp, for any arguments/flags that do not take a value.
        Flags should not include the '--' prefix.
    """

    if match_prefix is None:
        match_prefix = camera_prefix

    left, right = img_pair

    print(f"Running dem_mosaic on {fn_dem}.")
    p = subprocess.Popen(['dem_mosaic', '--dem-blur-sigma', '5', fn_dem,
                          '-o', 'tmp_warp_dem.tif'])
    p.wait()

    print('Creating hillshade of warped DEM')
    p = subprocess.Popen(['gdaldem', 'hillshade', '-multidirectional', '-compute_edges',
                          'tmp_warp_dem.tif', 'tmp_warp_dem_hs.tif'])
    p.wait()

    print('Cropping/masking reference DEM to valid area of warped DEM.')
    mapprojected_footprint('tmp_warp_dem_hs.tif').to_file('tmp_mask.gpkg')
    tmp_masked = data.crop_mask_dem(fn_ref,
                                    'tmp_mask.gpkg',
                                    buff=1000,
                                    use_rect=False)
    tmp_masked.to_file('tmp_dem.tif')

    print("Running dem_mosaic on cropped/masked reference DEM.")
    p = subprocess.Popen(['dem_mosaic', '--dem-blur-sigma', '5', 'tmp_dem.tif',
                          '-o', 'tmp_blur.tif'])
    p.wait()

    print('Creating hillshade of cropped/masked reference DEM')
    p = subprocess.Popen(['gdaldem', 'hillshade', '-multidirectional', '-compute_edges',
                          'tmp_blur.tif', 'tmp_blur_hs.tif'])
    p.wait()

    if not ps_kwargs:
        ps_kwargs = {'stereo-algorithm': 'asp_mgm', 'subpixel-mode': 9, 'ip-per-tile': 1000}
    elif 'stereo-algorithm' not in ps_kwargs:
        ps_kwargs['stereo-algorithm'] = 'asp_mgm'

    ps_cl_args = ['parallel_stereo', '--correlator-mode']
    for kwarg in ps_kwargs:
        ps_cl_args.extend(['--' + kwarg, str(ps_kwargs[kwarg])])

    for arg in ps_args:
        ps_cl_args.append('--' + arg)

    ps_cl_args.extend(['tmp_warp_dem_hs.tif', 'tmp_blur_hs.tif', warp_prefix])

    p = subprocess.Popen(ps_cl_args)
    p.wait()

    # now, call dem2gcp on the output
    gcp_cl_args = ['dem2gcp', '--warped-dem', 'tmp_warp_dem_hs.tif', '--ref-dem', 'tmp_dem.tif',
                   '--warped-to-ref-disparity', f"{warp_prefix}-F.tif",
                   '--left-image', left,
                   '--right-image', right,
                   '--left-camera', f"{camera_prefix}-{left.replace('.tif', '')}.tsai",
                   '--right-camera', f"{camera_prefix}-{right.replace('.tif', '')}.tsai"]

    if use_clean:
        match_ext = '-clean.match'
    else:
        match_ext = '.match'

    matchfile = f"{match_prefix}-{left.replace('.tif', '')}__{right.replace('.tif', '')}{match_ext}"
    if Path(matchfile).exists():
        gcp_cl_args.extend(['--match-file', matchfile,])
    else:
        matchfile = f"{match_prefix}-{right.replace('.tif', '')}__{left.replace('.tif', '')}{match_ext}"
        if Path(matchfile).exists():
            gcp_cl_args.extend(['--match-file', matchfile,])
        else:
            raise FileNotFoundError(f"Unable to find match file in {match_prefix} for {left}, {right}")

    gcp_cl_args.extend(['--output-gcp', fn_gcp])

    if not gcp_kwargs:
        gcp_kwargs = {'max-num-gcp': 20000}

    if 'gcp-sigma' not in gcp_kwargs:
        dem = gu.Raster('tmp_dem.tif')
        gcp_kwargs['gcp-sigma'] = dem.res[0] / 4

    for kwarg in gcp_kwargs:
        gcp_cl_args.extend(['--' + kwarg, str(gcp_kwargs[kwarg])])

    p = subprocess.Popen(gcp_cl_args)
    p.wait()

    tmp_files = ['tmp_dem.tif', 'tmp_blur.tif', 'tmp_blur_hs.tif',
                 'tmp_warp_dem.tif', 'tmp_warp_dem_hs.tif', 'tmp_mask.gpkg']
    log_files = glob('tmp_warp_dem.tif*.txt')
    log_files += glob('tmp_blur.tif*.txt')

    for fn_tmp in tmp_files + log_files:
        os.remove(fn_tmp)

    if fn_glacmask is not None or fn_landmask is not None:
        mask_asp_dem_gcps(fn_gcp, fn_inc=fn_landmask, fn_exc=fn_glacmask)


def gcps_from_ortho(fn_img: Union[str, Path],
                    fn_ref: Union[str, Path],
                    fn_dem: Union[str, Path],
                    fn_out: Union[str, Path],
                    out_pre: str = 'gcp_ortho/run',
                    fn_mapproj: Union[str, Path, None] = None,
                    fn_footprints: Union[str, Path, None] = None,
                    fn_landmask: Union[str, Path, None] = None,
                    fn_glacmask: Union[str, Path, None] = None,
                    use_asp: bool = False,
                    kwargs: dict = {},
                    args: list = [],
                    spacing: int = 200,
                    thresh: int = 100,
                    matching_kwargs: dict = {}) -> None:
    """
    Find GCPs for an image by matching between the image (or the mapprojected version of the image) and a
    reference ortho/satellite image, using ASP's gcp_gen.

    If use_asp = False (the default), uses template matching between the (normally mapprojected image) and the
    reference image. Otherwise, uses ASP interestpoint matching.

    :param fn_img: the image to find GCPs for
    :param fn_ref: the reference ortho/satellite image
    :param fn_dem: the reference DEM to use for the GCP elevation
    :param fn_out: the output filename for the .gcp file.
    :param out_pre: the output prefix to use for the outputs from gcp_gen
    :param fn_mapproj: the (optional) filename for the mapprojected image to use.
    :param fn_footprints: the (optional) footprint file to use to crop/mask the reference image before matching.
    :param fn_landmask: the (optional) filename for an inclusion (i.e., land/stable terrain) mask for matching.
    :param fn_glacmask: the (optional) filename for an exclusion (i.e., unstable terrain) mask for matching.
    :param use_asp: whether to use ASP's matching, or use template matching from spymicmac.matching.
    :param kwargs: optional keyword parameters for gcp_gen, for any arguments/flags that take a value.
        Keys should not include the '--' prefix - for example, use 'gcp-sigma' to set the GCP uncertainty
        to use, rather than '--gcp-sigma'.
    :param args: optional flags to pass to gcp_gen, for any arguments/flags that do not take a value.
        Flags should not include the '--' prefix.
    :param spacing: the spacing (in pixels) to use for finding GCPs, if using spymicmac.matching
    :param thresh: the residual threshold (in pixels) to use for RANSAC, if using spymicmac.matching
    :param matching_kwargs: additional kwargs to pass to spymicmac.matching.find_matches.
    """
    clean_name = fn_img.split('OIS-Reech_')[-1].split('.tif')[0]

    if use_asp:
        do_mask = fn_mapproj is not None or fn_footprints is not None
        if do_mask:
            print('Cropping/masking reference orthoimage to image footprint.')
            if fn_mapproj is not None:
                mapprojected_footprint(fn_mapproj).to_file('tmp_mask.gpkg')
            else:
                footprints = gpd.read_file(fn_footprints)
                footprints[footprints['ID'] == clean_name].to_file('tmp_mask.gpkg')

            tmp_masked = data.crop_mask_dem(fn_ort,
                                            'tmp_mask.gpkg',
                                            buff=1000,
                                            use_rect=False)
            if tmp_masked.nodata is None:
                tmp_masked.set_nodata(0)

            tmp_masked.to_file('tmp_ortho.tif')
            fn_ref = 'tmp_ortho.tif'
        else:
            print('Using uncropped, unmasked orthoimage.')

    elif 'match-file' in kwargs:
        print(f"Using match file: {kwargs['match-file']}")

    else:
        assert fn_mapproj is not None, "map-projected image must be provided if not using ASP matching."
        gcps = register.register_ortho(fn_mapproj,
                                       fn_ref,
                                       spacing = spacing,
                                       fn_landmask = fn_landmask,
                                       fn_glacmask = fn_glacmask,
                                       thresh = thresh,
                                       matching_kwargs = matching_kwargs)

        ortho_match = gcps[['match_j', 'match_i']].copy().rename(columns={'match_j': 'x', 'match_i': 'y'})
        ortho_match[['ix', 'iy']] = np.ceil(ortho_match[['x', 'y']]).astype(int)

        ref_match = gcps[['full_j', 'full_i']].copy().rename(columns={'full_j': 'x', 'full_i': 'y'})
        ref_match[['ix', 'iy']] = np.ceil(ref_match[['x', 'y']]).astype(int)

        for cn in ['orientation', 'scale', 'interest', 'polarity', 'octave', 'scale_lv', 'num_descr']:
            ortho_match[cn] = 0
            ref_match[cn] = 0

        ortho_match['scale'] = 1
        ref_match['scale'] = 1

        match_pts = pd.concat([ortho_match, ref_match], ignore_index=False).reset_index(drop=True)

        with open("tmp_match.txt", 'w') as f:
            print(f"{len(gcps)} {len(gcps)}", file=f)

            for ind in range(len(match_pts)):
                print(match_pts.loc[[ind]].to_string(header=False, index=False), file=f)

        left_fn = os.path.splitext(os.path.basename(fn_mapproj))[0]
        right_fn = os.path.splitext(os.path.basename(fn_ref))[0]

        parse_args = ['parse_match_file.py', 'tmp_match.txt', f"{out_pre}-{left_fn}__{right_fn}.match", '-rev']
        p = subprocess.Popen(parse_args)
        p.wait()

        #kwargs['match-file'] = f"{fn_mapproj.replace('.tif', '')}__ref_ortho_match.txt"

    cl_args = ['gcp_gen',
               '--camera-image', fn_img,
               '--ortho-image', fn_ref,
               '--dem', fn_dem,
               '--output-gcp', fn_out,
               '--output-prefix', out_pre]

    if fn_mapproj is not None:
        kwargs['mapproj-image'] = fn_mapproj

    if 'gcp-sigma' not in kwargs:
        ortho = gu.Raster(fn_ref)
        kwargs['gcp-sigma'] = ortho.res[0]

    for arg in args:
        cl_args.append('--' + arg)

    for kwarg in kwargs:
        cl_args.extend(['--' + kwarg, str(kwargs[kwarg])])

    print(cl_args)

    p = subprocess.Popen(cl_args)
    p.wait()

    tmp_files = ['tmp_ortho.tif', 'tmp_mask.gpkg']
    for fn_tmp in tmp_files:
        if os.path.exists(fn_tmp):
            os.remove(fn_tmp)


def merge_gcps(left: Union[str, Path, gpd.GeoDataFrame],
               right: Union[str, Path, gpd.GeoDataFrame],
               crs: Union[CRS, str, int, None],
               tol: float = 1e-3,
               nimg: bool = 2,
               no_singles: bool = False) -> gpd.GeoDataFrame:
    """
    Merge two sets of GCPs together by using sjoin_nearest to determine duplicated GCPs.

    :param left: the "left" set of GCPs to use in the join; can be a parsed GCP file or a filename.
    :param right: the "right" set of GCPs to use in the join; can be a parsed GCP file or a filename.
    :param crs: the (projected) CRS to use for joining the GCPs to the residual pointmap, using sjoin_nearest
    :param tol: the tolerance used to determine duplicated GCPs
    :param nimg: the number to use for the suffix of the repeated columns (e.g., image_name, pixel_x, pixel_y, etc.).
        If only 2 sets of GCPs are joined, use 2. If a third is being joined to a previously merged pair, use 3, and
        so on.
    :param no_singles: return only GCPs seen by more than one image. Defaults to using all GCPs.
    :return: the merged set of GCPs.
    """
    right_cols = ['image_name.1', 'pixel_x.1', 'pixel_y.1', 'sx_px.1', 'sy_px.1', 'geometry']
    new_cols = [c.replace('.1', f".{nimg}") for c in right_cols]

    if isinstance(left, (str, Path)):
        left = _parse_gcp(left).to_crs(crs)
    elif isinstance(left, gpd.GeoDataFrame):
        left = left.to_crs(crs)
    else:
        raise TypeError(f"must be a filename/path or a GeoDataFrame: {left}")

    if isinstance(right, (str, Path)):
        right = _parse_gcp(right).to_crs(crs)
    elif isinstance(right, gpd.GeoDataFrame):
        right = right.to_crs(crs)
    else:
        raise TypeError(f"must be a filename/path or a GeoDataFrame: {right}")

    if no_singles:
        how = 'inner'
    else:
        how = 'left'

    merged = left.sjoin_nearest(right[right_cols].rename(columns=dict(zip(right_cols, new_cols))),
                                max_distance=tol, how=how)

    if not no_singles:
        if 'id_right' in merged.columns:
            id_col = 'id_right'
        else:
            id_col = 'id'

        missing = ~right.index.isin(merged[id_col])
        merged = pd.concat([merged, right.loc[missing]], ignore_index=True)

    return merged.drop(columns=[id_col]).reset_index(drop=True)


def parse_pointmap(fn_csv: Union[str, Path], gcps_only: bool = True) -> gpd.GeoDataFrame:
    """
    Parse a residual pointmap CSV, output by bundle_adjust.

    :param fn_csv: the filename of the CSV to parse (e.g., ba/all-final_residuals_pointmap.csv')
    :param gcps_only: return only the residuals from GCPs
    :return: the parsed pointmap file as a GeoDataFrame with CRS EPSG:4326
    """
    pointmap = pd.read_csv(fn_csv, header=2, names=['lon', 'lat', 'elev', 'res', 'source'])
    pointmap['source'] = pointmap['source'].str.split('#', expand=True)[1].str.strip()

    if gcps_only:
        pointmap = pointmap.loc[pointmap['source'] == 'GCP'].copy()

    return gpd.GeoDataFrame(pointmap, geometry=gpd.points_from_xy(x=pointmap['lon'], y=pointmap['lat'], crs=4326))


def filter_gcps_pointmap(fn_gcp: Union[str, Path],
                         fn_pointmap: Union[str, Path],
                         crs: Union[CRS, str, int, None],
                         nfact: int = 2,
                         thresh: Union[float, None] = None,
                         use_neighbors: bool = False,
                         num_neighbors: int = 10,
                         max_dist: float = 1e4) -> gpd.GeoDataFrame:
    """
    Filter GCPs based on their residual, as read from the residuals pointmap CSV output by bundle_adjust.

    Note that this only works if bundle_adjust was called with --fix-gcp-xyz!

    :param fn_gcp: the filename of the GCPs
    :param fn_pointmap: the filename of the pointmap CSV to parse (e.g., ba/all-final_residuals_pointmap.csv')
    :param crs: the (projected) CRS to use for joining the GCPs to the residual pointmap, using sjoin_nearest
    :param nfact: the multiple of the NMAD to use as a filter. Values more than this times the NMAD away from the
        median residual will be discarded.
    :param thresh: the threshold value to use as a filter. If specified, uses this value directly; otherwise,
        threshold is calculated as a multiple of the NMAD value.
    :param use_neighbors: compare the residual of each point to its nearest neighbors, rather than the statistics of
        the entire datset
    :param num_neighbors: the number of nearest neighbors to use
    :param max_dist: the maximum distance (in the units of the projected crs) to use for neighbor-based filtering
    :return: the filtered GCP geodataframe
    """
    pointmap = parse_pointmap(fn_pointmap)
    gcps = _parse_gcp(fn_gcp)

    joined = gcps[['geometry']].to_crs(crs).sjoin_nearest(pointmap.to_crs(crs), distance_col='dist')
    joined = joined[joined['dist'] < 1e-3].copy()

    if use_neighbors:
        xy = np.hstack([joined.geometry.x.values.reshape(-1, 1), joined.geometry.y.values.reshape(-1, 1)])
        tree = BallTree(xy, leaf_size=15)

        for ind in joined.index:
            this_xy = np.array([joined.loc[ind].geometry.x, joined.loc[ind].geometry.y]).reshape(1, -1)
            distances, indices = tree.query(this_xy, k=num_neighbors + 1) # get the num_neighbors closest + self

            indices = indices[np.logical_and(distances > 0, distances < max_dist)]
            joined.loc[ind, 'res_diff'] = (joined.loc[ind, 'res'] - joined.loc[joined.index[indices], 'res']).median()

    if thresh is not None:
        if not use_neighbors:
            valid = joined['res'] < thresh
        else:
            valid = joined['res_diff'].abs() < thresh
    else:
        if not use_neighbors:
            valid = np.abs(joined['res'] - joined['res'].median()) < nfact * register.nmad(joined['res'])
        else:
            valid = joined['res_diff'].abs() < nfact * register.nmad(joined['res_diff'])

    return gcps.loc[gcps.index.isin(joined.loc[valid].index)]
