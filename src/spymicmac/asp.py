"""
spymicmac.asp is a collection of tools for interfacing with Ames Stereo Pipeline
"""
import os
from pathlib import Path
import subprocess
import numpy as np
import pyproj
import geoutils as gu
import pandas as pd
import geopandas as gpd
from osgeo import gdal
from shapely.ops import split, orient
from shapely.geometry import LineString, Point, Polygon
from . import data, declass, micmac
from typing import Union


gdal.UseExceptions()


def _isaft(fn_img: str) -> bool:
    return os.path.splitext(fn_img)[0][-4] == 'A'


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

    with open(out_name, 'w') as f:
        print('VERSION_4', file=f)
        print('OPTICAL_BAR', file=f)

        print(f'image_size = {width} {height}', file=f)
        print(f'image_center = {cx} {cy}', file=f)

        print(f'pitch = {scan_res}', file=f)
        print(f'f = {params["f"]}', file=f)
        print(f'scan_time = {params["scan_time"]}', file=f)

        if _isaft(fn_img):
            print(f'forward_tilt = {-params["tilt"]}', file=f)
        else:
            print(f'forward_tilt = {params["tilt"]}', file=f)

        if fprint is not None:
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
        print(f"motion_compensation_factor = {params['motion_comp']}", file=f)

        if use_3d_vel:
            print('scan_dir = right', file=f)
        else:
            if _isaft(fn_img):
                print('scan_dir = left', file=f)
            else:
                print('scan_dir = right', file=f)

        # new format
        if use_3d_vel:
            print('velocity = 0 0 0', file=f)

def cam_from_footprint(fn_img: str, flavor: str, scan_res: float, fn_dem: Union[str, Path],
                       north_up: bool=True, footprints: gpd.GeoDataFrame=None, mean_el: Union[float, int]=1000):
    """
    Generate a camera (.tsai) file from an image footprint.

    :param fn_img: the filename of the image to generate a camera for.
    :param flavor: what type of camera the image came from - currently either KH4 or KH9
    :param scan_res: the scanning resolution of the image
    :param fn_dem: the filename of the reference DEM
    :param north_up: whether the top of the image corresponds to North or not
    :param footprints: a GeoDataFrame containing image footprints and an ID field with image names. If not
        provided, will attempt to download from USGS.
    :param mean_el: the mean surface elevation covered by the image. If None, uses DEM and footprint to
        calculate the value.
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
        optical_bar_cam(fn_img, flavor, 'samp_aft.tsai', fprint, scan_res=scan_res, mean_el=mean_el)
        fn_samp = 'samp_aft.tsai'
    else:
        optical_bar_cam(fn_img, flavor, 'samp_for.tsai', fprint, scan_res=scan_res, mean_el=mean_el)
        fn_samp = 'samp_for.tsai'

    coords = _stanrogers(fprint, north_up)

    cl_args = ['cam_gen', '--sample-file', fn_samp, '--camera-type', 'opticalbar',
               '--lon-lat-values', '  '.join([f'{x} {y}' for x, y in coords]), fn_img,
               '--reference-dem', fn_dem, '--refine-camera', '-o', fn_img.replace('.tif', '.tsai')]

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


def bundle_adjust_from_gcp(fn_img: str, fn_cam: str, fn_out: str, fn_gcp: str) -> None:
    """
    Use an ASP GCP file to refine a camera position

    :param str fn_img: the filename of the image to generate a camera for
    :param str fn_cam: the camera filename to use for refinement
    :param str fn_out: the output folder and prefix to write the updated camera to
    :param str fn_gcp: the GCP filename
    """
    cl_args = ['bundle_adjust', fn_img, fn_cam, fn_gcp, '-t', 'opticalbar', '--inline-adjustments',
               '--num-passes', '1', '--camera-weight', '0', '--ip-detect-method', '1', '-o', fn_out,
               '--max-iterations', '30', '--fix-gcp-xyz']

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


def write_asp_gcp(fn_gcp: Union[str, Path], gcp_df: gpd.GeoDataFrame,
                  gcp_list: Union[None, list] = None, imlist: Union[None, list] = None,
                  scale: int = 1, singles: bool = True, meas: Union[None, pd.DataFrame] = None,
                  headers: Union[None, list[str]] = None) -> None:
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
    """
    with open(fn_gcp, 'w') as f:
        if gcp_list is not None:
            for gcp in gcp_list:
                _gcp = gcp_df.loc[gcp]
                lon, lat = _gcp.geometry.x, _gcp.geometry.y

                out_gcp = ','.join([gcp.strip('GCP'), str(lat), str(lon), str(_gcp.elevation), '1.0', '1.0', '1.0'])

                if not singles:
                    if all([gcp in meas.loc[meas.image == img]['name'].values for img in imlist]):
                        for img in sorted(imlist):
                            row, col = meas.loc[(meas.image == img) & (meas.name == gcp), ['i', 'j']].values[0]
                            out_gcp += ',' + ','.join([img, str(col / scale), str(row / scale), '1.0', '1.0'])
                        print(out_gcp, file=f)
                else:
                    for img in sorted(imlist):
                        try:
                            row, col = meas.loc[(meas.image == img) & (meas.name == gcp), ['i', 'j']].values[0]
                            out_gcp += ',' + ','.join([img, str(col / scale), str(row / scale), '1.0', '1.0'])
                        except IndexError as e:
                            continue
                    print(out_gcp, file=f)
        else:
            if headers is not None:
                for header in headers:
                    print(header, file=f)

            for gcp in gcp_df.drop(columns=['geometry']).itertuples():
                print(' '.join([str(c) for c in list(gcp)]), file=f)


def _parse_gcp(fn_gcp):
    cols = ['id', 'lat', 'lon', 'height_above_datum', 'sigma_x', 'sigma_y', 'sigma_z']
    img_headers = ['image_name', 'pixel_x', 'pixel_y', 'sx_px', 'sy_px']


    with open(fn_gcp, 'r') as f:
        all_lines = f.readlines()

        crs_wkt = all_lines[0].strip().replace('# ', '')

        delim = ' '
        ncols = len(all_lines[2].strip().split(delim))

        if ncols < 2:
            delim = ','
            ncols = len(all_lines[2].strip().split(delim))

    nimg = 0
    while len(cols) < ncols:
        nimg += 1
        cols += ['.'.join([c, str(nimg)]) for c in img_headers]

    cols += ['blank']
    df = pd.read_csv(fn_gcp, skiprows=2, delimiter=delim, names=cols).set_index('id')
    del df['blank']

    return gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat, crs='epsg:4326'))


def mask_asp_dem_gcps(fn_gcp: Union[str, Path],
                      fn_inc: Union[None, str, Path, list[str], list[Path]] = None,
                      fn_exc: Union[None, str, Path, list[str], list[Path]] = None) -> None:
    """
    Filter/mask GCPs from a .gcp file, using inclusion/exclusion masks. Inclusion masks indicate areas where points
    inside should be included as stable terrain; exclusion masks indicate areas where points inside should be excluded.
    At least one of fn_inc, fn_exc must be set in order to filter.

    TODO: implement multiple masks (WIP)
    :param fn_gcp: the .gcp file to filter
    :param fn_inc: the inclusion mask(s) to use, pointing to a vector file format
    :param fn_exc: the exclusion mask(s) to use, pointing to a vector file format
    """
    assert not all(el is None for el in [fn_inc, fn_exc]), "at least one mask has to be set!"

    gcps = _parse_gcp(fn_gcp)
    with open(fn_gcp, 'r') as f:
        headers = [l.strip() for l in f.readlines() if l.strip()[0] == '#']

    if fn_inc is not None:
        inc_mask = gpd.read_file(fn_inc)
        within_inc = inc_mask.sindex.query(gcps.geometry, predicate='intersects')[0]
    else:
        within_inc = gcps.index[~gcps.index.isin(gcps.index)]

    if fn_exc is not None:
        exc_mask = gpd.read_file(fn_exc)
        within_exc = exc_mask.sindex.query(gcps.geometry, predicate='intersects')[0]
    else:
        within_exc = gcps.index

    valid = gcps.index[(gcps.index.isin(within_inc)) & (~gcps.index.isin(within_exc))]

    write_asp_gcp(fn_gcp, gcps.loc[valid], headers=headers)
