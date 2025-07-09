"""
spymicmac.asp is a collection of tools for interfacing with Ames Stereo Pipeline
"""
import os
from pathlib import Path
import subprocess
import numpy as np
import pyproj
import geoutils as gu
import geopandas as gpd
from osgeo import gdal
from shapely.ops import split, orient
from shapely.geometry import LineString, Point, Polygon
from . import data, declass, micmac
from typing import Union


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
                    mean_el: Union[float, int] = 1000) -> None:
    """
    Generate a sample ASP camera file for a KH-4 Optical Bar camera.

    :param fn_img: the filename of the image. Used to read the image size, and determine whether the image is from
        the aft or forward camera.
    :param flavor: what type of camera the image came from - currently either KH4 or KH9
    :param out_name: the filename to write the camera file to.
    :param fprint: an optional image, footprint used to estimate the initial camera position
    :param scan_res: the image scanning resolution, in m per pixel (e.g., 7 microns -> 7.0e-6)
    :param mean_el: the mean elevation covered by the image
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
        params['motion_comp'] = 1

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

        print(f'speed = {params["speed"]}', file=f)
        print('mean_earth_radius = 6371000', file=f)
        # need a better value than this
        print(f"mean_surface_elevation = {mean_el}", file=f)
        print(f"motion_compensation_factor = {params['motion_comp']}", file=f)

        if _isaft(fn_img):
            print('scan_dir = left', file=f)
        else:
            print('scan_dir = right', file=f)


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

    with open(outname, 'w') as f:
        for gcp in gcp_list:
            _gcp = gcps.loc[gcp]
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
