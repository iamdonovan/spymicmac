"""
spymicmac.asp is a collection of tools for interfacing with Ames Stereo Pipeline
"""
import os
import subprocess
import numpy as np
import pyproj
from osgeo import gdal
from shapely.ops import split, orient
from shapely.geometry import LineString, Point, Polygon
from spymicmac import data


sample_params = {
    'KH4': {'f': 0.61, 'tilt': np.deg2rad(15), 'scan_time': 0.5, 'speed': 7700},
    'KH9': {'f': 1.5, 'tilt': np.deg2rad(10), 'scan_time': 0.7, 'speed': 8000}
}

usgs_datasets = {
    'KH4': 'corona2',
    'KH9': 'declassiii'
}


def _isaft(fn_img):
    return os.path.splitext(fn_img)[0][-4] == 'A'


def _parse_cam(fn_cam):
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


def _init_center(fprint):
    cx, cy = fprint.centroid.x, fprint.centroid.y
    alt = 180000  # very rough estimated altitude of 180 km

    geocent = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    geodet = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

    transformer = pyproj.Transformer.from_proj(geodet, geocent)
    x, y, z = transformer.transform(cx, cy, alt)

    return x, y, z


def optical_bar_cam(fn_img, flavor, out_name, fprint=None, scan_res=7e-6):
    """
    Generate a sample ASP camera file for a KH-4 Optical Bar camera.

    :param str fn_img: the filename of the image. Used to read the image size, and determine whether the image is from
        the aft or forward camera.
    :param str flavor: what type of camera the image came from - currently either KH4 or KH9
    :param str out_name: the filename to write the camera file to.
    :param Polygon fprint: an optional image, footprint used to estimate the initial camera position
    :param float scan_res: the image scanning resolution, in m per pixel (e.g., 7 microns -> 7.0e-6)
    """
    assert flavor in sample_params.keys(), "flavor must be one of {}".format(list(sample_params.keys()))
    ds = gdal.Open(fn_img)
    width, height = ds.RasterXSize, ds.RasterYSize
    cx, cy = width / 2, height / 2

    ds = None  # close the image

    params = sample_params[flavor]

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
        print('mean_surface_elevation = 4000', file=f)
        print('motion_compensation_factor = 1.0', file=f)

        if _isaft(fn_img):
            print('scan_dir = left', file=f)
        else:
            print('scan_dir = right', file=f)


def cam_from_footprint(fn_img, flavor, scan_res, fn_dem, north_up=True, footprints=None):
    """
    Generate a camera (.tsai) file from an image footprint.

    :param str fn_img: the filename of the image to generate a camera for.
    :param str flavor: what type of camera the image came from - currently either KH4 or KH9
    :param float scan_res: the scanning resolution of the image
    :param str fn_dem: the filename of the reference DEM
    :param bool north_up: whether the top of the image corresponds to North or not (default: True)
    :param GeoDataFrame footprints: a GeoDataFrame containing image footprints and an ID field with image names. If not
        provided, will attempt to download from USGS.
    :return:
    """
    clean_name = fn_img.split('OIS-Reech_')[-1].split('.tif')[0]

    # now, get the image footprint, and use ul_corner to get the ul, ur, lr, ll coordinates
    if footprints is None:
        footprints = data.get_usgs_footprints([clean_name], dataset=usgs_datasets[flavor])
        fprint = footprints.loc[0, 'geometry']
    else:
        fprint = footprints.loc[footprints['ID'] == clean_name, 'geometry'].values[0]

    if _isaft(fn_img):
        optical_bar_cam(fn_img, flavor, 'samp_aft.tsai', fprint, scan_res=scan_res)
        fn_samp = 'samp_aft.tsai'
    else:
        optical_bar_cam(fn_img, flavor, 'samp_for.tsai', fprint, scan_res=scan_res)
        fn_samp = 'samp_for.tsai'

    coords = _stanrogers(fprint, north_up)

    cl_args = ['cam_gen', '--sample-file', fn_samp, '--camera-type', 'opticalbar',
               '--lon-lat-values', '  '.join([f'{x} {y}' for x, y in coords]), fn_img,
               '--reference-dem', fn_dem, '--refine-camera', '-o', fn_img.replace('.tif', '.tsai')]

    print(cl_args)

    p = subprocess.Popen(cl_args)
    p.wait()

    os.remove(fn_samp)


# return the upper left, upper right, lower right, lower left coordinates for an image
def _stanrogers(fprint, north_up):

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
    if north_up:
        lower, upper = split(inner, horizontal).geoms
    else:
        # check this with an actual south up image!
        upper, lower = split(inner, horizontal).geoms

    upper_left, upper_right = split(upper, vertical).geoms
    lower_right, lower_left = split(lower, vertical).geoms

    fx, fy = fprint.exterior.xy
    vertices = np.array([Point(x, y) for x, y in zip(fx[:-1], fy[:-1])])

    # upper left, upper right, lower right, lower left
    ul = vertices[[upper_left.contains(pt) for pt in vertices]][0]
    ur = vertices[[upper_right.contains(pt) for pt in vertices]][0]
    lr = vertices[[lower_right.contains(pt) for pt in vertices]][0]
    ll = vertices[[lower_left.contains(pt) for pt in vertices]][0]

    return (ul.x, ul.y), (ur.x, ur.y), (lr.x, lr.y), (ll.x, ll.y)


def bundle_adjust_from_gcp(fn_img, fn_cam, fn_out, fn_gcp):
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
