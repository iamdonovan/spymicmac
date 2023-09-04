"""
spymicmac.data is a collection of tools for handling external datasets
"""
import os
import sys
import urllib
import netrc
import zipfile
import tarfile
import geopandas as gpd
import numpy as np
from pathlib import Path
from glob import glob
import pyproj
from osgeo import gdal
from shapely.geometry.polygon import Polygon
from usgs import api, USGSAuthExpiredError
from pybob.GeoImg import GeoImg


def _check_data_dir():
    if not _data_dir().exists():
        os.makedirs(_data_dir(), exist_ok=True)


def _data_dir():
    return Path(sys.prefix, 'share', 'spymicmac')


def _get_login_creds():
    """
    Read a user's .netrc file and return the login credentials.
    """

    # first, check whether .netrc exists in the home directory
    if os.path.exists(os.path.expanduser('~/.netrc')):
        return netrc.netrc(os.path.expanduser('~/.netrc'))
    # next, check whether _netrc exists
    elif os.path.exists(os.path.expanduser('~/_netrc')):
        return netrc.netrc(os.path.expanduser('~/_netrc'))
    # if we don't have ~/.netrc or ~/_netrc, then we can't log in.
    else:
        raise FileExistsError("Please ensure that you have a .netrc file stored in your home directory.")


def _authenticate():
    """
    Use the credentials stored in the user's .netrc file to authenticate the user on earthexplorer.usgs.gov
    """
    creds = _get_login_creds()
    user, _, pwd = creds.authenticators('earthexplorer.usgs.gov')

    try:
        login = api.login(user, pwd)
    except USGSAuthExpiredError as e:
        print('API key has expired. Attempting to remove .usgs from home directory.')
        os.remove(os.path.expanduser('~/.usgs'))

        login = api.login(user, pwd)

    del user, pwd, creds

    return login


def _read_coords(result):
    """
    Parse a search result returned from USGS and create a list of coordinates for the image footprint.
    """
    corner_names = ['NW', 'NE', 'SE', 'SW']
    corner_fields = [d for d in result['metadataFields'] if 'Corner' in d['fieldName'] and 'dec' in d['fieldName']]
    corner_dict = dict()

    for field in corner_fields:
        corner_dict[field['fieldName']] = float(field['value'])
    coords = []

    for corner in corner_names:
        coords.append((corner_dict['{} Corner Long dec'.format(corner)],
                       corner_dict['{} Corner Lat dec'.format(corner)]))

    return coords


def get_usgs_footprints(imlist, dataset='DECLASSII'):
    """
    Search for a list of images on USGS Earth Explorer. Note that the image names should be the USGS entity ID (e.g.,
    AR5840034159994 rather than AR5840034159994.tif).

    :param list imlist: a list of image names.
    :param str dataset: the USGS dataset name to search (default: DECLASSII).

    :return: **footprints** (*GeoDataFrame*) -- a GeoDataFrame of image footprints.
    """
    # air photos: 'AERIAL_COMBIN'
    gdf = gpd.GeoDataFrame()

    login = _authenticate()

    if login['errorCode'] is not None:
        print('Error logging in to USGS EROS.')
        raise login['error']
    else:
        search_results = api.scene_metadata(dataset, entity_id=imlist)
        for ii, result in enumerate(search_results['data']):
            poly = Polygon(result['spatialCoverage']['coordinates'][0])

            gdf.loc[ii, 'geometry'] = poly
            gdf.loc[ii, 'ID'] = result['entityId']

        gdf.set_crs(epsg=4326, inplace=True)

        return gdf


def landsat_to_gdf(results):
    """
    Convert USGS Landsat search results to a GeoDataFrame

    :param list results: a list of results (i.e., the 'data' value returned by api.scene_metadata)
    :return: **meta_gdf** (*geopandas.GeoDataFrame*) - a GeoDataFrame of search results
    """
    meta_gdf = gpd.GeoDataFrame()

    for ind, res in enumerate(results):
        keys = [f['fieldName'] for f in res['metadata']]
        values = [f['value'] for f in res['metadata']]
        metadata = dict(zip(keys, values))

        meta_gdf.loc[ind, 'ID'] = metadata['Landsat Product Identifier L1']
        meta_gdf.loc[ind, 'geometry'] = Polygon(res['spatialCoverage']['coordinates'][0])
        meta_gdf.loc[ind, 'acq_date'] = metadata['Date Acquired']
        meta_gdf.loc[ind, 'land_cc'] = float(metadata['Land Cloud Cover'])
        meta_gdf.loc[ind, 'scene_cc'] = float(metadata['Scene Cloud Cover L1'])
        meta_gdf.loc[ind, 'day_night'] = metadata['Day/Night Indicator']
        # meta_gdf.loc[ind, 'rmse_x'] = float(metadata['Geometric RMSE Model X'])
        # meta_gdf.loc[ind, 'rmse_y'] = float(metadata['Geometric RMSE Model Y'])

    return meta_gdf.set_crs(epsg=4326)


def _clean_imlist(imlist, globstr):
    if imlist is None:
        imlist = glob(globstr)
        imlist.sort()

    return [im.split('OIS-Reech_')[-1].split('.tif')[0] for im in imlist]


def download_cop30_vrt(imlist=None, footprints=None, imgsource='DECLASSII', globstr='OIS*.tif'):
    """
    Create a VRT using Copernicus 30m DSM tiles that intersect image footprints. Creates Copernicus_DSM.vrt using files
    downloaded to cop30_dem/ within the current directory.

    :param list imlist: a list of image filenames. If None, uses globstr to search for images in the current directory.
    :param GeoDataFrame footprints: a GeoDataFrame of image footprints. If None, uses spymicmac.usgs.get_usgs_footprints
        to download footprints based on imlist.
    :param str imgsource: the EE Dataset name for the images (default: DECLASSII)
    :param str globstr: the search string to use to find images in the current directory.
    """

    clean_imlist = _clean_imlist(imlist, globstr)

    if footprints is None:
        footprints = get_usgs_footprints(clean_imlist, dataset=imgsource)

    # now, get the envelope
    xmin, ymin, xmax, ymax = footprints.unary_union.bounds

    lat_min = int(np.floor(ymin))
    lat_max = int(np.ceil(ymax))

    lon_min = int(np.floor(xmin))
    lon_max = int(np.ceil(xmax))

    lats = np.arange(lat_min, lat_max)
    lons = np.arange(lon_min, lon_max)

    Lons, Lats = np.meshgrid(lons, lats)

    pairs = list(zip(Lons.flatten(), Lats.flatten()))

    tiles = []
    for pair in pairs:
        tiles.append(_format_cop30(_lat_prefix(pair[1]) + '{:02d}'.format(abs(pair[1])),
                                   _lon_prefix(pair[0]) + '{:03d}'.format(abs(pair[0]))))

    # now, download the tiles using boto3
    os.makedirs('cop30_dem', exist_ok=True)

    for tile in tiles:
        this_url = '/'.join(['https://copernicus-dem-30m.s3.amazonaws.com', tile, tile + '.tif'])
        try:
            urllib.request.urlretrieve(this_url, os.path.join('cop30_dem', tile + '.tif'))
        except urllib.error.HTTPError:
            print(f'No tile found for {tile}')

    filelist = glob(os.path.join('cop30_dem', '*DEM.tif'))
    out_vrt = gdal.BuildVRT('Copernicus_DSM.vrt', filelist, srcNodata=0)
    out_vrt = None


def _lon_prefix(lon):
    if lon < 0:
        return 'W'
    else:
        return 'E'


def _lat_prefix(lat):
    if lat < 0:
        return 'S'
    else:
        return 'N'


def _format_cop30(lat, lon):
    return f'Copernicus_DSM_COG_10_{lat}_00_{lon}_00_DEM'


def to_wgs84_ellipsoid(fn_dem):
    """
    Convert a DEM to WGS84 ellipsoid heights, using the EGM08 Geoid.

    :param str fn_dem: the filename of the DEM to convert.
    """
    proj_data = pyproj.datadir.get_data_dir()

    if not Path(proj_data, 'egm08_25.gtx').exists():
        print('Downloading egm08_25.gtx from osgeo.org')
        this_url = 'https://download.osgeo.org/proj/vdatum/egm08_25/egm08_25.gtx'
        urllib.request.urlretrieve(this_url, Path(proj_data, 'egm08_25.gtx'))

    dem = GeoImg(fn_dem)
    geoid = GeoImg(str(Path(proj_data, 'egm08_25.gtx'))).reproject(dem)

    ell = dem.copy(new_raster=(dem.img + geoid.img))
    ell.write(os.path.splitext(fn_dem)[0] + '_ell.tif')


def download_arcticdem_mosaic(imlist=None, footprints=None, imgsource='DECLASSII', globstr='OIS*.tif', res='2m'):
    """
    Download the ArcticDEM v3.0 Mosaic tiles that intersect image footprints. Downloads .tar.gz files to
    arctic_dem, extracts each DEM, and creates ArcticDEM.vrt in the current directory.

    :param list imlist: a list of image filenames. If None, uses globstr to search for images in the current directory.
    :param GeoDataFrame footprints: a GeoDataFrame of image footprints. If None, uses spymicmac.usgs.get_usgs_footprints
        to download footprints based on imlist.
    :param str imgsource: the EE Dataset name for the images (default: DECLASSII)
    :param str globstr: the search string to use to find images in the current directory.
    :param str res: the DEM resolution to download. Options are 2m, 10m, or 32m (default: 2m)
    """
    assert res in ['2m', '10m', '32m'], "res must be one of 2m, 10m, or 32m"

    clean_imlist = _clean_imlist(imlist, globstr)

    if footprints is None:
        footprints = get_usgs_footprints(clean_imlist, dataset=imgsource)

    os.makedirs('arctic_dem', exist_ok=True)

    arcticdem_tiles = _arcticdem_shp(res=res)
    intersects = arcticdem_tiles.intersects(footprints.to_crs(arcticdem_tiles.crs).unary_union)

    selection = arcticdem_tiles.loc[intersects].reset_index(drop=True)
    for ind, row in selection.iterrows():
        this_path = Path('arctic_dem', row['dem_id'] + '.tar.gz')
        print('Downloading', row['dem_id'], f'({ind+1}/{selection.shape[0]})')
        urllib.request.urlretrieve(row['fileurl'], this_path)

    tarlist = []
    for tarball in tarlist:
        _unpack_adem(tarball)

    filelist = glob(os.path.join('arctic_dem', '*reg_dem.tif'))
    out_vrt = gdal.BuildVRT('ArcticDEM.vrt', filelist)
    out_vrt = None


def _arcticdem_shp(res='2m'):
    _check_data_dir()

    # latest version is 4.1 - may need to update with future releases
    fn_shp = Path(_data_dir(), f'ArcticDEM_Mosaic_Index_v4_1_{res}.shp')

    if not fn_shp.exists():
        print('Downloading latest ArcticDEM Mosaic Tile Index from data.pgc.umn.edu')
        zip_url = 'https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/indexes/ArcticDEM_Strip_Index_latest_shp.zip'
        zip_path = Path(_data_dir(), 'ArcticDEM_Mosaic_Index_latest_shp.zip')
        urllib.request.urlretrieve(zip_url, zip_path)

        with zipfile.ZipFile(Path(_data_dir(), 'ArcticDEM_Mosaic_Index_latest_shp.zip'), 'r') as zip_ref:
            zip_ref.extractall(_data_dir())

        os.remove(zip_path)

    return gpd.read_file(fn_shp)


def _unpack_adem(tarball):
    with tarfile.open(Path('arctic_dem', tarball), 'r') as tfile:
        dem = tfile.getmember(tarball.replace('.tar.gz', '_reg_dem.tif'))
        dem.name = Path('arctic_dem', dem.name)  # will extract to arctic_dem
        tfile.extract(dem)
