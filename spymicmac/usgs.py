"""
sPyMicMac.usgs is a collection of tools for interfacing with USGS Earth Explorer.
"""
import os
import netrc
import geopandas as gpd
from shapely.geometry.polygon import Polygon
from usgs import api, USGSAuthExpiredError


def get_login_creds():
    """
    Read a user's .netrc file and return the login credentials.

    :return:
        - **creds** -- the netrc.netrc credentials.
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


def authenticate():
    """
    Use the credentials stored in the user's .netrc file to authenticate the user on earthexplorer.usgs.gov

    :return:
        -- **login** (*dict*) -- the response from the login attempt
    """
    creds = get_login_creds()
    user, _, pwd = creds.authenticators('earthexplorer.usgs.gov')

    try:
        login = api.login(user, pwd)
    except USGSAuthExpiredError as e:
        print('API key has expired. Attempting to remove .usgs from home directory.')
        os.remove(os.path.expanduser('~/.usgs'))

        login = api.login(user, pwd)

    del user, pwd, creds

    return login


def read_coords(result):
    """
    Parse a search result returned from USGS and create a list of coordinates for the image footprint.

    :param dict result: the USGS search result
    :return:
        -- **coords** (*list*) -- a list of coordinates
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

    :return:
        -- **gdf** (*GeoDataFrame*) -- a GeoDataFrame of image footprints.
    """
    # air photos: 'AERIAL_COMBIN'
    gdf = gpd.GeoDataFrame()

    login = authenticate()

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
    :return:
        - **meta_gdf** (*geopandas.GeoDataFrame*) - a GeoDataFrame of search results
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

