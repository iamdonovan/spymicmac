"""
sPyMicMac.usgs is a collection of tools for interfacing with USGS Earth Explorer.
"""
import os
import netrc
import geopandas as gpd
from shapely.geometry.polygon import Polygon
from usgs import api


def get_login_creds():
    return netrc.netrc(os.path.expanduser('~/.netrc'))


def read_coords(result):
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


def get_usgs_footprints(imglist, dataset='DECLASSII'):
    # air photos: 'AERIAL_COMBIN'
    creds = get_login_creds()

    gdf = gpd.GeoDataFrame()

    user, _, pwd = creds.authenticators('earthexplorer.usgs.gov')
    login = api.login(user, pwd)
    del user, pwd, creds

    if login['errorCode'] is not None:
        print('Error logging in to USGS EROS.')
        raise login['error']
    else:
        search_results = api.metadata(dataset,
                                      node='EE',
                                      entityids=imglist)
        for i, result in enumerate(search_results['data']):
            # coords = result['spatialFootprint']['coordinates'][0]
            coords = read_coords(result)
            poly = Polygon(coords)

            gdf.loc[i, 'geometry'] = poly
            gdf.loc[i, 'ID'] = result['entityId']

        gdf.crs = {'init': 'epsg:4326'}

        return gdf
