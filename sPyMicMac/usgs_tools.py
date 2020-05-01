import os
import netrc
import geopandas as gpd
from shapely.geometry.polygon import Polygon
from usgs import api


def get_login_creds():
    return netrc.netrc(os.path.expanduser('~/.netrc'))


def get_usgs_footprints(imglist, dataset='DECLASSII'):
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
            coords = result['spatialFootprint']['coordinates'][0]
            poly = Polygon(coords)

            gdf.loc[i, 'geometry'] = poly
            gdf.loc[i, 'ID'] = result['entityId']

        gdf.crs = {'init': 'epsg:4326'}

        return gdf
