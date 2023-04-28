"""
spymicmac.ee_tools is a collection of tools for getting GEE data for processing KH-9 Hexagon imagery.
"""
import ee


ee.Initialize()


def get_landsat_collections(filt_geom=None, sensors='all', tier='all'):
    """
    Return a collection of raw Landsat images filtered by an (optional) geometry, given a list of sensor names and
    tier levels. Default behavior is to select all sensors and all tiers.

    Valid sensor names: LC08, LE07, LT05, LT04, LM05, LM04, LM03, LM02, LM01

    Valid tier levels: T1, T2. Check Landsat docs for details.

    :param filt_geom: A geometry object to filter images by. See EE docs for more details.
    :param sensors: Landsat sensor(s) to select images from.
    :param tier: Processing tier(s) to select.
    :return:
    """
    assert tier.lower() in ['all', 't1', 't2'], "tier is not one of 'all', 't1', 't2': {}".format(tier)
    if tier == 'all':
        tiers = ['T1', 'T2']
    else:
        tiers = list(tier.upper())

    allSensors = ['LC08', 'LE07', 'LT05', 'LT04', 'LM05', 'LM04', 'LM03', 'LM02', 'LM01']
    if isinstance(sensors, str) and sensors.lower() == 'all':
        getSensors = allSensors
    else:
        getSensors = []
        if isinstance(sensors, list):
            for sens in sensors:
                assert sens.upper() in allSensors, "{} is not a recognized Landsat sensor name".format(sens)
                getSensors.append(sens)
        else:
            assert isinstance(sensors, str), "sensors is not a list of strings, or a single string: {}".format(sensors)
            assert sensors.upper() in allSensors, "{} is not a recognized Landsat sensor name".format(sensors)
            getSensors = [sensors]

    outColls = []
    for sens in getSensors:
        for t in tiers:
            if filt_geom is not None:
                outColls.append(ee.ImageCollection('LANDSAT/{}/C01/{}'.format(sens, t)).filterBounds(filt_geom))
            else:
                outColls.append(ee.ImageCollection('LANDSAT/{}/C01/{}'.format(sens, t)))

    if len(outColls) == 1:
        return outColls[0]
    else:
        _out = outColls[0]
        for coll in outColls[1:]:
            _out.merge(coll)
        return _out


def get_srtm(ver='usgs', filt_geom=None):
    '''
    Return the SRTM DEM, clipped to an optional geometry.

    :param ver: SRTM version to select, one of 'usgs', 'cgiar'
    :param filt_geom: A geometry object to clip images to

    :return:
    '''
    assert ver.lower() in ['usgs', 'cgiar'], "srtm version is not one of 'usgs', 'cgiar': {}".format(ver)

    verDict = {'usgs': 'USGS/SRTMGL1_003', 'cgiar': 'CGIAR/SRTM90_V4'}
    _ver = verDict[ver]

    if filt_geom is not None:
        return ee.Image(_ver).clip(filt_geom)
    else:
        return ee.Image(_ver)


def get_alos_dem(filt_geom=None):
    '''
    Return the ALOS DEM, clipped to an optional geometry.

    :param filt_geom: A geometry object to clip images to

    :return:
    '''

    if filt_geom is not None:
        return ee.Image('JAXA/ALOS/AW3D30/V2_2').select('AVE_DSM').clip(filt_geom)
    else:
        return ee.Image('JAXA/ALOS/AW3D30/V2_2').select('AVE_DSM')



def get_s2_collection():
    pass


def get_cloudfree_mosaic():
    pass
