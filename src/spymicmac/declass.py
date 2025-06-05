"""
spymicmac.declass is a module designed for handling metadata related to USGS declassified images.
"""
import os
import numpy as np
from typing import Union


# default initial parameters for panoramic cameras
sample_params = {
    'KH4': {'f': 0.61, 'tilt': np.deg2rad(15), 'scan_time': 0.36, 'speed': 7700},
    'KH9': {'f': 1.5, 'tilt': np.deg2rad(10), 'scan_time': 0.7, 'speed': 8000}
}


# dataset names for searching from USGS EarthExplorer API
usgs_datasets = {
    'KH4': 'corona2',
    'KH9': 'declassiii'
}


def _mission(fn_img: str) -> int:
    """
    Extract the mission number from a USGS declassified dataset filename.
    """
    if 'OIS-Reech_' in fn_img:
        fn_img = fn_img.split('OIS-Reech_')[-1]
    fn_img = os.path.splitext(fn_img)[0]

    return int(fn_img.split('-')[0][-4:])


def _mission_to_camera(num: int) -> Union[str, None]:
    """
    Given a mission number extracted from a USGS declassified dataset filename, return the camera type associated
    with that number.

    """
    if num in range(9031, 9063):
        return 'KH4'
    elif num in range(1001, 1053):
        return 'KH4A'
    elif num in range(1101, 1118):
        return 'KH4B'
    elif num > 1200:
        return 'KH9'
    else:
        return None


def get_declass_camera(fn_img: str) -> Union[str, None]:
    """
    Given a USGS Entity ID for a declassified image, returns the camera the image corresponds to, based on the
    mission number extracted from the Entity ID.

    Examples (Entity ID -> Mission Num -> Camera):
        - DS1022-2075DA023 -> 1022 -> KH4A
        - DS1105-1122DA131 -> 1105 -> KH4B
        - DZB1214-500252L001001 -> 1214 -> KH9

    :param fn_img: the image filename or USGS Entity ID
    :returns: the mission/camera designation for the given Entity ID, or None if the mission number does not match
    """
    mission = _mission(fn_img)
    return _mission_to_camera(mission)
