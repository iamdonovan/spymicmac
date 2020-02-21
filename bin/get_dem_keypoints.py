import os
from pybob.GeoImg import GeoImg 
import cv2 
import numpy as np 
import matplotlib.pyplot as plt 
from pybob.image_tools import hillshade 
from skimage.io import imsave, imread 
import geopandas as gpd 
from glob import glob
from pymmaster.mmaster_tools import get_track_angle
from shapely.affinity import rotate
from shapely.ops import cascaded_union


def _argparser():
    pass


footprints = gpd.read_file('/uio/kant/geo-geohyd-u1/robertwm/data/elevations/hexagon/footprints/HexagonFootprints.shp')
ext_dem = GeoImg('')

imlist = glob('OIS-Reech*.tif')
imgs = [os.path.splitext(im)[0].split('OIS-Reech_')[-1] for im in imlist]

fprints = footprints[footprints.ID.isin(imgs)].copy()
fprints.to_crs({'init': 'epsg:{}'.format(ext_dem.epsg)}, inplace=True)
fprints.sort_values('ID', inplace=True)

overlap = cascaded_union(fprints.geometry.values[1:-1])
outline = overlap.buffer(1000).minimum_rotated_rectangle

angle = get_track_angle(outline, 0)

# TODO: crop ext_dem to outline using gdal/ogr

rot_dem = imutils.rotate_bound(crop_dem.img, -angle)
# TODO: check how I did the rotating/cropping in hexagon_gcps script - once we have points,
# have to rotate them back to reality and make sure we have the right coordinates

# now, divide the relative dem into chunks, select keypoints
chunk = rel_dem.img[:500, :500]
# try both hillshade and highpass_filter on chunk
# remember that hillshade needs to somehow scale

orb = cv2.ORB_create()
kp = orb.detect(chunk.astype(np.uint8), None)
# find best matches in kp using cv2
# save best matches

# check - find residuals from campari, remove anything over X ground units? check residual files.
