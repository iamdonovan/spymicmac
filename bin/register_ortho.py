#!/usr/bin/env python
import argparse
import os
import subprocess
import shutil
from itertools import chain
import numpy as np
import matplotlib.pyplot as plt
import cv2
import gdal
import pandas as pd
import geopandas as gpd
from PIL import Image
import lxml.etree as etree
import lxml.builder as builder
import xml.etree.ElementTree as ET
from glob import glob
from shapely.geometry.point import Point
from skimage.io import imread, imsave
from skimage import exposure
from skimage.feature import peak_local_max
from skimage.measure import ransac
from skimage.transform import EuclideanTransform, SimilarityTransform, AffineTransform, warp, estimate_transform
from pybob.bob_tools import mkdir_p
from pybob.ddem_tools import nmad
from pybob.image_tools import create_mask_from_shapefile
from pybob.GeoImg import GeoImg
import sPyMicMac.image_tools as imtools
import sPyMicMac.micmac_tools as mmtools
from sPyMicMac.usgs_tools import get_usgs_footprints


def run_bascule(in_gcps, outdir, img_pattern, sub, ori):
    subprocess.Popen(['mm3d', 'GCPBascule', img_pattern, ori,
                      'TerrainRelAuto{}'.format(sub),
                      os.path.join(outdir, 'AutoGCPs{}.xml'.format(sub)),
                      os.path.join(outdir, 'AutoMeasures{}-S2D.xml'.format(sub))]).wait()

    out_gcps = mmtools.get_bascule_residuals(os.path.join('Ori-TerrainRelAuto{}'.format(sub),
                                                          'Result-GCP-Bascule.xml'), in_gcps)
    return out_gcps


def run_campari(in_gcps, outdir, img_pattern, sub, dx, ortho_res):
    subprocess.Popen(['mm3d', 'Campari', img_pattern,
                      'TerrainRelAuto{}'.format(sub),
                      'TerrainFirstPass{}'.format(sub),
                      'GCP=[{},{},{},{}]'.format(os.path.join(outdir, 'AutoGCPs{}.xml'.format(sub)),
                                                 np.abs(dx),
                                                 os.path.join(outdir, 'AutoMeasures{}-S2D.xml'.format(sub)),
                                                 np.abs(dx / ortho_res)),
                      'SH=Homol', 'AllFree=1']).wait()

    out_gcps = mmtools.get_campari_residuals('Ori-TerrainFirstPass{}/Residus.xml'.format(subscript), in_gcps)
    out_gcps.dropna(inplace=True)  # sometimes, campari can return no information for a gcp
    return out_gcps


def save_gcps(in_gcps, outdir, utmstr, sub):
    in_gcps.to_file(os.path.join(outdir, 'AutoGCPs{}.shp'.format(sub)))
    mmtools.write_auto_gcps(in_gcps, sub, outdir, utmstr)
    subprocess.Popen(['mm3d', 'GCPConvert', 'AppInFile',
                      os.path.join(outdir, 'AutoGCPs{}.txt'.format(sub))]).wait()

    auto_root = ET.parse(os.path.join(outdir, 'AutoMeasures{}-S2D.xml'.format(sub))).getroot()
    for im in auto_root.findall('MesureAppuiFlottant1Im'):
        for pt in im.findall('OneMesureAF1I'):
            if pt.find('NamePt').text not in gcps.id.values:
                im.remove(pt)

    # save AutoMeasures
    out_xml = ET.ElementTree(auto_root)
    out_xml.write(os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript)),
                  encoding="utf-8", xml_declaration=True)


def sliding_window_filter(img_shape, pts_df, winsize, stepsize=None, mindist=2000, how='residual', is_ascending=True):
    if stepsize is None:
        stepsize = winsize / 2

    _out_inds = []
    _out_pts = []

    for x_ind in np.arange(stepsize, img_shape[1], winsize):
        for y_ind in np.arange(stepsize, img_shape[0], winsize):
            min_x = x_ind - winsize / 2
            max_x = x_ind + winsize / 2
            min_y = y_ind - winsize / 2
            max_y = y_ind + winsize / 2
            samp_ = pts_df.loc[np.logical_and.reduce([pts_df.orig_i > min_x,
                                                      pts_df.orig_i < max_x,
                                                      pts_df.orig_j > min_y,
                                                      pts_df.orig_j < max_y])].copy()
            if samp_.shape[0] == 0:
                continue
            # only take the above-average z_corr values
            samp_ = samp_[samp_.z_corr >= samp_.z_corr.quantile(0.5)]
            # make sure we get the best residual
            samp_.sort_values(how, ascending=is_ascending, inplace=True)
            if len(_out_inds) == 0:
                best_ind = samp_.index[0]
                best_pt = Point(samp_.loc[best_ind, ['orig_j', 'orig_i']].values)

                _out_inds.append(best_ind)
                _out_pts.append(best_pt)
            else:
                for _ind, _row in samp_.iterrows():
                    this_pt = Point(_row[['orig_j', 'orig_i']].values)
                    this_min_dist = np.array([this_pt.distance(pt) for pt in _out_pts]).min()
                    if this_min_dist > mindist:
                        _out_inds.append(_ind)
                        _out_pts.append(this_pt)

    return np.array(_out_inds)


def get_img_grid(shape, spacing): 
    _i = np.arange(0, shape[0], spacing) 
    _j = np.arange(0, shape[1], spacing) 
 
    J, I = np.meshgrid(_j, _i) 
    return np.array([I.reshape(-1), J.reshape(-1)]).T


def _argparser():
    _parser = argparse.ArgumentParser(description="", formatter_class=argparse.RawDescriptionHelpFormatter)

    _parser.add_argument('ortho', action='store', type=str, help='non-referenced orthophoto mosaic')
    _parser.add_argument('master', action='store', type=str, help='georeferenced satellite image')
    _parser.add_argument('dem', action='store', type=str, help='dem')
    _parser.add_argument('rel_dem', action='store', type=str, help='relative dem corresponding to ortho')
    _parser.add_argument('-glacmask', action='store', type=str, default=None,
                         help='path to shapefile of glacier outlines')
    _parser.add_argument('-landmask', action='store', type=str, default=None,
                         help='path to shapefile of land outlines')
    _parser.add_argument('-footprints', action='store', type=str, default=None,
                         help='approximate outline of hexagon images to crop satellite image to')
    _parser.add_argument('-im_subset', action='store', type=str, default=None, nargs='+',
                         help='subset of raw images to work with (default all)')
    _parser.add_argument('-corr_thresh', action='store', type=float, default=0.5,
                         help='minimum correlation value to use for accepting a match.')
    _parser.add_argument('-tfm_pts', action='store', type=str, default=None,
                         help='CSV containing set of 4-5+ points to estimate a rough transform, in the form I,J,X,Y.')
    _parser.add_argument('-b', '--block', action='store', type=str, default=None,
                         help='Block number to use if multiple image blocks exist in directory.')
    _parser.add_argument('-ori', action='store', type=str, default='Relative',
                         help='name of orientation directory (after Ori-) [Relative]')
    _parser.add_argument('-init_res', action='store', type=int, default=400,
                         help='initial resolution to get rough transformation [400 m gsd]')
    _parser.add_argument('-ortho_res', action='store', type=float, default=8,
                         help='approx. ground sampling distance (pixel resolution) of ortho image. [8 m]')
    _parser.add_argument('-imgsource', action='store', type=str, default='DECLASSII',
                         help='USGS dataset name for images [DECLASSII]')
    _parser.add_argument('-density', action='store', type=int, default=200,
                         help='pixel spacing to look for GCPs [200]')
    return _parser


# def main():
parser = _argparser()
args = parser.parse_args()

plt.ion()
print('start.')

out_dir = 'auto_gcps'
mkdir_p(out_dir)
if args.block is not None:
    subscript = '_block{}'.format(args.block)
else:
    subscript = ''

if args.im_subset is None:
    imlist = glob('OIS*.tif')
    match_pattern = 'OIS.*.tif'
else:
    imlist = args.im_subset
    match_pattern = mmtools.get_match_pattern(imlist)

mst = GeoImg(args.master)

# assumes wgs84 utm zones (326XX for N, 327XX for S)
# there's probably a much better way to get this info...
epsg_str = str(mst.epsg)
hemi_dict = {'6': 'N', '7': 'S'}
utm_str = epsg_str[-2:] + hemi_dict[epsg_str[2]]

ortho = imread(args.ortho)
ortho_ = Image.fromarray(ortho)

lowres_ = int(args.init_res / args.ortho_res)

ortho_lowres = np.array(ortho_.resize((np.array(ortho_.size) / lowres_).astype(int), Image.LANCZOS))

mst_lowres = mst.resample(args.init_res, method=gdal.GRA_Lanczos)
mst_eq = imtools.stretch_image(mst_lowres.img, scale=(0.05, 0.95))

if args.footprints is not None:
    fmask, fprint = imtools.get_footprint_mask(args.footprints, mst_lowres, imlist, fprint_out=True)
else:
    clean_imlist = [im.split('OIS-Reech_')[-1].split('.tif')[0] for im in imlist]
    print('Attempting to get image footprints from USGS EarthExplorer.')
    _fprints = get_usgs_footprints(clean_imlist, dataset=args.imgsource)
    fmask, fprint = imtools.get_footprint_mask(_fprints, mst_lowres, imlist, fprint_out=True)

lowres_mask = imtools.make_binary_mask(mst_lowres.img, erode=3, mask_value=np.nan)
if args.landmask is not None:
    lmask = create_mask_from_shapefile(mst_lowres, args.landmask)
    lowres_mask[~lmask] = 0
if args.glacmask is not None:
    gmask = create_mask_from_shapefile(mst_lowres, args.glacmask)
    lowres_mask[gmask] = 0

lowres_mask[~fmask] = 0

if args.imgsource == 'DECLASSII':
    M = imtools.transform_from_fprint(ortho_lowres, mst_lowres, fprint)
else:
    ortho_mask = 255 * np.ones(ortho_lowres.shape, dtype=np.uint8)
    ortho_mask[ortho_lowres == 0] = 0

    kp, des, matches = imtools.get_matches(ortho_lowres.astype(np.uint8),
                                           mst_eq.astype(np.uint8),
                                           mask1=ortho_mask,
                                           mask2=lowres_mask,
                                           dense=True)
    print('{} matches found.'.format(len(matches)))

    src_pts = np.array([kp[0][m.queryIdx].pt for m in matches])
    dst_pts = np.array([kp[1][m.trainIdx].pt for m in matches])

    M, inliers = ransac((dst_pts, src_pts), EuclideanTransform,
                        min_samples=10, residual_threshold=25, max_trials=10000)
    print('{} matches used for initial transformation'.format(np.count_nonzero(inliers)))

rough_tfm = warp(ortho_lowres, M, output_shape=mst_lowres.img.shape, preserve_range=True)

lowres_mask[gmask] = 255

rough_gcps = imtools.find_grid_matches(rough_tfm, mst_lowres, lowres_mask, M, spacing=10, srcwin=25, dstwin=200)

best = np.logical_and.reduce((rough_gcps.pk_corr > rough_gcps.pk_corr.quantile(0.75),
                              rough_gcps.z_corr > rough_gcps.z_corr.quantile(0.5),
                              rough_tfm[rough_gcps.search_i, rough_gcps.search_j] > 0))

Minit, inliers = ransac((rough_gcps.loc[best, ['search_j', 'search_i']].values,
                         rough_gcps.loc[best, ['orig_j', 'orig_i']].values),
                        AffineTransform, min_samples=10, residual_threshold=10, max_trials=5000)
print('{} valid matches used for initial transformation'.format(np.count_nonzero(inliers)))

ortho_tfm = warp(ortho_lowres, Minit, output_shape=mst_lowres.img.shape, preserve_range=True, order=5)

dst_scale = mst_lowres.dx / mst.dx

fig, ax = plt.subplots(1, 2, figsize=(7, 5))
ax[0].imshow(ortho_tfm, cmap='gray')
ax[1].imshow(mst_eq, cmap='gray')

plt.savefig('initial_transformation.png', dpi=200, bbox_inches='tight')

i_ = np.arange(0, ortho_lowres.shape[0], 10)
j_ = np.arange(0, ortho_lowres.shape[1], 10)

I, J = np.meshgrid(i_, j_)

src_grd = np.array(list(zip(J.reshape(-1, 1), I.reshape(-1, 1)))).reshape(-1, 2)

dst_tfm = []
for pt in src_grd:
    dst_tfm.append(Minit.inverse(pt) * dst_scale)
dst_tfm = np.array(dst_tfm).reshape(-1, 2)

Minit_full, _ = ransac((dst_tfm, lowres_ * src_grd), AffineTransform, min_samples=3,
                       residual_threshold=1, max_trials=1000)
rough_tfm = warp(ortho, Minit_full, output_shape=mst.img.shape, preserve_range=True)

mask_full = imtools.make_binary_mask(mst.img, mask_value=np.nan)
if args.landmask is not None:
    lm = create_mask_from_shapefile(mst, args.landmask)
    mask_full[~lm] = 0
if args.glacmask is not None:
    gm = create_mask_from_shapefile(mst, args.glacmask)
    mask_full[gm] = 0
mask_full[rough_tfm == 0] = 0

# for each of these pairs (src, dst), find the precise subpixel match (or not...)
gcps = imtools.find_grid_matches(rough_tfm, mst, mask_full, Minit_full, spacing=args.density, dstwin=400)

xy = np.array([mst.ij2xy((pt[1], pt[0])) for pt in gcps[['search_j', 'search_i']].values]).reshape(-1, 2)
gcps['geometry'] = [Point(pt) for pt in xy]

# gcps = gcps[gcps.z_corr > gcps.z_corr.quantile(0.5)]

fn_tfw = args.ortho.replace('.tif', '.tfw')
with open(fn_tfw, 'r') as f:
    gt = [float(l.strip()) for l in f.readlines()]

gcps['rel_x'] = gt[4] + gcps['orig_j'].values * gt[0]  # need the original image coordinates
gcps['rel_y'] = gt[5] + gcps['orig_i'].values * gt[3]
gcps['elevation'] = 0
gcps.crs = mst.proj4

gcps.dropna(inplace=True)

print('loading dems')
dem = GeoImg(args.dem, dtype=np.float32)
rel_dem = GeoImg(args.rel_dem)

gcps.crs = {'init': 'epsg:{}'.format(mst.epsg)}
for i, row in gcps.to_crs(crs=dem.proj4).iterrows():
    gcps.loc[i, 'elevation'] = dem.raster_points([(row.geometry.x, row.geometry.y)], nsize=3, mode='linear')
    gcps.loc[i, 'el_rel'] = rel_dem.raster_points([(row.rel_x, row.rel_y)], nsize=3, mode='linear')

# drop any gcps where we don't have a DEM value or a valid match
gcps.dropna(inplace=True)

# run ransac to find the matches between the transformed image and the master image make a coherent transformation
# residual_threshold is 10 pixels to allow for some local distortions, but get rid of the big blunders
Mref, inliers_ref = ransac((gcps[['search_j', 'search_i']].values, gcps[['match_j', 'match_i']].values),
                           AffineTransform, min_samples=6, residual_threshold=20, max_trials=5000)

gcps_orig = gcps.copy()

gcps = gcps[inliers_ref]

out = sliding_window_filter([ortho.shape[1], ortho.shape[0]], gcps,
                            1000, mindist=500, how='z_corr', is_ascending=False)

gcps = gcps.loc[out]

Mfin, inliers_fin = ransac((gcps[['search_j', 'search_i']].values, gcps[['match_j', 'match_i']].values),
                           AffineTransform, min_samples=6, residual_threshold=20, max_trials=5000)
gcps = gcps[inliers_fin]

print('{} valid matches found'.format(gcps.shape[0]))

gcps.index = range(gcps.shape[0])  # make sure index corresponds to row we're writing out
gcps['id'] = ['GCP{}'.format(i) for i in range(gcps.shape[0])]
gcps.to_file(os.path.join(out_dir, 'AutoGCPs{}.shp'.format(subscript)))

print('writing AutoGCPs.txt')
mmtools.write_auto_gcps(gcps, subscript, out_dir, utm_str)

print('converting AutoGCPs.txt to AutoGCPs.xml')
subprocess.Popen(['mm3d', 'GCPConvert', 'AppInFile',
                  os.path.join(out_dir, 'AutoGCPs{}.txt'.format(subscript))]).wait()

print('writing AutoMeasures.txt')
mmtools.write_auto_mesures(gcps, subscript, out_dir)

print('running get_autogcp_locations.sh to get rough image locations for each point')
subprocess.Popen(['get_autogcp_locations.sh', 'Ori-{}'.format(args.ori),
                  os.path.join(out_dir, 'AutoMeasures{}.txt'.format(subscript))] + imlist).wait()

# print('searching for points in non-orthorectified images')
print('finding image measures')
mmtools.write_image_mesures(imlist, out_dir, subscript)

print('running mm3d GCPBascule to estimate terrain errors')
gcps = run_bascule(gcps, out_dir, match_pattern, subscript, args.ori)
gcps['res_dist'] = np.sqrt(gcps.xres**2 + gcps.yres**2)

gcps = gcps[gcps.res_dist < 5 * nmad(gcps.res_dist)]

# if args.imgsource != 'DECLASSII':
#     out = sliding_window_filter([ortho.shape[1], ortho.shape[0]], gcps, 2000, mindist=500, how='res_dist')
#     gcps = gcps.loc[out]

save_gcps(gcps, out_dir, utm_str, subscript)
gcps = run_bascule(gcps, out_dir, match_pattern, subscript, args.ori)
gcps['res_dist'] = np.sqrt(gcps.xres**2 + gcps.yres**2)

gcps = run_campari(gcps, out_dir, match_pattern, subscript, mst.dx, args.ortho_res)
gcps['camp_dist'] = np.sqrt(gcps.camp_xres**2 + gcps.camp_yres**2)

# if args.imgsource != 'DECLASSII':
#     out = sliding_window_filter([ortho.shape[1], ortho.shape[0]], gcps, 2000, mindist=1000, how='camp_dist')
#    gcps = gcps.loc[out]
niter = 0
while all([np.any(gcps.camp_res > 5 * nmad(gcps.camp_res)),
           np.any(gcps.camp_dist > 5 * nmad(gcps.camp_dist)),
           gcps.camp_res.max() > 2]) and niter < 3:
    gcps = gcps[np.logical_and.reduce((gcps.camp_res < 5 * nmad(gcps.camp_res),
                                       gcps.camp_res < gcps.camp_res.max(),
                                       gcps.z_corr > gcps.z_corr.min()))]
    save_gcps(gcps, out_dir, utm_str, subscript)
    gcps = run_bascule(gcps, out_dir, match_pattern, subscript, args.ori)
    gcps['res_dist'] = np.sqrt(gcps.xres ** 2 + gcps.yres ** 2)
    gcps = run_campari(gcps, out_dir, match_pattern, subscript, mst.dx, args.ortho_res)
    gcps['camp_dist'] = np.sqrt(gcps.camp_xres**2 + gcps.camp_yres**2)
    niter += 1

# if args.imgsource != 'DECLASSII':
#     out = sliding_window_filter([ortho.shape[1], ortho.shape[0]], gcps, 2000, mindist=3000, how='camp_dist')
#     gcps = gcps.loc[out]
#
#     save_gcps(gcps, out_dir, utm_str, subscript)
#     gcps = run_bascule(gcps, out_dir, match_pattern, subscript, args.ori)
#     gcps['res_dist'] = np.sqrt(gcps.xres ** 2 + gcps.yres ** 2)
#     gcps = run_campari(gcps, out_dir, match_pattern, subscript, mst.dx, args.ortho_res)
#     gcps['camp_dist'] = np.sqrt(gcps.camp_xres ** 2 + gcps.camp_yres ** 2)

# final write of gcps to disk.
gcps.to_file(os.path.join(out_dir, 'AutoGCPs{}.shp'.format(subscript)))

fig1 = plt.figure(figsize=(7, 5))
plt.imshow(ortho[::5, ::5], cmap='gray', extent=[0, ortho.shape[1], ortho.shape[0], 0])
plt.plot(gcps.orig_j, gcps.orig_i, 'r+')
plt.quiver(gcps.orig_j, gcps.orig_i, gcps.camp_xres, gcps.camp_yres, color='r')

plt.savefig(os.path.join(out_dir, 'relative_gcps{}.png'.format(subscript)), bbox_inches='tight', dpi=200)
plt.close(fig1)

fig2 = mst.display(sfact=10, fig=plt.figure(figsize=(7, 5)))
plt.plot(gcps.geometry.x, gcps.geometry.y, 'r+')
plt.savefig(os.path.join(out_dir, 'world_gcps{}.png'.format(subscript)), bbox_inches='tight', dpi=200)
plt.close(fig2)

print('cleaning up.')
# remove Auto-im.tif.txt, NoDist-im.tif.txt, etc. Ori-Relative-NoDist
shutil.rmtree('Ori-{}-NoDist/'.format(args.ori))

for txtfile in glob('Auto-OIS*.tif.txt') + \
               glob('NoDist-OIS*.tif.txt'):
    os.remove(txtfile)
print('end.')
