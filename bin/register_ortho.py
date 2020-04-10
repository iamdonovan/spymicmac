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
from glob import glob
from shapely.geometry.point import Point
from skimage.io import imread
from skimage import exposure
from skimage.feature import peak_local_max
from skimage.measure import ransac
from skimage.transform import EuclideanTransform, AffineTransform, warp
from pybob.bob_tools import mkdir_p
from pybob.image_tools import create_mask_from_shapefile
from pybob.GeoImg import GeoImg
import sPyMicMac.image_tools as imtools
import sPyMicMac.micmac_tools as mmtools


def sliding_window_filter(img_shape, pts_df, winsize, stepsize=None, mindist=2000):
    if stepsize is None:
        stepsize = winsize / 2

    out_inds = []
    out_pts = []

    for x_ind in np.arange(stepsize, img_shape[1], winsize):
        for y_ind in np.arange(stepsize, img_shape[0], winsize):
            min_x = x_ind - winsize / 2
            max_x = x_ind + winsize / 2
            min_y = y_ind - winsize / 2
            max_y = y_ind + winsize / 2
            samp_ = pts_df.loc[np.logical_and.reduce([pts_df.i > min_x,
                                                      pts_df.i < max_x,
                                                      pts_df.j > min_y,
                                                      pts_df.j < max_y])].copy()
            if samp_.shape[0] == 0:
                continue
            samp_.sort_values('residual', ascending=True, inplace=True)
            if len(out_inds) == 0:
                best_ind = samp_.index[0]
                best_pt = Point(samp_.loc[best_ind, ['j', 'i']].values)

                out_inds.append(best_ind)
                out_pts.append(best_pt)
            else:
                for ind, row in samp_.iterrows():
                    this_pt = Point(row[['j', 'i']].values)
                    this_min_dist = np.array([this_pt.distance(pt) for pt in out_pts]).min()
                    if this_min_dist > mindist:
                        out_inds.append(ind)
                        out_pts.append(this_pt)

    return np.array(out_inds)


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
    _parser.add_argument('-footprint', action='store', type=str, default=None,
                         help='approximate outline of hexagon images to crop satellite image to')
    _parser.add_argument('-im_subset', action='store', type=str, default=None, nargs='+',
                         help='subset of raw images to work with (default all)')
    _parser.add_argument('-corr_thresh', action='store', type=float, default=0.5,
                         help='minimum correlation value to use for accepting a match.')
    _parser.add_argument('-tfm_pts', action='store', type=str, default=None,
                         help='CSV containing set of 4-5 points to estimate a rough transform, in the form I,J,X,Y.')
    _parser.add_argument('-gcps', action='store', type=str, default=None,
                         help='Shapefile (or CSV) containing pre-selected control points.')
    _parser.add_argument('-b', '--block', action='store', type=str, default=None,
                         help='Block number to use if multiple image blocks exist in directory.')
    _parser.add_argument('-rsize', action='store', type=int, default=50,
                         help='half-size of reference search window [50 pixels]')
    _parser.add_argument('-tsize', action='store', type=int, default=600,
                         help='half-size of search window [600 pixels]')
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
else:
    imlist = args.im_subset

mst = GeoImg(args.master)

# assumes wgs84 utm zones (326XX for N, 327XX for S)
# there's probably a much better way to get this info...
epsg_str = str(mst.epsg)
hemi_dict = {'6': 'N', '7': 'S'}
utm_str = epsg_str[-2:] + hemi_dict[epsg_str[2]]

ortho = imread(args.ortho)
ortho_ = Image.fromarray(ortho)
ortho_lowres = np.array(ortho_.resize((np.array(ortho_.size) / 50).astype(int), Image.LANCZOS))

ortho_tfm, Minit, tfm_pts = imtools.get_rough_geotransform(ortho_lowres, mst, pRes=400, landmask=args.landmask)

mst_lowres = mst.resample(400, method=gdal.GRA_NearestNeighbour)
mst_eq = (255 * exposure.equalize_adapthist(mst_lowres.img.astype(np.uint16), clip_limit=0.03)).astype(np.uint8)
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

Minit_full, _ = ransac((dst_tfm, 50 * src_grd), AffineTransform, min_samples=3,
                       residual_threshold=1, max_trials=1000)

rough_tfm = warp(ortho, Minit_full, output_shape=mst.img.shape, preserve_range=True)

mask = 255 * np.ones(mst.img.shape, dtype=np.uint8)
lm = create_mask_from_shapefile(mst, args.landmask, buffer=200)
gm = create_mask_from_shapefile(mst, args.glacmask, buffer=200)

mask[np.isnan(mst.img)] = 0
mask[rough_tfm == 0] = 0
mask[~lm] = 0
mask[gm] = 0

# for each of these pairs (src, dst), find the precise subpixel match (or not...)
match_pts = []
z_corrs = []
peak_corrs = []
res_imgs = []

jj = np.arange(0, rough_tfm.shape[1], 200)
ii = np.arange(0, rough_tfm.shape[0], 200)

# I, J = np.meshgrid(ii, jj)
# search_pts = np.array(list(zip(J.reshape(-1, 1), I.reshape(-1, 1)))).reshape(-1, 2)
search_pts = []

for _i in ii:
    for _j in jj:
        search_pts.append((_j, _i))
        # for pt in search_pts:
        # if mask[pt[1], pt[0]] == 0:
        if mask[_i, _j] == 0:
            match_pts.append((-1, -1))
            z_corrs.append(np.nan)
            peak_corrs.append(np.nan)
            res_imgs.append(np.nan)
            continue

        try:
            # testchip, _, _ = imtools.make_template(rough_tfm, (pt[1], pt[0]), 40)
            # dst_chip, _, _ = imtools.make_template(mst.img, (pt[1], pt[0]), 200)
            testchip, _, _ = imtools.make_template(rough_tfm, (_i, _j), 40)
            dst_chip, _, _ = imtools.make_template(mst.img, (_i, _j), 400)

            dst_chip[np.isnan(dst_chip)] = 0

            test = np.ma.masked_values(imtools.highpass_filter(testchip), 0)
            dest = np.ma.masked_values(imtools.highpass_filter(dst_chip), 0)

            corr_res, this_i, this_j = imtools.find_gcp_match(dest.astype(np.float32), test.astype(np.float32))
            peak_corr = cv2.minMaxLoc(corr_res)[1]

            pks = peak_local_max(corr_res, min_distance=5, num_peaks=2)
            this_z_corrs = []
            for pk in pks:
                max_ = corr_res[pk[0], pk[1]]
                this_z_corrs.append((max_ - corr_res.mean()) / corr_res.std())
            dz_corr = max(this_z_corrs) / min(this_z_corrs)
            z_corr = max(this_z_corrs)

            # if the correlation peak is very high, or very unique, add it as a match
            # out_i, out_j = this_i - 200 + pt[1], this_j - 200 + pt[0]
            out_i, out_j = this_i - 400 + _i, this_j - 400 + _j
            z_corrs.append(z_corr)
            peak_corrs.append(peak_corr)
            match_pts.append([out_j, out_i])
            res_imgs.append(corr_res)
        except:
            match_pts.append((-1, -1))
            z_corrs.append(np.nan)
            peak_corrs.append(np.nan)
            res_imgs.append(np.nan)

search_pts = np.array(search_pts)
# have to find the initial, which means back-transforming search_pts with Minit_full
_src = np.dot(Minit_full.params, np.hstack([search_pts,
                                            np.ones(search_pts[:, 0].shape).reshape(-1, 1)]).T).T[:, :2]
_dst = np.array(match_pts)

xy = np.array([mst.ij2xy((pt[1], pt[0])) for pt in _src]).reshape(-1, 2)

fn_tfw = args.ortho.replace('.tif', '.tfw')
with open(fn_tfw, 'r') as f:
    gt = [float(l.strip()) for l in f.readlines()]

gcps = gpd.GeoDataFrame()
gcps['geometry'] = [Point(pt) for pt in xy]
gcps['pk_corr'] = peak_corrs
gcps['z_corr'] = z_corrs
gcps['match_j'] = _dst[:, 0]  # points matched in master image
gcps['match_i'] = _dst[:, 1]
gcps['orig_j'] = _src[:, 0]  # this should be the back-transformed search_pts
gcps['orig_i'] = _src[:, 1]
gcps['search_j'] = search_pts[:, 0]
gcps['search_i'] = search_pts[:, 1]
gcps['rel_x'] = gt[4] + _src[:, 0] * gt[0]  # need the original image coordinates
gcps['rel_y'] = gt[5] + _src[:, 1] * gt[3]
gcps['dj'] = gcps['search_j'] - gcps['match_j']
gcps['di'] = gcps['search_i'] - gcps['match_i']
gcps['elevation'] = 0
gcps.crs = mst.proj4

gcps.dropna(inplace=True)

# run ransac to find the matches between the transformed image and the master image make a coherent transformation
# residual_threshold is 10 pixels to allow for some local distortions, but get rid of the big blunders
Mfin, inliers_fin = ransac((gcps[['match_j', 'match_i']].values, gcps[['search_j', 'search_i']].values),
                           AffineTransform,
                           min_samples=10, residual_threshold=10, max_trials=1000)
gcps = gcps[inliers_fin]

print('loading dems')
dem = GeoImg(args.dem)
rel_dem = GeoImg(args.rel_dem)

gcps.crs = {'init': 'epsg:{}'.format(mst.epsg)}
for i, row in gcps.to_crs(crs=dem.proj4).iterrows():
    gcps.loc[i, 'elevation'] = dem.raster_points([(row.geometry.x, row.geometry.y)], nsize=9, mode='cubic')
    gcps.loc[i, 'el_rel'] = rel_dem.raster_points([(row.rel_x, row.rel_y)], nsize=9, mode='cubic')

# drop any gcps where we don't have a DEM value or a valid match
gcps.dropna(inplace=True)

print('{} valid matches found'.format(np.count_nonzero(inliers_fin)))

# ref_src = gcps[['j', 'i']].values
# ref_dst = gcps[['mj', 'mi']].values
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
subprocess.Popen(['get_autogcp_locations.sh', 'Ori-Relative',
                  os.path.join(out_dir, 'AutoMeasures{}.txt'.format(subscript))] + imlist).wait()

print('searching for points in non-orthorectified images')
E = builder.ElementMaker()
MesureSet = E.SetOfMesureAppuisFlottants()
for im in imlist:
    print(im)
    img = imread(im)
    maxi, maxj = img.shape

    impts = pd.read_csv('Auto-{}.txt'.format(im), sep=' ', names=['j', 'i'])
    impts_nodist = pd.read_csv('NoDist-{}.txt'.format(im), sep=' ', names=['j', 'i'])

    this_im_mes = E.MesureAppuiFlottant1Im(E.NameIm(im))
    for ind, row in impts.iterrows():
        in_im = 0 < row.j < maxj and 0 < row.i < maxi
        in_nd = -200 < impts_nodist.j[ind] < maxj + 200 and -200 < impts_nodist.i[ind] < maxi + 200
        if in_im and in_nd:
            try:
                testchip, _, _ = imtools.make_template(ortho, (gcps.orig_i[ind], gcps.orig_j[ind]), 20)
                dst_chip, _, _ = imtools.make_template(img, (row.i, row.j), 200)

                # test = imtools.highpass_filter(testchip)
                # dest = imtools.highpass_filter(dst_chip)
                test = np.ma.masked_values(testchip, 0)
                dest = np.ma.masked_values(dst_chip, 0)
                corr_res, this_i, this_j = imtools.find_gcp_match(dest.astype(np.float32), test.astype(np.float32))
                # z_corr = (corr_res.max() - corr_res.min()) / corr_res.std()
                peak_corr = cv2.minMaxLoc(corr_res)[1]

                # pks = peak_local_max(corr_res, min_distance=5, num_peaks=2)
                # this_z_corrs = []
                # for pk in pks:
                #     max_ = corr_res[pk[0], pk[1]]
                #     this_z_corrs.append((max_ - corr_res.mean()) / corr_res.std())
                # dz_corr = max(this_z_corrs) / min(this_z_corrs)
                # z_corr = max(this_z_corrs)

                # if args.corr_thresh < peak_corr or (z_corr > 5 and dz_corr > 1.5):
                if args.corr_thresh < peak_corr:
                    out_i, out_j = this_i - 200 + int(row.i), this_j - 200 + int(row.j)
                    this_mes = E.OneMesureAF1I(E.NamePt('GCP{}'.format(ind)),
                                               E.PtIm('{} {}'.format(out_j, out_i)))
                    this_im_mes.append(this_mes)
            except:
                continue
    MesureSet.append(this_im_mes)

tree = etree.ElementTree(MesureSet)
tree.write(os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript)),
           pretty_print=True, xml_declaration=True, encoding="utf-8")

print('running mm3d GCPBascule to estimate terrain errors')
subprocess.Popen(
    ['mm3d', 'GCPBascule', mmtools.get_match_pattern(imlist), 'Relative', 'TerrainRelAuto{}'.format(subscript),
     os.path.join(out_dir, 'AutoGCPs{}.xml'.format(subscript)),
     os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript))]).wait()

gcps = mmtools.get_bascule_residuals(os.path.join('Ori-TerrainRelAuto{}'.format(subscript), 'Result-GCP-Bascule.xml'),
                                     gcps)

# inds = gcps.residual < 5 * nmad(gcps.residual)
# gcps = gcps.loc[inds]

# gcps.to_file(os.path.join(out_dir, 'AutoGCPs{}.shp'.format(subscript)))
# write_auto_gcps(gcps, subscript, out_dir, utm_str)
# subprocess.Popen(['mm3d', 'GCPConvert', 'AppInFile',
#                   os.path.join(out_dir, 'AutoGCPs{}.txt'.format(subscript))]).wait()
#
# auto_root = ET.parse(os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript))).getroot()
# for im in auto_root.findall('MesureAppuiFlottant1Im'):
#     for pt in im.findall('OneMesureAF1I'):
#         if pt.find('NamePt').text not in gcps.id.values:
#             im.remove(pt)
# # # # save AutoMeasures
# out_xml = ET.ElementTree(auto_root)
# out_xml.write(os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript)),
#               encoding="utf-8", xml_declaration=True)
#
# subprocess.Popen(['mm3d', 'GCPBascule', get_match_pattern(imlist), 'Relative',
#                   'TerrainRelAuto{}'.format(subscript),
#                   os.path.join(out_dir, 'AutoGCPs{}.xml'.format(subscript)),
#                   os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript))]).wait()
#
# gcps = get_bascule_residuals(os.path.join('Ori-TerrainRelAuto{}'.format(subscript), 'Result-GCP-Bascule.xml'), gcps)
#
# out_inds = sliding_window_filter([ortho.shape[1], ortho.shape[0]], gcps, 1500, mindist=800)
# gcps = gcps.loc[out_inds]
# residuals = gcps.residual
# gcp_names = gcps.id
# #
# # out = np.abs(gcps.residual) > 5 * nmad(gcps.residual)
# # gcps = gcps[gcps.id.isin(gcp_names)]
# # out_list = gcp_names[out]
# #
# # gcps.loc[gcps.id.isin(out_list), 'z_corr'] = np.nan
# # gcps.dropna(inplace=True)
# #
# gcps.to_file(os.path.join(out_dir, 'AutoGCPs{}.shp'.format(subscript)))
# write_auto_gcps(gcps, subscript, out_dir, utm_str)
# subprocess.Popen(['mm3d', 'GCPConvert', 'AppInFile',
#                   os.path.join(out_dir, 'AutoGCPs{}.txt'.format(subscript))]).wait()
#
# auto_root = ET.parse(os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript))).getroot()
# for im in auto_root.findall('MesureAppuiFlottant1Im'):
#     for pt in im.findall('OneMesureAF1I'):
#         if pt.find('NamePt').text not in gcps.id.values:
#             im.remove(pt)
# # # save AutoMeasures
# out_xml = ET.ElementTree(auto_root)
# out_xml.write(os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript)),
#               encoding="utf-8", xml_declaration=True)
#
# subprocess.Popen(['mm3d', 'GCPBascule', get_match_pattern(imlist), 'Relative',
#                   'TerrainRelAuto{}'.format(subscript),
#                   os.path.join(out_dir, 'AutoGCPs{}.xml'.format(subscript)),
#                   os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript))]).wait()
#
# subprocess.Popen(['mm3d', 'Campari', get_match_pattern(imlist),
#                   'TerrainRelAuto{}'.format(subscript),
#                   'TerrainFinal{}'.format(subscript),
#                   'GCP=[{},5,{},2]'.format(os.path.join(out_dir, 'AutoGCPs{}.xml'.format(subscript)),
#                                            os.path.join(out_dir, 'AutoMeasures{}-S2D.xml'.format(subscript))),
#                   'SH=Homol', 'AllFree=1']).wait()
#
# gcps = get_bascule_residuals(os.path.join('Ori-TerrainRelAuto{}'.format(subscript), 'Result-GCP-Bascule.xml'), gcps)
# gcps = get_campari_residuals('Ori-TerrainFinal_block0/Residus.xml', gcps)
#
# # gcps = gcps[gcps.camp_res < 5*nmad(gcps.camp_res)]
#
# fig1 = plt.figure(figsize=(7, 5))
# plt.imshow(ortho[::5, ::5], cmap='gray', extent=[0, ortho.shape[1], ortho.shape[0], 0])
# plt.plot(gcps.j, gcps.i, 'r+')
# plt.quiver(gcps.j, gcps.i, gcps.xres, gcps.yres, color='r', scale=1000)
#
# plt.savefig(os.path.join(out_dir, 'relative_gcps{}.png'.format(subscript)), bbox_inches='tight', dpi=200)
# plt.close(fig1)
#
# fig2 = lsat_fullres.display(sfact=10, fig=plt.figure(figsize=(7,5)))
# plt.plot(gcps.geometry.x, gcps.geometry.y, 'r+')
# plt.savefig(os.path.join(out_dir, 'world_gcps{}.png'.format(subscript)), bbox_inches='tight', dpi=200)
# plt.close(fig2)
#
# print('cleaning up.')
# # remove Auto-im.tif.txt, NoDist-im.tif.txt, etc. Ori-Relative-NoDist
# shutil.rmtree('Ori-Relative-NoDist/')
#
# for txtfile in glob('Auto-OIS*.tif.txt') + \
#                glob('NoDist-OIS*.tif.txt'):
#     os.remove(txtfile)
print('end.')
