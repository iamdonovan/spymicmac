"""
spymicmac.orientation is a collection of tools for working with image orientation using micmac.
"""
import os
import subprocess
from collections import defaultdict
import numpy as np
import pandas as pd
import geopandas as gpd
import lxml.etree as etree
import lxml.builder as builder
import xml.etree.ElementTree as ET
from glob import glob
from scipy.interpolate import LinearNDInterpolator
from shapely.geometry.point import Point
from shapely.geometry import LineString, MultiPoint
from skimage.transform import AffineTransform
from skimage.measure import ransac
from pybob.ddem_tools import nmad
from . import register, micmac


def combine_block_measures(blocks, meas_out='AutoMeasures', gcp_out='AutoGCPs',
                           fn_mes='AutoMeasures_block', fn_gcp='AutoGCPs_block', dirname='auto_gcps'):
    """
    Combine GCPs and Measures files from multiple sub-blocks into a single file.

    :param list blocks: a list of the sub-block numbers to combine
    :param str meas_out: the output filename for the Measures file (no extension). (default: AutoMeasures)
    :param str gcp_out: the output filename for the GCP file (no extension). (default: AutoGCPs)
    :param str fn_mes: the name pattern of the measures files to combine (default: AutoMeasures_block)
    :param str fn_gcp: the name pattern of the GCP files to combine (default: AutoGCPs_block)
    :param str dirname: the output directory where the files are saved (default: auto_gcps)
    """
    ngcp = 0 # keep track of the number of GCPs

    mes_dicts = list()
    gcp_dicts = list()
    gcp_shps = list()

    for b in blocks:
        # load dirname/AutoMeasures_block{b}-S2D.xml
        this_root = ET.parse(os.path.join(dirname, fn_mes + '{}-S2D.xml'.format(b))).getroot()

        this_mes_dict, this_gcp_dict = micmac.rename_gcps(this_root, ngcp=ngcp)
        mes_dicts.append(this_mes_dict)
        gcp_dicts.append(this_gcp_dict)

        ngcp += len(this_gcp_dict)
        # load dirname/AutoGCPs_block{b}.shp
        this_gcp = gpd.read_file(os.path.join(dirname, fn_gcp + '{}.shp'.format(b)))

        for ii, row in this_gcp.iterrows():
            this_gcp.loc[ii, 'id'] = this_gcp_dict[row['id']]

        gcp_shps.append(this_gcp)

    out_gcp = gpd.GeoDataFrame(pd.concat(gcp_shps, ignore_index=True))
    out_gcp.sort_values('id', ignore_index=True, inplace=True)

    out_gcp.set_crs(gcp_shps[0].crs, inplace=True)
    out_gcp.to_file(os.path.join(dirname, gcp_out + '.shp'))

    micmac.write_auto_gcps(out_gcp, '', dirname, register.get_utm_str(out_gcp.crs.to_epsg), outname=gcp_out)

    echo = subprocess.Popen('echo', stdout=subprocess.PIPE)
    p = subprocess.Popen(['mm3d', 'GCPConvert', 'AppInFile',
                          os.path.join(dirname, gcp_out + '.txt')], stdin=echo.stdout)
    p.wait()

    # now have to combine mes_dicts based on key names
    comb_dict = defaultdict(list)

    for d in mes_dicts:
        for im, mes in d.items():
            comb_dict[im].append(mes)

    # now have to create new AutoMeasures file from comb_dict
    E = builder.ElementMaker()
    MesureSet = E.SetOfMesureAppuisFlottants()

    # have to write combined measures to a single file
    for im, mes in comb_dict.items():
        this_im_mes = E.MesureAppuiFlottant1Im(E.NameIm(im))
        for m in mes[0]:
            this_mes = E.OneMesureAF1I(E.NamePt(m.find('NamePt').text), E.PtIm(m.find('PtIm').text))
            this_im_mes.append(this_mes)

        MesureSet.append(this_im_mes)

    tree = etree.ElementTree(MesureSet)
    tree.write(os.path.join(dirname, meas_out + '-S2D.xml'), pretty_print=True, xml_declaration=True, encoding="utf-8")


def block_orientation(blocks, meas_out='AutoMeasures', gcp_out='AutoGCPs',
                      fn_mes='AutoMeasures_block', fn_gcp='AutoGCPs_block', dirname='auto_gcps',
                      rel_ori='Relative', outori='TerrainFinal', homol='Homol',
                      ref_dx=15, ortho_res=8, allfree=True, max_iter=1):
    """
    Combine GCPs, Measures files, and Ori directories from multiple sub-blocks into a single file and orientation.

    :param list blocks: a list of the sub-block numbers to combine
    :param str meas_out: the output filename for the Measures file (no extension). (default: AutoMeasures)
    :param str gcp_out: the output filename for the GCP file (no extension). (default: AutoGCPs)
    :param str fn_mes: the name pattern of the measures files to combine (default: AutoMeasures_block)
    :param str fn_gcp: the name pattern of the GCP files to combine (default: AutoGCPs_block)
    :param str dirname: the output directory where the files are saved (default: auto_gcps)
    :param str rel_ori: the name of the relative orientation to input to GCPBascule (default: Relative -> Ori-Relative)
    :param str outori: the output orientation from Campari (default: TerrainFinal -> Ori-TerrainFinal)
    :param str homol: the Homologue directory to use (default: Homol)
    :param int|float ref_dx: the pixel resolution of the reference image, in meters. (default: 15)
    :param int|float ortho_res: the pixel resolution of the orthoimage being used, in meters. (default: 8)
    :param bool allfree: run Campari with AllFree=1 (True), or AllFree=0 (False). (default: True)
    :param int max_iter: the maximum number of iterations to run. (default: 1)
    :return: **gcps** (*GeoDataFrame*) -- the combined GCPs output from spymicmac.micmac.iterate_campari
    """
    combine_block_measures(blocks, meas_out=meas_out, gcp_out=gcp_out, fn_mes=fn_mes, fn_gcp=fn_gcp, dirname=dirname)

    gcps = gpd.read_file(os.path.join(dirname, gcp_out + '.shp'))

    gcps = micmac.iterate_campari(gcps, dirname, "OIS.*tif", '', ref_dx, ortho_res, fn_gcp=gcp_out,
                                  fn_meas=meas_out, rel_ori=rel_ori, outori=outori, homol=homol,
                                  allfree=allfree, max_iter=max_iter)
    return gcps


######################################################################################################################
# orientation tools - used for visualizing, manipulating camera orientation files
######################################################################################################################
def load_orientation(fn_img, ori):
    """
    Read camera position and rotation information from an Orientation xml file.

    :param str fn_img: the name of the image to read the orientation file for.
    :param str ori: the name of the orientation directory (e.g., Ori-Relative).
    :return:
        - **centre** (*list*) -- the camera position (x, y, z)
        - **l1** (*list*) -- the L1 orientation parameters
        - **l2** (*list*) -- the L2 orientation parameters
        - **l3** (*list*) -- the L3 orientation parameters
        - **prof** (*float*) -- the 'Profondeur' value from the xml file.
        - **altisol** (*float*) -- the 'AltiSol' value from the xml file.

    """
    ori_root = ET.parse(os.path.join(ori, 'Orientation-{}.xml'.format(fn_img))).getroot()
    if ori_root.tag != 'OrientationConique':
        ori_coniq = ori_root.find('OrientationConique')
    else:
        ori_coniq = ori_root
    centre = [float(p) for p in ori_coniq.find('Externe').find('Centre').text.split()]
    rotMat = ori_coniq.find('Externe').find('ParamRotation').find('CodageMatr')

    if rotMat is not None:
        l1 = [float(p) for p in rotMat.find('L1').text.split()]
        l2 = [float(p) for p in rotMat.find('L2').text.split()]
        l3 = [float(p) for p in rotMat.find('L3').text.split()]
    else:
        l1 = [np.nan, np.nan, np.nan]
        l2 = [np.nan, np.nan, np.nan]
        l3 = [np.nan, np.nan, np.nan]

    prof = ori_coniq.find('Externe').find('Profondeur')
    if prof is not None:
        prof = prof.text
    else:
        prof = np.nan

    altisol = ori_coniq.find('Externe').find('AltiSol')
    if altisol is not None:
        altisol = altisol.text
    else:
        altisol = np.nan

    return centre, l1, l2, l3, prof, altisol


def load_all_orientation(ori, imlist=None):
    """
    Load all of the orientation parameters for a set of images from a given directory.

    :param str ori: the orientation directory to read
    :param list imlist: the images to load. If not set, loads all orientation files from the given directory.
    :return: **ori_df** (*pandas.DataFrame*) -- a DataFrame containing the orientation parameters for each image
    """
    df = pd.DataFrame()
    points = []

    if imlist is None:
        imlist = [os.path.basename(g).split('Orientation-')[1].split('.xml')[0] for g in
                  glob(os.path.join(ori, '*.tif.xml'))]
        imlist.sort()

    for i, fn_img in enumerate(imlist):
        centre, l1, l2, l3, prof, altisol = load_orientation(fn_img, ori)

        df.loc[i, 'name'] = fn_img
        points.append(Point(centre[0], centre[1], centre[2]))
        df.loc[i, 'x'] = centre[0]
        df.loc[i, 'y'] = centre[1]
        df.loc[i, 'z'] = centre[2]

        df.loc[i, 'l11'] = l1[0]
        df.loc[i, 'l12'] = l1[1]
        df.loc[i, 'l13'] = l1[2]

        df.loc[i, 'l21'] = l2[0]
        df.loc[i, 'l22'] = l2[1]
        df.loc[i, 'l23'] = l2[2]

        df.loc[i, 'l31'] = l3[0]
        df.loc[i, 'l32'] = l3[1]
        df.loc[i, 'l33'] = l3[2]

        df.loc[i, 'profondeur'] = prof
        df.loc[i, 'altisol'] = altisol

    df['geometry'] = points
    return df


def extend_line(df, first, last):
    """
    Extend a flightline using existing camera positions.

    :param GeoDataFrame df: a GeoDataFrame containing the camera positions and image names
    :param str first: the name of the image to start interpolating from.
    :param str last: the name of the image to end interpolating at.
    :return: **outpt** (*shapely.Point*) -- the new point along the flightline.
    """
    firstImg = df.loc[df.name.str.contains(first), 'geometry'].values[0]
    lastImg = df.loc[df.name.str.contains(last), 'geometry'].values[0]
    dx = firstImg.x - lastImg.x
    dy = firstImg.y - lastImg.y
    dz = firstImg.z - lastImg.z
    outpt = Point(firstImg.x + dx, firstImg.y + dy, firstImg.z + dz)

    return outpt


def interp_line(df, first, last, nimgs=None, pos=None):
    """
    Interpolate camera positions along a flightline.

    :param GeoDataFrame df: a GeoDataFrame containing the camera positions and image names
    :param str first: the name of the image to start interpolating from.
    :param str last: the name of the image to end interpolating at.
    :param int nimgs: the number of images to interpolate (default: calculated based on the image numbers)
    :param int pos: which image position to return (default: all images between first and last)
    :return: **ptList** (*list*) -- a list containing the interpolated camera positions (or, a tuple of the requested
      position).
    """
    if nimgs is None:
        nimgs = np.abs(int(last) - int(first)).astype(int)
    firstImg = df.loc[df.name.str.contains(first), 'geometry'].values[0]
    lastImg = df.loc[df.name.str.contains(last), 'geometry'].values[0]
    flightLine = LineString([firstImg, lastImg])
    avgDist = flightLine.length / nimgs
    if pos is not None:
        return flightLine.interpolate(pos * avgDist)
    else:
        ptList = []
        for i in range(1, nimgs):
            ptList.append((i, flightLine.interpolate(i * avgDist)))
        return ptList


def update_center(fn_img, ori, new_center):
    """
    Update the camera position in an Orientation file.

    :param str fn_img: the name of the image to update the orientation for (e.g., 'OIS-Reech_ARCSEA000590122.tif')
    :param str ori: the name of the orientation directory (e.g., 'Ori-Relative')
    :param list new_center: a list of the new camera position [x, y, z]
    """
    ori_root = ET.parse(os.path.join(ori, 'Orientation-{}.xml'.format(fn_img))).getroot()
    if ori_root.tag != 'OrientationConique':
        ori_coniq = ori_root.find('OrientationConique')
    else:
        ori_coniq = ori_root

    ori_coniq.find('Externe').find('Centre').text = ' '.join([str(f) for f in new_center])

    tree = ET.ElementTree(ori_root)
    tree.write(os.path.join(ori, 'Orientation-{}.xml'.format(fn_img)),
               encoding="utf-8", xml_declaration=True)


def update_pose(fn_img, ori, new_rot):
    """
    Update the camera pose (rotation matrix) in an Orientation file.

    :param str fn_img: the name of the image to update the orientation for (e.g., 'OIS-Reech_ARCSEA000590122.tif')
    :param str ori: the name of the orientation directory (e.g., 'Ori-Relative')
    :param array-like new_rot: the new 3x3 rotation matrix
    """
    ori_root = ET.parse(os.path.join(ori, 'Orientation-{}.xml'.format(fn_img))).getroot()
    if ori_root.tag != 'OrientationConique':
        ori_coniq = ori_root.find('OrientationConique')
    else:
        ori_coniq = ori_root

    for ii, row in enumerate(new_rot):
        ori_coniq.find('Externe').find('ParamRotation')\
            .find('CodageMatr').find('L{}'.format(ii+1)).text = ' '.join([str(f) for f in row])

    tree = ET.ElementTree(ori_root)
    tree.write(os.path.join(ori, 'Orientation-{}.xml'.format(fn_img)),
               encoding="utf-8", xml_declaration=True)


def update_params(fn_img, ori, profondeur, altisol):
    """
    Update the profondeur and altisol parameters in an orientation file.

    :param str fn_img: the name of the image to update the orientation for (e.g., 'OIS-Reech_ARCSEA000590122.tif')
    :param str ori: the name of the orientation directory (e.g., 'Ori-Relative')
    :param float profondeur: the new profondeur value
    :param float altisol: the new altisol value
    """
    ori_root = ET.parse(os.path.join(ori, 'Orientation-{}.xml'.format(fn_img))).getroot()
    if ori_root.tag != 'OrientationConique':
        ori_coniq = ori_root.find('OrientationConique')
    else:
        ori_coniq = ori_root

    ori_coniq.find('Externe').find('AltiSol').text = str(altisol)
    ori_coniq.find('Externe').find('Profondeur').text = str(profondeur)

    tree = ET.ElementTree(ori_root)
    tree.write(os.path.join(ori, 'Orientation-{}.xml'.format(fn_img)),
               encoding="utf-8", xml_declaration=True)


def fix_orientation(cameras, ori_df, ori, nsig=4):
    """
    Correct erroneous Tapas camera positions using an estimated affine transformation between the absolute camera
    locations and the relative locations read from the orientation directory.

    Once the positions have been updated, you should re-run Tapas using the InOri set to the directory; e.g., if you
    have updated Ori-Relative, you should run:

        mm3d Tapas RadialBasic "OIS.*tif" InOri=Relative Out=Relative LibFoc=0

    :param pandas.DataFrame cameras: A DataFrame containing camera positions (x, y, z) and a 'name' column that contains
        the image names.
    :param pandas.DataFrame ori_df: A DataFrame output from sPyMicMac.micmac.load_all_orientations, or that contains
        a 'name' column and camera positions in relative space (x, y, z)
    :param str ori: the Orientation directory to update (e.g., Ori-Relative)
    :param int|float nsig: the number of normalized absolute deviations from the median residual value to consider
        a camera an outlier (default: 4)
    """
    join = cameras.set_index('name').join(ori_df.set_index('name'), lsuffix='abs', rsuffix='rel')

    model = AffineTransform()
    model.estimate(join.dropna()[['xabs', 'yabs']].values, join.dropna()[['xrel', 'yrel']].values)

    res = model.residuals(join[['xabs', 'yabs']].values, join[['xrel', 'yrel']].values)

    outliers = res - np.nanmedian(res) > nsig * nmad(res)
    if np.count_nonzero(outliers) > 0:
        interp = LinearNDInterpolator(join.loc[~outliers, ['xrel', 'yrel']].values, join.loc[~outliers, 'zrel'])
        print('found {} outliers using nsig={}'.format(np.count_nonzero(outliers), nsig))
        for name, row in join[outliers].iterrows():
            new_x, new_y = model(row[['xabs', 'yabs']].values)[0]
            new_z = interp(new_x, new_y)

            print('new location for {}: {}, {}, {}'.format(name, new_x, new_y, new_z))
            print('writing new Orientation file for {}'.format(name))
            update_center(name, ori, [new_x, new_y, new_z])


def transform_centers(img_gt, ref, imlist, footprints, ori, imgeom=True):
    """
    Use the camera centers in relative space provided by MicMac Orientation files, along with camera footprints,
    to estimate a transformation between the relative coordinate system and the absolute coordinate system.

    :param array-like img_gt: the image GeoTransform (as read from a TFW file)
    :param GeoImg ref: the reference image to use to determine the output image shape
    :param list imlist: a list of the images that were used for the relative orthophoto
    :param GeoDataFrame footprints: the (approximate) image footprints - the centroid will be used for the absolute
        camera positions.
    :param str ori: name of orientation directory
    :param bool imgeom: calculate a transformation between image ij locations (True)
    :return:
        - **model** (*AffineTransform*) -- the estimated Affine Transformation between relative and absolute space
        - **inliers** (*array-like*) -- a list of the inliers returned by skimage.measure.ransac
        - **join** (*GeoDataFrame*) -- the joined image footprints and relative orientation files
    """

    rel_ori = load_all_orientation(ori, imlist=imlist)

    footprints = footprints.to_crs(epsg=ref.epsg).copy()

    for i, row in footprints.iterrows():
        footprints.loc[i, 'x'] = row['geometry'].centroid.x
        footprints.loc[i, 'y'] = row['geometry'].centroid.y
        footprints.loc[i, 'name'] = 'OIS-Reech_' + row['ID'] + '.tif'

    join = footprints.set_index('name').join(rel_ori.set_index('name'), lsuffix='abs', rsuffix='rel')
    join.dropna(inplace=True)

    if join.shape[0] > 3:
        width_ratio = _point_spread(rel_ori.geometry)
        # if the points are very linear, we want to add a point to keep the transformation from being too sheared
        if width_ratio > 10:
            ind1, ind2 = _find_add([Point(row.xrel, row.yrel) for ii, row in join.iterrows()])

            ref_pts = np.concatenate([join[['xabs', 'yabs']].values,
                                      _get_points([Point(join['xabs'].values[ind1], join['yabs'].values[ind1]),
                                                   Point(join['xabs'].values[ind2], join['yabs'].values[ind2])])])
            rel_pts = np.concatenate([join[['xrel', 'yrel']].values,
                                      _get_points([Point(join['xrel'].values[ind1], join['yrel'].values[ind1]),
                                                   Point(join['xrel'].values[ind2], join['yrel'].values[ind2])])])
        else:
            ref_pts = join[['xabs', 'yabs']].values
            rel_pts = join[['xrel', 'yrel']].values

    else:
        # if we only have 2 points, we add two (midpoint, perpendicular to midpoint) using _get_points()
        # this ensures that we can actually get an affine transformation
        ref_pts = _get_points([Point(row.xabs, row.yabs) for ii, row in join.iterrows()])
        rel_pts = _get_points([Point(row.xrel, row.yrel) for ii, row in join.iterrows()])

    if imgeom:
        model, inliers = transform_points(ref, ref_pts, img_gt, rel_pts)
    else:
        model, inliers = _transform(ref_pts, rel_pts)

    print('{} inliers for center transformation'.format(np.count_nonzero(inliers)))
    return model, inliers, join


def transform_points(ref, ref_pts, rel_gt, rel_pts):
    """
    Given x,y points and two "geo"-referenced images, finds an affine transformation between the two images.

    :param GeoImg ref: the reference image
    :param np.array ref_pts: an Mx2 array of the x,y points in the reference image
    :param array-like rel_gt: the "geo" transform for the second image, as read from a .tfw file.
    :param np.array rel_pts: an Mx2 array of the x,y points in the second image.
    :return:
        - **model** (*AffineTransform*) -- the estimated Affine Transformation between relative and absolute space
        - **inliers** (*array-like*) -- a list of the inliers returned by skimage.measure.ransac
    """
    ref_ij = np.array([ref.xy2ij(pt) for pt in ref_pts])
    rel_ij = np.array([((pt[0] - rel_gt[4]) / rel_gt[0],
                        (pt[1] - rel_gt[5]) / rel_gt[3]) for pt in rel_pts])

    model, inliers = _transform(ref_ij[:, ::-1], rel_ij)

    return model, inliers


def _transform(ref_pts, rel_pts):
    if ref_pts.shape[0] > 3:
        model, inliers = ransac((ref_pts, rel_pts), AffineTransform, min_samples=3,
                                residual_threshold=100, max_trials=5000)
    else:
        model = AffineTransform()
        model.estimate(ref_pts, rel_pts)
        residuals = model.residuals(ref_pts, rel_pts)
        inliers = residuals <= 1

    return model, inliers


def _find_add(pts):
    # find two points that are (a) far from the center of a distribution, and (b) far from each other
    cent = MultiPoint(pts).centroid
    cdist = [cent.distance(pt) for pt in pts]

    # find the point furthest from the centroid
    pt1 = pts[np.argmax(cdist)]

    # find the point furthest from that point
    pdist = [pt1.distance(pt) for pt in pts]

    return np.argmax(cdist), np.argmax(pdist)


def _point_spread(pts):
    # get the ratio of the length to the width of the minimum rotated rectangle covering a set of points
    rect = MultiPoint(pts).minimum_rotated_rectangle
    verts = [Point(pt) for pt in list(zip(rect.boundary.xy[0], rect.boundary.xy[1]))]
    dists = [verts[0].distance(pt) for pt in verts[1:]]

    dists.remove(max(dists))  # remove the longest - this is a diagonal
    dists.remove(min(dists))  # remove the shortest - this is the same point

    return max(dists) / min(dists)


def _get_points(centers):
    pt1 = centers[0]  # the first point
    pt2 = centers[1]  # the second point

    line = LineString([pt1, pt2])  # form a line between point 1, point 2
    norm = _norm_vector(line)  # get the normal vector to the line

    pt12 = line.centroid  # get the midpoint of the line
    # get a point perpendicular to the line at a distance of line.length from the midpoint
    endpt = Point(pt12.x + line.length * norm[0], pt12.y + line.length * norm[1])

    pts = [(p.x, p.y) for p in [pt1, pt2, pt12, endpt]]

    return np.array(pts)


def _norm_vector(line):
    x, y = line.xy
    a = np.array([x[-1] - x[0], y[-1] - y[0]])
    b = np.array([a[1], -a[0]])
    b = b / np.linalg.norm(a)
    return b
