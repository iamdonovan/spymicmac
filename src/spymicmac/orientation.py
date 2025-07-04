"""
spymicmac.orientation is a collection of tools for working with image orientation using micmac.
"""
import os
from pathlib import Path
import subprocess
from collections import defaultdict
import numpy as np
import scipy
import pandas as pd
import geopandas as gpd
import geoutils as gu
import matplotlib.pyplot as plt
import lxml.etree as etree
import lxml.builder as builder
import xml.etree.ElementTree as ET
from glob import glob
from scipy.interpolate import LinearNDInterpolator
from shapely.geometry.point import Point
from shapely.geometry import LineString, MultiPoint
from skimage.transform import AffineTransform
from skimage.measure import ransac
from . import register, micmac
from typing import Union
from numpy.typing import NDArray
from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib


def plot_camera_centers(ori: str, ax: Union[Axes3D, None] = None) -> Axes3D:
    """
    Plot camera center locations in a 3D axis using matplotlib.

    :param ori: the name of the orientation directory (e.g., Ori-Relative).
    :param ax: an existing 3d matplotlib axis. If None, one will be created.
    :return: **ax** -- a matplotlib axis with the camera centers plotted
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

    ori_df = load_all_orientation(ori)

    ax.plot(ori_df.x, ori_df.y, ori_df.z, 's')

    return ax


def combine_block_measures(blocks: list, meas_out: str = 'AutoMeasures', gcp_out: str = 'AutoGCPs',
                           fn_mes: str = 'AutoMeasures_block', fn_gcp: str = 'AutoGCPs_block',
                           dirname: str ='auto_gcps', share_gcps: bool = False) -> None:
    """
    Combine GCPs and Measures files from multiple sub-blocks into a single file.

    :param blocks: a list of the sub-block numbers to combine
    :param meas_out: the output filename for the Measures file (no extension).
    :param gcp_out: the output filename for the GCP file (no extension).
    :param fn_mes: the name pattern of the measures files to combine.
    :param fn_gcp: the name pattern of the GCP files to combine.
    :param dirname: the output directory where the files are saved.
    :param share_gcps: GCPs are shared between blocks.
    """
    ngcp = 0 # keep track of the number of GCPs

    mes_dicts = list()
    gcp_dicts = list()
    gcp_shps = list()

    for b in blocks:
        # load dirname/AutoMeasures_block{b}-S2D.xml
        this_root = ET.parse(Path(dirname, fn_mes + f"{b}-S2D.xml")).getroot()
        this_gcp = gpd.read_file(Path(dirname, fn_gcp + f"{b}.shp"))

        if share_gcps:
            this_mes_dict = dict()
            this_gcp_dict = dict()

            for im in this_root.findall('MesureAppuiFlottant1Im'):
                this_name = im.find('NameIm').text
                these_mes = im.findall('OneMesureAF1I')

                this_mes_dict[this_name] = these_mes
        else:
            this_mes_dict, this_gcp_dict = micmac.rename_gcps(this_root, ngcp=ngcp)

            ngcp += len(this_gcp_dict)
            # load dirname/AutoGCPs_block{b}.shp

            for ii, row in this_gcp.iterrows():
                this_gcp.loc[ii, 'id'] = this_gcp_dict[row['id']]

        gcp_shps.append(this_gcp)
        mes_dicts.append(this_mes_dict)
        gcp_dicts.append(this_gcp_dict)

    out_gcp = gpd.GeoDataFrame(pd.concat(gcp_shps, ignore_index=True))

    if share_gcps:
        out_gcp = out_gcp[['id', 'elevation', 'geometry']].drop_duplicates(subset='id')

    out_gcp.sort_values('id', ignore_index=True, inplace=True)

    out_gcp.set_crs(gcp_shps[0].crs, inplace=True)
    out_gcp.to_file(Path(dirname, gcp_out + '.shp'))

    micmac.write_auto_gcps(out_gcp, '', dirname, register._get_utm_str(out_gcp.crs.to_epsg), outname=gcp_out)

    echo = subprocess.Popen('echo', stdout=subprocess.PIPE)
    p = subprocess.Popen(['mm3d', 'GCPConvert', 'AppInFile',
                          Path(dirname, gcp_out + '.txt')], stdin=echo.stdout)
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
    tree.write(Path(dirname, meas_out + '-S2D.xml'), pretty_print=True, xml_declaration=True, encoding="utf-8")


def block_orientation(blocks: list, meas_out: str = 'AutoMeasures', gcp_out: str = 'AutoGCPs',
                      fn_mes: str = 'AutoMeasures_block', fn_gcp: str = 'AutoGCPs_block', dirname: str = 'auto_gcps',
                      rel_ori: str = 'Relative', outori: str = 'TerrainFinal', homol: str = 'Homol',
                      ref_dx: Union[int, float] = 15, ortho_res: Union[int, float] = 8, allfree: bool = True,
                      max_iter: int = 1, share_gcps: bool = False) -> gpd.GeoDataFrame:
    """
    Combine GCPs, Measures files, and Ori directories from multiple sub-blocks into a single file and orientation.
    After combining blocks into a single orientation, runs mm3d Campari to refine the orientation and shared
    camera parameters.

    :param blocks: a list of the sub-block numbers to combine
    :param meas_out: the output filename for the Measures file (no extension).
    :param gcp_out: the output filename for the GCP file (no extension).
    :param fn_mes: the name pattern of the measures files to combine.
    :param fn_gcp: the name pattern of the GCP files to combine.
    :param dirname: the output directory where the files are saved.
    :param rel_ori: the name of the relative orientation to input to GCPBascule (default: Relative -> Ori-Relative)
    :param outori: the output orientation from Campari (default: TerrainFinal -> Ori-TerrainFinal).
    :param homol: the Homologue directory to use.
    :param ref_dx: the pixel resolution of the reference image, in meters.
    :param ortho_res: the pixel resolution of the orthoimage being used, in meters.
    :param allfree: run Campari with AllFree=1 (True), or AllFree=0 (False).
    :param max_iter: the maximum number of iterations to run.
    :param share_gcps: GCPs are shared between blocks
    :return: **gcps** -- the combined GCPs output from spymicmac.micmac.iterate_campari
    """
    combine_block_measures(blocks, meas_out=meas_out, gcp_out=gcp_out,
                           fn_mes=fn_mes, fn_gcp=fn_gcp, dirname=dirname, share_gcps=share_gcps)

    gcps = gpd.read_file(Path(dirname, gcp_out + '.shp'))

    gcps = micmac.iterate_campari(gcps, dirname, "OIS.*tif", '', ref_dx, ortho_res, fn_gcp=gcp_out,
                                  fn_meas=meas_out, rel_ori=rel_ori, outori=outori, homol=homol,
                                  allfree=allfree, max_iter=max_iter)
    return gcps


######################################################################################################################
# orientation tools - used for visualizing, manipulating camera orientation files
######################################################################################################################
def load_orientation(fn_img: str, ori: str) -> tuple[list, list, list, list, float, float]:
    """
    Read camera position and rotation information from an Orientation xml file.

    :param fn_img: the name of the image to read the orientation file for.
    :param ori: the name of the orientation directory (e.g., Ori-Relative).
    :return:
        - **centre** -- the camera position (x, y, z)
        - **l1** -- the L1 orientation parameters
        - **l2** -- the L2 orientation parameters
        - **l3** -- the L3 orientation parameters
        - **prof** -- the 'Profondeur' value from the xml file.
        - **altisol** -- the 'AltiSol' value from the xml file.

    """
    ori_root = ET.parse(Path(ori, f"Orientation-{fn_img}.xml")).getroot()
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


def load_all_orientation(ori: str, imlist: Union[list, None] = None) -> pd.DataFrame:
    """
    Load all of the orientation parameters for a set of images from a given directory.

    :param ori: the orientation directory to read
    :param imlist: the images to load. If not set, loads all orientation files from the given directory.
    :return: **ori_df** -- a DataFrame containing the orientation parameters for each image
    """
    df = pd.DataFrame()
    points = []

    if imlist is None:
        imlist = [os.path.basename(g).split('Orientation-')[1].split('.xml')[0] for g in
                  glob(Path(ori, 'Orientation*.xml'))]
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
    return gpd.GeoDataFrame(df)


def extend_line(df: pd.DataFrame, first: str, last: str) -> Point:
    """
    Extend a flightline using existing camera positions.

    :param df: a GeoDataFrame containing the camera positions and image names
    :param first: the name of the image to start interpolating from.
    :param last: the name of the image to end interpolating at.
    :return: **outpt** -- the new point along the flightline.
    """
    firstImg = df.loc[df.name.str.contains(first), 'geometry'].values[0]
    lastImg = df.loc[df.name.str.contains(last), 'geometry'].values[0]
    dx = firstImg.x - lastImg.x
    dy = firstImg.y - lastImg.y
    dz = firstImg.z - lastImg.z
    outpt = Point(firstImg.x + dx, firstImg.y + dy, firstImg.z + dz)

    return outpt


def interp_line(df: pd.DataFrame, first: str, last: str,
                nimgs: Union[int, None] = None, pos: Union[int, None] = None) -> Union[list, Point]:
    """
    Interpolate camera positions along a flightline.

    :param GeoDataFrame df: a GeoDataFrame containing the camera positions and image names
    :param first: the name of the image to start interpolating from.
    :param last: the name of the image to end interpolating at.
    :param nimgs: the number of images to interpolate (default: calculated based on the image numbers)
    :param pos: which image position to return (default: all images between first and last)
    :return: **ptList** -- a list containing the interpolated camera positions (or, a tuple of the requested
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


def update_center(fn_img: str, ori: str, new_center: list[float, float, float]) -> None:
    """
    Update the camera position in an Orientation file.

    :param fn_img: the name of the image to update the orientation for (e.g., 'OIS-Reech_ARCSEA000590122.tif')
    :param ori: the name of the orientation directory (e.g., 'Ori-Relative')
    :param new_center: a list of the new camera position [x, y, z]
    """
    ori_root = ET.parse(Path(ori, f"Orientation-{fn_img}.xml")).getroot()
    if ori_root.tag != 'OrientationConique':
        ori_coniq = ori_root.find('OrientationConique')
    else:
        ori_coniq = ori_root

    ori_coniq.find('Externe').find('Centre').text = ' '.join([str(f) for f in new_center])

    tree = ET.ElementTree(ori_root)
    tree.write(Path(ori, f"Orientation-{fn_img}.xml"),
               encoding="utf-8", xml_declaration=True)


def update_pose(fn_img: str, ori: str, new_rot: NDArray) -> None:
    """
    Update the camera pose (rotation matrix) in an Orientation file.

    :param fn_img: the name of the image to update the orientation for (e.g., 'OIS-Reech_ARCSEA000590122.tif')
    :param ori: the name of the orientation directory (e.g., 'Ori-Relative')
    :param new_rot: the new 3x3 rotation matrix
    """
    ori_root = ET.parse(Path(ori, f"Orientation-{fn_img}.xml")).getroot()
    if ori_root.tag != 'OrientationConique':
        ori_coniq = ori_root.find('OrientationConique')
    else:
        ori_coniq = ori_root

    for ii, row in enumerate(new_rot):
        ori_coniq.find('Externe').find('ParamRotation')\
            .find('CodageMatr').find(f"L{ii+1}").text = ' '.join([str(f) for f in row])

    tree = ET.ElementTree(ori_root)
    tree.write(Path(ori, f"Orientation-{fn_img}.xml"),
               encoding="utf-8", xml_declaration=True)


def update_params(fn_img: str, ori: str, profondeur: float, altisol: float) -> None:
    """
    Update the profondeur and altisol parameters in an orientation file.

    :param fn_img: the name of the image to update the orientation for (e.g., 'OIS-Reech_ARCSEA000590122.tif')
    :param ori: the name of the orientation directory (e.g., 'Ori-Relative')
    :param profondeur: the new profondeur value
    :param altisol: the new altisol value
    """
    ori_root = ET.parse(Path(ori, f"Orientation-{fn_img}.xml")).getroot()
    if ori_root.tag != 'OrientationConique':
        ori_coniq = ori_root.find('OrientationConique')
    else:
        ori_coniq = ori_root

    ori_coniq.find('Externe').find('AltiSol').text = str(altisol)
    ori_coniq.find('Externe').find('Profondeur').text = str(profondeur)

    tree = ET.ElementTree(ori_root)
    tree.write(Path(ori, f"Orientation-{fn_img}.xml"),
               encoding="utf-8", xml_declaration=True)


def write_orientation(ori_df: pd.DataFrame, dir_ori: str, calfile: Union[str, Path],
                      known_conv: str = 'eConvApero_DistM2C'):
    """
    Write orientation xml files for a set of images

    :param pd.DataFrame ori_df: a pandas DataFrame like the kind output by load_all_orientation()
    :param str dir_ori: the name of the output orientation directory (e.g., Ori-Relative)
    :param str calfile: the path to the calibration file for this orientation
    :param str known_conv: the name of a conversion used by MicMac
    """

    os.makedirs(dir_ori, exist_ok=True)

    for _, row in ori_df.iterrows():
        # TODO: allow for multiple cameras by finding the "right" calibration file
        write_ind_ori(fn_img=row['name'],
                      center=list(row[['x', 'y', 'z']]),
                      codage_mat=np.array([row[['l11', 'l12', 'l13']].values,
                                           row[['l21', 'l22', 'l23']].values,
                                           row[['l31', 'l32', 'l33']].values]),
                      profondeur=row['profondeur'],
                      altisol=row['altisol'],
                      dir_ori=dir_ori,
                      calfile=calfile,
                      known_conv=known_conv
                      )


def write_ind_ori(fn_img: str, center: list, codage_mat: NDArray, profondeur: float, altisol: float, dir_ori: str,
                  calfile: Union[str, Path], known_conv: str = 'eConvApero_DistM2C'):
    """
    Write an orientation xml file for an individual image.

    :param fn_img: the name of the image
    :param center: the camera center position (x, y, z)
    :param codage_mat: the rotation matrix for the image
    :param profondeur: the camera depth parameter
    :param altisol: the camera altisol parameter
    :param dir_ori: the name of the output orientation directory (e.g., Ori-Relative)
    :param calfile: the filename of the camera calibration file for this image
    :param known_conv: the name of a conversion used by MicMac
    """
    fn_out = Path(dir_ori, f"Orientation-{fn_img}.xml")

    E = builder.ElementMaker()

    OrientationConique = E.OrientationConique(
        E.OrIntImaM2C(E.I00('0 0'), E.V10('1 0'), E.V01('0 1')),
        E.TypeProj('eProjStenope'),
        E.ZoneUtileInPixel('true'),
        E.FileInterne('/'.join(['.', calfile])),

        E.Externe(
            E.AltiSol(altisol),
            E.Profondeur(profondeur),
            E.Time('-1e+30'),
            E.KnownConv(known_conv),
            E.Centre(' '.join([str(p) for p in center])),
            E.IncCentre('1 1 1'),

            E.ParamRotation(
                E.CodageMatr(
                    E.L1(' '.join([str(r) for r in codage_mat[0]])),
                    E.L2(' '.join([str(r) for r in codage_mat[1]])),
                    E.L3(' '.join([str(r) for r in codage_mat[2]])),
                )
            )
        ),
        E.ConvOri(E.KnownConv(known_conv))
    )

    tree = etree.ElementTree(OrientationConique)
    tree.write(fn_out, pretty_print=True, xml_declaration=True, encoding="utf-8")


def fix_orientation(cameras: pd.DataFrame, ori_df: pd.DataFrame, ori: str, nsig: Union[int, float] = 4) -> None:
    """
    Correct erroneous Tapas camera positions using an estimated affine transformation between the absolute camera
    locations and the relative locations read from the orientation directory.

    Once the positions have been updated, you should re-run Tapas using the InOri set to the directory; e.g., if you
    have updated Ori-Relative, you should run:

        mm3d Tapas RadialBasic "OIS.*tif" InOri=Relative Out=Relative LibFoc=0

    :param cameras: A DataFrame containing camera positions (x, y, z) and a 'name' column that contains
        the image names.
    :param ori_df: A DataFrame output from sPyMicMac.micmac.load_all_orientations, or that contains
        a 'name' column and camera positions in relative space (x, y, z)
    :param ori: the Orientation directory to update (e.g., Ori-Relative)
    :param nsig: the number of normalized absolute deviations from the median residual value to consider
        a camera an outlier
    """
    join = cameras.set_index('name').join(ori_df.set_index('name'), lsuffix='abs', rsuffix='rel')
    join.dropna(subset=['xabs', 'yabs', 'zabs', 'xrel', 'yrel', 'zrel'], inplace=True)

    model = AffineTransform()
    est = model.estimate(join[['xabs', 'yabs']].values, join[['xrel', 'yrel']].values)

    if not est:
        print('Unable to estimate transformation. Trying with RANSAC.')
        model, inliers = ransac((join[['xabs', 'yabs']].values, join[['xrel', 'yrel']].values), AffineTransform,
                                min_samples=10, residual_threshold=10, max_trials=10000)
        print(f"transformation found with {np.count_nonzero(inliers)} inliers")

    res = model.residuals(join[['xabs', 'yabs']].values, join[['xrel', 'yrel']].values)
    outliers = np.abs(res - np.nanmedian(res)) > nsig * register.nmad(res)

    if np.count_nonzero(outliers) > 0:
        interp = LinearNDInterpolator(join.loc[~outliers, ['xrel', 'yrel']].values, join.loc[~outliers, 'zrel'])
        print(f"found {np.count_nonzero(outliers)} outliers using nsig={nsig}")
        for name, row in join[outliers].iterrows():
            new_x, new_y = model(row[['xabs', 'yabs']].values)[0]
            new_z = interp(new_x, new_y)

            print(f"new location for {name}: {new_x}, {new_y}, {new_z}")
            print(f"writing new Orientation file for {name}")
            update_center(name, ori, [new_x, new_y, new_z])


def transform_centers(rel: gu.Raster, ref: gu.Raster, imlist: list, footprints: gpd.GeoDataFrame,
                      ori: str, imgeom: bool = True) -> tuple[AffineTransform, NDArray, gpd.GeoDataFrame]:
    """
    Use the camera centers in relative space provided by MicMac Orientation files, along with camera centers or
    footprints, to estimate a transformation between the relative coordinate system and the absolute coordinate system.

    :param rel: the relative image
    :param ref: the reference image to use to determine the output image shape
    :param imlist: a list of the images that were used for the relative orthophoto
    :param footprints: the (approximate) image footprints or camera centers. If geom_type is Polygon, the
        centroid will be used for the absolute camera positions.
    :param ori: name of orientation directory
    :param imgeom: calculate a transformation between image ij locations, rather than real-world coordinates
    :return:
        - **model** -- the estimated Affine Transformation between relative and absolute space
        - **inliers** -- a list of the inliers returned by skimage.measure.ransac
        - **join** -- the joined image footprints and relative orientation files
    """

    rel_ori = load_all_orientation(ori, imlist=imlist)

    footprints = footprints.to_crs(crs=ref.crs).copy()
    if all(footprints.geom_type == 'Polygon'):
        footprints['x'] = footprints.geometry.centroid.x
        footprints['y'] = footprints.geometry.centroid.y
    elif all(footprints.geom_type == 'Point'):
        footprints['x'] = footprints.geometry.x
        footprints['y'] = footprints.geometry.y
    else:
        raise ValueError("footprint geometry contains mixed types - please ensure that only Point or Polygon is used.")

    footprints['name'] = 'OIS-Reech_' + footprints['ID'] + '.tif'

    join = footprints.set_index('name').join(rel_ori.set_index('name'), lsuffix='abs', rsuffix='rel')
    join.dropna(subset='geometryrel', inplace=True)

    if join.shape[0] > 3:
        width_ratio = _point_spread(rel_ori.geometry)
        # if the points are very linear, we want to add a point to keep the transformation from being too sheared
        if width_ratio > 10:
            ind1, ind2 = _find_add([Point(row.xrel, row.yrel) for row in join.itertuples()])

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
        ref_pts = _get_points([Point(row.xabs, row.yabs) for row in join.itertuples()])
        rel_pts = _get_points([Point(row.xrel, row.yrel) for row in join.itertuples()])

    if imgeom:
        model, inliers = transform_points(ref, ref_pts, rel, rel_pts)
    else:
        model, inliers = _transform(ref_pts, rel_pts)

    print(f"{np.count_nonzero(inliers)} inliers for center transformation")
    return model, inliers, join


def transform_points(ref: gu.Raster, ref_pts: NDArray,
                     rel: gu.Raster, rel_pts: NDArray) -> tuple[AffineTransform, NDArray]:
    """
    Given x,y points and two "geo"-referenced images, finds an affine transformation between the two images.

    :param ref: the reference image
    :param ref_pts: an Mx2 array of the x,y points in the reference image
    :param rel: the second image
    :param rel_pts: an Mx2 array of the x,y points in the second image.
    :return:
        - **model** -- the estimated Affine Transformation between relative and absolute space
        - **inliers** -- an array of the inliers returned by skimage.measure.ransac
    """
    ref_ij = np.array(ref.xy2ij(ref_pts[:, 0], ref_pts[:, 1])).T

    rel_ij = np.array(rel.xy2ij(rel_pts[:, 0], rel_pts[:, 1])).T
    # rel_ij = np.array([((pt[0] - rel_gt[4]) / rel_gt[0],
    #                     (pt[1] - rel_gt[5]) / rel_gt[3]) for pt in rel_pts])

    model, inliers = _transform(ref_ij[:, ::-1], rel_ij[:, ::-1])

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
    # get two points perpendicular to the line, centered on the midpoint
    endpts = [Point(pt12.x - 0.5 * line.length * norm[0],
                    pt12.y - 0.5 * line.length * norm[1]),
              Point(pt12.x + 0.5 * line.length * norm[0],
                    pt12.y + 0.5 * line.length * norm[1])]


    pts = [(p.x, p.y) for p in [pt1, pt2, pt12] + endpts]

    return np.array(pts)


def _norm_vector(line):
    x, y = line.xy
    a = np.array([x[-1] - x[0], y[-1] - y[0]])
    b = np.array([a[1], -a[0]])
    b = b / np.linalg.norm(a)
    return b


def scale_intrinsics(fn_cam: Union[str, Path], scale: float) -> dict:
    """
    Parse an AutoCal xml file into a dictionary of intrinsic parameters, then scale those parameters - for example,
    to go from pixel values (default for MicMac) to mm.

    :param fn_cam: the name of the AutoCal xml file to parse.
    :param scale: the scale to use to convert the intrinsic values
    :return: **cam_dict** -- a dictionary of intrinsic parameter values for the camera
    """
    cam_dict = micmac.load_cam_xml(fn_cam)

    cam_dict['cdist'] = np.array([float(p) for p in cam_dict['cdist'].split(' ')]) * scale
    cam_dict['size'] = np.array([float(p) for p in cam_dict['size'].split(' ')]) * scale

    for ii in range(1, 4):
        cam_dict[f"K{ii}"] = float(cam_dict[f"K{ii}"]) * scale

    if 'K4' in cam_dict:
        cam_dict['K4'] = float(cam_dict['K4']) * scale

    if 'K5' in cam_dict:
        cam_dict['K5'] = float(cam_dict['K5']) * scale

    for pp in ['P1', 'P2']:
        if pp in cam_dict:
            cam_dict[pp] = float(cam_dict[pp]) * scale

    for pp in ['b1', 'b2']:
        if pp in cam_dict:
            cam_dict[pp] = float(cam_dict[pp]) * scale

    return cam_dict


def radial_distortion(params: dict, spacing: int) -> tuple[NDArray, NDArray, NDArray, NDArray]:
    """
    Apply a radial distortion model to a grid of pixel values, to produce a grid of un-distorted x,y locations.

    Uses values of {K1, K2, K3, ...} to apply the following formulas:

        Dx = Cx + (1 + K1 r**2 + K2 r**4 + K3 r**6 + ... Kn r**(2n)) * (x - Cx)
        Dy = Cy + (1 + K1 r**2 + K2 r**4 + K3 r**6 + ... Kn r**(2n)) * (y - Cy)

    Where Cx, Cy are the x,y center of distortion for the camera, Ki are the radial distortion coefficients, and
    x, y are the coordinates in the image space. Dx, Dy are the un-distorted locations that would be produced by an
    ideal (pinhole) camera.

    :param params: a dictionary of intrinsic parameter values for the camera
    :param spacing: the grid spacing to use for computing the (x, y) locations
    :returns:
        - **Dx**, **Dy** - the un-distorted x and y coordinates
        - **x**, **y** - the coordinates in the image space
    """
    sx, sy = params['size']
    cx, cy = params['cdist']

    xx, yy = np.meshgrid(np.arange(0, sx + 1, spacing),
                         np.arange(0, sy + 1, spacing))

    du = xx - cx
    dv = yy - cy

    rr = np.sqrt(du ** 2 + dv ** 2)

    order = len([kk for kk in params.keys() if 'K' in kk])
    radpart = 1
    for nn in range(1, order+1):
        radpart += params[f"K{nn}"] * rr ** (2 * nn)

    xdist = cx + radpart * du
    ydist = cy + radpart * dv

    return xdist, ydist, xx, yy


def standard_distortion(params: dict, spacing: int) -> tuple[NDArray, NDArray, NDArray, NDArray]:
    """
    Apply a standard distortion model to a grid of pixel values, to produce a grid of un-distorted x,y locations
    (see, for example, MicMac manual §15.2 "Distorsion specification", and eq. 15.13). The model consists of
    a radial distortion part Dr, a decentric part (P1, P2), and an affine part (b1, b2).

    Using the given parameters, the following formulas are applied:

        Dx = Dr,x + (2(x - Cx)**2 + r**2) * P1 + (2(x - Cx)(y - Cy)) * P2 + (x - Cx) * b1 + (y - Cy) * b2
        Dy = Dr,y + (2(x - Cx)(y - Cy)) * P1 + (2(y - Cy)**2 + r**2) * P2

    Where Dr are the x,y components of radial distortion, x, y are the coordinates in the image space, and r is the
    distance from the center of distortion. Dx, Dy are the un-distorted locations that would be produced by an ideal
    (pinhole) camera.

    :param params: a dictionary of intrinsic parameter values for the camera
    :param spacing: the grid spacing to use for computing the (x, y) locations
    :returns:
        - **Dx**, **Dy** - the un-distorted x and y coordinates
        - **x**, **y** - the coordinates in the image space
    """
    xdist, ydist, xx, yy = radial_distortion(params, spacing)

    cx, cy = params['cdist']

    du = xx - cx
    dv = yy - cy

    rr = np.sqrt(du**2 + dv**2)

    xdist += (2 * du**2 + rr**2) * params['P1'] + (2 * du * dv) * params['P2'] + params['b1'] * du + params['b2'] * dv
    ydist += (2 * du * dv ) * params['P1'] + (2 * dv**2 + rr**2) * params['P2']

    return xdist, ydist, xx, yy


def plot_lens_distortion(fn_cam: Union[str, Path],
                         spacing: int, scale: float = 1.0) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    """
    Plot a lens distortion curve (D% = d / r), where d is the offset caused by the lens distortion, and r is the
    radial distance from the center of distortion.

    :param fn_cam: the name of the AutoCal xml file to parse.
    :param spacing: the grid spacing to use for computing the (x, y) locations
    :param scale: the scale to use to convert the intrinsic values
    :returns: **fig**, **ax** - the Figure and Axes objects containing the plot.
    """
    cam_params = scale_intrinsics(fn_cam, scale)

    for pp in ['P1', 'P2', 'b1', 'b2']:
        if pp not in cam_params.keys():
            cam_params[pp] = 0.

    xdist, ydist, xx, yy = standard_distortion(cam_params, spacing)

    cx, cy = cam_params['cdist']
    dd = np.sqrt((xdist - xx)**2 + (ydist - yy)**2)
    rr = np.sqrt((xx - cx)**2 + (yy - cy)**2)

    fdist = scipy.interpolate.interp1d(rr.flatten(), 100 * dd.flatten() / rr.flatten())
    _rr = np.linspace(rr.min(), rr.max(), 1000)

    fig, ax = plt.subplots(1, 1, figsize=(6, 3))

    ax.plot(_rr, fdist(_rr), 'k')
    ax.set_ylabel('D (%)')
    ax.set_xlabel('radial distance')

    return fig, ax
