"""
spymicmac.micmac is a collection of tools for interfacing with MicMac
"""
import os
import subprocess
import shutil
import numpy as np
from osgeo import gdal
import pandas as pd
import lxml.etree as etree
import lxml.builder as builder
import difflib
import xml.etree.ElementTree as ET
from glob import glob
from scipy.interpolate import LinearNDInterpolator
from shapely.geometry.point import Point
from shapely.geometry import LineString
from shapely.strtree import STRtree
from skimage.io import imread
from skimage.measure import ransac
from skimage.transform import AffineTransform
from pybob.bob_tools import mkdir_p
from pybob.GeoImg import GeoImg
from pybob.ddem_tools import nmad
from spymicmac.usgs import get_usgs_footprints


######################################################################################################################
# MicMac interfaces - write xml files for MicMac to read
######################################################################################################################
def write_neighbour_images(imlist, fprints=None, nameField='ID', prefix='OIS-Reech_', fileExt='.tif',
                           dataset='AERIAL_COMBIN'):
    """
    Using a list of images and a collection of image footprints, return a list of potential image pairs for processing
    with Tapioca.

    :param list imlist: a list of (original) image names to use (e.g., without 'OIS-Reech\_')
    :param GeoDataFrame fprints: a vector dataset of footprint polygons. If not provided, will attempt to download
        metadata from USGS for the images.
    :param str nameField: the field in fprints table that contains the image name
    :param str prefix: the prefix attached to the image name read by Tapioca (default: 'OIS-Reech\_')
    :param str fileExt: the file extension for the images read by Tapioca (default: .tif)
    :param dataset: the USGS dataset name to search if no footprints are provided (default: AERIAL_COMBIN)
    """
    E = builder.ElementMaker()
    NamedRel = E.SauvegardeNamedRel()

    if fprints is None:
        fprints = get_usgs_footprints(imlist, dataset=dataset)
    else:
        fprints = fprints[fprints[nameField].isin(imlist)]

    s = STRtree([f for f in fprints['geometry'].values])

    for i, row in fprints.iterrows():
        fn = row[nameField]
        fp = row['geometry']

        print(fn)

        res = s.query(fp)
        intersects = [c for c in res if fp.intersection(c).area > 0]
        fnames = [fprints[nameField][fprints['geometry'] == c].values[0] for c in intersects]
        try:
            fnames.remove(fn)
        except ValueError:
            pass

        for f in fnames:
            this_pair = E.Cple(' '.join([prefix + fn + fileExt,
                                         prefix + f + fileExt]))
            NamedRel.append(this_pair)

    tree = etree.ElementTree(NamedRel)
    tree.write('FileImagesNeighbour.xml', pretty_print=True, xml_declaration=True, encoding="utf-8")


def get_gcp_meas(im_name, meas_name, in_dir, E, nodist=None, gcp_name='GCP'):
    """
    Create an lxml.builder.ElementMaker object with a GCP name and the image (row, pixel) location.

    :param str im_name: the image name to write the GCP location for.
    :param str meas_name: the name of the file to read the point locations from.
    :param str in_dir: the name of the directory where the images and measures files are located.
    :param lxml.builder.ElementMaker E: an ElementMaker object for writing to the xml file.
    :param str nodist: the name of the directory
    :param str gcp_name: the prefix (e.g., GCP0, GCP1, etc.) for the GCP name (default: GCP).
    :return:
        - **this_im_meas** (*lxml.builder.ElementMaker*) -- an ElementMaker object with the GCP location in the image.
    """
    im = gdal.Open(os.path.sep.join([in_dir, im_name]))
    maxj = im.RasterXSize
    maxi = im.RasterYSize

    impts = pd.read_csv(os.path.join(in_dir, meas_name), sep=' ', names=['j', 'i'])
    if nodist is not None:
        impts_nodist = pd.read_csv(os.path.join(in_dir, nodist), sep=' ', names=['j', 'i'])

    this_im_mes = E.MesureAppuiFlottant1Im(E.NameIm(im_name))
    for ind, row in impts.iterrows():
        in_im = 0 < row.j < maxj and 0 < row.i < maxi
        if nodist is not None:
            in_nd = -200 < impts_nodist.j[ind]+200 < maxj and -200 < impts_nodist.i[ind] < maxi+200
            in_im = in_im and in_nd
        if in_im:
            this_mes = E.OneMesureAF1I(
                                E.NamePt('{}{}'.format(gcp_name, ind+1)),
                                E.PtIm('{} {}'.format(row['j'], row['i']))
                            )
            this_im_mes.append(this_mes)
    return this_im_mes


def get_im_meas(gcps, E):
    """
    Populate an lxml.builder.ElementMaker object with GCP image locations, for writing to xml files.

    :param pandas.DataFrame gcps: a DataFrame with the GCPs to find image locations for.
    :param lxml.builder.ElementMaker E: an ElementMaker object for writing to the xml file.
    :return:
        - **pt_els** (*list*) -- a list of ElementMaker objects corresponding to each GCP image location.
    """
    pt_els = []
    for ind, row in gcps.iterrows():
        this_mes = E.OneMesureAF1I(
                        E.NamePt(row['gcp']),
                        E.PtIm('{} {}'.format(row['im_col'], row['im_row']))
                        )
        pt_els.append(this_mes)
    return pt_els


def parse_im_meas(fn_meas):
    """
    Read an xml file with GCP image locations into a pandas DataFrame.

    :param fn_meas: the name of the measures file to read.
    :return:
        - **gcp_df** (*pandas.DataFrame*) -- a DataFrame with gcp names and image locations.
    """
    gcp_df = pd.DataFrame()
    root = ET.parse(fn_meas).getroot()
    measures = root.findall('MesureAppuiFlottant1Im')[0]
    for i, mes in enumerate(measures.findall('OneMesureAF1I')):
        gcp_df.loc[i, 'name'] = mes.find('NamePt').text
        pt = mes.find('PtIm').text.split()
        gcp_df.loc[i, 'i'] = float(pt[1])
        gcp_df.loc[i, 'j'] = float(pt[0])
    return gcp_df


def generate_measures_files(joined=False):
    """
    Create id_fiducial.txt, MeasuresCamera.xml, and Tmp-SL-Glob.xml files for KH-9 Hexagon images.

    :param bool joined: generate files for joined scene (220x460 mm) instead of half (220x230mm)
    """
    i_list = np.arange(22, -1, -1)
    if not joined:
        j_list = np.arange(0, 24)
    else:
        j_list = np.arange(0, 47)

    J, I = np.meshgrid(np.arange(0, j_list.size), np.arange(0, i_list.size))
    gcp_names = list(zip(I[0, :], J[0, :]))
    for i in range(1, i_list.size):
        gcp_names.extend(list(zip(I[i, :], J[i, :])))

    JJ, II = np.meshgrid(np.round(j_list).astype(int), np.round(i_list).astype(int))
    ij = list(zip(II[0, :], JJ[0, :]))
    for i in np.arange(1, i_list.size):
        ij.extend(list(zip(II[i, :], JJ[i, :])))
    ij = 10 * np.array(ij)

    E = builder.ElementMaker()
    ImMes = E.MesureAppuiFlottant1Im(E.NameIm('Glob'))
    SpGlob = E.SetPointGlob()

    gcp_df = pd.DataFrame()
    with open('id_fiducial.txt', 'w') as f:
        for i, ind in enumerate(gcp_names):
            row, col = ind
            gcp_name = 'GCP_{}_{}'.format(row, col)
            gcp_df.loc[i, 'gcp'] = gcp_name
            gcp_df.loc[i, 'im_row'] = ij[i, 0]
            gcp_df.loc[i, 'im_col'] = ij[i, 1]

            pt_glob = E.PointGlob(E.Type('eNSM_Pts'),
                                  E.Name(gcp_name),
                                  E.LargeurFlou('0'),
                                  E.NumAuto('0'),
                                  E.SzRech('-1'))
            SpGlob.append(pt_glob)
            print(gcp_name, file=f)

    tree = etree.ElementTree(SpGlob)
    tree.write('Tmp-SL-Glob.xml', pretty_print=True, xml_declaration=True, encoding="utf-8")

    pt_els = get_im_meas(gcp_df, E)
    for p in pt_els:
        ImMes.append(p)

    outxml = E.SetOfMesureAppuisFlottants(ImMes)
    tree = etree.ElementTree(outxml)
    tree.write('MeasuresCamera.xml', pretty_print=True, xml_declaration=True, encoding="utf-8")

    with open('id_fiducial.txt', 'w') as f:
        for gcp in gcp_names:
            row, col = gcp
            print('GCP_{}_{}'.format(row, col), file=f)


def get_match_pattern(imlist):
    """
    Given a list of image names, return a match pattern that can be passed to MicMac command line functions.

    :param list imlist: a list of image names.
    :return:
        - **pattern** (*str*) -- a match pattern (e.g., "OIS.*tif") that can be passed to MicMac functions.
    """
    matches = []
    for i, this_im in enumerate(imlist[:-1]):
        for im in imlist[i + 1:]:
            matches.extend(list(difflib.SequenceMatcher(None, this_im, im).get_matching_blocks()))

    good_matches = set([m for m in matches if m.size > 0 and m.a == m.b])
    start_inds = set([m.a for m in good_matches])
    ind_lengths = [(ind, min([m.size for m in good_matches if m.a == ind])) for ind in start_inds]
    ind_lengths.sort()

    first, last = ind_lengths[0][1], ind_lengths[1][0]
    return imlist[0][:first] + '(' + '|'.join([im[first:last] for im in imlist]) + ')' + imlist[0][last:]


def write_auto_mesures(gcps, sub, outdir, outname='AutoMeasures'):
    """
    Write a file with GCP locations in relaive space (x, y, z) to use with get_autogcp_locations.sh

    :param pandas.DataFrame gcps: a DataFrame with the GCPs to save.
    :param str sub: the name of the block, if multiple blocks are being used (e.g., '_block1'). If not, use ''.
    :param str outdir: the output directory to save the files to.
    :param str outname: the base name of the file to create (default: AutoMeasures).
    """
    with open(os.path.join(outdir, '{}{}.txt'.format(outname, sub)), 'w') as f:
        for i, row in gcps.iterrows():
            print('{} {} {}'.format(row.rel_x, row.rel_y, row.el_rel), file=f)


def get_valid_image_points(shape, pts, pts_nodist):
    """
    Find which image points are located within an image based on the size of the image.

    :param shape: the shape of the image (rows, columns) to determine valid points for.
    :param pandas.DataFrame pts: a DataFrame containing point locations (i, j)
    :param pandas.DataFrame pts_nodist: a DataFrame containing point locations (i, j) calculated using no camera distortion.
    :return:
        - **valid_pts** (*array-like*) -- an array of the points that are located within the image shape.
    """
    maxi, maxj = shape

    in_im = np.logical_and.reduce((0 < pts.j, pts.j < maxj,
                                   0 < pts.i, pts.i < maxi))
    in_nd = np.logical_and.reduce((-200 < pts_nodist.j, pts_nodist.j < maxj + 200,
                                   -200 < pts_nodist.i, pts_nodist.i < maxi + 200))

    return np.logical_and(in_im, in_nd)


def write_image_mesures(imlist, gcps, outdir='.', sub='', ort_dir='Ortho-MEC-Relative'):
    """
    Create a Measures-S2D.xml file (row, pixel) for each GCP in each image from a list of image names.

    :param list imlist: a list of image names.
    :param pandas.DataFrame gcps: a DataFrame of GCPs.
    :param str outdir: the output directory to save the files to.
    :param str sub: the name of the block, if multiple blocks are being used (e.g., '_block1').
    :param str ort_dir: the Ortho-MEC directory where the images are located.
    """
    E = builder.ElementMaker()
    MesureSet = E.SetOfMesureAppuisFlottants()

    for im in imlist:
        print(im)
        img_geo = GeoImg(os.path.join(ort_dir, 'Ort_' + im))
        impts = pd.read_csv('Auto-{}.txt'.format(im), sep=' ', names=['j', 'i'])
        # impts_nodist = pd.read_csv('NoDist-{}.txt'.format(im), sep=' ', names=['j', 'i'])

        xmin, xmax, ymin, ymax = img_geo.find_valid_bbox()

        # valid_pts = get_valid_image_points(img.shape, impts, impts_nodist)
        valid_pts = np.logical_and.reduce([xmin <= gcps.rel_x, gcps.rel_x < xmax,
                                           ymin <= gcps.rel_y, gcps.rel_y < ymax])

        if np.count_nonzero(valid_pts) == 0:
            continue

        this_im_mes = E.MesureAppuiFlottant1Im(E.NameIm(im))

        for i, (ind, row) in enumerate(impts[valid_pts].iterrows()):
            this_mes = E.OneMesureAF1I(E.NamePt(gcps.iloc[ind]['id']),
                                       E.PtIm('{} {}'.format(row.j, row.i)))
            this_im_mes.append(this_mes)

        MesureSet.append(this_im_mes)

    tree = etree.ElementTree(MesureSet)
    tree.write(os.path.join(outdir, 'AutoMeasures{}-S2D.xml'.format(sub)),
               pretty_print=True, xml_declaration=True, encoding="utf-8")


def write_auto_gcps(gcp_df, sub, outdir, utm_zone, outname='AutoGCPs'):
    """
    Write GCP name, x, y, and z information to a text file to use with mm3d GCPConvert.

    :param pandas.DataFrame gcp_df: a DataFrame with the GCPs to save.
    :param str sub: the name of the block, if multiple blocks are being used (e.g., '_block1'). If not, use ''.
    :param str outdir: the output directory to save the files to.
    :param str utm_zone: the UTM zone name (e.g., 8N).
    :param str outname: the name to use for the GCPs file (default: AutoGCPs.txt)
    """
    with open(os.path.join(outdir, '{}{}.txt'.format(outname, sub)), 'w') as f:
        # print('#F= N X Y Z Ix Iy Iz', file=f)
        print('#F= N X Y Z', file=f)
        print('#Here the coordinates are in UTM {} X=Easting Y=Northing Z=Altitude'.format(utm_zone), file=f)
        for i, row in gcp_df.iterrows():
            # print('{} {} {} {} {} {} {}'.format(row.id, row.geometry.x, row.geometry.y, row.elevation,
            #                                        5/row.z_corr, 5/row.z_corr, 1), file=f)
            print('{} {} {} {}'.format(row.id, row.geometry.x, row.geometry.y, row.elevation), file=f)


def get_bascule_residuals(fn_basc, gcp_df):
    """
    Read a given GCPBascule residual file, and add the residuals to a DataFrame with GCP information.

    :param str fn_basc: the GCPBascule xml file to read the residuals from.
    :param pandas.DataFrame gcp_df: a DataFrame with the GCPs to read the residuals for.
    :return:
        - **gcp_df** (*pandas.DataFrame*) -- the input GCPs with the Bascule residuals added.
    """
    root = ET.parse(fn_basc).getroot()
    gcp_res = root.findall('Residus')
    gcp_names = np.array([res.find('Name').text for res in gcp_res])
    # residuals = np.array([float(res.find('Dist').text) for res in gcp_res])
    x_res = np.array([float(res.find('Offset').text.split()[0]) for res in gcp_res])
    y_res = np.array([float(res.find('Offset').text.split()[1]) for res in gcp_res])
    z_res = np.array([float(res.find('Offset').text.split()[2]) for res in gcp_res])
    dist = np.array([float(res.find('Dist').text) for res in gcp_res])

    for data_ in zip(gcp_names, x_res):
        gcp_df.loc[gcp_df.id == data_[0], 'xres'] = data_[1]

    for data_ in zip(gcp_names, y_res):
        gcp_df.loc[gcp_df.id == data_[0], 'yres'] = data_[1]

    for data_ in zip(gcp_names, z_res):
        gcp_df.loc[gcp_df.id == data_[0], 'zres'] = data_[1]

    for data_ in zip(gcp_names, dist):
        gcp_df.loc[gcp_df.id == data_[0], 'residual'] = data_[1]
    # gcp_df['residual'] = np.sqrt(gcp_df['xres'].values**2 + gcp_df['yres'].values**2)

    return gcp_df


def get_campari_residuals(fn_resids, gcp_df):
    """
    Read a given Campari residual file, and add the residuals to a DataFrame with GCP information.

    :param fn_resids: the Campari residual xml file to read.
    :param pandas.DataFrame gcp_df: a DataFrame with the GCPs to read the residuals for.
    :return:
        - **gcp_df** (*pandas.DataFrame*) -- the input GCPs with the Campari residuals added.
    """
    camp_root = ET.parse(fn_resids).getroot()

    last_iter = camp_root.findall('Iters')[-1].findall('OneAppui')
    camp_gcp_names = [a.find('Name').text for a in last_iter]

    err_max = []
    for a in last_iter:
        try:
            err_max.append(float(a.find('EcartImMax').text))
        except AttributeError:
            err_max.append(np.nan)

    camp_x = []
    camp_y = []
    camp_z = []

    for a in last_iter:
        try:
            camp_x.append(float(a.find('EcartFaiscTerrain').text.split()[0]))
            camp_y.append(float(a.find('EcartFaiscTerrain').text.split()[1]))
            camp_z.append(float(a.find('EcartFaiscTerrain').text.split()[2]))
        except AttributeError:
            camp_x.append(np.nan)
            camp_y.append(np.nan)
            camp_z.append(np.nan)

    for data_ in zip(camp_gcp_names, err_max):
        gcp_df.loc[gcp_df.id == data_[0], 'camp_res'] = data_[1]

    for data_ in zip(camp_gcp_names, camp_x):
        gcp_df.loc[gcp_df.id == data_[0], 'camp_xres'] = data_[1]

    for data_ in zip(camp_gcp_names, camp_y):
        gcp_df.loc[gcp_df.id == data_[0], 'camp_yres'] = data_[1]

    for data_ in zip(camp_gcp_names, camp_z):
        gcp_df.loc[gcp_df.id == data_[0], 'camp_zres'] = data_[1]

    return gcp_df


def move_bad_tapas(ori):
    """
    Read residual files output from Tapas (or Campari, GCPBascule), and move images with a NaN residual.

    :param str ori: the orientation directory to read the residuals file from (e.g., 'Ori-Relative').
    """
    root = ET.parse(os.path.join(ori, 'Residus.xml')).getroot()
    res_df = pd.DataFrame()

    nimgs = [len(a.findall('OneIm')) for a in root.findall('Iters')]

    nmax = max(nimgs)

    ind = np.where(np.array(nimgs) == nmax)[0].min()

    res_df['name'] = [a.find('Name').text for a in root.findall('Iters')[ind].findall('OneIm')]
    res_df['residual'] = [float(a.find('Residual').text) for a in root.findall('Iters')[ind].findall('OneIm')]
    res_df['pct_ok'] = [float(a.find('PercOk').text) for a in root.findall('Iters')[ind].findall('OneIm')]
    res_df['npts'] = [float(a.find('NbPts').text) for a in root.findall('Iters')[ind].findall('OneIm')]

    mkdir_p('bad')
    for im in res_df['name'][np.isnan(res_df.residual)]:
        print('{} -> bad/{}'.format(im, im))
        shutil.move(im, 'bad')


def run_bascule(in_gcps, outdir, img_pattern, sub, ori):
    """
    Interface for running mm3d GCPBascule and reading the residuals from the resulting xml file.

    :param pandas.DataFrame in_gcps: a DataFrame with the GCPs that are being input to Campari.
    :param str outdir: the output directory where the AutoGCPs.xml file is saved.
    :param str img_pattern: the match pattern for the images being input to Campari (e.g., "OIS.*tif")
    :param str sub: the name of the block, if multiple blocks are being used (e.g., '_block1'). If not, use ''.
    :param str ori: the name of the orientation directory (e.g., Ori-Relative).
    :return:
        - **out_gcps** (*pandas.DataFrame*) -- the input gcps with the updated Campari residuals.
    """
    echo = subprocess.Popen('echo', stdout=subprocess.PIPE)
    p = subprocess.Popen(['mm3d', 'GCPBascule', img_pattern, ori,
                          'TerrainRelAuto{}'.format(sub),
                          os.path.join(outdir, 'AutoGCPs{}.xml'.format(sub)),
                          os.path.join(outdir, 'AutoMeasures{}-S2D.xml'.format(sub))], stdin=echo.stdout)
    p.wait()

    out_gcps = get_bascule_residuals(os.path.join('Ori-TerrainRelAuto{}'.format(sub),
                                                  'Result-GCP-Bascule.xml'), in_gcps)
    return out_gcps


def run_campari(in_gcps, outdir, img_pattern, sub, dx, ortho_res, allfree=True):
    """
    Interface for running mm3d Campari and reading the residuals from the residual xml file.

    :param pandas.DataFrame in_gcps: a DataFrame with the GCPs that are being input to Campari.
    :param str outdir: the output directory where the AutoGCPs.xml file is saved.
    :param str img_pattern: the match pattern for the images being input to Campari (e.g., "OIS.*tif")
    :param str sub: the name of the block, if multiple blocks are being used (e.g., '_block1'). If not, use ''.
    :param int|float dx: the pixel resolution of the reference image.
    :param int|float ortho_res: the pixel resolution of the orthoimage being used.
    :param bool allfree: run Campari with AllFree=1 (True), or AllFree=0 (False). (default: True)
    :return:
        - **out_gcps** (*pandas.DataFrame*) -- the input gcps with the updated Campari residuals.
    """
    echo = subprocess.Popen('echo', stdout=subprocess.PIPE)
    if allfree:
        allfree_tag = 1
    else:
        allfree_tag = 0

    p = subprocess.Popen(['mm3d', 'Campari', img_pattern,
                          'TerrainRelAuto{}'.format(sub),
                          'TerrainFirstPass{}'.format(sub),
                          'GCP=[{},{},{},{}]'.format(os.path.join(outdir, 'AutoGCPs{}.xml'.format(sub)),
                                                     np.abs(dx),
                                                     os.path.join(outdir, 'AutoMeasures{}-S2D.xml'.format(sub)),
                                                     np.abs(dx / ortho_res)),
                          'SH=Homol', 'AllFree={}'.format(allfree_tag)], stdin=echo.stdout)
    p.wait()

    out_gcps = get_campari_residuals('Ori-TerrainFirstPass{}/Residus.xml'.format(sub), in_gcps)
    # out_gcps.dropna(inplace=True)  # sometimes, campari can return no information for a gcp
    return out_gcps


def save_gcps(in_gcps, outdir, utmstr, sub):
    """
    Save a GeoDataFrame of GCP information to shapefile, txt, and xml formats.

    After running, the following new files will be created:

        - outdir/AutoGCPs.shp (+ associated files)
        - outdir/AutoGCPs.txt
        - outdir/AutoGCPs.xml (output from mm3d GCPConvert)
        - outdir/AutoMeasures.xml (a file with image locations for each GCP)

    :param GeoDataFrame in_gcps: the gcps GeoDataFrame to save
    :param str outdir: the output directory to save the files to
    :param str utmstr: a UTM string generated by register.get_utm_str()
    :param str sub: the name of the block, if multiple blocks are being used (e.g., '_block1'). If not, use ''.

    """
    in_gcps.to_file(os.path.join(outdir, 'AutoGCPs{}.shp'.format(sub)))
    write_auto_gcps(in_gcps, sub, outdir, utmstr)
    echo = subprocess.Popen('echo', stdout=subprocess.PIPE)
    p = subprocess.Popen(['mm3d', 'GCPConvert', 'AppInFile',
                          os.path.join(outdir, 'AutoGCPs{}.txt'.format(sub))], stdin=echo.stdout)
    p.wait()

    auto_root = ET.parse(os.path.join(outdir, 'AutoMeasures{}-S2D.xml'.format(sub))).getroot()
    for im in auto_root.findall('MesureAppuiFlottant1Im'):
        for pt in im.findall('OneMesureAF1I'):
            if pt.find('NamePt').text not in in_gcps.id.values:
                im.remove(pt)

    # save AutoMeasures
    out_xml = ET.ElementTree(auto_root)
    out_xml.write(os.path.join(outdir, 'AutoMeasures{}-S2D.xml'.format(sub)),
                  encoding="utf-8", xml_declaration=True)


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
    :return:
        - **df** (*pandas.DataFrame*) -- a DataFrame containing the orientation parameters for each image
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
    :return:
        - **outpt** (*shapely.Point*) -- the new point along the flightline.
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
    :return:
        - **ptList** (*list*) -- a list containing the interpolated camera positions (or, a tuple of the requested position).
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


def fix_orientation(cameras, ori_df, ori, nsig=4):
    """
    Correct erroneous Tapas camera positions using an estimated affine transformation between the absolute camera locations
    and the relative locations read from the orientation directory.

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
    model.estimate(join[['xabs', 'yabs']].values, join[['xrel', 'yrel']].values)

    res = model.residuals(join[['xabs', 'yabs']].values, join[['xrel', 'yrel']].values)

    outliers = res - np.median(res) > nsig * nmad(res)
    if np.count_nonzero(outliers) > 0:
        interp = LinearNDInterpolator(join.loc[~outliers, ['xrel', 'yrel']].values, join.loc[~outliers, 'zrel'])
        print('found {} outliers using nsig={}'.format(np.count_nonzero(outliers), nsig))
        for name, row in join[outliers].iterrows():
            new_x, new_y = model(row[['xabs', 'yabs']].values)[0]
            new_z = interp(new_x, new_y)

            print('new location for {}: {}, {}, {}'.format(name, new_x, new_y, new_z))
            print('writing new Orientation file for {}'.format(name))
            update_center(name, ori, [new_x, new_y, new_z])


######################################################################################################################
#
######################################################################################################################
def dem_to_text(fn_dem, fn_out='dem_pts.txt', spacing=100):
    """
    Write elevations from a DEM raster to a text file for use in mm3d PostProc Banana.

    :param str fn_dem: the filename of the DEM to read.
    :param str fn_out: the name of the text file to write out (default: dem_pts.txt)
    :param int spacing: the pixel spacing of the DEM to write (default: every 100 pixels)
    """
    if isinstance(fn_dem, GeoImg):
        dem = fn_dem
    else:
        dem = GeoImg(fn_dem)
    x, y = dem.xy()

    z = dem.img[::spacing, ::spacing].flatten()
    x = x[::spacing, ::spacing].flatten()
    y = y[::spacing, ::spacing].flatten()

    x = x[np.isfinite(z)]
    y = y[np.isfinite(z)]
    z = z[np.isfinite(z)]

    with open(fn_out, 'w') as f:
        for pt in list(zip(x, y, z)):
            print(pt[0], pt[1], pt[2], file=f)

