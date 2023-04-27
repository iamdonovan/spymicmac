"""
spymicmac.micmac is a collection of tools for interfacing with MicMac
"""
import os
import re
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
from shapely.strtree import STRtree
from skimage.io import imread, imsave
from pybob.GeoImg import GeoImg
from pybob.ddem_tools import nmad
from pybob.image_tools import create_mask_from_shapefile
from spymicmac import data, register


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
        fprints = data.get_usgs_footprints(imlist, dataset=dataset)
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


def write_xml(fn_img, fn_mask='./MEC-Malt/Masq_STD-MALT_DeZoom1.tif', fn_xml=None, geomname='eGeomMNTEuclid'):
    """
    Given a GDAL dataset, create a MicMac xml worldfile.

    :param str fn_img: the filename of the image.
    :param str fn_mask: the filename of the mask file (default: ./MEC-Malt/Masq_STD-MALT_DeZoom1.tif)
    :param str fn_xml: the filename of the xml file to create (default: fn_img + '.xml')
    :param str geomname: the MicMac Geometry name to use (default: eGeomMNTEuclid)
    """
    ds = gdal.Open(fn_img)
    ext = os.path.splitext(fn_img)[-1]
    ulx, dx, _, uly, _, dy = ds.GetGeoTransform()

    E = builder.ElementMaker()
    FileOriMnt = E.FileOriMnt
    NameFileMnt = E.NameFileMnt
    NameFileMasque = E.NameFileMasque
    NombrePixels = E.NombrePixels
    OriginePlani = E.OriginePlani
    ResolutionPlani = E.ResolutionPlani
    OrigineAlti = E.OrigineAlti
    ResolutionAlti = E.ResolutionAlti
    Geometrie = E.Geometrie

    outxml = FileOriMnt(
        NameFileMnt(fn_img),
        NameFileMasque(fn_mask),
        NombrePixels(' '.join([str(ds.RasterXSize), str(ds.RasterYSize)])),
        OriginePlani(' '.join([str(ulx), str(uly)])),
        ResolutionPlani(' '.join([str(dx), str(dy)])),
        OrigineAlti('0'),
        ResolutionAlti('1'),
        Geometrie(geomname)
    )

    tree = etree.ElementTree(outxml)
    if fn_xml is None:
        fn_xml = fn_img.replace(ext, '.xml')

    tree.write(fn_xml, pretty_print=True,
               xml_declaration=False, encoding="utf-8")


def get_gcp_meas(im_name, meas_name, in_dir, E, nodist=None, gcp_name='GCP'):
    """
    Create an lxml.builder.ElementMaker object with a GCP name and the image (row, pixel) location.

    :param str im_name: the image name to write the GCP location for.
    :param str meas_name: the name of the file to read the point locations from.
    :param str in_dir: the name of the directory where the images and measures files are located.
    :param lxml.builder.ElementMaker E: an ElementMaker object for writing to the xml file.
    :param str nodist: the name of the directory
    :param str gcp_name: the prefix (e.g., GCP0, GCP1, etc.) for the GCP name (default: GCP).
    :return: **this_im_meas** (*lxml.builder.ElementMaker*) -- an ElementMaker object with the GCP location in
      the image.
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
    :return: **pt_els** (*list*) -- a list of ElementMaker objects corresponding to each GCP image location.
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
    :return: **gcp_df** (*pandas.DataFrame*) -- a DataFrame with gcp names and image locations.
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


def create_localchantier_xml(name='KH9MC', short_name='KH-9 Hexagon Mapping Camera', film_size=(460, 220),
                             pattern='.*', focal=304.8, add_sfs=False):
    """
    Create a MicMac-LocalChantierDescripteur.xml file for a given camera. Default is the KH-9 Hexagon Mapping Camera.

    :param str name: The name to use for the camera [KH9MC]
    :param str short_name: A short description of the camera [KH-9 Hexagon Mapping Camera]
    :param array-like film_size: the film size (width, height) in mm [460, 220]
    :param str pattern: the matching pattern to use for the images [.*]
    :param float focal: the nominal focal length, in mm [304.8]
    :param bool add_sfs: use SFS to help find tie points in low-contrast images [False]
    """
    E = builder.ElementMaker()

    chantier = E.ChantierDescripteur(
        E.LocCamDataBase(
            E.CameraEntry(
                E.Name(name),
                E.SzCaptMm('{} {}'.format(film_size[0], film_size[1])),
                E.ShortName(short_name)
            )
        ),
        E.KeyedNamesAssociations(
            E.Calcs(
                E.Arrite('1 1'),
                E.Direct(
                    E.PatternTransform(pattern),
                    E.CalcName(name)
                )
            ),
            E.Key('NKS-Assoc-STD-CAM')
        ),
        E.KeyedNamesAssociations(
            E.Calcs(
                E.Arrite('1 1'),
                E.Direct(
                    E.PatternTransform('OIS.*'),
                    E.CalcName('{}'.format(focal))
                )
            ),
            E.Key('NKS-Assoc-STD-FOC')
        )
    )

    if add_sfs:
        chantier.append(
            E.KeyedNamesAssociations(
                E.Calcs(
                    E.Arrite('1 1'),
                    E.Direct(
                        E.PatternTransform(pattern),
                        E.CalcName('SFS'),
                    )
                ),
                E.Key('NKS-Assoc-SFS')
            )
        )

    outxml = E.Global(chantier)
    tree = etree.ElementTree(outxml)
    tree.write('MicMac-LocalChantierDescripteur.xml', pretty_print=True, xml_declaration=True, encoding="utf-8")


def init_autocal(imsize=(32200, 15400), framesize=(460, 220), foc=304.8, camname='KH9MC'):
    """
    Create an AutoCal xml file for use in the Tapas step. Default values are for KH-9 Hexagon Mapping Camera.

    When calling mm3d Tapas, be sure to use "InCal=Init":

        mm3d Tapas RadialBasic "OIS.*tif" InCal=Init Out=Relative LibFoc=0

    The name of the file changes based on the focal length and camera name. Using the default values of
    foc=304.8 and camname='KH9MC' creates the following file in Ori-Init:

        AutoCal_Foc-KH9MC_304800.xml

    :param array-like imsize: the size of the image (width, height) in pixels (default: (32200, 15400))
    :param array-like framesize: the size of the image (width, height) in mm (default: (460, 220))
    :param float foc: nominal focal length, in mm (default: 304.8)
    :param int camname: the camera short name to use (default: KH9MC)
    """
    os.makedirs('Ori-Init', exist_ok=True)

    scale = np.mean([imsize[0] / framesize[0], imsize[1] / framesize[1]])

    pp = np.array(imsize) / 2

    E = builder.ElementMaker()
    outxml = E.ExportAPERO(
        E.CalibrationInternConique(
            E.KnownConv('eConvApero_DistM2C'),
                E.PP('{} {}'.format(pp[0], pp[1])),
                E.F('{}'.format(scale * foc)),
                E.SzIm('{} {}'.format(imsize[0], imsize[1])),
                E.CalibDistortion(
                    E.ModRad(
                        E.CDist('{} {}'.format(pp[0], pp[1])),
                        E.CoeffDist('2.74e-11'),
                        E.CoeffDist('-1.13e-21'),
                        E.CoeffDist('4.01e-29'),
                        E.CoeffDist('1.28e-38'),
                        E.CoeffDist('-4.32e-46'),
                        E.CeoffDistInv('2.74e-11'),
                        E.CeoffDistInv('3.88e-21'),
                        E.CeoffDistInv('-4.79e-29'),
                        E.CeoffDistInv('3.82e-38'),
                        E.CeoffDistInv('2.13e-46'),
                        E.CeoffDistInv('4.47e-55'),
                        E.PPaEqPPs('true')
                    )
                )
        )
    )

    tree = etree.ElementTree(outxml)
    tree.write(os.path.join('Ori-Init', 'AutoCal_Foc-{}_{}.xml'.format(int(foc*1000), camname)),
               pretty_print=True, xml_declaration=True, encoding="utf-8")


def parse_localchantier(fn_chant):
    pass


def get_match_pattern(imlist):
    """
    Given a list of image names, return a match pattern that can be passed to MicMac command line functions.

    :param list imlist: a list of image names.
    :return: **pattern** (*str*) -- a match pattern (e.g., "OIS.*tif") that can be passed to MicMac functions.
    """
    imlist.sort()

    first_ind = len(imlist[0])
    last_ind = len(imlist[0])

    for fn_img in imlist[1:]:
        diff = list(difflib.ndiff(imlist[0], fn_img))
        first_ind = min(first_ind, diff.index(next(d for d in diff if len(d.strip()) > 1)))
        last_ind = min(last_ind, diff[::-1].index(next(d for d in diff if len(d.strip()) > 1)))

    last_ind = len(diff) - last_ind # because we reversed diff to find the last ind

    first = imlist[0][:first_ind]
    last = imlist[0][last_ind:]

    middle = ''.join(['(', '|'.join([fn_img[first_ind:last_ind] for fn_img in imlist]), ')'])

    return first + middle + last


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
    :return: **valid_pts** (*array-like*) -- an array of the points that are located within the image shape.
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


def remove_measure(fn_meas, name):
    """
    Remove all instances of a given measure from an xml file.

    :param str fn_meas: the xml file (e.g., AutoMeasures-S2D.xml)
    :param str name: the measurement name (e.g., GCP0)
    """
    root = ET.parse(fn_meas).getroot()
    for im in root.findall('MesureAppuiFlottant1Im'):
        for mes in im.findall('OneMesureAF1I'):
            if mes.find('NamePt').text == name:
                im.remove(mes)

    tree = ET.ElementTree(root)
    tree.write(fn_meas, encoding="utf-8", xml_declaration=True)


def rename_gcps(root, ngcp=0):
    """
    Rename all GCPs in order of their appearance in an (opened) xml file.

    :param xml.etree.ElementTree.Element root: the root element of an xml tree
    :param int ngcp: the number to start counting from (defaults to 0)

    :return:
        - **mes_dict** (*dict*) -- a dict containing image, gcp key/value pairs
        - **gcp_dict** (*dict*) -- a dict containing old/new gcp name key/value pairs
    """
    mes_dict = dict()
    gcp_dict = dict()

    for im in root.findall('MesureAppuiFlottant1Im'):
        this_name = im.find('NameIm').text
        these_mes = im.findall('OneMesureAF1I')

        for ii, mes in enumerate(these_mes):
            old_name = mes.find('NamePt').text

            if mes.find('NamePt').text not in gcp_dict.keys():
                gcp_dict[old_name] = 'GCP{}'.format(ngcp)
                ngcp += 1

            mes.find('NamePt').text = gcp_dict[old_name]
            these_mes[ii] = mes

        mes_dict[this_name] = these_mes

    return mes_dict, gcp_dict


def get_bascule_residuals(fn_basc, gcp_df):
    """
    Read a given GCPBascule residual file, and add the residuals to a DataFrame with GCP information.

    :param str fn_basc: the GCPBascule xml file to read the residuals from.
    :param pandas.DataFrame gcp_df: a DataFrame with the GCPs to read the residuals for.
    :return: **gcp_df** (*pandas.DataFrame*) -- the input GCPs with the Bascule residuals added.
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
    :return: **gcp_df** (*pandas.DataFrame*) -- the input GCPs with the Campari residuals added.
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

    os.makedirs('bad', exist_ok=True)
    for im in res_df['name'][np.isnan(res_df.residual)]:
        print('{} -> bad/{}'.format(im, im))
        shutil.move(im, 'bad')


def tapioca(img_pattern='OIS.*tif', res_low=400, res_high=1200):
    """
    Run mm3d Tapioca MulScale

    :param str img_pattern: The image pattern to pass to Tapioca (default: OIS.*tif)
    :param int res_low: the size of the largest image axis, in pixels, for low-resolution matching (default: 400)
    :param int res_high: the size of the largest image axis, in pixels, for high-resolution matching (default: 1200)
    """
    echo = subprocess.Popen('echo', stdout=subprocess.PIPE)
    p = subprocess.Popen(['mm3d', 'Tapioca', 'MulScale', img_pattern,
                          str(res_low), str(res_high)], stdin=echo.stdout)
    return p.wait()


def tapas(cam_model, ori_out, img_pattern='OIS.*tif', in_cal=None, lib_foc=True, lib_pp=True, lib_cd=True):
    """
    Run mm3d Tapas with a given camera calibration model.

    Some basic camera calibration models for air photos:
        - RadialBasic
        - RadialStd
        - RadialExtended
        - FraserBasic
        - Fraser

    See MicMac docs for a full list/explanation of the camera models.

    :param str cam_model: the camera calibration model to use.
    :param str ori_out: the output orientation. Will create a directory, Ori-{ori_out}, with camera parameter files.
    :param str img_pattern: the image pattern to pass to Tapas (default: OIS.*tif)
    :param str in_cal: an input calibration model to refine (default: None)
    :param bool lib_foc: allow the focal length to be calibrated (default: True)
    :param bool lib_pp: allow the principal point to be calibrated (default: True)
    :param bool lib_cd: allow the center of distortion to be calibrated (default: True)
    """
    echo = subprocess.Popen('echo', stdout=subprocess.PIPE)

    if in_cal is not None:
        p = subprocess.Popen(['mm3d', 'Tapas', cam_model, img_pattern, 'InCal=' + in_cal,
                              'LibFoc={}'.format(int(lib_foc)), 'LibPP={}'.format(int(lib_pp)),
                              'LibCD={}'.format(int(lib_cd)), 'Out=' + ori_out], stdin=echo.stdout)
    else:
        p = subprocess.Popen(['mm3d', 'Tapas', cam_model, img_pattern,
                              'LibFoc={}'.format(int(lib_foc)), 'LibPP={}'.format(int(lib_pp)),
                              'LibCD={}'.format(int(lib_cd)), 'Out=' + ori_out], stdin=echo.stdout)
    return p.wait()


def apericloud(ori, img_pattern='OIS.*tif'):
    """
    Run mm3d AperiCloud to create a point cloud layer

    :param str ori: the input orientation to use
    :param str img_pattern: the image pattern to pass to AperiCloud (default: OIS.*tif)
    """
    echo = subprocess.Popen('echo', stdout=subprocess.PIPE)
    p = subprocess.Popen(['mm3d', 'AperiCloud', img_pattern, ori], stdin=echo.stdout)
    return p.wait()


def malt(imlist, ori, zoomf=1, zoomi=None, dirmec='MEC-Malt', seed_img=None, seed_xml=None):
    """
    Run mm3d Malt Ortho.

    :param str|iterable imlist: either a match pattern (e.g., OIS.*tif) or an iterable object of image filenames.
    :param str ori: the orientation directory to use for Malt.
    :param int zoomf: the final Zoom level to use (default: 1)
    :param int zoomi: the initial Zoom level to use (default: not set)
    :param str dirmec: the output MEC directory to create (default: MEC-Malt)
    :param str seed_img: a DEM to pass to Malt as DEMInitImg. Note that if seed_img is set, seed_xml
        must also be set. (default: not used)
    :param str seed_xml: an XML file corresponding to the seed_img (default: not used)
    """
    echo = subprocess.Popen('echo', stdout=subprocess.PIPE)

    if type(imlist) is str:
        matchstr = imlist
    else:
        try:
            matchstr = '|'.join(imlist)
        except TypeError as te:
            raise TypeError(f"imlist is not iterable: {imlist}")

    args = ['mm3d', 'Malt', 'Ortho', matchstr, ori, 'DirMEC={}'.format(dirmec),
            'NbVI=2', 'ZoomF={}'.format(zoomf), 'DefCor=0', 'CostTrans=1', 'EZA=1']

    if zoomi is not None:
        args.append('ZoomI={}'.format(zoomi))

    if seed_img is not None:
        assert seed_xml is not None
        args.append('DEMInitImg=' + seed_img)
        args.append('DEMInitXML=' + seed_xml)

    p = subprocess.Popen(args, stdin=echo.stdout)

    p.wait()


def tawny(dirmec, radiomegal=False):
    """
    Run mm3d Tawny to create an orthomosaic.

    :param str dirmec: the MEC directory to use
    :param bool radiomegal: run Tawny with RadiomEgal=1 (default: False)
    """
    echo = subprocess.Popen('echo', stdout=subprocess.PIPE)

    p = subprocess.Popen(['mm3d', 'Tawny', 'Ortho-{}'.format(dirmec), 'Out=Orthophotomosaic.tif',
                          'RadiomEgal={}'.format(int(radiomegal))], stdin=echo.stdout)
    p.wait()


def block_malt(imlist, nimg=3, ori='Relative', zoomf=8):
    """
    Run mm3d Malt Ortho and mm3d Tawny on successive blocks of images.

    :param iterable imlist: an iterable object of image filenames.
    :param int nimg: the number of images to use in a block (default: 3)
    :param str ori: the name of the orientation directory (e.g., Ori-Relative). (default: Relative)
    :param int zoomf: the final Zoom level to use (default: 8)
    """
    dirmec = 'MEC-' + ori

    inds = range(0, len(imlist) - (nimg - 1), nimg - 1)
    if len(inds) == 1 and len(imlist) > nimg:
        inds = [0, 1]

    for block, ind in enumerate(inds):
        print(imlist[ind:ind + nimg])

        malt(imlist[ind:ind + nimg], ori, dirmec='{}_block{}'.format(dirmec, block), zoomf=zoomf)

        tawny('{}_block{}'.format(dirmec, block))


def bascule(in_gcps, outdir, img_pattern, sub, ori, outori='TerrainRelAuto',
            fn_gcp='AutoGCPs', fn_meas='AutoMeasures'):
    """
    Interface for running mm3d GCPBascule and reading the residuals from the resulting xml file.

    :param pandas.DataFrame in_gcps: a DataFrame with the GCPs that are being input to Campari.
    :param str outdir: the output directory where the AutoGCPs.xml file is saved.
    :param str img_pattern: the match pattern for the images being input to Campari (e.g., "OIS.*tif")
    :param str sub: the name of the block, if multiple blocks are being used (e.g., '_block1'). If not, use ''.
    :param str ori: the name of the orientation directory (e.g., Ori-Relative).
    :param str outori: the name of the output orientation directory (default: TerrainRelAuto).
    :param str fn_gcp: the filename pattern for the GCP file. The file that will be loaded will be
        fn_gcp + sub + '.xml' (e.g., default: AutoGCPs -> AutoGCPs_block0.xml)
    :param str fn_meas: the filename pattern for the measures file. The file that will be loaded will be
        fn_meas + sub + '-S2D.xml' (e.g., default: AutoMeasures -> AutoMeasures_block0-S2D.xml)
    :return: **out_gcps** (*pandas.DataFrame*) -- the input gcps with the updated Bascule residuals.
    """
    fn_gcp = fn_gcp + sub + '.xml'
    fn_meas = fn_meas + sub + '-S2D.xml'

    echo = subprocess.Popen('echo', stdout=subprocess.PIPE)
    p = subprocess.Popen(['mm3d', 'GCPBascule', img_pattern, ori,
                          outori + sub,
                          os.path.join(outdir, fn_gcp),
                          os.path.join(outdir, fn_meas)], stdin=echo.stdout)
    p.wait()

    out_gcps = get_bascule_residuals(os.path.join('Ori-{}{}'.format(outori, sub),
                                                  'Result-GCP-Bascule.xml'), in_gcps)
    return out_gcps


def campari(in_gcps, outdir, img_pattern, sub, dx, ortho_res, allfree=True,
            fn_gcp='AutoGCPs', fn_meas='AutoMeasures', inori='TerrainRelAuto',
            outori='TerrainFinal', homol='Homol'):
    """
    Interface for running mm3d Campari and reading the residuals from the residual xml file.

    :param pandas.DataFrame in_gcps: a DataFrame with the GCPs that are being input to Campari.
    :param str outdir: the output directory where the AutoGCPs.xml file is saved.
    :param str img_pattern: the match pattern for the images being input to Campari (e.g., "OIS.*tif")
    :param str sub: the name of the block, if multiple blocks are being used (e.g., '_block1'). If not, use ''.
    :param int|float dx: the pixel resolution of the reference image.
    :param int|float ortho_res: the pixel resolution of the orthoimage being used.
    :param bool allfree: run Campari with AllFree=1 (True), or AllFree=0 (False). (default: True)
    :param str fn_gcp: the filename pattern for the GCP file. The file that will be loaded will be
        fn_gcp + sub + '.xml' (e.g., default: AutoGCPs -> AutoGCPs_block0.xml)
    :param str fn_meas: the filename pattern for the measures file. The file that will be loaded will be
        fn_meas + sub + '-S2D.xml' (e.g., default: AutoMeasures -> AutoMeasures_block0-S2D.xml)
    :param str inori: the input orientation to Campari (default: Ori-TerrainRelAuto -> TerrainRelAuto)
    :param str outori: the output orientation from Campari (default: Ori-TerrainFinal -> TerrainFinal)
    :param str homol: the Homologue directory to use (default: Homol)
    :return: **out_gcps** (*pandas.DataFrame*) -- the input gcps with the updated Campari residuals.
    """
    echo = subprocess.Popen('echo', stdout=subprocess.PIPE)

    fn_gcp = fn_gcp + sub + '.xml'
    fn_meas = fn_meas + sub + '-S2D.xml'

    p = subprocess.Popen(['mm3d', 'Campari', img_pattern,
                          inori + sub,
                          outori + sub,
                          'GCP=[{},{},{},{}]'.format(os.path.join(outdir, fn_gcp),
                                                     np.abs(dx),
                                                     os.path.join(outdir, fn_meas),
                                                     np.abs(dx / ortho_res)),
                          'SH={}'.format(homol),
                          'AllFree={}'.format(int(allfree))], stdin=echo.stdout)
    p.wait()

    out_gcps = get_campari_residuals('Ori-{}/Residus.xml'.format(outori + sub), in_gcps)
    # out_gcps.dropna(inplace=True)  # sometimes, campari can return no information for a gcp
    return out_gcps


def remove_worst_mesures(fn_meas, ori):
    """
    Remove outlier measures from an xml file, given the output from Campari.

    :param str fn_meas: the filename for the measures file.
    :param str ori: the orientation directory output from Campari (e.g., Ori-TerrainFinal -> TerrainFinal)
    """
    camp_root = ET.parse('Ori-{}/Residus.xml'.format(ori)).getroot()
    auto_root = ET.parse(fn_meas).getroot()

    last_iter = camp_root.findall('Iters')[-1].findall('OneAppui')

    resids_df = pd.DataFrame()

    for ii, appui in enumerate(last_iter):
        resids_df.loc[ii, 'id'] = appui.find('Name').text
        if appui.find('EcartImMoy') is not None:
            resids_df.loc[ii, 'errmoy'] = float(appui.find('EcartImMoy').text)
        if appui.find('EcartImMax') is not None:
            resids_df.loc[ii, 'errmax'] = float(appui.find('EcartImMax').text)
        if appui.find('NameImMax') is not None:
            resids_df.loc[ii, 'immax'] = appui.find('NameImMax').text

    bad_meas = np.abs(resids_df.errmax - resids_df.errmax.median()) > nmad(resids_df.errmax)
    bad_resids = resids_df[bad_meas].copy()

    for im in auto_root.findall('MesureAppuiFlottant1Im'):
        if im.find('NameIm').text in list(bad_resids.immax):
            these_mes = bad_resids[bad_resids.immax == im.find('NameIm').text]
            for pt in im.findall('OneMesureAF1I'):
                if pt.find('NamePt').text in these_mes.id.values:
                    im.remove(pt)

    out_xml = ET.ElementTree(auto_root)
    out_xml.write(fn_meas, encoding="utf-8", xml_declaration=True)


def iterate_campari(gcps, out_dir, match_pattern, subscript, dx, ortho_res, fn_gcp='AutoGCPs', fn_meas='AutoMeasures',
                    rel_ori='Relative', inori='TerrainRelAuto', outori='TerrainFinal', homol='Homol',
                    allfree=True, max_iter=5):
    """
    Run Campari iteratively, refining the orientation by removing outlier GCPs and Measures, based on their fit to the
    estimated camera model.

    :param pandas.DataFrame gcps: a DataFrame with the GCPs that are being input to Campari.
    :param str out_dir: the output directory where the GCP and Measures files are located.
    :param str match_pattern: the match pattern for the images being input to Campari (e.g., "OIS.*tif")
    :param str subscript: the name of the block, if multiple blocks are being used (e.g., '_block1'). If not, use ''.
    :param int|float dx: the pixel resolution of the reference image.
    :param int|float ortho_res: the pixel resolution of the orthoimage being used.
    :param str fn_gcp: the filename pattern for the GCP file. The file that will be loaded will be
        fn_gcp + sub + '.xml' (e.g., default: AutoGCPs -> AutoGCPs_block0.xml)
    :param str fn_meas: the filename pattern for the measures file. The file that will be loaded will be
        fn_meas + sub + '-S2D.xml' (e.g., default: AutoMeasures -> AutoMeasures_block0-S2D.xml)
    :param str rel_ori: the name of the relative orientation to input to GCPBascule (default: Relative -> Ori-Relative + sub)
    :param str inori: the input orientation to Campari (default: Ori-TerrainRelAuto -> TerrainRelAuto)
    :param str outori: the output orientation from Campari (default: Ori-TerrainFinal -> TerrainFinal)
    :param str homol: the Homologue directory to use (default: Homol)
    :param bool allfree: run Campari with AllFree=1 (True), or AllFree=0 (False). (default: True)
    :param int max_iter: the maximum number of iterations to run. (default: 5)
    :return: **gcps** (*pandas.DataFrame*) -- the gcps with updated residuals after the iterative process.
    """
    niter = 0

    gcps = bascule(gcps, out_dir, match_pattern, subscript, rel_ori, fn_gcp=fn_gcp, fn_meas=fn_meas, outori=inori)

    gcps['res_dist'] = np.sqrt(gcps.xres ** 2 + gcps.yres ** 2)

    gcps = campari(gcps, out_dir, match_pattern, subscript, dx, ortho_res,
                   inori=inori, outori=outori, fn_gcp=fn_gcp, fn_meas=fn_meas,
                   allfree=allfree)

    gcps['camp_dist'] = np.sqrt(gcps.camp_xres ** 2 + gcps.camp_yres ** 2)

    while any([np.any(np.abs(gcps.camp_res - gcps.camp_res.median()) > 2 * nmad(gcps.camp_res)),
               np.any(np.abs(gcps.camp_dist - gcps.camp_dist.median()) > 2 * nmad(gcps.camp_dist)),
               gcps.camp_res.max() > 2]) and niter <= max_iter:
        valid_inds = np.logical_and.reduce((np.abs(gcps.camp_res - gcps.camp_res.median()) < 2 * nmad(gcps.camp_res),
                                            gcps.camp_res < gcps.camp_res.max()))
        if np.count_nonzero(valid_inds) < 10:
            break

        gcps = gcps.loc[valid_inds]
        save_gcps(gcps, out_dir, register.get_utm_str(gcps.crs.to_epsg), subscript, fn_gcp=fn_gcp, fn_meas=fn_meas)

        gcps = bascule(gcps, out_dir, match_pattern, subscript, rel_ori, fn_gcp=fn_gcp,
                       fn_meas=fn_meas, outori=inori)
        gcps['res_dist'] = np.sqrt(gcps.xres ** 2 + gcps.yres ** 2)

        gcps = campari(gcps, out_dir, match_pattern, subscript, dx, ortho_res,
                       inori=inori, outori=outori, fn_gcp=fn_gcp, fn_meas=fn_meas,
                       allfree=allfree, homol=homol)

        gcps['camp_dist'] = np.sqrt(gcps.camp_xres ** 2 + gcps.camp_yres ** 2)
        niter += 1

    return gcps


def mask_invalid_els(dir_mec, fn_dem, fn_mask, ori, match_pattern='OIS.*tif', zoomf=1):
    """
    Mask invalid elevations (e.g., water) in a DEM, then re-run the final step of mm3d Malt Ortho to make nicer
    orthophotos.

    :param str dir_mec: the MEC directory (e.g., MEC-Malt) to use
    :param str fn_dem: the filename of the reference DEM
    :param str fn_mask: filename for the mask vector file
    :param str ori: the orientation directory used to run Malt
    :param str match_pattern: the match pattern used to
    :param int zoomf: the final zoom level to run Malt at (default: ZoomF=1)
    """
    zlist = glob(os.path.join(dir_mec, 'Z*.tif'))
    zlist.sort()

    etapes = [int(f.split('_')[1].replace('Num', '')) for f in zlist]

    ind = np.argmax(etapes)
    etape0 = max(etapes)

    fn_auto = os.path.join(dir_mec, 'AutoMask_STD-MALT_Num_{}.tif'.format(etape0 - 1))

    print(fn_auto)
    print(zlist[ind])

    dem = GeoImg(zlist[ind])

    automask = GeoImg(os.path.join(dir_mec, 'AutoMask_STD-MALT_Num_{}.tif'.format(etape0 - 1)))

    shutil.copy(zlist[ind].replace('tif', 'tfw'),
                fn_auto.replace('tif', 'tfw'))

    ref_dem = GeoImg(fn_dem).reproject(dem)

    mask = create_mask_from_shapefile(dem, fn_mask)

    dem.img[mask] = ref_dem.img[mask]
    automask.img[mask] = 1

    automask.write(fn_auto)
    dem.write(zlist[ind])

    echo = subprocess.Popen('echo', stdout=subprocess.PIPE)
    p = subprocess.Popen(['mm3d', 'Malt', 'Ortho', match_pattern, ori, 'DirMEC={}'.format(dir_mec),
                          'NbVI=2', 'MasqImGlob=filtre.tif', 'ZoomF={}'.format(zoomf),
                          'DefCor=1', 'CostTrans=1', 'EZA=1', 'DoMEC=0', 'DoOrtho=1', 'Etape0={}'.format(etape0)],
                         stdin=echo.stdout)
    p.wait()


def save_gcps(in_gcps, outdir, utmstr, sub, fn_gcp='AutoGCPs', fn_meas='AutoMeasures'):
    """
    Save a GeoDataFrame of GCP information to shapefile, txt, and xml formats.

    After running, the following new files will be created:

        - outdir/fn_gcp.shp (+ associated files)
        - outdir/fn_gcp.txt
        - outdir/fn_gcp.xml (output from mm3d GCPConvert)
        - outdir/fn_meas.xml (a file with image locations for each GCP)

    :param GeoDataFrame in_gcps: the gcps GeoDataFrame to save
    :param str outdir: the output directory to save the files to
    :param str utmstr: a UTM string generated by register.get_utm_str()
    :param str sub: the name of the block, if multiple blocks are being used (e.g., '_block1'). If not, use ''.
    :param str fn_gcp: the filename pattern for the GCP file. The file that will be loaded will be
        fn_gcp + sub + '.xml' (e.g., default: AutoGCPs -> AutoGCPs_block0.xml)
    :param str fn_meas: the filename pattern for the measures file. The file that will be loaded will be
        fn_meas + sub + '-S2D.xml' (e.g., default: AutoMeasures -> AutoMeasures_block0-S2D.xml)
    """
    in_gcps.to_file(os.path.join(outdir, fn_gcp + sub + '.shp'))
    write_auto_gcps(in_gcps, sub, outdir, utmstr, outname=fn_gcp)
    echo = subprocess.Popen('echo', stdout=subprocess.PIPE)
    p = subprocess.Popen(['mm3d', 'GCPConvert', 'AppInFile',
                          os.path.join(outdir, fn_gcp + sub + '.txt')], stdin=echo.stdout)
    p.wait()

    auto_root = ET.parse(os.path.join(outdir, fn_meas + sub + '-S2D.xml')).getroot()
    for im in auto_root.findall('MesureAppuiFlottant1Im'):
        for pt in im.findall('OneMesureAF1I'):
            if pt.find('NamePt').text not in in_gcps.id.values:
                im.remove(pt)

    # save AutoMeasures
    out_xml = ET.ElementTree(auto_root)
    out_xml.write(os.path.join(outdir, fn_meas + sub + '-S2D.xml'),
                  encoding="utf-8", xml_declaration=True)


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


# copied from pymmaster
def mosaic_micmac_tiles(filename, dirname='.'):
    """
    Re-stitch images tiled by MicMac.

    :param str filename: MicMac filename to mosaic together
    :param str dirname: Directory containing images to Mosaic (default: .)
    """
    filelist = glob(os.path.sep.join([dirname, '{}_Tile*'.format(filename)]))

    tiled = arrange_tiles(filelist, filename, dirname)
    I, J = tiled.shape

    arr_cols = []
    for j in range(J):
        arr_cols.append(np.concatenate(tiled[:, j], axis=0))

    img = np.concatenate(arr_cols, axis=1)

    imsave(os.path.sep.join([dirname, '{}.tif'.format(filename)]), img)


def arrange_tiles(flist, filename, dirname='.'):
    tmp_inds = [os.path.splitext(f)[0].split('Tile_')[-1].split('_') for f in flist]
    arr_inds = np.array([[int(a) for a in ind] for ind in tmp_inds])
    nrows = arr_inds[:, 1].max() + 1
    ncols = arr_inds[:, 0].max() + 1
    img_arr = np.array(np.zeros((nrows, ncols)), dtype='object')
    for i in range(nrows):
        for j in range(ncols):
            img_arr[i, j] = imread(os.path.sep.join([dirname, '{}_Tile_{}_{}.tif'.format(filename, j, i)]))
    return img_arr


def post_process(projstr, out_name, dirmec, do_ortho=True):
    """
    Apply georeferencing and masking to the final DEM and Correlation images (optionally, the orthomosaic as well).

    Output files are written as follows:
        - DEM: post_processed/{out_name}_Z.tif
        - Hillshade: post_processed/{out_name}_HS.tif
        - Correlation: post_processed/{out_name}_CORR.tif
        - Orthomosaic: post_processed/{out_name}_Ortho.tif

    :param str projstr: A string corresponding to the DEM's CRS that GDAL can use to georeference the rasters.
    :param str out_name: The name that the output files should have.
    :param str dirmec: The MEC directory to process files from (e.g., MEC-Malt)
    :param bool do_ortho: Post-process the orthomosaic in Ortho-{dirmec}, as well. Assumes that you have run
        mm3d Tawny with Out=Orthophotomosaic first.
    """
    # TODO: re-implement this with geoutils/xdem instead of subprocess calls
    os.makedirs('post_processed', exist_ok=True)

    # first, the stuff in MEC
    dem_list = sorted(glob('Z_Num*STD-MALT.tif', root_dir=dirmec))
    level = int(re.findall(r'\d+', dem_list[-1].split('_')[1])[0])
    zoomf = int(re.findall(r'\d+', dem_list[-1].split('_')[2])[0])

    shutil.copy(os.path.join(dirmec, f'Z_Num{level}_DeZoom{zoomf}_STD-MALT.tfw'),
                os.path.join(dirmec, f'Correl_STD-MALT_Num_{level-1}.tfw'))

    shutil.copy(os.path.join(dirmec, f'Z_Num{level}_DeZoom{zoomf}_STD-MALT.tfw'),
                os.path.join(dirmec, f'AutoMask_STD-MALT_Num_{level-1}.tfw'))

    subprocess.Popen(['gdal_translate', '-a_nodata', '0', '-a_srs', projstr,
                      os.path.join(dirmec, f'Correl_STD-MALT_Num_{level-1}.tif'),
                      'tmp_corr.tif']).wait()

    subprocess.Popen(['gdal_translate', '-a_srs', projstr,
                      os.path.join(dirmec, f'Z_Num{level}_DeZoom{zoomf}_STD-MALT.tif'),
                      'tmp_geo.tif']).wait()

    subprocess.Popen(['gdal_translate', '-a_nodata', '0', '-a_srs', projstr,
                      os.path.join(dirmec, f'AutoMask_STD-MALT_Num_{level-1}.tif'),
                      'tmp_mask.tif']).wait()

    subprocess.Popen(['gdal_calc.py', '--quiet', '-A', 'tmp_mask.tif', '-B', 'tmp_geo.tif',
                      '--outfile={}'.format(os.path.join('post_processed', f'{out_name}_Z.tif')),
                      '--calc="B*(A>0)"', '--NoDataValue=-9999']).wait()

    subprocess.Popen(['gdaldem', 'hillshade', os.path.join('post_processed', f'{out_name}_Z.tif'),
                      os.path.join('post_processed', f'{out_name}_HS.tif')]).wait()

    subprocess.Popen(['gdal_calc.py', '--quiet', '-A', 'tmp_corr.tif',
                      '--outfile={}'.format(os.path.join('post_processed', f'{out_name}_CORR.tif')),
                      '--calc="((A.astype(float)-127)/128)*100"', '--NoDataValue=-9999']).wait()

    # clean up the temporary files
    os.remove('tmp_geo.tif')
    os.remove('tmp_corr.tif')

    # now, mask the ortho image(s)
    if do_ortho:
        ortho = os.path.join('Ortho-' + dirmec, 'Orthophotomosaic.tif')

        subprocess.Popen(['gdal_translate', '-a_nodata', '0', '-a_srs', projstr, ortho, 'tmp_ortho.tif']).wait()

        # TODO: re-size the mask to fit the ortho image, if needed
        subprocess.Popen(['gdal_calc.py', '--quiet', '-A', 'tmp_mask.tif', '-B', 'tmp_ortho.tif',
                          '--outfile={}'.format(os.path.join('post_processed', f'{out_name}_Ortho.tif')),
                          '--calc="B*(A>0)"', '--NoDataValue=0']).wait()
        os.remove('tmp_ortho.tif')

    os.remove('tmp_mask.tif')


# converted from bash script
def get_autogcp_locations(ori, meas_file, imlist):
    """
    Find location of automatically-detected control points in individual images using mm3d XYZ2Im.

    :param str ori: The orientation directory name (e.g., Ori-Relative)
    :param str meas_file: The Measures file to find image locations for
    :param list imlist: a list of image names
    """
    nodist = '-'.join([ori, 'NoDist'])

    # copy the orientation directory to a new, "nodist" directory
    shutil.copytree(ori, nodist, dirs_exist_ok=True)

    autocals = glob('AutoCal*.xml', root_dir=nodist)
    for autocal in autocals:
        _remove_distortion_coeffs(os.path.join(nodist, autocal))

    for im in imlist:
        _update_autocal(nodist, im)

        p = subprocess.Popen(['mm3d', 'XYZ2Im', os.path.join(nodist, f'Orientation-{im}.xml'),
                              meas_file, f'NoDist-{im}.txt'])
        p.wait()

        p = subprocess.Popen(['mm3d', 'XYZ2Im', os.path.join(ori, f'Orientation-{im}.xml'),
                              meas_file, f'Auto-{im}.txt'])
        p.wait()


def _remove_distortion_coeffs(fn_xml):
    root = ET.parse(fn_xml).getroot()

    dist_coeffs = root.find('CalibrationInternConique').find('CalibDistortion').find('ModRad').findall('CoeffDist')
    inv_coeffs = root.find('CalibrationInternConique').find('CalibDistortion').find('ModRad').findall('CoeffDistInv')

    for coeff in dist_coeffs + inv_coeffs:
        coeff.text = '0.0'

    tree = ET.ElementTree(root)
    tree.write(fn_xml, encoding="utf-8", xml_declaration=True)


def _update_autocal(ori, im):
    fn_xml = os.path.join(ori, f'Orientation-{im}.xml')
    root = ET.parse(fn_xml).getroot()

    old_autocal = root.find('OrientationConique').find('FileInterne').text
    old_autocal = os.path.join(ori, os.path.basename(old_autocal))

    root.find('OrientationConique').find('FileInterne').text = old_autocal

    tree = ET.ElementTree(root)
    tree.write(fn_xml, encoding="utf-8", xml_declaration=True)
