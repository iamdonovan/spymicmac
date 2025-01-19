"""
spymicmac.micmac is a collection of tools for interfacing with MicMac
"""
import os
from pathlib import Path
import sys
import re
import subprocess
import shutil
import PIL
from itertools import combinations
import numpy as np
from osgeo import gdal
import pandas as pd
import geopandas as gpd
import lxml.etree as etree
import lxml.builder as builder
import difflib
import xml.etree.ElementTree as ET
from glob import glob
from shapely.strtree import STRtree
from shapely.geometry import LineString, MultiPoint
from skimage.io import imread, imsave
from skimage.transform import AffineTransform, SimilarityTransform
from skimage.measure import ransac
import geoutils as gu
from spymicmac import data, matching, register


######################################################################################################################
# MicMac interfaces - write xml files for MicMac to read
######################################################################################################################
def write_neighbour_images(imlist=None, fprints=None, name_field='ID', prefix='OIS-Reech_', file_ext='.tif',
                           dataset='AERIAL_COMBIN', from_homol=False, img_pattern='OIS*.tif', dir_homol='Homol'):
    """
    Write an xml file containing image pairs for processing with Tapioca, using either image footprints or a homologue
    directory.

    :param list imlist: a list of (original) image names to use (e.g., without 'OIS-Reech\_')
    :param GeoDataFrame fprints: a vector dataset of footprint polygons. If not provided, will attempt to download
        metadata from USGS for the images.
    :param str name_field: the field in fprints table that contains the image name
    :param str prefix: the prefix attached to the image name read by Tapioca (default: 'OIS-Reech\_')
    :param str file_ext: the file extension for the images read by Tapioca (default: .tif)
    :param str dataset: the USGS dataset name to search if no footprints are provided (default: AERIAL_COMBIN)
    :param bool from_homol: get a list of pairs based on homologue files (default: False)
    :param str img_pattern: the image pattern to pass to glob to get a list of filenames (default: OIS*.tif)
    :param str dir_homol: the directory where the homologue files are (default: Homol)
    """
    E = builder.ElementMaker()
    NamedRel = E.SauvegardeNamedRel()

    if from_homol:
        pairs = pairs_from_homol(img_pattern, dir_homol)
    else:
        if imlist is None:
            imlist = sorted([fn.strip(prefix).strip(file_ext) for fn in glob(img_pattern)])

        pairs = pairs_from_footprints(imlist=imlist, fprints=fprints, name_field=name_field, prefix=prefix,
                                      file_ext=file_ext, dataset=dataset)

    for pair in pairs:
        this_pair = E.Cple(' '.join(pair))
        NamedRel.append(this_pair)

    tree = etree.ElementTree(NamedRel)
    tree.write('FileImagesNeighbour.xml', pretty_print=True, xml_declaration=True, encoding="utf-8")


def pairs_from_footprints(imlist, fprints=None, name_field='ID', prefix='OIS-Reech_', file_ext='.tif',
                          dataset='AERIAL_COMBIN'):
    """
    Using a list of images and a collection of image footprints, return a list of potential image pairs for processing
    with Tapioca.

    :param list imlist: a list of (original) image names to use (e.g., without 'OIS-Reech_')
    :param GeoDataFrame fprints: a vector dataset of footprint polygons. If not provided, will attempt to download
        metadata from USGS for the images.
    :param str name_field: the field in fprints table that contains the image name
    :param str prefix: the prefix attached to the image name read by Tapioca (default: 'OIS-Reech_')
    :param str file_ext: the file extension for the images read by Tapioca (default: .tif)
    :param str dataset: the USGS dataset name to search if no footprints are provided (default: AERIAL_COMBIN)
    :return: **pairs** (*list*) -- a list of tuples representing image pairs
    """

    if fprints is None:
        fprints = data.get_usgs_footprints(imlist, dataset=dataset)
    else:
        fprints = fprints[fprints[name_field].isin(imlist)]

    fprints.reset_index(inplace=True)  # do this to ensure that strtree indices are correct
    s = STRtree([f for f in fprints['geometry'].values])

    all_pairs = []

    for ind, row in fprints.iterrows():
        fn = row[name_field]
        fp = row['geometry']

        print(fn)

        res = s.query(fp)
        intersects = [fprints.loc[c, 'geometry'] for c in res if fp.intersection(fprints.loc[c, 'geometry']).area > 0]
        fnames = [fprints[name_field][fprints['geometry'] == c].values[0] for c in intersects]
        try:
            fnames.remove(fn)
        except ValueError:
            pass

        all_pairs += [(prefix + fn + file_ext, prefix + fn_match + file_ext) for fn_match in fnames]

    return all_pairs


def _get_pairs(fn_img, dir_homol):
    datlist = sorted(glob('*.dat', root_dir=Path(dir_homol, f"Pastis{fn_img}")))
    return sorted([fn.strip('.dat') for fn in datlist])


def pairs_from_homol(img_pattern='OIS*.tif', dir_homol='Homol'):
    """
    Get a list of image pairs based on homologue files.

    :param str img_pattern: the image pattern to pass to glob to get a list of filenames (default: OIS*.tif)
    :param str dir_homol: the directory where the homologue files are (default: Homol)
    :return: **pairs** (*list*) -- a list of tuples representing image pairs
    """
    imlist = sorted(glob(img_pattern))

    all_pairs = []

    for fn_img in imlist:
        pairs = _get_pairs(fn_img, dir_homol)
        all_pairs += zip(len(pairs) * [fn_img], pairs)

    return all_pairs


def _get_dat_sizes(fn_img, dir_homol):
    datlist = sorted(glob('*.dat', root_dir=Path(dir_homol, f"Pastis{fn_img}")))
    sizes = np.array([os.stat(Path(dir_homol, f"Pastis{fn_img}", fn)).st_size for fn in datlist])
    percs = sizes / sizes.sum()

    size_df = pd.DataFrame({'image': fn_img, 'filename': datlist, 'size': sizes, 'percentage': percs})

    size_df['cumulative'] = size_df.sort_values('percentage', ascending=False)['percentage'].cumsum()

    return size_df


def clean_homol(img_pattern='OIS*.tif', dir_homol='Homol', min_size=None, remove_asymmetric=False, return_df=False):
    """
    Remove spurious homologue files based on a threshold file size.

    :param str img_pattern: the image pattern to pass to glob to get a list of filenames (default: OIS*.tif)
    :param str dir_homol: the directory where the homologue files are (default: Homol)
    :param int min_size: the size, in bytes, to use as a threshold for removing file (default: calculated from all files)
    :param bool remove_asymmetric: remove asymmetric homologue files (pairs where only one image in the pair "sees" the
        other one) (default: False)
    :param bool return_df: return a DataFrame of all homologue files, rather than removing them (default: False)
    """
    imlist = sorted(glob(img_pattern))

    dat_sizes = pd.concat([_get_dat_sizes(fn, dir_homol) for fn in imlist], ignore_index=True)
    dat_sizes['symmetric'] = True

    for ind, row in dat_sizes.iterrows():
        dat_sizes.loc[ind, 'symmetric'] = row['image'] + '.dat' in dat_sizes.loc[
            dat_sizes['image'] == row['filename'].split('.dat')[0], 'filename'].to_list()

    if return_df:
        return dat_sizes

    for fn_img, sizes in dat_sizes.groupby('image'):
        if min_size is None:
            ind = (sizes.sort_values('cumulative')['cumulative'] > 0.95).idxmax()
            this_cutoff = min(250, sizes.loc[ind, 'size'] + 1)
        else:
            this_cutoff = min_size

        for fn in sizes.loc[sizes['size'] < this_cutoff, 'filename']:
            os.remove(Path(dir_homol, f"Pastis{fn_img}", fn))

        if remove_asymmetric:
            for fn in sizes.loc[(~sizes['symmetric']) & (sizes['size'] > this_cutoff), 'filename']:
                os.remove(Path(dir_homol, f"Pastis{fn_img}", fn))


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


def get_im_meas(points_df, E, name='gcp', x='im_col', y='im_row'):
    """
    Populate an lxml.builder.ElementMaker object with GCP image locations, for writing to xml files.

    :param pandas.DataFrame points_df: a DataFrame with the points to find image locations for.
    :param lxml.builder.ElementMaker E: an ElementMaker object for writing to the xml file.
    :param str name: the column name in points_df corresponding to the point name [gcp]
    :param str x: the column name in points_df corresponding to the image x location [im_col]
    :param str y: the column name in points_df corresponding to the image y location [im_row]
    :return: **pt_els** (*list*) -- a list of ElementMaker objects corresponding to each GCP image location.
    """
    pt_els = []
    for ind, row in points_df.iterrows():
        this_mes = E.OneMesureAF1I(
                        E.NamePt(row[name]),
                        E.PtIm('{} {}'.format(row[x], row[y]))
                        )
        pt_els.append(this_mes)
    return pt_els


def parse_im_meas(fn_meas):
    """
    Read an xml file with image locations into a pandas DataFrame.

    :param fn_meas: the name of the measures file to read.
    :return: **gcp_df** (*pandas.DataFrame*) -- a DataFrame with gcp names and image locations.
    """
    root = ET.parse(fn_meas).getroot()
    if root.tag == 'MesureAppuiFlottant1Im':
        measures = root
    else:
        measures = root.findall('MesureAppuiFlottant1Im')

    meas_df = pd.DataFrame()

    if type(measures) == list:
        for img_mes in measures:
            this_df = pd.DataFrame()
            for ind, mes in enumerate(img_mes.findall('OneMesureAF1I')):
                this_df.loc[ind, 'image'] = img_mes.find('NameIm').text
                this_df.loc[ind, 'name'] = mes.find('NamePt').text
                pt = mes.find('PtIm').text.split()
                this_df.loc[ind, 'i'] = float(pt[1])
                this_df.loc[ind, 'j'] = float(pt[0])

            meas_df = pd.concat([meas_df, this_df], ignore_index=True)
    else:
        for ind, mes in enumerate(measures.findall('OneMesureAF1I')):
            meas_df.loc[ind, 'name'] = mes.find('NamePt').text
            pt = mes.find('PtIm').text.split()
            meas_df.loc[ind, 'i'] = float(pt[1])
            meas_df.loc[ind, 'j'] = float(pt[0])

    return meas_df


def write_measures_im(meas_df, fn_img):
    """
    Create a MeasuresIm xml file for an image.

    :param DataFrame meas_df: a DataFrame of image measures, with [gcp, im_col, im_row] columns
    :param str fn_img: the filename of the image.
    :return:
    """
    os.makedirs('Ori-InterneScan', exist_ok=True)

    # write the measures
    E = builder.ElementMaker()
    ImMes = E.MesureAppuiFlottant1Im(E.NameIm(fn_img))

    pt_els = get_im_meas(meas_df, E)
    for p in pt_els:
        ImMes.append(p)

    outxml = E.SetOfMesureAppuisFlottants(ImMes)

    tree = etree.ElementTree(outxml)
    tree.write(os.path.join('Ori-InterneScan', 'MeasuresIm-' + fn_img + '.xml'), pretty_print=True,
               xml_declaration=True, encoding="utf-8")


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


def create_measurescamera_xml(fn_csv, ori='InterneScan', translate=False, name='gcp', x='im_col', y='im_row'):
    """
    Create a MeasuresCamera.xml file from a csv of fiducial marker locations.

    :param str|DataFrame fn_csv: the filename of the CSV file, or a pandas DataFrame.
    :param str ori: the Ori directory to write the MeasuresCamera.xml file to. Defaults to (Ori-)InterneScan.
    :param bool translate: translate coordinates so that the origin is the upper left corner, rather than the principal
        point
    :param str name: the column name in the csv file corresponding to the point name [gcp]
    :param str x: the column name in the csv file corresponding to the image x location [im_col]
    :param str y: the column name in the csv file corresponding to the image y location [im_row]
    """
    assert type(fn_csv) in [pd.core.frame.DataFrame, str], "fn_csv must be one of [str, DataFrame]"
    if isinstance(fn_csv, str):
        fids = pd.read_csv(fn_csv)
    else:
        fids = fn_csv

    # if coordinates are relative to the principal point,
    # convert them to be relative to the upper left corner
    if translate:
        fids[x] = fids[x] - fids[x].min()
        fids[y] = -fids[y] - min(-fids[y])

    E = builder.ElementMaker()
    ImMes = E.MesureAppuiFlottant1Im(E.NameIm('Glob'))

    pt_els = get_im_meas(fids, E, name=name, x=x, y=y)
    for p in pt_els:
        ImMes.append(p)

    outxml = E.SetOfMesureAppuisFlottants(ImMes)
    tree = etree.ElementTree(outxml)

    os.makedirs(f'Ori-{ori}', exist_ok=True)
    tree.write(os.path.join(f'Ori-{ori}', 'MeasuresCamera.xml'), pretty_print=True,
               xml_declaration=True, encoding="utf-8")


def estimate_measures_camera(approx, pairs, ori='InterneScan', scan_res=2.5e-5, how='mean', write_xml=True):
    """
    Use a set of located fiducial markers to create a MeasuresCamera file using the average location of each fiducial
    marker.

    :param approx: a dataframe of approximate fiducial marker locations.
    :param pairs: a list of pairs of co-linear fiducial markers
    :param str ori: The Ori- directory containing the MeasuresIm files (default: InterneScan)
    :param float scan_res: the scanning resolution of the images in m (default: 2.5e-5; 25 µm)
    :param str how: what average to use for the output locations. Must be one of [mean, median].
    :param bool write_xml: write the MeasuresCamera.xml file in addition to a CSV (default: True)
    """
    assert how in ['mean', 'median'], "how must be one of [mean, median]"

    meas_list = sorted(glob('MeasuresIm*.xml', root_dir=f'Ori-{ori}'))

    all_meas = []

    for fn_meas in meas_list:
        meas = parse_im_meas(os.path.join(f'Ori-{ori}', fn_meas)).set_index('name')

        collinear = [LineString(meas.loc[p, ['j', 'i']].values) for p in pairs]

        for ind, pair in enumerate(pairs):
            meas.loc[pair, ['collim_dist']] = collinear[ind].length

        scale = np.mean([c.length for c in collinear])

        meas['j'] = meas['j'] / scale
        meas['i'] = meas['i'] / scale

        ppx, ppy = _meas_center(meas, pairs)
        meas['j'] -= ppx
        meas['i'] -= ppy

        model = AffineTransform()
        joined = meas.join(approx.set_index('name'), lsuffix='_img', rsuffix='_cam')
        model.estimate(joined[['j_img', 'i_img']].values, joined[['j_cam', 'i_cam']].values)

        meas['resid'] = model.residuals(joined[['j_img', 'i_img']].values, joined[['j_cam', 'i_cam']].values)

        noscale = AffineTransform(translation=model.translation, rotation=model.rotation, shear=model.shear)
        rot = noscale(meas[['j', 'i']].values)
        # rot = model(meas[['j', 'i']].values)

        meas['j'] = rot[:, 0]
        meas['i'] = rot[:, 1]

        all_meas.append(meas.reset_index())

    all_meas = pd.concat(all_meas, ignore_index=True)

    if how == 'mean':
        avg_meas = all_meas.groupby('name').mean(numeric_only=True)
    else:
        avg_meas = all_meas.groupby('name').median(numeric_only=True)

    joined = avg_meas.join(approx.set_index('name'), lsuffix='_img', rsuffix='_cam')
    model = AffineTransform()
    model.estimate(joined[['j_img', 'i_img']].values, joined[['j_cam', 'i_cam']].values)

    noscale = AffineTransform(translation=model.translation, rotation=model.rotation)

    rot = noscale(meas[['j', 'i']].values)

    avg_meas['j'] = rot[:, 0]
    avg_meas['i'] = rot[:, 1]

    avg_meas['j'] -= avg_meas['j'].min()
    avg_meas['i'] -= avg_meas['i'].min()

    scale = avg_meas['collim_dist'].mean()

    avg_meas['j'] *= scale * scan_res * 1000  # convert from m to mm
    avg_meas['i'] *= scale * scan_res * 1000  # convert from m to mm

    all_meas.to_csv('AllMeasures.csv', index=False)
    avg_meas.to_csv('AverageMeasures.csv')

    if write_xml:
        create_measurescamera_xml('AverageMeasures.csv', ori=ori, translate=False, name='name', x='j', y='i')


def _meas_center(meas, pairs):
    collims = [LineString(meas.loc[p, ['j', 'i']].values) for p in pairs]
    pp = MultiPoint([a.intersection(b) for a, b in list(combinations(collims, 2))]).centroid

    return pp.x, pp.y


def generate_multicam_csv(patterns=None, prefix='OIS-Reech_', fn_out='camera_defs.csv',
                          name='', short_name='', film_size='', focal=''):
    """
    Create a CSV file with camera parameters than can be read by create_localchantier_xml() to use images acquired by
    multiple cameras.

    Can be used to create a blank CSV template to be filled out manually, or generated using the optional function
    arguments.

    :param patterns: a list of filename patterns corresponding to each camera [None]
    :param str prefix: an optional prefix to add to the matching patterns [OIS-Reech_]
    :param str fn_out: the name of the CSV file to create [camera_defs.csv]
    :param name: the name to give each camera. Must be unique.
    :param short_name: the "short name" description of each camera. Does not need to be unique.
    :param film_size: the size (width, height in mm) of the frame for each camera. Can be a list of tuples or a str.
    :param focal: the focal length of each camera, in mm.
    """
    cameras = pd.DataFrame()

    if patterns is None:
        cameras['pattern'] = ''
    else:
        patterns = [p + '.*' for p in patterns if '.*' not in p]
        patterns = [prefix + p for p in patterns if prefix not in p]

        cameras['pattern'] = patterns

    cameras['name'] = name
    cameras['short_name'] = short_name

    if not isinstance(film_size, str):
        film_size = [','.join([str(p) for p in pp]) for pp in film_size]

    cameras['film_size'] = film_size
    cameras['focal'] = focal

    cameras.to_csv(fn_out, index=False)


def create_localchantier_xml(name='KH9MC', short_name='KH-9 Hexagon Mapping Camera', film_size=(460, 220),
                             pattern='.*', focal=304.8, add_sfs=False, cam_csv=None):
    """
    Create a MicMac-LocalChantierDescripteur.xml file for a given camera. Default is the KH-9 Hexagon Mapping Camera.

    :param str name: The name to use for the camera [KH9MC]
    :param str short_name: A short description of the camera [KH-9 Hexagon Mapping Camera]
    :param array-like film_size: the film size (width, height) in mm [460, 220]
    :param str pattern: the matching pattern to use for the images [.*]
    :param float focal: the nominal focal length, in mm [304.8]
    :param bool add_sfs: use SFS to help find tie points in low-contrast images [False]
    :param str cam_csv: the CSV file containing parameters for multiple cameras [None]
    """
    E = builder.ElementMaker()

    chantier = E.ChantierDescripteur()
    cam_db = E.LocCamDataBase()
    cam_assocs = E.KeyedNamesAssociations()
    foc_assocs = E.KeyedNamesAssociations()

    if cam_csv is not None:
        cameras = pd.read_csv(cam_csv)
        for ind, cam in cameras.iterrows():
            width, height = [p.strip() for p in cam['film_size'].split(',')]

            cam_db.append(
                E.CameraEntry(
                    E.Name(cam['name']),
                    E.SzCaptMm(f"{width} {height}"),
                    E.ShortName(cam['short_name'])
                )
            )

            cam_assocs.append(
                E.Calcs(
                    E.Arrite('1 1'),
                    E.Direct(
                        E.PatternTransform(cam['pattern']),
                        E.CalcName(cam['name'])
                    )
                )
            )

            foc_assocs.append(
                E.Calcs(
                    E.Arrite('1 1'),
                    E.Direct(
                        E.PatternTransform(cam['pattern']),
                        E.CalcName(f"{cam['focal']}")
                    )
                )
            )

    else:
        width, height = film_size
        cam_db.append(
            E.CameraEntry(
                E.Name(name),
                E.SzCaptMm(f"{width} {height}"),
                E.ShortName(short_name)
            )
        )

        cam_assocs.append(
            E.Calcs(
                E.Arrite('1 1'),
                E.Direct(
                    E.PatternTransform(pattern),
                    E.CalcName(name)
                )
            )
        )

        foc_assocs.append(
            E.Calcs(
                E.Arrite('1 1'),
                E.Direct(
                    E.PatternTransform(pattern),
                    E.CalcName(f"{focal}")
                )
            )
        )

    cam_assocs.append(
        E.Key('NKS-Assoc-STD-CAM')
    )

    foc_assocs.append(
        E.Key('NKS-Assoc-STD-FOC')
    )

    for item in [cam_db, cam_assocs, foc_assocs]:
        chantier.append(item)

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


def write_image_mesures(imlist, gcps, outdir='.', sub='', ort_dir='Ortho-MEC-Relative', outname='AutoMeasures'):
    """
    Create a Measures-S2D.xml file (row, pixel) for each GCP in each image from a list of image names.

    :param list imlist: a list of image names.
    :param pandas.DataFrame gcps: a DataFrame of GCPs.
    :param str outdir: the output directory to save the files to.
    :param str sub: the name of the block, if multiple blocks are being used (e.g., '_block1').
    :param str ort_dir: the Ortho-MEC directory where the images are located (default: Ortho-MEC-Relative)
    :param str outname: the base name of the file to write (default: AutoMeasures)
    """
    E = builder.ElementMaker()
    MesureSet = E.SetOfMesureAppuisFlottants()

    for im in imlist:
        print(im)
        img_geo = gu.Raster(os.path.join(ort_dir, 'Ort_' + im))
        impts = pd.read_csv('Auto-{}.txt'.format(im), sep=' ', names=['j', 'i'])
        # impts_nodist = pd.read_csv('NoDist-{}.txt'.format(im), sep=' ', names=['j', 'i'])

        footprint = (img_geo > 0).polygonize().ds.union_all()
        valid_pts = footprint.contains(gpd.points_from_xy(gcps.rel_x, gcps.rel_y))

        # valid_pts = get_valid_image_points(img.shape, impts, impts_nodist)
        # valid_pts = np.logical_and.reduce([xmin <= gcps.rel_x, gcps.rel_x < xmax,
        #                                    ymin <= gcps.rel_y, gcps.rel_y < ymax])

        if np.count_nonzero(valid_pts) == 0:
            continue

        this_im_mes = E.MesureAppuiFlottant1Im(E.NameIm(im))

        for i, (ind, row) in enumerate(impts[valid_pts].iterrows()):
            this_mes = E.OneMesureAF1I(E.NamePt(gcps.iloc[ind]['id']),
                                       E.PtIm('{} {}'.format(row.j, row.i)))
            this_im_mes.append(this_mes)

        MesureSet.append(this_im_mes)

    tree = etree.ElementTree(MesureSet)
    tree.write(os.path.join(outdir, '{}{}-S2D.xml'.format(outname, sub)),
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
    camp_dist = []

    for a in last_iter:
        try:
            camp_x.append(float(a.find('EcartFaiscTerrain').text.split()[0]))
            camp_y.append(float(a.find('EcartFaiscTerrain').text.split()[1]))
            camp_z.append(float(a.find('EcartFaiscTerrain').text.split()[2]))
            camp_dist.append(float(a.find('DistFaiscTerrain').text))
        except AttributeError:
            camp_x.append(np.nan)
            camp_y.append(np.nan)
            camp_z.append(np.nan)
            camp_dist.append(np.nan)

    for data_ in zip(camp_gcp_names, err_max):
        gcp_df.loc[gcp_df.id == data_[0], 'camp_res'] = data_[1]

    for data_ in zip(camp_gcp_names, camp_x):
        gcp_df.loc[gcp_df.id == data_[0], 'camp_xres'] = data_[1]

    for data_ in zip(camp_gcp_names, camp_y):
        gcp_df.loc[gcp_df.id == data_[0], 'camp_yres'] = data_[1]

    for data_ in zip(camp_gcp_names, camp_z):
        gcp_df.loc[gcp_df.id == data_[0], 'camp_zres'] = data_[1]

    for data_ in zip(camp_gcp_names, camp_dist):
        gcp_df.loc[gcp_df.id == data_[0], 'camp_dist'] = data_[1]

    return gcp_df


def get_tapas_residuals(ori):
    """
    Read the image residuals output from Tapas.

    :param str ori: the name of the Ori directory to read the residuals from (e.g., 'Relative' for Ori-Relative)
    :return: img_df (DataFrame) -- a DataFrame with image names and residuals
    """
    root = ET.parse(os.path.join(f'Ori-{ori}', 'Residus.xml'))
    last = root.findall('Iters')[-1]

    img_df = pd.DataFrame()
    for ind, img in enumerate(last.findall('OneIm')):
        img_df.loc[ind, 'name'] = img.find('Name').text
        img_df.loc[ind, 'res'] = float(img.find('Residual').text)
        img_df.loc[ind, 'perc_ok'] = float(img.find('PercOk').text)
        img_df.loc[ind, 'npts'] = int(img.find('NbPts').text)
        img_df.loc[ind, 'nmult'] = int(img.find('NbPtsMul').text)

    return img_df


def find_empty_homol(imlist=None, dir_homol='Homol', pattern='OIS*.tif'):
    """
    Search through a Homol directory to find images without any matches, then move them to a new directory called
    'EmptyMatch'

    :param list|None imlist: a list of images in the current directory. If None, uses pattern to find images.
    :param str dir_homol: the Homol directory to search in (default: Homol)
    :param str pattern: the search pattern to use to find images (default: OIS*.tif)
    """

    if imlist is None:
        imlist = glob(pattern)

    empty = [fn_img for fn_img in imlist if len(_get_homol(fn_img, dir_homol)) == 0]

    os.makedirs('EmptyMatch', exist_ok=True)

    for fn_img in empty:
        print(f'{fn_img} -> EmptyMatch/{fn_img}')
        shutil.move(fn_img, 'EmptyMatch')


def _get_homol(fn_img, dir_homol='Homol'):
    if not os.path.exists(os.path.join(dir_homol, 'Pastis' + fn_img)):
        return []
    else:
        return sorted([h.split('.dat')[0] for h in glob('*.dat', root_dir=os.path.join(dir_homol, 'Pastis' + fn_img))])


# adapted from the fantastic answer provided by
# user matias-thayer and edited by user redbeam_
# at https://stackoverflow.com/a/50639220
def _get_connected_block(img, seen, hdict):
    result = []
    imgs = set([img])

    while imgs:
        img = imgs.pop()
        seen.add(img)
        imgs.update(set(hdict[img]) - seen)
        result.append(img)

    return result, seen


# adapted from the fantastic answer provided by
# user matias-thayer and edited by user redbeam_
# at https://stackoverflow.com/a/50639220
def find_connected_blocks(pattern='OIS*.tif', dir_homol='Homol'):
    """
    Find connected blocks of images.

    :param str pattern: the search pattern to use to get image names
    :param str dir_homol: the Homologue directory to use to determine what images are connected
    :return: blocks -- a list containing lists of connected blocks of images
    """
    imlist = glob(pattern)
    homols = [[fn for fn in _get_homol(fn_img, dir_homol) if fn in imlist] for fn_img in imlist]

    hdict = dict(zip(imlist, homols))

    seen = set()
    blocks = []

    for img in hdict:
        if img not in seen:
            block, seen = _get_connected_block(img, seen, hdict)
            blocks.append(sorted(block))

    return blocks


def separate_blocks(pattern='OIS*.tif', dir_homol='Homol', min_size=2):
    """
    Based on homologous points, find connected blocks of images and then separate the files into sub-folders.
    Moves files from {dir_homol} and Pastis, along with the image files.

    :param str pattern: the search pattern to use to get image names (default: OIS*.tif)
    :param str dir_homol: the Homologue directory to use to determine what images are connected (default: Homol)
    :param int min_size: the minimum number of images to be considered a block (default: 2)
    """
    # get connected blocks
    blocks = find_connected_blocks(pattern, dir_homol)

    # find single unconnected images
    singles = sorted([im for imgs in blocks for im in imgs if len(imgs) < min_size])
    blocks = [b for b in blocks if len(b) >= min_size]

    # make directory
    os.makedirs('singles', exist_ok=True)

    for fn_sin in singles:
        shutil.move(fn_sin, 'singles')

    for num, block in enumerate(blocks):

        os.makedirs(f"Block{num}", exist_ok=True)
        os.makedirs(Path(f"Block{num}", dir_homol), exist_ok=True)
        os.makedirs(Path(f"Block{num}", 'Pastis'), exist_ok=True)

        for fn_img in block:
            # move homol files
            if os.path.exists(Path(dir_homol, f"Pastis{fn_img}")):
                shutil.move(Path(dir_homol, f"Pastis{fn_img}"), Path(f"Block{num}", dir_homol, f"Pastis{fn_img}"))

            # move pastis files
            for fn_pas in glob(f"*{fn_img}*", root_dir='Pastis'):
                shutil.move(Path('Pastis', fn_pas), Path(f"Block{num}", 'Pastis', fn_pas))

            # move image
            if os.path.exists(fn_img):
                shutil.move(fn_img, f"Block{num}")

        # copy xml files if they exist
        if os.path.exists('MicMac-LocalChantierDescripteur.xml'):
            shutil.copy('MicMac-LocalChantierDescripteur.xml', f"Block{num}")


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


def _generate_glob(fn_ids):
    E = builder.ElementMaker()
    ImMes = E.MesureAppuiFlottant1Im(E.NameIm('Glob'))
    SpGlob = E.SetPointGlob()

    gcp_df = pd.DataFrame()
    with open(fn_ids, 'r') as f:
        fids = [l.strip() for l in f.readlines()]

    for fid in fids:
        pt_glob = E.PointGlob(E.Type('eNSM_Pts'),
                              E.Name(fid),
                              E.LargeurFlou('0'),
                              E.NumAuto('0'),
                              E.SzRech('-1'))
        SpGlob.append(pt_glob)

    tree = etree.ElementTree(SpGlob)
    tree.write('Tmp-SL-Glob.xml', pretty_print=True, xml_declaration=True, encoding="utf-8")


def batch_saisie_fids(imlist, flavor='qt', fn_cam=None, clean=True, gamma=None):
    """
    Run SaisieAppuisInit to locate the fiducial markers for a given list of images.

    :param list imlist: the list of image filenames.
    :param str flavor: which version of SaisieAppuisInit to run. Must be one of [qt, og] (default: qt)
    :param str fn_cam: the filename for the MeasuresCamera.xml file (default: Ori-InterneScan/MeasuresCamera.xml)
    :param bool clean: remove any image files in Tmp-SaisieAppuis
    :param float gamma: Gamma adjustment value for Saisie (default: 1.0)
    """
    assert flavor in ['qt', 'og'], "flavor must be one of [qt, og]"

    os.makedirs('Ori-InterneScan', exist_ok=True)
    os.makedirs('Tmp-SaisieAppuis', exist_ok=True)

    if fn_cam is None:
        fn_cam = os.path.join('Ori-InterneScan', 'MeasuresCamera.xml')

    if os.path.exists(fn_cam):
        measures_cam = parse_im_meas(fn_cam)
        with open('id_fiducials.txt', 'w') as f:
            for fid in measures_cam['name']:
                print(fid, file=f)
        _generate_glob('id_fiducials.txt')
    else:
        try:
            _generate_glob('id_fiducials.txt')
        except FileNotFoundError as e:
            raise FileNotFoundError('id_fiducials.txt not found. Please specify fn_cam, '
                                    'or ensure that id_fiducials.txt exists in the current directory.')

    if flavor == 'qt':
        saisie = 'SaisieAppuisInitQT'
    else:
        saisie = 'SaisieAppuisInit'

    for fn_img in imlist:
        if os.path.exists(os.path.join('Ori-InterneScan', f'MeasuresIm-{fn_img}.xml')):
            if clean:
                tmplist = glob('*' + fn_img + '*', root_dir='Tmp-SaisieAppuis')
                for fn_tmp in tmplist:
                    os.remove(os.path.join('Tmp-SaisieAppuis', fn_tmp))

            shutil.copy(os.path.join('Ori-InterneScan', f'MeasuresIm-{fn_img}.xml'),
                        f'MeasuresIm-{fn_img}-S2D.xml')

            shutil.copy('Tmp-SL-Glob.xml',
                        os.path.join('Tmp-SaisieAppuis', f'Tmp-SL-Glob-MeasuresIm-{fn_img}.xml'))

        saisie_args = ['mm3d', saisie, fn_img, 'NONE', 'id_fiducials.txt', f'MeasuresIm-{fn_img}.xml']

        if gamma is not None:
            saisie_args.append(f"Gama={gamma}")

        p = subprocess.Popen(saisie_args)
        p.wait()

        shutil.move(f'MeasuresIm-{fn_img}-S2D.xml', os.path.join('Ori-InterneScan', f'MeasuresIm-{fn_img}.xml'))
        os.remove(f'MeasuresIm-{fn_img}-S3D.xml')

    os.remove('Tmp-SL-Glob.xml')


def tapioca(img_pattern='OIS.*tif', res_low=400, res_high=1200, fn_neighbours=None):
    """
    Run mm3d Tapioca

    :param str img_pattern: The image pattern to pass to Tapioca (default: OIS.*tif)
    :param int res_low: the size of the largest image axis, in pixels, for low-resolution matching (default: 400)
    :param int res_high: the size of the largest image axis, in pixels, for high-resolution matching (default: 1200)
    :param str fn_neighbours:
    """
    if os.name == 'nt':
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE, shell=True)
    else:
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE)

    if fn_neighbours is None:
        args = ['mm3d', 'Tapioca', 'MulScale', img_pattern, str(res_low), str(res_high)]
    else:
        args = ['mm3d', 'Tapioca', 'File', fn_neighbours, str(res_high)]

    p = subprocess.Popen(args, stdin=echo.stdout)

    return p.wait()


def martini(img_pattern='OIS.*tif', in_ori=None, ori_out=None):
    """
    Run mm3d Martini, which provides a quick way to orient images without solving for camera parameters.

    :param str img_pattern: The image pattern to pass to Martini (default: OIS.*tif)
    :param str in_ori: the orientation directory to use to initialize the calibration (default: None)
    :param str ori_out: the name of the output orientation directory (default: Martini)
    """
    if os.name == 'nt':
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE, shell=True)
    else:
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE)

    args = ['mm3d', 'Martini', img_pattern]

    if in_ori is not None:
        args.append(f"InOri={in_ori}")

    if ori_out is not None:
        args.append(f"OriOut={ori_out}")

    p = subprocess.Popen(args, stdin=echo.stdout)

    return p.wait()


def tapas(cam_model, ori_out=None, img_pattern='OIS.*tif', in_cal=None, in_ori=None,
          lib_foc=True, lib_pp=True, lib_cd=True):
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
    :param str in_ori: a set of orientations to initialize the calibration (default: None)
    :param bool lib_foc: allow the focal length to be calibrated (default: True)
    :param bool lib_pp: allow the principal point to be calibrated (default: True)
    :param bool lib_cd: allow the center of distortion to be calibrated (default: True)
    """
    if os.name == 'nt':
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE, shell=True)
    else:
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE)

    args = ['mm3d', 'Tapas', cam_model, img_pattern,
            'LibFoc={}'.format(int(lib_foc)), 'LibPP={}'.format(int(lib_pp)),
            'LibCD={}'.format(int(lib_cd))]

    if ori_out is not None:
        args.append('Out=' + ori_out)

    if in_cal is not None:
        args.append('InCal=' + in_cal)

    if in_ori is not None:
        args.append('InOri=' + in_ori)

    p = subprocess.Popen(args, stdin=echo.stdout)

    return p.wait()


def apericloud(ori, img_pattern='OIS.*tif', fn_out=None, with_points=True):
    """
    Run mm3d AperiCloud to create a point cloud layer

    :param str ori: the input orientation to use
    :param str img_pattern: the image pattern to pass to AperiCloud (default: OIS.*tif)
    :param str fn_out: the output filename (default: AperiCloud_{ori}.ply)
    :param bool with_points: display the point cloud (default: True)
    """
    if os.name == 'nt':
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE, shell=True)
    else:
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE)

    args = ['mm3d', 'AperiCloud', img_pattern, ori]

    if fn_out is not None:
        args.append(f"Out={fn_out}")

    if not with_points:
        args.append(f"WithPoints=0")

    p = subprocess.Popen(args, stdin=echo.stdout)

    return p.wait()


def malt(imlist, ori, zoomf=1, zoomi=None, dirmec='MEC-Malt', seed_img=None, seed_xml=None,
         resol_terr=None, resol_ort=None, cost_trans=None, szw=None, regul=None, do_ortho=True, do_mec=True):
    """
    Run mm3d Malt Ortho.

    :param str|iterable imlist: either a match pattern (e.g., OIS.*tif) or an iterable object of image filenames.
    :param str ori: the orientation directory to use for Malt.
    :param int zoomf: the final Zoom level to use (default: 1)
    :param int zoomi: the initial Zoom level to use (default: not set)
    :param str dirmec: the output MEC directory to create (default: MEC-Malt)
    :param str seed_img: a DEM to pass to Malt as DEMInitImg. Note that if seed_img is set, seed_xml
        must also be set. If used, it is recommended to set zoomi to be approximately equal to the DEM resolution -
        i.e., if the ortho resolution is 5 m and the seed DEM is 20 m, ZoomI should be 4. (default: not used)
    :param str seed_xml: an XML file corresponding to the seed_img (default: not used)
    :param float resol_terr: the resolution of the output DEM, in ground units (default: computed by mm3d)
    :param float resol_ort: the resolution of the ortho images, relative to the output DEM - e.g., resol_ort=1 means
        the DEM and Orthoimage have the same resolution (default: 2.0)
    :param float cost_trans: cost to change from correlation to decorrelation (default: 2.0)
    :param int szw: the half-size of the correlation window to use - e.g., szw=1 means a 3x3 correlation window.
        (default: 2)
    :param float regul: the regularization factor to use. Lower values mean higher potential variability between
        adjacent pixels, higher values (up to 1) mean smoother outputs (default: 0.05)
    :param bool do_ortho: whether to generate the orthoimages (default: True)
    :param bool do_mec: whether to generate an output DEM (default: True)
    """
    if os.name == 'nt':
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE, shell=True)
    else:
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE)

    if type(imlist) is str:
        matchstr = imlist
    else:
        try:
            matchstr = '|'.join(imlist)
        except TypeError as te:
            raise TypeError(f"imlist is not iterable: {imlist}")

    args = ['mm3d', 'Malt', 'Ortho', matchstr, ori, f'DirMEC={dirmec}',
            'NbVI=2', f'ZoomF={zoomf}', 'DefCor=0', 'EZA=1',
            f'DoOrtho={int(do_ortho)}', f'DoMEC={int(do_mec)}']

    if zoomi is not None:
        args.append(f'ZoomI={zoomi}')

    if seed_img is not None:
        assert seed_xml is not None
        args.append(f'DEMInitIMG={seed_img}')
        args.append(f'DEMInitXML={seed_xml}')

    if resol_terr is not None:
        args.append(f'ResolTerrain={resol_terr}')

    if resol_ort is not None:
        args.append(f'ResolOrtho={resol_ort}')

    if cost_trans is not None:
        args.append(f'CostTrans={cost_trans}')

    if szw is not None:
        args.append(f'SzW={szw}')

    if regul is not None:
        args.append(f'Regul={regul}')

    p = subprocess.Popen(args, stdin=echo.stdout)

    return p.wait()


def tawny(dirmec, radiomegal=False):
    """
    Run mm3d Tawny to create an orthomosaic.

    :param str dirmec: the MEC directory to use
    :param bool radiomegal: run Tawny with RadiomEgal=1 (default: False)
    """
    if os.name == 'nt':
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE, shell=True)
    else:
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE)

    p = subprocess.Popen(['mm3d', 'Tawny', 'Ortho-{}'.format(dirmec), 'Out=Orthophotomosaic.tif',
                          'RadiomEgal={}'.format(int(radiomegal))], stdin=echo.stdout)
    return p.wait()


def block_malt(imlist, nimg=3, ori='Relative', zoomf=8):
    """
    Run mm3d Malt Ortho and mm3d Tawny on successive blocks of images.

    :param iterable imlist: an iterable object of image filenames, or an iterable object of lists of image filenames
    :param int nimg: the number of images to use in a block (default: 3)
    :param str ori: the name of the orientation directory (e.g., Ori-Relative). (default: Relative)
    :param int zoomf: the final Zoom level to use (default: 8)
    """
    dirmec = 'MEC-' + ori

    if type(imlist[0]) is str:
        inds = range(0, len(imlist) - (nimg - 1), nimg - 1)
        if len(inds) == 1 and len(imlist) > nimg:
            inds = [0, 1]
        blocklist = [imlist[ind:ind + nimg] for ind in inds]
    else:
        blocklist = imlist

    for block, imgs in enumerate(blocklist):
        print(imgs)

        malt(imgs, ori, dirmec='{}_block{}'.format(dirmec, block), zoomf=zoomf)

        tawny('{}_block{}'.format(dirmec, block))

        mosaic_micmac_tiles('Orthophotomosaic', '{}_block{}'.format(dirmec, block))


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

    if os.name == 'nt':
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE, shell=True)
    else:
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
    if os.name == 'nt':
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE, shell=True)
    else:
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE)

    fn_gcp = fn_gcp + sub + '.xml'
    fn_meas = fn_meas + sub + '-S2D.xml'

    p = subprocess.Popen(['mm3d', 'Campari', img_pattern,
                          inori + sub,
                          outori + sub,
                          'GCP=[{},{},{},{}]'.format(os.path.join(outdir, fn_gcp),
                                                     np.abs(dx) / 4,  # should be correct within 0.25 pixel
                                                     os.path.join(outdir, fn_meas),
                                                     0.5),  # best balance for distortion
                          'SH={}'.format(homol),
                          'AllFree={}'.format(int(allfree))], stdin=echo.stdout)
    p.wait()

    out_gcps = get_campari_residuals('Ori-{}/Residus.xml'.format(outori + sub), in_gcps)
    # out_gcps.dropna(inplace=True)  # sometimes, campari can return no information for a gcp
    return out_gcps


def checkpoints(img_pattern, ori, fn_cp, fn_meas, fn_resids=None, ret_df=True):
    """
    Interface to run GCPCtrl to calculate checkpoint residuals for a given Orientation.

    :param str img_pattern: the match pattern for the images being input to Campari (e.g., "OIS.*tif")
    :param str ori: the full name of the orientation directory to use (e.g., Ori-TerrainFinal)
    :param str fn_cp: the filename of the CPs.xml file to use
    :param str fn_meas: the filename of the CP Measures.xml file to use
    :param str fn_resids: the (optional) filename to write the residuals for each checkpoint to
    :param bool ret_df: return a DataFrame with the residuals for each checkpoint
    :return:
    """
    args = ['mm3d', 'GCPCtrl', img_pattern, ori, fn_cp, fn_meas]
    if fn_resids is not None:
        args += [f'OutTxt={fn_resids}']

    p = subprocess.Popen(args)
    p.wait()

    if fn_resids is not None and ret_df:
        return pd.read_csv(str(fn_resids) + '_RollCtrl.txt', delimiter=r'\s+', names=['id', 'xres', 'yres', 'zres'])


def banana(fn_dem, fn_ref, deg=2, dZthresh=200., fn_mask=None, spacing=100):
    """
    Interface for running mm3d Postproc Banana, for computing a polynomial correction to a "banana" or dome.

    :param str fn_dem: the filename of the input DEM to correct
    :param str fn_ref: the filename of the reference DEM to use for the correction
    :param int deg: the degree of the polynomial correction (0 - 3, default: 2)
    :param float dZthresh: the threshold elevation difference between the reference and input DEMs (default: 200)
    :param str fn_mask: an (optional) exclusion mask to use for the reference DEM
    :param int spacing: the pixel spacing of the DEM to write (default: every 100 pixels)
    """
    assert deg in range(0, 4), "Polynomial degree limited to 0, 1, 2, or 3"

    dem_ext = os.path.splitext(fn_ref)[-1]
    fn_txt = fn_ref.replace(dem_ext, '.txt')
    dem_to_text(fn_ref, fn_out=fn_txt, spacing=spacing, fn_mask=fn_mask)

    # first, run gdal_translate to make a TFW file for the input dem
    p = subprocess.Popen(['gdal_translate', fn_dem, 'tmp.tif', '-co', 'TFW=YES'])
    p.wait()

    shutil.move('tmp.tfw', fn_dem.replace(os.path.splitext(fn_dem)[-1], '.tfw'))
    os.remove('tmp.tif')

    p = subprocess.Popen(['mm3d', 'Postproc', 'Banana',
                          fn_dem,
                          f'ListGCPs={fn_txt}',
                          f'DegPoly={deg}',
                          f'dZthresh={dZthresh}'])

    p.wait()


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

    bad_meas = np.abs(resids_df.errmax - resids_df.errmax.median()) > register.nmad(resids_df.errmax)
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

    gcps['camp_xy'] = np.sqrt(gcps.camp_xres ** 2 + gcps.camp_yres ** 2)

    while any([np.any(np.abs(gcps.camp_res - gcps.camp_res.median()) > 2 * register.nmad(gcps.camp_res)),
               np.any(np.abs(gcps.camp_dist - gcps.camp_dist.median()) > 2 * register.nmad(gcps.camp_dist)),
               gcps.camp_res.max() > 2]) and niter <= max_iter:
        valid_inds = np.logical_and.reduce((np.abs(gcps.camp_dist - gcps.camp_dist.median()) < 2 * register.nmad(gcps.camp_dist),
                                            gcps.camp_res < gcps.camp_res.max()))
        if np.count_nonzero(valid_inds) < 10:
            break

        gcps = gcps.loc[valid_inds]
        save_gcps(gcps, out_dir, register._get_utm_str(gcps.crs.to_epsg), subscript, fn_gcp=fn_gcp, fn_meas=fn_meas)

        gcps = bascule(gcps, out_dir, match_pattern, subscript, rel_ori, fn_gcp=fn_gcp,
                       fn_meas=fn_meas, outori=inori)
        gcps['res_dist'] = np.sqrt(gcps.xres ** 2 + gcps.yres ** 2)

        gcps = campari(gcps, out_dir, match_pattern, subscript, dx, ortho_res,
                       inori=inori, outori=outori, fn_gcp=fn_gcp, fn_meas=fn_meas,
                       allfree=allfree, homol=homol)

        gcps['camp_xy'] = np.sqrt(gcps.camp_xres ** 2 + gcps.camp_yres ** 2)
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
    zlist = glob('Z*.tif', root_dir=dir_mec)
    zlist.sort()

    etapes = [int(f.split('_')[1].replace('Num', '')) for f in zlist]

    ind = np.argmax(etapes)
    etape0 = max(etapes)

    fn_auto = os.path.join(dir_mec, 'AutoMask_STD-MALT_Num_{}.tif'.format(etape0 - 1))

    print(fn_auto)
    print(zlist[ind])

    dem = gu.Raster(os.path.join(dir_mec, zlist[ind]))

    automask = gu.Raster(os.path.join(dir_mec, 'AutoMask_STD-MALT_Num_{}.tif'.format(etape0 - 1)))

    shutil.copy(zlist[ind].replace('tif', 'tfw'),
                fn_auto.replace('tif', 'tfw'))

    # TODO: dem needs to have CRS set
    ref_dem = gu.Raster(fn_dem).reproject(dem)

    mask = gu.Vector(fn_mask).create_mask(dem).data

    dem.data[mask] = ref_dem.data[mask]
    automask.data[mask] = 1

    automask.save(fn_auto)
    dem.save(zlist[ind])

    if os.name == 'nt':
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE, shell=True)
    else:
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

    if os.name == 'nt':
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE, shell=True)
    else:
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


def dem_to_text(fn_dem, fn_out='dem_pts.txt', spacing=100, fn_mask=None):
    """
    Write elevations from a DEM raster to a text file for use in mm3d PostProc Banana.

    :param str fn_dem: the filename of the DEM to read.
    :param str fn_out: the name of the text file to write out (default: dem_pts.txt)
    :param int spacing: the pixel spacing of the DEM to write (default: every 100 pixels)
    """
    if isinstance(fn_dem, gu.Raster):
        dem = fn_dem
    else:
        dem = gu.Raster(fn_dem)

    if fn_mask is not None:
        mask = gu.Vector(fn_mask).create_mask(dem).data
        dem.data[mask] = np.nan

    x, y = dem.coords()

    z = dem.data[::spacing, ::spacing].flatten()
    x = x[::spacing, ::spacing].flatten()
    y = y[::spacing, ::spacing].flatten()

    x = x[np.logical_and(np.isfinite(z), ~z.mask)]
    y = y[np.logical_and(np.isfinite(z), ~z.mask)]
    z = z[np.logical_and(np.isfinite(z), ~z.mask)]

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


def _gdal_calc():
    if os.name == 'nt':
        # if we're on windows, call gdal_calc.py with the currently active python
        return ['python', os.path.join(sys.prefix, 'Scripts', 'gdal_calc.py')]
    else:
        # if we're not on windows, call gdal_calc.py as a shell script
        return ['gdal_calc.py']


def post_process(projstr, out_name, dirmec, do_ortho=True, ind_ortho=False):
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
    :param bool ind_ortho: apply a mask to each individual ortho image (default: false)
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

    subprocess.Popen(_gdal_calc() + ['--quiet', '-A', 'tmp_mask.tif', '-B', 'tmp_geo.tif',
                      '--outfile={}'.format(os.path.join('post_processed', f'{out_name}_Z.tif')),
                      '--calc="B*(A>0)"', '--NoDataValue=-9999']).wait()

    subprocess.Popen(['gdaldem', 'hillshade', os.path.join('post_processed', f'{out_name}_Z.tif'),
                      os.path.join('post_processed', f'{out_name}_HS.tif')]).wait()

    subprocess.Popen(_gdal_calc() + ['--quiet', '-A', 'tmp_corr.tif',
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
        subprocess.Popen(_gdal_calc() + ['--quiet', '-A', 'tmp_mask.tif', '-B', 'tmp_ortho.tif',
                          '--outfile={}'.format(os.path.join('post_processed', f'{out_name}_Ortho.tif')),
                          '--calc="B*(A>0)"', '--NoDataValue=0']).wait()
        os.remove('tmp_ortho.tif')

    os.remove('tmp_mask.tif')

    if ind_ortho:
        imlist = sorted(glob('OIS*.tif'))
        for fn_img in imlist:
            _mask_ortho(fn_img, out_name, dirmec, projstr)


def _mask_ortho(fn_img, out_name, dirmec, projstr):
    fn_ortho = os.path.join('-'.join(['Ortho', dirmec]), '_'.join(['Ort', fn_img]))
    fn_incid = os.path.join('-'.join(['Ortho', dirmec]), '_'.join(['Incid', fn_img]))
    fn_mask = os.path.join('-'.join(['Ortho', dirmec]), '_'.join(['Mask', fn_img]))

    shutil.copy(fn_ortho.replace('tif', 'tfw'), fn_mask.replace('tif', 'tfw'))

    mask = imread(fn_incid) < 1
    oy, ox = imread(fn_ortho).shape

    _mask = PIL.Image.fromarray(mask)
    mask = np.array(_mask.resize((ox, oy)))
    imsave(fn_mask, 255 * mask.astype(np.uint8))

    subprocess.Popen(['gdal_translate', '-a_srs', projstr, fn_ortho, 'tmp_ortho.tif']).wait()
    subprocess.Popen(['gdal_translate', '-a_srs', projstr, fn_mask, 'tmp_mask.tif']).wait()

    subprocess.Popen(_gdal_calc() + ['--quiet', '-A', 'tmp_ortho.tif', '-B', 'tmp_mask.tif',
                     '--outfile={}'.format(os.path.join('post_processed', f'{out_name}_{fn_img}')),
                     '--calc="A*(B>0)"', '--NoDataValue=0', '--type', 'Byte']).wait()

    os.remove('tmp_mask.tif')
    os.remove('tmp_ortho.tif')


# converted from bash script
def get_autogcp_locations(ori, meas_file, imlist):
    """
    Find location of automatically-detected control points in individual images using mm3d XYZ2Im.

    :param str ori: The orientation directory name (e.g., Ori-Relative)
    :param str meas_file: The Measures file to find image locations for
    :param list imlist: a list of image names
    """
    # nodist = '-'.join([ori, 'NoDist'])

    # copy the orientation directory to a new, "nodist" directory
    # shutil.copytree(ori, nodist, dirs_exist_ok=True)

    # autocals = glob('AutoCal*.xml', root_dir=nodist)
    # for autocal in autocals:
    #    _remove_distortion_coeffs(os.path.join(nodist, autocal))

    for im in imlist:
        # _update_autocal(nodist, im)

        # p = subprocess.Popen(['mm3d', 'XYZ2Im', os.path.join(nodist, f'Orientation-{im}.xml'),
        #                       meas_file, f'NoDist-{im}.txt'])
        # p.wait()

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


def init_git():
    """
    Initialize a git repository in the current working directory.
    """

    if shutil.which('git') is not None:
        # initialize an empty git repository with a main branch
        p = subprocess.Popen(['git', 'init'])
        p.wait()

        # copy the .gitignore file to the current directory
        _gitignore()

    else:
        # not sure if EnvironmentError is the best choice here
        raise EnvironmentError("unable to find git using shutil.which - please ensure that git is installed.")


def _gitignore():
    # section headers for the .gitignore file
    sect_headers = ['tar files', 'image files', 'point clouds', 'directories', 'xml files', 'log/txt files']

    # patterns to ignore for each section
    ignore = [['*.tar.gz', '*.tgz'],
              ['*.tif'],
              ['*.ply'],
              ['Homol*/', 'MEC-*/', 'Ortho-MEC-*/', 'Pastis/', 'Tmp-MM-Dir/'],
              ['SauvApero.xml', 'Tmp-SL-Glob.xml'],
              ['mm3d-LogFile.txt', 'WarnApero.txt']]

    ignore_dict = dict(zip(sect_headers, ignore))

    with open('.gitignore', 'w') as f:
        for sect in sect_headers:
            print(f'# {sect}', file=f)
            for ig in ignore_dict[sect]:
                print(ig, file=f)
            print('\n', file=f)


def meas_to_asp_gcp(fn_gcp, fn_meas, imlist, outname=None, scale=1):
    """
    Convert image measures stored in a micmac xml file to an ASP .gcp file format.

    :param str fn_gcp: the filename of the shapefile with the GCP coordinates
    :param str fn_meas: the filename of the xml file with the image measures
    :param list imlist: the image(s) to write point locations for
    :param str outname: the name of the output filename to create
    :param int scale: the factor by which to scale the image point locations
    """
    if outname is None:
        outname = fn_meas.replace('.xml', '.gcp')

    gcps = gpd.read_file(fn_gcp).to_crs(crs='epsg:4326').set_index('id')
    meas = parse_im_meas(fn_meas)

    meas = meas.loc[meas['image'].isin(imlist)]

    gcp_list = sorted(meas.name.unique())

    with open(outname, 'w') as f:
        for gcp in gcp_list:
            if all([gcp in meas.loc[meas.image == img]['name'].values for img in imlist]):
                _gcp = gcps.loc[gcp]
                lon, lat = _gcp.geometry.x, _gcp.geometry.y

                out_gcp = ','.join([gcp.strip('GCP'), str(lat), str(lon), str(_gcp.elevation), '1.0', '1.0', '1.0'])

                for img in sorted(imlist):
                    row, col = meas.loc[(meas.image == img) & (meas.name == gcp), ['i', 'j']].values[0]
                    out_gcp += ',' + ','.join([img, str(col / scale), str(row / scale), '1.0', '1.0'])

                print(out_gcp, file=f)
