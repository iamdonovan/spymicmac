"""
spymicmac.micmac is a collection of tools for interfacing with MicMac commands, as well as for working with the various
XML and other files created/used by MicMac.
"""
import os
from pathlib import Path
import sys
import re
import subprocess
import shutil
import PIL
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
from skimage.transform import AffineTransform
import geoutils as gu
from . import data, register, resample, matching
from typing import Union
from numpy.typing import NDArray


######################################################################################################################
# MicMac interfaces - write xml files for MicMac to read
######################################################################################################################
def write_neighbour_images(imlist: Union[list, None] = None,
                           footprints: Union[str, Path, gpd.GeoDataFrame, None] = None,
                           name_field: str = 'ID', prefix: str = 'OIS-Reech_', file_ext: str = '.tif',
                           dataset: str = 'AERIAL_COMBIN', from_homol: bool = False,
                           img_pattern: str = 'OIS*.tif', dir_homol: str = 'Homol') -> None:
    """
    Write an xml file containing image pairs for processing with Tapioca, using either image footprints or a homologue
    directory to determine pairs of overlapping images.

    :param imlist: a list of (original) image names to use (e.g., without 'OIS-Reech\\_')
    :param footprints: an optional filename for a vector dataset of image footprints or a GeoDataFrame of image
        footprints. If None, uses spymicmac.data.get_usgs_footprints to download footprints based on imlist.
    :param name_field: the field in fprints table that contains the image name
    :param prefix: the prefix attached to the image name read by Tapioca
    :param file_ext: the file extension for the images read by Tapioca
    :param dataset: the USGS dataset name to search if no footprints are provided
    :param from_homol: get a list of pairs based on homologue files
    :param img_pattern: the image pattern to pass to glob to get a list of filenames
    :param dir_homol: the directory where the homologue files are
    """
    E = builder.ElementMaker()
    NamedRel = E.SauvegardeNamedRel()

    if from_homol:
        pairs = pairs_from_homol(img_pattern, dir_homol)
    else:
        if imlist is None:
            imlist = sorted([fn.strip(prefix).strip(file_ext) for fn in glob(img_pattern)])

        pairs = pairs_from_footprints(imlist=imlist, footprints=footprints, name_field=name_field, prefix=prefix,
                                      file_ext=file_ext, dataset=dataset)

    for pair in pairs:
        this_pair = E.Cple(' '.join(pair))
        NamedRel.append(this_pair)

    tree = etree.ElementTree(NamedRel)
    tree.write('FileImagesNeighbour.xml', pretty_print=True, xml_declaration=True, encoding="utf-8")


def pairs_from_footprints(imlist: list,
                          footprints: Union[str, Path, gpd.GeoDataFrame, None] = None,
                          name_field: str = 'ID', prefix: str = 'OIS-Reech_', file_ext: str = '.tif',
                          dataset='AERIAL_COMBIN') -> list:
    """
    Using a list of images and a collection of image footprints, return a list of potential image pairs for processing
    with Tapioca.

    :param imlist: a list of (original) image names to use (e.g., without 'OIS-Reech\\_')
    :param footprints: an optional filename for a vector dataset of image footprints or a GeoDataFrame of image
        footprints. If None, uses spymicmac.data.get_usgs_footprints to download footprints based on imlist.
    :param name_field: the field in fprints table that contains the image name
    :param prefix: the prefix attached to the image name read by Tapioca
    :param file_ext: the file extension for the images read by Tapioca
    :param dataset: the USGS dataset name to search if no footprints are provided
    :return: **pairs** -- a list of tuples representing image pairs
    """

    if footprints is None:
        footprints = data.get_usgs_footprints(imlist, dataset=dataset)
    elif isinstance(footprints, (str, Path)):
        footprints = gpd.read_file(footprints)
        footprints = footprints[footprints[name_field].isin(imlist)]
    else:
        footprints = footprints[footprints[name_field].isin(imlist)]

    footprints.reset_index(inplace=True)  # do this to ensure that strtree indices are correct
    s = STRtree([f for f in footprints['geometry'].values])

    all_pairs = []

    for ind, row in footprints.iterrows():
        fn = row[name_field]
        fp = row['geometry']

        print(fn)

        res = s.query(fp)
        intersects = [footprints.loc[c, 'geometry'] for c in res if fp.intersection(footprints.loc[c, 'geometry']).area > 0]
        fnames = [footprints[name_field][footprints['geometry'] == c].values[0] for c in intersects]
        try:
            fnames.remove(fn)
        except ValueError:
            pass

        all_pairs += [(prefix + fn + file_ext, prefix + fn_match + file_ext) for fn_match in fnames]

    return all_pairs


def _get_pairs(fn_img: str, dir_homol: Union[str, Path]) -> list:
    datlist = sorted(glob('*.dat', root_dir=Path(dir_homol, f"Pastis{fn_img}")))
    return sorted([fn.strip('.dat') for fn in datlist])


def pairs_from_homol(img_pattern: str = 'OIS*.tif', dir_homol: str = 'Homol') -> list:
    """
    Get a list of image pairs based on homologue files.

    :param img_pattern: the image pattern to pass to glob to get a list of filenames
    :param dir_homol: the directory where the homologue files are
    :return: **pairs** -- a list of tuples representing image pairs
    """
    imlist = sorted(glob(img_pattern))

    all_pairs = []

    for fn_img in imlist:
        pairs = _get_pairs(fn_img, dir_homol)
        all_pairs += zip(len(pairs) * [fn_img], pairs)

    return all_pairs


def _get_dat_sizes(fn_img: str, dir_homol: str) -> pd.DataFrame:
    datlist = sorted(glob('*.dat', root_dir=Path(dir_homol, f"Pastis{fn_img}")))
    sizes = np.array([os.stat(Path(dir_homol, f"Pastis{fn_img}", fn)).st_size for fn in datlist])
    percs = sizes / sizes.sum()

    size_df = pd.DataFrame({'image': fn_img, 'filename': datlist, 'size': sizes, 'percentage': percs})

    size_df['cumulative'] = size_df.sort_values('percentage', ascending=False)['percentage'].cumsum()

    return size_df


def clean_homol(img_pattern: str = 'OIS*.tif', dir_homol: str = 'Homol',
                min_size: Union[int, None] = None, remove_asymmetric: bool = False,
                return_df: bool = False) -> Union[None, pd.DataFrame]:
    """
    Remove spurious homologue files based on a threshold file size.

    :param img_pattern: the image pattern to pass to glob to get a list of filenames
    :param dir_homol: the directory where the homologue files are
    :param min_size: the size, in bytes, to use as a threshold for removing file (default: calculated from all files)
    :param remove_asymmetric: remove asymmetric homologue files (pairs where only one image in the pair "sees" the
        other one)
    :param return_df: return a DataFrame of all homologue files, rather than removing them
    :returns: a DataFrame of homologue files and associated information
    """
    imlist = sorted(glob(img_pattern))

    dat_sizes = pd.concat([_get_dat_sizes(fn, dir_homol) for fn in imlist], ignore_index=True)
    dat_sizes['symmetric'] = True

    for ind, row in dat_sizes.iterrows():
        dat_sizes.loc[ind, 'symmetric'] = row['image'] + '.dat' in dat_sizes.loc[
            dat_sizes['image'] == row['filename'].split('.dat')[0], 'filename'].to_list()

    if not return_df:
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

        return None
    else:
        return dat_sizes


def write_xml(fn_img: str, fn_mask: str = './MEC-Malt/Masq_STD-MALT_DeZoom1.tif',
              fn_xml: Union[None, str] = None, geomname: str = 'eGeomMNTEuclid') -> None:
    """
    Given a GDAL dataset, create a MicMac xml worldfile.

    :param fn_img: the filename of the image.
    :param fn_mask: the filename of the mask file
    :param fn_xml: the filename of the xml file to create (default: fn_img + '.xml')
    :param geomname: the MicMac Geometry name to use
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


def get_gcp_meas(im_name: str, meas_name: str, in_dir: str,
                 E: builder.ElementMaker,
                 nodist: Union[None, str] = None, gcp_name: str = 'GCP') -> builder.ElementMaker:
    """
    Create an lxml.builder.ElementMaker object with a GCP name and the image (row, pixel) location.

    :param im_name: the image name to write the GCP location for.
    :param meas_name: the name of the file to read the point locations from.
    :param in_dir: the name of the directory where the images and measures files are located.
    :param E: an ElementMaker object for writing to the xml file.
    :param nodist: the name of the directory
    :param gcp_name: the prefix (e.g., GCP0, GCP1, etc.) for the GCP name
    :return: **this_im_meas** -- an ElementMaker object with the GCP location in the image.
    """
    im = gdal.Open(Path(in_dir, im_name))
    maxj = im.RasterXSize
    maxi = im.RasterYSize

    impts = pd.read_csv(Path(in_dir, meas_name), sep=' ', names=['j', 'i'])
    if nodist is not None:
        impts_nodist = pd.read_csv(Path(in_dir, nodist), sep=' ', names=['j', 'i'])

    this_im_mes = E.MesureAppuiFlottant1Im(E.NameIm(im_name))
    for ind, row in impts.iterrows():
        in_im = 0 < row.j < maxj and 0 < row.i < maxi
        if nodist is not None:
            in_nd = -200 < impts_nodist.loc[ind, 'j'] + 200 < maxj and -200 < impts_nodist.loc[ind, 'i'] < maxi + 200
            in_im = in_im and in_nd
        if in_im:
            this_mes = E.OneMesureAF1I(
                                E.NamePt(f"{gcp_name}{ind+1}"),
                                E.PtIm(f"{row['j']} {row['i']}")
                            )
            this_im_mes.append(this_mes)
    return this_im_mes


def get_im_meas(points_df: pd.DataFrame, E: builder.ElementMaker,
                name: str = 'gcp', x: str = 'im_col', y: str = 'im_row') -> list:
    """
    Populate an lxml.builder.ElementMaker object with GCP image locations, for writing to xml files.

    :param points_df: a DataFrame with the points to find image locations for.
    :param E: an ElementMaker object for writing to the xml file.
    :param name: the column name in points_df corresponding to the point name
    :param x: the column name in points_df corresponding to the image x location
    :param y: the column name in points_df corresponding to the image y location
    :return: **pt_els** -- a list of ElementMaker objects corresponding to each GCP image location.
    """
    pt_els = []
    for row in points_df.itertuples():
        this_mes = E.OneMesureAF1I(
                        E.NamePt(getattr(row, name)),
                        E.PtIm(f"{getattr(row, x)} {getattr(row, y)}")
                        )
        pt_els.append(this_mes)
    return pt_els


def parse_im_meas(fn_meas: Union[str, Path]) -> pd.DataFrame:
    """
    Read an xml file with image locations into a pandas DataFrame.

    :param fn_meas: the name of the measures file to read.
    :return: **gcp_df** -- a DataFrame with gcp names and image locations.
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


def write_measures_im(meas_df: pd.DataFrame, fn_img: str) -> None:
    """
    Create a MeasuresIm xml file for an image.

    :param meas_df: a DataFrame of image measures, with [gcp, im_col, im_row] columns
    :param fn_img: the filename of the image.
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
    tree.write(Path('Ori-InterneScan', 'MeasuresIm-' + fn_img + '.xml'), pretty_print=True,
               xml_declaration=True, encoding="utf-8")


def generate_measures_files(joined: bool = True) -> None:
    """
    Create id_fiducial.txt, MeasuresCamera.xml, and Tmp-SL-Glob.xml files for KH-9 Hexagon mapping camera images.

    :param joined: generate files for joined scene (220x460 mm) instead of half (220x230mm)
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
            gcp_name = f"GCP_{row}_{col}"
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
            print(f"GCP_{row}_{col}", file=f)


def create_measurescamera_xml(fn_csv: Union[str, Path, pd.DataFrame],
                              ori: str = 'InterneScan', translate: bool = False,
                              name: str = 'gcp', x: str = 'im_col', y: str = 'im_row') -> None:
    """
    Create a MeasuresCamera.xml file from a csv of fiducial marker locations.

    :param fn_csv: the filename of the CSV file, or a pandas DataFrame.
    :param ori: the Ori directory to write the MeasuresCamera.xml file to. Defaults to (Ori-)InterneScan.
    :param translate: translate coordinates so that the origin is the upper left corner,
        rather than the principal point
    :param name: the column name in the csv file corresponding to the point name
    :param x: the column name in the csv file corresponding to the image x location
    :param y: the column name in the csv file corresponding to the image y location
    """
    assert type(fn_csv) in [pd.core.frame.DataFrame, str, Path], "fn_csv must be one of [str, DataFrame, Path]"
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
    tree.write(Path(f'Ori-{ori}', 'MeasuresCamera.xml'), pretty_print=True,
               xml_declaration=True, encoding="utf-8")


def estimate_measures_camera(approx: Union[pd.DataFrame, dict], pairs: list[tuple],
                             ori: str = 'InterneScan', scan_res: float = 2.5e-5,
                             how: str = 'mean', write_xml: bool = True) -> None:
    """
    Use a set of located fiducial markers to create a MeasuresCamera file using the average location of each fiducial
    marker.

    :param approx: a dataframe of approximate fiducial marker locations, or a dict of fiducial marker names and
        their angle with respect to the principal point. i.e., a mid-side marker on the right side of the frame should
        have an angle of 0, the fiducial marker in the upper right-hand corner should have an angle of 45° (pi / 4),
        a mid-side marker on the top of the frame should have an angle of 90° (pi / 2), and so on.
    :param pairs: a list of pairs of co-linear fiducial markers
    :param ori: The Ori- directory containing the MeasuresIm files
    :param scan_res: the scanning resolution of the images in m
    :param how: what average to use for the output locations. Must be one of [mean, median].
    :param write_xml: write the MeasuresCamera.xml file in addition to a CSV
    """
    assert how in ['mean', 'median'], "how must be one of [mean, median]"

    if isinstance(approx, dict):
        if max(approx.values()) > 2 * np.pi:
            angles = pd.Series(approx).apply(np.deg2rad)
        else:
            angles = pd.Series(approx)
    else:
        ppx, ppy = matching._meas_center(approx, pairs)
        approx['j'] -= ppx
        approx['i'] -= ppy

        angles = np.arctan2(approx['i'], approx['j'])
        angles[angles < 0] += 2 * np.pi

    meas_list = sorted(glob('MeasuresIm*.xml', root_dir=f'Ori-{ori}'))

    all_meas = []
    all_angles = []

    for fn_meas in meas_list:
        meas = parse_im_meas(Path(f'Ori-{ori}', fn_meas)).set_index('name')

        meas, rot = matching.tfm_measures(meas, pairs, angles)

        all_meas.append(meas.reset_index())
        all_angles.append(rot)

    all_meas = pd.concat(all_meas, ignore_index=True)

    all_meas['j'] /= all_meas['collim_dist']
    all_meas['i'] /= all_meas['collim_dist']

    if how == 'mean':
        avg_meas = all_meas.groupby('name').mean(numeric_only=True)
    else:
        avg_meas = all_meas.groupby('name').median(numeric_only=True)

    avg_meas['angle'] = np.arctan2(avg_meas['i'], avg_meas['j'])

    avg_meas['j'] *= avg_meas['collim_dist'] * scan_res * 1000 # convert from m to mm
    avg_meas['i'] *= avg_meas['collim_dist'] * scan_res * 1000 # convert from m to mm

    all_meas.to_csv('AllMeasures.csv', index=False)
    avg_meas.to_csv('AverageMeasures.csv')

    if write_xml:
        create_measurescamera_xml('AverageMeasures.csv', ori=ori, translate=True, name='name', x='j', y='i')


def generate_multicam_csv(patterns: Union[list, None] = None, prefix: str = 'OIS-Reech_',
                          fn_out: str = 'camera_defs.csv', name: Union[str, list] = '',
                          short_name: Union[str, list] = '', film_width: Union[str, list] = '',
                          film_height: Union[str, list] = '', focal: Union[str, float, list] = '') -> None:
    """
    Create a CSV file with camera parameters than can be read by create_localchantier_xml() to use images acquired by
    multiple cameras.

    Can be used to create a blank CSV template to be filled out manually, or generated using the optional function
    arguments.

    :param patterns: a list of filename patterns corresponding to each camera
    :param prefix: an optional prefix to add to the matching patterns
    :param fn_out: the name of the CSV file to create
    :param name: the name to give each camera. Must be unique.
    :param short_name: the "short name" description of each camera. Does not need to be unique.
    :param film_width: the width in mm of the frame for each camera.
    :param film_height: the height in mm of the frame for each camera.
    :param focal: the focal length of each camera, in mm.
    """
    cameras = pd.DataFrame()

    if patterns is None:
        cameras['pattern'] = ''
    else:
        patterns = [p + '.*' if '.*' not in p else p for p in patterns]
        patterns = [prefix + p if prefix not in p else p for p in patterns]

        cameras['pattern'] = patterns

    cameras['name'] = name
    cameras['short_name'] = short_name

    cameras['width'] = film_width
    cameras['height'] = film_height

    cameras['focal'] = focal

    cameras.to_csv(fn_out, index=False)


def create_localchantier_xml(name: str = 'KH9MC', short_name: str = 'KH-9 Hexagon Mapping Camera',
                             film_size: tuple[Union[int, float], Union[int, float]] = (460, 220),
                             pattern: str = '.*', focal: float = 304.8, add_sfs: bool = False,
                             cam_csv: Union[str, Path, None] = None) -> None:
    """
    Create a MicMac-LocalChantierDescripteur.xml file for a given camera. Default is the KH-9 Hexagon Mapping Camera.

    :param name: The name to use for the camera
    :param short_name: A short description of the camera
    :param film_size: the film size (width, height) in mm
    :param pattern: the matching pattern to use for the images
    :param focal: the nominal focal length, in mm
    :param add_sfs: use SFS to help find tie points in low-contrast images
    :param cam_csv: the CSV file containing parameters for multiple cameras
    """
    E = builder.ElementMaker()

    chantier = E.ChantierDescripteur()
    cam_db = E.LocCamDataBase()
    cam_assocs = E.KeyedNamesAssociations()
    foc_assocs = E.KeyedNamesAssociations()

    if cam_csv is not None:
        cameras = pd.read_csv(cam_csv)
        for ind, cam in cameras.iterrows():
            cam_db.append(
                E.CameraEntry(
                    E.Name(cam['name']),
                    E.SzCaptMm(f"{cam['width']} {cam['height']}"),
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


def init_autocal(imsize: tuple[int, int] = (32200, 15400),
                 framesize: tuple[Union[int, float], Union[int, float]] = (460, 220),
                 focal: float = 304.8, camname: str = 'KH9MC') -> None:
    """
    Create an AutoCal xml file for use in the Tapas step. Default values are for KH-9 Hexagon Mapping Camera.

    When calling mm3d Tapas, be sure to use "InCal=Init":

        mm3d Tapas RadialBasic "OIS.*tif" InCal=Init Out=Relative LibFoc=0

    The name of the file changes based on the focal length and camera name. Using the default values of
    foc=304.8 and camname='KH9MC' creates the following file in Ori-Init:

        AutoCal_Foc-KH9MC_304800.xml

    :param imsize: the size of the image (width, height) in pixels
    :param framesize: the size of the image (width, height) in mm
    :param focal: nominal focal length, in mm
    :param camname: the camera short name to use
    """
    os.makedirs('Ori-Init', exist_ok=True)

    scale = np.mean([imsize[0] / framesize[0], imsize[1] / framesize[1]])

    pp = np.array(imsize) / 2

    E = builder.ElementMaker()
    outxml = E.ExportAPERO(
        E.CalibrationInternConique(
            E.KnownConv('eConvApero_DistM2C'),
                E.PP(f"{pp[0]} {pp[1]}"),
                E.F(f"{scale * focal}"),
                E.SzIm(f"{imsize[0]} {imsize[1]}"),
                E.CalibDistortion(
                    E.ModRad(
                        E.CDist(f"{pp[0]} {pp[1]}"),
                        E.CoeffDist('2.74e-11'),
                        E.CoeffDist('-1.13e-21'),
                        E.CoeffDist('4.01e-29'),
                        E.CoeffDist('1.28e-38'),
                        E.CoeffDist('-4.32e-46'),
                        E.PPaEqPPs('true')
                    )
                )
        )
    )

    tree = etree.ElementTree(outxml)
    tree.write(Path('Ori-Init', f"AutoCal_Foc-{int(focal*1000)}_{camname}.xml"),
               pretty_print=True, xml_declaration=True, encoding="utf-8")


def load_cam_xml(fn_cam: Union[str, Path]) -> dict:
    """
    Parse an AutoCal xml file into a dictionary of intrinsic parameters.

    - pp: the principal point of the camera, in image coordinates
    - focal: the focal length of the camera
    - size: the size of the image (width height)
    - cdist: the image coordinates of the center of distortion
    - K*: the coefficients for the radial component of the distortion model, K1, K2, K3, etc.

    If the camera model used has both affine and decentric distortion parameters, these will also be included as
    [P1, P2] for the affine and [b1, b2] for the decentric distortion, respectively.

    Currently tested/working on outputs from Tapas Radial[Basic, Std, Extended] and Fraser[Basic].

    :param fn_cam: the name of the AutoCal xml file to parse.
    :return: **cam_dict** -- a dictionary of intrinsic parameter values for the camera
    """
    root = ET.parse(fn_cam).getroot()
    dist_model = root.find('CalibrationInternConique').find('CalibDistortion')

    cam_dict = {'pp': root.find('CalibrationInternConique').find('PP').text,
                'focal': root.find('CalibrationInternConique').find('F').text,
                'size': root.find('CalibrationInternConique').find('SzIm').text}

    if dist_model.find('ModRad') is not None:
        dist_model = dist_model.find('ModRad')
        cam_dict['cdist'] = dist_model.find('CDist').text
        for ind, coef in enumerate(dist_model.findall('CoeffDist')):
            cam_dict[f"K{ind+1}"] = coef.text

    elif dist_model.find('ModPhgrStd') is not None:
        dist_model = dist_model.find('ModPhgrStd')
        rad_part = dist_model.find('RadialePart')

        cam_dict['cdist'] = rad_part.find('CDist').text
        for ind, coef in enumerate(rad_part.findall('CoeffDist')):
            cam_dict[f"K{ind+1}"] = coef.text

        for param in ['P1', 'P2', 'b1', 'b2']:
            if dist_model.find(param) is not None:
                cam_dict[param] = dist_model.find(param).text
            else:
                cam_dict[param] = 0.0

    else:
        raise NotImplementedError("Camera model is not yet implemented...")

    return cam_dict


def write_cam_xml(fn_xml: Union[str, Path], cam_dict: dict, fraser: bool = True) -> None:
    """
    Write a camera xml file.

    :param fn_xml: the name of the xml file to write
    :param cam_dict: a dictionary containing camera parameters, as read by micmac.load_cam_xml()
    :param fraser: whether to add decentric and affine parameters to the xml
    """
    E = builder.ElementMaker()

    rad_coefs = [p for p in cam_dict.keys() if 'K' in p]

    if fraser:
        for param in ['P1', 'P2', 'b1', 'b2']:
            if param not in cam_dict.keys():
                cam_dict[param] = '0.0'

        rad_part = E.RadialePart(E.CDist(cam_dict['cdist']))
        for coef in rad_coefs:
            rad_part.append(E.CoeffDist(cam_dict[coef]))
        rad_part.append(E.PPaEqPPs('true'))

        dist_model = E.CalibDistortion(
            E.ModPhgrStd(
                rad_part,
                E.P1(cam_dict['P1']),
                E.P2(cam_dict['P2']),
                E.b1(cam_dict['b1']),
                E.b2(cam_dict['b2']),
            )
        )

    else:
        rad_part = E.ModRad(E.CDist(cam_dict['cdist']))
        for coef in rad_coefs:
            rad_part.append(E.CoeffDist(cam_dict[coef]))

        dist_model = E.CalibDistortion(
            rad_part
        )

    outxml = E.ExportAPERO(
        E.CalibrationInternConique(
            E.KnownConv('eConvApero_DistM2C'),
                E.PP(cam_dict['pp']),
                E.F(cam_dict['focal']),
                E.SzIm(cam_dict['size']),
                dist_model
        )
    )

    tree = etree.ElementTree(outxml)
    tree.write(fn_xml, pretty_print=True, xml_declaration=True, encoding="utf-8")


def get_match_pattern(imlist: list) -> str:
    """
    Given a list of image names, return a match pattern that can be passed to MicMac command line functions.

    :param imlist: a list of image names.
    :return: **pattern** -- a match pattern (e.g., "OIS.*tif") that can be passed to MicMac functions.
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


def write_auto_mesures(gcps: Union[pd.DataFrame, gpd.GeoDataFrame],
                       sub: str, outdir: Union[str, Path], outname: str = 'AutoMeasures') -> None:
    """
    Write a file with GCP locations in relaive space (x, y, z) to use with get_autogcp_locations.sh

    :param pandas.DataFrame gcps: a DataFrame with the GCPs to save.
    :param sub: the name of the block, if multiple blocks are being used (e.g., '_block1'). If not, use ''.
    :param outdir: the output directory to save the files to.
    :param outname: the base name of the file to create.
    """
    with open(Path(outdir, f"{outname}{sub}.txt"), 'w') as f:
        for row in gcps.itertuples():
            print(f"{row.rel_x} {row.rel_y} {row.el_rel}", file=f)


def get_valid_image_points(shape: tuple[int, int], pts: pd.DataFrame, pts_nodist: pd.DataFrame) -> NDArray:
    """
    Find which image points are located within an image based on the size of the image.

    :param shape: the shape of the image (rows, columns) to determine valid points for.
    :param pts: a DataFrame containing point locations (i, j)
    :param pts_nodist: a DataFrame containing point locations (i, j) calculated using no camera distortion.
    :return: **valid_pts** -- an array of the points that are located within the image shape.
    """
    maxi, maxj = shape

    in_im = np.logical_and.reduce((0 < pts.j, pts.j < maxj,
                                   0 < pts.i, pts.i < maxi))
    in_nd = np.logical_and.reduce((-200 < pts_nodist.j, pts_nodist.j < maxj + 200,
                                   -200 < pts_nodist.i, pts_nodist.i < maxi + 200))

    return np.logical_and(in_im, in_nd)


def write_image_mesures(imlist: list, gcps: Union[pd.DataFrame, gpd.GeoDataFrame],
                        outdir: str = '.', sub: str = '', ort_dir: str = 'Ortho-MEC-Relative',
                        outname: str = 'AutoMeasures') -> None:
    """
    Create a Measures-S2D.xml file (row, pixel) for each GCP in each image from a list of image names.

    :param imlist: a list of image names.
    :param gcps: a DataFrame of GCPs.
    :param outdir: the output directory to save the files to.
    :param sub: the name of the block, if multiple blocks are being used (e.g., '_block1').
    :param ort_dir: the Ortho-MEC directory where the images are located
    :param outname: the base name of the file to write
    """
    E = builder.ElementMaker()
    MesureSet = E.SetOfMesureAppuisFlottants()

    for im in imlist:
        print(im)
        ort_img = gu.Raster(Path(ort_dir, f"Ort_{im}"))
        dx, _, xmin, _, dy, ymin, _, _, _ = ort_img.transform
        ort_img = gu.Raster.from_array(resample.downsample(ort_img.data, fact=10),
                                       (10 * dx, 0, xmin, 0, 10 * dy, ymin), None)

        footprint = (ort_img > 0).polygonize().ds.union_all()
        valid = footprint.contains(gpd.points_from_xy(gcps.rel_x, gcps.rel_y))

        impts = pd.read_csv(f"Auto-{im}.txt", sep=' ', names=['j', 'i'])

        if np.count_nonzero(valid) == 0:
            continue

        this_im_mes = E.MesureAppuiFlottant1Im(E.NameIm(im))

        for ii, row in enumerate(impts[valid].itertuples()):
            this_mes = E.OneMesureAF1I(E.NamePt(gcps.iloc[row.Index]['id']),
                                       E.PtIm(f"{row.j} {row.i}"))
            this_im_mes.append(this_mes)

        MesureSet.append(this_im_mes)

    tree = etree.ElementTree(MesureSet)
    tree.write(Path(outdir, f"{outname}{sub}-S2D.xml"),
               pretty_print=True, xml_declaration=True, encoding="utf-8")


def write_auto_gcps(gcp_df: Union[pd.DataFrame, gpd.GeoDataFrame],
                    sub: str, outdir: str, utm_zone: str, outname: str = 'AutoGCPs') -> None:
    """
    Write GCP name, x, y, and z information to a text file to use with mm3d GCPConvert.

    :param pandas.DataFrame gcp_df: a DataFrame with the GCPs to save.
    :param sub: the name of the block, if multiple blocks are being used (e.g., '_block1'). If not, use ''.
    :param outdir: the output directory to save the files to.
    :param utm_zone: the UTM zone name (e.g., 8N).
    :param outname: the base name to use for the GCPs file
    """
    with open(Path(outdir, f"{outname}{sub}.txt"), 'w') as f:
        # print('#F= N X Y Z Ix Iy Iz', file=f)
        print('#F= N X Y Z', file=f)
        print(f"#Here the coordinates are in UTM {utm_zone} X=Easting Y=Northing Z=Altitude", file=f)
        for row in gcp_df.itertuples():
            print(f"{row.id} {row.geometry.x} {row.geometry.y} {row.elevation}", file=f)


def remove_measure(fn_meas: Union[str, Path], name: str) -> None:
    """
    Remove all instances of a given measure from an xml file.

    :param fn_meas: the xml file (e.g., AutoMeasures-S2D.xml)
    :param name: the measurement name (e.g., GCP0)
    """
    root = ET.parse(fn_meas).getroot()
    for im in root.findall('MesureAppuiFlottant1Im'):
        for mes in im.findall('OneMesureAF1I'):
            if mes.find('NamePt').text == name:
                im.remove(mes)

    tree = ET.ElementTree(root)
    tree.write(fn_meas, encoding="utf-8", xml_declaration=True)


def rename_gcps(root: ET.Element, ngcp: int = 0) -> tuple[dict, dict]:
    """
    Rename all GCPs in order of their appearance in an (opened) xml file.

    :param root: the root element of an xml tree
    :param ngcp: the number to start counting from

    :return:
        - **mes_dict** -- a dict containing image, gcp key/value pairs
        - **gcp_dict** -- a dict containing old/new gcp name key/value pairs
    """
    mes_dict = dict()
    gcp_dict = dict()

    for im in root.findall('MesureAppuiFlottant1Im'):
        this_name = im.find('NameIm').text
        these_mes = im.findall('OneMesureAF1I')

        for ii, mes in enumerate(these_mes):
            old_name = mes.find('NamePt').text

            if mes.find('NamePt').text not in gcp_dict.keys():
                gcp_dict[old_name] = f"GCP{ngcp}"
                ngcp += 1

            mes.find('NamePt').text = gcp_dict[old_name]
            these_mes[ii] = mes

        mes_dict[this_name] = these_mes

    return mes_dict, gcp_dict


def get_bascule_residuals(fn_basc: Union[str, Path], gcp_df: pd.DataFrame) -> pd.DataFrame:
    """
    Read a given GCPBascule residual file, and add the residuals to a DataFrame with GCP information.

    :param fn_basc: the GCPBascule xml file to read the residuals from.
    :param gcp_df: a DataFrame with the GCPs to read the residuals for.
    :return: **gcp_df** -- the input GCPs with the Bascule residuals added.
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


def get_campari_residuals(fn_resids: Union[str, Path], gcp_df: pd.DataFrame) -> pd.DataFrame:
    """
    Read a given Campari residual file, and add the residuals to a DataFrame with GCP information.

    :param fn_resids: the Campari residual xml file to read.
    :param gcp_df: a DataFrame with the GCPs to read the residuals for.
    :return: **gcp_df** -- the input GCPs with the Campari residuals added.
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


def get_tapas_residuals(ori: str) -> pd.DataFrame:
    """
    Read the image residuals output from Tapas.

    :param ori: the name of the Ori directory to read the residuals from (e.g., 'Relative' for Ori-Relative)
    :return: **img_df** -- a DataFrame with image names and residuals
    """
    root = ET.parse(Path(f'Ori-{ori}', 'Residus.xml'))
    last = root.findall('Iters')[-1]

    img_df = pd.DataFrame()
    for ind, img in enumerate(last.findall('OneIm')):
        img_df.loc[ind, 'name'] = img.find('Name').text
        img_df.loc[ind, 'res'] = float(img.find('Residual').text)
        img_df.loc[ind, 'perc_ok'] = float(img.find('PercOk').text)
        img_df.loc[ind, 'npts'] = int(img.find('NbPts').text)
        img_df.loc[ind, 'nmult'] = int(img.find('NbPtsMul').text)

    return img_df


def find_empty_homol(imlist: Union[list, None] = None, dir_homol: str = 'Homol',
                     pattern: str = 'OIS*.tif') -> None:
    """
    Search through a Homol directory to find images without any matches, then move them to a new directory called
    'EmptyMatch'

    :param imlist: a list of images in the current directory. If None, uses pattern to find images.
    :param dir_homol: the Homol directory to search in
    :param pattern: the search pattern to use to find images
    """

    if imlist is None:
        imlist = glob(pattern)

    empty = [fn_img for fn_img in imlist if len(_get_homol(fn_img, dir_homol)) == 0]

    os.makedirs('EmptyMatch', exist_ok=True)

    for fn_img in empty:
        print(f'{fn_img} -> EmptyMatch/{fn_img}')
        shutil.move(fn_img, 'EmptyMatch')


def _get_homol(fn_img: str, dir_homol: Union[str, Path] = 'Homol') -> list:
    if not os.path.exists(Path(dir_homol, 'Pastis' + fn_img)):
        return []
    else:
        return sorted([h.split('.dat')[0] for h in glob('*.dat', root_dir=Path(dir_homol, 'Pastis' + fn_img))])


# adapted from the fantastic answer provided by
# user matias-thayer and edited by user redbeam_
# at https://stackoverflow.com/a/50639220
def _get_connected_block(img: Union[list, str], seen: set, hdict: dict) -> tuple[list, set]:
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
def find_connected_blocks(pattern: str = 'OIS*.tif', dir_homol: Union[str, Path] = 'Homol') -> list:
    """
    Find connected blocks of images.

    :param pattern: the search pattern to use to get image names
    :param dir_homol: the Homologue directory to use to determine what images are connected
    :return: **blocks** -- a list containing lists of connected blocks of images
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


def separate_blocks(pattern: str = 'OIS*.tif', dir_homol: Union[str, Path] = 'Homol', min_size: int = 2) -> None:
    """
    Based on homologous points, find connected blocks of images and then separate the files into sub-folders.
    Moves files from {dir_homol} and Pastis, along with the image files.

    :param pattern: the search pattern to use to get image names
    :param dir_homol: the Homologue directory to use to determine what images are connected
    :param min_size: the minimum number of images to be considered a block
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


def move_bad_tapas(ori: str) -> None:
    """
    Read residual files output from Tapas (or Campari, GCPBascule), and move images with a NaN residual.

    :param ori: the orientation directory to read the residuals file from (e.g., 'Ori-Relative').
    """
    root = ET.parse(Path(ori, 'Residus.xml')).getroot()
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
        print(f"{im} -> bad/{im}")
        shutil.move(im, 'bad')


def _generate_glob(fn_ids: Union[str, Path]) -> None:
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


def batch_saisie_fids(imlist: list, flavor: str = 'qt', fn_cam: Union[None, str] = None,
                      clean: bool = True, gamma: Union[None, float] = None) -> None:
    """
    Run SaisieAppuisInit to locate the fiducial markers for a given list of images.

    :param imlist: the list of image filenames.
    :param flavor: which version of SaisieAppuisInit to run. Must be one of [qt, og]
    :param fn_cam: the filename for the MeasuresCamera.xml file
    :param clean: remove any image files in Tmp-SaisieAppuis
    :param gamma: Gamma adjustment value for Saisie
    """
    assert flavor in ['qt', 'og'], "flavor must be one of [qt, og]"

    os.makedirs('Ori-InterneScan', exist_ok=True)
    os.makedirs('Tmp-SaisieAppuis', exist_ok=True)

    if fn_cam is None:
        fn_cam = Path('Ori-InterneScan', 'MeasuresCamera.xml')

    if os.path.exists(fn_cam):
        measures_cam = parse_im_meas(fn_cam)
        with open('id_fiducial.txt', 'w') as f:
            for fid in measures_cam['name']:
                print(fid, file=f)
        _generate_glob('id_fiducial.txt')
    else:
        try:
            _generate_glob('id_fiducial.txt')
        except FileNotFoundError as e:
            raise FileNotFoundError('id_fiducial.txt not found. Please specify fn_cam, '
                                    'or ensure that id_fiducial.txt exists in the current directory.')

    if flavor == 'qt':
        saisie = 'SaisieAppuisInitQT'
    else:
        saisie = 'SaisieAppuisInit'

    for fn_img in imlist:
        if os.path.exists(Path('Ori-InterneScan', f'MeasuresIm-{fn_img}.xml')):
            if clean:
                tmplist = glob('*' + fn_img + '*', root_dir='Tmp-SaisieAppuis')
                for fn_tmp in tmplist:
                    os.remove(Path('Tmp-SaisieAppuis', fn_tmp))

            shutil.copy(Path('Ori-InterneScan', f'MeasuresIm-{fn_img}.xml'),
                        f'MeasuresIm-{fn_img}-S2D.xml')

            shutil.copy('Tmp-SL-Glob.xml',
                        Path('Tmp-SaisieAppuis', f'Tmp-SL-Glob-MeasuresIm-{fn_img}.xml'))

        saisie_args = ['mm3d', saisie, fn_img, 'NONE', 'id_fiducial.txt', f'MeasuresIm-{fn_img}.xml']

        if gamma is not None:
            saisie_args.append(f"Gama={gamma}")

        p = subprocess.Popen(saisie_args)
        p.wait()

        shutil.move(f'MeasuresIm-{fn_img}-S2D.xml', Path('Ori-InterneScan', f'MeasuresIm-{fn_img}.xml'))
        os.remove(f'MeasuresIm-{fn_img}-S3D.xml')

    os.remove('Tmp-SL-Glob.xml')


def tapioca(img_pattern: str = 'OIS.*tif', res_low: int = 400, res_high: int = 1200,
            fn_neighbours: Union[str, Path, None] = None) -> None:
    """
    Run mm3d Tapioca to find image tie points.

    :param img_pattern: The image pattern to pass to Tapioca (default: OIS.*tif)
    :param res_low: the size of the largest image axis, in pixels, for low-resolution matching
    :param res_high: the size of the largest image axis, in pixels, for high-resolution matching
    :param fn_neighbours: filename for an optional XML file containing image pairs
    """
    if os.name == 'nt':
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE, shell=True)
    else:
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE)

    if fn_neighbours is None:
        args = ['mm3d', 'Tapioca', 'MulScale', img_pattern, str(res_low), str(res_high)]
    else:
        args = ['mm3d', 'Tapioca', 'File', str(fn_neighbours), str(res_high)]

    p = subprocess.Popen(args, stdin=echo.stdout)

    return p.wait()


def martini(img_pattern: str = 'OIS.*tif', in_ori: Union[None, str] = None, ori_out: Union[None, str] = None,
            quick: bool = True) -> None:
    """
    Run mm3d Martini, which provides a quick way to orient images without solving for camera parameters.

    :param img_pattern: The image pattern to pass to Martini
    :param in_ori: the orientation directory to use to initialize the calibration
    :param ori_out: the name of the output orientation directory
    :param quick: run Martini in "quick" mode
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

    args.append(f"Quick={int(quick)}")

    p = subprocess.Popen(args, stdin=echo.stdout)

    return p.wait()


def tapas(cam_model: str, ori_out: Union[str, None] = None, img_pattern: str = 'OIS.*tif',
          in_cal: Union[str, None] = None, in_ori: Union[str, None] = None, lib_foc: bool = True,
          lib_pp: bool = True, lib_cd: bool = True) -> None:
    """
    Run mm3d Tapas with a given camera calibration model.

    Some basic camera calibration models for air photos:
        - RadialBasic
        - RadialStd
        - RadialExtended
        - FraserBasic
        - Fraser

    See MicMac docs for a full list/explanation of the camera models.

    :param cam_model: the camera calibration model to use.
    :param ori_out: the output orientation. Will create a directory, Ori-{ori_out}, with camera parameter files.
    :param img_pattern: the image pattern to pass to Tapas
    :param in_cal: an input calibration model to refine
    :param in_ori: a set of orientations to initialize the calibration
    :param lib_foc: allow the focal length to be calibrated
    :param lib_pp: allow the principal point to be calibrated
    :param lib_cd: allow the center of distortion to be calibrated
    """
    if os.name == 'nt':
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE, shell=True)
    else:
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE)

    args = ['mm3d', 'Tapas', cam_model, img_pattern,
            f"LibFoc={int(lib_foc)}", f"LibPP={int(lib_pp)}",
            f"LibCD={int(lib_cd)}"]

    if ori_out is not None:
        args.append('Out=' + ori_out)

    if in_cal is not None:
        args.append('InCal=' + in_cal)

    if in_ori is not None:
        args.append('InOri=' + in_ori)

    p = subprocess.Popen(args, stdin=echo.stdout)

    return p.wait()


def apericloud(ori: str, img_pattern: str = 'OIS.*tif',
               fn_out: Union[str, None] = None, with_points: bool = True) -> None:
    """
    Run mm3d AperiCloud to create a point cloud layer

    :param ori: the input orientation to use
    :param img_pattern: the image pattern to pass to AperiCloud
    :param fn_out: the output filename (default: AperiCloud_{ori}.ply)
    :param with_points: display the point cloud
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


def malt(imlist: Union[str, list], ori: str, zoomf: int = 1, zoomi: Union[None, int] = None,
         dirmec: str = 'MEC-Malt', seed_img: Union[str, Path, None] = None, seed_xml: Union[str, Path, None] = None,
         resol_terr: Union[float, int, None] = None, resol_ort: Union[float, int, None] = None,
         cost_trans: Union[float, int, None] = None, szw: Union[int, None] = None,
         regul: Union[float, None] = None, do_ortho: bool = True, do_mec: bool = True) -> None:
    """
    Run mm3d Malt Ortho.

    :param imlist: either a match pattern (e.g., OIS.*tif) or an iterable object of image filenames.
    :param ori: the orientation directory to use for Malt.
    :param zoomf: the final Zoom level to use
    :param zoomi: the initial Zoom level to use (default: not set)
    :param dirmec: the output MEC directory to create
    :param seed_img: a DEM to pass to Malt as DEMInitImg. Note that if seed_img is set, seed_xml
        must also be set. If used, it is recommended to set zoomi to be approximately equal to the DEM resolution -
        i.e., if the ortho resolution is 5 m and the seed DEM is 20 m, ZoomI should be 4.
    :param seed_xml: an XML file corresponding to the seed_img
    :param resol_terr: the resolution of the output DEM, in ground units (default: computed by mm3d)
    :param resol_ort: the resolution of the ortho images, relative to the output DEM - e.g., resol_ort=1 means
        the DEM and Orthoimage have the same resolution (default: 2.0)
    :param cost_trans: cost to change from correlation to decorrelation (default: 2.0)
    :param szw: the half-size of the correlation window to use - e.g., szw=1 means a 3x3 correlation window.
    :param regul: the regularization factor to use. Lower values mean higher potential variability between
        adjacent pixels, higher values (up to 1) mean smoother outputs
    :param do_ortho: whether to generate the orthoimages
    :param do_mec: whether to generate an output DEM
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


def tawny(dirmec: str, radiomegal: bool = False) -> None:
    """
    Run mm3d Tawny to create an orthomosaic.

    :param dirmec: the MEC directory to use
    :param radiomegal: run Tawny with RadiomEgal=1
    """
    if os.name == 'nt':
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE, shell=True)
    else:
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE)

    p = subprocess.Popen(['mm3d', 'Tawny', f"Ortho-{dirmec}", 'Out=Orthophotomosaic.tif',
                          f"RadiomEgal={int(radiomegal)}"], stdin=echo.stdout)
    return p.wait()


def block_malt(imlist: list, nimg: int = 3, ori: str = 'Relative', zoomf: int = 8) -> None:
    """
    Run mm3d Malt Ortho and mm3d Tawny on successive blocks of images.

    :param imlist: an iterable object of image filenames, or an iterable object of lists of image filenames
    :param nimg: the number of images to use in a block
    :param ori: the name of the orientation directory (e.g., Ori-Relative).
    :param zoomf: the final Zoom level to use
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

        malt(imgs, ori, dirmec=f"{dirmec}_block{block}", zoomf=zoomf)

        tawny(f"{dirmec}_block{block}")

        mosaic_micmac_tiles('Orthophotomosaic', f"{dirmec}_block{block}")


def bascule(in_gcps: pd.DataFrame, outdir: str, img_pattern: str, sub: str, ori: str,
            outori: str = 'TerrainRelAuto', fn_gcp: str = 'AutoGCPs', fn_meas: str ='AutoMeasures') -> pd.DataFrame:
    """
    Interface for running mm3d GCPBascule and reading the residuals from the resulting xml file.

    :param in_gcps: a DataFrame with the GCPs that are being input to Campari.
    :param outdir: the output directory where the AutoGCPs.xml file is saved.
    :param img_pattern: the match pattern for the images being input to Campari (e.g., "OIS.*tif")
    :param sub: the name of the block, if multiple blocks are being used (e.g., '_block1'). If not, use ''.
    :param ori: the name of the orientation directory (e.g., Ori-Relative).
    :param outori: the name of the output orientation directory.
    :param fn_gcp: the filename pattern for the GCP file. The file that will be loaded will be
        fn_gcp + sub + '.xml' (e.g., default: AutoGCPs -> AutoGCPs_block0.xml)
    :param fn_meas: the filename pattern for the measures file. The file that will be loaded will be
        fn_meas + sub + '-S2D.xml' (e.g., default: AutoMeasures -> AutoMeasures_block0-S2D.xml)
    :return: **out_gcps** -- the input gcps with the updated Bascule residuals.
    """
    fn_gcp = fn_gcp + sub + '.xml'
    fn_meas = fn_meas + sub + '-S2D.xml'

    if os.name == 'nt':
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE, shell=True)
    else:
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE)

    p = subprocess.Popen(['mm3d', 'GCPBascule', img_pattern, ori,
                          outori + sub,
                          Path(outdir, fn_gcp),
                          Path(outdir, fn_meas)], stdin=echo.stdout)
    p.wait()

    out_gcps = get_bascule_residuals(Path(f"Ori-{outori}{sub}",
                                                  'Result-GCP-Bascule.xml'), in_gcps)
    return out_gcps


def campari(in_gcps: pd.DataFrame, outdir: str, img_pattern: str, sub: str, dx: Union[int, float],
            ortho_res: Union[int, float], allfree: bool = True, fn_gcp: str = 'AutoGCPs',
            fn_meas: str = 'AutoMeasures', inori: str = 'TerrainRelAuto',
            outori: str = 'TerrainFinal', homol: str = 'Homol') -> pd.DataFrame:
    """
    Interface for running mm3d Campari and reading the residuals from the residual xml file.

    :param in_gcps: a DataFrame with the GCPs that are being input to Campari.
    :param outdir: the output directory where the AutoGCPs.xml file is saved.
    :param img_pattern: the match pattern for the images being input to Campari (e.g., "OIS.*tif")
    :param sub: the name of the block, if multiple blocks are being used (e.g., '_block1'). If not, use ''.
    :param dx: the pixel resolution of the reference image.
    :param ortho_res: the pixel resolution of the orthoimage being used.
    :param allfree: run Campari with AllFree=1 (True), meaning that all camera parameters will be optimized,
        or AllFree=0 (False), meaning that only the orientation will be optimized.
    :param fn_gcp: the filename pattern for the GCP file. The file that will be loaded will be
        fn_gcp + sub + '.xml'
    :param fn_meas: the filename pattern for the measures file. The file that will be loaded will be
        fn_meas + sub + '-S2D.xml'
    :param inori: the input orientation to Campari
    :param outori: the output orientation from Campari
    :param homol: the Homologue directory to use
    :return: **out_gcps** -- the input gcps with the updated Campari residuals.
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
                          f"GCP=[{Path(outdir, fn_gcp)},{np.abs(dx) / 4},{Path(outdir, fn_meas)},{0.5}]",
                          f"SH={homol}",
                          f"AllFree={int(allfree)}"], stdin=echo.stdout)
    p.wait()

    out_gcps = get_campari_residuals(Path(f"Ori-{outori+sub}", "Residus.xml"), in_gcps)
    return out_gcps


def checkpoints(img_pattern: str, ori: str, fn_cp: Union[str, Path], fn_meas: Union[str, Path],
                fn_resids: Union[str, Path, None] = None, ret_df: bool = True) -> Union[None, pd.DataFrame]:
    """
    Interface to run GCPCtrl to calculate checkpoint residuals for a given Orientation.

    :param str img_pattern: the match pattern for the images being input to Campari (e.g., "OIS.*tif")
    :param str ori: the full name of the orientation directory to use (e.g., Ori-TerrainFinal)
    :param str fn_cp: the filename of the CPs.xml file to use
    :param str fn_meas: the filename of the CP Measures.xml file to use
    :param str fn_resids: the (optional) filename to write the residuals for each checkpoint to
    :param bool ret_df: return a DataFrame with the residuals for each checkpoint
    :return: **cp_resids** -- a DataFrame with residuals for each checkpoint
    """
    args = ['mm3d', 'GCPCtrl', img_pattern, ori, fn_cp, fn_meas]
    if fn_resids is not None:
        args += [f'OutTxt={fn_resids}']

    p = subprocess.Popen(args)
    p.wait()

    if fn_resids is not None and ret_df:
        return pd.read_csv(str(fn_resids) + '_RollCtrl.txt', delimiter=r'\s+', names=['id', 'xres', 'yres', 'zres'])


def banana(fn_dem: Union[str, Path], fn_ref: Union[str, Path], deg: int = 2,
           dz_thresh: float = 200., fn_mask: Union[str, Path, None] = None, spacing: int = 100) -> None:
    """
    Interface for running mm3d Postproc Banana, for computing a polynomial correction to a "banana" or dome.

    :param fn_dem: the filename of the input DEM to correct
    :param fn_ref: the filename of the reference DEM to use for the correction
    :param deg: the degree of the polynomial correction (0 - 3)
    :param dz_thresh: the threshold elevation difference between the reference and input DEMs
    :param fn_mask: an (optional) exclusion mask to use for the reference DEM
    :param spacing: the sample spacing of the DEM to write, in pixels
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
                          f'dZthresh={dz_thresh}'])

    p.wait()


def remove_worst_mesures(fn_meas: Union[str, Path], ori: str) -> None:
    """
    Remove outlier measures from an xml file, given the output from Campari.

    :param fn_meas: the filename for the measures file.
    :param ori: the orientation directory output from Campari (e.g., Ori-TerrainFinal -> TerrainFinal)
    """
    camp_root = ET.parse(Path(f"Ori-{ori}", 'Residus.xml')).getroot()
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


def iterate_campari(gcps: pd.DataFrame, out_dir: str, match_pattern: str, subscript: str, dx: Union[int, float],
                    ortho_res: Union[int, float], fn_gcp: str = 'AutoGCPs', fn_meas: str = 'AutoMeasures',
                    rel_ori: str = 'Relative', inori: str = 'TerrainRelAuto', outori: str = 'TerrainFinal',
                    homol: str = 'Homol', allfree: bool = True, max_iter: int = 5) -> pd.DataFrame:
    """
    Run Campari iteratively, refining the orientation by removing outlier GCPs and Measures, based on their fit to the
    estimated camera model.

    :param gcps: a DataFrame with the GCPs that are being input to Campari.
    :param out_dir: the output directory where the GCP and Measures files are located.
    :param match_pattern: the match pattern for the images being input to Campari (e.g., "OIS.*tif")
    :param subscript: the name of the block, if multiple blocks are being used (e.g., '_block1'). If not, use ''.
    :param dx: the pixel resolution of the reference image.
    :param ortho_res: the pixel resolution of the orthoimage being used.
    :param fn_gcp: the filename pattern for the GCP file. The file that will be loaded will be
        fn_gcp + sub + '.xml' (e.g., default: AutoGCPs -> AutoGCPs_block0.xml)
    :param fn_meas: the filename pattern for the measures file. The file that will be loaded will be
        fn_meas + sub + '-S2D.xml' (e.g., default: AutoMeasures -> AutoMeasures_block0-S2D.xml)
    :param rel_ori: the name of the relative orientation to input to GCPBascule
    :param inori: the input orientation to Campari
    :param outori: the output orientation from Campari
    :param homol: the Homologue directory to use
    :param allfree: run Campari with AllFree=1 (True), meaning that all camera parameters will be optimized,
        or AllFree=0 (False), meaning that only the orientation will be optimized.
    :param max_iter: the maximum number of iterations to run.
    :return: **gcps** -- the gcps with updated residuals after the iterative process.
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


def mask_invalid_els(dir_mec: str, fn_dem: str, fn_mask: Union[str, Path], ori: str,
                     match_pattern: str = 'OIS.*tif', zoomf: int = 1) -> None:
    """
    Mask invalid elevations (e.g., water) in a DEM, then re-run the final step of mm3d Malt Ortho to make nicer
    orthophotos.

    :param dir_mec: the MEC directory (e.g., MEC-Malt) to use
    :param fn_dem: the filename of the reference DEM
    :param fn_mask: filename for the mask vector file
    :param ori: the orientation directory used to run Malt
    :param match_pattern: the match pattern used to
    :param zoomf: the final zoom level to run Malt at
    """
    zlist = glob('Z*.tif', root_dir=dir_mec)
    zlist.sort()

    etapes = [int(f.split('_')[1].replace('Num', '')) for f in zlist]

    ind = np.argmax(etapes)
    etape0 = max(etapes)

    fn_auto = Path(dir_mec, f"AutoMask_STD-MALT_Num_{etape0 - 1}.tif")

    print(fn_auto)
    print(zlist[ind])

    dem = gu.Raster(Path(dir_mec, zlist[ind]))

    automask = gu.Raster(Path(dir_mec, f"AutoMask_STD-MALT_Num_{etape0 - 1}.tif"))

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

    p = subprocess.Popen(['mm3d', 'Malt', 'Ortho', match_pattern, ori, f"DirMEC={dir_mec}",
                          'NbVI=2', 'MasqImGlob=filtre.tif', f"ZoomF={zoomf}",
                          'DefCor=1', 'CostTrans=1', 'EZA=1', 'DoMEC=0', 'DoOrtho=1', f"Etape0={etape0}"],
                         stdin=echo.stdout)
    p.wait()


def save_gcps(in_gcps: gpd.GeoDataFrame, outdir: str, utmstr: str, sub: str,
              fn_gcp: str = 'AutoGCPs', fn_meas: str = 'AutoMeasures') -> None:
    """
    Save a GeoDataFrame of GCP information to shapefile, txt, and xml formats.

    After running, the following new files will be created:

        - outdir/fn_gcp.shp (+ associated files)
        - outdir/fn_gcp.txt
        - outdir/fn_gcp.xml (output from mm3d GCPConvert)
        - outdir/fn_meas.xml (a file with image locations for each GCP)

    :param in_gcps: the gcps GeoDataFrame to save
    :param outdir: the output directory to save the files to
    :param utmstr: a UTM string generated by register.get_utm_str()
    :param sub: the name of the block, if multiple blocks are being used (e.g., '_block1'). If not, use ''.
    :param fn_gcp: the filename pattern for the GCP file. The file that will be loaded will be
        fn_gcp + sub + '.xml' (e.g., default: AutoGCPs -> AutoGCPs_block0.xml)
    :param fn_meas: the filename pattern for the measures file. The file that will be loaded will be
        fn_meas + sub + '-S2D.xml' (e.g., default: AutoMeasures -> AutoMeasures_block0-S2D.xml)
    """
    in_gcps.to_file(Path(outdir, fn_gcp + sub + '.shp'))
    write_auto_gcps(in_gcps, sub, outdir, utmstr, outname=fn_gcp)

    if os.name == 'nt':
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE, shell=True)
    else:
        echo = subprocess.Popen('echo', stdout=subprocess.PIPE)

    p = subprocess.Popen(['mm3d', 'GCPConvert', 'AppInFile',
                          Path(outdir, fn_gcp + sub + '.txt')], stdin=echo.stdout)
    p.wait()

    auto_root = ET.parse(Path(outdir, fn_meas + sub + '-S2D.xml')).getroot()
    for im in auto_root.findall('MesureAppuiFlottant1Im'):
        for pt in im.findall('OneMesureAF1I'):
            if pt.find('NamePt').text not in in_gcps.id.values:
                im.remove(pt)

    # save AutoMeasures
    out_xml = ET.ElementTree(auto_root)
    out_xml.write(Path(outdir, fn_meas + sub + '-S2D.xml'),
                  encoding="utf-8", xml_declaration=True)


def dem_to_text(fn_dem: Union[str, Path], fn_out: str = 'dem_pts.txt',
                spacing: int = 100, fn_mask: Union[str, Path, gpd.GeoDataFrame, None] = None) -> None:
    """
    Write elevations from a DEM raster to a text file for use in mm3d PostProc Banana.

    :param fn_dem: the filename of the DEM to read.
    :param fn_out: the name of the text file to write out (default: dem_pts.txt)
    :param spacing: the pixel spacing of the DEM to write (default: every 100 pixels)
    :param fn_mask: an optional filename or GeoDataFrame containing areas to mask out of the DEM.
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
def mosaic_micmac_tiles(filename: str, dirname: Union[str, Path] = '.') -> None:
    """
    Re-stitch images tiled by MicMac.

    :param filename: MicMac filename to mosaic together (e.g., Orthophotomosaic)
    :param dirname: directory containing images to mosaic
    """
    filelist = glob(f"{filename}_Tile*", root_dir=dirname)
    if len(filelist) == 0:
        print(f"No tiles found for {Path(dirname, filename)}; exiting.")
        return

    tiled = arrange_tiles(filelist, filename, dirname)
    I, J = tiled.shape

    arr_cols = []
    for j in range(J):
        arr_cols.append(np.concatenate(tiled[:, j], axis=0))

    img = np.concatenate(arr_cols, axis=1)

    imsave(Path(dirname, f"{filename}.tif"), img)


def arrange_tiles(flist: list, filename: str, dirname: Union[str, Path] = '.') -> NDArray:
    tmp_inds = [os.path.splitext(f)[0].split('Tile_')[-1].split('_') for f in flist]
    arr_inds = np.array([[int(a) for a in ind] for ind in tmp_inds])
    nrows = arr_inds[:, 1].max() + 1
    ncols = arr_inds[:, 0].max() + 1
    img_arr = np.array(np.zeros((nrows, ncols)), dtype='object')
    for i in range(nrows):
        for j in range(ncols):
            img_arr[i, j] = imread(Path(dirname, f"{filename}_Tile_{j}_{i}.tif"))
    return img_arr


def _gdal_calc() -> list:
    if os.name == 'nt':
        # if we're on windows, call gdal_calc.py with the currently active python
        return ['python', Path(sys.prefix, 'Scripts', 'gdal_calc.py')]
    else:
        # if we're not on windows, call gdal_calc.py as a shell script
        return ['gdal_calc.py']


def post_process(projstr: Union[str, int], out_name: str, dirmec: str,
                 do_ortho: bool = True, ind_ortho: bool = False) -> None:
    """
    Apply georeferencing and masking to the final DEM and Correlation images (optionally, the orthomosaic as well).

    Output files are written as follows:
        - DEM: post_processed/{out_name}_Z.tif
        - Hillshade: post_processed/{out_name}_HS.tif
        - Correlation: post_processed/{out_name}_CORR.tif
        - Orthomosaic: post_processed/{out_name}_Ortho.tif

    :param projstr: A string corresponding to the DEM's CRS that GDAL can use to georeference the rasters, or
        an int corresponding to the EPSG code for the DEM's CRS.
    :param out_name: The name that the output files should have.
    :param dirmec: The MEC directory to process files from (e.g., MEC-Malt)
    :param do_ortho: Post-process the orthomosaic in Ortho-{dirmec}, as well. Assumes that you have run
        mm3d Tawny with Out=Orthophotomosaic first.
    :param ind_ortho: apply a mask to each individual ortho image
    """

    os.makedirs('post_processed', exist_ok=True)

    # first, the stuff in MEC
    dem_list = sorted(glob('Z_Num*STD-MALT.tif', root_dir=dirmec))
    level = int(re.findall(r'\d+', dem_list[-1].split('_')[1])[0])
    zoomf = int(re.findall(r'\d+', dem_list[-1].split('_')[2])[0])

    shutil.copy(Path(dirmec, f'Z_Num{level}_DeZoom{zoomf}_STD-MALT.tfw'),
                Path(dirmec, f'Correl_STD-MALT_Num_{level-1}.tfw'))

    shutil.copy(Path(dirmec, f'Z_Num{level}_DeZoom{zoomf}_STD-MALT.tfw'),
                Path(dirmec, f'AutoMask_STD-MALT_Num_{level-1}.tfw'))

    if _needs_mosaic(Path(dirmec, f"Correl_STD-MALT_Num_{level-1}.tif")):
        mosaic_micmac_tiles(f"Correl_STD-MALT_Num_{level-1}", dirmec)

    if _needs_mosaic(Path(dirmec, f"Z_Num{level}_DeZoom{zoomf}_STD-MALT.tif")):
        mosaic_micmac_tiles(f"Z_Num{level}_DeZoom{zoomf}_STD-MALT", dirmec)

    if isinstance(projstr, int):
        projstr = f"EPSG:{projstr}"

    subprocess.Popen(['gdal_translate', '-a_nodata', '0', '-a_srs', projstr,
                      Path(dirmec, f'Correl_STD-MALT_Num_{level-1}.tif'),
                      'tmp_corr.tif']).wait()

    subprocess.Popen(['gdal_translate', '-a_srs', projstr,
                      Path(dirmec, f'Z_Num{level}_DeZoom{zoomf}_STD-MALT.tif'),
                      'tmp_geo.tif']).wait()

    subprocess.Popen(['gdal_translate', '-a_nodata', '0', '-a_srs', projstr,
                      Path(dirmec, f'AutoMask_STD-MALT_Num_{level-1}.tif'),
                      'tmp_mask.tif']).wait()

    subprocess.Popen(_gdal_calc() + ['--quiet', '-A', 'tmp_mask.tif', '-B', 'tmp_geo.tif',
                      f"--outfile={Path('post_processed', f"{out_name}_Z.tif")}",
                      '--calc="B*(A>0)"', '--NoDataValue=-9999']).wait()

    subprocess.Popen(['gdaldem', 'hillshade', Path('post_processed', f'{out_name}_Z.tif'),
                      Path('post_processed', f'{out_name}_HS.tif')]).wait()

    subprocess.Popen(_gdal_calc() + ['--quiet', '-A', 'tmp_corr.tif',
                      f"--outfile={Path('post_processed', f"{out_name}_CORR.tif")}",
                      '--calc="((A.astype(float)-127)/128)*100"', '--NoDataValue=-9999']).wait()

    # clean up the temporary files
    os.remove('tmp_geo.tif')
    os.remove('tmp_corr.tif')
    os.remove('tmp_mask.tif')

    # now, mask the ortho image(s)
    if do_ortho:
        fn_ortho = Path('Ortho-' + dirmec, 'Orthophotomosaic.tif')

        if _needs_mosaic(fn_ortho):
            mosaic_micmac_tiles('Orthophotomosaic', 'Ortho-' + dirmec)

        subprocess.Popen(['gdal_translate', '-a_nodata', '0', '-a_srs', projstr, fn_ortho,
                          Path('post_processed', f'{out_name}_Ortho.tif')]).wait()

    if ind_ortho:
        imlist = sorted(glob('OIS*.tif'))
        for fn_img in imlist:
            fn_ort = Path('Ortho-' + dirmec, f"Ort_{fn_img}")

            if _needs_mosaic(fn_ort):
                mosaic_micmac_tiles(f"Ort_{os.path.splitext(fn_img)[0]}", 'Ortho-' + dirmec)

            _mask_ortho(fn_img, out_name, dirmec, projstr)

def _needs_mosaic(fn_img: str) -> bool:
    fn_tile = os.path.splitext(fn_img)[0] + '_Tile_0_0.tif'
    if not os.path.exists(fn_tile):
        return False
    # if the tile exists, we compare the file size of the full image to the size of the first tile
    return os.path.getsize(fn_img) < os.path.getsize(fn_tile)


def _mask_ortho(fn_img: str, out_name: str, dirmec: str, projstr: Union[str, int]) -> None:
    fn_ortho = Path('-'.join(['Ortho', dirmec]), '_'.join(['Ort', fn_img]))
    fn_incid = Path('-'.join(['Ortho', dirmec]), '_'.join(['Incid', fn_img]))
    fn_mask = Path('-'.join(['Ortho', dirmec]), '_'.join(['Mask', fn_img]))

    shutil.copy(str(fn_ortho).replace('tif', 'tfw'), str(fn_mask).replace('tif', 'tfw'))

    mask = imread(fn_incid) < 1
    oy, ox = imread(fn_ortho).shape

    _mask = PIL.Image.fromarray(mask)
    mask = np.array(_mask.resize((ox, oy)))
    imsave(fn_mask, 255 * mask.astype(np.uint8))

    if isinstance(projstr, int):
        projstr = f"EPSG:{projstr}"

    subprocess.Popen(['gdal_translate', '-a_srs', projstr, fn_ortho, 'tmp_ortho.tif']).wait()
    subprocess.Popen(['gdal_translate', '-a_srs', projstr, fn_mask, 'tmp_mask.tif']).wait()

    subprocess.Popen(_gdal_calc() + ['--quiet', '-A', 'tmp_ortho.tif', '-B', 'tmp_mask.tif',
                     f"--outfile={Path('post_processed', f"{out_name}_{fn_img}")}",
                     '--calc="A*(B>0)"', '--NoDataValue=0', '--type', 'Byte']).wait()

    os.remove('tmp_mask.tif')
    os.remove('tmp_ortho.tif')


# converted from bash script
def get_autogcp_locations(ori: str, meas_file: Union[str, Path], imlist: list) -> None:
    """
    Find location of automatically-detected control points in individual images using mm3d XYZ2Im.

    :param ori: The orientation directory name (e.g., Ori-Relative)
    :param meas_file: The Measures file to find image locations for
    :param imlist: a list of image names
    """
    # nodist = '-'.join([ori, 'NoDist'])

    # copy the orientation directory to a new, "nodist" directory
    # shutil.copytree(ori, nodist, dirs_exist_ok=True)

    # autocals = glob('AutoCal*.xml', root_dir=nodist)
    # for autocal in autocals:
    #    _remove_distortion_coeffs(Path(nodist, autocal))

    for im in imlist:
        # _update_autocal(nodist, im)

        # p = subprocess.Popen(['mm3d', 'XYZ2Im', Path(nodist, f'Orientation-{im}.xml'),
        #                       meas_file, f'NoDist-{im}.txt'])
        # p.wait()

        p = subprocess.Popen(['mm3d', 'XYZ2Im', Path(ori, f'Orientation-{im}.xml'),
                              meas_file, f'Auto-{im}.txt'])
        p.wait()


def _remove_distortion_coeffs(fn_xml: Union[str, Path]) -> None:
    root = ET.parse(fn_xml).getroot()

    dist_coeffs = root.find('CalibrationInternConique').find('CalibDistortion').find('ModRad').findall('CoeffDist')
    inv_coeffs = root.find('CalibrationInternConique').find('CalibDistortion').find('ModRad').findall('CoeffDistInv')

    for coeff in dist_coeffs + inv_coeffs:
        coeff.text = '0.0'

    tree = ET.ElementTree(root)
    tree.write(fn_xml, encoding="utf-8", xml_declaration=True)


def _update_autocal(ori: str, im: str) -> None:
    fn_xml = Path(ori, f'Orientation-{im}.xml')
    root = ET.parse(fn_xml).getroot()

    old_autocal = root.find('OrientationConique').find('FileInterne').text
    old_autocal = Path(ori, os.path.basename(old_autocal))

    root.find('OrientationConique').find('FileInterne').text = old_autocal

    tree = ET.ElementTree(root)
    tree.write(fn_xml, encoding="utf-8", xml_declaration=True)


def init_git() -> None:
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


def _gitignore() -> None:
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
