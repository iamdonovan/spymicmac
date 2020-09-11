"""
sPyMicMac.micmac_tools is a collection of tools for interfacing with MicMac
"""
import os
import shutil
import numpy as np
import gdal
import pandas as pd
import lxml.etree as etree
import lxml.builder as builder
import difflib
import xml.etree.ElementTree as ET
from shapely.strtree import STRtree
from sPyMicMac.usgs_tools import get_usgs_footprints


######################################################################################################################
# MicMac interfaces - write xml files for MicMac to read
######################################################################################################################
def write_neighbour_images(imlist, fprints=None, nameField='ID', prefix='OIS-Reech_', fileExt='.tif', **kwargs):
    E = builder.ElementMaker()
    NamedRel = E.SauvegardeNamedRel()

    if fprints is None:
        if 'dataset' not in kwargs:
            dset = 'AERIAL_COMBIN'
        else:
            dset = kwargs['dataset']
        fprints = get_usgs_footprints(imlist, dataset=dset)
    else:
        fprints = fprints[fprints[nameField].isin(imlist)]

    s = STRtree([f for f in fprints['geometry'].values])

    for i, row in fprints.iterrows():
        fn = row[nameField]
        fp = row['geometry']

        res = s.query(fp)
        intersects = [c for c in res if fp.intersection(c).area > 0]
        fnames = [fprints[nameField][fprints['geometry'] == c].values[0] for c in intersects]
        fnames.remove(fn)

        for f in fnames:
            this_pair = E.Cple(' '.join([prefix + fn + fileExt,
                                         prefix + f + fileExt]))
            NamedRel.append(this_pair)

    tree = etree.ElementTree(NamedRel)
    tree.write('FileImagesNeighbour.xml', pretty_print=True, xml_declaration=True, encoding="utf-8")


def get_gcp_meas(im_name, meas_name, in_dir, E, nodist=None, gcp_name='GCP'):
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
    pt_els = []
    for ind, row in gcps.iterrows():
        this_mes = E.OneMesureAF1I(
                        E.NamePt(row['gcp']),
                        E.PtIm('{} {}'.format(row['im_col'], row['im_row']))
                        )
        pt_els.append(this_mes)
    return pt_els


def generate_measures_files(joined=False):
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


def write_auto_mesures(gcp_df, subscript, out_dir, outname='AutoMeasures'):
    with open(os.path.join(out_dir, '{}{}.txt'.format(outname, subscript)), 'w') as f:
        for i, row in gcp_df.iterrows():
            print('{} {} {}'.format(row.rel_x, row.rel_y, row.el_rel), file=f)


def write_auto_gcps(gcp_df, subscript, out_dir, utm_zone, outname='AutoGCPs'):
    with open(os.path.join(out_dir, '{}{}.txt'.format(outname, subscript)), 'w') as f:
        # print('#F= N X Y Z Ix Iy Iz', file=f)
        print('#F= N X Y Z', file=f)
        print('#Here the coordinates are in UTM {} X=Easting Y=Northing Z=Altitude'.format(utm_zone), file=f)
        for i, row in gcp_df.iterrows():
            # print('{} {} {} {} {} {} {}'.format(row.id, row.geometry.x, row.geometry.y, row.elevation,
            #                                        5/row.z_corr, 5/row.z_corr, 1), file=f)
            print('{} {} {} {}'.format(row.id, row.geometry.x, row.geometry.y, row.elevation), file=f)


def get_bascule_residuals(fn_basc, gcp_df):
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
    camp_root = ET.parse(fn_resids).getroot()
    camp_gcp_names = [a.find('Name').text for a in camp_root.findall('Iters')[-1].findall('OneAppui')]
    err_max = [float(a.find('EcartImMax').text) for a in camp_root.findall('Iters')[-1].findall('OneAppui')]
    camp_x = [float(a.find('EcartFaiscTerrain').text.split()[0])
              for a in camp_root.findall('Iters')[-1].findall('OneAppui')]
    camp_y = [float(a.find('EcartFaiscTerrain').text.split()[1])
              for a in camp_root.findall('Iters')[-1].findall('OneAppui')]

    for data_ in zip(camp_gcp_names, err_max):
        gcp_df.loc[gcp_df.id == data_[0], 'camp_res'] = data_[1]

    for data_ in zip(camp_gcp_names, camp_x):
        gcp_df.loc[gcp_df.id == data_[0], 'camp_xres'] = data_[1]

    for data_ in zip(camp_gcp_names, camp_y):
        gcp_df.loc[gcp_df.id == data_[0], 'camp_yres'] = data_[1]


    return gcp_df


def move_bad_tapas(ori):
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
