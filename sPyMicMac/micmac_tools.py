"""
sPyMicMac.micmac_tools is a collection of tools for interfacing with MicMac
"""
from __future__ import print_function, division
import os
import numpy as np
import gdal
import pandas as pd
import lxml.etree as etree
import lxml.builder as builder


######################################################################################################################
# MicMac interfaces - write xml files for MicMac to read
######################################################################################################################
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


