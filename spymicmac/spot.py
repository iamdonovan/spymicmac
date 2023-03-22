import os
import xml.etree.ElementTree as ET
import lxml.etree as etree
import lxml.builder as builder
import numpy as np
import pandas as pd
from skimage.transform import AffineTransform
from shapely.geometry import Point, LineString
import pybob.landsat_tools as lt
from pybob.GeoImg import GeoImg


def rfm(X, Y, Z):

    if X.size == 1:
        return np.array([1, X, Y, Z,
                         X * Y, X * Z, Y * Z, X * X, Y * Y, Z * Z,
                         X * Y * Z, X * X * X, Y * Y * X, X * Z * Z,
                         X * X * Y, Y * Y * Y, Y * Z * Z, X * X * Z,
                         Y * Y * Z, Z * Z * Z]).flatten()
    else:
        return np.concatenate([np.ones(X.shape), X, Y, Z,
                               X * Y, X * Z, Y * Z, X * X, Y * Y, Z * Z,
                               X * Y * Z, X * X * X, Y * Y * X, X * Z * Z,
                               X * X * Y, Y * Y * Y, Y * Z * Z, X * X * Z,
                               Y * Y * Z, Z * Z * Z], axis=1)

# note: uses the same form that micmac uses:
# 1, X, Y, Z, Y*X, X*Z, Y*Z, X*X, Y*Y, Z*Z, X*Y*Z, X*X*X, Y*Y*X, X*Z*Z, X*X*Y, Y*Y*Y, Y*Z*Z, X*X*Z, Y*Y*Z, Z*Z*Z
def build_matrix(X, Y, Z, V):
    M_num = np.concatenate([np.ones(V.shape), X, Y, Z,
                            X * Y, X * Z, Y * Z, X * X, Y * Y, Z * Z,
                            X * Y * Z, X * X * X, Y * Y * X, X * Z * Z,
                            X * X * Y, Y * Y * Y, Y * Z * Z, X * X * Z,
                            Y * Y * Z, Z * Z * Z], axis=1)

    M_den = np.concatenate([X, Y, Z,
                            Y * X, X * Z, Y * Z, X * X, Y * Y, Z * Z,
                            X * Y * Z, X * X * X, Y * Y * X, X * Z * Z,
                            X * X * Y, Y * Y * Y, Y * Z * Z, X * X * Z,
                            Y * Y * Z, Z * Z * Z], axis=1)
    M = np.concatenate([M_num, -V * M_den], axis=1)

    return M


def weight_matrix(X, Y, Z, b):
    M = np.concatenate([np.ones(X.shape), X, Y, Z,
                        X * Y, X * Z, Y * Z, X * X, Y * Y, Z * Z,
                        X * Y * Z, X * X * X, Y * Y * X, X * Z * Z,
                        X * X * Y, Y * Y * Y, Y * Z * Z, X * X * Z,
                        Y * Y * Z, Z * Z * Z], axis=1)

    return np.diag(M.dot(b.T))


def build_deriv_matrix(X, Y, Z):
    # 1, X, Y, Z, XY, YZ, XX, YY, ZZ, XYZ, XXX, YYX, XZZ, XXY, YYY, YZZ, XXZ, YYZ, ZZZ
    # dx: 0, 1, 0, 0, Y, 0, 2X, 0, 0, YZ, 3XX, YY, ZZ, 2XY, 0, 0, 2XZ, 0, 0
    Mx = np.concatenate([np.zeros(X.shape), np.ones(X.shape), np.zeros(X.shape), np.zeros(X.shape),
                         Y, np.zeros(X.shape), 2 * X, np.zeros(X.shape), np.zeros(X.shape), Y * Z,
                         3 * X * X, Y * Y, Z * Z, 2 * X * Y, np.zeros(X.shape), np.zeros(X.shape),
                         2 * X * Z, np.zeros(X.shape), np.zeros(X.shape)], axis=1)

    # dy: 0, 0, 1, 0, X, Z, 0, 2Y, 0, XZ, 0, 2YX, 0, XX, 3YY, ZZ, 0, 2YZ, 0
    My = np.concatenate([np.zeros(X.shape), np.zeros(X.shape), np.ones(X.shape), np.zeros(X.shape),
                         X, Z, np.zeros(X.shape), 2 * Y, np.zeros(X.shape), X * Z, np.zeros(X.shape),
                         2 * Y * X, np.zeros(X.shape), X * X, 3 * Y * Y, Z * Z,
                         np.zeros(X.shape), 2 * Y * Z, np.zeros(X.shape)], axis=1)

    return np.concatenate([Mx, My], axis=1)


def find_min_max_el(pts, dem):
    ext = [pts.lon.min(), pts.lon.max(), pts.lat.min(), pts.lat.max()]
    cropped = dem.crop_to_extent(ext)
    return 100 * np.ceil(np.nanmax(cropped.img / 100)), 100 * np.floor(np.nanmin(cropped.img / 100))


def parse_meta(fn_meta):

    outdict = dict()

    meta_root = ET.parse(fn_meta).getroot()
    frames = meta_root.findall('Dataset_Frame')[0].findall('Vertex')

    lons = [float(vert.find('FRAME_LON').text) for vert in frames]
    lats = [float(vert.find('FRAME_LAT').text) for vert in frames]
    rows = [float(vert.find('FRAME_ROW').text) for vert in frames]
    cols = [float(vert.find('FRAME_COL').text) for vert in frames]
    init_pts = pd.DataFrame()

    init_pts['lon'] = lons
    init_pts['lat'] = lats
    init_pts['hgt'] = 0

    init_pts['row'] = rows
    init_pts['col'] = cols

    outdict['corners'] = init_pts

    outdict['center'] = meta_root.findall('Dataset_Frame')[0].find('Scene_Center')

    lcs = meta_root.find('Geoposition').find('Simplified_Location_Model').find('Direct_Location_Model')\
        .find('lc_List').findall('lc')
    outdict['line_coeffs'] = np.array([float(lc.text) for lc in lcs])

    pcs = meta_root.find('Geoposition').find('Simplified_Location_Model')\
        .find('Direct_Location_Model').find('pc_List').findall('pc')
    outdict['samp_coeffs'] = np.array([float(pc.text) for pc in pcs])

    return outdict


def rpc_grid(metadict, fn_dem=None, npts=20):
    corners = metadict['corners']

    _rows = np.linspace(corners.row.min(), corners.row.max(), npts)
    _cols = np.linspace(corners.col.min(), corners.col.max(), npts)

    CC, RR = np.meshgrid(_cols, _rows)
    grid = np.concatenate((CC.reshape(-1, 1), RR.reshape(-1, 1)), axis=1)

    line = grid[:, 1]
    samp = grid[:, 0]

    X = np.array([np.ones(line.shape), line, samp, samp * line, line * line, samp * samp])
    grid_lat = (metadict['samp_coeffs'].reshape(-1, 1) * X).sum(axis=0)
    grid_lon = (metadict['line_coeffs'].reshape(-1, 1) * X).sum(axis=0)

    out_pts = pd.DataFrame()

    if fn_dem is not None: # if we have a dem, use the min/max elevation from that
        dem = GeoImg(fn_dem)
        min_el, max_el = find_min_max_el(corners, dem)
    else: # default to min/max elevation for the whole earth
        min_el = -100
        max_el = 8800

    els = np.linspace(min_el, max_el + 1, npts)

    for el in els:
        these_pts = pd.DataFrame()
        these_pts['lat'] = grid_lat
        these_pts['lon'] = grid_lon
        these_pts['row'] = line
        these_pts['col'] = samp
        these_pts['hgt'] = el

        out_pts = pd.concat([out_pts, these_pts], ignore_index=True)

    return out_pts


def get_rescaling_factors(gridpts, metadict):

    valid_dict = dict()
    valid_dict['min_row'] = metadict['corners'].row.min()
    valid_dict['min_col'] = metadict['corners'].col.min()

    valid_dict['max_row'] = metadict['corners'].row.max()
    valid_dict['max_col'] = metadict['corners'].col.max()

    valid_dict['row_off'] = np.mean([valid_dict['max_row'], valid_dict['min_row']])
    valid_dict['row_scale'] = np.mean([valid_dict['max_row'], valid_dict['min_row']])

    valid_dict['col_off'] = np.mean([valid_dict['max_col'], valid_dict['min_col']])
    valid_dict['col_scale'] = np.mean([valid_dict['max_col'], valid_dict['min_col']])

    ul = metadict['corners'].loc[np.logical_and(metadict['corners'].row == valid_dict['min_row'],
                                                metadict['corners'].col == valid_dict['min_col'])]
    ur = metadict['corners'].loc[np.logical_and(metadict['corners'].row == valid_dict['min_row'],
                                                metadict['corners'].col == valid_dict['max_col'])]

    ll = metadict['corners'].loc[np.logical_and(metadict['corners'].row == valid_dict['max_row'],
                                                metadict['corners'].col == valid_dict['min_col'])]
    lr = metadict['corners'].loc[np.logical_and(metadict['corners'].row == valid_dict['max_row'],
                                                metadict['corners'].col == valid_dict['max_col'])]

    valid_dict['last_lon'] = np.mean([ur.lon, lr.lon])
    valid_dict['first_lon'] = np.mean([ul.lon, ll.lon])

    valid_dict['lon_off'] = gridpts.lon.mean()
    valid_dict['lon_scale'] = max(np.abs(gridpts.lon.max() - valid_dict['lon_off']),
                                  np.abs(valid_dict['lon_off'] - gridpts.lon.min()))

    valid_dict['first_lat'] = np.mean([ul.lat, ur.lat])
    valid_dict['last_lat'] = np.mean([ll.lat, lr.lat])

    valid_dict['lat_off'] = gridpts.lat.mean()
    valid_dict['lat_scale'] = max(np.abs(gridpts.lat.max() - valid_dict['lat_off']),
                                  np.abs(valid_dict['lat_off'] - gridpts.lat.min()))

    valid_dict['hgt_off'] = gridpts.hgt.mean()
    valid_dict['hgt_scale'] = max(np.abs(gridpts.hgt.max() - valid_dict['hgt_off']),
                                  np.abs(valid_dict['hgt_off'] - gridpts.hgt.min()))

    return valid_dict


def rescale_grid(gridpts, metadict):

    rescaled = gridpts.copy()

    valid_dict = get_rescaling_factors(gridpts, metadict)

    rescaled['row'] = (gridpts['row'] - valid_dict['row_off']) / valid_dict['row_scale']
    rescaled['col'] = (gridpts['col'] - valid_dict['col_off']) / valid_dict['col_scale']

    rescaled['lat'] = (gridpts['lat'] - valid_dict['lat_off']) / valid_dict['lat_scale']
    rescaled['lon'] = (gridpts['lon'] - valid_dict['lon_off']) / valid_dict['lon_scale']

    rescaled['hgt'] = (gridpts['hgt'] - valid_dict['hgt_off']) / valid_dict['hgt_scale']

    return rescaled


def solve_rpcs(rescaled):

    resdict = dict()

    resdict['inv_samp'] = _rpc_solver(rescaled['col'].values.reshape(-1, 1),
                                      rescaled['row'].values.reshape(-1, 1),
                                      rescaled['hgt'].values.reshape(-1, 1),
                                      rescaled['lon'].values.reshape(-1, 1))

    resdict['inv_line'] = _rpc_solver(rescaled['col'].values.reshape(-1, 1),
                                      rescaled['row'].values.reshape(-1, 1),
                                      rescaled['hgt'].values.reshape(-1, 1),
                                      rescaled['lat'].values.reshape(-1, 1))

    resdict['samp'] = _rpc_solver(rescaled['lon'].values.reshape(-1, 1),
                                  rescaled['lat'].values.reshape(-1, 1),
                                  rescaled['hgt'].values.reshape(-1, 1),
                                  rescaled['col'].values.reshape(-1, 1))

    resdict['line'] = _rpc_solver(rescaled['lon'].values.reshape(-1, 1),
                                  rescaled['lat'].values.reshape(-1, 1),
                                  rescaled['hgt'].values.reshape(-1, 1),
                                  rescaled['row'].values.reshape(-1, 1))

    return resdict


def _rpc_solver(X, Y, Z, V):
    Mat = build_matrix(X, Y, Z, V)

    init_soln = np.linalg.solve(Mat.T.dot(Mat), Mat.T.dot(V))

    for ii in range(10):
        B = np.concatenate([[1], init_soln[20:].reshape(-1)])
        W = weight_matrix(X, Y, Z, B)

        this_soln = np.linalg.solve(Mat.T.dot(W*W).dot(Mat), Mat.T.dot(W*W).dot(V))
        init_soln = this_soln

    return this_soln


def _split_result(result):
    num = result[:20].flatten()
    den = np.concatenate([np.array([[1]]), result[20:]]).flatten()

    return num, den


def write_rpc_xml(outname, valid_dict, resdict):

    samp_num, samp_den = _split_result(resdict['samp'])
    line_num, line_den = _split_result(resdict['line'])

    lat_num, lat_den = _split_result(resdict['inv_samp'])
    lon_num, lon_den = _split_result(resdict['inv_line'])

    E = builder.ElementMaker()

    SampNumCoeff = E.SAMP_NUM_COEFF()
    SampDenCoeff = E.SAMP_DEN_COEFF()

    InvSampNumCoeff = E.SAMP_NUM_COEFF()
    InvSampDenCoeff = E.SAMP_DEN_COEFF()

    LineNumCoeff = E.LINE_NUM_COEFF()
    LineDenCoeff = E.LINE_DEN_COEFF()

    InvLineNumCoeff = E.LINE_NUM_COEFF()
    InvLineDenCoeff = E.LINE_DEN_COEFF()

    for ii in range(1, 21):
        SampNumCoeff.append(E('SAMP_NUM_COEFF_{}'.format(ii), '{}'.format(samp_num[ii-1])))
        SampDenCoeff.append(E('SAMP_DEN_COEFF_{}'.format(ii), '{}'.format(samp_den[ii-1])))

        InvSampNumCoeff.append(E('SAMP_NUM_COEFF_{}'.format(ii), '{}'.format(lon_num[ii-1])))
        InvSampDenCoeff.append(E('SAMP_DEN_COEFF_{}'.format(ii), '{}'.format(lon_den[ii-1])))

        LineNumCoeff.append(E('LINE_NUM_COEFF_{}'.format(ii), '{}'.format(line_num[ii-1])))
        LineDenCoeff.append(E('LINE_DEN_COEFF_{}'.format(ii), '{}'.format(line_den[ii-1])))

        InvLineNumCoeff.append(E('LINE_NUM_COEFF_{}'.format(ii), '{}'.format(lat_num[ii-1])))
        InvLineDenCoeff.append(E('LINE_DEN_COEFF_{}'.format(ii), '{}'.format(lat_den[ii-1])))

    xml_rpc = E.Xml_RPC(
        E.METADATA_FORMAT('DIMAP'),
        E.METADATA_VERSION('2.0'),

        E.Direct_Model(
            SampNumCoeff,
            SampDenCoeff,
            LineNumCoeff,
            LineDenCoeff
        ),

        E.Inverse_Model(
            InvSampNumCoeff,
            InvSampDenCoeff,
            InvLineNumCoeff,
            InvLineDenCoeff
        ),

        E.RFM_Validity(
            E.FIRST_ROW('{}'.format(valid_dict['min_row'])),
            E.FIRST_COL('{}'.format(valid_dict['min_col'])),
            E.LAST_ROW('{}'.format(valid_dict['max_row'])),
            E.LAST_COL('{}'.format(valid_dict['max_col'])),
            E.FIRST_LON('{}'.format(valid_dict['first_lon'])),
            E.FIRST_LAT('{}'.format(valid_dict['first_lat'])),
            E.LAST_LON('{}'.format(valid_dict['last_lon'])),
            E.LAST_LAT('{}'.format(valid_dict['last_lat'])),
            E.LONG_SCALE('{}'.format(valid_dict['lon_scale'])),
            E.LONG_OFF('{}'.format(valid_dict['lon_off'])),
            E.LAT_SCALE('{}'.format(valid_dict['lat_scale'])),
            E.LAT_OFF('{}'.format(valid_dict['lat_off'])),
            E.HEIGHT_SCALE('{}'.format(valid_dict['hgt_scale'])),
            E.HEIGHT_OFF('{}'.format(valid_dict['hgt_off'])),
            E.SAMP_SCALE('{}'.format(valid_dict['col_scale'])),
            E.SAMP_OFF('{}'.format(valid_dict['col_off'])),
            E.LINE_SCALE('{}'.format(valid_dict['row_scale'])),
            E.LINE_OFF('{}'.format(valid_dict['row_off']))
        )
    )

    outxml = E.Global(xml_rpc)
    tree = etree.ElementTree(outxml)
    tree.write(outname, pretty_print=True, xml_declaration=True, encoding="utf-8")


def read_asp(fn_asp):

    meta_root = ET.parse(fn_asp).getroot()
    rpc_meta = meta_root.find('RPB').find('IMAGE')

    meta = dict()
    meta['min_row'] = int(1)
    meta['min_col'] = int(1)
    meta['max_row'] = int(12000)
    meta['max_col'] = int(12000)
    meta['first_lon'] = float(1)
    meta['first_lat'] = float(1)
    meta['last_lon'] = float(1)
    meta['last_lat'] = float(1)
    meta['lon_scale'] = float(rpc_meta.find('LONGSCALE').text)
    meta['lon_off'] = float(rpc_meta.find('LONGOFFSET').text)
    meta['lat_scale'] = float(rpc_meta.find('LATSCALE').text)
    meta['lat_off'] = float(rpc_meta.find('LATOFFSET').text)
    meta['hgt_scale'] = float(rpc_meta.find('HEIGHTSCALE').text)
    meta['hgt_off'] = float(rpc_meta.find('HEIGHTOFFSET').text)
    meta['col_scale'] = float(rpc_meta.find('SAMPSCALE').text)
    meta['col_off'] = float(rpc_meta.find('SAMPOFFSET').text)
    meta['row_scale'] = float(rpc_meta.find('LINESCALE').text)
    meta['row_off'] = float(rpc_meta.find('LINEOFFSET').text)

    linenum = np.array([float(p) for p in rpc_meta.find('LINENUMCOEFList').find('LINENUMCOEF').text.split()])
    lineden = np.array([float(p) for p in rpc_meta.find('LINEDENCOEFList').find('LINEDENCOEF').text.split()])

    meta['line'] = np.concatenate([linenum, lineden[1:]]).reshape(-1, 1)

    sampnum = np.array([float(p) for p in rpc_meta.find('SAMPNUMCOEFList').find('SAMPNUMCOEF').text.split()])
    sampden = np.array([float(p) for p in rpc_meta.find('SAMPDENCOEFList').find('SAMPDENCOEF').text.split()])

    meta['samp'] = np.concatenate([sampnum, sampden[1:]]).reshape(-1, 1)

    meta['inv_line'] = np.ones(meta['samp'].shape)
    meta['inv_samp'] = np.ones(meta['samp'].shape)

    return meta


def asp2micmac(fn_asp, fn_micmac):
    """

    :param fn_asp:
    :param fn_micmac:
    :return:
    """
    asp_meta = read_asp(fn_asp)
    write_rpc_xml(fn_micmac, asp_meta, asp_meta)
