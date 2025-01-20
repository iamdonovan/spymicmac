#!/usr/bin/env python
import argparse
from spymicmac.register import register_relative


def _argparser():
    helpstr = "Register a relative DEM or orthoimage to a reference DEM and/or orthorectified image."

    _parser = argparse.ArgumentParser(description=helpstr,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    _parser.add_argument('dirmec', action='store', type=str,
                         help='the name of the MEC directory to read the relative DEM from (e.g., MEC-Relative)')
    _parser.add_argument('fn_dem', action='store', type=str, help='path to reference DEM')
    _parser.add_argument('-ort', '--fn_ortho', action='store', type=str, default=None,
                         help='path to relative orthoimage (optional)')
    _parser.add_argument('-ref', '--fn_ref', action='store', type=str, default=None,
                         help='path to reference orthorectified image (optional)')
    _parser.add_argument('-glacmask', action='store', type=str, default=None,
                         help='path to shapefile of glacier outlines (i.e., an exclusion mask)')
    _parser.add_argument('-landmask', action='store', type=str, default=None,
                         help='path to shapefile of land outlines (i.e., an inclusion mask)')
    _parser.add_argument('-footprints', action='store', type=str, default=None,
                         help='path to shapefile of image outlines. If not set, will attempt to download from USGS.')
    _parser.add_argument('-im_subset', action='store', type=str, default=None, nargs='+',
                         help='subset of raw images to work with (default all)')
    _parser.add_argument('-b', '--block_num', action='store', type=str, default=None,
                         help='Block number to use if multiple image blocks exist in directory.')
    _parser.add_argument('--subscript', action='store', type=str, default=None,
                         help='Optional subscript to add to filenames.')
    _parser.add_argument('-ori', action='store', type=str, default='Relative',
                         help='name of orientation directory (after Ori-) [Relative]')
    _parser.add_argument('-ortho_res', action='store', type=float, default=8,
                         help='approx. ground sampling distance (pixel resolution) of ortho image. [8 m]')
    _parser.add_argument('-imgsource', action='store', type=str, default='DECLASSII',
                         help='USGS dataset name for images [DECLASSII]')
    _parser.add_argument('-density', action='store', type=int, default=200,
                         help='pixel spacing to look for GCPs [200]')
    _parser.add_argument('-no_allfree', action='store_false',
                         help='run Campari with AllFree set to False')
    _parser.add_argument('-useortho', action='store_true',
                         help='use the orthomosaic in Ortho-{dirmec} rather than the DEM [False]. '
                              'If fn_ortho is set, uses that file instead.')
    _parser.add_argument('-max_iter', action='store', type=int, default=5,
                         help='the maximum number of Campari iterations to run [5]')
    _parser.add_argument('-use_cps', action='store_true',
                         help='split the GCPs into GCPs and CPs, to quantify the uncertainty of the '
                              'camera model [False]')
    _parser.add_argument('-cp_frac', type=float, default=0.2,
                         help='the fraction of GCPs to use when splitting into GCPs and CPs [0.2]')
    _parser.add_argument('-o', '--use_orb', action='store_true',
                         help='use skimage.feature.ORB to identify GCP locations in the reference image '
                              '(default: use regular grid for matching)')
    _parser.add_argument('-fn_gcps', action='store', type=str, default=None,
                         help='(optional) shapefile or CSV of GCP coordinates to use. Column names should be '
                              '[(name | id), (z | elevation), x, y]. If CSV is used, x,y should have the same '
                              'CRS as the reference image.')
    return _parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    register_relative(args.dirmec, args.fn_dem,
                      fn_ortho=args.fn_ortho,
                      fn_ref=args.fn_ref,
                      glacmask=args.glacmask,
                      landmask=args.landmask,
                      footprints=args.footprints,
                      im_subset=args.im_subset,
                      block_num=args.block_num,
                      subscript=args.subscript,
                      ori=args.ori,
                      ortho_res=args.ortho_res,
                      imgsource=args.imgsource,
                      density=args.density,
                      allfree=args.no_allfree,
                      useortho=args.useortho,
                      max_iter=args.max_iter,
                      use_cps=args.use_cps,
                      cp_frac=args.cp_frac,
                      use_orb=args.use_orb,
                      fn_gcps=args.fn_gcps)


if __name__ == "__main__":
    main()
