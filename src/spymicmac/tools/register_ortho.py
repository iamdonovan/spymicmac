#!/usr/bin/env python
import argparse
from spymicmac.register import register_ortho


def _argparser():
    _parser = argparse.ArgumentParser(description="Register a relative orthoimage and DEM to a reference orthorectified image and DEM.",
                                      formatter_class=argparse.RawDescriptionHelpFormatter)

    _parser.add_argument('fn_ortho', action='store', type=str, help='non-referenced orthophoto mosaic')
    _parser.add_argument('fn_ref', action='store', type=str, help='georeferenced satellite image')
    _parser.add_argument('fn_reldem', action='store', type=str, help='relative dem corresponding to ortho')
    _parser.add_argument('fn_dem', action='store', type=str, help='dem')
    _parser.add_argument('-glacmask', action='store', type=str, default=None,
                         help='path to shapefile of glacier outlines (i.e., an exclusion mask)')
    _parser.add_argument('-landmask', action='store', type=str, default=None,
                         help='path to shapefile of land outlines (i.e., an inclusion mask)')
    _parser.add_argument('-footprints', action='store', type=str, default=None,
                         help='path to shapefile of image outlines. If not set, will download from USGS.')
    _parser.add_argument('-im_subset', action='store', type=str, default=None, nargs='+',
                         help='subset of raw images to work with (default all)')
    _parser.add_argument('-b', '--block_num', action='store', type=str, default=None,
                         help='Block number to use if multiple image blocks exist in directory.')
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
    return _parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    register_ortho(args.fn_ortho, args.fn_ref, args.fn_reldem, args.fn_dem,
                   glacmask=args.glacmask,
                   landmask=args.landmask,
                   footprints=args.footprints,
                   im_subset=args.im_subset,
                   block_num=args.block_num,
                   ori=args.ori,
                   ortho_res=args.ortho_res,
                   imgsource=args.imgsource,
                   density=args.density,
                   allfree=args.no_allfree)


if __name__ == "__main__":
    main()
