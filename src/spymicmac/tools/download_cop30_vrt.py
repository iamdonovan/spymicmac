import argparse
import geopandas as gpd
from spymicmac.data import download_cop30_vrt


def _argparser():
    helpstr = """
    Create a VRT using Copernicus 30m DSM tiles that intersect image footprints. Creates Copernicus_DSM.vrt 
    using files downloaded to cop30_dem/ within the current directory.
    """
    _parser = argparse.ArgumentParser(description=helpstr,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    _parser.add_argument('-imlist', action='store', type=str, nargs='+',
                         help='image(s) to use for geographic extent. If not set, '
                              'will search for images of form OIS*.tif')
    _parser.add_argument('-footprints', action='store', type=str,
                         help='filename for image footprints. By default, '
                              'downloads footprints from USGS Earth Explorer.')
    _parser.add_argument('-imgsource', action='store', type=str, default='DECLASSII',
                         help='the EE Dataset name for the images (default: DECLASSII)')

    return _parser


def main():
    parser = _argparser()
    args = parser.parse_args()
    print(args)

    if args.footprints is not None:
        footprints = gpd.read_file(args.footprints)
    else:
        footprints = None

    download_cop30_vrt(imlist=args.imlist,
                       footprints=footprints,
                       imgsource=args.imgsource)


if __name__ == "__main__":
    main()
