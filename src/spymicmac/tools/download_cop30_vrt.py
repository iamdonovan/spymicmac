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
                         help='filename for image footprints. If not set, uses '
                              'spymicmac.data.get_usgs_footprints to download footprints'
                              'based on imlist.')
    _parser.add_argument('-imgsource', action='store', type=str, default='DECLASSII',
                         help='the EarthExplorer Dataset name for the images (default: DECLASSII)')
    _parser.add_argument('-g', '--globstr', action='store', type=str, default='OIS.*tif',
                         help='the search string to use to find images in the '
                              'current directory. (default: OIS.*tif)')
    _parser.add_argument('-crs', action='store', type=int, default=None,
                         help='the epsg code for the CRS to re-project the mosaic to. '
                              'If not set, uses EPSG:4326 (WGS84 Lat/Lon).')
    _parser.add_argument('--geoid_height', action='store_true',
                         help='do not convert the elevations from height above EGM2008 Geoid '
                              'to height above WGS84 Ellipsoid.')

    return _parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    if args.footprints is not None:
        footprints = gpd.read_file(args.footprints)
    else:
        footprints = None

    download_cop30_vrt(
        imlist=args.imlist,
        footprints=footprints,
        imgsource=args.imgsource,
        globstr=args.globstr,
        crs=args.crs,
        to_ellipsoid=(not args.geoid_height)
    )


if __name__ == "__main__":
    main()
