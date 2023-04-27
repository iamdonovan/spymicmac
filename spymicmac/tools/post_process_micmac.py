import argparse
from spymicmac.micmac import post_process


def _argparser():
    helpstr = """
    """
    parser = argparse.ArgumentParser(description=helpstr,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('projstr', action='store', type=str,
                        help="A string corresponding to the DEM's CRS that GDAL can use to georeference the rasters.")
    parser.add_argument('out_name', action='store', type=str,
                        help='The name that the output files should have.')
    parser.add_argument('dirmec', action='store', type=str,
                        help='The MEC directory to process files from (e.g., MEC-Malt)')
    parser.add_argument('--do_ortho', action='store_true',
                        help='Post-process the orthomosaic in Ortho-{dirmec}, as well. Assumes that you have run'
                             'mm3d Tawny with Out=Orthophotomosaic first.')
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    post_process(**vars(args))


if __name__ == "__main__":
    main()

