import argparse
from spymicmac.micmac import post_process


def _argparser():
    helpstr = """
    Apply georeferencing and masking to the final DEM and Correlation images (optionally, the orthomosaic as well).

    Output files are written as follows:
        - DEM: post_processed/{out_name}_Z.tif
        - Hillshade: post_processed/{out_name}_HS.tif
        - Correlation: post_processed/{out_name}_CORR.tif
        - Orthomosaic: post_processed/{out_name}_Ortho.tif
    """
    parser = argparse.ArgumentParser(description=helpstr,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('projstr', action='store', type=str,
                        help="A string corresponding to the DEM's CRS that GDAL can use "
                             "to georeference the rasters, or an EPSG code.")
    parser.add_argument('out_name', action='store', type=str,
                        help='The name that the output files should have.')
    parser.add_argument('dirmec', action='store', type=str,
                        help='The MEC directory to process files from (e.g., MEC-Malt)')
    parser.add_argument('--do_ortho', action='store_true',
                        help='Post-process the orthomosaic in Ortho-{dirmec}, as well. Assumes '
                             'that you have run mm3d Tawny with Out=Orthophotomosaic first.')
    parser.add_argument('--ind_ortho', action='store_true',
                        help='Post-process the individual orthophotos in  in Ortho-{dirmec}, as well.')
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    # if projstr is an epsg code, try to convert to int; otherwise, keep as str
    try:
        args.projstr = int(args.projstr)
    except ValueError as e:
        pass

    post_process(**vars(args))


if __name__ == "__main__":
    main()

