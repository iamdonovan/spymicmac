import argparse
from spymicmac.micmac import write_xml


def _argparser():
    _parser = argparse.ArgumentParser(description="Given a GDAL dataset, create a MicMac xml worldfile.",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    _parser.add_argument('filename', action='store', type=str,
                        help='the filename of the image.')
    _parser.add_argument('-m', '--mask', action='store', type=str, default='./MEC-Malt/Masq_STD-MALT_DeZoom1.tif',
                        help='Path to mask file [./MEC-Malt/Masq_STD-MALT_DeZoom1.tif]')
    _parser.add_argument('-g', '--geom', action='store', type=str, default='eGeomMNTEuclid',
                        help='MicMac Geometry name [eGeomMNTEuclid]')
    return _parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    write_xml(args.filename, args.mask, args.geom)


if __name__ == "__main__":
    main()
