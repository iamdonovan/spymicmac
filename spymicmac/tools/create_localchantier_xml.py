#!/usr/bin/env python
import argparse
from spymicmac.micmac import create_localchantier_xml


def _argparser():
    parser = argparse.ArgumentParser(description="Create a MicMac-LocalChantierDescripteur.xml file",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-n', '--name', action='store', type=str, default='KH9MC',
                        help='Camera name to use. [Default: KH9MC]')
    parser.add_argument('-s', '--short_name', action='store', type=str, default='KH-9 Hexagon Mapping Camera',
                        help='Short description of camera. [Default: KH-9 Hexagon Mapping Camera]')
    parser.add_argument('--film_size', action='store', type=float, nargs=2, default=(460, 220),
                        help='Film size (width, height), in mm. [Default: 460, 220]')
    parser.add_argument('-p', '--pattern', action='store', type=str, default='.*',
                        help='Camera name to use. [Default: .*]')
    parser.add_argument('-f', '--focal', action='store', type=float, default=304.8,
                        help='Camera focal length (in mm). [Default: 304.8]')
    parser.add_argument('--add_sfs', action='store_true',
                        help='Use SFS for tie point matching [Default: False].')
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    create_localchantier_xml(**vars(args))


if __name__ == "__main__":
    main()

