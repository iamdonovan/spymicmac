#!/usr/bin/env python
import argparse
from spymicmac.micmac import create_localchantier_xml


def _argparser():
    parser = argparse.ArgumentParser(description="Create a MicMac-LocalChantierDescripteur.xml file for a given camera. "
                                                 "Default is the KH-9 Hexagon Mapping Camera.",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-n', '--name', action='store', type=str, default='KH9MC',
                        help='Camera name to use. (default: KH9MC).')
    parser.add_argument('-s', '--short_name', action='store', type=str, default='KH-9 Hexagon Mapping Camera',
                        help='Short description of camera. (default: KH-9 Hexagon Mapping Camera).')
    parser.add_argument('--film_size', action='store', type=float, nargs=2, default=(460, 220),
                        help='Film size (width, height), in mm. (default: 460 220).')
    parser.add_argument('-p', '--pattern', action='store', type=str, default='.*',
                        help='Camera name to use. (default: .*).')
    parser.add_argument('-f', '--focal', action='store', type=float, default=304.8,
                        help='Camera focal length in mm (default: 304.8).')
    parser.add_argument('--add_sfs', action='store_true',
                        help='Use SFS for tie point matching (default: False).')
    parser.add_argument('--cam_csv', action='store', type=str, default=None,
                        help='a CSV file containing parameters for multiple cameras (default: None).')
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    create_localchantier_xml(**vars(args))


if __name__ == "__main__":
    main()

