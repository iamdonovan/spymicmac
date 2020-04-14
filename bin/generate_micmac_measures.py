#!/usr/bin/env python
import argparse
import sPyMicMac.micmac_tools as mt


def _argparser():
    parser = argparse.ArgumentParser(description="", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-joined', action='store_true',
                        help='generate files for joined scene (220x460 mm) instead of half (220x230mm)')
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()
    mt.generate_measures_files(joined=args.joined)


if __name__ == "__main__":
    main()
