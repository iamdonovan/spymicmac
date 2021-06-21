#!/usr/bin/env python
import argparse
import spymicmac.micmac as mt


def _argparser():
    parser = argparse.ArgumentParser(description="Create id_fiducial.txt, MeasuresCamera.xml, and Tmp-SL-Glob.xml files for KH-9 Hexagon images.",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-joined', action='store_true',
                        help='generate files for joined scene (220x460 mm) instead of half (220x230mm)')
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()
    mt.generate_measures_files(joined=args.joined)


if __name__ == "__main__":
    main()
