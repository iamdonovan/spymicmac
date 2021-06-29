#!/usr/bin/env python
import argparse
import numpy as np
from spymicmac.image import find_reseau_grid


def _argparser():
    _parser = argparse.ArgumentParser(description="Find Reseau marks in a scanned KH-9 Hexagon image.",
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    _parser.add_argument('fn_img', action='store', type=str, help='Image to find Reseau grid in.')
    _parser.add_argument('-csize', action='store', type=int, default=361, help='Reseau mark template size [361 pixels]')
    _parser.add_argument('-tsize', action='store', type=int, default=300, help='half-size of search window [300 pixels]')
    _parser.add_argument('-nproc', action='store', type=int, default=1,
                         help='Number of processors to use [1].')
    return _parser


def main():
    np.seterr(divide='ignore', invalid='ignore')
    parser = _argparser()
    args = parser.parse_args()

    find_reseau_grid(args.fn_img,
                     csize=args.csize,
                     tsize=args.tsize)


if __name__ == "__main__":
    main()
