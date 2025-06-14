#!/usr/bin/env python
import argparse
from glob import glob
from spymicmac.image import join_hexagon


def _argparser():
    parser = argparse.ArgumentParser(description="Join multiple parts of a scanned image.",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-p', '--pattern', action='store', type=str, default='DZB',
                        help='Match pattern for images (default: DZB)')
    parser.add_argument('-o', '--overlap', action='store', type=int, default=2000,
                        help='the overlap, in pixels, between the image parts (default: 2000)')
    parser.add_argument('-k', '--block_size', action='store', type=int, default=None,
                        help='the number of rows each sub-block should cover. Defaults to overlap value.')
    parser.add_argument('-b', '--blend', action='store_true',
                        help='apply a linear blend between the two scanned halves.')
    parser.add_argument('-r', '--reversed', action='store_true',
                        help='parts are in reversed order (i.e., part b is the left part, part a is the right part)')
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    imlist = [f.split('_a.tif')[0] for f in glob(f"{args.pattern}*a.tif")]
    imlist.sort()

    for im in imlist:
        print(im)
        join_hexagon(im,
                     overlap=args.overlap,
                     block_size=args.block_size,
                     blend=args.blend,
                     is_reversed=args.reversed)


if __name__ == "__main__":
    main()
