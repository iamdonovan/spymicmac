#!/usr/bin/env python
import argparse
from glob import glob
from spymicmac.image import join_hexagon


def _argparser():
    parser = argparse.ArgumentParser(description="Join parts of a scanned image",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-p', '--pattern', action='store', type=str, default='DZB',
                        help='Match pattern for images [DZB]')
    parser.add_argument('-o', '--overlap', action='store', type=int, default=2000,
                        help='overlap search width between two images [2000]')
    parser.add_argument('-k', '--block_size', action='store', type=int, default=None,
                        help='the number of rows each sub-block should cover. Defaults to overlap value.')
    parser.add_argument('-b', '--blend', action='store_true',
                        help='Blend across image halves to prevent a sharp line at edge.')
    parser.add_argument('-r', '--reversed', action='store_true',
                        help='parts are in reversed order (i.e., part b is the left part, part a is the right part)')
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    imlist = [f.split('_a.tif')[0] for f in glob('{}*a.tif'.format(args.pattern))]
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
