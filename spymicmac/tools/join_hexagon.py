#!/usr/bin/env python
import argparse
from glob import glob
from spymicmac.image import join_hexagon


def _argparser():
    parser = argparse.ArgumentParser(description="Join two halves of a scanned KH-9 Hexagon image (or four parts of a scanned KH-4 Corona image).",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-o', '--overlap', action='store', type=int, default=2000,
                        help='overlap search width between two images [2000]')
    parser.add_argument('-p', '--pattern', action='store', type=str, default='DZB',
                       help='Match pattern for images [DZB]')
    parser.add_argument('-b', '--blend', action='store_true',
                        help='Blend across image halves to prevent a sharp line at edge.')
    parser.add_argument('-c', '--corona', action='store_true',
                        help='Images are Corona KH-4/KH-4A scans (i.e., there are 4 parts).')
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    imlist = [f.split('_a.tif')[0] for f in glob('{}*a.tif'.format(args.pattern))]
    for im in imlist:
        print(im)
        join_hexagon(im,
                     overlap=args.overlap,
                     blend=args.blend,
                     corona=args.corona)


if __name__ == "__main__":
    main()
