#!/usr/bin/env python
import os
from glob import glob
import argparse
import pyvips
from sPyMicMac.image_tools import join_halves


def _argparser():
    parser = argparse.ArgumentParser(description="", formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('img', action='store', type=str, nargs='+', help='split image(s) to join and balance')
    parser.add_argument('-overlap', action='store', type=int, nargs=4,
                        metavar=('c_left', 'r_left', 'c_right', 'r_right'), default=[33145, 16000, 144, 16000],
                        help=('overlap point between images. Default is central overlap point, ' \
                              'assuming 233x224 mm images output from mm3d ReSampFid with 0.007 mm/px resolution. ' \
                              'Given as [col_left, row_left, col_right, row_right].'))
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    for im in args.img:
        imname = im.rsplit('_', 1)[0]
        print(imname)

        join_halves(imname, args.overlap)

        img = pyvips.Image.new_from_file('{}.tif'.format(imname), memory=True)
        img_bal = img.hist_equal()
        img_bal.write_to_file('{}.tif'.format(imname))


if __name__ == "__main__":
    main()
