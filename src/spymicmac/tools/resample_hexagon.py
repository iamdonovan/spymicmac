#!/usr/bin/env python
import argparse
import multiprocessing as mp
from spymicmac.resample import resample_hex


def batch_wrapper(argsin):
    print(argsin['fn_img'])
    resample_hex(**argsin)


def _argparser():
    parser = argparse.ArgumentParser(description="Use a piecewise affine transformation to resample KH-9 images"
                                                 " using the reseau grid",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('img', action='store', type=str, nargs='+', help='Image(s) to resample from.')
    parser.add_argument('-s', '--scale', action='store', type=int, default=70,
                        help='The scale of the resampled image, in pixels per mm. [Default: 70]')
    parser.add_argument('-o', '--ori', action='store', type=str, default='InterneScan',
                        help='The Ori directory that contains both MeasuresCamera.xml and MeasuresIm for each image. '
                             '[Default: InterneScan]')
    parser.add_argument('-n', '--nproc', type=int, default=1, help='number of sub-processes to use [Default: 1].')
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    if args.nproc > 1 and len(args.img) > 1:
        pool = mp.Pool(args.nproc, maxtasksperchild=1)

        arg_dict = {'scale': args.scale, 'ori': args.ori}
        pool_args = [{'fn_img': fn_img} for fn_img in args.img]

        for d in pool_args:
            d.update(arg_dict)

        pool.map(batch_wrapper, pool_args, chunksize=1)
        pool.close()
        pool.join()
    else:
        for fn_img in args.img:
            print(fn_img)
            resample_hex(fn_img, scale=args.scale, ori=args.ori)


if __name__ == "__main__":
    main()
