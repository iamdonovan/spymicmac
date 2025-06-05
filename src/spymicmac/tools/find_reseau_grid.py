#!/usr/bin/env python
import argparse
import multiprocessing as mp
import numpy as np
from spymicmac.matching import find_reseau_grid


def batch_wrapper(argsin):
    find_reseau_grid(**argsin)


def _argparser():
    _parser = argparse.ArgumentParser(description="Find the locations of the Réseau marks in scanned KH-9 images. "
                                                  "Locations are saved to Ori-InterneScan/MeasuresIm-{fn_img}.xml.",
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    _parser.add_argument('fn_img', action='store', type=str, nargs='+',
                         help='Image(s) to find Réseau marks in.')
    _parser.add_argument('-csize', action='store', type=int, default=361,
                         help='the size of the cross template (default: 361 -> 361x361).')
    _parser.add_argument('-n', '--nproc', type=int, default=1,
                         help='number of sub-processes to use (default: 1).')
    return _parser


def main():
    np.seterr(divide='ignore', invalid='ignore')
    parser = _argparser()
    args = parser.parse_args()

    if args.nproc > 1 and len(args.img) > 1:
        pool = mp.Pool(args.nproc, maxtasksperchild=1)

        arg_dict = {'csize': args.csize}
        pool_args = [{'fn_img': fn_img} for fn_img in args.img]

        for d in pool_args:
            d.update(arg_dict)

        pool.map(batch_wrapper, pool_args, chunksize=1)
        pool.close()
        pool.join()
    else:
        for fn_img in args.img:
            find_reseau_grid(fn_img, csize=args.csize)


if __name__ == "__main__":
    main()
