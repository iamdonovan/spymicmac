#!/usr/bin/env python
import argparse
from spymicmac.micmac import move_bad_tapas


def _argparser():
    parser = argparse.ArgumentParser(description="Read Tapas residuals and move images with NaN residual.",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('ori', action='store', type=str, help='Orientation directory to read.')

    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    move_bad_tapas(args.ori)


if __name__ == "__main__":
    main()
