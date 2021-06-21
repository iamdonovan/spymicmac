#!/usr/bin/env python
import argparse
from spymicmac.image import remove_crosses


def _argparser():
    parser = argparse.ArgumentParser(description="Remove Reseau marks from KH-9 image(s).",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('img', action='store', type=str, nargs='+', help='Image(s) to remove crosses from.')
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    for im in args.img:
        print(im)
        remove_crosses(im)


if __name__ == "__main__":
    main()
