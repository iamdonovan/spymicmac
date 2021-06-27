#!/usr/bin/env python
import argparse
from spymicmac.micmac import remove_measure


def _argparser():
    parser = argparse.ArgumentParser(description="Remove GCP(s) from a Measures xml file.",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fn_meas', action='store', type=str, nargs='+', help='xml file to remove GCP(s) from.')
    parser.add_argument('gcp', action='store', type=str, nargs='+', help='GCP name(s) to remove from <fn_meas>.')
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    for gcp in args.gcps:
        remove_measure(args.fn_meas, gcp)


if __name__ == "__main__":
    main()
