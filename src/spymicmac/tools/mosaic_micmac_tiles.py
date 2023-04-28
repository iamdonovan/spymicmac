import argparse
from spymicmac.micmac import mosaic_micmac_tiles


def _argparser():
    parser = argparse.ArgumentParser(description="Re-stitch images tiled by MicMac.",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename', action='store', type=str,
                        help="MicMac filename to mosaic together")
    parser.add_argument('-imgdir', action='store', type=str, default='.',
                        help="Directory containing images to Mosaic (default: .)")
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    mosaic_micmac_tiles(args.filename, dirname=args.imgdir)


if __name__ == "__main__":
    main()
