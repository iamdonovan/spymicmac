import argparse
from spymicmac.preprocessing import preprocess_pan


def _argparser():
    helpstr = """
    Run pre-processing steps for declassified panoramic images (e.g., KH-4 and KH-9 Panoramic Camera). By default, runs 
    all steps (equivalent to "--steps all"):

    - extract: extracts images from tar files (skips if no tar files are found)
    - join: joins scanned image halves
    - filter: use a 1-sigma gaussian filter to smooth the images before cropping
    - crop: rotate and crop images to remove image frame
    - balance: use contrast-limited adaptive histogram equalization (clahe) to improve contrast in the image

    To run steps individually, use the --steps flag with the corresponding step name(s). For example, to only run the
    'join' and 'crop' steps:

    preprocess_pan <flavor> --steps join crop <additional arguments>

    """
    parser = argparse.ArgumentParser(description=helpstr,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('flavor', action='store', type=str,
                        help='The camera type (KH4 or KH9)')
    parser.add_argument('--steps', action='store', type=str, nargs='+', default='all',
                        help='The pre-processing steps to run.')
    parser.add_argument('--skip', action='store', type=str, nargs='+', default='none',
                        help='The pre-processing steps to skip.')
    parser.add_argument('--tar_ext', action='store', type=str, default='.tgz',
                        help='Extension for tar files (default: .tgz)')
    parser.add_argument('-r', '--is_reversed', action='store_true',
                        help='Order of image parts is reversed (a is the right side of the image). (default: False)')
    parser.add_argument('-o', '--overlap', action='store', type=int, default=None,
                        help='The amount of overlap between image parts to use to search for matches. Default depends'
                             'on flavor: KH4 -> 8000, KH9 -> 1000')
    parser.add_argument('-k', '--block_size', action='store', type=int, default=2000,
                    help='the number of rows each search sub-block should cover [2000].')
    parser.add_argument('-b', '--blend', action='store_true',
                        help='Blend across image halves to prevent a sharp line at edge.')
    parser.add_argument('-m', '--marker_size', action='store', type=int, default=31,
                        help='The size of the wagon wheel markers to identify in the image (default: 31 px)')
    parser.add_argument('-s', '--factor', action='store', type=int, default=None,
                        help='The number by which to divide the image width and height to scale the image '
                             '(default: do not scale)')
    parser.add_argument('--clip_limit', action='store', type=float, default=0.005,
                        help='Clipping limit, for contrast-limited adaptive histogram equalization. (default: 0.005)')
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    preprocess_pan(**vars(args))


if __name__ == "__main__":
    main()
