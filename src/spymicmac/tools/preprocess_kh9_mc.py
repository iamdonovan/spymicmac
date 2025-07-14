import argparse
import multiprocessing as mp
from spymicmac.preprocessing import preprocess_kh9_mc


def _argparser():
    helpstr = """
    Run pre-processing steps for KH-9 Hexagon Mapping Camera images. By default, runs all steps (equivalent to 
    "--steps all"):
    
    - extract: extracts images from tar files (skips if no tar files are found)
    - join: joins scanned image halves
    - reseau: finds reseau marker locations in the joined image
    - erase: erases reseau markers from image
    - filter: use a 1-sigma gaussian filter to smooth the images before resampling
    - resample: resamples images to common size using the reseau marker locations
    - balance: use contrast-limited adaptive histogram equalization (clahe) to improve contrast in the image
    - tapioca: calls mm3d Tapioca MulScale to find tie points
    - schnaps: calls mm3d Schnaps to clean/filter tie points    
    - tapas: calls mm3d Tapas to calibrate camera model, find relative image orientation
    - aperi: calls mm3d AperiCloud to create point cloud using calibrated camera model

    Additional optional steps can be included using the 'option' argument:

    - filter: use a 1-sigma gaussian filter to smooth the images before resampling. Done before resampling the images.
    - balance: use contrast-limited adaptive histogram equalization (clahe) to improve contrast in the image. Done
        after resampling the images.
    - schnaps: calls mm3d Schnaps to clean/filter tie points. Done after calling Tapioca and before calling Tapas.

    To run steps individually, use the --steps flag with the corresponding step name(s). For example, to only run the
    'reseau' and 'erase' steps:
    
        preprocess_kh9_mc --steps reseau erase <additional arguments>

    Similarly, to run the 'reseau', 'erase', and 'balance' steps:
    
        preprocess_kh9_mc --steps reseau erase --option balance <additional arguments>

    """
    parser = argparse.ArgumentParser(description=helpstr,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--steps', action='store', type=str, nargs='+', default='all',
                        help='The default pre-processing steps to run (default: all).')
    parser.add_argument('--skip', action='store', type=str, nargs='+', default='none',
                        help='The pre-processing steps to skip (default: none).')
    parser.add_argument('--option', action='store', type=str, nargs='+', default='none',
                        help='The optional pre-processing steps to skip (default: none).')
    parser.add_argument('-n', '--nproc', type=str, default=1,
                        help='The number of sub-processes to use - either an integer value, or max. '
                             'If max, uses mp.cpu_count() to determine the total number of processors '
                             'available (default: 1).')
    parser.add_argument('--add_sfs', action='store_true',
                        help='use SFS to help find tie points in low-contrast images')
    parser.add_argument('--cam_csv', action='store', type=str, default='camera_defs.csv',
                        help='Name of the CSV file containing camera information (default: camera_defs.csv)')
    parser.add_argument('--tar_ext', action='store', type=str, default='.tgz',
                        help='Extension for tar files (default: .tgz).')
    parser.add_argument('-r', '--is_reversed', action='store_true',
                        help='Order of image parts is reversed (a is the right side of the image). (default: False)')
    parser.add_argument('-o', '--overlap', action='store', type=int, default=2000,
                        help='The amount of overlap between image parts to use to search for matches. (default: 2000).')
    parser.add_argument('-k', '--block_size', action='store', type=int, default=2000,
                    help='the number of rows each search sub-block should cover (default: 2000).')
    parser.add_argument('-b', '--blend', action='store_true',
                        help='Blend across image halves to prevent a sharp line at edge.')
    parser.add_argument('-s', '--scale', action='store', type=int, default=70,
                        help='The scale of the resampled images, in pixels per mm (default: 70).')
    parser.add_argument('--clip_limit', action='store', type=float, default=0.005,
                        help='Clipping limit, for contrast-limited adaptive histogram equalization. (default: 0.005)')
    parser.add_argument('--res_low', action='store', type=int, default=400,
                        help='the size of the largest image axis, in pixels, '
                             'for low-resolution matching with Tapioca (default: 400).')
    parser.add_argument('--res_high', action='store', type=int, default=1200,
                        help='the size of the largest image axis, in pixels, '
                             'for low-resolution matching with Tapioca (default: 1200).')
    parser.add_argument('--camera_model', action='store', type=str, default='RadialExtended',
                        help='The camera calibration model to use for Tapas (default: RadialExtended)')
    parser.add_argument('--ori', action='store', type=str, default='Relative',
                        help='The output Ori directory to create using Tapas (default: Relative)')
    parser.add_argument('--init_cal', action='store', type=str, default='Init',
                        help='The initial calibration Ori to use for Tapas (default: Init)')
    parser.add_argument('--lib_foc', action='store_true',
                        help='Use LibFoc=1 for mm3d Tapas')
    parser.add_argument('--lib_pp', action='store_true',
                        help='Use LibPP=1 for mm3d Tapas')
    parser.add_argument('--lib_cd', action='store_true',
                        help='Use LibCD=1 for mm3d Tapas')
    parser.add_argument('--add_params', action='store_true',
                        help='Add decentric and affine parameters to the camera model')

    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    try:
        args.nproc = int(args.nproc)
    except ValueError as e:
        assert args.nproc == 'max', f"-nproc must either be max or a integer value: {args.nproc}"
        args.nproc = mp.cpu_count()

    preprocess_kh9_mc(**vars(args))


if __name__ == "__main__":
    main()
