import os
import argparse
import shutil
import tarfile
from glob import glob
from spymicmac.image import join_hexagon, find_reseau_grid, resample_hex, remove_crosses
from spymicmac import micmac


def extract(tar_ext='.tgz'):
    # make a directory to move the tarballs to
    os.makedirs('tarballs', exist_ok=True)

    tarlist = glob('*' + tar_ext)
    tarlist.sort()

    if len(tarlist) > 0:
        print('Extracting images from tar files.')
        for tarball in tarlist:
            print(tarball)
            with tarfile.open(tarball, 'r') as tfile:
                tfile.extractall('.')
            shutil.move(tarball, 'tarballs')
    else:
        print('No tar files found, skipping.')
        # we still want to make sure that we have a list of images to work with.
        tarlist = list(set([os.path.splitext(fn)[0].split('_')[0] for fn in glob('DZB*.tif')]))
        tarlist.sort()

    return tarlist


def _argparser():
    helpstr = """
    Run pre-processing steps for KH-9 Hexagon Mapping Camera images. By default, runs all steps (equivalent to 
    "--steps all"):
    
    - extract: extracts images from tar files (skips if no tar files are found)
    - join: joins scanned image halves
    - reseau: finds reseau marker locations in the joined image
    - erase: erases reseau markers from image
    - resample: resamples images to common size using the reseau marker locations
    - tapioca: calls mm3d Tapioca MulScale to find tie points
    - tapas: calls mm3d Tapas to calibrate camera model, find relative image orientation
    - aperi: calls mm3d AperiCloud to create point cloud using calibrated camera model

    To run steps individually, use the --steps flag with the corresponding step name(s). For example, to only run the
    'reseau' and 'erase' steps:
    
    preprocess_kh9 --steps reseau erase <additional arguments>

    """
    parser = argparse.ArgumentParser(description=helpstr,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--steps', action='store', type=str, nargs='+', default='all',
                        help='The pre-processing steps to run.')
    parser.add_argument('--tar_ext', action='store', type=str, default='.tgz',
                        help='Extension for tar files (default: .tgz)')
    parser.add_argument('-s', '--scale', action='store', type=int, default=70,
                        help='The scale of the resampled images, in pixels per mm. (default: 70)')
    parser.add_argument('-b', '--blend', action='store_true',
                        help='Blend across image halves to prevent a sharp line at edge.')
    parser.add_argument('--add_sfs', action='store_true',
                        help='use SFS to help find tie points in low-contrast images [False]')
    parser.add_argument('--res_low', action='store', type=int, default=400,
                        help='the size of the largest image axis, in pixels, '
                             'for low-resolution matching with Tapioca (default: 400)')
    parser.add_argument('--res_high', action='store', type=int, default=1200,
                        help='the size of the largest image axis, in pixels, '
                             'for low-resolution matching with Tapioca (default: 1200)')
    parser.add_argument('--camera_model', action='store', type=str, default='RadialExtended',
                        help='The camera calibration model to use for Tapas (default: RadialExtended)')
    parser.add_argument('--ori', action='store', type=str, default='Relative',
                        help='The output Ori directory to create using Tapas (default: Relative)')
    parser.add_argument('--init_cal', action='store', type=str, default='Init',
                        help='The initial calibration Ori to use for Tapas (default: Init)')
    parser.add_argument('--lib_foc', action='store_true',
                        help='Use LibFoc=1 for mm3d Tapas (default: LibFoc=0)')
    parser.add_argument('--lib_pp', action='store_true',
                        help='Use LibPP=1 for mm3d Tapas (default: LibPP=0)')
    parser.add_argument('--lib_cd', action='store_true',
                        help='Use LibCD=1 for mm3d Tapas (default: LibCD=0)')
    parser.add_argument('-n', '--nproc', type=int, default=1, help='number of sub-processes to use (default: 1).')
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    proc_steps = ['extract', 'join', 'reseau', 'erase', 'resample', 'tapioca', 'tapas', 'aperi']

    # check what steps we're running
    if args.steps == 'all':
        do_step = [True] * len(proc_steps)
    else:
        do_step = [step in args.steps for step in proc_steps]
    do = dict(zip(proc_steps, do_step))

    # create the necessary xml files
    micmac.init_autocal()
    micmac.create_localchantier_xml(add_sfs=args.add_sfs)
    micmac.generate_measures_files(joined=True)

    os.makedirs('Ori-InterneScan', exist_ok=True)
    if not os.path.isfile(os.path.join('Ori-InterneScan', 'MeasuresCamera.xml')):
        shutil.move('MeasuresCamera.xml', 'Ori-InterneScan')

    # now, do all the steps we were asked to do, in order
    if do['extract']:
        imlist = extract(args.tar_ext)
        imlist = [tfile.split(args.tar_ext)[0] for tfile in imlist]
    else:
        imlist = list(set([os.path.splitext(fn)[0].split('_')[0] for fn in glob('DZB*.tif')]))
        imlist.sort()

    if do['join']:
        print('Joining scanned image halves.')
        os.makedirs('halves', exist_ok=True)

        for fn_img in imlist:
            print(fn_img)
            join_hexagon(fn_img, blend=args.blend)
            shutil.move(fn_img + '_a.tif', 'halves')
            shutil.move(fn_img + '_b.tif', 'halves')

    if do['reseau']:
        print('Finding Reseau marks in images.')
        for fn_img in imlist:
            find_reseau_grid(fn_img + '.tif')

    if do['erase']:
        print('Erasing Reseau marks from images.')
        for fn_img in imlist:
            print(fn_img)
            remove_crosses(fn_img + '.tif', nproc=args.nproc)

    if do['resample']:
        os.makedirs('Orig', exist_ok=True)
        # now, resample the images
        for fn_img in imlist:
            print('Resampling {}'.format(fn_img))
            resample_hex(fn_img + '.tif', scale=70)
            shutil.move(fn_img + '.tif', 'Orig')

    # run tapioca
    if do['tapioca']:
        exit_code = micmac.tapioca(res_low=args.res_low, res_high=args.res_high)
        if exit_code != 0:
            raise RuntimeError('Error in mm3d Tapioca - check Tapioca output for details.')

    # run tapas
    if do['tapas']:
        exit_code = micmac.tapas(args.camera_model, args.ori, lib_foc=args.lib_foc,
                                 lib_pp=args.lib_pp, lib_cd=args.lib_cd)
        if exit_code != 0:
            raise RuntimeError('Error in mm3d Tapas - check Tapas output for details.')

    # run apericloud
    if do['aperi']:
        exit_code = micmac.apericloud('Relative')
        if exit_code != 0:
            raise RuntimeError('Error in mm3d AperiCloud - check AperiCloud output for details.')


if __name__ == "__main__":
    main()
