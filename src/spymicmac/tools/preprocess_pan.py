import os
import argparse
import shutil
import tarfile
import numpy as np
import pandas as pd
from glob import glob
from skimage import io, exposure, filters
from spymicmac.image import join_hexagon, get_parts_list
from spymicmac.preprocessing import extract
from spymicmac.resample import crop_panoramic


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
    parser.add_argument('-m', '--marker_size', action='store', type=int, default=31,
                        help='The size of the wagon wheel markers to identify in the image (default: 31 px)')
    parser.add_argument('-s', '--factor', action='store', type=int, default=None,
                        help='The number by which to divide the image width and height to scale the image '
                             '(default: do not scale)')
    parser.add_argument('-b', '--blend', action='store_true',
                        help='Blend across image halves to prevent a sharp line at edge.')
    parser.add_argument('--clip_limit', action='store', type=float, default=0.005,
                        help='Clipping limit, for contrast-limited adaptive histogram equalization. (default: 0.005)')
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    assert args.flavor in ['KH4', 'KH9'], "flavor must be one of [KH4, KH9]"
    pattern_dict = {'KH4': 'DS', 'KH9': 'D3C'}
    patt = pattern_dict[args.flavor]

    proc_steps = ['extract', 'join', 'filter', 'crop', 'balance']

    # check what steps we're running
    if args.steps == 'all':
        do_step = [True] * len(proc_steps)
    else:
        do_step = [step in args.steps for step in proc_steps]
    do = dict(zip(proc_steps, do_step))

    if args.skip != 'none':
        for step in args.skip:
            if step in do.keys():
                do.update({step: False})

    # now, do all the steps we were asked to do, in order
    if do['extract']:
        imlist = extract(args.tar_ext)
        imlist = [tfile.split(args.tar_ext)[0] for tfile in imlist]
    elif any([do['join'], do['filter'], do['crop']]):
        imlist = list(set([os.path.splitext(fn)[0].split('_')[0] for fn in glob(f'{patt}*.tif')]))
        imlist.sort()
    else:
        imlist = [os.path.splitext(fn.split('OIS-Reech_')[1])[0] for fn in glob('OIS*.tif')]
        imlist.sort()

    print(imlist)

    if do['join']:
        print('Joining scanned image parts.')
        os.makedirs('parts', exist_ok=True)

        for fn_img in imlist:
            print(fn_img)
            parts_list = get_parts_list(fn_img)
            if len(parts_list) < 2:
                continue

            join_hexagon(fn_img, overlap=8000, block_size=2000, blend=args.blend)
            for fn in parts_list:
                shutil.move(f'{fn_img}_{fn}.tif', 'parts')

    if do['filter']:
        print('Filtering images with a 1-sigma Gaussian Filter')
        for fn_img in imlist:
            print(fn_img)
            img = io.imread(fn_img + '.tif')
            filt = filters.gaussian(img, sigma=1, preserve_range=True).astype(np.uint8)
            io.imsave(fn_img + '.tif', filt)

    if do['crop']:
        print('Rotating and cropping images')
        os.makedirs('Orig', exist_ok=True)

        img_params = pd.DataFrame(columns=['fn_img', 'left', 'right', 'top', 'bot', 'angle'])

        for fn_img in imlist:
            print(fn_img)
            border, angle = crop_panoramic(fn_img + '.tif', args.flavor, marker_size=args.marker_size,
                                           fact=args.factor, return_vals=True)
            left, right, top, bot = border

            shutil.move(fn_img + '.tif', 'Orig')

            row = pd.Series()
            row['fn_img'] = fn_img
            row['left'] = int(left)
            row['right'] = int(right)
            row['top'] = int(top)
            row['bot'] = int(bot)
            row['angle'] = angle

            img_params.loc[len(img_params)] = row
            img_params.to_csv('crop_parameters.csv', index=False)

    if do['balance']:
        print('Using CLAHE to balance image contrast')
        for fn_img in imlist:
            print('OIS-Reech_' + fn_img)
            img = io.imread('OIS-Reech_' + fn_img + '.tif')
            img_adj = 255 * exposure.equalize_adapthist(img, clip_limit=args.clip_limit)
            io.imsave('OIS-Reech_' + fn_img + '.tif', img_adj.astype(np.uint8))


if __name__ == "__main__":
    main()
