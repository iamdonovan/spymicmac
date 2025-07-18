"""
spymicmac.preprocessing is a collection of tools for preprocessing images
"""
import os
from pathlib import Path
import multiprocessing as mp
import shutil
import tarfile
import numpy as np
import pandas as pd
from glob import glob
from skimage import io, filters, exposure
from . import micmac, resample, matching, image
from typing import Union

from .matching import find_reseau_grid


def initialize_kh9_mc(add_sfs: bool = False, cam_csv: str = 'camera_defs.csv', overwrite: bool = False) -> None:
    """
    Initialize the following files needed for processing KH-9 MC images, if they do not already exist:

    - MicMac-LocalChantierDescripteur.xml
    - Ori-Init/AutoCal_Foc-304800_KH9MC.xml
    - Ori-InterneScan/MeasuresCamera.xml

    If multiple "cameras" (i.e., images acquired from different missions) are being used, these should be
    defined in a cam_csv file, as generated by micmac.generate_multicam_csv().

    :param add_sfs: use SFS to help find tie points in low-contrast images
    :param cam_csv: the name of the CSV file containing camera information
    :param overwrite: overwrite existing files
    """
    # first, check whether we have multiple cameras
    is_multicam = os.path.exists(cam_csv)

    if is_multicam:
        cameras = pd.read_csv(cam_csv)

        if not os.path.exists('MicMac-LocalChantierDescripteur.xml'):
            micmac.create_localchantier_xml(add_sfs=add_sfs, cam_csv=cam_csv)

        for cam in cameras.itertuples():
            foc = int(cam.focal * 1000)
            if not os.path.exists(Path('Ori-Init', f"AutoCal_Foc-{foc}_{cam.name}.xml")):
                micmac.init_autocal(framesize=(cam.width, cam.height), foc=cam.focal, camname=cam.name)

    else:
        if not os.path.exists(Path('Ori-Init', 'AutoCal_Foc-304800_KH9MC.xml')) or overwrite:
            micmac.init_autocal()

        if not os.path.exists('MicMac-LocalChantierDescripteur.xml') or overwrite:
            micmac.create_localchantier_xml(add_sfs=add_sfs)

    if not os.path.isfile(Path('Ori-InterneScan', 'MeasuresCamera.xml')) or overwrite:
        os.makedirs('Ori-InterneScan', exist_ok=True)
        micmac.generate_measures_files(joined=True)
        shutil.move('MeasuresCamera.xml', 'Ori-InterneScan')


def extract(tar_ext: str = '.tgz') -> list:
    """
    Extract image parts from tar files in the current directory.

    :param tar_ext: the extension of the files to look for.
    :returns: **tarlist** -- the list of tarfiles found in the current directory.
    """
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
        tarlist = list(set([os.path.splitext(fn)[0].split('_')[0] for fn in glob('D*.tif')]))
        tarlist.sort()

    return tarlist


def check_reseau() -> bool:
    """
    Check whether all KH-9 Mapping Camera images in the current directory have a corresponding Measures{fn_img}.xml
    file in Ori-InterneScan.

    :return: True if the files exist; False if they do not exist
    """
    imlist = glob('DZB*.tif')
    measlist = glob('MeasuresIm*.tif.xml', dir_fd='Ori-InterneScan')

    # if there are no DZB*.tif, but there are MeasuresIm files, we've probably done this step
    if len(imlist) == 0 and len(measlist) > 0:
        return True
    else:
        # check that each image file has a matching measuresim file
        return all([any([im in meas for meas in measlist]) for im in imlist])


def _reseau_wrapper(argsin: dict) -> None:
    find_reseau_grid(**argsin)


def batch_reseau(imlist: list, args) -> None:

    pool = mp.Pool(args['nproc'], maxtasksperchild=1)

    pool_args = [{'fn_img': fn_img + '.tif'} for fn_img in imlist]
    pool.map(_reseau_wrapper, pool_args, chunksize=1)

    pool.close()
    pool.join()


def _resample_wrapper(argsin: dict) -> None:
    print(argsin['fn_img'])
    resample.resample_hex(**argsin)

def batch_resample(imlist: list, args) -> None:
    pool = mp.Pool(args['nproc'], maxtasksperchild=1)

    arg_dict = {'scale': args['scale']}
    pool_args = [{'fn_img': fn_img + '.tif'} for fn_img in imlist]

    for d in pool_args:
        d.update(arg_dict)

    pool.map(_resample_wrapper, pool_args, chunksize=1)
    pool.close()
    pool.join()


def _handle_steps(proc_steps, steps, skips, opt_steps=[], option=None):

    if type(steps) is str:
        assert steps in set(proc_steps + ['all']), f"steps must be one of {['all'] + proc_steps}"
    else:
        for step in steps:
            assert step in proc_steps, f"{step} not recognized"

    if type(skips) is str:
        assert skips in set(proc_steps + ['none']), f"skipped steps must be one of {['none'] + proc_steps}"
    else:
        for step in skips:
            assert step in proc_steps, f"{step} not recognized"

    if steps == 'all':
        do_step = [True] * len(proc_steps)
    else:
        if type(steps) is not list: steps = [ steps ]
        do_step = [step in steps for step in proc_steps]

    if type(option) is str:
        assert option in set(opt_steps + ['all', 'none']), f"option must be one of {['all', 'none'] + opt_steps}"
    else:
        for opt in option:
            assert opt in opt_steps, f"{opt} not recognized"

    if option == 'all':
        do_step += [True] * len(opt_steps)
    elif option == 'none':
        do_step += [False] * len(opt_steps)
    else:
        if type(option) is not list: option = [ option ]
        do_step += [step in option for step in opt_steps]

    do = dict(zip(proc_steps + opt_steps, do_step))

    if skips != 'none':
        if type(skips) is not list: skips = [ skips ]
        for step in skips:
            if step in do.keys():
                do.update({step: False})

    return do


def preprocess_kh9_mc(steps: Union[str, list] = 'all', skip: Union[str, list] = 'none',
                      option: Union[str, list] = 'none', nproc: Union[int, str] = 1,
                      add_sfs: bool = False, cam_csv: str = 'camera_defs.csv', tar_ext: str = '.tgz',
                      is_reversed: bool = False, overlap: int = 2000, block_size: int = 2000, blend: bool = False,
                      scale: int = 70, clip_limit: float = 0.005, res_low: int = 400, res_high: int = 1200,
                      camera_model: str = 'RadialExtended', ori: str = 'Relative', init_cal: str = 'Init',
                      lib_foc: bool = False, lib_pp: bool = False, lib_cd: bool = False, add_params: bool = False,
                      ) -> None:
    """
    Run pre-processing steps for KH-9 Hexagon Mapping Camera images. By default, runs all steps (equivalent
    to steps='all'):

    - extract: extracts images from tar files (skips if no tar files are found)
    - join: joins scanned image halves
    - reseau: finds reseau marker locations in the joined image
    - erase: erases reseau markers from image
    - resample: resamples images to common size using the reseau marker locations
    - tapioca: calls mm3d Tapioca MulScale to find tie points
    - tapas: calls mm3d Tapas to calibrate camera model, find relative image orientation
    - aperi: calls mm3d AperiCloud to create point cloud using calibrated camera model

    Additional optional steps can be included using the 'option' argument:

    - filter: use a 1-sigma gaussian filter to smooth the images before resampling. Done before resampling the images.
    destripe: remove horizontal/vertical scanner-induced stripes from images
    - balance: use contrast-limited adaptive histogram equalization (clahe) to improve contrast in the image. Done
        after resampling the images.
    - schnaps: calls mm3d Schnaps to clean/filter tie points. Done after calling Tapioca and before calling Tapas.

    To run steps individually, pass them to the steps argument as a list. For example, to only run the 'reseau' and
    'erase' steps:

        preprocess_kh9_mc(steps=['reseau', 'erase'])

    Similarly, to run the 'reseau', 'erase', and 'balance' steps:

        preprocess_kh9_mc(steps=['reseau', 'erase'], option='balance')

    :param steps: The default pre-processing steps to run
    :param skip: The pre-processing steps to skip
    :param option: The optional pre-processing steps to run
    :param nproc: The number of sub-processes to use - either an integer value, or 'max'. If 'max',
        uses mp.cpu_count() to determine the total number of processors available.
    :param add_sfs: use SFS to help find tie points in low-contrast images
    :param cam_csv: Name of the CSV file containing camera information
    :param tar_ext: Extension for tar files
    :param is_reversed: Order of image parts is reversed (a is the right side of the image)
    :param overlap: The amount of overlap between image parts to use to search for matches
    :param block_size: The number of rows each search sub-block should cover
    :param blend: Blend across image halves to prevent a sharp line at edge
    :param scale: The scale of the resampled images, in pixels per mm
    :param clip_limit: Clipping limit, for contrast-limited adaptive histogram equalization
    :param res_low: The size of the largest image axis, in pixels, for low-resolution matching with Tapioca
    :param res_high: The size of the largest image axis, in pixels, for high-resolution matching with Tapioca
    :param camera_model: The camera calibration model to use for Tapas
    :param ori: The output Ori directory to create using Tapas
    :param init_cal: The initial calibration Ori to use for Tapas
    :param lib_foc: Use LibFoc=1 for mm3d Tapas
    :param lib_pp: Use LibPP=1 for mm3d Tapas
    :param lib_cd: Use LibCD=1 for mm3d Tapas
    :param add_params: Add decentric and affine parameters to the camera model
    """
    assert isinstance(nproc, int) or nproc == 'max', f"nproc must be an integer or 'max': {nproc}"

    if nproc == 'max':
        nproc = mp.cpu_count()
        print(f"Using {nproc} processors for steps that use multiprocessing.")

    proc_steps = ['extract', 'join', 'reseau', 'erase', 'resample', 'tapioca', 'tapas', 'aperi']
    opt_steps = ['filter', 'destripe', 'balance', 'schnaps']

    do = _handle_steps(proc_steps, steps, skip, opt_steps=opt_steps, option=option)

    # initialize the xml files needed
    initialize_kh9_mc(add_sfs=add_sfs, cam_csv=cam_csv)

    # now, do the requested steps, in order
    if do['extract']:
        imlist = extract(tar_ext)
        imlist = [tfile.split(tar_ext)[0] for tfile in imlist]
    elif any([do['join'], do['reseau'], do['erase'], do['resample']]):
        imlist = list(set([os.path.splitext(fn)[0].split('_')[0] for fn in glob('DZB*.tif')]))
        imlist.sort()
    else:
        imlist = [os.path.splitext(fn.split('OIS-Reech_')[1])[0] for fn in glob('OIS*.tif')]
        imlist.sort()

    print(imlist)

    if do['join']:
        print('Joining scanned image halves.')
        os.makedirs('halves', exist_ok=True)

        half_list = [fn.split('_a.tif')[0] for fn in glob('DZB*_a.tif')]
        half_list.sort()

        if len(half_list) == 0:
            print('No image halves found, skipping.')
        else:
            for fn_img in half_list:
                print(fn_img)

                join_args = {'overlap': overlap,
                             'block_size': block_size,
                             'blend': blend,
                             'is_reversed': is_reversed}

                image.join_hexagon(fn_img, **join_args)
                shutil.move(fn_img + '_a.tif', 'halves')
                shutil.move(fn_img + '_b.tif', 'halves')

    if do['reseau']:
        # if we're doing all steps, check that we need to; if we explicitly asked to
        # do this step, then do it.
        if (steps == 'all' and not check_reseau()) or ('reseau' in steps):
            print('Finding Reseau marks in images.')

            if nproc > 1 and len(imlist) > 1:
                batch_reseau(imlist, locals())
            else:
                for fn_img in imlist:
                    matching.find_reseau_grid(fn_img + '.tif')

        else:
            print('All images have reseau marks found. To re-run, explicitly call with --steps reseau.')

    if do['erase']:
        orig_imlist = glob('DZB*.tif', dir_fd='original')
        measlist = glob('MeasuresIm*.tif.xml', dir_fd='Ori-InterneScan')
        if len(orig_imlist) == len(measlist):
            print('Reseau marks have already been erased, skipping.')

        else:
            print('Erasing Reseau marks from images.')
            for fn_img in imlist:
                print(fn_img)
                matching.remove_crosses(fn_img + '.tif', nproc=nproc)

    if do['filter']:
        print('Filtering images with a 1-sigma Gaussian Filter')
        for fn_img in imlist:
            print(fn_img)
            img = io.imread(fn_img + '.tif')
            filt = filters.gaussian(img, sigma=1, preserve_range=True).astype(np.uint8)
            io.imsave(fn_img + '.tif', filt)

    if do['resample']:
        os.makedirs('Orig', exist_ok=True)
        # now, resample the images
        if nproc > 1 and len(imlist) > 1:
            batch_resample(imlist, locals())
        else:
            for fn_img in imlist:
                print(f'Resampling {fn_img}')
                resample.resample_hex(fn_img + '.tif', scale=scale)

        for fn_img in imlist:
            shutil.move(fn_img + '.tif', 'Orig')

    if do['destripe'] and do['balance']:
        print('Removing scanner stripes by subtracting row/column median values.')
        for fn_img in imlist:
            print('OIS-Reech_' + fn_img)
            img = io.imread('OIS-Reech_' + fn_img + '.tif')
            img_adj = image.high_low_subtract(img)
            img_adj = image.balance_image(img_adj, clip_limit=clip_limit)
            io.imsave('OIS-Reech_' + fn_img + '.tif', img_adj.astype(np.uint8))
    elif do['destripe']:
        print('Removing scanner stripes by subtracting row/column median values.')
        for fn_img in imlist:
            print('OIS-Reech_' + fn_img)
            img = io.imread('OIS-Reech_' + fn_img + '.tif')
            img_adj = image.high_low_subtract(img)
            io.imsave('OIS-Reech_' + fn_img + '.tif', img_adj.astype(np.uint8))
    elif do['balance']:
        print('Using CLAHE to balance image contrast')
        for fn_img in imlist:
            print('OIS-Reech_' + fn_img)
            img = io.imread('OIS-Reech_' + fn_img + '.tif')
            img_adj = image.balance_image(img, clip_limit=clip_limit)
            io.imsave('OIS-Reech_' + fn_img + '.tif', img_adj.astype(np.uint8))

    # run tapioca
    if do['tapioca']:
        exit_code = micmac.tapioca(res_low=res_low, res_high=res_high)
        if exit_code != 0:
            raise RuntimeError('Error in mm3d Tapioca - check Tapioca output for details.')

    # run schnaps
    if do['schnaps']:
        exit_code = micmac.schnaps("OIS.*tif")
        if exit_code != 0:
            raise RuntimeError('Error in mm3d Schnaps - check Schnaps output for details.')

    # run tapas
    if do['tapas']:
        tapas_kwargs = {'in_cal': init_cal, 'lib_foc': lib_foc, 'lib_pp': lib_pp, 'lib_cd': lib_cd}
        if do['schnaps']:
            tapas_kwargs['dir_homol'] = 'Homol_mini'

        exit_code = micmac.tapas(camera_model, ori, **tapas_kwargs)
        if exit_code != 0:
            raise RuntimeError('Error in mm3d Tapas - check Tapas output for details.')

    if add_params:
        cam_xmls = glob('AutoCal*.xml', root_dir='Ori-Relative')
        for fn_cam in cam_xmls:
            cam = micmac.load_cam_xml(Path('Ori-Relative', fn_cam))
            micmac.write_cam_xml(Path('Ori-Relative', fn_cam), cam)

    # run apericloud
    if do['aperi']:
        exit_code = micmac.apericloud('Relative')
        if exit_code != 0:
            raise RuntimeError('Error in mm3d AperiCloud - check AperiCloud output for details.')


def preprocess_pan(flavor: str, steps: Union[str, list] = 'all', skip: Union[str, list] = 'none',
                   tar_ext: str = '.tgz', is_reversed: bool = False, overlap: int = 2000, block_size: int = 2000,
                   blend: bool = False, marker_size: int = 31, factor: Union[int, None] = None,
                   clip_limit: float = 0.005) -> None:
    """
    Run pre-processing steps for declassified panoramic images (e.g., KH-4 and KH-9 Panoramic Camera). By default, runs
    all steps (equivalent to steps='all'):

    - extract: extracts images from tar files (skips if no tar files are found)
    - join: joins scanned image halves
    - filter: use a 1-sigma gaussian filter to smooth the images before cropping
    - crop: rotate and crop images to remove image frame
    - balance: use contrast-limited adaptive histogram equalization (clahe) to improve contrast in the image

    To run steps individually, pass them to the steps argument as a list. For example, to only run the 'join' and
    'crop' steps:

        preprocess_pan(steps=['join', 'crop'])

    :param flavor: The camera type (KH4 or KH9)
    :param steps: The pre-processing steps to run
    :param skip: The pre-processing steps to skip
    :param tar_ext: Extension for tar files
    :param is_reversed: Order of image parts is reversed (a is the right side of the image)
    :param overlap: The amount of overlap between image parts to use to search for matches
    :param block_size: The number of rows each search sub-block should cover
    :param blend: Blend across image halves to prevent a sharp line at edge
    :param marker_size: The size of the wagon wheel or rail markers to identify in the image
    :param factor: The number by which to divide the image width and height to scale the image (default: do not scale)
    :param clip_limit: Clipping limit, for contrast-limited adaptive histogram equalization
    """
    assert flavor in ['KH4', 'KH9'], "flavor must be one of [KH4, KH9]"

    pattern_dict = {'KH4': 'DS', 'KH9': 'D3C'}
    patt = pattern_dict[flavor]

    proc_steps = ['extract', 'join', 'filter', 'crop', 'balance']

    # check what steps we're running
    do = _handle_steps(proc_steps, steps, skip)

    # now, do all the steps we were asked to do, in order
    if do['extract']:
        imlist = extract(tar_ext)
        imlist = [tfile.split(tar_ext)[0] for tfile in imlist]
        if len(imlist) == 0:
            imlist = list(set([os.path.splitext(fn)[0].split('_')[0] for fn in glob(f'{patt}*.tif')]))
            imlist.sort()
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
            parts_list = image.get_parts_list(fn_img)
            if len(parts_list) < 2:
                continue

            join_args = {'overlap': overlap,
                         'block_size': block_size,
                         'blend': blend,
                         'is_reversed': is_reversed}

            if flavor == 'KH4' and overlap is None:
                join_args.update({'overlap': 8000})
            elif flavor == 'KH9' and overlap is None:
                join_args.update({'overlap': 1000})

            image.join_hexagon(fn_img, **join_args)

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
            border, angle = resample.crop_panoramic(fn_img + '.tif', flavor, marker_size=marker_size,
                                                    fact=factor, return_vals=True)
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
            img_adj = 255 * exposure.equalize_adapthist(img, clip_limit=clip_limit)
            io.imsave('OIS-Reech_' + fn_img + '.tif', img_adj.astype(np.uint8))

