import os
import argparse
from glob import glob
from skimage.io import imread, imsave
from pybob.bob_tools import mkdir_p
from spymicmac.image import balance_image


def _argparser():
    _parser = argparse.ArgumentParser(description="Apply Contrast-limited Adaptive Histogram Equalization (CLAHE) to all re-sampled images in current directory.",
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    return _parser


def main():
    mkdir_p('balanced')

    imlist = glob('OIS*.tif')
    imlist.sort()

    for im in imlist:
        print(im)
        img_filt = balance_image(imread(im))
        imsave(os.path.join('balanced', im), img_filt)


if __name__ == "__main__":
    main()
