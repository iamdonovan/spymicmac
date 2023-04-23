import os
import argparse
from glob import glob
from skimage.io import imread, imsave
from spymicmac.image import balance_image


def _argparser():
    _parser = argparse.ArgumentParser(description="Apply Contrast-limited Adaptive Histogram Equalization (CLAHE) to all re-sampled images in current directory.",
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    return _parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    os.makedirs('balanced', exist_ok=True)

    imlist = glob('OIS*.tif')
    imlist.sort()

    for im in imlist:
        print(im)
        img_filt = balance_image(imread(im))
        imsave(os.path.join('balanced', im), img_filt)


if __name__ == "__main__":
    main()
