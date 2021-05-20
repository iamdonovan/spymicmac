import os
from glob import glob
from skimage.io import imread, imsave
from skimage.exposure import equalize_adapthist
from skimage.filters import median
from skimage.morphology import disk
import numpy as np
from pybob.bob_tools import mkdir_p

mkdir_p('balanced')

imlist = glob('OIS*.tif')
imlist.sort()

for im in imlist:
    print(im)
    img = imread(im)
    img_eq = (255 * equalize_adapthist(img)).astype(np.uint8)
    img_filt = median(img_eq, selem=disk(1))
    imsave(os.path.join('balanced', im), img_filt)
