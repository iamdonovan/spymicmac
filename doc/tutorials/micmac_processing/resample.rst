re-sampling the images
======================

After you have found each of the fiducial marks in each image and generated a MeasuresIm file for each image,
either using :doc:`../../sPyMicMac/scripts/find_reseau_shifts` for KH-9 images, or by hand/using ``mm3d Kugelhupf``
for historic aerial photographs, you can run ``ReSampFid``:
::

    mm3d ReSampFid <Pattern> <Resolution>

where *<Resolution>* is the pixel size in mm. For example, if you are using KH-9 images from Earth Explorer, you would
run the following command to re-sample the images to 14 microns (0.014 mm):
::

    mm3d ReSampFid "DZB.*tif" 0.014

The re-sampled images will have OIS-Reech\_ appended to the filename:
::

    AR5840034159994.tif -> OIS-Reech_AR5840034159994.tif

These are the images that you will use for the remaining steps - you might want to create a new folder to place the
original images.

The next step is to find tie points using ``Tapioca``.
