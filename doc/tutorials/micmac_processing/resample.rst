re-sampling the images
======================

After you have found each of the fiducial marks in each image and generated a MeasuresIm file for each image,
either using :doc:`../../spymicmac/scripts/find_reseau_grid` for KH-9 images, or by hand/using ``mm3d Kugelhupf``
for historic aerial photographs, you can run ``ReSampFid``:
::

    *****************************
    *  Help for Elise Arg main  *
    *****************************
    Mandatory unnamed args :
      * string :: {Pattern image}
      * REAL :: {Resolution of scan, mm/pix}
    Named args :
      * [Name=BoxCh] Box2dr :: {Box in Chambre (generally in mm, [xmin,ymin,xmax,ymax])}
      * [Name=Kern] INT :: {Kernel of interpol,0 Bilin, 1 Bicub, other SinC (fix size of apodisation window), Def=5}
      * [Name=AttrMasq] string :: {Atribut for masq toto-> toto_AttrMasq.tif, NONE if unused, Def=NONE}
      * [Name=ExpAff] bool :: {Export the affine transformation}

For example, if you are using KH-9 images from Earth Explorer, you would run the following command to re-sample
the images to 14 microns (0.014 mm):
::

    mm3d ReSampFid "DZB.*tif" 0.014

.. note::
    There is also a shell script, :doc:`../../spymicmac/scripts/resample_hexagon`, that will leave a 1 mm overlap
    between the two image halves, to aid in joining them together.


The re-sampled images will have OIS-Reech\_ appended to the filename:
::

    AR5840034159994.tif -> OIS-Reech_AR5840034159994.tif

These are the images that you will use for the remaining steps - you might want to create a new folder to place the
original images.

The next step is to find tie points using ``Tapioca``.
