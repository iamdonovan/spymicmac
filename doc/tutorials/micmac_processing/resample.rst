re-sampling the images
======================

air photos
------------

After you have found each of the fiducial marks in each image and generated a MeasuresIm file for each image, either by
hand or using ``mm3d Kugelhupf``, you can run ``ReSampFid``:

.. code-block:: text

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

For example, to re-sample the images to 14 microns (0.014 mm):

.. code-block:: sh

    mm3d ReSampFid "AR5.*tif" 0.014


The re-sampled images will have OIS-Reech\_ appended to the filename:

.. code-block:: text

    AR5840034159994.tif -> OIS-Reech_AR5840034159994.tif

These are the images that you will use for the remaining steps - you might want to create a new folder to place the
original images.

kh-9 hexagon mapping camera
-----------------------------

To resample KH-9 Hexagon Mapping Camera images, use either :py:meth:`spymicmac.resample.resample_hex` or
:doc:`../../spymicmac/scripts/resample_hexagon`.


next step
----------

The next step is to find tie points using ``Tapioca``.
