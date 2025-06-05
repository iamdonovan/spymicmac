tl;dr guide
=======================

This page is designed to provide a "bare bones" user guide to spymicmac/MicMac for processing KH-9 Hexagon mapping
camera images or scanned air photos. It is presented primarily as a list of commands to run; if you would like a
longer explanation of the different steps, see the additional guides :doc:`here <index>`.

.. note::

    If even this is too verbose for your tastes, you can view the sample "one big script" :doc:`here <script>`.

setup
-----

Pre-processing and processing with MicMac requires a few files to be present in the directory - I recommend reading
the :doc:`setup tutorial <tutorials/setup>` first.

.. note::

    If you are working with KH-9 Hexagon mapping camera images, you can skip to the pre-processing section below, as
    the convenience tool :py:meth:`spymicmac.preprocessing.preprocess_kh9_mc` (command-line tool:
    :doc:`spymicmac/scripts/preprocess_kh9_mc`) will generate these files for the KH-9 Hexagon mapping camera.


.. note::

    The rest of this guide assumes that you have the following directory structure:

    .. code-block:: text

        project
        ├── id_fiducial.txt
        ├── Img1.tif
        ├── Img2.tif
        ...
        ├── MicMac-LocalChantierDescripteur.xml
        ├── Ori-InterneScan
        │   └── MeasuresCamera.xml

    where ``Img{n}.tif`` are the original scanned images.

.. note::

    You will need a minimum of two images that overlap at least in part, because this is stereo photogrammetry.


pre-processing
--------------

Before DEM and Orthophoto processing can happen, scanned images need to be *pre-processed*:

- geometrically, so that they have a common format
- radiometrically, to improve contrast and reduce noise

How these steps are done using ``spymicmac`` depends on what type of images you are working with.

KH-9 mapping camera
^^^^^^^^^^^^^^^^^^^

Pre-processing KH-9 mapping camera images requires a number of different steps.

The convenience function :py:meth:`spymicmac.preprocessing.preprocess_kh9_mc` (command-line tool:
:doc:`spymicmac/scripts/preprocess_kh9_mc`) will perform all of these steps in order. All you need to do is run the
command/function from the same directory where you have the ``.tgz`` (or ``.tar.gz``) files.

Otherwise, the required steps are implemented in the following functions/modules, which should be run in order:

- joining image scans together (:py:meth:`spymicmac.image.join_hexagon`)
- identifying réseau marks (:py:meth:`spymicmac.matching.find_reseau_grid`)
- resampling the images to a common format using the réseau mark locations (:py:meth:`spymicmac.resample.resample_hex`)

The re-sampled images will have **OIS-Reech\\_** appended to the filename, e.g.:

    DZB1214-500206L002001.tif -> OIS-Reech_DZB1214-500206L002001.tif

These are the images that you will use for the remaining steps - you might want to create a new folder to place the
original images.

Optional steps include:

- cross removal (removing/erasing réseau markers) - this should be done before resampling the images
  (:py:meth:`spymicmac.matching.remove_crosses`)
- radiometric pre-processing (e.g., de-noising with a gaussian filter or other contrast enhancements using
  :py:mod:`spymicmac.image`)


panoramic cameras
^^^^^^^^^^^^^^^^^

As with the KH-9 mapping camera images, declassified panoramic camera images can be pre-processed using the
convenience function :py:meth:`spymicmac.preprocessing.preprocess_pan` (command-line tool:
:doc:`spymicmac/scripts/preprocess_pan`). All you need to do is run the command/function from the same directory where
you have the ``.tgz`` (or ``.tar.gz``) files.

The required steps are:

- joining image scans together (:py:meth:`spymicmac.image.join_hexagon`)
- resampling the images using the image frame/border (:py:meth:`spymicmac.resample.crop_panoramic`)

Optional steps include:

- radiometric pre-processing (de-noising with a gaussian filter, other contrast enhancements using :py:mod:`spymicmac.image`)

aerial photos
^^^^^^^^^^^^^

Before processing scanned aerial images, you need to resample them to a common geometry using
`mm3d ReSampFid <https://micmac.ensg.eu/index.php/ReSampFid>`_. The easiest way to do this is using the fiducial
markers which are usually visible around the frame of the image.

You can manually identify the location of the fiducial markers using :py:meth:`spymicmac.micmac.batch_saisie_fids`:

.. code-block:: python

    from spymicmac.micmac import batch_saisie_fids
    from glob import glob

    imlist = glob('*.tif')
    batch_saisie_fids(imlist, flavor='qt')

If you know the type of camera you are using, you can also use :py:mod:`spymicmac.matching` to try to automatically
locate the fiducial markers in each image.

After you have identified the fiducial marker locations, use ``mm3d ReSampFid`` to re-sample the images. For example,
to resample your images to a resolution of 14 microns (0.014 mm per pixel):

.. code-block:: sh

    mm3d ReSampFid ".*tif" 0.014

The re-sampled images will have **OIS-Reech\\_** appended to the filename, e.g.:

    AR5840034159994.tif -> OIS-Reech_AR5840034159994.tif

These are the images that you will use for the remaining steps - you might want to create a new folder to place the
original images.


relative geometry
-------------------

.. note::

    Panoramic camera processing with MicMac is not currently supported; however, these can be processed by following
    the `example guides <https://stereopipeline.readthedocs.io/en/latest/examples.html>`__ from Ames Stereo Pipeline
    (ASP).

    :py:mod:`spymicmac.asp` has a number of functions that can be used to help set this processing up.

.. note::

    If you have used :py:meth:`spymicmac.preprocessing.preprocess_kh9_mc` or :doc:`spymicmac/scripts/preprocess_kh9_mc`,
    the tie points and camera calibration and orientation steps will all be run, and you can skip to computing the
    relative DEM and orthophoto.

Once the images have been re-sampled, the rest of the workflow is largely the same for both KH-9 mapping camera images
and scanned air photos.

tie points
^^^^^^^^^^

First, you need to compute tie points for your re-sampled images, using ``mm3d Tapioca``
(or :py:meth:`spymicmac.micmac.tapioca`):

.. code-block:: sh

    mm3d Tapioca MulScale "OIS.*tif" 400 1200

.. code-block:: python

    from spymicmac import micmac
    micmac.tapioca('OIS.*tif', res_low=400, res_high=1200)


camera calibration and orientation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next, you need to compute the relative orientation of the images, and calibrate the intrinsic camera parameters using
``mm3d Tapas`` (or :py:meth:`spymicmac.micmac.tapas`):

.. code-block:: sh

    mm3d Tapas FraserBasic "OIS.*tif" Out=Relative

.. code-block:: python

    from spymicmac import micmac
    micmac.tapas(
        cam_model='FraserBasic',
        ori_out='Relative',
        img_pattern='OIS.*tif'
    )

.. tip::

    If you have a large number of images, you might want to initialize the calibration using a small number of
    good-quality images, before working on the whole set.

.. tip::

    Using ``mm3d Martini`` (:py:meth:`spymicmac.micmac.martini`) can help by initializing the image orientation.

.. note::

    ``Tapas`` has a number of different camera calibration models available. You might wish to experiment with these
    to find one that works well.

visualizing the orientation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can visualize the orientation of the images using ``mm3d AperiCloud`` (or :py:meth:`spymicmac.micmac.aperi`):

.. code-block:: sh

    mm3d AperiCloud "OIS.*tif" Relative

.. code-block:: python

    from spymicmac import micmac
    micmac.apericloud('Relative', 'OIS.*tif')

dem and orthomosaic
^^^^^^^^^^^^^^^^^^^

If the orientation looks fine (no obvious camera outliers, point cloud looks roughly similar to the study area), you
can process a relative DEM/orthophotos using ``mm3d Malt`` (or :py:meth:`spymicmac.micmac.malt`):

.. code-block:: sh

    mm3d Malt Ortho "OIS.*tif" Relative DirMEC=MEC-Relative NbVI=2 ZoomF=4 DefCor=0 CostTrans=4 EZA=1 SzW=3 Regul=0.1

.. code-block:: python

    from spymicmac import micmac
    micmac.malt('OIS.*tif', 'Relative',
        zoomf=4,
        dirmec='MEC-Relative',
        cost_trans=4,
        szw=3,
        regul=0.1
    )

``mm3d Malt`` only orthorectifies the individual images; to create an orthomosaic, use the ``mm3d Tawny``
command (or :py:meth:`spymicmac.micmac.tawny`):

.. code-block:: sh

    mm3d Tawny Ortho-MEC-Relative Out=Orthophotomosaic.tif RadiomEgal=0

.. code-block:: python

    from spymicmac import micmac
    micmac.tawny('MEC-Relative', radiomegal=False)

.. note::

    If the image is very large, you may need to run :py:meth:`spymicmac.micmac.mosaic_micmac_tiles` (or
    :doc:`spymicmac/scripts/mosaic_micmac_tiles`) to combine the tiles into a single image:

    For the orthophoto, run the following from the project directory:

    .. code-block:: sh

        mosaic_micmac_tiles -filename Orthophotomosaic -imgdir Ortho-MEC-Relative

    .. code-block:: python

        from spymicmac import micmac
        micmac.mosaic_micmac_tiles('Orthophotomosaic', dirname='Ortho-MEC-Relative')

    For the relative DEM, this will probably not be needed.

You can view the DEM or orthomosaic in a GIS software - if it looks okay, you can move on to the next steps. If not,
you might need to work on the camera calibration or orientation.

absolute geometry
-------------------

Now that you have a DEM/orthophotos in relative geometry, you can use a reference image to transform them to an
absolute geometry, and compute the final DEM/orthophotos.

registering the image
^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

    This step requires a DEM to find control points and estimate the absolute orientation of the images. In the
    examples below, I assume that you have named this file ``DEM.tif``.

The main function to use here is :py:meth:`spymicmac.register.register_relative` (or the command-line tool:
:doc:`../../spymicmac/scripts/register_relative`):

.. code-block:: sh

    register_relative MEC-Malt DEM.tif

.. code-block:: python

    register.register_relative(
        'MEC-Relative',
        'DEM.tif',
    )

.. note::

    If you have a file containing the image footprints, use the ``-footprints`` flag; otherwise, they will automatically
    be loaded from ``Footprints.gpkg`` (if it exists), or they will be downloaded from USGS Earth Explorer:

    .. code-block:: sh

        register_relative MEC-Malt DEM.tif -footprints Footprints.gpkg

    .. code-block:: python

        register.register_relative(
            'MEC-Relative',
            'DEM.tif',
            footprints='Footprints.gpkg'
        )

    The file should contain a polygon for each image, with the name of the original image (minus the file extension)
    included in an ``ID`` column.

dem and orthomosaic
^^^^^^^^^^^^^^^^^^^^^

After the absolute orientation has been estimated by registering the image, run ``mm3d Malt`` again with this new
orientation to extract the final DEM and orthophoto:

.. code-block:: sh

    mm3d Malt Ortho "OIS.*tif" TerrainFinal DirMEC=MEC-Malt NbVI=2 ZoomF=1 DefCor=0 CostTrans=4 EZA=1 SzW=3 Regul=0.1

.. code-block:: python

    from spymicmac import micmac
    micmac.malt('OIS.*tif', 'TerrainFinal',
        zoomf=1,
        dirmec='MEC-Malt',
        cost_trans=4,
        szw=3,
        regul=0.1
    )

To generate the orthomosaic, run the following command (or, using python):

.. code-block:: sh

    mm3d Tawny Ortho-MEC-Malt Out=Orthophotomosaic.tif RadiomEgal=0

.. code-block:: python

    from spymicmac import micmac
    micmac.tawny('MEC-Relative', radiomegal=False)


post-processing
---------------

As a final step, run :py:meth:`spymicmac.micmac.post_process` (:doc:`spymicmac/scripts/post_process_micmac`) to apply
the AutoMask to the DEM and Orthomosaic:

.. code-block:: sh

    post_process_micmac epsg_code desired_prefix MEC-Malt --do_ortho

.. code-block:: python

    from spymicmac import micmac
    micmac.post_process(epsg_code, desired_prefix, 'MEC-Malt', do_ortho=True)

And that's it. You should now have an orthorectified KH-9 (or air photo) mosaic, and a DEM. Enjoy.
