tl;dr guide
=======================

This page is designed to provide a "bare bones" user guide to spymicmac/MicMac for processing KH-9 Hexagon Mapping
Camera images or scanned air photos. It is presented primarily as a list of commands to run; if you would like a
longer explanation of the different steps, see the additional guides :doc:`here <index>`.

KH-9 processing
-----------------

To begin, run :doc:`../../spymicmac/scripts/generate_micmac_measures` to generate ``id_fiducial.txt`` and
``MeasuresCamera.xml``.

Next, create a new directory, ``Ori-InterneScan``, and move ``MeasuresCamera.xml`` into it.

Finally, run :py:meth:`spymicmac.micmac.create_localchantier_xml` to create ``MicMac-LocalChantierDescripteur.xml``,
or :download:`download <../examples/MicMac-LocalChantierDescripteur.xml>` the example provided.

.. note::

    The rest of this guide assumes that you have the following directory structure:

    .. code-block:: text

        project
        ├── id_fiducial.txt
        ├── Img1_a.tif
        ├── Img1_b.tif
        ├── Img2_a.tif
        ├── Img2_b.tif
        ...
        ├── MicMac-LocalChantierDescripteur.xml
        ├── Ori-InterneScan
        │   └── MeasuresCamera.xml

    where ``Img1_a.tif`` and ``Img1_b.tif`` are the two halves of one of the KH-9 images, ``Img2_a.tif`` and
    ``Img2_b.tif`` are two halves of the next image in the acquisition, and so on. For stereo processing, you need
    a minimum of two images.


joining image halves
^^^^^^^^^^^^^^^^^^^^^

To join image halves use the command line tool :doc:`../../spymicmac/scripts/join_hexagon`:

.. code-block:: sh

    join_hexagon

This will search for all files that fit the pattern ``DZB*N001_(a|b).tif``, and join them into a single image,
``DZB*N001.tif``. After joining the images, create a new directory, ``halves``, and move the image halves into it.

reseau grid
^^^^^^^^^^^^

Next, use :doc:`../../spymicmac/scripts/find_reseau_grid` to find the location of the 1081 Reseau marks in the image:

.. code-block:: sh

    find_reseau_grid DZB*.tif

This will create a file, ``Ori-InterneScan/MeasuresIm-DZB*N001.tif.xml``, for each of the images in the directory. It
will also create an image in ``match_imgs`` that shows the estimated distortion pattern.

removing crosses (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To erase the Reseau marks from each image, run :doc:`../../spymicmac/scripts/remove_crosses`:

.. code-block:: sh

    remove_crosses DZB*.tif

This will move the original images to a new directory, ``original``, and save the erased images to the working
directory.

resampling
^^^^^^^^^^^

After this step, you can use :doc:`../../spymicmac/scripts/resample_hexagon` to resample the images to the same size,
and correct the distortion using the Reseau marks:

.. code-block:: sh

    resample_hexagon DZB*.tif -s SCALE

where scale is the scale of the resampled image, in pixels per mm (default is 70). This will create a new file,
``OIS-Reech_DZB*N001.tif`` for each of the original images.


general micmac processing
---------------------------

Once the images have been re-sampled, the rest of the workflow is largely the same for both KH-9 images and scanned
air photos.

tapioca
^^^^^^^^^

Once the images are re-sampled, run ``mm3d Tapioca MulScale`` to find tie points in downscaled versions of the images:

.. code-block:: sh

    mm3d Tapioca MulScale "OIS.*tif" 400 1200

In the above command, the numbers after "OIS.*tif" are the size, in pixels, of the longest dimension of the downscaled
image.


tapas
^^^^^^

To find the relative orientation of the images, and calibrate the camera parameters, use ``mm3d Tapas``:

.. code-block:: sh

    mm3d Tapas RadialBasic "OIS.*tif" Out=Relative LibFoc=0

The ``LibFoc=0`` flag will keep the focal length fixed at the value provided in ``MicMac-LocalChantierDescripteur.xml``.


relative dem and orthomosaic
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a relative orthophoto and DEM using ``mm3d Malt``:

.. code-block:: sh

    mm3d Malt Ortho "OIS.*tif" Relative DirMEC=MEC-Relative NbVI=2 ZoomF=8 DefCor=0 CostTrans=1 EZA=1

If you have used an image mask, run the following command instead:

.. code-block:: sh

    mm3d Malt Ortho "OIS.*tif" Relative DirMEC=MEC-Relative NbVI=2 MasqImGlob=filtre.tif ZoomF=8 DefCor=0 CostTrans=1 EZA=1

By default, ``mm3d Malt`` only orthorectifies the individual images; to create an orthomosaic, use the ``mm3d Tawny``
command:

.. code-block:: sh

    mm3d Tawny Ortho-MEC-Malt Out=Orthophotomosaic.tif RadiomEgal=0

.. note::

    If the image is very large, you may need to run ``mosaic_micmac_tiles.py`` to combine the tiles into a single
    image. For the relative DEM, this will probably not be needed.

    For the orthophoto, run the following from the ``Ortho-MEC-Relative`` directory:

    .. code-block:: sh

        mosaic_micmac_tiles.py -filename Orthophoto


registering the image
^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

    This step requires a DEM and an orthoimage to find control points and estimate the absolute orientation of the
    images. In the examples below, I assume that these are named ``DEM.tif`` and ``Landsat.tif``, respectively.


The main command to run here is :doc:`../../spymicmac/scripts/register_ortho`:

.. code-block:: sh

    register_ortho Ortho-MEC-Malt/Orthophotomosaic.tif Landsat.tif MEC-Malt/Z_Num6_DeZoom8_STD-MALT.tif DEM.tif

.. note::

    If you have a shapefile of the image footprints, use the ``-footprints`` flag; otherwise, they will be downloaded from
    USGS Earth Explorer:

    .. code-block:: sh

        register_ortho Ortho-MEC-Malt/Orthophotomosaic.tif Landsat.tif MEC-Malt/Z_Num6_DeZoom8_STD-MALT.tif DEM.tif -footprints Footprints.shp

    The shapefile should have a polygon for each image, with the name of the original image (minus the file extension)
    included in an ``ID`` column.

dem and orthomosaic
^^^^^^^^^^^^^^^^^^^^^

After the absolute orientation has been estimated by registering the image, run ``mm3d Malt`` again with this new
orientation to extract the final DEM and orthophoto:

.. code-block:: sh

    mm3d Malt Ortho "OIS.*tif" TerrainFirstPass DirMEC=MEC-Malt NbVI=2 ZoomF=1 DefCor=0 CostTrans=1 EZA=1

If you have used an image mask, run the following command instead:

.. code-block:: sh

    mm3d Malt Ortho "OIS.*tif" TerrainFirstPass DirMEC=MEC-Malt NbVI=2 MasqImGlob=filtre.tif ZoomF=1 DefCor=0 CostTrans=1 EZA=1

To generate an orthomosaic, run the following command:

.. code-block:: sh

    mm3d Tawny Ortho-MEC-Malt Out=Orthophotomosaic.tif RadiomEgal=0

.. note::

    If the image is very large, you may need to run ``mosaic_micmac_tiles.py`` to combine the tiles into a single
    image. For the DEM, run the following from the ``MEC-Relative`` directory:

    .. code-block:: sh

        mosaic_micmac_tiles.py

    And for the orthophoto, run the following from the ``Ortho-MEC-Relative`` directory:

    .. code-block:: sh

        mosaic_micmac_tiles.py -filename Orthophoto


You can also run :doc:`../../spymicmac/scripts/post_process_micmac` to apply the AutoMask to the DEM and
Orthomosaic, and georeference the correlation mask.

And that's it. You should now have an orthorectified KH-9 (or air photo) mosaic, and a DEM. Enjoy.
