computing tie points
====================

The basic command for computing tie points is `Tapioca <https://micmac.ensg.eu/index.php/Tapioca>`_. If you don't
have very many images, you can use ``Tapioca`` to find the tie points by matching all pairs of images:

.. code-block:: sh

    mm3d Tapioca MulScale "OIS.*tif" 400 1200

kh-9 hexagon mapping camera
-----------------------------

When working with KH-9 images, the single ``Tapioca`` call will usually be sufficient, and you can move on to the
next step: :doc:`tapas`.

air photos
-----------

the neighbours file
.....................

With a large number of images, this will be a very slow process. If you have vector data (e.g., a shapefile) of
the image footprints, you can use :py:meth:`spymicmac.micmac.write_neighbour_images` to narrow down the number of
image pairs where ``Tapioca`` will search for pairs. This will create a file, ``FileImagesNeighbour.xml``, that specifies
which images overlap based on their footprints.

For more modern images where more precise location information is available, you can also use the ``OriConvert`` tool:

.. code-block:: sh

    mm3d OriConvert OriTxtInFile GPS_sommets.txt Sommets NameCple=FileImagesNeighbour.xml

Then, you can run ``Tapioca`` using the ``File`` option:

.. code-block:: sh

    mm3d Tapioca File FileImagesNeighbour.xml 1200

creating a mask
.....................

You can also create a mask to filter out tie points that are created due to the presence of fiducial marks or
inscriptions on the images. The basic command for this is:

.. code-block:: sh

    mm3d SaisieMasqQT "OIS-Reech_<Img>.tif"

A slightly more detailed instructional video can be found `here <https://youtu.be/xOHEkKiiRnM>`_. Once you have created
the mask, be sure to rename the file:

.. code-block:: sh

    mv OIS-Reech_<Img>_Masq.tif filtre.tif

filtering tie points
.....................

Once you have created a mask, you can use it to filter tie points:

.. code-block:: sh

    mm3d HomolFilterMasq "OIS.*tif" GlobalMasq=filtre.tif

Once this is done, you can move on to computing the relative orientation using ``Tapas``.
