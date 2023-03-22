computing the absolute orthophoto and DEM
=========================================

Once you have successfully run :py:meth:`spymicmac.register.register_ortho`, you will have an orientation folder,
``Ori-TerrainFirstPass`` [#]_. This is the folder to pass to ``Malt``:

.. code-block:: sh

    mm3d Malt Ortho "OIS.*tif" TerrainFirstPass DirMEC=MEC-Malt NbVI=2 MasqImGlob=filtre.tif ZoomF=1 DefCor=0 CostTrans=1 EZA=1

Just like with the :doc:`relative` step, once this command finishes, you will have two new directories: ``MEC-Malt``
and ``Ortho-MEC-Malt``. The DEM and associated correlation masks are found in ``MEC-Malt``, while the
orthophotos are found in ``Ortho-MEC-Malt``.

Depending on the scale of your images and how many you are using, you will probably need to mosaic the tiles
of the DEM, and possibly the correlation masks [#]_:

.. code-block:: sh

    cd MEC-Malt
    mosaic_micmac_tiles.py -filename Z_Num9_DeZoom1_STD-Malt
    mosaic_micmac_tiles.py -filename Correl_STD-MALT_Num_8

Once the tiles have been mosaicked, you can also run **post_process_micmac.sh** to apply the AutoMask and add a CRS
to the DEM.

creating the orthomosaic using Tawny
------------------------------------
Just like with the :doc:`relative` step, the orthoimages are not mosaicked - they are just the individual images
orthorectified using the extracted DEM - you'll need to run ``Tawny`` again to mosaic the images:

.. code-block:: sh

    mm3d Tawny Ortho-MEC-Malt Out=Orthophotomosaic.tif RadiomEgal=0

Again, we are using ``RadiomEgal=0`` to use the images as-is, rather than attempting to balance the radiometry (as this
can lead to undesirable results). Finally, you might need to re-combine the image tiles, depending on how large they
are:

.. code-block:: sh

    cd Ortho-MEC-Malt
    mosaic_micmac_tiles.py -filename Orthophotomosaic

At this point, you should have a finished DEM and orthomosaic. You may want to check the accuracy of your DEM by
co-registering it to a DEM of known quality. You may also wish to remove residual
`doming effects <https://doi.org/10.5194/isprs-annals-V-3-2020-375-2020>`_ using ``mm3d PostProc Banana``.

You can also run :doc:`../../spymicmac/scripts/post_process_micmac` to apply the AutoMask to the DEM and
georeference the correlation mask:

.. code-block:: sh

    cd MEC-Malt
    post_process_micmac.sh -z "8 +north" -n Block1


.. [#] If you are running :py:meth:`spymicmac.register.register_ortho` with ``block`` set to a value (e.g., ``-b 1``), this
    will be ``Ori-TerrainFirstPass_block1``. You should also change ``DirMEC`` to a different name, (e.g., ``MEC-Malt_block1``),
    otherwise it will be overwritten with each new block that you run.

.. [#] Note that in some cases with very large image blocks, the final steps will actually be ``Z_Num10_DeZoom1_STD-MALT``
    and ``Correl_STD-MALT_Num_9`` - it's probably a good idea to check this!
