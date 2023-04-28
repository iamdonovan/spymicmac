automatically finding control points
====================================
At this point, you should have successfully run ``mm3d Malt`` and ``mm3d Tawny`` to generate a relative orthophoto
and DEM. You should also have run :py:meth:`spymicmac.micmac.mosaic_micmac_tiles` (or
:doc:`../../spymicmac/scripts/mosaic_micmac_tiles`) if needed - otherwise, :py:meth:`spymicmac.register.register_ortho`
(or :doc:`../../spymicmac/scripts/register_ortho`) will most likely fail.

In order to run :py:meth:`spymicmac.register.register_ortho`, you will need a number of files, detailed in the section
below. After that, this document will describe the process that :py:meth:`spymicmac.register.register_ortho`
uses to find GCPs and iteratively refine the orientation.

necessary files
----------------

reference orthoimage and DEM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
At the risk of stating the obvious, the reference orthoimage and DEM should cover your study area. The reference
orthoimage can be a high-resolution orthomosaic, or it can be a comparatively low-resolution satellite image.
:py:meth:`spymicmac.register.register_ortho` has been tested on both Landsat-8 and Sentinel-2 images, with
satisfactory results for both.

In general, the results will depend in part on the accuracy of the control points - so if your reference orthoimage
is highly accurate and of a similar resolution to your air photos, the more accurate the result will be.

Similar rules apply for the reference DEM - the more accurate the DEM, the better the final result.

image footprints
^^^^^^^^^^^^^^^^^^
This should be a vector dataset (e.g., shapefile, geopackage - anything that can be read by
`geopandas.read_file <https://geopandas.org/docs/reference/api/geopandas.read_file.html>`_). The footprints do not have
to be highly accurate - most of the routines in :py:meth:`spymicmac.register` have been developed using USGS
datasets/metadata, which are only approximate footprints.

The main use for the footprints in :py:meth:`spymicmac.register.register_ortho` is in the call to
:py:meth:`spymicmac.orientation.transform_centers`, which uses
`RANSAC <https://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.ransac>`_
to estimate an affine transformation between the footprint centroids and the relative camera centers estimated
by ``mm3d Tapas``.

As long as the distribution of the footprints/centroids is approximately correct, this step
should work according to plan.

.. note::
    Some `USGS Aerial Photo Single Frames <https://doi.org/10.5066/F7610XKM>`_ (as well as the KH-9 Hexagon images)
    have footprints that are incorrectly placed, or the images have been scanned but not georeferenced. In this
    case, you may need to create your own footprints. It is most likely worth checking the metadata for these images
    before you start working.

exclusion and inclusion masks (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Finally, in areas with lots of water or terrain that has changed substantially (e.g., glaciers),
:py:meth:`spymicmac.register.register_ortho` can use both exclusion and inclusion masks to avoid searching for
matches on unstable terrain. Just like with the footprints, these files should be any data format that can be
read by `geopandas.read_file <https://geopandas.org/docs/reference/api/geopandas.read_file.html>`_.

The **exclusion** mask should be polygons of glaciers, landslides, or other terrain that should be *excluded* from
the search. Areas covered by this mask will be excluded from the search.

The **inclusion** mask should be polygons of land - any terrain that should be *included* in the search. Areas that
are not covered by this mask will be excluded from the search.

relative to absolute transformation
------------------------------------
The first step in :py:meth:`spymicmac.register.register_ortho` to use :py:meth:`spymicmac.orientation.transform_centers`
to transform between the relative and absolute spaces, using the centroids of the footprint polygons and the camera
positions estimated by ``mm3d Tapas``.

.. image:: ../../img/relative_ply.png
    :width: 600
    :align: center
    :alt: diagram showing the transformation of the relative orthophoto to absolute space, using the camera positions

|br| Because the footprints are most likely approximate, especially for historic datasets, this step uses
`RANSAC <https://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.ransac>`_ with a fairly large
residual threshold. The goal is to create a rough transformation of the relative orthophoto that can be used for
the gridded template matching step.

gridded template matching
--------------------------
Once the relative orthophoto has been roughly transformed to absolute space,
:py:meth:`spymicmac.register.register_ortho` find matches between the orthophoto and the reference image using
:py:meth:`spymicmac.matching.find_grid_matches`. The size of each search window is set by ``dstwin``, and the templates
(of size 121x121 pixels) are taken from a grid with spacing determined by the ``density`` parameter.

Each template and search image are first run through :py:meth:`spymicmac.image.highpass_filter`, to help minimize
radiometric differences between the two images (and maximizing the high-frequency variation). After that, the
template and search image are passed to `OpenCV matchTemplate <https://docs.opencv.org/4.5.2/d4/dc6/tutorial_py_template_matching.html>`_,
and the best match is found using normalized cross-correlation.

The correlation value of each potential match is then compared to the standard deviation of all of the correlation
values from the search image. This value (``z_corr``) is then used to filter out poor matches later on, as higher
quality matches are more likely to represent larger departures from the background correlation value:

.. image:: ../../img/correlation_match.png
    :width: 600
    :align: center
    :alt: a comparison of (a) the template, (b) the search space (with match indicated by a red plus), and (c) the correlation between the template and search image

|br|

iterative outlier removal
--------------------------
After the potential matches are found, a number of filtering steps are used to refine the results. First, any matches
where the DEM does not have a value are removed. Then, an affine transformation between the relative orthoimage
and reference orthoimage locations is estimated using RANSAC, to help remove obvious blunders.

Next, `mm3d GCPBascule <https://micmac.ensg.eu/index.php/GCPBascule>`_ is called, which transforms the camera locations
to the absolute space. The residuals for each GCP are then calculated, and outliers more than 2 normalized median
absolute deviations (NMAD) from the median residual value are discarded, and ``mm3d GCPBascule`` is called again.

This is followed by a call to `mm3d Campari <https://micmac.ensg.eu/index.php/Campari>`_ using
:py:meth:`spymicmac.micmac.campari`, and again residuals more than 2 NMAD from the median residual value are discarded.

After this, this process (``mm3d GCPBascule`` -> ``mm3d Campari`` -> outlier removal) is run up to 4 more times,
until there are no further outliers found.

final result
-------------
Once the outliers have been removed, the final GCP locations are stored in a number of files:

    - auto_gcps/AutoGCPs.xml
    - auto_gcps/AutoGCPs.txt
    - auto_gcps/AutoGCPs.shp (+ other files)
    - AutoMeasures.xml -- the GCP locations in each of the individual images

If there are still problematic GCPs, you can manually delete them from ``AutoMeasures.xml`` and re-run
``mm3d GCPBascule`` and ``mm3d Campari``.

The next step will be to run `mm3d Malt <https://micmac.ensg.eu/index.php/Malt>`_ using the ``Ori-TerrainFirstPass``
directory, to produce the absolute orthophoto and DEM.
