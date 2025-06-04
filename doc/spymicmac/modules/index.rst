modules
=========

``spymicmac`` contains the following modules, organized roughly by what they are used for.

photogrammetry software
------------------------

The following modules are used for interfacing with different photogrammetric processing software. In spite of the
name, ``spymicmac`` can also be used to help prepare files for using
`Ames Stereo Pipeline <https://stereopipeline.readthedocs.io/en/latest/introduction.html>`__, especially for processing
panoramic camera images.

- :py:mod:`spymicmac.micmac`: for running ``mm3d`` commands and reading/writing XML files used by MicMac.
- :py:mod:`spymicmac.register`: for finding GCPs to register images to reference datasets
- :py:mod:`spymicmac.orientation`: for reading MicMac orientation directories and transforming coordinate systems
- :py:mod:`spymicmac.asp`: for preparing camera files for ASP.


declassified images
-------------------

The following modules are primarily used for working with declassified datasets:

- :py:mod:`spymicmac.preprocessing`: for pre-processing KH-9 Hexagon mapping camera and declassified panoramic camera
  files.
- :py:mod:`spymicmac.declass`: for various metadata information about declassified datasets


reference data
--------------

The following modules are used for preparing/accessing reference datasets:

- :py:mod:`spymicmac.data`: for accessing reference datasets, including USGS footprints, Copernicus 30m DSM tiles,
  and the Polar Geospatial Center DEMs, ArcticDEM and REMA.
- :py:mod:`spymicmac.ee_tools`: (in progress) for accessing reference images and datasets using Google Earth Engine.


image processing
----------------

The following modules are primarily used for image processing, resampling, and feature matching:

- :py:mod:`spymicmac.image`: for radiometric processing and image mosaicking
- :py:mod:`spymicmac.resample`: for resampling images
- :py:mod:`spymicmac.matching`: for matching templates or keypoints between images



.. toctree::
   :glob:
   :titlesonly:
   :hidden:

   *