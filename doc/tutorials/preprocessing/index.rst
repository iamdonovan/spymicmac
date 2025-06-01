image pre-processing
=======================

This section describes a few of the different pre-processing tools available for historical scanned images, which
are broken into two categories:

- :doc:`geometric <geometric/index>` tools, for re-sampling images so that they have the same shape. This makes it easier to
  estimate the intrinsic camera parameters (e.g., principal point, radial distortion, affine/decentric parameters).
- :doc:`radiometric <radiometric>` tools, which include steps such as de-striping, contrast enhancement, and
  de-noising scanned images. These steps can improve the end result, as well as the performance of the matching
  algorithm(s) used to extract elevation.

.. toctree::
   :maxdepth: 2
   :hidden:

   geometric/index
   radiometric
