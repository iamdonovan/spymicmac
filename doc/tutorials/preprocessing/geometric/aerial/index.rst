aerial images
=============

Geometric pre-processing for aerial images typically involves identifying the location of 4-8 "fiducial" markers
(fixed points whose location is known) in each image.

For a general workflow of identifying the location of fiducial markers, see :doc:`this page <fiducials>`.

For example camera models and fiducial marker patterns, see :doc:`this page <cameras>`. Note that this page is not (yet)
exhaustive, and it is possible that there are other fiducial marker patterns for the camera model you are working with.

Finally, for cases where you are unable to use the original fiducial markers, or do not have a calibration report to
use, you can :doc:`estimate fiducial marker <estimating>` locations to generate the file(s) needed for geometric
resampling.

.. toctree::
   :glob:
   :maxdepth: 1
   :hidden:

   fiducials
   cameras
   estimating
