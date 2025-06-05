tutorials
=========

Here, you will find information about the workflow for using spymicmac to process KH-9 mapping camera and historic air
photos.

The basic workflow looks like this:

.. image:: ../img/workflow.png
    :width: 720
    :align: center
    :alt: a flow diagram showing the steps involved in spymicmac processing

|br|

This can be divided into four main steps: initial setup, pre-processing, relative geometry, and absolute geometry:

- :doc:`initial setup <setup>` covers the different files that are necessary for processing using MicMac
- :doc:`pre-processing <preprocessing/index>` covers all of the steps taken to get the images ready for processing, and
  includes both :doc:`geometric <preprocessing/geometric/index>` pre-processing and
  :doc:`radiometric <preprocessing/geometric/index>` pre-processing.
- :doc:`relative <relative/index>` geometry covers all of the steps taken to estimate the intrinsic camera parameters
  (e.g., principal point, radial distortion, affine/decentric correction) and the relative external orientation of
  the images.
- :doc:`absolute <absolute/index>` geometry covers all of the steps taken to process the images in an absolute
  ("real-world") geometry: registration and bundle block adjustment (BBA), dense matching and orthorectification.

If you don't want to read more about the different processing steps, there is also a :doc:`../tealdeer` which runs through
the different tools and commands used in the workflow with minimal explanation. :)

.. toctree::
   :glob:
   :maxdepth: 2
   :hidden:

   setup
   preprocessing/index
   relative/index
   absolute/index
