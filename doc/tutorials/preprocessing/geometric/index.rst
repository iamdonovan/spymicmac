geometric pre-processing
========================

Whether working with scanned :doc:`declassified <declass/index>` or :doc:`aerial <aerial/index>`, images, some form
of geometric pre-processing is typically required. Most of the time, this involves identifying the location of
known/fixed points in each image, followed by some form of :doc:`resampling <resampling>`, either a geometric
transformation or cropping, so that the images have a consistent geometry.

.. toctree::
   :glob:
   :maxdepth: 1
   :hidden:

   declass/index
   aerial/index
   resampling
