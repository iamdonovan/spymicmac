absolute geometry
=================

After computing the relative orientation of the images and the relative DEM and orthomosaic, the next step is to
register the images to reference DEM or orthoimage, compute the bundle block adjustment, and compute the absolute
DEM and orthophoto.

- First, :py:meth:`spymicmac.register.register_relative` is used to :doc:`register <register>` the images to the
  reference DEM or orthoimage and compute the bundle block adjustment to refine the camera parameters and orientation.
- If using a large number of images, it might be preferable to process and register the relative DEM/orthomosaic in
  smaller blocks. If this is done, the orientation of each :doc:`block <blocks>` should be re-combined before computing
  the absolute DEM and orthophotos.
- Finally, :doc:`mm3d Malt <malt>` is again used, this time to compute the absolute DEM and orthphotos. As before,
  ``mm3d Tawny`` can be used to compute the orthomosaic.

.. toctree::
   :glob:
   :maxdepth: 1
   :hidden:

   register
   blocks
   malt
