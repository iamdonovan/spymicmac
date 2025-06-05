relative geometry
=================

After pre-processing the images, the next step is to calibrate the camera model(s) by estimating the intrinsic
parameters for each camera, and estimating a relative orientation of the different images.

- First, :doc:`mm3d Tapioca <tapioca>` is used to compute tie points between the images
- Then, :doc:`mm3d Tapas <tapas>` is used to calibrate the camera parameters and compute the initial relative
  orientation of the images. ``mm3d Martini`` can be used first, to help initialize the orientation and improve the
  overall ``Tapas`` result.
- Finally, :doc:`mm3d Malt <malt>` is used to compute the relative DEM and orthophotos. ``mm3d Tawny`` can also be used
  to compute the relative orthomosaic.

.. toctree::
   :glob:
   :maxdepth: 1
   :hidden:

   tapioca
   tapas
   malt
