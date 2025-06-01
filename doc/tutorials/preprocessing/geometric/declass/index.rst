declassified images
===================

Scanned declassified images, provided by the USGS through `EarthExplorer <https://earthexplorer.usgs.gov/>`__,
typically need to be :doc:`joined <joining>` together, owing to the large film size. After joining, they need to be
resampled to a common format, in exactly the same fashion as scanned aerial photographs.

However, because information about the fiducial markers for most of these cameras is either unknown, still
classified, or otherwise unavailable, we typically use alternative methods for resampling these images:

- for KH-9 Hexagon mapping camera images, we use the grid of :doc:`réseau markers <reseau>` to re-sample the images to
  a common geometry. This also has the benefit of potentially correcting film distortion, as the réseau markers were
  on a glass plate exposed to the film, and should therefore be fixed points.
- for :doc:`panoramic <panoramic>` cameras such as KH-4/A/B ("Corona") or the KH-9 panoramic camera, which did not have
  réseau marks or fiducial markers, we use the image frame to re-sample the images to as common a format as possible.

.. toctree::
   :glob:
   :maxdepth: 1
   :hidden:

   joining
   reseau
   panoramic
