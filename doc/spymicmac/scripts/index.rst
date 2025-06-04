command line tools
===================

``spymicmac`` comes with the following command line tools, organized roughly by what they are used for:

setting up micmac files
-----------------------

- :doc:`generate_measures_files`: for creating the various input files needed for processing KH-9 Hexagon mapping
  camera images.
- :doc:`create_localchantier_xml`: for creating a MicMac-LocalChantierDescripteur.xml file.
- :doc:`write_micmac_xml`: creates an XML file for a geotiff. Used for preparing a seed DEM to pass to ``mm3d Malt``,
  for example.

pre-processing
--------------

- :doc:`preprocess_kh9_mc`: for pre-processing KH-9 Hexagon mapping camera images from .tgz files to a
  relative orientation
- :doc:`preprocess_pan`: for pre-processing panoramic camera images from .tgz files to cropped images.

individual steps
^^^^^^^^^^^^^^^^

These tools do individual steps as part of pre-processing:

- :doc:`join_hexagon`: for joining scanned parts of an image together.
- :doc:`find_reseau_grid`: for finding réseau marks in a KH-9 Hexagon mapping camera image.
- :doc:`remove_crosses`: for removing réseau marks in a KH-9 Hexagon mapping camera image.
- :doc:`resample_hexagon`: for resampling a KH-9 Hexagon mapping camera image using the réseau marks found previously.

image registration and gcps
---------------------------

The following tools perform steps related to automatically finding GCPs:

- :doc:`register_relative`: for automatically finding GCPs, given a relative DEM/orthophoto and a reference
  DEM/orthoimage.
- :doc:`download_cop30_vrt`: for downloading and creating a VRT of Copernicus 30m DSM tiles.
- :doc:`block_orientation`: combines GCPs and Measures files from multiple blocks into a single orientation, then
  runs :py:meth:`spymicmac.micmac.iterate_campari` to refine the orientation + camera parameters.
- :doc:`combine_auto_measures`: used for combining the outputs of ``mm3d XYZ2Im`` into a single xml file.
- :doc:`remove_measures`: used for removing named GCPs from a Measures xml file - for example, if the GCP is badly
  placed.


image processing
----------------

- :doc:`balance_images`: for applying CLAHE to all images in the current directory.


post-processing
---------------

:doc:`post_process_micmac`: for masking and georeferencing the DEM and orthophoto(s) or orthomosaic.
:doc:`mosaic_micmac_tiles`: for mosaicking the tiled outputs from Malt.


.. toctree::
   :glob:
   :titlesonly:
   :hidden:

   *