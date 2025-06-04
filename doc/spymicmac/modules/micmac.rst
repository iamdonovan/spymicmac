spymicmac.micmac
==================

.. automodule:: spymicmac.micmac


setup
-----

These methods are used to set up the different files needed to run MicMac.

.. autofunction:: spymicmac.micmac.create_localchantier_xml

.. autofunction:: spymicmac.micmac.generate_multicam_csv

.. autofunction:: spymicmac.micmac.create_measurescamera_xml

.. autofunction:: spymicmac.micmac.generate_measures_files

.. autofunction:: spymicmac.micmac.init_autocal


micmac interface
----------------

These methods are used to interface with their corresponding ``mm3d`` commands, based on where they fall in the
:doc:`../../tutorials/index` processing workflow.

- **pre-processing**: the steps taken to get the images ready for processing
- **relative**: the steps taken to estimate the intrinsic camera parameters (e.g., principal point, radial distortion,
  affine/decentric correction) and the relative external orientation of the images
- **absolute**: the steps taken to process the images in an absolute (“real-world”) geometry: registration and bundle
  block adjustment (BBA), dense matching and orthorectification
- **post-processing**: the steps taken to correct or adjust the results after processing

pre-processing
^^^^^^^^^^^^^^

.. autofunction:: spymicmac.micmac.batch_saisie_fids

.. autofunction:: spymicmac.micmac.estimate_measures_camera


finding tie points
^^^^^^^^^^^^^^^^^^

.. autofunction:: spymicmac.micmac.tapioca

.. autofunction:: spymicmac.micmac.write_neighbour_images

.. autofunction:: spymicmac.micmac.pairs_from_footprints


relative geometry
^^^^^^^^^^^^^^^^^

The following methods are used for working with images in a relative geometry:

.. autofunction:: spymicmac.micmac.martini

.. autofunction:: spymicmac.micmac.tapas


absolute geometry
^^^^^^^^^^^^^^^^^

.. autofunction:: spymicmac.micmac.bascule

.. autofunction:: spymicmac.micmac.campari

.. autofunction:: spymicmac.micmac.iterate_campari

.. autofunction:: spymicmac.micmac.checkpoints

.. autofunction:: spymicmac.micmac.get_autogcp_locations


post-processing
^^^^^^^^^^^^^^^

.. autofunction:: spymicmac.micmac.post_process

.. autofunction:: spymicmac.micmac.mosaic_micmac_tiles

.. autofunction:: spymicmac.micmac.banana


general
^^^^^^^

.. autofunction:: spymicmac.micmac.apericloud

.. autofunction:: spymicmac.micmac.malt

.. autofunction:: spymicmac.micmac.tawny


homologous points
------------------

The following methods are used for filtering homologous point files, or for finding connected blocks of images
based on their shared tie points.

.. autofunction:: spymicmac.micmac.clean_homol

.. autofunction:: spymicmac.micmac.find_empty_homol

.. autofunction:: spymicmac.micmac.find_connected_blocks

.. autofunction:: spymicmac.micmac.separate_blocks

.. autofunction:: spymicmac.micmac.pairs_from_homol


measures files
--------------

These methods are used for working with various Measures XML files.

.. autofunction:: spymicmac.micmac.parse_im_meas

.. autofunction:: spymicmac.micmac.remove_measure

.. autofunction:: spymicmac.micmac.remove_worst_mesures


camera xml files
----------------

.. autofunction:: spymicmac.micmac.load_cam_xml

.. autofunction:: spymicmac.micmac.write_cam_xml


version control
---------------

In case you wish to keep track of changes to the various XML files created by MicMac, especially the camera
calibration and orientation files, the following method will initialize an empty ``git`` repository in the current
directory, and populate the ``.gitignore`` file with common patterns to ignore (e.g., ``*.tif``).

.. autofunction:: spymicmac.micmac.init_git