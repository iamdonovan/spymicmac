initial setup files
===================

There are a number of files that you will need to process your images using MicMac, which are detailed below.

MicMac-LocalChantierDescripteur.xml
------------------------------------

.. note::

    The following explains the use of :py:meth:`spymicmac.micmac.create_localchantier_xml` or
    :doc:`../../../../spymicmac/scripts/create_localchantier_xml` to generate this file.

    Alternatively, you can download and edit the example file linked below.

To create the ``MicMac-LocalChantierDescripteur.xml`` file, you need the following information:

- **name** - This is the name that is used in the ``AutoCal.xml`` file that describes the camera parameters, so it
  should be as short as possible and without spaces.
- **short_name** - a (longer) description of the camera
- **film_size** - the film size [width, height] in mm. This should be the size of the image frame after resampling
  using the fiducial markers - for example, ``xmax - xmin``, ``ymax - ymin`` for the maximum/minimum x,y position of the
  fiducial markers.
- **pattern** - the file matching pattern to use for the images. Use ``'.*'`` to match all images in the directory.
- **focal** - the nominal focal length for the camera, in mm.
- **add_sfs** - whether to use SFS matching to help find tie points in low-contrast images.
- **cam_csv** - an optional name of a CSV file that defines multiple cameras, if attempting to process images from
  multiple cameras at once (see more information below).

MeasuresCamera.xml
-------------------
.. warning::

    This file must be placed in a directory called ``Ori-InterneScan`` in order for MicMac to be able to use it.

Ideally, you will have a camera calibration report, that will tell you the location
of the different fiducial markers in the image geometry. Note that using ``ReSampFid`` **requires** a file,
``Ori-InterneScan/MeasuresCamera.xml``, that tells MicMac what the location of each fiducial mark is.

.. note::
    The image coordinates are defined with the origin in the upper left corner, rather than the center
    of the image used by most calibration files. You can translate from one system to the other with the following:

    .. code-block:: text

        xp = x - min(x)
        yp = (-y) - min(-y)

If you do not have a calibration report for your particular camera, you can have a look at some
:doc:`preprocessing/geometric/aerial/cameras` for approximate locations of fiducial markers.

Rather than editing the ``MeasuresCamera.xml`` file with the fiducial marker locations, you can also put the fiducial
marker locations into a CSV file, then use :py:meth:`spymicmac.micmac.create_measurescamera_xml` to create
``Ori-InterneScan/MeasuresCamera.xml``.


id_fiducial.txt
----------------

This is just a plain text file, with the "names" of the different fiducial marks:

.. code-block:: text

    P1
    P2
    P3

... and so on.

.. note::

    The names in the file should match the names written in ``MeasuresCamera.xml``.

file structure
----------------
Before starting, your file structure should look something like this:

.. code-block:: text

    project
    ├── id_fiducial.txt
    ├── Img1.tif
    ├── Img2.tif
    ...
    ├── MicMac-LocalChantierDescripteur.xml
    ├── Ori-InterneScan
    │   └── MeasuresCamera.xml

Once you have this set up, you can work on the preprocessing steps.

example micmac files
--------------------

Below, you can download some example files to get started with your MicMac project. Note you will
most likely need to modify them as explained above.

- :download:`id_fiducial.txt <../examples/id_fiducial.txt>`
- :download:`MeasuresCamera.xml <../examples/MeasuresCamera.xml>`
- :download:`MicMac-LocalChantierDescripteur.xml <../examples/MicMac-LocalChantierDescripteur.xml>`

multiple cameras
----------------

``spymicmac`` also provides some tools for setting up multiple cameras for processing at the same time.
:py:meth:`spymicmac.micmac.generate_multicam_csv` will create a CSV based on multiple input parameters. This
filename can then be passed to :py:meth:`spymicmac.micmac.create_localchantier_xml` to create
``MicMac-LocalChantierDescripteur.xml`` with multiple cameras defined. Tools such as ``Tapas`` will then
calibrate the intrinsic parameters for each camera separately.

.. note::

    Geometric processing for the images from each camera should be done separately. Once the images have been
    resampled to a common format, you can proceed with steps like ``Tapioca`` and ``Tapas`` with all of the
    resampled images in the same directory.
