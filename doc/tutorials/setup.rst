initial setup
==============

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
marker locations into a CSV file like this:

.. code-block:: text

    gcp,im_col,im_row
    P1,-105.99,-105.994
    P2,106.002,106.004
    P3,-105.996,106.003
    P4,106.012,-105.994
    P5,-109.998,0.005
    P6,109.997,0.009
    P7,0.004,109.996
    P8,0.009,-109.995

Then, use :py:meth:`spymicmac.micmac.create_measurescamera_xml` to create ``Ori-InterneScan/MeasuresCamera.xml``:

.. code-block:: python

    from spymicmac import micmac
    micmac.create_measurescamera_xml(fn_csv, translate=True)

This will generate the following XML file:

.. code-block:: text

    <?xml version='1.0' encoding='UTF-8'?>
    <SetOfMesureAppuisFlottants>
      <MesureAppuiFlottant1Im>
        <NameIm>Glob</NameIm>
        <OneMesureAF1I>
          <NamePt>P1</NamePt>
          <PtIm>4.00800000000001 215.99</PtIm>
        </OneMesureAF1I>
        <OneMesureAF1I>
          <NamePt>P2</NamePt>
          <PtIm>216.0 3.9919999999999902</PtIm>
        </OneMesureAF1I>
        <OneMesureAF1I>
          <NamePt>P3</NamePt>
          <PtIm>4.0020000000000095 3.992999999999995</PtIm>
        </OneMesureAF1I>
        <OneMesureAF1I>
          <NamePt>P4</NamePt>
          <PtIm>216.01 215.99</PtIm>
        </OneMesureAF1I>
        <OneMesureAF1I>
          <NamePt>P5</NamePt>
          <PtIm>0.0 109.991</PtIm>
        </OneMesureAF1I>
        <OneMesureAF1I>
          <NamePt>P6</NamePt>
          <PtIm>219.995 109.987</PtIm>
        </OneMesureAF1I>
        <OneMesureAF1I>
          <NamePt>P7</NamePt>
          <PtIm>110.00200000000001 0.0</PtIm>
        </OneMesureAF1I>
        <OneMesureAF1I>
          <NamePt>P8</NamePt>
          <PtIm>110.007 219.99099999999999</PtIm>
        </OneMesureAF1I>
      </MesureAppuiFlottant1Im>
    </SetOfMesureAppuisFlottants>

.. _id_fid:

id_fiducial.txt
----------------

This is just a plain text file, with the "names" of the different fiducial marks:

.. code-block:: text

    P1
    P2
    P3

... and so on. You can easily generate this from the command line:

.. code-block:: sh

    for nn in {1..8}; do echo P$nn >> id_fiducial.txt

Alternatively, in python:

.. code-block:: python

    with open('id_fiducial.txt', 'w') as f:
        for nn in range(1, 9):
            print(f"P{nn}", file=f)

Making sure to adjust the number of fiducials based on your camera.

.. note::

    The names in the file must match the names written in ``MeasuresCamera.xml``.

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

.. _multicam:

multiple cameras
----------------

``spymicmac`` also provides some tools for setting up multiple cameras for processing at the same time.
:py:meth:`spymicmac.micmac.generate_multicam_csv` will create a CSV based on multiple input parameters. You can
either generate the CSV directly:

.. code-block:: python

    from spymicmac import micmac
    micmac.generate_multicam_csv(
        patterns = ['AR1VDYL.*', 'ARH61.*', 'AR58.*', 'AR1QM.*'],
        fn_out = 'camera_defs.csv',
        name = ['AR1VDYL', 'ARH61', 'AR58', 'AR1QM'],
        short_name = ['Wild RC8', 'Zeiss RMK A 15/23', 'Wild RC10', 'Hurd + Metrogon'],
        film_width = [220.000, 225.969, 219.995, 217.750],
        film_height = [220.009, 225.992, 219.992, 217.750],
        focal = [152.22, 152.335, 302.281, 151.1]
    )

.. note::

    Alternatively, you can run this with no arguments to generate a blank CSV template that you can edit:

    .. code-block:: python

        from spymicmac import micmac
        micmac.generate_multicam_csv()

Once you have a CSV with multiple camera definitions, run :py:meth:`spymicmac.micmac.create_localchantier_xml`, passing
the filename of your CSV with the ``cam_csv`` argument:

.. code-block:: python

    from spymicmac import micmac
    micmac.create_localchantier_xml(cam_csv='camera_defs.csv', add_sfs=True)

This will then create the ``MicMac-LocalChantierDescripteur.xml`` file with multiple cameras defined:

.. code-block:: text

    <?xml version='1.0' encoding='UTF-8'?>
    <Global>
      <ChantierDescripteur>
        <LocCamDataBase>
          <CameraEntry>
            <Name>AR1VDYL</Name>
            <SzCaptMm>220.000 220.009</SzCaptMm>
            <ShortName>Wild RC8</ShortName>
          </CameraEntry>
          <CameraEntry>
            <Name>ARH61</Name>
            <SzCaptMm>225.969 225.992</SzCaptMm>
            <ShortName>Zeiss RMK A 15/23</ShortName>
          </CameraEntry>
          <CameraEntry>
            <Name>AR58</Name>
            <SzCaptMm>219.995 219.992</SzCaptMm>
            <ShortName>Wild RC10</ShortName>
          </CameraEntry>
          <CameraEntry>
            <Name>AR1QM</Name>
            <SzCaptMm>217.7500 217.750</SzCaptMm>
            <ShortName>Hurd + Metrogon</ShortName>
          </CameraEntry>
        </LocCamDataBase>
        <KeyedNamesAssociations>
          <Calcs>
            <Arrite>1 1</Arrite>
            <Direct>
              <PatternTransform>OIS-Reech_AR1VDYL.*</PatternTransform>
              <CalcName>AR1VDYL</CalcName>
            </Direct>
          </Calcs>
          <Calcs>
            <Arrite>1 1</Arrite>
            <Direct>
              <PatternTransform>OIS-Reech_ARH61.*</PatternTransform>
              <CalcName>ARH61</CalcName>
            </Direct>
          </Calcs>
          <Calcs>
            <Arrite>1 1</Arrite>
            <Direct>
              <PatternTransform>OIS-Reech_AR58.*</PatternTransform>
              <CalcName>AR58</CalcName>
            </Direct>
          </Calcs>
          <Calcs>
            <Arrite>1 1</Arrite>
            <Direct>
              <PatternTransform>OIS-Reech_AR1QM.*</PatternTransform>
              <CalcName>AR1QM</CalcName>
            </Direct>
          </Calcs>
          <Key>NKS-Assoc-STD-CAM</Key>
        </KeyedNamesAssociations>
        <KeyedNamesAssociations>
          <Calcs>
            <Arrite>1 1</Arrite>
            <Direct>
              <PatternTransform>OIS-Reech_AR1VDYL.*</PatternTransform>
              <CalcName>152.22</CalcName>
            </Direct>
          </Calcs>
          <Calcs>
            <Arrite>1 1</Arrite>
            <Direct>
              <PatternTransform>OIS-Reech_ARH61.*</PatternTransform>
              <CalcName>152.335</CalcName>
            </Direct>
          </Calcs>
          <Calcs>
            <Arrite>1 1</Arrite>
            <Direct>
              <PatternTransform>OIS-Reech_AR58.*</PatternTransform>
              <CalcName>302.281</CalcName>
            </Direct>
          </Calcs>
          <Calcs>
            <Arrite>1 1</Arrite>
            <Direct>
              <PatternTransform>OIS-Reech_AR1QM.*</PatternTransform>
              <CalcName>151.1</CalcName>
            </Direct>
          </Calcs>
          <Key>NKS-Assoc-STD-FOC</Key>
        </KeyedNamesAssociations>
      </ChantierDescripteur>
    </Global>

Tools such as ``Tapas`` will then calibrate the intrinsic parameters for each camera separately.

.. note::

    If you are processing images from multiple cameras simultaneously, you should still do the geometric processing
    for the images from each camera separately. Once the images have been resampled, you can proceed with steps like
    ``Tapioca`` and ``Tapas`` with all of the resampled images in the same directory.
