micmac processing files
================================
There are a number of files that you will need to process your images using MicMac, which are detailed below. If you
are working with KH-9 images, :doc:`../../spymicmac/scripts/generate_micmac_measures` will automatically create the
necessary files.

You can use the provided :doc:`../../examples/index` as a starting point, though note that you will probably need to
make some modifications as detailed below.

MicMac-LocalChantierDescripteur.xml
------------------------------------
You will need to set the image size (in mm) by changing the values in ``SzCaptMm`` under ``CameraEntry`` (line 8 of the
example). You will also need to set the image matching pattern (``PatternTransform``, line 29 of the example) and focal
length in mm (``CalcName``, line 30) for each set of images - if, for example, different focal lengths were used.

To improve tie point density and matching, especially in low-contrast images, you can try changing ``CalcName`` (line 18)
to ``SFS``, and ``Key`` (line 21) to ``NKS-Assoc-SFS``:
::

    <KeyedNamesAssociations>
        <Calcs>
             <Arrite>  1 1 </Arrite>
             <Direct>
                   <PatternTransform> .*  </PatternTransform> <!-- Regular expressions of the group of images with the following camera model -->
                   <CalcName> SFS </CalcName> <!-- Name of the camera for these images -->
             </Direct>
         </Calcs>
         <Key> NKS-Assoc-SFS </Key>
    </KeyedNamesAssociations>

MeasuresCamera.xml
-------------------
.. note::
    This file must be placed in a directory called ``Ori-InterneScan``

Ideally, you will have a camera calibration report, that will tell you the location
of the different fiducial markers in the image geometry. Note that using ``ReSampFid`` **requires** a file,
``Ori-InterneScan/MeasuresCamera.xml``, that tells MicMac what the location of each fiducial mark is.

.. note::
    The image coordinates are defined with the origin in the upper left corner, rather than the center
    of the image used by most calibration files. You can translate from one system to the other with the following:
    ::

        xp = x - min(x)
        yp = (-y) - min(-y)

id_fiducial.txt
----------------
This is just a plain text file, with the "names" of the different fiducial marks:
::

    P1
    P2
    P3

... and so on.

.. note::
    The names in the file should match the names written in ``MeasuresCamera.xml``.

file structure
----------------
Before starting, your file structure should look something like this:
::

    project
    ├── id_fiducial.txt
    ├── Img1.tif
    ├── Img2.tif
    ...
    ├── MicMac-LocalChantierDescripteur.xml
    ├── Ori-InterneScan
    │   └── MeasuresCamera.xml

Once you have this set up, you can work on the preprocessing steps.
