generate_measures_files
=========================

``generate_measures_files`` calls :py:meth:`spymicmac.micmac.generate_measures_files` to create:

- ``id_fiducial.txt``,
- ``MeasuresCamera.xml``, and
- ``Tmp-SL-Glob.xml``

files for KH-9 Hexagon mapping camera images.

.. warning::

    Once created, you should move ``MeasuresCamera.xml`` into the ``Ori-InterneScan`` directory so that it can be
    found by ``mm3d ReSampFid`` or :py:meth:`spymicmac.resample.resample_hex`.

.. note::
    To use the ``Tmp-SL-Glob.xml`` file, copy it into the ``Tmp-SaisieAppuis`` directory with the name of the image
    appended; e.g.,:

    .. code-block:: sh

        cp Tmp-SL-Glob.xml Tmp-SaisieAppuis/Tmp-SL-Glob-MeasuresIm-DZB1215-500425L002001.tif.xml

    You will also need MeasuresIm-DZB1215-500425L002001.tif-S2D.xml (or whatever the image name is) to exist in the directory:

    .. code-block:: sh

        cp Ori-InterneScan/MeasuresIm-DZB1215-500425L002001.tif.xml MeasuresIm-DZB1215-500425L002001.tif-S2D.xml

    and, you should remove any temporary files from ``Tmp-SaisieAppuis``:

    .. code-block:: sh

        rm Tmp-SaisieAppuis/Tmp-SL-Im-MeasuresIm-MeasuresIm-DZB1215-500425L002001.tif.xml.*


.. argparse::
   :module: spymicmac.tools.generate_measures_files
   :func: _argparser
   :prog: generate_measures_files
