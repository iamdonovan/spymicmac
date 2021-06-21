generate_micmac_measures.py
=================================

``generate_micmac_measures.py`` calls :py:meth:`spymicmac.micmac.generate_micmac_measures` to create ``id_fiducial.txt``,
``MeasuresCamera.xml``, and ``Tmp-SL-Glob.xml`` files for KH-9 Hexagon images.

.. note::
    To use the ``Tmp-SL-Glob.xml`` file, copy it into the ``Tmp-SaisieAppuis`` directory with the name of the image
    appended; e.g.,:
    ::

        cp Tmp-SL-Glob.xml Tmp-SaisieAppuis/Tmp-SL-Glob-MeasuresIm-DZB1215-500425L002001.tif.xml

    You will also need MeasuresIm-DZB1215-500425L002001.tif.xml (or whatever the image name is) to exist in the directory.


.. argparse::
   :filename: ../bin/generate_micmac_measures.py
   :func: _argparser
   :prog: generate_micmac_measures.py
