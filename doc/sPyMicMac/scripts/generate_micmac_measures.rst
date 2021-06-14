generate_micmac_measures.py
=================================

``generate_micmac_measures.py`` calls :py:meth:`sPyMicMac.micmac.generate_micmac_measures` to create ``id_fiducial.txt``,
``MeasuresCamera.xml``, and ``Tmp-SL-Glob.xml`` files (that can be used to visualize fiducial mark locations
using ``SaisieAppuisInit``) for KH-9 Hexagon images.

.. argparse::
   :filename: ../bin/generate_micmac_measures.py
   :func: _argparser
   :prog: generate_micmac_measures.py
