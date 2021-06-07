general workflow
=======================

The workflow and steps here are built partly from the steps outlined in the excellent
`Historical Orthoimage <https://micmac.ensg.eu/index.php/Historical_Orthoimage>`_ tutorial, written by Luc Girod. It
has been modified to include specific :doc:`../preprocessing/kh9_preprocessing`,
as well as optional :doc:`../preprocessing/image_preprocessing` steps for historic air photos.

The other main difference is the use of :doc:`../../sPyMicMac/modules/register` to find control points automatically,
using an orthorectified image and external DEM.

.. toctree::
   :maxdepth: 2

   ../preprocessing/index
   fiducials
   resample
   tapioca
   tapas
   relative
   register
   malt
