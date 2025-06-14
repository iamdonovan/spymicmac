example camera models
=======================

This page contains example diagrams for different camera models that have acquired historical aerial photos, including
what their fiducial marks look like, and the approximate coordinates of the fiducial marks that you can use to populate
the ``MeasuresCamera.xml`` file.

.. note::

    If you have a calibration report that corresponds to your specific images, **you should use that instead**.

    The information provided here is for those cases where a calibration report does not exist, or has been lost to time.

Fairchild Cameras
-----------------

"wing-type" fiducials
^^^^^^^^^^^^^^^^^^^^^

These cameras tend to have 4 mid-side fiducial markers with varying patterns. Some examples of these patterns are shown
here:



.. _fairchild k17:

F224, K17B (Metrogon lens)
""""""""""""""""""""""""""

.. image:: img/fairchild.png
    :width: 400
    :align: center
    :alt: a diagram of a Fairchild F224 K17B with fiducial markers labeled

|br|

+----+-----+-----+
|    | x   | y   |
+----+-----+-----+
| P1 | 1   | 113 |
+----+-----+-----+
| P2 | 225 | 113 |
+----+-----+-----+
| P3 | 113 | 1   |
+----+-----+-----+
| P4 | 113 | 225 |
+----+-----+-----+

.. note::

    These are the approximate coordinates for the corners of the fiducial marker, as these tend to be more stable
    than the tips of the "wings"

notch-type fiducials
^^^^^^^^^^^^^^^^^^^^

.. _fairchild t11d:

KC-1, KC-1B, T-11 (Metrogon Lens)
"""""""""""""""""""""""""""""""""

.. image:: img/fairchild_kc1b.png
    :width: 400
    :align: center
    :alt: a diagram of a Fairchild KC-1B with fiducial markers labeled

|br|

+----+-----+-------+
|    | x   | y     |
+----+-----+-------+
| P1 | 0   | 117.5 |
+----+-----+-------+
| P2 | 238 | 117.5 |
+----+-----+-------+
| P3 | 120 | 0     |
+----+-----+-------+
| P4 | 120 | 235   |
+----+-----+-------+

check-type fiducials
^^^^^^^^^^^^^^^^^^^^

.. _fairchild t11s:

T-11
""""

These cameras may have different fiducial patterns - for example, they may have a checked pattern like this:

.. image:: img/fairchild_t11.png
    :width: 400
    :align: center
    :alt: a diagram of a Fairchild T-11 with fiducial markers labeled

+----+-----+-------+
|    | x   | y     |
+----+-----+-------+
| P1 | 0   | 117.5 |
+----+-----+-------+
| P2 | 238 | 117.5 |
+----+-----+-------+
| P3 | 120 | 0     |
+----+-----+-------+
| P4 | 120 | 235   |
+----+-----+-------+


Wild Heerbrugg Cameras
----------------------

.. _wild rc5:

Wild RC5A, RC5, RC8 (Aviogon Lens)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These cameras tend to have four crosshair fidcuial marks in the corners, set inside of a rounded frame:

.. image:: img/wild_rc5.png
    :width: 400
    :align: center
    :alt: a diagram of a Wild RC5-type with fiducial markers labeled

|br|

+----+-----+-----+
|    | x   | y   |
+----+-----+-----+
| P1 | 212 | 0   |
+----+-----+-----+
| P2 | 212 | 212 |
+----+-----+-----+
| P3 | 0   | 212 |
+----+-----+-----+
| P4 | 0   | 0   |
+----+-----+-----+

.. _wild rc10:

Wild RC8, RC10 (Aviogon Lens)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The RC8 and RC10 camera might also come in a pattern with 8 fiducial markers: 4 corner markers as on the RC5-type, and
4 mid-side markers:

.. image:: img/wild_rc10.png
    :width: 400
    :align: center
    :alt: a diagram of a Wild RC10 with fiducial markers labeled

|br|

+----+-----+-----+
|    | x   | y   |
+----+-----+-----+
| P1 | 4   | 216 |
+----+-----+-----+
| P2 | 216 | 4   |
+----+-----+-----+
| P3 | 4   | 4   |
+----+-----+-----+
| P4 | 216 | 216 |
+----+-----+-----+
| P5 | 0   | 110 |
+----+-----+-----+
| P6 | 220 | 110 |
+----+-----+-----+
| P7 | 110 | 0   |
+----+-----+-----+
| P8 | 110 | 220 |
+----+-----+-----+

.. _zeiss midside:

Zeiss Cameras
-------------

Zeiss RMK 15/23 (Pleogon Lens)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: img/zeiss.png
    :width: 400
    :align: center
    :alt: a diagram of a Zeiss RMK with fiducial markers labeled

|br|

+----+-----+-----+
|    | x   | y   |
+----+-----+-----+
| P1 | 0   | 113 |
+----+-----+-----+
| P2 | 226 | 113 |
+----+-----+-----+
| P3 | 113 | 0   |
+----+-----+-----+
| P4 | 113 | 226 |
+----+-----+-----+

.. note::

    The coordinates above correspond to the center of the small dot near the tip of the fiducial marker.

.. _zeiss corner:

Zeiss RMK A 15/23 (Pleogon Lens)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Later versions of the Zeiss RMK camera used 8 fiducial markers: 4 mid-side markers, and 4 corner markers.

.. image:: img/zeiss_rmka.png
    :width: 400
    :align: center
    :alt: a diagram of a Zeiss RMK with fiducial markers labeled

|br|

+----+-------+-------+
|    | x     | y     |
+----+-------+-------+
| P1 | 9     | 217   |
+----+-------+-------+
| P2 | 217   | 9     |
+----+-------+-------+
| P3 | 9     | 9     |
+----+-------+-------+
| P4 | 217   | 217   |
+----+-------+-------+
| P5 | 0     | 113   |
+----+-------+-------+
| P6 | 226   | 113   |
+----+-------+-------+
| P7 | 113   | 0     |
+----+-------+-------+
| P8 | 113   | 226   |
+----+-------+-------+

.. note::

    The coordinates for P5-P8 above correspond to the center of the small dot near the tip of the fiducial marker.

