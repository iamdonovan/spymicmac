example fiducial marker patterns
================================

This page contains example fiducial marker patterns for different cameras that have acquired historical aerial photos,
including what their fiducial marks look like, and the approximate coordinates of the fiducial marks that you can use
to populate the ``MeasuresCamera.xml`` file.

.. note::

    The fiducial markers in the following examples are labeled as follows, with the assumption that the data strip /
    direction of flight marker is on the **left side** of the image:

    - **P1**: the lower left corner of the image
    - **P2**: the upper right corner of the image
    - **P3**: the upper left corner of the image
    - **P4**: the lower right corner of the image
    - **P5**: the middle of the left edge image
    - **P6**: the middle of the right edge of the image
    - **P7**: the middle of the top edge of the image
    - **P8**: the middle of the bottom edge of the image

a note on identification
------------------------

Aside from a single 1960 article\ [1]_, I have not managed to find much solid information about different fiducial
marker patterns used by different camera manufacturers in different cameras. Some patterns are relatively distinct,
but others are both very similar and widely used, making a definitive identification based on the fiducial pattern
alone challenging, especially for earlier cameras that did not expose this ancillary information on the image negative.

A good deal of the information on this page has been pieced together from looking at a very large number of images
and calibration reports downloaded from `EarthExplorer <https://earthexplorer.usgs.gov/>`__ and the USGS
`Camera Calibration Report Database <https://calval.cr.usgs.gov/cameracal/reports.html>`__. This guide is by no means
exhaustive, though it is built on a fairly large dataset of example images.

.. note::

    If you have a calibration report that corresponds to your specific images, **you should use that instead**.

    The information provided here is for those cases where a calibration report does not exist, or has been lost to time.

.. tip::

    For images acquired by the US Government or other agencies and hosted on `EarthExplorer <https://earthexplorer.usgs.gov/>`__,
    you may be able to find a calibration report by using the lens or camera serial number, if this was recorded on the
    camera's data strip (or, sometimes, written around the edges of the frame).

    With the serial number(s) and/or other identifying information (e.g., make/model), you can use the USGS
    `Camera Calibration Report Database <https://calval.cr.usgs.gov/cameracal/reports.html>`__ to search for the
    calibration report for your camera (or a similar camera).


"wing-type" fiducials
---------------------

A number of different camera manufacturers have used four "wing-shaped" mid-side fiducial markers of varying designs:

- Aero Service Corporation (also Aero/View cameras, or just "Aero")
- Airagon
- Mark Hurd Aerial Surveys, Inc.
- Park Aerial Surveys, Inc.
- Fairchild
- Wes Smith

Many of the cameras made/used by Mark Hurd Aerial Surveys, Inc., Park Aerial Surveys, Inc., or similar companies
do not appear to have model names - rather, they are identified only by a number. The lenses used in these cameras
tend to be made by manufacturers such as Bausch and Lomb (e.g., the Metrogon Lens), C. P. Goerz (e.g., the Aerotar
lens), or Zeiss (e.g., the Pleogon lens).

.. tip::

    Because:

    - the measured separation between P5-P6 and P7-P8 is somewhat variable for these images, even for the same camera
      manufacturer (and sometimes even the same camera...); and
    - the majority of these cameras were used before coordinates were regularly included in calibration reports,

    it is probably a better idea to use :py:meth:`estimate_measures_camera <spymicmac.micmac.estimate_measures_camera>`
    with a large sample of these images (and the known scanning resolution), rather than the approximate coordinates
    given below.

hollow-type wing
^^^^^^^^^^^^^^^^^

Cameras with this style of fiducial marker (four mid-side 'wing' markers that are at least semi-transparent) seem to
have been made by Mark Hurd Aerial Surveys, Inc. They do not have any particular model name, and are recorded with
single- or double-digit serial numbers in most calibration reports on file. They used a few different lens types,
mostly Bausch and Lomb Metrogon with a 6" (152 mm) focal length.

.. image:: img/hollow_wing.png
    :width: 500
    :align: center
    :alt: an example of a Mark Hurd camera with four transparent mid-side wing fiducial markers


flat wing (small)
^^^^^^^^^^^^^^^^^

This style is perhaps the most common of the wing-type (at least by volume in EarthExplorer), and seems to have been
used by a wide range of manufacturers:

- Fairchild (K-3B "modified", K-17, K-17 "modified", K-17B "modified")\ [2]_
- Aero Services (Aero, Aero/View)
- Airagon
- Western Aerial Contractors
- Park (Twinplex, unnamed)

Some, but not all, of these used a direction of flight indicator similar to the one shown below on P5. Some models
also had early data strip styles, with both serial number and calibrated focal length included in the image (the
example below, from a Park camera, has these on the left side of the frame).

.. image:: img/small_flat_wing.png
    :width: 500
    :align: center
    :alt: an example of an aerial camera with four "small flat" mid-side wing fiducial markers


flat wing (large)
^^^^^^^^^^^^^^^^^

A larger version of the "flat" wing was used in at least the Fairchild F-224 and K-17A cameras. In the lower-right
corner of the frame, there is sometimes an indicator with the focal length (when not overwritten by other markings
on the print, that is).

Note also the diamond-shaped marker on P5, which I am taking to be the "direction of flight marker" as the arrow-shaped
notches on both P5 and P6 point in opposing directions.

.. image:: img/large_flat_wing.png
    :width: 500
    :align: center
    :alt: an example of an aerial camera with four "large flat" mid-side wing fiducial markers

flat wing (tall)
^^^^^^^^^^^^^^^^^

A "tall" version of the "flat" wing was used in primarily Aero Services (Aero, Aero/View) cameras. Some versions of
these cameras have a small marker that includes the camera number on it, usually just below P5. Note also the
diamond-shaped marker on P5, which I am taking to be the "direction of flight marker".

.. image:: img/tall_flat_wing.png
    :width: 500
    :align: center
    :alt: an example of an aerial camera with four "tall flat" mid-side wing fiducial markers


flat wing (curved)
^^^^^^^^^^^^^^^^^^

This design was used primarily in cameras manufactured by Wes Smith (or just "Smith"). Note that marker P5 in the
example below has a "direction of flight" indicator, meant to point in the direction that the plane was flying (in
this case, towards the right of the frame):

.. image:: img/curved_flat_wing.png
    :width: 500
    :align: center
    :alt: a diagram of a camera with four "curved wing" mid-side fiducial markers


curved wing
^^^^^^^^^^^

This design was used in at least some versions of the K-17 "modified" camera. Note that three of the markers have
notches on either side, or are slightly inset from the rest of the frame - in the absence of any other identifying
marks, I am taking P5 to be the odd one out.

.. image:: img/curved_wing.png
    :width: 500
    :align: center
    :alt: a diagram of a camera with four "curved wing" mid-side fiducial markers


.. _fairchild k17:

sharp wing
^^^^^^^^^^

This design was used in at least the Fairchild K-17 and some "modified" versions of the K-17 camera. As with the
F-224 and K-17B cameras, there is a strip indicating the calibrated focal length in the lower right-hand corner of the
frame, as well as a "direction of flight" marker on P5:

.. image:: img/sharp_wing.png
    :width: 500
    :align: center
    :alt: a diagram of a camera with four "sharp wing" mid-side fiducial markers


Fairchild Cameras
-----------------

arrow-type fiducials
^^^^^^^^^^^^^^^^^^^^^

.. caution::

    I have only been able to identify a single example of this type of fiducial marker connected to a calibration
    report, and to notes made on the edges of the frame of an image from a separate survey.

    The reported camera in both cases was a Fairchild K-17, but I am hesitant to definitively state this. As with the
    wing-type examples shown above, I recommend using
    :py:meth:`estimate_measures_camera <spymicmac.micmac.estimate_measures_camera>` in the absence of any definitive
    marker coordinates.

These cameras have four mid-side fiducial markers. Both P5 and P6 are "arrows" that point in the direction of flight.
P5 is the tip of the arrow-shaped notch that is level with the frame, while P6 is the tip of the arrow-shaped notch
that is cut into the frame.

Many prints (or scans) of images of this type will cut off the point of the notch on P7 and P8, so it is probably
better to use the point where the vertical edge of the notch instersects the image frame, rather than the "point" of
the notch.

.. image:: img/arrow.png
    :width: 500
    :align: center
    :alt: a diagram of a fairchild k-17 "arrow" type fiducial marker pattern

notch-type fiducials
^^^^^^^^^^^^^^^^^^^^

.. _fairchild t11d:

T-11
"""""

These cameras use four mid-side fiducial markers. Note that the marker location is given by a small pinhole just
inside the frame of the image, rather than the notch-shaped marks in the frame. Note also the distance between P5 and
the principal point is longer than the distance between P6 and the principal point.

.. image:: img/fairchild_t11_notch.png
    :width: 500
    :align: center
    :alt: a diagram of a fairchild t-11 notch-type fiducial marker pattern with fiducial markers labeled

|br|


T-12
"""""

Similar to the T-11 pattern shown above, the T-12 also used corner fiducial markers. Note that the fiducial marker
for each of these is a small pinhole just inside the frame, rather than the notch-shaped marks in the frame.

.. image:: img/fairchild_t12_notch.png
    :width: 500
    :align: center
    :alt: a diagram of a fairchild t-12 notch-type fiducial marker pattern with fiducial markers labeled

|br| This camera also used cross-type fiducial markers instead of dots, with the same approximate measurements:

.. image:: img/fairchild_t12_cross.png
    :width: 500
    :align: center
    :alt: a diagram of a fairchild t-12 notch-type fiducial marker pattern with fiducial markers labeled

T-5
""""

Based on diagrams like this one, which are included in some calibration reports, I am taking the location of each
fiducial marker to be the corner formed by the horizontal/vertical edge of the notch, and the image frame:

.. image:: img/report_diagram.png
    :width: 700
    :align: center
    :alt: a diagram showing where the measurements were made for estimating the fiducial marker separation

|br| If the data strip is not visible on the image, you should be able to identify P5 as being a slightly larger than
the other three notches.

.. image:: img/t5_notch.png
    :width: 500
    :align: center
    :alt: a diagram of a Fairchild T-5 with fiducial markers labeled

KC-6A
""""""

Example coming soon...


checker-type fiducials
^^^^^^^^^^^^^^^^^^^^^^

.. _fairchild t11s:

T-11
""""

Some versions of the T-11 used checker-type fiducials like the one shown below. Note the direction of flight indicator
next to the P5 marker:

.. image:: img/fairchild_t11_check.png
    :width: 500
    :align: center
    :alt: a diagram of a Fairchild T-11 with fiducial markers labeled


Park Cameras
-------------

In addition to the wing-style fiducial marker pattern, Park Aerial Surveys, Inc. manufactured and used cameras with
eight crosshair-style fiducial markers (four corner, four mid-side).

These look very similar to the Wild Heerbrugg RC8 and RC10 crosshair-style cameras shown below, but are distinguished
by the corner fiducial markers being closer to the frame (further from the principal point) and set inside of a
smaller rounded frame:

.. image:: img/park_crosshair.png
    :width: 500
    :align: center
    :alt: an example of a Park crosshair-style camera with corner and mid-side fiducial markers labeled

|br| Some, but not all, of these cameras also displayed either the camera serial number and focal length, or just the
focal length, on the side of the frame near one of the fiducial markers.

Wild Heerbrugg Cameras
----------------------

.. _wild rc5:

corner-only fiducials
^^^^^^^^^^^^^^^^^^^^^

Earlier Wild cameras, such as the RC5 or RC5A, or earlier versions of the RC8, have four crosshair fidcuial marks in
the corners, set inside of a rounded frame:

.. image:: img/wild_corner.png
    :width: 500
    :align: center
    :alt: an example of a Wild corner with corner fiducial markers labeled

|br| Some versions also used an alternate corner marker:

.. image:: img/wild_light_corner.png
    :width: 500
    :align: center
    :alt: an example of a Wild camera with corner fiducial markers labeled

|br|

+----+-----+-----+
|    | x   | y   |
+----+-----+-----+
| P1 | 212 | 0   |
+----+-----+-----+
| P2 | 0   | 212 |
+----+-----+-----+
| P3 | 0   | 0   |
+----+-----+-----+
| P4 | 212 | 212 |
+----+-----+-----+

.. _wild rc10:

corner + mid-side fiducials
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Later versions of the RC8 camera, as well as the RC10 camera, used a pattern with eight fiducial markers: four corner
markers as on the corner type, and four mid-side markers.

These markers may be all larger crosses, as in the corner marker types (these were typically RC8 cameras):

.. image:: img/wild_all_cross.png
    :width: 500
    :align: center
    :alt: an example of a Wild camera with corner and mid-side fiducial markers labeled

|br| Or, the markers might be a mix of cross-style markers in the corner, and crosshair style on the mid-side markers
(again, these are typically RC8 cameras):

.. image:: img/wild_mid_crosshair.png
    :width: 500
    :align: center
    :alt: an example of a Wild camera with corner and mid-side fiducial markers labeled

|br| Or, they might be all crosshair style markers (typically RC10 cameras):

.. image:: img/wild_all_crosshair.png
    :width: 500
    :align: center
    :alt: an example of a Wild camera with corner and mid-side fiducial markers labeled

|br| No matter the style, the approximate location of these fiducial markers is largely the same:

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

Zeiss RMK Cameras
-----------------

.. note::

    The model name for the Zeiss RMK includes information about the lens and film type. For example:

    - **Zeiss RMK 15/23** means that the camera uses a Pleogon lens with a ~15 cm (~150 mm) focal length and 23 cm film.
    - **Zeiss RMK A 15/23** means that the camera uses a Pleogon A lens with a ~15 cm focal length and 23 cm film.
    - **Zeiss RMK A 21/23** means that the camera uses a Pleogon A lens with a ~21 cm focal length and 23 cm film.
    - **Zeiss RMK A 30/23** means that the camera uses a Pleogon A lens with a ~30 cm focal length and 23 cm film

    ... and so on. This pattern doesn't necessarily tell you the pattern of the fiducial markers, but it should at
    least give you a rough idea of the focal length and film size of the camera.

mid-side only
^^^^^^^^^^^^^

Earlier models tended to use only mid-side fiducial markers:

.. image:: img/zeiss_mid.png
    :width: 500
    :align: center
    :alt: an example image taken by a Zeiss RMK with mid-side fiducial markers labeled

|br|

+----+-----+-----+
|    | x   | y   |
+----+-----+-----+
| P5 | 0   | 113 |
+----+-----+-----+
| P6 | 226 | 113 |
+----+-----+-----+
| P7 | 113 | 0   |
+----+-----+-----+
| P8 | 113 | 226 |
+----+-----+-----+

.. note::

    The coordinates above correspond to the center of the small dot near the tip of the fiducial marker.

.. _zeiss corner:

corner fiducial markers
^^^^^^^^^^^^^^^^^^^^^^^^

Later(?) versions of the Zeiss RMK camera used eight fiducial markers: four mid-side markers, and four corner markers.
These came in (at least) two main styles: a "floating style":

.. image:: img/zeiss_corner_float.png
    :width: 500
    :align: center
    :alt: an example image taken by a Zeiss RMK with corner and mid-side fiducial markers labeled

|br| and a "fixed" style:

.. image:: img/zeiss_corner.png
    :width: 500
    :align: center
    :alt: an example image taken by a Zeiss RMK with corner and mid-side fiducial markers labeled

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

references
----------

.. [1] Fleming EA (1960) Recognition of Air Survey Lens Types. *The Canadian Surveyor* 15(**2**), 91–96.
             doi: `10.1139/tcs-1960-0027 <https://doi.org/10.1139/tcs-1960-0027>`__.

.. [2] In a number of the calibration reports for these "modified" cameras, the phrase "... which has been modified to
       meet the requirements of a precision camera" is included. Presumably, this modification is adding fiducial
       markers, which might explain why this camera type corresponds to so many different fiducial marker patterns.
