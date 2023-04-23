KH-9 pre-processing steps
==========================
There are a number of necessary pre-processing steps when working with KH-9 imagery. Distortions in the film,
caused by storage and other conditions, must be corrected to be able to accurately process the images. After
resampling, it is also possible to remove the Reseau marks from the images, in order to improve the final results.

Finally, because of the size of the film (approx. 9"x18"), the images are scanned in two halves that must be joined
together. And finally, it can also be helpful to improve contrast in the images, in order to help improve the
final results.

image joining
--------------
Because of the large size of the film, USGS scans the images in two halves with a small amount of overlap, as shown in
the example below.

.. image:: ../../img/half_a.png
    :width: 49%

.. image:: ../../img/half_b.png
    :width: 49%

|br| In ``spymicmac``, the function to join the images is :py:meth:`spymicmac.image.join_hexagon`, with a corresponding
command-line tool :doc:`../../spymicmac/scripts/join_hexagon`.

Normally, the scans are labelled 'a' and 'b', with 'a' corresponding to the left-hand scan, and 'b' corresponding to
the right-hand scan. This is what :py:meth:`spymicmac.image.join_hexagon` is expecting - that the overlap between the
two halves is the right-hand side of image 'a', and the left-hand side of image 'b'.

After calling ``join_hexagon``, the image should look something like this:

.. image:: ../../img/joined.png
    :width: 98%
    :align: center
    :alt: a re-sampled and joined KH-9 image showing Hofsjökull, Iceland

|br| As there is sometimes a difference in brightness between the two halves, :py:meth:`spymicmac.image.join_hexagon`
has the option to blend the two halves over the overlap by averaging the values from the two halves, starting from
100% of the value of image 'a', linearly increasing to 100% of the value of image 'b' at the end of the
overlapping part.


reseau field
-------------
To help correct some of the distortion in the images caused by film storage, :py:meth:`spymicmac.matching` includes
a routine to automatically find the Reseau markers in the image and use their locations to resample the images using
:py:meth:`spymicmac.resample.resample_hex`.

In the images below, you can see the difference between the expected location of each Reseau marker and the
automatically detected locations:

.. image:: ../../img/reseau_field.png
    :width: 98%
    :alt: a KH-9 image with the Reseau field warping shown

To run the routine, use either :py:meth:`spymicmac.matching.find_reseau_grid` or
:doc:`../../spymicmac/scripts/find_reseau_grid`. This will produce a ``MeasuresIm`` file that will be read by
:py:meth:`spymicmac.resample.resample_hex`.

.. note::
    Before running :py:meth:`spymicmac.resample.resample_hex`, you will also need to run
    :doc:`../../spymicmac/scripts/generate_micmac_measures` in order to generate the ``MeasuresCamera.xml`` file needed,
    then move ``MeasuresCamera.xml`` to the ``Ori-InterneScan`` directory in the correct folder.


cross removal
--------------

Once you have found the Reseau marks in each image half, you can "remove" the Reseau marks using either
:py:meth:`spymicmac.matching.remove_crosses` or :doc:`../../spymicmac/scripts/remove_crosses`.

.. image:: ../../img/fixed_cross.png
    :width: 600
    :align: center
    :alt: an image showing a Reseau mark on the left, and the Reseau mark erased on the right.

|br| After this step, you can use :doc:`../../spymicmac/scripts/resample_hexagon`.

.. note::

    Because ``mm3d ReSampFid`` calculates an affine transform based on the fiducial marker locations, it does not
    actually correct the image using the marker locations. For KH-9 Mapping Camera images, it's better to use
    :doc:`../../spymicmac/scripts/resample_hexagon`.

contrast enhancement
---------------------

Most of the scanned KH-9 images provided by USGS do not have issues with striping. However, they can still be
low contrast, and it can help to use either of :py:meth:`spymicmac.image.stretch_image` or
:py:meth:`spymicmac.image.contrast_enhance` for this.

For examples of these functions applied to a historical aerial image, see :ref:`contrast-enhancement`.
