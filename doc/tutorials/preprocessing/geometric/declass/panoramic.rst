panoramic cameras
=================

.. note::

    MicMac does not currently have an `optical bar <https://en.wikipedia.org/wiki/Optical_Bar_Camera>`__ camera
    available, which means that it is not currently possible to process KH-4/A/B or KH-9 panoramic camera imagery
    using MicMac.

    To process these images, you will need to install
    `NASA's Ames Stereo Pipeline (ASP) <https://stereopipeline.readthedocs.io/en/latest/introduction.html>`__, which
    can be used for further processing.

resampling using the frame
---------------------------

Unlike with the KH-9 Mapping Camera, for KH-4/A/B and KH-9 panoramic camera images, there are no réseau markers available
to help resample the images. Instead, :py:meth:`spymicmac.resample.crop_panoramic` can be used to first rotate the
image based on the detected locations of either rail marks (for KH-4 images) or "wagon wheel" marks (for KH-9 PC
images), though not all cameras will have these marks.

Then, the (rough) image border is detected, using :py:meth:`spymicmac.image.get_rough_frame`:

.. image:: img/pan_border.png
    :width: 98%
    :align: center
    :alt: a re-sampled and joined KH-4 image showing the image border outlined in red

|br| the image is then cropped to this border and, optionally, re-sampled to a smaller size:

.. image:: img/pan_cropped.png
    :width: 98%
    :align: center
    :alt: a re-sampled and joined KH-4 image with the original border removed

|br| In both cases, if the image is from the aft camera (as determined by the filename), then the image is rotated by
180° before being saved to the disk.

If you have a rough idea of where the image border is (using the left, right, top, and bottom coordinates), you can
also use :py:meth:`spymicmac.resample.crop_from_extent` to resample the image. This is especially useful for cases
where one part of the image is over dark water, for example.

next steps
-----------

Once the images have been resampled, and you have done any desired radiometric pre-processing, you can follow the
`ASP Examples <https://stereopipeline.readthedocs.io/en/latest/examples.html>`__ for the relevant image type in order
to georeference the images and extract a DEM. :py:mod:`spymicmac.asp` has a number of functions available for
generating the necessary files for ASP processing.
