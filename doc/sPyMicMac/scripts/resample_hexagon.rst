resample_hexagon.sh
=================================

``resample_hexagon.sh`` is a shell script that will run ``mm3d ReSampFid`` on the two scanned halves of a KH-9
Hexagon image, leaving a 2 mm overlap between the two images to aid in using :py:meth`sPyMicMac.image.join_hexagon`
to join the images together.
::

    Resample KH-9 halves using mm3d ReSampFid.
    usage: resamp_hexagon.sh -r 'RES'
        -r RES      : Scan resolution to re-sample to (default: 0.014 mm/pixel)
        -h          : displays this message and exits.
