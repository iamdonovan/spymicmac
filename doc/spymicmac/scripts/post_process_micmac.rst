post_process_micmac.sh
=================================

``post_process_micmac.sh`` is a shell script that re-names the MicMac-produced DEM and correlation mask, georeferences
the correlation mask, and applies the AutoMask (denoting areas visible in the images) to the DEM.

The results will be written in a directory one level above the files called "post_processed".
::

    Post-process outputs from MicMac into nice files to use.
    usage: PostProcessMicMac.sh -z 'UTMZONE' -n 'NAME' -l 'LEVEL' -h
        -z UTMZONE  : UTM Zone of area of interest. Takes form 'NN +north(south)'
        -n NAME     : Output name to use for the DEM and related files.
        -l LEVEL    : MicMac processing level to use (default 9)
        -d DIR      : Directory to process (default: current directory)
        -h          : displays this message and exits.
