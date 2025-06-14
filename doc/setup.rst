installation and setup
=======================

The following is a (non-exhaustive) set of instructions for getting setup to run spymicmac on your own machine. Note
that this can be an **extremely** computationally intensive process, so we don't really recommend trying to run this on your
personal laptop.

As this is a (non-exhaustive) set of instructions, it may not work 100% with your particular setup.
We are happy to try to provide guidance/support, **but we make no promises**.

installing MicMac
-----------------

Detailed installation instructions for MicMac on multiple platforms can be found `here <https://micmac.ensg.eu/index.php/Install/>`_,
but we've added a short summary to help guide through the process.

First, clone the MicMac repository to a folder on your computer (you can also do this online via github):

.. code-block:: sh

    /home/bob/software:~$ git clone https://github.com/micmacIGN/micmac.git
    ...
    /home/bob/software:~$ cd micmac
    /home/bob/software/micmac:~$ git fetch
    /home/bob/software/micmac:~$ git checkout IncludeALGLIB

This will clone the MicMac git repository to your machine, fetch the remote, and switch to the *IncludeALGLIB* branch.
Check the **README.md** (or **LISEZMOI.md**) file to install any dependencies, then:

.. code-block:: sh

    /home/bob/software/micmac:~$ mkdir build && cd build/
    /home/bob/software/micmac/build:~$ cmake .. -DWITH_QT5=1 -DWERROR=0 -DWITH_CCACHE=OFF
    ...
    /home/bob/software/micmac/build:~$ make install -j$n

where ``$n`` is the number of cores to compile MicMac with. The compiler flag ``-DWERROR=0`` is needed, as some of the dependencies
will throw warnings that will force the compiler to quit with errors if we don't turn it off.

Finally, make sure to add the MicMac bin directory (``/home/bob/software/micmac/bin`` in the above example)
to your ``$PATH`` environment variable, in order to be able to run MicMac. You can check that all dependencies are
installed by running the following:

.. code-block:: text

    /home/bob:~$ mm3d CheckDependencies
    git revision : v1.0.beta13-844-g21d990533

    byte order   : little-endian
    address size : 64 bits

    micmac directory : [/home/bob/software/micmac/]
    auxilary tools directory : [/home/bob/software/micmac/binaire-aux/linux/]

    --- Qt enabled : 5.9.5
        library path:  [/home/bob/miniconda3/envs/bobtools/plugins]

    make:  found (/usr/bin/make)
    exiftool:  found (/usr/bin/exiftool)
    exiv2:  found (/usr/bin/exiv2)
    convert:  found (/usr/bin/convert)
    proj:  found (/usr/bin/proj)
    cs2cs:  found (/usr/bin/cs2cs

In a nutshell, the basic idea is: clone the MicMac git repository, then build the source code. Simple!

installing spymicmac
--------------------
spymicmac is available in a number of ways - either installing from source or packaged via PyPI or conda-forge.

via pip
^^^^^^^^

As of version 0.1, spymicmac is available via PyPI. To install the latest packaged version into your python environment,
simply run:

.. code-block:: sh

    pip install spymicmac

via conda-forge
^^^^^^^^^^^^^^^

As of version 0.1.1, spymicmac is available via conda-forge. To install the latest version, run:

.. code-block:: sh

    conda install -c conda-forge spymicmac


from source
^^^^^^^^^^^

To get started, clone the repository, then navigate to the directory where the repository is downloaded:

.. code-block:: sh

    git clone https://github.com/iamdonovan/spymicmac.git

optional: Preparing a python environment
""""""""""""""""""""""""""""""""""""""""
If you like, you can set up a dedicated python environment for your spymicmac needs. This can be handy, in case any
packages required by spymicmac clash with packages in your default environment. Our personal preference
is `conda <https://docs.conda.io/en/latest/>`_, but your preferences may differ.

The git repository has a file, environment.yml, which provides a working environment for spymicmac and conda.
Once you have conda installed, simply run:

.. code-block:: sh

    conda env create -f environment.yml

This will create a new conda environment, called spymicmac, which will have all of the various python packages
necessary to run spymicmac. To activate the new environment, type:

.. code-block:: sh

    conda activate spymicmac

And you should be ready to go. Note that you will have to activate this environment any time you wish to run
spymicmac scripts and tools, if it is not already activated in your terminal.

installing via pip
""""""""""""""""""
Once you have the environment prepared (or not), run pip from inside the ``spymicmac`` directory:

.. code-block:: sh

    pip install .

Alternatively, you can install a development version, which allows you to make changes to the code (either via git updates
or your own tinkering) without having to re-install each time. To install a development version, use the ``-e`` option:

.. code-block:: sh

    pip install -e .

checking the installation
^^^^^^^^^^^^^^^^^^^^^^^^^

Assuming that you haven't run into any errors, you should be set up. You can verify this by running:

.. code-block:: sh

    register_relative -h

From the command line. You should see the following output (or something very similar):

.. code-block:: text

    usage: register_relative [-h] [-ort FN_ORTHO] [-ref FN_REF] [-glacmask GLACMASK] [-landmask LANDMASK]
                             [-footprints FOOTPRINTS] [-im_subset IM_SUBSET [IM_SUBSET ...]] [-b BLOCK_NUM]
                             [--subscript SUBSCRIPT] [-ori ORI] [-ortho_res ORTHO_RES] [-imgsource IMGSOURCE]
                             [-strategy STRATEGY] [-density DENSITY] [-no_allfree] [-useortho] [-max_iter MAX_ITER]
                             [-use_cps] [-cp_frac CP_FRAC] [-o] [-fn_gcps FN_GCPS]
                             dirmec fn_dem

    Register a relative DEM or orthoimage to a reference DEM and/or orthorectified image.

    positional arguments:
      dirmec                the name of the MEC directory to read the relative DEM from (e.g., MEC-Relative)
      fn_dem                path to reference DEM

    options:
      -h, --help            show this help message and exit
      -ort FN_ORTHO, --fn_ortho FN_ORTHO
                            path to relative orthoimage (optional)
      -ref FN_REF, --fn_ref FN_REF
                            path to reference orthorectified image (optional)
      -glacmask GLACMASK    path to shapefile of glacier outlines (i.e., an exclusion mask)
      -landmask LANDMASK    path to shapefile of land outlines (i.e., an inclusion mask)
      -footprints FOOTPRINTS
                            path to shapefile of image outlines. If not set, will attempt to download from USGS.
      -im_subset IM_SUBSET [IM_SUBSET ...]
                            subset of raw images to work with (default all)
      -b BLOCK_NUM, --block_num BLOCK_NUM
                            Block number to use if multiple image blocks exist in directory.
      --subscript SUBSCRIPT
                            Optional subscript to add to filenames.
      -ori ORI              name of orientation directory (after Ori-) [Relative]
      -ortho_res ORTHO_RES  approx. ground sampling distance (pixel resolution) of ortho image. [8 m]
      -imgsource IMGSOURCE  USGS dataset name for images [DECLASSII]
      -strategy STRATEGY    strategy for generating GCPs. Must be one of 'grid', 'random', or 'chebyshev' [grid]
      -density DENSITY      pixel spacing to look for GCPs [200]
      -no_allfree           run Campari with AllFree set to False
      -useortho             use the orthomosaic in Ortho-{dirmec} rather than the DEM [False]. If fn_ortho is set,
                            uses that file instead.
      -max_iter MAX_ITER    the maximum number of Campari iterations to run [5]
      -use_cps              split the GCPs into GCPs and CPs, to quantify the uncertainty of the camera model [False]
      -cp_frac CP_FRAC      the fraction of GCPs to use when splitting into GCPs and CPs [0.2]
      -o, --use_orb         use skimage.feature.ORB to identify GCP locations in the reference image (default: use
                            regular grid for matching)
      -fn_gcps FN_GCPS      (optional) shapefile or CSV of GCP coordinates to use. Column names should be [(name | id),
                            (z | elevation), x, y]. If CSV is used, x,y should have the same CRS as the reference image.

.. _usgs_setup:

using the USGS M2M API
----------------------

:py:mod:`spymicmac.data` is set up to provide a way to search for image footprints using the USGS M2M API, through
the `usgs <http://kapadia.github.io/usgs/>`__ python package.

.. note::

    Because of recent changes to the USGS API, you will need to install ``usgs`` from github
    (https://github.com/kapadia/usgs), rather than from PyPI.

In order for this to work, you will need to do the following things:

- create a free USGS EarthExplorer account: https://earthexplorer.usgs.gov/
- add your EarthExplorer login credentials to a ``.netrc`` (or ``_netrc``) file in your ``$home`` directory.
- create a USGS M2M access token (https://m2m.cr.usgs.gov/), and save the token to a file, ``.usgs_token``, in your
  ``$home`` directory.

Once you have done this, you can check access:

.. code-block:: python

    from spymicmac import data
    data._authenticate()

If this prints something like the following:

.. code-block:: text

    {'requestId': 2011843395,
     'version': 'stable',
     'data': 'eyJjaWQiOjMyNDc5MCwicyI6IjE3NDg5NDQwNDYiLCJyIjo2MTgsInAiOlsidXNlciJdfQ==',
     'errorCode': None,
     'errorMessage': None,
     'sessionId': 316673291}

You have successfully authenticated.

On the other hand, if you see something like:

.. code-block:: python

    Traceback (most recent call last):
      File "<string>", line 1, in <module>
      File "/home/bob/software/spymicmac/src/spymicmac/data.py", line 60, in _authenticate
        login = api.login(user, token)
                ^^^^^^^^^^^^^^^^^^^^^^
      File "/home/bob/software/usgs/usgs/api.py", line 159, in login
        _check_for_usgs_error(response)
      File "/home/bob/software/usgs/usgs/api.py", line 38, in _check_for_usgs_error
        raise USGSError('%s: %s' % (error_code, error))
    usgs.USGSError: AUTH_INVALID: User credential verification failed

it means that you have not successfully authenticated, and will need to check that your credentials are correct and
in the right place.
