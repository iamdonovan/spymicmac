Installation and Setup
=======================

The following is a (non-exhaustive) set of instructions for getting setup to run spymicmac on your own machine. Note
that this can be an **extremely** computationally intensive process, so we don't really recommend trying to run this on your
personal laptop.

As this is a (non-exhaustive) set of instructions, it may not work 100% with your particular setup.
We are happy to try to provide guidance/support, **but we make no promises**.

Installing MicMac
#################

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

Installing spymicmac
#####################
spymicmac is available in a number of ways - either installing from source or packaged via PyPI or conda-forge.

via PyPI
------------
As of version 0.1, spymicmac is available via PyPI. To install the latest packaged version into your python environment,
simply run:

.. code-block:: sh

    pip install spymicmac

via conda-forge
-----------------
As of version 0.1.1, spymicmac is available via conda-forge. To install the latest version, run:

.. code-block:: sh

    conda install -c conda-forge spymicmac


from source
-------------
To get started, clone the repository, then navigate to the directory where the repository is downloaded:

.. code-block:: sh

    git clone https://github.com/iamdonovan/spymicmac.git

Optional: Preparing a python environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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

Installing via pip
^^^^^^^^^^^^^^^^^^^^
Once you have the environment prepared (or not), run pip from inside the ``spymicmac`` directory:

.. code-block:: sh

    pip install .

Alternatively, you can install a development version, which allows you to make changes to the code (either via git updates
or your own tinkering) without having to re-install each time. To install a development version, use the ``-e`` option:

.. code-block:: sh

    pip install -e .

Checking the installation
--------------------------
Assuming that you haven't run into any errors, you should be set up. You can verify this by running:

.. code-block:: sh

    register_ortho -h

From the command line. You should see the following output (or something very similar):

.. code-block:: text

    usage: register_ortho [-h] [-glacmask GLACMASK] [-landmask LANDMASK] [-footprints FOOTPRINTS]
                          [-im_subset IM_SUBSET [IM_SUBSET ...]] [-b BLOCK_NUM] [-ori ORI]
                          [-ortho_res ORTHO_RES] [-imgsource IMGSOURCE] [-density DENSITY]
                          fn_ortho fn_ref fn_dem fn_reldem

    Register a relative orthoimage and DEM to a reference orthorectified image and DEM.

    positional arguments:
      fn_ortho              non-referenced orthophoto mosaic
      fn_ref                georeferenced satellite image
      fn_dem                dem
      fn_reldem             relative dem corresponding to ortho

    optional arguments:
      -h, --help            show this help message and exit
      -glacmask GLACMASK    path to shapefile of glacier outlines (i.e., an exclusion mask)
      -landmask LANDMASK    path to shapefile of land outlines (i.e., an inclusion mask)
      -footprints FOOTPRINTS
                            path to shapefile of image outlines. If not set, will download from USGS.
      -im_subset IM_SUBSET [IM_SUBSET ...]
                            subset of raw images to work with (default all)
      -b BLOCK_NUM, --block_num BLOCK_NUM
                            Block number to use if multiple image blocks exist in directory.
      -ori ORI              name of orientation directory (after Ori-) [Relative]
      -ortho_res ORTHO_RES  approx. ground sampling distance (pixel resolution) of ortho image. [8 m]
      -imgsource IMGSOURCE  USGS dataset name for images [DECLASSII]
      -density DENSITY      pixel spacing to look for GCPs [200]
      -no_allfree           run Campari with AllFree set to False


