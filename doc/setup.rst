Installation and Setup
=======================

The following is a (non-exhaustive) set of instructions for getting setup to run MMASTER on your own machine. Note
that this can be an **extremely** computationally intensive process, so we don't really recommend trying to run this on your
personal laptop.

As this is a (non-exhaustive) set of instructions, it may not work 100% with your particular setup.
We are happy to try to provide guidance/support, **but we make no promises**.

Installing MicMac
#################

Detailed installation instructions for MicMac on multiple platforms can be found `here <https://micmac.ensg.eu/index.php/Install/>`_,
but we've added a short summary to help guide through the process.

First, clone the MicMac repository to a folder on your computer (you can also do this online via github):
::

    /home/bob/software:~$ git clone https://github.com/micmacIGN/micmac.git
    ...
    /home/bob/software:~$ cd micmac
    /home/bob/software/micmac:~$ git fetch
    /home/bob/software/micmac:~$ git checkout IncludeALGLIB

This will clone the MicMac git repository to your machine, fetch the remote, and switch to the *IncludeALGLIB* branch.
Check the **README.md** (or **LISEZMOI.md**) file to install any dependencies, then:
::

    /home/bob/software/micmac:~$ mkdir build && cd build/
    /home/bob/software/micmac/build:~$ cmake .. -DWITH_QT5=1 -DWERROR=0 -DWITH_CCACHE=OFF
    ...
    /home/bob/software/micmac/build:~$ make install -j$n

where $n is the number of cores to compile MicMac with. The compiler flag **-DWERROR=0** is needed, as some of the dependencies
will throw warnings that will force the compiler to quit with errors if we don't turn it off.

Finally, make sure to add the MicMac bin directory (/home/bob/software/micmac/bin in the above example) to your $PATH
environment variable, in order to be able to run MicMac. You can check that all dependencies are installed by running
the following:
::

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

Optional: Preparing a python environment
########################################
If you like, you can set up a dedicated python environment for your sPyMicMac needs. This can be handy, in case any
packages required by sPyMicMac clash with packages in your default environment. Our personal preference is `conda <https://docs.conda.io/en/latest/>`_,
but your preferences may differ.

The git repository has a file, environment.yml, which provides a working environment for sPyMicMac and conda.
Once you have conda installed, simply run:
::

    conda env create -f environment.yml

This will create a new conda environment, called spymicmac, which will have all of the various python packages
necessary to run sPyMicMac. To activate the new environment, type:
::

    conda activate spymicmac

And you should be ready to go. Note that you will have to activate this environment any time you wish to run
sPyMicMac scripts and tools, if it is not already activated in your terminal.

Installing sPyMicMac
####################
Next, use **pip** to install the scripts and python modules:
::

    pip install -e sPyMicMac

from the repository folder. Note: the *-e* allows you to make changes to the code (for example, from git updates
or through your own tinkering), that will then be updated within your python install. If you run *pip install*
without this option, it will install a static version of the package, and any changes/updates will have to be
explictly re-installed.

Assuming that you haven't run into any errors, you should be set up. You can verify this by running:
::

    register_ortho.py -h

From the command line (in a non-Windows environment; *Windows instructions coming soon-ish*).

You should see the following output (or something very similar):

