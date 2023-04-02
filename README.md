[![Documentation Status](https://readthedocs.org/projects/spymicmac/badge/?version=latest)](https://spymicmac.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/spymicmac.svg)](https://badge.fury.io/py/spymicmac)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/spymicmac/badges/version.svg)](https://anaconda.org/conda-forge/spymicmac)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/spymicmac/badges/license.svg)](https://anaconda.org/conda-forge/spymicmac)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/spymicmac/badges/platforms.svg)](https://anaconda.org/conda-forge/spymicmac)


# sPyMicMac
A python package to help in the processing of KH-9 and historical aerial imagery using
[MicMac](https://micmac.ensg.eu/index.php/Accueil).

## Installing MicMac
General instructions for installing MicMac can be found via the [Installation](https://micmac.ensg.eu/index.php/Install)
tutorial found on the MicMac wiki. For specific instructions for installing the version of
MicMac used in developing sPyMicMac, check out the 
[readthedocs page](https://spymicmac.readthedocs.io/en/latest/setup.html) for sPyMicMac. 

## Installing sPyMicMac from PyPI
As of version 0.1, sPyMicMac is available via PyPI. To install the latest packaged version (0.1.1), simply run:

```sh
pip install spymicmac
```

to install sPyMicMac in your environment.

## Installing sPyMicMac from conda-forge
As of version 0.1.1, sPyMicMac is available via conda-forge. To install the latest packaged version (0.1.1), 
simply run:

```sh
conda install -c conda-forge spymicmac
```

to install spymicmac in your environment.

## Installing the latest version of sPyMicMac

To install the latest version of sPyMicMac, first clone the GitHub repository:

```sh
git clone git@github.com:iamdonovan/sPyMicMac.git
```

### Environment setup

Using [conda](https://docs.conda.io/en/latest/) and the environment.yml file from this repository:

```sh
cd path/to/repository
conda env create -f environment.yml
```

Once the environment has been created successfully (this can take some time, so grab a coffee), activate it:

```sh
conda activate spymicmac
```
To install this version of sPyMicMac as a static package, run ``pip`` without the ``-e`` flag:

```sh
# install the development version
pip install .
```

To install a development version of sPyMicMac, you can use the ``-e`` flag:
```sh
# install the development version in editing mode
pip install -e .
```

## Basic usage

For more detailed information on using sPyMicMac, see the [readthedocs page](https://spymicmac.readthedocs.io).
```python
import spymicmac.register as rt
```
