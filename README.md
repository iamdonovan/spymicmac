[![Documentation Status](https://readthedocs.org/projects/spymicmac/badge/?version=latest)](https://spymicmac.readthedocs.io/en/latest/?badge=latest)

# sPyMicMac
A python package to help in the processing of KH-9 and historical aerial imagery using
[MicMac](https://micmac.ensg.eu/index.php/Accueil).

## Installing MicMac
General instructions for installing MicMac can be found via the [Installation](https://micmac.ensg.eu/index.php/Install)
tutorial found on the MicMac wiki. For specific instructions for installing the version of
MicMac used in developing sPyMicMac, check out the 
[readthedocs page](https://spymicmac.readthedocs.io/en/latest/setup.html) for sPyMicMac. 

## Installing sPyMicMac

```sh
# Clone the repository
git clone git@github.com:iamdonovan/sPyMicMac.git

# install the development verion in editing mode
pip install -e [path2folder/sPyMicMac]
```

## Environment setup

Using [conda](https://docs.conda.io/en/latest/) and the environment.yml file from this repository:

```sh
conda env create -f environment.yml
```

## Basic usage

For more detailed information on using sPyMicMac, see the [readthedocs page](https://spymicmac.readthedocs.io).
```python
import sPyMicMac.register as rt
```
