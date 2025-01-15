"""
spymicmac is a python package to help in the processing of historical aerial imagery using MicMac
"""
from importlib.metadata import version
from . import asp
from . import data
from . import declass
from . import image
from . import matching
from . import micmac
from . import orientation
from . import register
from . import resample
# from . import ee_tools


__version__ = version(__name__)

__all__ = [
    'asp',
    'data',
    'declass',
    'image',
    'matching',
    'micmac',
    'orientation',
    'preprocessing',
    'register',
    'resample',
    '__version__'
]
