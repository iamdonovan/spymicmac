"""
spymicmac is a python package to help in the processing of historical aerial imagery using MicMac
"""
from importlib.metadata import version
from . import data
from . import image
from . import micmac
from . import matching
from . import orientation
from . import register
from . import resample
# from . import ee_tools


__version__ = version(__name__)

__all__ = [
    'data',
    'image',
    'micmac',
    'matching',
    'orientation',
    'register',
    'resample',
    '__version__'
]
