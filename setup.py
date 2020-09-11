from setuptools import setup

setup(name='sPyMicMac',
      version='0.1',
      description='a python package for processing KH-9 imagery using MicMac',
      url='http://github.com/iamdonovan/sPyMicMac',
      author='Bob McNabb',
      author_email='robertmcnabb@gmail.com',
      license='GPL-3.0',
      packages=['sPyMicMac'],
      install_requires=['numpy', 'scipy', 'matplotlib', 'fiona', 'pyvips', 'lxml', 'pybob>=0.25',
                        'shapely', 'opencv-python', 'pandas', 'geopandas', 'pymmaster',
                        'scikit-image', 'gdal', 'h5py', 'pyproj', 'llc', 'numba', 'descartes',
                        'sphinx-argparse', 'earthengine-api', 'pyasn1', 'usgs'],
      scripts=['bin/combine_auto_measures.py', 'bin/get_autogcp_locations.sh', 'bin/find_reseau_shifts.py',
               'bin/join_balance_images.py', 'bin/join_hexagon_halves.py', 'bin/remove_crosses.py',
               'bin/register_ortho.py', 'bin/generate_micmac_measures.py', 'bin/move_bad_tapas.py'],
      zip_safe=False)
