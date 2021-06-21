from setuptools import setup

setup(name='spymicmac',
      version='0.1',
      description='a python package for processing KH-9 imagery using MicMac',
      url='http://github.com/iamdonovan/sPyMicMac',
      author='Bob McNabb',
      author_email='robertmcnabb@gmail.com',
      license='GPL-3.0',
      packages=['spymicmac'],
      install_requires=['numpy', 'scipy', 'matplotlib', 'lxml', 'pybob>0.25',
                        'shapely', 'opencv-python', 'pandas', 'geopandas', 'pymmaster>0.1',
                        'scikit-image>=0.18', 'gdal', 'llc',
                        'sphinx-argparse', 'earthengine-api', 'pyasn1', 'usgs'],
      scripts=['bin/balance_images.py', 'bin/combine_auto_measures.py', 'bin/find_reseau_grid.py',
               'bin/generate_micmac_measures.py', 'bin/get_autogcp_locations.sh',
               'bin/join_hexagon.py', 'bin/post_process_micmac.sh', 'bin/register_ortho.py',
               'bin/remove_crosses.py', 'bin/resample_hexagon.sh'],
      zip_safe=False)
