from setuptools import setup

setup(name='spymicmac',
      version='0.1.1',
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
      scripts=['bin/get_autogcp_locations.sh',
               'bin/post_process_micmac.sh',
               'bin/resample_hexagon.sh'],
      entry_points={
          'console_scripts': [
              'balance_images = spymicmac.tools.balance_images:main',
              'combine_auto_measures = spymicmac.tools.combine_auto_measures:main',
              'find_reseau_grid = spymicmac.tools.find_reseau_grid:main',
              'generate_micmac_measures = spymicmac.tools.generate_micmac_measures:main',
              'join_hexagon = spymicmac.tools.join_hexagon:main',
              'move_bad_tapas = spymicmac.tools.move_bad_tapas:main',
              'register_ortho = spymicmac.tools.register_ortho:main',
              'remove_crosses = spymicmac.tools.remove_crosses:main'
          ],
      },
      zip_safe=False)
