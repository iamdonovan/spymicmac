from pathlib import Path
from setuptools import find_packages, setup

readme = Path(__file__).parent / 'README.md'

setup(name='spymicmac',
      version='0.1.1',
      description='a python package for processing KH-9 and historical aerial imagery using MicMac',
      long_description=readme.read_text(),
      long_description_content_type='text/markdown',
      url='https://github.com/iamdonovan/sPyMicMac',
      doc_url='https://spymicmac.readthedocs.io/',
      author='Bob McNabb',
      author_email='robertmcnabb@gmail.com',
      maintainer='iamdonovan',
      license='GPL-3.0',
      license_file='LICENSE',
      include_package_data=True,
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
          'Programming Language :: Unix Shell',
          'Topic :: Scientific/Engineering :: GIS',
          'Topic :: Scientific/Engineering :: Image Processing'
      ],
      python_requires='>=3.7',
      install_requires=['numpy', 'scipy', 'matplotlib', 'lxml', 'pybob>0.25',
                        'shapely', 'opencv-python', 'pandas', 'geopandas', 'pymmaster>0.1',
                        'scikit-image>=0.18', 'gdal>=3.2.0',
                        'sphinx-argparse', 'earthengine-api', 'pyasn1', 'usgs'],
      packages=find_packages(),
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
              'remove_crosses = spymicmac.tools.remove_crosses:main',
              'remove_measures = spymicmac.tools.remove_measures:main'
          ],
      },
      zip_safe=False)
