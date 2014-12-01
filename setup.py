# coding: utf-8
import os
import sys
from setuptools import setup
from setuptools.command.test import test as TestCommand

setup(name='gisrastertools',
      version='0.1',
      description='Tools for working with Geographical Information System Rasters',
      url='http://github.com/ozak/gis-raster-tools',
      author='Ömer Özak',
      author_email='omer@omerozak.com',
      license='GPLv3',
      #package_dir={'': 'src'},
      packages=['gisrastertools'],
      install_requires=[
          'shapely',
          'numpy',
          'GDAL',
          'docopt',
          'pandas', 
          'pyproj',
          'scikit-image',
      ],
      classifiers=[
          "Development Status :: 1 - Planning",
          "Topic :: Utilities",
          "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering :: GIS',
      ],
      zip_safe=False)
