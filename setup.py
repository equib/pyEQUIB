#!/usr/bin/env python

#
# Setup script for pyequib
#

import os
import codecs
try:
      from setuptools import setup
except ImportError:
      from distutils.core import setup

import pyequib


with codecs.open('README.rst', 'r', 'utf-8') as fd:
    setup(name='pyequib',
          version=pyequib.__version__,
          description = 'pyEQUIB: Python package for equilibrium atomic populations and line emissivities',
          long_description=fd.read(),
          author='Ashkbiz Danehkar',
          author_email='ashkbiz.danehkar@students.mq.edu.au',
          url='http://physics.mq.edu.au/~ashkbiz/pyequib/',
          download_url = 'https://github.com/equib/pyEQUIB',
          keywords = ['pyEQUIB', 'EQUIB', 'equilibrium', 'emissivities', 'atomic', 'astronomy', 'physics', 'astrophysics'],
          license='http://www.gnu.org/licenses/gpl.html',
          platforms=['any'],
          packages=['pyequib'],
          package_data={'pyequib': ['*.txt', 'text/*.readme', 'atomic-data/*.dat']},
          install_requires=['numpy','scipy'],
         )

