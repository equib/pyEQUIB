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
          description = 'pyequib: Python Package for Plasma Diagnostics and Abundance Analysis',
          long_description=fd.read(),
          author='Ashkbiz Danehkar',
          author_email='ashkbiz.danehkar@students.mq.edu.au',
          url='https://equib.github.io/pyEQUIB/',
          download_url = 'https://github.com/equib/pyEQUIB',
          keywords = ['pyEQUIB', 'EQUIB', 'equilibrium', 'emissivities', 'plasma diagnostics', 'abundance analysis', 'astronomy', 'physics', 'astrophysics'],
          license='http://www.gnu.org/licenses/gpl.html',
          platforms=['any'],
          packages=['pyequib'],
          package_data={'pyequib': ['*.txt', 'text/*.readme', 'atomic-data/*.dat', 'line-list/*.dat']},
          install_requires=['numpy','scipy','atomneb'],
         )

