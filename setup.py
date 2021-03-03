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

#import pyequib

with codecs.open('README.rst', 'r', 'utf-8') as fd:
    setup(name='pyequib',
          version="0.4.2", #pyequib.__version__,
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
          #package_data={'pyequib': ['*.txt', 'text/*.readme', 'atomic-data/*.dat', 'line-list/*.dat']},
          data_files = [('', ['LICENSE', 
                              'atomic-data/cel_diagnostic.dat', 
                              'atomic-data/cel_diagnostic2.dat',
                              'atomic-data/cel_list.dat',
                              'atomic-data/ext_list.dat',
                              'atomic-data/hei_diagnostic.dat',
                              'atomic-data/orl_diagnostic.dat',
                              'atomic-data/orl_list.dat',
                              'line-list/cel_list_doc.dat'])],
          install_requires=['numpy','scipy','atomneb'],
         )
