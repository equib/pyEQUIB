"""
pyEQUIB - python package for atomic level populations and line emissivities in statistical equilibrium
"""

__all__ = []
from .version import __version__

import sys
from numpy.version import version as numpy_version

if sys.version_info[0:2] < (2, 6):
    log_.warn('pyEQUIB requires Python version >= 2.6, but it is version {0}'.format(sys.version_info), calling='pyEQUIB')
try:
    if [int(n) for n in (numpy_version.split('.')[:3])] < [1, 5, 1] :
        log_.warn('pyEQUIB Numpy version >= 1.5.1, but it is version {0}'.format(numpy_version), calling='pyEQUIB')
except:
    log_.warn('Cannot find Numpy version {0}, report the bug'.format(numpy_version), calling='pyEQUIB')
    
from pyequib import ext
from pyequib import cel
from pyequib import orl

