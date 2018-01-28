=======
pyEQUIB
=======

.. image:: https://img.shields.io/pypi/v/pyequib.svg?style=flat
    :target: https://pypi.python.org/pypi/pyequib/
    :alt: PyPI Version
    
.. image:: https://travis-ci.org/equib/pyEQUIB.svg?branch=master
    :target: https://travis-ci.org/equib/pyEQUIB
    :alt: Build Status

.. image:: https://ci.appveyor.com/api/projects/status/b3gw6vgf8s0vu8nv?svg=true
    :target: https://ci.appveyor.com/project/danehkar/pyequib
    :alt: Build Status

.. image:: https://coveralls.io/repos/github/equib/pyEQUIB/badge.svg?
    :target: https://coveralls.io/github/equib/pyEQUIB?branch=master
    :alt: Coverage Status

.. image:: http://mybinder.org/badge.svg
    :target: http://mybinder.org/repo/equib/pyequib
    :alt: Binder

.. image:: https://img.shields.io/aur/license/yaourt.svg
    :target: https://github.com/equib/pyEQUIB/blob/master/LICENSE
    :alt: GitHub license

.. image:: https://img.shields.io/badge/python-2.7%2C%203.4%2C%203.5-blue.svg
    :alt: Support Python versions 2.7, 3.4 and 3.5

Description
============

The **pyEQUIB** library is a collection of `Python <https://www.python.org/>`_ programs developed to calculate atomic level populations and line emissivities in statistical equilibrium in multi-level atoms for different physical conditions of the stratification layers where the chemical elements are ionized. This library includes the Python implementation of the program `EQUIB <http://adsabs.harvard.edu/abs/2016ascl.soft03005H>`_, which was originally written in FORTRAN by `Howarth & Adams (1981) <http://adsabs.harvard.edu/abs/1981ucl..rept.....H>`_, and was recently converted to Python. It also includes the Python implementation of deredden functions from `STSDAS IRAF Package <http://www.stsci.edu/institute/software_hardware/stsdas>`_, hydrogen emissivities (`Osterbrock & Ferland 2006 <http://adsabs.harvard.edu/abs/2006agna.book.....O>`_; `Storey & Hummer 1995 <http://adsabs.harvard.edu/abs/1995yCat.6064....0S>`_), helium emissivities (`Smits 1996 <http://adsabs.harvard.edu/abs/1996MNRAS.278..683S>`_; `Porter et al. 2012 <http://adsabs.harvard.edu/abs/2012MNRAS.425L..28P>`_). It uses the Python implementation of the `MIDAS <http://www.eso.org/~ohainaut/ccd/midas.html>`_ scripts by X. W. Liu to read heavy element recombination emissivities of C II (`Davey et al. 2000  <http://adsabs.harvard.edu/abs/2000A%26AS..142...85D>`_), N II (`Escalante & Victor 1990 <http://adsabs.harvard.edu/abs/1990ApJS...73..513E>`_), O II (`Storey 1994 <http://adsabs.harvard.edu/abs/1994A%26A...282..999S>`_; `Liu et al. 1995 <http://adsabs.harvard.edu/abs/1995MNRAS.272..369L>`_), Ne II (`Kisielius et al. 1998 <http://adsabs.harvard.edu/abs/1998A%26AS..133..257K>`_), and C III and N III (`Pequignot et al. 1991 <http://adsabs.harvard.edu/abs/1991A%26A...251..680P>`_).

History of codes can be found `here <https://physics.mq.edu.au/~ashkbiz/proequib/history/>`_.

Website: `physics.mq.edu.au/~ashkbiz/pyequib <https://physics.mq.edu.au/~ashkbiz/pyequib/>`_

Installation
============

To install the last version, all you should need to do is

.. code-block::

    $ python setup.py install

To install the stable version, you can use the preferred installer program (pip):

.. code-block::

    $ pip install pyequib


References
==========
A Danehkar, Q A Parker & W Steffen, `AJ, 151, 38, 2016 <http://adsabs.harvard.edu/abs/2016AJ....151...38D>`_

A Danehkar, H Todt, B Ercolano & A Y Kniazev, `MNRAS, 439, 3605, 2014 <http://adsabs.harvard.edu/abs/2014MNRAS.439.3605D>`_

A Danehkar, Q A Parker & B Ercolano, `MNRAS, 434, 1513, 2013 <http://adsabs.harvard.edu/abs/2013MNRAS.434.1513D>`_

A Danehkar, PhD Thesis, `Macquarie University, 2014 <http://adsabs.harvard.edu/abs/2014PhDT........76D>`_
