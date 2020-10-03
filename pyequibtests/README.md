## pyEQUIB
[![PyPI version](https://badge.fury.io/py/pyequib.svg)](https://badge.fury.io/py/pyequib)
[![Build Status](https://travis-ci.org/equib/pyEQUIB.svg?branch=master)](https://travis-ci.org/equib/pyEQUIB)
[![Coverage Status](https://coveralls.io/repos/github/equib/pyEQUIB/badge.svg?)](https://coveralls.io/github/equib/pyEQUIB?branch=master)
[![GitHub license](https://img.shields.io/aur/license/yaourt.svg)](https://github.com/equib/pyEQUIB/blob/master/LICENSE)

Python package for atomic level populations and line emissivities in statistical equilibrium

### Description
**pyEQUIB** is a collection of [Python](https://www.python.org/) programs developed to perform plasma diagnostics and abundance analysis using emission line fluxes measured in ionzed nebulae. It uses the [AtomNeb Python Package](https://github.com/atomneb/AtomNeb-py) to read collision strengths and transition probabilities for collisionally excited lines (CEL), and recombination coefficients for recombination lines (RL). This Python package can be used to determine interstellar extinctions, electron temperatures, electron densities, and ionic abundances from the measured fluxes of emission lines.

### Installation
To install the last version, all you should need to do is

    python setup.py install

To install the stable version, you can use the preferred installer program (pip):

    pip install pyequib

### References

* Danehkar, A. (2018). proEQUIB: IDL Library for Plasma Diagnostics and Abundance Analysis. *J. Open Source Softw.*, **3**, 899. doi:[10.21105/joss.00899](https://doi.org/10.21105/joss.00899) ads:[2018JOSS....3..899D](https://ui.adsabs.harvard.edu/abs/2018JOSS....3..899D).

