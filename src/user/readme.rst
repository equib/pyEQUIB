Introduction
============

**pyEQUIB** package is a collection of `Python <https://www.python.org/>`_ programs developed to perform plasma diagnostics and abundance analysis using emission line fluxes measured in ionzed nebulae. It uses the `AtomNeb Python Package <https://github.com/atomneb/AtomNeb-py>`_ to read collision strengths and transition probabilities for collisionally excited lines (CEL), and recombination coefficients for recombination lines (RL). This Python package can be used to determine interstellar extinctions, electron temperatures, electron densities, and ionic abundances from the measured fluxes of emission lines. It mainly contains the follwing API functions written purely in Python: 

Collisional Excitation Unit
---------------------------

**API functions for collisionally excited lines (CEL)** have been developed based on the algorithm of the FORTRAN program `EQUIB <http://adsabs.harvard.edu/abs/2016ascl.soft03005H>`_ written in FORTRAN by `Howarth & Adams (1981) <http://adsabs.harvard.edu/abs/1981ucl..rept.....H>`_. The program EQUIB calculates atomic level populations and line emissivities in statistical equilibrium in multi-level atoms for different physical conditions of the stratification layers where the chemical elements are ionized. Using the Python implementation of the program `EQUIB <http://adsabs.harvard.edu/abs/2016ascl.soft03005H>`_, electron temperatures, electron densities, and ionic abundances are determined from the measured fluxes of collisionally excited lines.

Recombination Unit
------------------

**API functions for recombination lines (RL)** have been developed based on the algorithm of the recombination scripts by X. W. Liu and Y. Zhang included in the FORTRAN program `MOCASSIN <https://github.com/mocassin/MOCASSIN-2.0>`_. These API functiosn are used to determine ionic abundances from recombination lines for some heavy element ions.
 
Reddening Unit
--------------

**API functions for reddening and extinctions** have been developed according to the methods of the reddening law functions from `STSDAS IRAF Package <http://www.stsci.edu/institute/software_hardware/stsdas>`_, which are used to obtain interstellar extinctions and deredden measured fluxes based on different reddening laws.

