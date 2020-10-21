---
title: "pyEQUIB: Python Package for Plasma Diagnostics and Abundance Analysis"
tags:
  - python
  - astrophysics
  - gaseous nebulae
  - plasma diagnostics
  - abundance analysis
authors:
  - name: Ashkbiz Danehkar
    orcid: 0000-0003-4552-5997
    affiliation: "1, 2, 3"
affiliations:
 - name: Department of Physics and Astronomy, Macquarie University, Sydney, NSW 2109, Australia
   index: 1
 - name: Harvard-Smithsonian Center for Astrophysics, 60 Garden Street, Cambridge, MA 02138, USA 
   index: 2
 - name: Department of Astronomy, University of Michigan, 1085 S. University Avenue, Ann Arbor, MI 48109, USA 
   index: 3
date: 20 October 2020
bibliography: paper.bib
---

# Introduction

The interpretation of emission-line spectra from ionized gases are crucial in studying planetary nebulae and H II regions in our Galaxy and other galaxies. Physical properties of the interstellar medium such as the electron temperature and concentration can be determined from nebular emission lines. Moreover, chemical elements can be probed using the collisionally excited lines (CEL) and recombination lines (RL) emitted from ionized gases in the interstellar medium [see e.g. @Danehkar:2013; @Danehkar:2016; @Danehkar:2018a]. The interstellar extinction can also be derived from the observed Balmer decrements. 

# Statement of Need

The spectral analysis of nebular emission lines with the Python programming language requires a package that can calculate emissivities for collisional excitation and recombination processes with the given electron temperature and density, and estimate the  reddening factor for the given interstellar extinction law and physical condition. The previous package `proEQUIB` [@Danehkar:2018b] requires the Interactive Data Language (IDL) compiler. It was necessary to make a similar package purely developed in Python that can be widely used by astronomers and researchers, who prefer to conduct the data analysis of nebular spectra with the high-level programming language Python that is freely available for all operating systems.

# Description

`pyEQUIB` is a pure Python open-source package containing several application programming interface (API) functions that can be employed for plasma diagnostics and abundance analysis of nebular emission lines. This package is a Python implementation of the IDL library `proEQUIB` [@Danehkar:2018b] that is coupled to the IDL library `AtomNeb` [@Danehkar:2019]. The collisional excitation and recombination units of this package need to have the energy levels, collision strengths, transition probabilities, and recombination coefficients, which can be retrieved from the Python package `AtomNeb` for _Atomic Data of Ionized Nebulae_ [@Danehkar:2020]. The API functions of this package can be used to deduce the electron temperature, electron concentration, chemical elements from CELs and Rls, and the interstellar extinction from the Balmer decrements emitted from ionized gaseous nebulae. This package has three units, namely:

- the _collisional excitation unit_ that contains API functions for the interpretation of CELs. This unit was initially developed based on the algorithm of the FORTRAN program `EQUIB` [@Howarth:1981; @Howarth:2016], which was also incorporated into the nebular empirical analysis tool `NEAT` [@Wesson:2012]. 
The program `EQUIB` computes atomic level populations and emissivities for multi-level atoms in the equilibrium condition using collision strengths and transition probabilities for the specified temperature and density. These API functions are useful for deriving the electron temperature, electron density, and chemical abundances from CELs emitted from ionized gases in planetary nebulae and H II regions.

- the _recombination unit_ that contains API functions for the interpretation of RLs. This unit was originally developed based on the algorithm of the recombination scripts written by X. W. Liu and Y. Zhang used by the FORTRAN photoionization program `MOCASSIN` [@Ercolano:2003; @Ercolano:2005] and the FORTRAN spectral analysis program `NEAT` [@Wesson:2012], which calculates emissivities using recombination coefficients for the given physical condition. These API functions are valuable to deduce the chemical abundances from RLs emitted from gaseous nebulae.

- the _reddening unit_ that contains API functions for the extinction calculation and reddening correction. This unit was developed based on the reddening IRAF scripts included in the Space Telescope Science Data Analysis System, former `STSDAS` IRAF package [@Bushouse:1994; @Shaw:1994]. These API functions can be used to calculate the interstellar extinction from the observed Balmer decrements for different reddening laws, and perform the reddening correction on the observed fluxes of nebular emission lines. 

`pyEQUIB` can be simply used by astronomers, who are familiar with the high-level, general-purpose programming language Python. The previous IDL version `proEQUIB` [@Danehkar:2018b] has been used for studies of gaseous nebulae [@Danehkar:2016; @Danehkar:2018a]. This package requires the Python packages `NumPy` [@Walt:2011; @Harris:2020], `SciPy` [@Virtanen:2020], and `AtomNeb` [@Danehkar:2020]. The API functions of this package can be used to analyze emission-line spectra from planetary nebulae and H II regions, as well as extragalactic sources [see e.g. @Danehkar:2014; @Danehkar:2014b; @Danehkar:2016]. This package is released under the GNU General Public License. The source code is publicly available on the GitHub platform. The latest version of this package can be installed directly from its repository on the GitHub, and its stable version from the Python Package Index (PyPi) via ``pip install pyequib`` or alternatively from the Anaconda Python package distributor via ``conda install -c conda-forge pyequib``. The online documentation, tutorials and examples are available on its GitHub page (https://equib.github.io/pyEQUIB/) and its Read the Docs documentation page (https://pyequib.readthedocs.io/).

# Acknowledgements

The author acknowledges the support of Research Excellence Scholarship from Macquarie University.

# References
